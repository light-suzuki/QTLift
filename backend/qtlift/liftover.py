from __future__ import annotations

import gzip
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
from pathlib import Path

from .models import Interval


def _fingerprint(path: Path) -> str:
    stat = path.stat()
    return f"{path.resolve()}|{stat.st_size}|{stat.st_mtime_ns}"


def cache_path(source_fasta: str | Path, target_fasta: str | Path, cache_root: str | Path) -> Path:
    source, target = Path(source_fasta), Path(target_fasta)
    key = hashlib.sha256((_fingerprint(source) + "\n" + _fingerprint(target)).encode()).hexdigest()[:20]
    return Path(cache_root) / f"{source.stem}__to__{target.stem}__{key}.paf.gz"


def _chunk_query_fasta(source: Path, cache_root: Path, chunk_size: int = 20_000_000,
                       distro: str | None = None) -> Path:
    """Split very large query contigs so minimap2 cannot exhaust WSL memory on one chromosome."""
    key = hashlib.sha256(_fingerprint(source).encode()).hexdigest()[:16]
    output = cache_root / f"{source.stem}__chunks-{chunk_size}__{key}.fa"
    if output.is_file() and output.stat().st_size > 0:
        return output
    tmp = output.with_suffix(output.suffix + f".{os.getpid()}.tmp")
    if shutil.which("wsl.exe"):
        command = ["wsl.exe"]
        if distro:
            command += ["-d", distro]
        helper = Path(__file__).resolve().parents[2] / "scripts" / "chunk_fasta.py"
        command += ["python3", _wsl_path(helper, distro), _wsl_path(source, distro), _wsl_path(tmp, distro), str(chunk_size)]
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if result.returncode:
            tmp.unlink(missing_ok=True)
            raise RuntimeError(result.stderr.decode("utf-8", errors="replace").strip() or "FASTA chunking failed")
        tmp.replace(output)
        return output
    opener = gzip.open if source.suffix.lower() == ".gz" else open
    with opener(source, "rt", encoding="utf-8") as src, tmp.open("wb") as dst:
        name, sequence, offset = None, bytearray(), 0
        def flush_chunks(final: bool = False) -> None:
            nonlocal sequence, offset
            while len(sequence) >= chunk_size or (final and sequence):
                size = min(chunk_size, len(sequence))
                dst.write(f">QTLIFT|{name}|{offset}\n".encode())
                chunk = bytes(sequence[:size])
                dst.write(chunk + b"\n")
                del sequence[:size]
                offset += size
        for line in src:
            if line.startswith(">"):
                if name is not None:
                    flush_chunks(True)
                name, sequence, offset = line[1:].split()[0], bytearray(), 0
            else:
                sequence.extend(line.strip().upper().encode("ascii"))
                flush_chunks(False)
        if name is not None:
            flush_chunks(True)
    tmp.replace(output)
    return output


def _wsl_path(path: Path, distro: str | None = None) -> str:
    raw = str(path)
    unc_match = re.match(r"^\\\\wsl(?:\.localhost)?\\([^\\]+)\\(.+)$", raw, re.IGNORECASE)
    if unc_match:
        requested_distro, linux_path = unc_match.groups()
        if distro and requested_distro.casefold() != distro.casefold():
            raise RuntimeError(f"WSL path belongs to {requested_distro}, not {distro}")
        return "/" + linux_path.replace("\\", "/")
    drive_match = re.match(r"^([A-Za-z]):[\\/](.+)$", raw)
    if drive_match:
        drive, relative = drive_match.groups()
        return f"/mnt/{drive.lower()}/{relative.replace(chr(92), '/')}"
    command = ["wsl.exe"]
    if distro:
        command += ["-d", distro]
    command += ["wslpath", "-a", str(path.resolve())]
    result = subprocess.run(command, text=True, capture_output=True, timeout=30)
    if result.returncode:
        raise RuntimeError(result.stderr.strip() or f"wslpath failed for {path}")
    return result.stdout.strip()


def build_alignment_cache(source_fasta: str | Path, target_fasta: str | Path, cache_root: str | Path,
                          minimap2_path: str = "minimap2", distro: str | None = None,
                          threads: int | None = None) -> tuple[Path, bool]:
    """Build target-vs-source whole-genome PAF once. Returns (path, cache_hit)."""
    source, target = Path(source_fasta), Path(target_fasta)
    output = cache_path(source, target, cache_root)
    metadata = output.with_suffix(output.suffix + ".json")
    if output.is_file() and output.stat().st_size > 0 and metadata.is_file():
        return output, True
    output.parent.mkdir(parents=True, exist_ok=True)
    tmp = output.with_suffix(output.suffix + f".{os.getpid()}.tmp")
    thread_count = threads or max(1, min(16, (os.cpu_count() or 4) // 2))
    query_source = _chunk_query_fasta(source, output.parent, distro=distro) if source.stat().st_size > 500_000_000 else source
    target_source = _chunk_query_fasta(target, output.parent, distro=distro) if target.stat().st_size > 500_000_000 else target
    if shutil.which("wsl.exe"):
        src, dst, out = (_wsl_path(p, distro) for p in (query_source, target_source, tmp))
        command = ["wsl.exe"]
        if distro:
            command += ["-d", distro]
        shell_command = "set -o pipefail; " + " ".join(shlex.quote(x) for x in
            [minimap2_path, "-x", "asm20", "-I", "1G", "-K", "50M", "--secondary=no", "-c", "-t", str(thread_count), dst, src]) \
            + " | gzip -1 > " + shlex.quote(out)
        command += ["bash", "-lc", shell_command]
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        # subprocess.run(stdout=handle) writes to handle.fileno() — the raw file under the gzip
        # wrapper — so the ".paf.gz" ends up uncompressed and lift_interval later fails to read
        # it. Stream minimap2's stdout through the compressor instead, with stderr to a file so a
        # full stderr pipe cannot deadlock the single-reader copy.
        stderr_file = tmp.with_suffix(tmp.suffix + ".stderr")
        with gzip.open(tmp, "wb", compresslevel=1) as handle, stderr_file.open("wb") as errors:
            proc = subprocess.Popen([minimap2_path, "-x", "asm20", "-I", "1G", "-K", "50M", "--secondary=no", "-c", "-t", str(thread_count),
                                     str(target), str(source)], stdout=subprocess.PIPE, stderr=errors)
            shutil.copyfileobj(proc.stdout, handle)
            proc.stdout.close()
            proc.wait()
        result = subprocess.CompletedProcess(proc.args, proc.returncode, stderr=stderr_file.read_bytes())
        stderr_file.unlink(missing_ok=True)
    if result.returncode:
        tmp.unlink(missing_ok=True)
        stderr = result.stderr.decode("utf-8", errors="replace") if isinstance(result.stderr, bytes) else result.stderr
        raise RuntimeError((stderr or "minimap2 failed").strip())
    tmp.replace(output)
    metadata.write_text(json.dumps({"source": _fingerprint(source), "target": _fingerprint(target),
                                    "preset": "asm20", "index_batch": "1G", "query_batch": "50M", "threads": thread_count}, indent=2), encoding="utf-8")
    return output, False


def _project(record: list[str], position: int) -> int:
    qstart, qend = int(record[2]), int(record[3])
    tstart, tend = int(record[7]), int(record[8])
    fraction = min(1.0, max(0.0, (position - 1 - qstart) / max(1, qend - qstart)))
    if record[4] == "+":
        return round(tstart + fraction * (tend - tstart)) + 1
    return round(tend - fraction * (tend - tstart))


def _covered_length(rows: list[list[str]], start: int, end: int) -> int:
    """Length of [start, end] covered by the union of the rows' 1-based query spans."""
    spans = sorted((max(start, int(r[2]) + 1), min(end, int(r[3]))) for r in rows)
    total, reach = 0, start - 1
    for s, e in spans:
        s = max(s, reach + 1)
        if e >= s:
            total += e - s + 1
            reach = e
    return total


def _collinear(r1: list[str], r2: list[str], strand: str, max_gap: int) -> bool:
    """True when two same-contig, same-strand records (r1 before r2 in query order) are a
    single collinear alignment merely split apart, e.g. by the query FASTA chunk boundary."""
    src_gap = int(r2[2]) - int(r1[3])
    tgt_gap = int(r2[7]) - int(r1[8]) if strand == "+" else int(r1[7]) - int(r2[8])
    return -max_gap <= src_gap <= max_gap and -max_gap <= tgt_gap <= max_gap


def _best_run(rows: list[list[str]], strand: str, start: int, end: int, max_gap: int) -> list[list[str]]:
    """Split rows (one target contig + strand) into maximal collinear runs and return the run
    covering the most of the requested interval, so genuine breaks are not merged over."""
    rows = sorted(rows, key=lambda r: int(r[2]))
    runs: list[list[list[str]]] = []
    for r in rows:
        if runs and _collinear(runs[-1][-1], r, strand, max_gap):
            runs[-1].append(r)
        else:
            runs.append([r])
    return max(runs, key=lambda run: _covered_length(run, start, end))


def lift_interval(paf_gz: str | Path, source_contig: str, start: int, end: int,
                  target_contig: str | None = None, max_gap: int = 100_000) -> tuple[Interval | None, list[str]]:
    """Project a 1-based source interval onto the target assembly.

    Records are de-chunked, grouped by (target contig, strand), and the group covering most of
    the interval is chosen. Within it, adjacent collinear alignments — including the two halves
    an interval is split into when it crosses a query FASTA chunk boundary — are merged so the
    interval projects as one continuous span instead of being truncated to a single chunk.
    ``max_gap`` bounds how far apart (in source and target bases) records may sit and still be
    treated as collinear; genuine cross-contig, opposite-strand, or distant records are reported
    as ambiguous/partial rather than silently merged.
    """
    candidates: list[list[str]] = []
    with gzip.open(paf_gz, "rt", encoding="utf-8") as handle:
        for line in handle:
            row = line.rstrip().split("\t")
            if len(row) < 12:
                continue
            query_name, query_offset = row[0], 0
            if query_name.startswith("QTLIFT|"):
                _, query_name, raw_offset = query_name.split("|", 2)
                query_offset = int(raw_offset)
                row[0] = query_name
                row[2] = str(int(row[2]) + query_offset)
                row[3] = str(int(row[3]) + query_offset)
            target_name = row[5]
            if target_name.startswith("QTLIFT|"):
                _, target_name, raw_offset = target_name.split("|", 2)
                target_offset = int(raw_offset)
                row[5] = target_name
                row[7] = str(int(row[7]) + target_offset)
                row[8] = str(int(row[8]) + target_offset)
            if query_name != source_contig or (target_contig and target_name != target_contig):
                continue
            qstart, qend = int(row[2]) + 1, int(row[3])
            if max(0, min(end, qend) - max(start, qstart) + 1):
                candidates.append(row)
    if not candidates:
        return None, ["Liftover failed: the source interval has no cached whole-genome alignment."]
    interval_len = end - start + 1
    groups: dict[tuple[str, str], list[list[str]]] = {}
    for row in candidates:
        groups.setdefault((row[5], row[4]), []).append(row)
    ranked = sorted(groups.items(), key=lambda kv: _covered_length(kv[1], start, end), reverse=True)
    (target_name, strand), rows = ranked[0]
    run = _best_run(rows, strand, start, end, max_gap)
    covered_start = max(start, min(int(r[2]) for r in run) + 1)
    covered_end = min(end, max(int(r[3]) for r in run))

    def record_for(pos: int) -> list[str]:
        covering = [r for r in run if int(r[2]) + 1 <= pos <= int(r[3])]
        return covering[0] if covering else min(run, key=lambda r: min(abs(int(r[2]) + 1 - pos), abs(int(r[3]) - pos)))

    a, b = _project(record_for(covered_start), covered_start), _project(record_for(covered_end), covered_end)
    warnings: list[str] = []
    coverage = (covered_end - covered_start + 1) / interval_len
    if len(run) > 1:
        warnings.append(f"Liftover merged {len(run)} adjacent alignment blocks (e.g. across a chunk boundary) into one collinear interval.")
    if sum(1 for _, g in ranked[1:] if _covered_length(g, start, end) >= 0.2 * interval_len):
        warnings.append("Liftover evidence is split across target contigs or orientations; the largest collinear block was reported and the mapping is ambiguous.")
    elif _covered_length(rows, start, end) > _covered_length(run, start, end) + 0.05 * interval_len:
        warnings.append("Liftover alignments on the target contig are not collinear across the interval; only the largest collinear block was projected (partial).")
    if coverage < 0.8:
        warnings.append(f"Liftover alignment covers only {coverage:.1%} of the source interval.")
    mapq = min(int(r[11]) for r in run)
    if mapq < 20:
        warnings.append(f"Liftover alignment has low mapping quality ({mapq}).")
    return Interval(target_name, min(a, b), max(a, b), strand, "liftover"), warnings
