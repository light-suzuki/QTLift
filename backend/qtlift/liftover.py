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
        with gzip.open(tmp, "wb", compresslevel=1) as handle:
            result = subprocess.run([minimap2_path, "-x", "asm20", "-I", "1G", "-K", "50M", "--secondary=no", "-c", "-t", str(thread_count),
                                     str(target), str(source)], stdout=handle, stderr=subprocess.PIPE)
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


def lift_interval(paf_gz: str | Path, source_contig: str, start: int, end: int,
                  target_contig: str | None = None) -> tuple[Interval | None, list[str]]:
    """Project a 1-based source interval through the best overlapping primary PAF alignment."""
    candidates: list[tuple[int, int, list[str]]] = []
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
            overlap = max(0, min(end, qend) - max(start, qstart) + 1)
            if overlap:
                mapq = int(row[11])
                candidates.append((overlap, mapq, row))
    if not candidates:
        return None, ["Liftover failed: the source interval has no cached whole-genome alignment."]
    _, mapq, best = max(candidates, key=lambda item: (item[0], item[1], int(item[2][9])))
    qstart, qend = int(best[2]) + 1, int(best[3])
    covered_start, covered_end = max(start, qstart), min(end, qend)
    a, b = _project(best, covered_start), _project(best, covered_end)
    warnings: list[str] = []
    coverage = (covered_end - covered_start + 1) / (end - start + 1)
    if coverage < 0.8:
        warnings.append(f"Liftover alignment covers only {coverage:.1%} of the source interval.")
    if mapq < 20:
        warnings.append(f"Liftover alignment has low mapping quality ({mapq}).")
    return Interval(best[5], min(a, b), max(a, b), best[4], "liftover"), warnings
