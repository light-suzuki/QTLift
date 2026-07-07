from __future__ import annotations

import gzip
import json
import os
import re
from functools import lru_cache
from pathlib import Path
from urllib.parse import unquote

from Bio import SeqIO

from .models import Gene

FASTA_SUFFIXES = (".fa", ".fasta", ".fna")
GFF_SUFFIXES = (".gff3", ".gff")
LOCAL_LIBRARY_TOKEN = "@this-pc"
LOCAL_CONFIG = Path(os.environ.get("QTLIFT_LOCAL_GENOMES", Path(__file__).resolve().parents[2] / "data" / "local_genomes.json"))


def _base_suffix(path: Path) -> str:
    name = path.name.lower()
    if name.endswith(".gz"):
        name = name[:-3]
    return next((s for s in FASTA_SUFFIXES + GFF_SUFFIXES if name.endswith(s)), "")


def detect_genomes(root: str | Path) -> list[dict]:
    if str(root) == LOCAL_LIBRARY_TOKEN:
        return detect_configured_genomes()
    root = Path(root).expanduser().resolve()
    if not root.is_dir():
        raise ValueError(f"Genome root is not a directory: {root}")
    result = []
    for directory in sorted(p for p in root.iterdir() if p.is_dir()):
        files = [p for p in directory.iterdir() if p.is_file()]
        fasta = next((p for p in files if _base_suffix(p) in FASTA_SUFFIXES), None)
        gff = next((p for p in files if _base_suffix(p) in GFF_SUFFIXES), None)
        contigs = fasta_lengths(fasta) if fasta else []
        result.append({"name": directory.name, "path": str(directory), "fasta": str(fasta) if fasta else None,
                       "gff": str(gff) if gff else None, "fasta_status": "ready" if fasta else "missing",
                       "gff_status": "ready" if gff else "missing", "contigs": contigs,
                       "gene_count": count_genes(gff) if gff else 0})
    return result


@lru_cache(maxsize=1)
def detect_configured_genomes() -> list[dict]:
    if not LOCAL_CONFIG.is_file():
        raise ValueError(f"Local genome configuration is missing: {LOCAL_CONFIG}")
    profiles = json.loads(LOCAL_CONFIG.read_text(encoding="utf-8"))
    result = []
    for profile in profiles:
        fasta = Path(profile["fasta"])
        gff = Path(profile["gff"])
        result.append({
            "name": profile["name"], "path": str(fasta.parent),
            "fasta": str(fasta), "gff": str(gff),
            "fasta_status": "ready" if fasta.is_file() else "missing",
            "gff_status": "ready" if gff.is_file() else "missing",
            "contigs": fasta_lengths(fasta) if fasta.is_file() else [],
            "gene_count": count_genes(gff) if gff.is_file() else 0,
            "source": "this-pc",
            "gff_contig_aliases": profile.get("gff_contig_aliases", {}),
            "blast_db": profile.get("blast_db"),
        })
    return result


def _open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.name.lower().endswith(".gz") else path.open("r", encoding="utf-8")


def fasta_lengths(path: str | Path) -> list[dict]:
    path = Path(path)
    fai = Path(str(path) + ".fai")
    if fai.is_file():
        rows = []
        with fai.open("r", encoding="utf-8") as handle:
            for line in handle:
                cols = line.rstrip().split("\t")
                if len(cols) >= 2:
                    rows.append({"name": cols[0], "length": int(cols[1])})
        return rows
    with _open_text(path) as handle:
        return [{"name": record.id, "length": len(record.seq)} for record in SeqIO.parse(handle, "fasta")]


@lru_cache(maxsize=32)
def count_genes(path: str | Path) -> int:
    with _open_text(Path(path)) as handle:
        return sum(1 for line in handle if not line.startswith("#") and "\tgene\t" in line)


def parse_attrs(text: str) -> dict[str, str]:
    return {k: unquote(v.strip('"')) for part in text.strip().split(";") if part and "=" in part for k, v in [part.split("=", 1)]}


def genes_in_interval(path: str | Path, contig: str, start: int, end: int, include_outside: bool = False) -> list[Gene]:
    validate_interval(start, end)
    genes: dict[str, Gene] = {}
    cds_rows: list[tuple[str, int, int]] = []
    with _open_text(Path(path)) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip().split("\t")
            if len(cols) != 9 or cols[0] != contig:
                continue
            feature, fstart, fend, strand, attrs = cols[2], int(cols[3]), int(cols[4]), cols[6], parse_attrs(cols[8])
            overlaps = fend >= start and fstart <= end
            if feature == "gene" and overlaps:
                gid = attrs.get("ID") or attrs.get("Name") or f"gene_{fstart}_{fend}"
                genes[gid] = Gene(gid, contig, fstart, fend, strand)
            elif feature == "CDS" and overlaps:
                parent = (attrs.get("Parent") or "").split(",")[0]
                cds_rows.append((parent, attrs.get("Name", ""), fstart, fend))
    for parent, name, fstart, fend in cds_rows:
        # CDS Parent is usually the transcript id (gene id + a transcript suffix such as
        # ".1", ".t1", "-T1"); resolve it back to the gene by trying the raw parent, the
        # parent with a trailing transcript suffix stripped, then the CDS Name attribute.
        stripped = re.sub(r"(?:[.]t\d+|-T\d+|[.]\d+|[.-]mRNA[._\d]*)$", "", parent, flags=re.IGNORECASE)
        key = next((c for c in (parent, stripped, name) if c and c in genes), None)
        if key:
            genes[key].cds.append((fstart, fend))
    return sorted(genes.values(), key=lambda g: (g.start, g.end))


def validate_interval(start: int, end: int, length: int | None = None, peak: int | None = None) -> None:
    if start < 1 or end < start:
        raise ValueError("Coordinates must be 1-based and start <= end")
    if length and end > length:
        raise ValueError(f"Interval end {end} exceeds contig length {length}")
    if peak is not None and not start <= peak <= end:
        raise ValueError("Peak must fall inside the source interval")


def sequence_slice(path: str | Path, contig: str, start: int, end: int) -> str:
    path = Path(path)
    if path.name.lower().endswith(".gz"):
        with _open_text(path) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if record.id == contig:
                    return str(record.seq[start - 1:end]).upper()
    elif Path(str(path) + ".fai").is_file():
        fai = Path(str(path) + ".fai")
        row = None
        with fai.open("r", encoding="utf-8") as handle:
            for line in handle:
                cols = line.rstrip().split("\t")
                if cols[0] == contig:
                    row = cols
                    break
        if row is None:
            raise KeyError(contig)
        length, offset, line_bases, line_width = map(int, row[1:5])
        start = max(1, start); end = min(length, end)
        start0, end0 = start - 1, end
        byte_start = offset + (start0 // line_bases) * line_width + (start0 % line_bases)
        newline_bytes = line_width - line_bases
        byte_count = (end0 - start0) + ((end0 // line_bases) - (start0 // line_bases) + 2) * newline_bytes
        with path.open("rb") as handle:
            handle.seek(byte_start)
            raw = handle.read(byte_count)
        return raw.replace(b"\n", b"").replace(b"\r", b"")[:end0-start0].decode("ascii").upper()
    else:
        index = SeqIO.index(str(path), "fasta")
        try:
            if contig not in index:
                raise KeyError(contig)
            return str(index[contig].seq[start - 1:end]).upper()
        finally:
            index.close()
    raise KeyError(contig)
