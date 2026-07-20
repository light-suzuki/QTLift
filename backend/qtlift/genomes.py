from __future__ import annotations

import gzip
import json
import os
from functools import lru_cache
from pathlib import Path
from urllib.parse import unquote

from Bio import SeqIO
from Bio.Seq import Seq

from .models import Gene

# Feature types that are never a coding transcript and never carry CDS children; skipping
# them keeps the per-contig feature graph small on large annotations.
LEAF_FEATURE_TYPES = frozenset({"exon", "five_prime_utr", "three_prime_utr", "utr", "intron"})
# Attribute conventions used across GFF3 dialects to flag the representative transcript.
CANONICAL_TAGS = frozenset({"ensembl_canonical", "canonical", "mane_select", "mane_plus_clinical"})
CANONICAL_KEYS = ("canonical", "is_canonical", "representative", "canonical_transcript")

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


def _parent_ids(col9: str) -> list[str]:
    # Per the GFF3 spec multiple parents are comma-separated and literal commas inside a
    # value are percent-encoded, so split on the raw commas before unquoting each id.
    for part in col9.split(";"):
        if part.strip().startswith("Parent="):
            raw = part.split("=", 1)[1].strip().strip('"')
            return [unquote(v) for v in raw.split(",") if v]
    return []


def _load_contig_features(path: str | Path, contig: str):
    """Read one contig into a feature graph: genes, CDS grouped by their parent id, and the
    parents/attributes of every id-bearing intermediate feature (candidate transcripts)."""
    genes: dict[str, tuple[int, int, str]] = {}
    cds_by_parent: dict[str, list[tuple[int, int]]] = {}
    tx_parents: dict[str, list[str]] = {}
    tx_attrs: dict[str, dict[str, str]] = {}
    with _open_text(Path(path)) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9 or cols[0] != contig:
                continue
            ftype, fstart, fend, strand = cols[2], int(cols[3]), int(cols[4]), cols[6]
            attrs = parse_attrs(cols[8])
            fid = attrs.get("ID")
            if ftype == "gene":
                gid = fid or attrs.get("Name") or f"gene_{fstart}_{fend}"
                genes[gid] = (fstart, fend, strand)
            elif ftype == "CDS":
                for parent in _parent_ids(cols[8]) or ([fid] if fid else []):
                    cds_by_parent.setdefault(parent, []).append((fstart, fend))
            elif fid and ftype.lower() not in LEAF_FEATURE_TYPES:
                tx_parents[fid] = _parent_ids(cols[8])
                tx_attrs[fid] = attrs
    return genes, cds_by_parent, tx_parents, tx_attrs


def _is_canonical(attrs: dict[str, str]) -> bool:
    tags = attrs.get("tag", "")
    if any(t.strip().lower() in CANONICAL_TAGS for t in tags.split(",")):
        return True
    return any(attrs.get(k, "").strip().lower() in {"1", "true", "yes"} for k in CANONICAL_KEYS)


def _representative_cds(gid: str, cds_by_parent, tx_parents, tx_attrs) -> tuple[str | None, list[tuple[int, int]]]:
    """Pick one transcript for a gene through real gene->transcript->CDS relationships.
    Preference: an explicitly canonical/representative transcript, else the longest valid
    CDS, with transcript-id order as a deterministic tie-break. Never mixes isoforms."""
    candidates: list[tuple[str, list[tuple[int, int]], dict[str, str]]] = []
    for tid, parents in tx_parents.items():
        if gid in parents and tid in cds_by_parent:
            candidates.append((tid, sorted(cds_by_parent[tid]), tx_attrs.get(tid, {})))
    if gid in cds_by_parent:  # CDS parented directly to the gene (no transcript layer)
        candidates.append((gid, sorted(cds_by_parent[gid]), {}))
    if not candidates:
        return None, []
    canonical = [c for c in candidates if _is_canonical(c[2])]
    pool = canonical or candidates
    length = lambda cds: sum(e - s + 1 for s, e in cds)
    tid, cds, _ = min(pool, key=lambda c: (-length(c[1]), c[0]))
    return tid, cds


def genes_in_interval(path: str | Path, contig: str, start: int, end: int, include_outside: bool = False) -> list[Gene]:
    """Select genes for a source interval and resolve one complete CDS model per gene.

    ``include_outside`` follows the anchor contract: when False only genes fully contained in
    the interval are returned; when True genes overlapping either edge are also returned. In
    both cases the complete CDS of the chosen transcript is loaded, so edge genes keep every
    CDS segment instead of being truncated to the part that falls inside the interval.
    """
    validate_interval(start, end)
    genes_meta, cds_by_parent, tx_parents, tx_attrs = _load_contig_features(path, contig)
    result: list[Gene] = []
    for gid, (gstart, gend, strand) in genes_meta.items():
        keep = (gend >= start and gstart <= end) if include_outside else (gstart >= start and gend <= end)
        if not keep:
            continue
        gene = Gene(gid, contig, gstart, gend, strand)
        tid, cds = _representative_cds(gid, cds_by_parent, tx_parents, tx_attrs)
        if cds:
            gene.cds, gene.transcript_id, gene.sequence_source = cds, tid, "cds"
        result.append(gene)
    return sorted(result, key=lambda g: (g.start, g.end))


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


def anchor_sequence(fasta_path: str | Path, gene: Gene) -> str:
    """Assemble a gene's anchor sequence from its complete CDS model (the whole-gene span is
    the documented fallback), reverse-complemented for minus-strand genes so a multi-exon CDS
    is read in biological 5'->3' order rather than plain genomic order."""
    ranges = sorted(gene.cds) if gene.cds else [(gene.start, gene.end)]
    seq = "".join(sequence_slice(fasta_path, gene.contig, a, b) for a, b in ranges)
    if gene.strand == "-":
        seq = str(Seq(seq).reverse_complement())
    return seq
