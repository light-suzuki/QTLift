from __future__ import annotations

from collections import Counter

from Bio.Seq import Seq

from .genomes import sequence_slice
from .models import Gene, Hit, Interval, Params


def select_anchors(genes: list[Gene], start: int, end: int, peak: int | None, params: Params) -> list[Gene]:
    if not genes:
        return []
    genes = sorted(genes, key=lambda g: g.start)
    selected: dict[str, Gene] = {}
    def add(gene, role):
        if gene.id not in selected:
            gene.role = role
            selected[gene.id] = gene
    # The two interval boundaries and the peak gene are always kept so orientation is
    # anchored. Internal anchors are then driven by a minimum genomic separation
    # (anchor_spacing_mb): walk genes left-to-right and keep one only when it lies at least
    # `spacing` from the previously kept anchor. There is no upper cap — a wider interval
    # simply yields proportionally more anchors.
    add(genes[0], "left-edge")
    add(genes[-1], "right-edge")
    center = peak if peak is not None else (start + end) // 2
    add(min(genes, key=lambda g: abs((g.start + g.end)//2 - center)), "peak")
    span = max(1, end - start)

    def greedy(separation: float) -> list[Gene]:
        picks: list[Gene] = []
        last: float | None = None
        for gene in genes:
            mid = (gene.start + gene.end) // 2
            if last is None or mid - last >= separation:
                picks.append(gene)
                last = mid
        return picks

    spacing = params.anchor_spacing_mb * 1_000_000 if params.anchor_spacing_mb and params.anchor_spacing_mb > 0 else span
    picks = greedy(spacing)
    # Guarantee a floor so short or gene-sparse intervals still get a usable result: if the
    # spacing yields fewer than min_anchors, spread that many genes evenly across the interval
    # (bounded by how many genes actually exist).
    target = min(params.min_anchors, len(genes))
    if len(picks) < target:
        step = (len(genes) - 1) / (target - 1) if target > 1 else 0
        picks = [genes[round(i * step)] for i in range(target)]
    for gene in picks:
        add(gene, "internal")
    return sorted(selected.values(), key=lambda g: g.start)


def exact_hits(query_id: str, sequence: str, target_fasta: str, contigs: list[dict], source_start: int | None = None, target_contig: str | None = None) -> list[Hit]:
    hits: list[Hit] = []
    reverse = str(Seq(sequence).reverse_complement())
    for item in contigs:
        if target_contig and item["name"] != target_contig:
            continue
        target = sequence_slice(target_fasta, item["name"], 1, item["length"])
        for strand, needle in (("+", sequence), ("-", reverse)):
            pos = target.find(needle)
            while pos >= 0:
                hits.append(Hit(query_id, item["name"], pos+1, pos+len(needle), strand, 100.0, 100.0, method="exact-fallback", source_start=source_start))
                pos = target.find(needle, pos+1)
    for hit in hits: hit.hit_count = len(hits)
    return hits


def evaluate_synteny(hits: list[Hit]) -> tuple[str, Interval | None, list[str]]:
    usable = [h for h in hits if h.hit_count == 1]
    if len(usable) < 2:
        return "failed", None, ["Too few unique anchor hits to infer synteny."]
    contigs = Counter(h.contig for h in usable)
    major, count = contigs.most_common(1)[0]
    if count < len(usable) * 0.6:
        return "split", None, ["Anchor hits are split across target contigs."]
    major_hits = sorted((h for h in usable if h.contig == major), key=lambda h: h.source_start or 0)
    coords = [(h.start+h.end)//2 for h in major_hits]
    inc = sum(b > a for a,b in zip(coords, coords[1:])); dec = sum(b < a for a,b in zip(coords, coords[1:]))
    total = max(1, len(coords)-1)
    if inc/total >= .8: state, strand = "forward", "+"
    elif dec/total >= .8: state, strand = "reverse", "-"
    else: state, strand = "partial", "."
    warnings = [] if state in ("forward", "reverse") else ["Anchor order is only partially collinear."]
    return state, Interval(major, min(h.start for h in major_hits), max(h.end for h in major_hits), strand, "synteny"), warnings


def marker_interval(hits: list[Hit]) -> tuple[Interval | None, list[str]]:
    unique = [h for h in hits if h.hit_count == 1]
    if len(unique) < 2:
        return None, ["Fewer than two uniquely mapped markers; marker interval unavailable."]
    major, count = Counter(h.contig for h in unique).most_common(1)[0]
    if count < 2:
        return None, ["Unique marker hits occur on different target contigs."]
    rows = [h for h in unique if h.contig == major]
    return Interval(major, min(h.start for h in rows), max(h.end for h in rows), ".", "markers"), []


def reconcile_orientation(intervals: list[Interval]) -> dict:
    """Reconcile the strand of a set of evidence intervals. An informative +/- strand always
    beats an uninformative '.', informative strands that disagree are flagged as a conflict, and
    the result is '.' only when no evidence establishes orientation."""
    informative = [x for x in intervals if x.strand in ("+", "-")]
    strands = {x.strand for x in informative}
    conflict = len(strands) > 1
    return {
        "strand": "." if conflict or not strands else next(iter(strands)),
        "orientation_evidence": sorted({x.evidence for x in informative if x.evidence}),
        "uninformative_orientation_evidence": sorted({x.evidence for x in intervals if x.strand not in ("+", "-") and x.evidence}),
        "conflict": conflict,
    }


def orientation_audit(synteny: Interval | None, marker: Interval | None, liftover: Interval | None) -> dict:
    """Auditable record of how the reported orientation was determined across evidence classes."""
    return reconcile_orientation([x for x in (synteny, marker, liftover) if x])


def score_confidence(synteny_state: str, synteny: Interval | None, marker: Interval | None, liftover: Interval | None, anchor_hits: list[Hit]) -> tuple[str, list[str], list[str]]:
    reasons, warnings = [], []
    evidence = [x for x in (synteny, marker, liftover) if x]
    oriented = [x for x in evidence if x.strand in ("+", "-")]
    unique = sum(h.hit_count == 1 for h in anchor_hits)
    if synteny_state in ("split", "failed"):
        return "Manual check", [f"Synteny state is {synteny_state}."], ["Mapping requires manual review."]
    agree = False
    if len(evidence) >= 2:
        contigs = {x.contig for x in evidence}
        agree = len(contigs) == 1 and max(x.start for x in evidence) <= min(x.end for x in evidence)
        if not agree:
            warnings.append("Marker/synteny/liftover intervals are inconsistent.")
            return "Manual check", ["Independent evidence intervals disagree and must be reviewed separately."], warnings
        if len({x.strand for x in oriented}) > 1:
            detail = ", ".join(f"{x.evidence or 'evidence'} {x.strand}" for x in oriented)
            warnings.append(f"Conflicting orientation evidence: {detail}; intervals retained separately.")
            return "Manual check", ["Independent evidence overlaps but disagrees on orientation and cannot be combined."], warnings
    if evidence and not oriented:
        warnings.append("Orientation is unresolved: no strand-informative evidence; the interval strand is reported as '.'.")
    if synteny_state in ("forward", "reverse") and unique >= 4 and len(evidence) >= 2 and agree:
        confidence = "High"; reasons.append("At least two evidence classes agree with four or more unique collinear anchors.")
    elif synteny and unique >= 3 and synteny_state in ("forward", "reverse", "partial"):
        confidence = "Medium"; reasons.append("A coherent anchor interval is supported, but independent evidence is limited.")
    elif evidence:
        confidence = "Low"; reasons.append("Only sparse or partially consistent evidence is available.")
    else:
        confidence = "Manual check"; reasons.append("No defensible target interval was produced.")
    return confidence, reasons, warnings


def combine_intervals(synteny: Interval | None, marker: Interval | None, liftover: Interval | None) -> list[Interval]:
    available = [x for x in (liftover, marker, synteny) if x]
    grouped: dict[str, list[Interval]] = {}
    for item in available: grouped.setdefault(item.contig, []).append(item)
    result: list[Interval] = []
    for contig, rows in grouped.items():
        if len(rows) == 1:
            result.append(rows[0]); continue
        overlap_start, overlap_end = max(x.start for x in rows), min(x.end for x in rows)
        orient = reconcile_orientation(rows)
        # Merge only when the intervals overlap AND their informative orientations agree; an
        # uninformative '.' never overwrites a +/- strand, and a real orientation conflict keeps
        # the evidence intervals separate instead of flattening it into one.
        if overlap_start <= overlap_end and not orient["conflict"]:
            result.append(Interval(contig, overlap_start, overlap_end, orient["strand"], "+".join(x.evidence for x in rows)))
        else:
            result.extend(rows)
    # The pipeline takes candidates[0] as the final interval, so order deterministically —
    # most evidence classes first, then widest span, then contig/start — rather than leaving it
    # to dict-insertion order.
    result.sort(key=lambda iv: (-len(iv.evidence.split("+")) if iv.evidence else 0, -(iv.end - iv.start), iv.contig, iv.start))
    return result
