from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Literal

Confidence = Literal["High", "Medium", "Low", "Manual check"]


@dataclass(slots=True)
class Gene:
    id: str
    contig: str
    start: int
    end: int
    strand: str = "+"
    cds: list[tuple[int, int]] = field(default_factory=list)
    role: str = "internal"
    # Auditability for anchor construction: the transcript that supplied the CDS model and
    # whether the anchor sequence comes from that CDS ("cds") or a whole-gene fallback ("gene").
    transcript_id: str | None = None
    sequence_source: str = "gene"


@dataclass(slots=True)
class Hit:
    query_id: str
    contig: str
    start: int
    end: int
    strand: str
    identity: float
    coverage: float
    hit_count: int = 1
    method: str = "exact-fallback"
    source_start: int | None = None
    query_start: int | None = None
    query_end: int | None = None
    query_length: int | None = None


@dataclass(slots=True)
class Interval:
    contig: str
    start: int
    end: int
    strand: str = "+"
    evidence: str = ""

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(slots=True)
class Params:
    # Anchor selection is spacing-driven: `anchor_spacing_mb` is the minimum genomic
    # separation between BLASTed anchors (the main knob, no upper cap on count), and
    # `min_anchors` is the floor that guarantees a usable result on short/gene-sparse
    # intervals by shrinking the separation until at least this many genes are chosen.
    anchor_spacing_mb: float = 0.25
    min_anchors: int = 6
    marker_flank_length: int = 250
    min_identity: float = 90.0
    min_coverage: float = 70.0
    multi_hit_threshold: int = 5
    include_outside_edge_genes: bool = True

    @classmethod
    def preset(cls, name: str) -> "Params":
        values = {
            "Fast": cls(anchor_spacing_mb=0.5, min_anchors=4, marker_flank_length=150, min_coverage=65, include_outside_edge_genes=False),
            "Standard": cls(),
            "Precise": cls(anchor_spacing_mb=0.1, min_anchors=10, marker_flank_length=500, min_identity=92, min_coverage=80, multi_hit_threshold=3),
            "Auto": cls(anchor_spacing_mb=0.5, min_anchors=6, marker_flank_length=250),
        }
        return values.get(name, cls())
