from __future__ import annotations

import re
from dataclasses import dataclass, asdict


@dataclass(slots=True)
class Marker:
    role: str
    name: str
    sequence: str | None = None
    contig: str | None = None
    start: int | None = None
    end: int | None = None

    def as_dict(self): return asdict(self)


DNA = re.compile(r"^[ACGTRYSWKMBDHVN]+$", re.I)


def parse_marker(role: str, text: str) -> tuple[Marker | None, list[str]]:
    # Markers are DNA sequences only: FASTA (>header + bases) or raw bases. Coordinate and
    # name-only markers are no longer accepted.
    text = text.strip()
    if not text:
        return None, []
    if text.startswith(">"):
        lines = [x.strip() for x in text.splitlines() if x.strip()]
        header = lines[0][1:].split()
        name, seq = (header[0] if header else role), "".join(lines[1:]).replace(" ", "").upper()
        if not seq or not DNA.fullmatch(seq):
            raise ValueError(f"Invalid FASTA marker: {role}")
        return Marker(role, name, seq), short_warning(role, seq)
    compact = re.sub(r"\s+", "", text).upper()
    if DNA.fullmatch(compact) and len(compact) >= 4:
        return Marker(role, role, compact), short_warning(role, compact)
    return None, [f"Marker '{role}' must be a DNA sequence (FASTA or raw bases); it was ignored."]


def short_warning(role: str, sequence: str) -> list[str]:
    return [f"Marker {role} is very short ({len(sequence)} bp); hits may be ambiguous."] if len(sequence) < 20 else []

