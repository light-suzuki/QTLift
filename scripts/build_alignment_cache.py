from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

from qtlift.liftover import build_alignment_cache


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a reusable QTLift whole-genome minimap2 PAF cache.")
    parser.add_argument("source_fasta")
    parser.add_argument("target_fasta")
    parser.add_argument("--cache-root", default=str(ROOT / "output" / "alignment_cache"))
    parser.add_argument("--minimap2", default="minimap2")
    parser.add_argument("--wsl-distro")
    parser.add_argument("--threads", type=int)
    args = parser.parse_args()
    path, cache_hit = build_alignment_cache(args.source_fasta, args.target_fasta, args.cache_root,
                                            args.minimap2, args.wsl_distro, args.threads)
    print(f"cache={'reused' if cache_hit else 'created'}\t{path}")


if __name__ == "__main__":
    main()
