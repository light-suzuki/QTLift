"""End-to-end integration tests that drive the *real* external tools (blastn, minimap2).

These exercise the code paths that unit tests deliberately stub out — the historically
bug-prone BLAST/liftover subsystems (#15, #16, #17). They run against small self-contained
synthetic sequences so they are deterministic and need no reference data.

Tool selection is automatic: a native binary on PATH (the CI runner installs
``ncbi-blast+`` and ``minimap2``) is used directly; on a Windows dev box the same code
reaches the tools through WSL. When neither is available the test skips, so the default
unit-test run is unaffected.
"""
import random
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

from qtlift.blast import blast_many
from qtlift.liftover import build_alignment_cache, lift_interval

# "windows" here just means "run the local blastn binary directly" (not via wsl.exe); on the
# Linux CI runner that is the native blastn, on a Windows dev box we fall back to WSL blastn.
BLAST_BACKEND = "windows" if shutil.which("blastn") else ("wsl" if shutil.which("wsl.exe") else None)
HAVE_MINIMAP2 = bool(shutil.which("minimap2") or shutil.which("wsl.exe"))


def _random_seq(seed: int, length: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


class IntegrationTests(unittest.TestCase):
    @unittest.skipUnless(BLAST_BACKEND, "requires blastn (native binary or via WSL)")
    def test_blast_backend_maps_a_known_sequence(self):
        motif = _random_seq(7, 300)
        flank = _random_seq(8, 500)
        with tempfile.TemporaryDirectory() as d:
            target = Path(d) / "target.fa"
            target.write_text(f">chrT\n{flank}{motif}{flank}\n", encoding="ascii")
            hits = blast_many([("q", motif, 1)], str(target), BLAST_BACKEND,
                              min_identity=95, min_coverage=90).get("q", [])
        self.assertTrue(hits, "blastn returned no hit for an exact 300 bp query")
        self.assertEqual(hits[0].contig, "chrT")
        self.assertGreaterEqual(hits[0].identity, 95)
        # The motif was embedded after the first 500 bp flank (1-based position 501).
        self.assertLess(abs(hits[0].start - 501), 5)

    @unittest.skipUnless(HAVE_MINIMAP2, "requires minimap2 (native binary or via WSL)")
    def test_minimap2_liftover_projects_interval(self):
        seq = _random_seq(11, 60_000)
        with tempfile.TemporaryDirectory() as d:
            src = Path(d) / "src.fa"
            tgt = Path(d) / "tgt.fa"
            src.write_text(f">chr1\n{seq}\n", encoding="ascii")
            tgt.write_text(f">chr1\n{seq}\n", encoding="ascii")  # identical -> collinear alignment
            cache = Path(d) / "cache"
            cache.mkdir()
            paf, cache_hit = build_alignment_cache(src, tgt, cache)
            self.assertFalse(cache_hit)
            interval, _warnings = lift_interval(paf, "chr1", 20_000, 30_000)
        self.assertIsNotNone(interval, "no liftover interval from identical 60 kb sequences")
        self.assertEqual(interval.contig, "chr1")
        # Identical sequences: the interval must project onto essentially itself.
        self.assertLess(abs(interval.start - 20_000), 300)
        self.assertLess(abs(interval.end - 30_000), 300)


if __name__ == "__main__":
    unittest.main()
