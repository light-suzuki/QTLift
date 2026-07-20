import gzip
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

from qtlift.liftover import lift_interval


class LiftoverTests(unittest.TestCase):
    def test_liftover_from_cached_paf(self):
        with tempfile.TemporaryDirectory() as directory:
            paf = Path(directory) / "pair.paf.gz"
            with gzip.open(paf, "wt") as handle:
                handle.write(
                    "Chr1\t1000\t99\t900\t+\tTarget1\t2000\t499\t1300\t780\t801\t60"
                    "\ttp:A:P\tcg:Z:801M\n"
                )
            interval, warnings = lift_interval(paf, "Chr1", 200, 800, "Target1")
            self.assertEqual(
                (interval.contig, interval.start, interval.end, interval.strand),
                ("Target1", 600, 1200, "+"),
            )
            self.assertFalse(warnings)

    def test_liftover_merges_across_chunk_boundary(self):
        with tempfile.TemporaryDirectory() as directory:
            paf = Path(directory) / "pair.paf.gz"
            with gzip.open(paf, "wt") as handle:
                handle.write(
                    "QTLIFT|Chr1|0\t20000000\t0\t20000000\t+\tTarget1\t40000000"
                    "\t0\t20000000\t20000000\t20000000\t60\n"
                )
                handle.write(
                    "QTLIFT|Chr1|20000000\t20000000\t0\t20000000\t+\tTarget1\t40000000"
                    "\t20000000\t40000000\t20000000\t20000000\t60\n"
                )
            interval, _warnings = lift_interval(
                paf, "Chr1", 19_500_000, 20_500_000, "Target1"
            )
            self.assertEqual(
                (interval.contig, interval.start, interval.end, interval.strand),
                ("Target1", 19_500_000, 20_500_000, "+"),
            )
            self.assertEqual(interval.end - interval.start + 1, 1_000_001)

    def test_liftover_opposite_strand_chunks_stay_ambiguous(self):
        with tempfile.TemporaryDirectory() as directory:
            paf = Path(directory) / "pair.paf.gz"
            with gzip.open(paf, "wt") as handle:
                handle.write(
                    "QTLIFT|Chr1|0\t20000000\t0\t20000000\t+\tTarget1\t40000000"
                    "\t0\t20000000\t20000000\t20000000\t60\n"
                )
                handle.write(
                    "QTLIFT|Chr1|20000000\t20000000\t0\t20000000\t-\tTarget1\t40000000"
                    "\t20000000\t40000000\t20000000\t20000000\t60\n"
                )
            interval, warnings = lift_interval(
                paf, "Chr1", 19_500_000, 20_500_000, "Target1"
            )
            self.assertIn("ambiguous", " ".join(warnings).lower())
            self.assertLess(interval.end - interval.start + 1, 1_000_001)

    def test_liftover_distant_same_strand_blocks_not_merged(self):
        with tempfile.TemporaryDirectory() as directory:
            paf = Path(directory) / "pair.paf.gz"
            with gzip.open(paf, "wt") as handle:
                handle.write(
                    "QTLIFT|Chr1|0\t20000000\t0\t20000000\t+\tTarget1\t200000000"
                    "\t0\t20000000\t20000000\t20000000\t60\n"
                )
                handle.write(
                    "QTLIFT|Chr1|20000000\t20000000\t0\t20000000\t+\tTarget1\t200000000"
                    "\t100000000\t120000000\t20000000\t20000000\t60\n"
                )
            interval, warnings = lift_interval(
                paf, "Chr1", 19_500_000, 20_500_000, "Target1"
            )
            self.assertLess(interval.end - interval.start + 1, 1_000_001)
            self.assertIn("partial", " ".join(warnings).lower())


if __name__ == "__main__":
    unittest.main()
