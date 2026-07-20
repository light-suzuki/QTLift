import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

from qtlift.analysis import combine_intervals
from qtlift.models import Interval


class CombineOrderingTests(unittest.TestCase):
    """combine_intervals returns a deterministic order so pipeline's candidates[0] is stable (#33)."""

    def test_most_supported_interval_ranks_first(self):
        syn = Interval("chr1", 100, 400, "+", "synteny")
        marker = Interval("chr1", 150, 350, "+", "markers")
        lift = Interval("chr2", 100, 200, "+", "liftover")
        result = combine_intervals(syn, marker, lift)
        # The chr1 synteny+marker overlap merges into a two-class interval, which outranks the
        # single-class chr2 liftover interval.
        self.assertEqual(result[0].contig, "chr1")
        self.assertIn("+", result[0].evidence)

    def test_order_is_input_position_independent(self):
        a = Interval("chr3", 100, 200, "+", "synteny")
        b = Interval("chr1", 100, 500, "+", "liftover")
        r1 = [(iv.contig, iv.start, iv.end) for iv in combine_intervals(a, None, b)]
        r2 = [(iv.contig, iv.start, iv.end) for iv in combine_intervals(None, a, b)]
        self.assertEqual(r1, r2)
        # Two single-class intervals: the wider one (chr1, 401 bp) ranks before chr3 (101 bp).
        self.assertEqual(r1[0][0], "chr1")


if __name__ == "__main__":
    unittest.main()
