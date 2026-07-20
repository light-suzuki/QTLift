import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

from qtlift.markers import parse_marker


class MarkerParsingTests(unittest.TestCase):
    def test_marker_parsing(self):
        marker, warnings = parse_marker("left", ">m1\nACGTACGTACGTACGTACGT")
        self.assertEqual(marker.name, "m1")
        self.assertFalse(warnings)

        marker, _warnings = parse_marker("peak", "ACGTACGTACGTACGTACGT")
        self.assertEqual(marker.sequence, "ACGTACGTACGTACGTACGT")

        marker, warnings = parse_marker("peak", "mk=Chr1:10-20")
        self.assertIsNone(marker)
        self.assertTrue(warnings)

        marker, warnings = parse_marker("right", "SSR-12")
        self.assertIsNone(marker)
        self.assertTrue(warnings)


if __name__ == "__main__":
    unittest.main()
