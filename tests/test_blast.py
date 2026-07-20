import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

from qtlift.blast import _collapse_loci, _find_blast_entry
from qtlift.models import Hit


class BlastTests(unittest.TestCase):
    def test_blast_contig_entry_resolution(self):
        metadata = "gnl|BL_ORD_ID|0\tchr1\ngnl|BL_ORD_ID|1\tchr2 description\n"
        self.assertEqual(_find_blast_entry(metadata, "chr2"), "gnl|BL_ORD_ID|1")
        self.assertIsNone(_find_blast_entry(metadata, "chr20"))

    def test_hsp_coverage_uses_query_union_and_separates_paralogs(self):
        hsps = [
            Hit(
                "q", "chr1", 100, 149, "+", 99, 50, method="blastn",
                query_start=1, query_end=50, query_length=100,
            ),
            Hit(
                "q", "chr1", 200, 249, "+", 97, 50, method="blastn",
                query_start=51, query_end=100, query_length=100,
            ),
        ]
        loci = _collapse_loci(hsps)
        self.assertEqual(len(loci), 1)
        self.assertEqual(loci[0].coverage, 100)

        paralog = Hit(
            "q", "chr1", 50_000, 50_099, "+", 96, 100, method="blastn",
            query_start=1, query_end=100, query_length=100,
        )
        self.assertEqual(len(_collapse_loci(hsps + [paralog])), 2)
        self.assertTrue(all(hit.coverage < 70 for hit in hsps))
        self.assertGreaterEqual(loci[0].coverage, 70)


if __name__ == "__main__":
    unittest.main()
