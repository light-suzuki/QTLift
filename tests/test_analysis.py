import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

from qtlift.analysis import (
    anchor_contig_distribution,
    combine_intervals,
    evaluate_synteny,
    orientation_audit,
    score_confidence,
    select_anchors,
)
from qtlift.models import Gene, Hit, Interval, Params


class SyntenyMajorContigTests(unittest.TestCase):
    """Unique anchors outside the adopted major contig must be visible and must not
    overstate confidence (#51)."""

    def _forward(self, contig, count, start=0):
        return [Hit(f"{contig}{i}", contig, start + i * 100, start + i * 100 + 10, "+", 100, 100, 1, "x", i * 90) for i in range(1, count + 1)]

    def test_minority_anchors_warn_and_are_tracked(self):
        hits = self._forward("A", 3) + self._forward("B", 2, start=300)
        state, interval, warnings = evaluate_synteny(hits)
        self.assertEqual(state, "forward")
        self.assertEqual(interval.contig, "A")
        self.assertTrue(any("outside the major contig" in w for w in warnings))
        self.assertIn("B", " ".join(warnings))
        distribution = anchor_contig_distribution(hits)
        self.assertEqual(distribution["major_contig"], "A")
        self.assertEqual(distribution["major_count"], 3)
        self.assertEqual(distribution["total"], 5)
        self.assertEqual(distribution["major_fraction"], 0.6)
        self.assertEqual(len(distribution["minority"]), 2)
        self.assertTrue(all(x["contig"] == "B" for x in distribution["minority"]))

    def test_minority_fraction_caps_high(self):
        hits = self._forward("A", 3) + self._forward("B", 2, start=300)
        _state, synteny, _warnings = evaluate_synteny(hits)
        marker = Interval("A", 100, 500)
        confidence, _reasons, _warnings = score_confidence("forward", synteny, marker, None, hits, major_fraction=0.6)
        self.assertNotEqual(confidence, "High")
        confidence = score_confidence("forward", synteny, marker, None, hits, major_fraction=0.8)[0]
        self.assertEqual(confidence, "High")

    def test_all_anchors_same_contig_unchanged(self):
        hits = self._forward("A", 5)
        state, _interval, warnings = evaluate_synteny(hits)
        self.assertEqual(state, "forward")
        self.assertEqual(warnings, [])
        distribution = anchor_contig_distribution(hits)
        self.assertEqual(distribution["major_fraction"], 1.0)
        self.assertEqual(distribution["minority"], [])

    def test_four_of_five_warns_but_stays_high(self):
        hits = self._forward("A", 4) + self._forward("B", 1, start=400)
        _state, synteny, warnings = evaluate_synteny(hits)
        self.assertTrue(any("outside the major contig" in w for w in warnings))
        confidence = score_confidence("forward", synteny, Interval("A", 100, 500), None, hits, major_fraction=0.8)[0]
        self.assertEqual(confidence, "High")


class AnalysisTests(unittest.TestCase):
    def test_anchor_selection(self):
        genes = [Gene(f"g{i}", "c", i * 100, i * 100 + 20) for i in range(1, 11)]
        anchors = select_anchors(genes, 100, 1020, 500, Params())
        self.assertLessEqual(len(anchors), 16)
        self.assertIn("left-edge", {x.role for x in anchors})

    def test_synteny_forward_reverse_split(self):
        forward = [
            Hit(str(i), "A", i * 100, i * 100 + 10, "+", 100, 100, 1, "x", i * 90)
            for i in range(1, 6)
        ]
        state, interval, _warnings = evaluate_synteny(forward)
        self.assertEqual(state, "forward")
        self.assertEqual(interval.contig, "A")

        reverse = [
            Hit(str(i), "A", 1000 - i * 100, 1010 - i * 100, "-", 100, 100, 1, "x", i * 90)
            for i in range(1, 6)
        ]
        self.assertEqual(evaluate_synteny(reverse)[0], "reverse")

    def test_confidence(self):
        interval = Interval("A", 100, 500)
        hits = [Hit(str(i), "A", i, i + 1, "+", 100, 100) for i in range(5)]
        confidence, _reasons, _warnings = score_confidence(
            "forward", interval, Interval("A", 150, 450), None, hits
        )
        self.assertEqual(confidence, "High")

    def test_disagreeing_evidence_stays_separate_and_requires_review(self):
        synteny = Interval("chr1", 100, 200, evidence="synteny")
        marker = Interval("chr1", 500, 600, evidence="markers")
        hits = [
            Hit(str(i), "chr1", 100 + i, 101 + i, "+", 99, 100, source_start=i)
            for i in range(5)
        ]
        confidence, _reasons, warnings = score_confidence(
            "forward", synteny, marker, None, hits
        )
        self.assertEqual(confidence, "Manual check")
        self.assertEqual(len(combine_intervals(synteny, marker, None)), 2)
        self.assertIn("inconsistent", " ".join(warnings))


class OrientationTests(unittest.TestCase):
    """Preserve/validate interval orientation when reconciling independent evidence (#19)."""

    def _hits(self, count, strand="+"):
        return [
            Hit(str(i), "chr1", 100 + i, 101 + i, strand, 99, 100, source_start=i)
            for i in range(count)
        ]

    def test_reverse_synteny_plus_marker_stays_reverse(self):
        synteny = Interval("chr1", 100, 300, "-", "synteny")
        marker = Interval("chr1", 150, 250, ".", "markers")
        result = combine_intervals(synteny, marker, None)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].strand, "-")

    def test_forward_synteny_plus_marker_stays_forward(self):
        synteny = Interval("chr1", 100, 300, "+", "synteny")
        marker = Interval("chr1", 150, 250, ".", "markers")
        result = combine_intervals(synteny, marker, None)
        self.assertEqual(result[0].strand, "+")

    def test_marker_only_unresolved_orientation(self):
        marker = Interval("chr1", 150, 250, ".", "markers")
        self.assertEqual(combine_intervals(None, marker, None)[0].strand, ".")
        _confidence, _reasons, warnings = score_confidence(
            "partial", None, marker, None, self._hits(5, ".")
        )
        self.assertIn("orientation", " ".join(warnings).lower())

    def test_opposite_synteny_and_liftover_conflict_not_high(self):
        synteny = Interval("chr1", 100, 300, "-", "synteny")
        liftover = Interval("chr1", 120, 280, "+", "liftover")
        confidence, _reasons, warnings = score_confidence(
            "reverse", synteny, None, liftover, self._hits(6)
        )
        self.assertEqual(confidence, "Manual check")
        self.assertIn("orientation", " ".join(warnings).lower())
        self.assertEqual(len(combine_intervals(synteny, None, liftover)), 2)

    def test_agreeing_informative_evidence_combines(self):
        synteny = Interval("chr1", 100, 300, "+", "synteny")
        liftover = Interval("chr1", 120, 280, "+", "liftover")
        confidence, _reasons, _warnings = score_confidence(
            "forward", synteny, None, liftover, self._hits(6)
        )
        self.assertEqual(confidence, "High")
        result = combine_intervals(synteny, None, liftover)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].strand, "+")

    def test_non_overlapping_evidence_stays_separate(self):
        synteny = Interval("chr1", 100, 200, "+", "synteny")
        marker = Interval("chr1", 500, 600, ".", "markers")
        self.assertEqual(len(combine_intervals(synteny, marker, None)), 2)

    def test_orientation_audit_records_provenance(self):
        synteny = Interval("chr1", 100, 300, "-", "synteny")
        marker = Interval("chr1", 150, 250, ".", "markers")
        audit = orientation_audit(synteny, marker, None)
        self.assertEqual(audit["strand"], "-")
        self.assertEqual(audit["orientation_evidence"], ["synteny"])
        self.assertEqual(audit["uninformative_orientation_evidence"], ["markers"])
        self.assertFalse(audit["conflict"])


if __name__ == "__main__":
    unittest.main()
