import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

from qtlift.reporting import write_outputs


class ReportingAuditTests(unittest.TestCase):
    """anchors.tsv provenance columns and orientation provenance in the report (#32)."""

    def _summary(self):
        return {
            "job_id": "qtlift-abcdef0123", "name": "audit", "confidence": "Medium",
            "source_label": "B chr1:1-2", "final_label": "A chr1:2-3",
            "warnings": [], "reasons": ["r"], "evidence": {},
            "orientation": {"strand": "-", "orientation_evidence": ["synteny"],
                            "uninformative_orientation_evidence": ["markers"], "conflict": False},
            "anchors": [{"id": "g1", "contig": "chr1", "start": 100, "end": 200, "strand": "-",
                         "role": "peak", "cds": [(100, 200)], "transcript_id": "g1.t2",
                         "sequence_source": "cds"}],
        }

    def test_anchors_tsv_carries_transcript_provenance(self):
        with tempfile.TemporaryDirectory() as d:
            files = write_outputs(d, self._summary(), [], [], {})
            self.assertIn("anchors.tsv", files)
            tsv = (Path(d) / "anchors.tsv").read_text(encoding="utf-8")
        header, row = tsv.splitlines()[0], tsv.splitlines()[1]
        self.assertIn("transcript_id", header)
        self.assertIn("sequence_source", header)
        self.assertIn("g1.t2", row)
        self.assertIn("cds", row)

    def test_report_shows_orientation_provenance(self):
        with tempfile.TemporaryDirectory() as d:
            write_outputs(d, self._summary(), [], [], {})
            report = (Path(d) / "report.html").read_text(encoding="utf-8")
        self.assertIn("Orientation provenance", report)
        self.assertIn("synteny", report)
        self.assertIn("markers", report)

    def test_missing_audit_fields_do_not_break_writing(self):
        # A minimal summary (no anchors/orientation) must still produce the report.
        minimal = {"job_id": "qtlift-0000000000", "confidence": "Low", "source_label": "B c:1-2",
                   "final_label": "A c:2-3", "warnings": [], "reasons": [], "evidence": {}}
        with tempfile.TemporaryDirectory() as d:
            files = write_outputs(d, minimal, [], [], {})
            self.assertIn("anchors.tsv", files)
            self.assertIn("report.html", files)


if __name__ == "__main__":
    unittest.main()
