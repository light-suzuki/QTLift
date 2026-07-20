import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

from qtlift.reporting import write_outputs


class ReportingOutputTests(unittest.TestCase):
    def test_output_writing(self):
        with tempfile.TemporaryDirectory() as directory:
            summary = {
                "job_id": "x",
                "confidence": "Low",
                "source_label": "B c:1-2",
                "final_label": "A c:2-3",
                "warnings": ["w"],
                "reasons": ["r"],
                "evidence": {},
            }
            files = write_outputs(directory, summary, [], [], {})
            self.assertIn("report.html", files)
            self.assertIn("summary.json", files)


if __name__ == "__main__":
    unittest.main()
