import shutil
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

from qtlift.pipeline import run_job


class PipelineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from scripts.create_sample_data import main

        main()
        cls.root = ROOT / "sample_data" / "genomes"

    @unittest.skipUnless(
        shutil.which("wsl.exe"),
        "asserts the WSL BLAST backend; unavailable on non-Windows CI",
    )
    def test_sample_pipeline(self):
        from scripts.create_sample_data import motif

        with tempfile.TemporaryDirectory() as directory:
            payload = {
                "job_id": "test",
                "genome_root": str(self.root),
                "target_ref": "RefA",
                "source_ref": "RefB",
                "contig": "Chr1",
                "start": 100,
                "end": 850,
                "peak": 450,
                "name": "test",
                "preset": "Standard",
                "markers": {
                    "left": motif(1),
                    "peak": motif(4),
                    "right": motif(7),
                },
            }
            result = run_job(payload, directory)
            self.assertEqual(result["synteny_state"], "forward")
            self.assertIn(result["confidence"], ("High", "Medium"))
            self.assertTrue((Path(directory) / "test" / "report.html").exists())
            self.assertIn("Liftover", " ".join(result["warnings"]))
            self.assertEqual(result["effective_backend"], "wsl")


if __name__ == "__main__":
    unittest.main()
