import shutil
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

from qtlift.pipeline import run_job


class AmbiguousReferenceTests(unittest.TestCase):
    """A reference folder with several FASTA/GFF files must block analysis (#47)."""

    def _genome(self, directory: Path, name: str, fasta_names: list[str], gff: bool = True):
        ref = directory / name
        ref.mkdir()
        for i, fasta_name in enumerate(fasta_names):
            (ref / fasta_name).write_text(">Chr1\n" + "ACGT" * 500 + f"\n{i}", encoding="ascii")
        if gff:
            (ref / "ann.gff3").write_text("##gff-version 3\nChr1\t.\tgene\t1\t100\t.\t+\t.\tID=g1\n", encoding="ascii")

    def _payload(self, root: Path, **overrides):
        payload = {"job_id": "amb", "genome_root": str(root), "target_ref": "RefA", "source_ref": "RefB",
                   "contig": "Chr1", "start": 1, "end": 500, "peak": 250}
        payload.update(overrides)
        return payload

    def test_multiple_source_fasta_blocks_job(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self._genome(root, "RefA", ["a.fa"])
            self._genome(root, "RefB", ["b1.fa", "b2.fa"])
            with self.assertRaises(ValueError) as cm:
                run_job(self._payload(root), directory)
            message = str(cm.exception)
            self.assertIn("RefB", message)
            self.assertIn("FASTA", message)
            self.assertIn("b2.fa", message)

    def test_multiple_source_gff_blocks_job(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self._genome(root, "RefA", ["a.fa"])
            self._genome(root, "RefB", ["b.fa"])
            ref_b = root / "RefB"
            (ref_b / "ann2.gff").write_text("##gff-version 3\nChr1\t.\tgene\t1\t100\t.\t+\t.\tID=g2\n", encoding="ascii")
            with self.assertRaises(ValueError) as cm:
                run_job(self._payload(root), directory)
            message = str(cm.exception)
            self.assertIn("RefB", message)
            self.assertIn("GFF", message)
            self.assertIn("ann2.gff", message)


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
