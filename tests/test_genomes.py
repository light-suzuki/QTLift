import sys
import tempfile
import unittest
from pathlib import Path

from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

from qtlift.genomes import anchor_sequence, detect_genomes, genes_in_interval, validate_interval
from qtlift.models import Gene


class GenomeLibraryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from scripts.create_sample_data import main

        main()
        cls.root = ROOT / "sample_data" / "genomes"

    def test_detection(self):
        rows = detect_genomes(self.root)
        self.assertEqual([x["name"] for x in rows], ["RefA", "RefB"])
        self.assertTrue(all(x["gene_count"] == 8 for x in rows))

    def test_interval_validation(self):
        validate_interval(1, 10, 20, 5)
        with self.assertRaises(ValueError):
            validate_interval(20, 10)
        with self.assertRaises(ValueError):
            validate_interval(1, 30, 20)

    def test_gene_extraction(self):
        genes = genes_in_interval(self.root / "RefB" / "refB.gff3", "Chr1", 100, 850)
        self.assertGreaterEqual(len(genes), 6)
        self.assertTrue(genes[0].cds)


class Gff3ParserTests(unittest.TestCase):
    """Species-agnostic GFF3 gene->transcript->CDS resolution and include_outside contract."""

    def _gff(self, directory, rows):
        path = Path(directory) / "ann.gff3"
        path.write_text(
            "##gff-version 3\n"
            + "\n".join("\t".join(str(column) for column in row) for row in rows)
            + "\n",
            encoding="utf-8",
        )
        return path

    def test_edge_gene_include_outside_contract(self):
        with tempfile.TemporaryDirectory() as directory:
            gff = self._gff(directory, [
                ("Chr1", ".", "gene", 50, 250, ".", "+", ".", "ID=g1"),
                ("Chr1", ".", "CDS", 50, 70, ".", "+", "0", "Parent=g1"),
                ("Chr1", ".", "CDS", 100, 120, ".", "+", "0", "Parent=g1"),
                ("Chr1", ".", "CDS", 200, 220, ".", "+", "0", "Parent=g1"),
            ])
            self.assertEqual(genes_in_interval(gff, "Chr1", 90, 180, include_outside=False), [])
            included = genes_in_interval(gff, "Chr1", 90, 180, include_outside=True)
            self.assertEqual([x.id for x in included], ["g1"])
            self.assertEqual(included[0].cds, [(50, 70), (100, 120), (200, 220)])

    def test_contained_gene_identical_in_both_modes(self):
        with tempfile.TemporaryDirectory() as directory:
            gff = self._gff(directory, [
                ("Chr1", ".", "gene", 100, 200, ".", "+", ".", "ID=g1"),
                ("Chr1", ".", "CDS", 100, 130, ".", "+", "0", "Parent=g1"),
                ("Chr1", ".", "CDS", 160, 200, ".", "+", "0", "Parent=g1"),
            ])
            contained = genes_in_interval(gff, "Chr1", 50, 300, include_outside=False)
            included = genes_in_interval(gff, "Chr1", 50, 300, include_outside=True)
            self.assertEqual([x.id for x in contained], ["g1"])
            self.assertEqual([x.id for x in included], ["g1"])
            self.assertEqual(contained[0].cds, included[0].cds)
            self.assertEqual(contained[0].cds, [(100, 130), (160, 200)])

    def test_direct_gene_to_cds(self):
        with tempfile.TemporaryDirectory() as directory:
            gff = self._gff(directory, [
                ("Chr1", ".", "gene", 100, 200, ".", "+", ".", "ID=g1"),
                ("Chr1", ".", "CDS", 100, 200, ".", "+", "0", "Parent=g1"),
            ])
            genes = genes_in_interval(gff, "Chr1", 1, 300, include_outside=True)
            self.assertEqual(genes[0].cds, [(100, 200)])
            self.assertEqual(genes[0].sequence_source, "cds")

    def test_implicit_transcript_suffix_parent(self):
        with tempfile.TemporaryDirectory() as directory:
            for suffix in (".t1", ".1", "-T1"):
                gff = self._gff(directory, [
                    ("Chr1", ".", "gene", 100, 300, ".", "+", ".", "ID=g1"),
                    ("Chr1", ".", "CDS", 100, 150, ".", "+", "0", f"Parent=g1{suffix}"),
                    ("Chr1", ".", "CDS", 200, 300, ".", "+", "0", f"Parent=g1{suffix}"),
                ])
                genes = genes_in_interval(gff, "Chr1", 1, 400, include_outside=True)
                self.assertEqual(genes[0].cds, [(100, 150), (200, 300)], suffix)
                self.assertEqual(genes[0].sequence_source, "cds", suffix)

    def test_unrelated_transcript_ids_resolve(self):
        with tempfile.TemporaryDirectory() as directory:
            gff = self._gff(directory, [
                ("Chr1", ".", "gene", 100, 200, ".", "+", ".", "ID=gene-LOC123"),
                ("Chr1", ".", "mRNA", 100, 200, ".", "+", ".", "ID=rna-XM_001;Parent=gene-LOC123"),
                ("Chr1", ".", "CDS", 100, 130, ".", "+", "0", "Parent=rna-XM_001"),
                ("Chr1", ".", "CDS", 160, 200, ".", "+", "0", "Parent=rna-XM_001"),
            ])
            genes = genes_in_interval(gff, "Chr1", 1, 300, include_outside=True)
            self.assertEqual([x.id for x in genes], ["gene-LOC123"])
            self.assertEqual(genes[0].cds, [(100, 130), (160, 200)])
            self.assertEqual(genes[0].transcript_id, "rna-XM_001")

    def test_two_isoforms_pick_longest_no_concatenation(self):
        with tempfile.TemporaryDirectory() as directory:
            gff = self._gff(directory, [
                ("Chr1", ".", "gene", 100, 300, ".", "+", ".", "ID=g1"),
                ("Chr1", ".", "mRNA", 100, 300, ".", "+", ".", "ID=g1.t1;Parent=g1"),
                ("Chr1", ".", "CDS", 100, 150, ".", "+", "0", "Parent=g1.t1"),
                ("Chr1", ".", "mRNA", 100, 300, ".", "+", ".", "ID=g1.t2;Parent=g1"),
                ("Chr1", ".", "CDS", 100, 150, ".", "+", "0", "Parent=g1.t2"),
                ("Chr1", ".", "CDS", 200, 300, ".", "+", "0", "Parent=g1.t2"),
            ])
            genes = genes_in_interval(gff, "Chr1", 1, 400, include_outside=True)
            self.assertEqual(genes[0].transcript_id, "g1.t2")
            self.assertEqual(genes[0].cds, [(100, 150), (200, 300)])

    def test_canonical_tag_overrides_length(self):
        with tempfile.TemporaryDirectory() as directory:
            gff = self._gff(directory, [
                ("Chr1", ".", "gene", 100, 300, ".", "+", ".", "ID=g1"),
                ("Chr1", ".", "mRNA", 100, 300, ".", "+", ".", "ID=g1.long;Parent=g1"),
                ("Chr1", ".", "CDS", 100, 300, ".", "+", "0", "Parent=g1.long"),
                ("Chr1", ".", "mRNA", 100, 180, ".", "+", ".", "ID=g1.canon;Parent=g1;tag=Ensembl_canonical"),
                ("Chr1", ".", "CDS", 100, 180, ".", "+", "0", "Parent=g1.canon"),
            ])
            genes = genes_in_interval(gff, "Chr1", 1, 400, include_outside=True)
            self.assertEqual(genes[0].transcript_id, "g1.canon")
            self.assertEqual(genes[0].cds, [(100, 180)])

    def test_multiple_parents_and_percent_encoding(self):
        with tempfile.TemporaryDirectory() as directory:
            gff = self._gff(directory, [
                ("Chr1", ".", "gene", 100, 200, ".", "+", ".", "ID=gene%3A1"),
                ("Chr1", ".", "mRNA", 100, 200, ".", "+", ".", "ID=t1;Parent=gene%3A1"),
                ("Chr1", ".", "mRNA", 100, 200, ".", "+", ".", "ID=t2;Parent=gene%3A1"),
                ("Chr1", ".", "CDS", 100, 200, ".", "+", "0", "Parent=t1,t2"),
            ])
            genes = genes_in_interval(gff, "Chr1", 1, 300, include_outside=True)
            self.assertEqual([x.id for x in genes], ["gene:1"])
            self.assertEqual(genes[0].transcript_id, "t1")
            self.assertEqual(genes[0].cds, [(100, 200)])

    def test_no_cds_falls_back_to_whole_gene(self):
        with tempfile.TemporaryDirectory() as directory:
            gff = self._gff(directory, [
                ("Chr1", ".", "gene", 100, 200, ".", "+", ".", "ID=g1"),
                ("Chr1", ".", "mRNA", 100, 200, ".", "+", ".", "ID=g1.t1;Parent=g1"),
                ("Chr1", ".", "exon", 100, 200, ".", "+", ".", "ID=g1.e1;Parent=g1.t1"),
            ])
            genes = genes_in_interval(gff, "Chr1", 1, 300, include_outside=True)
            self.assertEqual(genes[0].sequence_source, "gene")
            self.assertIsNone(genes[0].transcript_id)
            self.assertEqual(genes[0].cds, [])

    def test_minus_strand_anchor_is_biological_order(self):
        with tempfile.TemporaryDirectory() as directory:
            sequence = "A" * 10 + "CCCCCGGGGG" + "T" * 10 + "ACGTACGTAC" + "A" * 10
            fasta = Path(directory) / "g.fa"
            fasta.write_text(">Chr1\n" + sequence + "\n", encoding="ascii")
            gene = Gene("g1", "Chr1", 11, 40, "-", cds=[(11, 20), (31, 40)])
            result = anchor_sequence(fasta, gene)
            exon1, exon2 = sequence[10:20], sequence[30:40]
            self.assertEqual(result, str(Seq(exon1 + exon2).reverse_complement()))
            self.assertTrue(result.startswith(str(Seq(exon2).reverse_complement())))


if __name__ == "__main__":
    unittest.main()
