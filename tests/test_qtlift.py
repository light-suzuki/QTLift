import gzip,json,shutil,sys,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/"backend"))
from qtlift.analysis import combine_intervals,evaluate_synteny,score_confidence,select_anchors
from qtlift.genomes import detect_genomes,genes_in_interval,validate_interval
from qtlift.markers import parse_marker
from qtlift.models import Gene,Hit,Interval,Params
from qtlift.pipeline import run_job
from qtlift.reporting import write_outputs
from qtlift.liftover import lift_interval
from qtlift.blast import _collapse_loci,_find_blast_entry

class QTLiftTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  from scripts.create_sample_data import main;main();cls.root=ROOT/"sample_data"/"genomes"
 def test_detection(self):
  rows=detect_genomes(self.root);self.assertEqual([x['name'] for x in rows],['RefA','RefB']);self.assertTrue(all(x['gene_count']==8 for x in rows))
 def test_interval_validation(self):
  validate_interval(1,10,20,5)
  with self.assertRaises(ValueError):validate_interval(20,10)
  with self.assertRaises(ValueError):validate_interval(1,30,20)
 def test_gene_extraction(self):
  genes=genes_in_interval(self.root/'RefB'/'refB.gff3','Chr1',100,850);self.assertGreaterEqual(len(genes),6);self.assertTrue(genes[0].cds)
 def test_marker_parsing(self):
  m,w=parse_marker('left','>m1\nACGTACGTACGTACGTACGT');self.assertEqual(m.name,'m1');self.assertFalse(w)
  m,w=parse_marker('peak','ACGTACGTACGTACGTACGT');self.assertEqual(m.sequence,'ACGTACGTACGTACGTACGT')
  m,w=parse_marker('peak','mk=Chr1:10-20');self.assertIsNone(m);self.assertTrue(w)
  m,w=parse_marker('right','SSR-12');self.assertIsNone(m);self.assertTrue(w)
 def test_anchor_selection(self):
  genes=[Gene(f'g{i}','c',i*100,i*100+20) for i in range(1,11)];a=select_anchors(genes,100,1020,500,Params());self.assertLessEqual(len(a),16);self.assertIn('left-edge',{x.role for x in a})
 def test_synteny_forward_reverse_split(self):
  hs=[Hit(str(i),'A',i*100,i*100+10,'+',100,100,1,'x',i*90) for i in range(1,6)];state,iv,w=evaluate_synteny(hs);self.assertEqual(state,'forward');self.assertEqual(iv.contig,'A')
  rev=[Hit(str(i),'A',1000-i*100,1010-i*100,'-',100,100,1,'x',i*90) for i in range(1,6)];self.assertEqual(evaluate_synteny(rev)[0],'reverse')
 def test_confidence(self):
  iv=Interval('A',100,500);hs=[Hit(str(i),'A',i,i+1,'+',100,100) for i in range(5)];c,r,w=score_confidence('forward',iv,Interval('A',150,450),None,hs);self.assertEqual(c,'High')
 def test_output_writing(self):
  with tempfile.TemporaryDirectory() as d:
   s={'job_id':'x','confidence':'Low','source_label':'B c:1-2','final_label':'A c:2-3','warnings':['w'],'reasons':['r'],'evidence':{}}
   files=write_outputs(d,s,[],[],{});self.assertIn('report.html',files);self.assertIn('summary.json',files)
 def test_liftover_from_cached_paf(self):
  with tempfile.TemporaryDirectory() as d:
   paf=Path(d)/'pair.paf.gz'
   with gzip.open(paf,'wt') as h:h.write('Chr1\t1000\t99\t900\t+\tTarget1\t2000\t499\t1300\t780\t801\t60\ttp:A:P\tcg:Z:801M\n')
   iv,w=lift_interval(paf,'Chr1',200,800,'Target1')
   self.assertEqual((iv.contig,iv.start,iv.end,iv.strand),('Target1',600,1200,'+'));self.assertFalse(w)
 def test_blast_contig_entry_resolution(self):
  metadata='gnl|BL_ORD_ID|0\tchr1\ngnl|BL_ORD_ID|1\tchr2 description\n'
  self.assertEqual(_find_blast_entry(metadata,'chr2'),'gnl|BL_ORD_ID|1')
  self.assertIsNone(_find_blast_entry(metadata,'chr20'))
 def test_hsp_coverage_uses_query_union_and_separates_paralogs(self):
  hs=[Hit('q','chr1',100,149,'+',99,50,method='blastn',query_start=1,query_end=50,query_length=100),
      Hit('q','chr1',200,249,'+',97,50,method='blastn',query_start=51,query_end=100,query_length=100)]
  loci=_collapse_loci(hs);self.assertEqual(len(loci),1);self.assertEqual(loci[0].coverage,100)
  paralog=Hit('q','chr1',50000,50099,'+',96,100,method='blastn',query_start=1,query_end=100,query_length=100)
  self.assertEqual(len(_collapse_loci(hs+[paralog])),2)
  self.assertTrue(all(h.coverage < 70 for h in hs));self.assertGreaterEqual(loci[0].coverage,70)
 def test_disagreeing_evidence_stays_separate_and_requires_review(self):
  syn=Interval('chr1',100,200,evidence='synteny');marker=Interval('chr1',500,600,evidence='markers')
  hs=[Hit(str(i),'chr1',100+i,101+i,'+',99,100,source_start=i) for i in range(5)]
  confidence,reasons,warnings=score_confidence('forward',syn,marker,None,hs)
  self.assertEqual(confidence,'Manual check');self.assertEqual(len(combine_intervals(syn,marker,None)),2)
  self.assertIn('inconsistent',' '.join(warnings))
 @unittest.skipUnless(shutil.which("wsl.exe"), "asserts the WSL BLAST backend; unavailable on non-Windows CI")
 def test_sample_pipeline(self):
  with tempfile.TemporaryDirectory() as d:
   from scripts.create_sample_data import motif
   p={"job_id":"test","genome_root":str(self.root),"target_ref":"RefA","source_ref":"RefB","contig":"Chr1","start":100,"end":850,"peak":450,"name":"test","preset":"Standard","markers":{"left":motif(1),"peak":motif(4),"right":motif(7)}}
   r=run_job(p,d);self.assertEqual(r['synteny_state'],'forward');self.assertIn(r['confidence'],('High','Medium'));self.assertTrue((Path(d)/'test'/'report.html').exists());self.assertIn('Liftover'," ".join(r['warnings']))
   self.assertEqual(r['effective_backend'],'wsl')
if __name__=='__main__':unittest.main()
