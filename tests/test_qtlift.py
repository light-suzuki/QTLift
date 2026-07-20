import gzip,json,shutil,sys,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/"backend"))
from qtlift.analysis import combine_intervals,evaluate_synteny,orientation_audit,score_confidence,select_anchors
from qtlift.genomes import anchor_sequence,detect_genomes,genes_in_interval,validate_interval
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
 def test_liftover_merges_across_chunk_boundary(self):
  with tempfile.TemporaryDirectory() as d:
   paf=Path(d)/'pair.paf.gz'
   with gzip.open(paf,'wt') as h:
    h.write('QTLIFT|Chr1|0\t20000000\t0\t20000000\t+\tTarget1\t40000000\t0\t20000000\t20000000\t20000000\t60\n')
    h.write('QTLIFT|Chr1|20000000\t20000000\t0\t20000000\t+\tTarget1\t40000000\t20000000\t40000000\t20000000\t20000000\t60\n')
   iv,w=lift_interval(paf,'Chr1',19500000,20500000,'Target1')
   self.assertEqual((iv.contig,iv.start,iv.end,iv.strand),('Target1',19500000,20500000,'+'));self.assertEqual(iv.end-iv.start+1,1000001)
 def test_liftover_opposite_strand_chunks_stay_ambiguous(self):
  with tempfile.TemporaryDirectory() as d:
   paf=Path(d)/'pair.paf.gz'
   with gzip.open(paf,'wt') as h:
    h.write('QTLIFT|Chr1|0\t20000000\t0\t20000000\t+\tTarget1\t40000000\t0\t20000000\t20000000\t20000000\t60\n')
    h.write('QTLIFT|Chr1|20000000\t20000000\t0\t20000000\t-\tTarget1\t40000000\t20000000\t40000000\t20000000\t20000000\t60\n')
   iv,w=lift_interval(paf,'Chr1',19500000,20500000,'Target1')
   self.assertIn('ambiguous',' '.join(w).lower());self.assertLess(iv.end-iv.start+1,1000001)
 def test_liftover_distant_same_strand_blocks_not_merged(self):
  with tempfile.TemporaryDirectory() as d:
   paf=Path(d)/'pair.paf.gz'
   with gzip.open(paf,'wt') as h:
    h.write('QTLIFT|Chr1|0\t20000000\t0\t20000000\t+\tTarget1\t200000000\t0\t20000000\t20000000\t20000000\t60\n')
    h.write('QTLIFT|Chr1|20000000\t20000000\t0\t20000000\t+\tTarget1\t200000000\t100000000\t120000000\t20000000\t20000000\t60\n')
   iv,w=lift_interval(paf,'Chr1',19500000,20500000,'Target1')
   self.assertLess(iv.end-iv.start+1,1000001);self.assertIn('partial',' '.join(w).lower())
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

class Gff3ParserTests(unittest.TestCase):
 """Species-agnostic GFF3 gene->transcript->CDS resolution and the include_outside contract (#18, #14)."""
 def _gff(self,d,rows):
  p=Path(d)/'ann.gff3';p.write_text("##gff-version 3\n"+"\n".join("\t".join(str(c) for c in r) for r in rows)+"\n",encoding="utf-8");return p
 def test_edge_gene_include_outside_contract(self):
  with tempfile.TemporaryDirectory() as d:
   g=self._gff(d,[("Chr1",".","gene",50,250,".","+",".","ID=g1"),("Chr1",".","CDS",50,70,".","+","0","Parent=g1"),("Chr1",".","CDS",100,120,".","+","0","Parent=g1"),("Chr1",".","CDS",200,220,".","+","0","Parent=g1")])
   self.assertEqual(genes_in_interval(g,"Chr1",90,180,include_outside=False),[])
   inc=genes_in_interval(g,"Chr1",90,180,include_outside=True)
   self.assertEqual([x.id for x in inc],["g1"]);self.assertEqual(inc[0].cds,[(50,70),(100,120),(200,220)])
 def test_contained_gene_identical_in_both_modes(self):
  with tempfile.TemporaryDirectory() as d:
   g=self._gff(d,[("Chr1",".","gene",100,200,".","+",".","ID=g1"),("Chr1",".","CDS",100,130,".","+","0","Parent=g1"),("Chr1",".","CDS",160,200,".","+","0","Parent=g1")])
   a=genes_in_interval(g,"Chr1",50,300,include_outside=False);b=genes_in_interval(g,"Chr1",50,300,include_outside=True)
   self.assertEqual([x.id for x in a],["g1"]);self.assertEqual([x.id for x in b],["g1"]);self.assertEqual(a[0].cds,b[0].cds);self.assertEqual(a[0].cds,[(100,130),(160,200)])
 def test_direct_gene_to_cds(self):
  with tempfile.TemporaryDirectory() as d:
   g=self._gff(d,[("Chr1",".","gene",100,200,".","+",".","ID=g1"),("Chr1",".","CDS",100,200,".","+","0","Parent=g1")])
   genes=genes_in_interval(g,"Chr1",1,300,include_outside=True)
   self.assertEqual(genes[0].cds,[(100,200)]);self.assertEqual(genes[0].sequence_source,"cds")
 def test_implicit_transcript_suffix_parent(self):
  # CDS Parent names a transcript id (g1.t1/g1.1/g1-T1) with no explicit transcript row (#23 review P1).
  with tempfile.TemporaryDirectory() as d:
   for suffix in (".t1",".1","-T1"):
    g=self._gff(d,[("Chr1",".","gene",100,300,".","+",".","ID=g1"),("Chr1",".","CDS",100,150,".","+","0",f"Parent=g1{suffix}"),("Chr1",".","CDS",200,300,".","+","0",f"Parent=g1{suffix}")])
    genes=genes_in_interval(g,"Chr1",1,400,include_outside=True)
    self.assertEqual(genes[0].cds,[(100,150),(200,300)],suffix);self.assertEqual(genes[0].sequence_source,"cds",suffix)
 def test_unrelated_transcript_ids_resolve(self):
  with tempfile.TemporaryDirectory() as d:
   g=self._gff(d,[("Chr1",".","gene",100,200,".","+",".","ID=gene-LOC123"),("Chr1",".","mRNA",100,200,".","+",".","ID=rna-XM_001;Parent=gene-LOC123"),("Chr1",".","CDS",100,130,".","+","0","Parent=rna-XM_001"),("Chr1",".","CDS",160,200,".","+","0","Parent=rna-XM_001")])
   genes=genes_in_interval(g,"Chr1",1,300,include_outside=True)
   self.assertEqual([x.id for x in genes],["gene-LOC123"]);self.assertEqual(genes[0].cds,[(100,130),(160,200)]);self.assertEqual(genes[0].transcript_id,"rna-XM_001")
 def test_two_isoforms_pick_longest_no_concatenation(self):
  with tempfile.TemporaryDirectory() as d:
   g=self._gff(d,[("Chr1",".","gene",100,300,".","+",".","ID=g1"),("Chr1",".","mRNA",100,300,".","+",".","ID=g1.t1;Parent=g1"),("Chr1",".","CDS",100,150,".","+","0","Parent=g1.t1"),("Chr1",".","mRNA",100,300,".","+",".","ID=g1.t2;Parent=g1"),("Chr1",".","CDS",100,150,".","+","0","Parent=g1.t2"),("Chr1",".","CDS",200,300,".","+","0","Parent=g1.t2")])
   genes=genes_in_interval(g,"Chr1",1,400,include_outside=True)
   self.assertEqual(genes[0].transcript_id,"g1.t2");self.assertEqual(genes[0].cds,[(100,150),(200,300)])
 def test_canonical_tag_overrides_length(self):
  with tempfile.TemporaryDirectory() as d:
   g=self._gff(d,[("Chr1",".","gene",100,300,".","+",".","ID=g1"),("Chr1",".","mRNA",100,300,".","+",".","ID=g1.long;Parent=g1"),("Chr1",".","CDS",100,300,".","+","0","Parent=g1.long"),("Chr1",".","mRNA",100,180,".","+",".","ID=g1.canon;Parent=g1;tag=Ensembl_canonical"),("Chr1",".","CDS",100,180,".","+","0","Parent=g1.canon")])
   genes=genes_in_interval(g,"Chr1",1,400,include_outside=True)
   self.assertEqual(genes[0].transcript_id,"g1.canon");self.assertEqual(genes[0].cds,[(100,180)])
 def test_multiple_parents_and_percent_encoding(self):
  with tempfile.TemporaryDirectory() as d:
   g=self._gff(d,[("Chr1",".","gene",100,200,".","+",".","ID=gene%3A1"),("Chr1",".","mRNA",100,200,".","+",".","ID=t1;Parent=gene%3A1"),("Chr1",".","mRNA",100,200,".","+",".","ID=t2;Parent=gene%3A1"),("Chr1",".","CDS",100,200,".","+","0","Parent=t1,t2")])
   genes=genes_in_interval(g,"Chr1",1,300,include_outside=True)
   self.assertEqual([x.id for x in genes],["gene:1"]);self.assertEqual(genes[0].transcript_id,"t1");self.assertEqual(genes[0].cds,[(100,200)])
 def test_no_cds_falls_back_to_whole_gene(self):
  with tempfile.TemporaryDirectory() as d:
   g=self._gff(d,[("Chr1",".","gene",100,200,".","+",".","ID=g1"),("Chr1",".","mRNA",100,200,".","+",".","ID=g1.t1;Parent=g1"),("Chr1",".","exon",100,200,".","+",".","ID=g1.e1;Parent=g1.t1")])
   genes=genes_in_interval(g,"Chr1",1,300,include_outside=True)
   self.assertEqual(genes[0].sequence_source,"gene");self.assertIsNone(genes[0].transcript_id);self.assertEqual(genes[0].cds,[])
 def test_minus_strand_anchor_is_biological_order(self):
  from Bio.Seq import Seq
  with tempfile.TemporaryDirectory() as d:
   seq="A"*10+"CCCCCGGGGG"+"T"*10+"ACGTACGTAC"+"A"*10
   fa=Path(d)/"g.fa";fa.write_text(">Chr1\n"+seq+"\n",encoding="ascii")
   gene=Gene("g1","Chr1",11,40,"-",cds=[(11,20),(31,40)])
   got=anchor_sequence(fa,gene)
   exon1,exon2=seq[10:20],seq[30:40]
   self.assertEqual(got,str(Seq(exon1+exon2).reverse_complement()));self.assertTrue(got.startswith(str(Seq(exon2).reverse_complement())))

class OrientationTests(unittest.TestCase):
 """Preserve/validate interval orientation when reconciling independent evidence (#19)."""
 def _hits(self,n,strand='+'):
  return [Hit(str(i),'chr1',100+i,101+i,strand,99,100,source_start=i) for i in range(n)]
 def test_reverse_synteny_plus_marker_stays_reverse(self):
  syn=Interval('chr1',100,300,'-','synteny');mk=Interval('chr1',150,250,'.','markers')
  res=combine_intervals(syn,mk,None);self.assertEqual(len(res),1);self.assertEqual(res[0].strand,'-')
 def test_forward_synteny_plus_marker_stays_forward(self):
  syn=Interval('chr1',100,300,'+','synteny');mk=Interval('chr1',150,250,'.','markers')
  res=combine_intervals(syn,mk,None);self.assertEqual(res[0].strand,'+')
 def test_marker_only_unresolved_orientation(self):
  mk=Interval('chr1',150,250,'.','markers')
  self.assertEqual(combine_intervals(None,mk,None)[0].strand,'.')
  c,r,w=score_confidence('partial',None,mk,None,self._hits(5,'.'))
  self.assertIn('orientation',' '.join(w).lower())
 def test_opposite_synteny_and_liftover_conflict_not_high(self):
  syn=Interval('chr1',100,300,'-','synteny');lift=Interval('chr1',120,280,'+','liftover')
  c,r,w=score_confidence('reverse',syn,None,lift,self._hits(6))
  self.assertEqual(c,'Manual check');self.assertIn('orientation',' '.join(w).lower())
  self.assertEqual(len(combine_intervals(syn,None,lift)),2)
 def test_agreeing_informative_evidence_combines(self):
  syn=Interval('chr1',100,300,'+','synteny');lift=Interval('chr1',120,280,'+','liftover')
  c,r,w=score_confidence('forward',syn,None,lift,self._hits(6));self.assertEqual(c,'High')
  res=combine_intervals(syn,None,lift);self.assertEqual(len(res),1);self.assertEqual(res[0].strand,'+')
 def test_non_overlapping_evidence_stays_separate(self):
  syn=Interval('chr1',100,200,'+','synteny');mk=Interval('chr1',500,600,'.','markers')
  self.assertEqual(len(combine_intervals(syn,mk,None)),2)
 def test_orientation_audit_records_provenance(self):
  syn=Interval('chr1',100,300,'-','synteny');mk=Interval('chr1',150,250,'.','markers')
  aud=orientation_audit(syn,mk,None)
  self.assertEqual(aud['strand'],'-');self.assertEqual(aud['orientation_evidence'],['synteny'])
  self.assertEqual(aud['uninformative_orientation_evidence'],['markers']);self.assertFalse(aud['conflict'])
if __name__=='__main__':unittest.main()
