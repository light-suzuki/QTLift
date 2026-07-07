import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/"backend"))
from qtlift.pipeline import run_job
from create_sample_data import motif
payload={"job_id":"sample-forward","genome_root":str(ROOT/"sample_data"/"genomes"),"target_ref":"RefA","source_ref":"RefB","contig":"Chr1","start":100,"end":850,"peak":450,"name":"Artificial forward-collinear QTL","preset":"Standard","markers":{"left_flanking":motif(1),"peak":motif(4),"right_flanking":motif(7)}}
result=run_job(payload,ROOT/"data"/"jobs");print(json.dumps({k:result[k] for k in ('job_id','confidence','source_label','final_label','synteny_state','warnings','files')},indent=2))
