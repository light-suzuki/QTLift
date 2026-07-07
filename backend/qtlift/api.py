from __future__ import annotations

import json
import os
import shutil
import re
from pathlib import Path

from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

from .genomes import LOCAL_CONFIG, detect_genomes
from .pipeline import run_job
from .tools import detect_tools
from .providers import PROVIDERS, provider_status
from .jobs import JobManager

ROOT = Path(__file__).resolve().parents[2]
DATA = Path(os.environ.get("QTLIFT_DATA", ROOT/"data")); JOBS = DATA/"jobs"; JOBS.mkdir(parents=True, exist_ok=True)
app = FastAPI(title="QTLift API", version="1.0.0")
manager = JobManager(JOBS)

class LibraryRequest(BaseModel): genome_root: str
class ContigRequest(BaseModel):
    genome_root: str; reference: str; query: str = ""; offset: int = Field(default=0, ge=0); limit: int = Field(default=100, ge=1, le=500)
class JobRequest(BaseModel):
    genome_root: str; target_ref: str; source_ref: str; contig: str; start: int = Field(ge=1); end: int = Field(ge=1)
    peak: int | None = None; name: str = Field(default="", max_length=200); preset: str = "Standard"
    params: dict = Field(default_factory=dict); markers: dict[str,str] = Field(default_factory=dict)
    tool_paths: dict = Field(default_factory=dict); mapping_backend: str = "auto"
    provider_options: dict = Field(default_factory=dict); all_genes: bool = False; target_contig: str = ""

@app.get("/api/health")
def health(): return {"status":"ok","version":"1.0.0","sample_root":str(ROOT/"sample_data"/"genomes"),"local_library_available":LOCAL_CONFIG.is_file()}

@app.get("/api/capabilities")
def capabilities(): return {"tools":detect_tools(),"providers":PROVIDERS,"provider_status":provider_status(),"liftover_enabled":os.environ.get("QTLIFT_ENABLE_LIFTOVER") == "1"}

@app.post("/api/genomes/scan")
def scan(req: LibraryRequest):
    try:
        genomes = detect_genomes(req.genome_root)
        compact = []
        for genome in genomes:
            row = dict(genome)
            row["contig_count"] = len(genome["contigs"])
            row["assembly_size"] = sum(item["length"] for item in genome["contigs"])
            row["contigs"] = genome["contigs"][:200]
            compact.append(row)
        return {"genomes": compact}
    except Exception as exc: raise HTTPException(400,str(exc)) from exc

@app.post("/api/genomes/contigs")
def contigs(req: ContigRequest):
    try:
        genome = next((x for x in detect_genomes(req.genome_root) if x["name"] == req.reference), None)
        if not genome: raise ValueError(f"Unknown reference: {req.reference}")
        query = req.query.casefold().strip()
        rows = [x for x in genome["contigs"] if not query or query in x["name"].casefold()]
        return {"contigs": rows[req.offset:req.offset+req.limit], "total": len(rows), "offset": req.offset, "limit": req.limit}
    except Exception as exc: raise HTTPException(400,str(exc)) from exc

@app.post("/api/jobs")
def create_job(req: JobRequest):
    try: return manager.submit(req.model_dump())
    except Exception as exc: raise HTTPException(400,str(exc)) from exc

@app.get("/api/jobs")
def list_jobs():
    rows=[]
    for p in JOBS.glob("*/summary.json"):
        try: rows.append(json.loads(p.read_text(encoding="utf-8")))
        except Exception: pass
    rows.sort(key=lambda r: r.get("created_at",""), reverse=True)
    return {"jobs":rows}

@app.get("/api/jobs/{job_id}")
def get_job(job_id: str):
    if not re.fullmatch(r"qtlift-[a-f0-9]{10}", job_id): raise HTTPException(404,"Job not found")
    try: return manager.get(job_id)
    except FileNotFoundError: raise HTTPException(404,"Job not found")

@app.post("/api/jobs/{job_id}/cancel")
def cancel_job(job_id: str):
    if not re.fullmatch(r"qtlift-[a-f0-9]{10}", job_id): raise HTTPException(404,"Job not found")
    try: return manager.cancel(job_id)
    except FileNotFoundError: raise HTTPException(404,"Job not found")

@app.delete("/api/jobs/{job_id}")
def delete_job(job_id: str):
    if not re.fullmatch(r"qtlift-[a-f0-9]{10}", job_id): raise HTTPException(404,"Job not found")
    base=(JOBS/job_id).resolve()
    if JOBS.resolve() not in base.parents or not base.is_dir(): raise HTTPException(404,"Job not found")
    row = manager.get(job_id)
    if row.get("status") in ("queued","running","cancelling"): raise HTTPException(409,"Cancel the running job before deleting it.")
    shutil.rmtree(base)
    return {"deleted":job_id}

@app.get("/api/jobs/{job_id}/files/{file_path:path}")
def download(job_id: str,file_path: str):
    if not re.fullmatch(r"qtlift-[a-f0-9]{10}", job_id): raise HTTPException(404,"Job not found")
    base=(JOBS/job_id).resolve(); path=(base/file_path).resolve()
    if base not in path.parents or not path.is_file(): raise HTTPException(404,"File not found")
    return FileResponse(path,filename=path.name)

dist=ROOT/"frontend"/"dist"
if dist.exists():
    app.mount("/assets",StaticFiles(directory=dist/"assets"),name="assets")
    @app.get("/{path:path}",include_in_schema=False)
    def spa(path:str): return FileResponse(dist/"index.html")
