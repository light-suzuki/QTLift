from __future__ import annotations

import os
import time
from threading import Event
from typing import Callable
from dataclasses import asdict
from datetime import datetime
from pathlib import Path
from uuid import uuid4

from .analysis import combine_intervals, evaluate_synteny, exact_hits, marker_interval, score_confidence, select_anchors
from .genomes import anchor_sequence, detect_genomes, genes_in_interval, sequence_slice, validate_interval
from .markers import parse_marker
from .models import Params
from .reporting import write_outputs
from .tools import detect_tools
from .providers import validate_provider
from .blast import blast_many, blast_subject
from .liftover import build_alignment_cache, lift_interval

# Auto maps every gene (all genes + orientation, the most accurate mode) when a source
# interval has at most this many, and falls back to spacing-driven anchor selection above it.
AUTO_ALL_GENES_MAX = 50


class JobCancelled(RuntimeError):
    pass


def run_job(payload: dict, jobs_root: str | Path, progress: Callable[[int, str], None] | None = None,
            cancel_event: Event | None = None) -> dict:
    _t0 = time.monotonic()
    def update(percent: int, stage: str) -> None:
        if cancel_event and cancel_event.is_set():
            raise JobCancelled("Job cancelled by user.")
        if progress:
            progress(percent, stage)
    update(2, "Loading genome library")
    libraries = {x["name"]: x for x in detect_genomes(payload["genome_root"])}
    target, source = libraries[payload["target_ref"]], libraries[payload["source_ref"]]
    start, end, peak = int(payload["start"]), int(payload["end"]), payload.get("peak")
    peak = int(peak) if peak not in (None, "") else None
    length = next(x["length"] for x in source["contigs"] if x["name"] == payload["contig"])
    validate_interval(start, end, length, peak)
    update(8, "Reading source annotation")
    params = Params.preset(payload.get("preset", "Standard"))
    if payload.get("params"): params = Params(**{**asdict(params), **payload["params"]})
    backend=payload.get("mapping_backend", "auto"); provider_options=payload.get("provider_options") or {}
    target_contig = (payload.get("target_contig") or "").strip() or None
    warnings: list[str] = validate_provider(backend, provider_options)
    tools = detect_tools(payload.get("tool_paths"))
    effective_backend = backend
    if backend == "auto":
        effective_backend = tools["blastn"].get("runtime") or ("wsl" if __import__('shutil').which('wsl.exe') else "exact")
    def map_query(qid,seq,source_start=None):
        if effective_backend in ("windows","wsl"):
            try: return blast_subject(qid,seq,target["fasta"],effective_backend,source_start,provider_options.get("wsl_distro"),target.get("blast_db"),target_contig,params.min_identity,params.min_coverage)
            except Exception as exc:
                if __import__('pathlib').Path(target["fasta"]).stat().st_size > 100_000_000:
                    warnings.append(f"{effective_backend} BLAST failed for {qid}: {exc}; exact fallback was disabled for this large target genome.")
                    return []
                warnings.append(f"{effective_backend} BLAST failed for {qid}: {exc}; exact fallback used.")
        return exact_hits(qid,seq,target["fasta"],target["contigs"],source_start,target_contig)
    gff_contig = source.get("gff_contig_aliases", {}).get(payload["contig"], payload["contig"])
    genes = genes_in_interval(source["gff"], gff_contig, start, end, params.include_outside_edge_genes)
    for gene in genes:
        gene.contig = payload["contig"]
    if len(genes) < 3: warnings.append(f"Too few genes in source interval ({len(genes)}).")
    preset_name = payload.get("preset", "Standard")
    if payload.get("all_genes"):
        anchors = genes
    elif preset_name == "Auto" and 0 < len(genes) <= AUTO_ALL_GENES_MAX:
        # Few enough genes to map every one: most accurate (all genes + orientation).
        anchors = genes
    else:
        anchors = select_anchors(genes, start, end, peak, params)
    if len(anchors) > 500:
        warnings.append(f"Selected {len(anchors)} genes to map; runtime may be long.")
    update(15, f"Selected {len(anchors)} anchor genes")
    fallback = [gene.id for gene in anchors if gene.sequence_source != "cds"]
    if fallback:
        warnings.append(f"No CDS model for {len(fallback)} anchor gene(s) ({', '.join(fallback[:5])}); whole-gene sequence was used.")
    anchor_hits = []
    anchor_queries = []
    for gene in anchors:
        seq = anchor_sequence(source["fasta"], gene)
        anchor_queries.append((gene.id, seq, gene.start))
    batched = None
    if effective_backend in ("windows", "wsl"):
        try:
            batched = blast_many(anchor_queries, target["fasta"], effective_backend, provider_options.get("wsl_distro"), target.get("blast_db"), target_contig, params.min_identity, params.min_coverage)
        except Exception as exc:
            warnings.append(f"{effective_backend} batched BLAST failed: {exc}; per-query fallback was used.")
    update(65, "Anchor mapping completed")
    for index, (gene, query) in enumerate(zip(anchors, anchor_queries)):
        update(65 + int(15 * (index / max(1, len(anchors)))), f"Evaluating anchor {gene.id}")
        hits = batched.get(gene.id, []) if batched is not None else map_query(*query)
        if len(hits) > params.multi_hit_threshold: warnings.append(f"Anchor {gene.id} has excessive multi-hits ({len(hits)}).")
        anchor_hits.extend(hits)
    synteny_state, synteny, synteny_warn = evaluate_synteny(anchor_hits); warnings += synteny_warn
    marker_hits = []
    for role, text in (payload.get("markers") or {}).items():
        marker, marker_warn = parse_marker(role, text); warnings += marker_warn
        if not marker: continue
        sequence = marker.sequence
        if not sequence and marker.contig and marker.start and marker.end:
            sequence = sequence_slice(source["fasta"], marker.contig, marker.start, marker.end)
        if sequence: marker_hits += map_query(marker.name, sequence, marker.start)
    marker_iv, marker_warn = marker_interval(marker_hits); warnings += marker_warn
    update(84, "Marker mapping completed")
    liftover = None
    liftover_cache = None
    liftover_enabled = os.environ.get("QTLIFT_ENABLE_LIFTOVER") == "1" or bool(provider_options.get("enable_liftover"))
    if not liftover_enabled:
        warnings.append("Liftover skipped: whole-genome minimap2 liftover is disabled.")
    elif not tools["minimap2"]["available"]:
        warnings.append("Liftover skipped: minimap2 is unavailable.")
    else:
        try:
            cache_root = Path(os.environ.get("QTLIFT_ALIGNMENT_CACHE", Path(jobs_root).parent / "alignment_cache"))
            liftover_cache, cache_hit = build_alignment_cache(source["fasta"], target["fasta"], cache_root,
                                                               tools["minimap2"]["path"], provider_options.get("wsl_distro"))
            liftover, lift_warn = lift_interval(liftover_cache, payload["contig"], start, end, target_contig)
            warnings += lift_warn
            if not cache_hit:
                warnings.append("Whole-genome minimap2 alignment cache was created for this reference pair.")
        except Exception as exc:
            warnings.append(f"Liftover failed: {exc}")
    if effective_backend == "exact":
        warnings.append("BLAST+ unavailable; exact-sequence fallback was used for sample-compatible mapping.")
    confidence, reasons, score_warn = score_confidence(synteny_state, synteny, marker_iv, liftover, anchor_hits); warnings += score_warn
    candidates = combine_intervals(synteny, marker_iv, liftover)
    final = None if confidence == "Manual check" and len(candidates) > 1 else (candidates[0] if candidates else None)
    if final and (final.end-final.start+1) > (end-start+1)*3: warnings.append("Target interval is much larger than the source interval.")
    job_id = payload.get("job_id") or f"qtlift-{uuid4().hex[:10]}"
    summary = {"job_id": job_id, "name": payload.get("name") or "Unnamed region", "status": "completed", "confidence": confidence,
               "created_at": payload.get("_created_at") or datetime.now().isoformat(timespec="seconds"), "duration_sec": round(time.monotonic()-_t0, 1), "progress": 100, "stage": "Completed", "target_contig": target_contig or "",
               "reasons": reasons, "warnings": list(dict.fromkeys(warnings)), "synteny_state": synteny_state,
               "source": {"reference": source["name"], "contig": payload["contig"], "start": start, "end": end, "peak": peak},
               "source_label": f"{source['name']} {payload['contig']}:{start:,}-{end:,}", "target_reference": target["name"],
               "candidates": [x.as_dict() for x in candidates], "final": final.as_dict() if final else None,
               "final_label": f"{target['name']} {final.contig}:{final.start:,}-{final.end:,}" if final else "Unavailable",
               "evidence": {"liftover": liftover.as_dict() if liftover else None, "markers": marker_iv.as_dict() if marker_iv else None, "synteny": synteny.as_dict() if synteny else None},
               "marker_hits": [asdict(x) for x in marker_hits], "anchor_hits": [asdict(x) for x in anchor_hits], "anchors": [asdict(x) for x in anchors],
               "params": asdict(params), "tools": tools, "mapping_backend": backend, "effective_backend": effective_backend,
               "liftover_cache": str(liftover_cache) if liftover_cache else None}
    job_dir = Path(jobs_root)/job_id
    update(95, "Writing reports")
    summary["files"] = write_outputs(job_dir, summary, summary["marker_hits"], summary["anchor_hits"], summary["params"], "Job completed\n")
    (job_dir/"summary.json").write_text(__import__('json').dumps(summary,indent=2,ensure_ascii=False),encoding='utf-8')
    return summary
