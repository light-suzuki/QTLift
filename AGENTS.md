# QTLift — implementation guide for contributors and AI agents

This document is written so that a human **or an AI coding agent** can pick up QTLift and extend it safely. Read it before making changes. `SPEC.md` is the product source of truth; this file is the engineering contract.

QTLift was itself built by AI agents (Claude + Codex) via "vibe coding", so it is designed to be agent-friendly: small modules, explicit data shapes, fast tests, and a hard rule against fabricating results.

---

## 0. The non-negotiable contract

1. **Never fabricate mapping evidence.** If a tool, index, alignment, or annotation is missing, emit a warning and degrade gracefully — do not invent coordinates. Absent liftover/alignment ⇒ no BED, no interval, an explicit warning.
2. **Lower confidence or force manual review** when evidence is sparse, split, ambiguous, or inconsistent. It is correct to return "Manual check" and no `final` interval.
3. **Coordinates:** UI/API are **1-based inclusive**; BED is **0-based half-open**. Keep this straight at every boundary.
4. **Privacy:** never commit real genomes, job data, machine-specific catalogs, or private paths. `data/`, `alignment_cache/`, `output/`, and `.venv/` are git-ignored. Keep it that way.
5. **The app must always start**, even with zero external tools installed. Missing tools disable only their evidence step.

---

## 1. Architecture

```
backend/qtlift/        FastAPI + analysis pipeline (Python 3.13)
  api.py               HTTP routes; JobManager wiring; serves the built frontend
  jobs.py              Async job queue (ThreadPoolExecutor), progress, cancel, restart-recovery
  pipeline.py          run_job(): the end-to-end orchestration (the heart)
  genomes.py           Genome discovery, FASTA (.fai) slicing, streaming GFF parsing
  analysis.py          Anchor selection, exact-match fallback, synteny, confidence, interval combine
  blast.py             blastn adapters (Windows / WSL, -db or -subject), HSP parsing, locus collapse
  liftover.py          minimap2 whole-genome PAF cache + interval projection
  markers.py           Marker parsing (DNA sequence only)
  models.py            Dataclasses: Gene, Hit, Interval, Params (+ presets)
  providers.py         Mapping-provider registry & validation
  tools.py             External-tool detection (PATH / WSL)
  reporting.py         Job artifacts: report.html, TSVs, BED, params.json, logs
frontend/src/          React + Vite (TypeScript), single-file main.tsx
  main.tsx             All pages: Library, Setup, Monitor, Results, How-it-works, Settings
  locales/{en,ja}.json Bilingual UI strings
tests/test_*.py        unittest suite split by backend module (fast, no network, artificial sample)
scripts/               create_sample_data.py, run_sample.py, chunk_fasta.py, build_alignment_cache.py, ...
sample_data/genomes/   Artificial RefA/RefB (safe to redistribute)
```

### Pipeline data flow (`pipeline.run_job`)

```
detect_genomes → validate_interval → Params.preset
   → genes_in_interval (GFF)          [source B interval]
   → select_anchors  (or all genes)
   → BLAST anchor CDS onto target A   (blast_many / blast_subject; exact fallback)
   → evaluate_synteny (order+orientation → forward/reverse/partial/split/failed)
   → map marker sequences (optional)  → marker_interval
   → liftover (optional, minimap2)    → lift_interval
   → score_confidence → combine_intervals → final
   → write_outputs (artifacts) → summary.json
```

`run_job` accepts a `progress(percent, stage)` callback and a `cancel_event`; `JobManager` runs it on a single-worker executor and persists progress to `summary.json` so the UI can poll `/api/jobs/{id}`.

---

## 2. Running & verifying

```powershell
# Backend tests (fast, deterministic; builds the artificial sample first)
$env:PYTHONPATH="$PWD\backend"
py -3.13 -m unittest discover -s tests -v

# Frontend build (also type-checks)
cd frontend; npm install; npm run build

# Full app
powershell -ExecutionPolicy Bypass -File .\start-qtlift.ps1   # http://127.0.0.1:8765
```

**Definition of done for any change:**
- `unittest` passes (full suite) and you added a test for new backend behavior.
- `npm run build` passes (TypeScript is part of the build; no type errors).
- The artificial sample still completes end-to-end.
- New UI strings exist in **both** `en.json` and `ja.json`.
- No private data or identifiers introduced (grep before committing).

---

## 3. Key invariants & gotchas (learned the hard way)

- **BLAST target names:** a `-db` search puts the real accession in `stitle` but a `BL_ORD_ID` in `sseqid`; a `-subject` search is the opposite (`sseqid` has the name, `stitle` is `N/A`). `blast._parse` requests both and prefers `stitle`, falling back to `sseqid`. Don't "simplify" this away.
- **GFF dialects:** CDS `Parent` is usually the transcript id (`<gene>.1`, `<gene>.t1`, `<gene>-T1`). `genes_in_interval` strips those suffixes and also tries the CDS `Name` to re-attach CDS to the gene. If CDS don't attach, anchors become intron-containing full-gene spans and BLAST multi-hits explode. Preserve/extend this resolver for new dialects rather than removing it.
- **Anchor selection** is spacing-driven (`anchor_spacing_mb` = minimum separation, no upper cap) with a `min_anchors` floor, plus the two boundary genes and the peak gene always kept. See the "How it works" tab / `select_anchors`.
- **Confidence ceiling:** synteny alone is one evidence class, so **without markers or liftover, confidence maxes out at Medium**. `score_confidence` needs ≥2 agreeing classes for High.
- **Contig-name matching:** the target-contig filter compares BLAST hit names against the library FASTA/GFF names. They must share a namespace; a `makeblastdb` built without `-parse_seqids` won't expose accessions.
- **WSL paths & memory:** on Windows, tools often run via `wsl.exe`. `liftover.py` chunks contigs >500 MB to avoid WSL OOM and un-chunks the PAF coordinates on projection. Respect the chunk header format `>QTLIFT|name|offset`.
- **Frontend fetch races:** prefer data already returned by `/api/genomes/scan` (it carries up to 200 contigs + true `contig_count`) over per-keystroke fetches; guard any async fetch against stale responses.

---

## 4. Common extension points

- **Add a preset:** add an entry to `Params.preset` (spacing + min_anchors + thresholds) and to the segmented control + `en/ja.json` in `main.tsx`. Document it in the How-it-works preset table.
- **Add a mapping provider / backend:** extend `providers.PROVIDERS` + `provider_status` + `validate_provider`, branch in `pipeline.map_query`, and add the option to the Setup backend dropdown.
- **Add an evidence class:** produce an `Interval` with a distinct `evidence` label, feed it into `score_confidence` and `combine_intervals`, surface it in the Results evidence track and `reporting.write_outputs`.
- **New annotation dialect:** extend the CDS→gene resolver in `genomes.genes_in_interval` and add a unit test with the new attribute style.
- **Non-Windows launcher:** `run.py` is cross-platform (uvicorn on 127.0.0.1:8765); add a shell script mirroring `start-qtlift.ps1` (venv, npm build, launch). Keep the Windows path first-class.

---

## 5. Style

- Match the surrounding code. Backend favors compact, explicit functions; the frontend is a single `main.tsx` with small components and a `t` localization object (`t.high === "高"` detects Japanese for prose blocks).
- Keep changes minimal and well-scoped. Add tests. Update `SPEC.md` if you change the product contract.
- Every user-facing string is bilingual (EN + JA).

Happy hacking — and remember rule #0: **surface uncertainty, never fake a result.**
