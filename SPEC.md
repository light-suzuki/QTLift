# QTLift Product Specification

QTLift is a Windows-first local FastAPI + React/Vite application that infers the corresponding canonical reference A interval for a QTL/candidate interval reported on source reference B.

## Functional requirements

- Genome root selection: each child directory is a reference containing `.fa`, `.fasta`, `.fna`, `.gff3`, or `.gff`, optionally gzipped.
- Library view: reference status, sequence names and lengths, gene counts.
- Mapping input: A/B references, optional target-A chromosome, B contig/start/end, optional peak/name, and a dynamic list of markers. Markers are DNA sequences only (FASTA or raw bases); both strands are searched.
- Presets: Auto, Fast, Standard, Precise, Manual, and All genes. Anchor selection is spacing-driven with a minimum-count floor and no upper cap; Manual exposes minimum spacing and minimum anchors; All genes maps every annotated gene overlapping the source interval.
- Analysis: disk-backed FASTA access, streaming GFF parsing, B genes/CDS, uncapped anchor selection, marker/anchor BLASTn mapping, anchor order/orientation evaluation, evidence reconciliation, confidence with reasons and warnings. Whole-genome minimap2 liftover is disabled by default and never fabricated. When enabled, liftover projects interval endpoints by linear interpolation across the chosen alignment (it does not consult the CIGAR), so its coordinates are **approximate** supporting evidence for broad intervals, not marker-grade breakpoints.
- GFF3 parsing follows the real `gene → transcript → CDS` hierarchy and picks one representative transcript per gene; each anchor records that `transcript_id` and its `sequence_source` (`cds`, or `gene` when no CDS model is available).
- States: forward, reverse/inversion-like, partial, split, ambiguous, failed.
- Evidence: liftover, marker, synteny/anchor intervals. The result also records an `orientation` audit — the reconciled `strand`, the `orientation_evidence` and `uninformative_orientation_evidence` classes, and whether they `conflict` — so a flattened or opposing orientation is never silently hidden.
- Confidence: High, Medium, Low, Manual check.
- Required job artifacts: `report.html`, `summary.tsv`, `summary.json`, `marker_hits.tsv`, `anchor_genes.tsv`, `anchors.tsv` (per-anchor transcript/sequence-source provenance), optional `liftover_interval.bed`, `warnings.txt`, `params.json`, and `logs/qtlift.log`.
- UI pages: Genome Library, Mapping Setup, Job Monitor, Results, Settings; Japanese and English locale files.
- Results include all evidence, tables, warnings, parameters, downloads, interval summary and anchor-order plot.
- External tools are configurable and independently optional. Their absence never prevents app startup.

## Scientific contract

Coordinates entered by users are 1-based inclusive. BED output is 0-based half-open. Every result records its evidence source, orientation provenance, and parameters. QTLift produces practical inferred intervals, not exact breakpoint claims. It must lower confidence or require manual review when evidence is sparse, ambiguous, split, or inconsistent — including when independent evidence overlaps but disagrees on orientation.

## Completion gates

The app launches locally, the artificial sample completes, required artifacts exist, backend tests pass, frontend builds, and the browser workflow is verified.
