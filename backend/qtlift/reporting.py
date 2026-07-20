from __future__ import annotations

import csv
import html
import json
from pathlib import Path


def _tsv(path: Path, rows: list[dict], fields: list[str]):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)


def write_outputs(job_dir: str | Path, summary: dict, marker_hits: list[dict], anchor_hits: list[dict], params: dict, log_text: str = "") -> list[str]:
    out = Path(job_dir); out.mkdir(parents=True, exist_ok=True); (out/"logs").mkdir(exist_ok=True)
    (out/"summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    flat = {"job_id": summary["job_id"], "name": summary.get("name", ""), "status": summary.get("status", "completed"),
            "confidence": summary["confidence"], "source": summary["source_label"], "target": summary.get("final_label", ""),
            "duration_sec": summary.get("duration_sec", ""), "warning_count": len(summary.get("warnings", []))}
    _tsv(out/"summary.tsv", [flat], list(flat))
    hit_fields = ["query_id","contig","start","end","strand","identity","coverage","hit_count","method","source_start","query_start","query_end","query_length"]
    _tsv(out/"marker_hits.tsv", marker_hits, hit_fields); _tsv(out/"anchor_genes.tsv", anchor_hits, hit_fields)
    (out/"warnings.txt").write_text("\n".join(summary["warnings"]) + "\n", encoding="utf-8")
    (out/"params.json").write_text(json.dumps(params, indent=2, ensure_ascii=False), encoding="utf-8")
    liftover = summary.get("evidence", {}).get("liftover")
    if liftover:
        (out/"liftover_interval.bed").write_text(
            f"{liftover['contig']}\t{liftover['start'] - 1}\t{liftover['end']}\t{summary['job_id']}\t0\t{liftover['strand']}\n",
            encoding="utf-8")
    (out/"logs"/"qtlift.log").write_text(log_text, encoding="utf-8")
    report = f'''<!doctype html><html lang="ja"><meta charset="utf-8"><title>QTLift {html.escape(summary.get('name') or summary['job_id'])}</title>
<style>body{{font:14px system-ui;max-width:1100px;margin:32px auto;color:#14253d}}table{{border-collapse:collapse;width:100%}}th,td{{border:1px solid #ccd6e2;padding:6px;text-align:left}}.confidence{{font-size:24px;color:#087f79}}.warn{{background:#fff5db;padding:8px}}</style>
<h1>{html.escape(summary.get('name') or 'QTLift解析レポート')}</h1><p><code>{html.escape(summary['job_id'])}</code></p><p class="confidence">{html.escape(summary['confidence'])}</p><p><b>ソース / Source:</b> {html.escape(summary['source_label'])}</p><p><b>最終ターゲット / Final target:</b> {html.escape(summary.get('final_label','Unavailable'))}</p>
<h2>根拠 / Evidence</h2><pre>{html.escape(json.dumps(summary['evidence'],indent=2,ensure_ascii=False))}</pre><h2>信頼度の理由 / Confidence rationale</h2><ul>{''.join('<li>'+html.escape(x)+'</li>' for x in summary['reasons'])}</ul>
<h2>警告 / Warnings</h2><div class="warn">{'<br>'.join(html.escape(x) for x in summary['warnings']) or 'なし / None'}</div><h2>パラメータ / Parameters</h2><pre>{html.escape(json.dumps(params,indent=2))}</pre></html>'''
    (out/"report.html").write_text(report, encoding="utf-8")
    return sorted(str(p.relative_to(out)).replace("\\", "/") for p in out.rglob("*") if p.is_file())
