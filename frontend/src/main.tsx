import React, { useEffect, useMemo, useState } from "react";
import { createRoot } from "react-dom/client";
import {
  BookOpen,
  FlaskConical,
  Gauge,
  Settings as SettingsIcon,
  Library,
  Download,
  Play,
  RefreshCw,
  TriangleAlert,
  HelpCircle,
  Trash2,
  Clock,
  ArrowRight,
} from "lucide-react";
import en from "./locales/en.json";
import ja from "./locales/ja.json";
import "./styles.css";
import "./manual.css";
type Genome = {
  name: string;
  fasta_status: string;
  gff_status: string;
  contigs: { name: string; length: number }[];
  contig_count?: number;
  assembly_size?: number;
  gene_count: number;
};
type Job = any;
const api = async (path: string, options?: RequestInit) => {
  const r = await fetch(path, options);
  const x = await r.json();
  if (!r.ok) throw Error(x.detail || "Request failed");
  return x;
};
function App() {
  const [lang, setLang] = useState<"ja" | "en">("ja"),
    [page, setPage] = useState("library"),
    [root, setRoot] = useState(localStorage.qtliftRoot || ""),
    [genomes, setGenomes] = useState<Genome[]>([]),
    [jobs, setJobs] = useState<Job[]>([]),
    [active, setActive] = useState<Job | null>(null),
    [health, setHealth] = useState<any>(null),
    [error, setError] = useState("");
  const t = lang === "ja" ? ja : en;
  const scanAt = async (selectedRoot: string) => {
    try {
      setError("");
      const x = await api("/api/genomes/scan", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ genome_root: selectedRoot }),
      });
      setGenomes(x.genomes);
      setRoot(selectedRoot);
      localStorage.qtliftRoot = selectedRoot;
    } catch (e: any) {
      setGenomes([]);
      setError(e.message);
    }
  };
  const scan = () => scanAt(root);
  const autoDetect = () => scanAt("@this-pc");
  const refresh = async () => {
    const x = await api("/api/jobs");
    setJobs(x.jobs);
    setActive((previous: any) => previous ? (x.jobs.find((j: any) => j.job_id === previous.job_id) || previous) : (x.jobs[0] || null));
  };
  useEffect(() => {
    Promise.all([api("/api/health"), api("/api/capabilities")]).then(([x, capabilities]) => {
      setHealth({...x, ...capabilities});
      const defaultRoot = x.local_library_available ? "@this-pc" : x.sample_root;
      const storedRoot = localStorage.qtliftRoot || defaultRoot;
      const selected = localStorage.qtliftLibraryV2 && (storedRoot !== "@this-pc" || x.local_library_available) ? storedRoot : defaultRoot;
      localStorage.qtliftLibraryV2 = "1";
      scanAt(selected);
    });
    refresh();
  }, []);
  useEffect(() => {
    const timer = window.setInterval(refresh, 2000);
    return () => window.clearInterval(timer);
  }, []);
  const del = async (id: string) => {
    await fetch(`/api/jobs/${id}`, { method: "DELETE" });
    if (active?.job_id === id) setActive(null);
    const x = await api("/api/jobs");
    setJobs(x.jobs);
  };
  const cancel = async (id: string) => {
    await api(`/api/jobs/${id}/cancel`, { method: "POST" });
    await refresh();
  };
  const nav = [
    ["library", Library, t.library],
    ["setup", FlaskConical, t.setup],
    ["monitor", Gauge, t.monitor],
    ["results", BookOpen, t.results],
    ["howto", HelpCircle, t.howto],
    ["settings", SettingsIcon, t.settings],
  ] as const;
  return (
    <div className="app">
      <aside>
        <div className="brand">
          <span>QT</span>Lift
        </div>
        <nav>
          {nav.map(([id, I, label]) => (
            <button
              className={page === id ? "active" : ""}
              onClick={() => setPage(id)}
            >
              <I size={18} />
              {label}
            </button>
          ))}
        </nav>
        <div className="asideFoot">
          <span className="pulse" />
          Local · 127.0.0.1
        </div>
      </aside>
      <main>
        <header>
          <div>
            <b>{nav.find((x) => x[0] === page)?.[2]}</b>
            <small>{t.subtitle}</small>
          </div>
          <div className="headerTools">
            <button onClick={() => setLang(lang === "ja" ? "en" : "ja")}>
              {lang === "ja" ? "EN" : "日本語"}
            </button>
            <span className="version">v1.0</span>
          </div>
        </header>
        <section className="content">
          {error && (
            <div className="error">
              <TriangleAlert size={18} />
              {error}
            </div>
          )}
          {page === "library" && (
            <LibraryPage {...{ root, setRoot, scan, autoDetect, genomes, health, t }} />
          )}
          {page === "setup" && (
            <SetupPage
              genomes={genomes}
              root={root}
              health={health}
              onDone={(j: any) => {
                setActive(j);
                refresh();
                setPage("monitor");
              }}
              t={t}
            />
          )}{" "}
          {page === "monitor" && (
            <Monitor
              jobs={jobs}
              refresh={refresh}
              del={del}
              cancel={cancel}
              setActive={(j: any) => {
                setActive(j);
                setPage("results");
              }}
              t={t}
            />
          )}{" "}
          {page === "results" && <Results job={active} t={t} />}{" "}
          {page === "howto" && <HowItWorks t={t} />}{" "}
          {page === "settings" && <Settings health={health} t={t} />}
        </section>
      </main>
    </div>
  );
}
function LibraryPage({ root, setRoot, scan, autoDetect, genomes, health, t }: any) {
  return (
    <>
      <div className="titleRow">
        <div>
          <h1>{t.library}</h1>
          <p>{t.libraryDesc}</p>
        </div>
        <div className="actions">
          <button onClick={scan}><RefreshCw size={16} />{t.scan}</button>
          {health?.local_library_available && <button className="primary" onClick={autoDetect}><Library size={16} />{t.autoDetect}</button>}
        </div>
      </div>
      <div className="path">
        <label>{t.root}</label>
        <input value={root} onChange={(e) => setRoot(e.target.value)} />
      </div>
      <div className="tableWrap">
        <table>
          <thead>
            <tr>
              <th>{t.reference}</th>
              <th>FASTA</th>
              <th>GFF</th>
              <th>{t.contigs}</th>
              <th>{t.genes}</th>
              <th>{t.assemblySize}</th>
            </tr>
          </thead>
          <tbody>
            {genomes.map((g: Genome) => (
              <tr>
                <td>
                  <b>{g.name}</b>
                </td>
                <td>
                  <Status ok={g.fasta_status === "ready"} t={t} />
                </td>
                <td>
                  <Status ok={g.gff_status === "ready"} t={t} />
                </td>
                <td>
                  {g.contig_count ?? g.contigs.length}
                  <small>
                    {g.contigs
                      .map((x) => x.name)
                      .slice(0, 2)
                      .join(", ")}
                  </small>
                </td>
                <td>{g.gene_count}</td>
                <td>
                  {(g.assembly_size ?? g.contigs.reduce((a, x) => a + x.length, 0)).toLocaleString()}{" "}
                  bp
                </td>
              </tr>
            ))}
          </tbody>
        </table>
        {!genomes.length && (
          <div className="empty">{t.selectRoot}</div>
        )}
      </div>
    </>
  );
}
const Status = ({ ok, t }: { ok: boolean; t: any }) => (
  <span className={"status " + (ok ? "ok" : "bad")}>
    <i />
    {ok ? t.ready : t.missing}
  </span>
);
const parseCoord = (s: string): number | undefined => {
  const digits = String(s).replace(/[^0-9]/g, "");
  if (!digits) return undefined;
  const n = parseInt(digits, 10);
  return Number.isFinite(n) ? n : undefined;
};
function SetupPage({ genomes, root, health, onDone, t }: any) {
  const minimapOk = !!health?.tools?.minimap2?.available;
  const [a, setA] = useState(genomes[0]?.name || "RefA"),
    [b, setB] = useState(genomes[1]?.name || "RefB"),
    [contig, setContig] = useState("Chr1"),
    [targetContig, setTargetContig] = useState(""),
    [start, setStart] = useState("100"),
    [end, setEnd] = useState("850"),
    [peak, setPeak] = useState("450"),
    [jobName, setJobName] = useState(""),
    [preset, setPreset] = useState("Auto"),
    [backend, setBackend] = useState("auto"),
    [liftover, setLiftover] = useState(false),
    [manual, setManual] = useState({ anchor_spacing_mb: 0.25, min_anchors: 6 }),
    [busy, setBusy] = useState(false),
    [err, setErr] = useState(""),
    [markers, setMarkers] = useState<{ key: string; label: string; value: string; fixed: boolean }[]>([
      { key: "left_flanking", label: t.leftFlanking, value: "", fixed: true },
      { key: "peak", label: t.peakMarker, value: "", fixed: true },
      { key: "right_flanking", label: t.rightFlanking, value: "", fixed: true },
    ]);
  const sourceGenome = genomes.find((g: Genome) => g.name === b);
  const targetGenome = genomes.find((g: Genome) => g.name === a);
  const contigLength = sourceGenome?.contigs?.find((x: any) => x.name === contig)?.length;
  useEffect(() => {
    setTargetContig("");
  }, [a]);
  const updateMarker = (i: number, patch: any) =>
    setMarkers(markers.map((m, j) => (j === i ? { ...m, ...patch } : m)));
  const addMarker = () =>
    setMarkers([...markers, { key: `marker_${markers.length + 1}`, label: "", value: "", fixed: false }]);
  const removeMarker = (i: number) => setMarkers(markers.filter((_, j) => j !== i));
  useEffect(() => {
    if (genomes.length && !genomes.some((g: Genome) => g.name === a)) setA(genomes[0].name);
    if (genomes.length && !genomes.some((g: Genome) => g.name === b)) setB((genomes[1] || genomes[0]).name);
  }, [genomes]);
  useEffect(() => {
    const first = sourceGenome?.contigs?.[0]?.name;
    if (first && !sourceGenome.contigs.some((x: any) => x.name === contig)) setContig(first);
  }, [b, genomes]);
  const run = async () => {
    setErr("");
    const startN = parseCoord(start), endN = parseCoord(end), peakN = parseCoord(peak);
    if (startN === undefined || endN === undefined) return setErr(t.coordRequired);
    if (endN < startN) return setErr(t.coordOrder);
    if (peakN !== undefined && (peakN < startN || peakN > endN)) return setErr(t.peakRange);
    const markerDict: Record<string, string> = {};
    markers.forEach((m, i) => {
      const key = (m.fixed ? m.key : m.label.trim() || `marker_${i + 1}`);
      if (m.value.trim()) markerDict[key] = m.value;
    });
    setBusy(true);
    try {
      onDone(
        await api("/api/jobs", {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({
            genome_root: root,
            target_ref: a,
            source_ref: b,
            contig,
            target_contig: targetContig,
            start: startN,
            end: endN,
            peak: peakN ?? null,
            name: jobName.trim() || `${b} ${contig}:${startN}-${endN}`,
            preset,
            params: preset === "Manual" ? manual : {},
            all_genes: preset === "All genes",
            markers: markerDict,
            mapping_backend: backend,
            provider_options: liftover && minimapOk ? { enable_liftover: true } : {},
          }),
        }),
      );
    } catch (e: any) {
      setErr(e.message);
    } finally {
      setBusy(false);
    }
  };
  return (
    <>
      <div className="titleRow">
        <div>
          <h1>{t.setup}</h1>
          <p>{t.setupDesc}</p>
        </div>
      </div>
      <section className="card">
        <h2 className="cardTitle"><span className="step">1</span>{t.stepReferences}</h2>
        <div className="cardGrid">
          <Field label={t.target}>
            <select value={a} onChange={(e) => setA(e.target.value)}>
              {genomes.map((g: Genome) => (<option>{g.name}</option>))}
            </select>
          </Field>
          <Field label={t.source}>
            <select value={b} onChange={(e) => setB(e.target.value)}>
              {genomes.map((g: Genome) => (<option>{g.name}</option>))}
            </select>
          </Field>
          <Field label={t.targetContig}>
            <ContigPicker root={root} genome={targetGenome} value={targetContig} onChange={setTargetContig} allowAll t={t} />
          </Field>
        </div>
      </section>

      <section className="card">
        <h2 className="cardTitle"><span className="step">2</span>{t.stepInterval}</h2>
        <div className="cardGrid">
          <Field label={t.jobName}>
            <input value={jobName} maxLength={200} placeholder={t.jobNamePlaceholder} onChange={(e) => setJobName(e.target.value)} />
          </Field>
          <Field label={t.contig}>
            <ContigPicker root={root} genome={sourceGenome} value={contig} onChange={setContig} t={t} />
          </Field>
          <Field label={t.start}>
            <input inputMode="numeric" value={start} placeholder="1" onChange={(e) => setStart(e.target.value)} />
          </Field>
          <Field label={t.end}>
            <input inputMode="numeric" value={end} placeholder={contigLength ? contigLength.toLocaleString() : ""} onChange={(e) => setEnd(e.target.value)} />
          </Field>
          <Field label={`${t.peak}（${t.optional}）`}>
            <input inputMode="numeric" value={peak} placeholder={t.optional} onChange={(e) => setPeak(e.target.value)} />
          </Field>
        </div>
        <p className="cardHint">{t.stepIntervalHint}</p>
      </section>

      <section className="card">
        <h2 className="cardTitle"><span className="step">3</span>{t.stepAnalysis}</h2>
        <div className="cardGrid">
          <Field label={t.backend}>
            <select value={backend} onChange={(e) => setBackend(e.target.value)}>
              <option value="auto">{t.auto}</option>
              <option value="windows">{t.windowsBlast}</option>
              <option value="wsl">{t.wslBlast}</option>
              <option value="exact">{t.exactFallback}</option>
            </select>
          </Field>
        </div>
        <Field label={t.preset}>
          <div className="segments">
            {[["Auto",t.autoScale],["Fast",t.fast],["Standard",t.standard],["Precise",t.precise],["Manual",t.manual],["All genes",t.allGenes]].map(([x,label]) => (
              <button className={preset === x ? "selected" : ""} onClick={() => setPreset(x)}>{label}</button>
            ))}
          </div>
        </Field>
        {preset === "Manual" && <div className="manualParams">
          <Field label={t.anchorSpacing}><input type="number" min="0.001" step="0.01" value={manual.anchor_spacing_mb} onChange={e=>setManual({...manual,anchor_spacing_mb:+e.target.value})}/></Field>
          <Field label={t.minimumAnchors}><input type="number" min="1" value={manual.min_anchors} onChange={e=>setManual({...manual,min_anchors:+e.target.value})}/></Field>
        </div>}
        {preset === "Manual" && <div className="notice">{t.spacingHelp}</div>}
        {preset === "All genes" && <div className="notice">{t.allGenesHelp}</div>}
        <label className={"toggle" + (minimapOk ? "" : " disabled")}>
          <input type="checkbox" checked={liftover && minimapOk} disabled={!minimapOk} onChange={(e) => setLiftover(e.target.checked)} />
          <span>{t.enableLiftover}</span>
        </label>
        <p className="cardHint subtle">{minimapOk ? t.liftoverHint : t.liftoverUnavailable}</p>
        <p className="cardHint subtle">{t.presetHintShort}</p>
      </section>

      <section className="card">
        <h2 className="cardTitle row">
          <span><span className="step">4</span>{t.markers} <small>{t.markersOptional}</small></span>
          <button className="ghost" onClick={addMarker}>+ {t.addMarker}</button>
        </h2>
        <div className="markerRows">
          {markers.map((m, i) => (
            <div className="markerRow">
              {m.fixed ? (<span className="markerLabel">{m.label}</span>) : (
                <input className="markerName" value={m.label} placeholder={t.markerName} onChange={(e) => updateMarker(i, { label: e.target.value })} />
              )}
              <textarea value={m.value} placeholder={t.markerHint} onChange={(e) => updateMarker(i, { value: e.target.value })} />
              {!m.fixed && (<button className="rowRemove" onClick={() => removeMarker(i)}>{t.removeRow}</button>)}
            </div>
          ))}
        </div>
      </section>
      {err && (
        <div className="error">
          <TriangleAlert size={18} />
          {err}
        </div>
      )}
      <button
        className="primary run"
        disabled={busy || !genomes.length}
        onClick={run}
      >
        <Play size={16} />
        {busy ? t.running : t.run}
      </button>
    </>
  );
}
function Field({ label, children }: any) {
  return (
    <label className="field">
      <span>{label}</span>
      {children}
    </label>
  );
}
function ContigPicker({ root, genome, value, onChange, allowAll, t }: any) {
  // The scanned genome already carries up to 200 contigs plus the true contig_count, so the
  // dropdown is populated directly (no fetch, no race). The /api/genomes/contigs search is
  // used only to reach contigs beyond the first 200 on scaffold-heavy assemblies.
  const baseContigs = genome?.contigs || [];
  const total = genome?.contig_count ?? baseContigs.length;
  const overflow = total > baseContigs.length;
  const [search, setSearch] = useState("");
  const [found, setFound] = useState<any[]>([]);
  useEffect(() => {
    if (!overflow || !search.trim() || !root || !genome?.name) { setFound([]); return; }
    let alive = true;
    const timer = window.setTimeout(() => {
      api("/api/genomes/contigs", {method:"POST",headers:{"Content-Type":"application/json"},body:JSON.stringify({genome_root:root,reference:genome.name,query:search,limit:300})})
        .then((x:any) => { if (alive) setFound(x.contigs || []); })
        .catch(() => { if (alive) setFound([]); });
    }, 200);
    return () => { alive = false; window.clearTimeout(timer); };
  }, [root, genome?.name, search, overflow]);
  const searching = overflow && !!search.trim();
  const shown = searching ? found : baseContigs;
  const known = value && shown.some((x:any) => x.name === value);
  return <div className="contigPicker">
    <select value={value} onChange={e => onChange(e.target.value)}>
      {allowAll && <option value="">{t.allContigs}</option>}
      {value && !known && <option value={value}>{value}</option>}
      {shown.map((x:any) => <option key={x.name} value={x.name}>{x.name} ({x.length.toLocaleString()} bp)</option>)}
    </select>
    {overflow && <input className="contigSearch" value={search} placeholder={t.contigSearch} onChange={e => setSearch(e.target.value)} />}
    <small>{overflow ? `${searching ? found.length : baseContigs.length} / ${total.toLocaleString()} ${t.contigMatches}` : `${total.toLocaleString()} ${t.contigMatches}`}</small>
  </div>;
}
function fmtDate(iso: string | undefined) {
  if (!iso) return "—";
  return iso.replace("T", " ").slice(0, 16);
}
function fmtBp(n: number) {
  if (n >= 1e6) return (n / 1e6).toFixed(n >= 1e7 ? 0 : 1) + " Mb";
  if (n >= 1e3) return Math.round(n / 1e3) + " kb";
  return n + " bp";
}
function fmtDur(s: number | undefined) {
  if (s == null) return "—";
  const m = Math.floor(s / 60), sec = Math.round(s % 60);
  return m ? `${m}m ${sec}s` : `${sec}s`;
}
function Monitor({ jobs, refresh, del, cancel, setActive, t }: any) {
  return (
    <>
      <div className="titleRow">
        <div>
          <h1>{t.monitor}</h1>
          <p>{t.monitorDesc}</p>
        </div>
        <button onClick={refresh}>
          <RefreshCw size={16} />
          {t.refresh}
        </button>
      </div>
      <div className="tableWrap">
        <table>
          <thead>
            <tr>
              <th>{t.job}</th>
              <th>{t.runAt}</th>
              <th>{t.sourceLabel}</th>
              <th>{t.target}</th>
              <th>{t.confidence}</th>
              <th>{t.status}</th>
              <th></th>
            </tr>
          </thead>
          <tbody>
            {jobs.map((j: any) => (
              <tr onClick={() => setActive(j)} className="jobRow">
                <td>
                  <b>{j.name}</b>
                  <small>{j.job_id}</small>
                </td>
                <td className="runAt"><Clock size={13} />{fmtDate(j.created_at)}</td>
                <td>{j.source_label}</td>
                <td>{j.final_label}</td>
                <td>
                  <Confidence value={j.confidence} t={t} />
                </td>
                <td><div className="progressCell"><span>{t["status_"+j.status] || j.status}</span><progress max="100" value={j.progress || 0}/><small>{j.stage}</small></div></td>
                <td>
                  {(j.status === "queued" || j.status === "running" || j.status === "cancelling") && <button className="rowRemove" title={t.cancelJob} onClick={(e)=>{e.stopPropagation();cancel(j.job_id)}}>{t.cancel}</button>}
                  <button
                    className="rowRemove danger"
                    title={t.deleteJob}
                    onClick={(e) => {
                      e.stopPropagation();
                      if (confirm(t.deleteConfirm)) del(j.job_id);
                    }}
                  >
                    <Trash2 size={14} />
                  </button>
                </td>
              </tr>
            ))}
          </tbody>
        </table>
        {!jobs.length && <div className="empty">{t.noJobs}</div>}
      </div>
    </>
  );
}
function Confidence({ value, t }: any) {
  const label = value === "High" ? t.high : value === "Medium" ? t.medium : value === "Low" ? t.low : t.manualCheck;
  return (
    <span
      className={"confidence " + String(value).toLowerCase().replace(" ", "-")}
    >
      {label}
    </span>
  );
}
function localizeMessage(message: string | undefined, t: any) {
  if (!message || t.high !== "高") return message || "";
  const fixed: Record<string, string> = {
    "At least two evidence classes agree with four or more unique collinear anchors.": "2種類以上の根拠が一致し、4個以上の一意なアンカーが共線性を支持しています。",
    "A coherent anchor interval is supported, but independent evidence is limited.": "整合したアンカー区間を支持しますが、独立した根拠が限定的です。",
    "Only sparse or partially consistent evidence is available.": "利用できる根拠が少ないか、一部のみ整合しています。",
    "No defensible target interval was produced.": "根拠をもって提示できるターゲット区間を作成できませんでした。",
    "Fewer than two uniquely mapped markers; marker interval unavailable.": "一意にマッピングできたマーカーが2個未満のため、マーカー区間を算出できません。",
    "Liftover skipped: minimap2 and/or paftools.js is unavailable.": "minimap2またはpaftools.jsを利用できないため、リフトオーバーを省略しました。",
    "Liftover not executed for this interval-only job; build or reuse a pairwise whole-genome alignment cache to enable BED liftover.": "この区間解析ではリフトオーバーを実行していません。BEDリフトオーバーには参照ペアの全ゲノムアラインメントキャッシュが必要です。",
    "Liftover skipped: whole-genome minimap2 liftover is disabled.": "全ゲノムminimap2リフトオーバーは無効のため省略しました。",
    "Independent evidence intervals disagree and must be reviewed separately.": "独立した根拠区間が一致しないため、候補を分けて手動確認してください。"
  };
  if (fixed[message]) return fixed[message];
  let match = message.match(/^Too few genes in source interval \((\d+)\)\.$/);
  if (match) return `ソース区間内の遺伝子が少なすぎます（${match[1]}個）。`;
  match = message.match(/^Anchor (.+) has excessive multi-hits \((\d+)\)\.$/);
  if (match) return `アンカー ${match[1]} の複数ヒットが多すぎます（${match[2]}件）。`;
  if (message.includes("exact fallback used")) return `外部マッピングに失敗したため完全一致フォールバックを使用しました: ${message}`;
  return message;
}
const SYNTENY_LABEL: Record<string, string> = { forward: "forward", reverse: "reverse", partial: "partial", split: "split", ambiguous: "ambiguous", failed: "failed" };
function Results({ job, t }: any) {
  if (!job) return <div className="empty big">{t.noJobs}</div>;
  if (job.status && job.status !== "completed") return <div className="jobPending"><h1>{job.name}</h1><p>{t["status_"+job.status] || job.status}</p><progress max="100" value={job.progress || 0}/><b>{job.stage}</b>{job.error && <div className="error">{job.error}</div>}</div>;
  const ev = job.evidence || {};
  const anchors = job.anchor_hits || [];
  const uniqA = anchors.filter((h: any) => h.hit_count === 1);
  const uniqM = (job.marker_hits || []).filter((h: any) => h.hit_count === 1);
  const fwd = uniqA.filter((h: any) => h.strand !== "-").length;
  const rev = uniqA.filter((h: any) => h.strand === "-").length;
  const syntenyKey = "synteny_" + job.synteny_state;
  return (
    <>
      <div className="resHead">
        <span>{t.job}: <b>{job.name}</b> <small>{job.job_id}</small></span>
        <span><Clock size={13} />{fmtDate(job.created_at)}</span>
        <span>{t.duration}: <b>{fmtDur(job.duration_sec)}</b></span>
        <span>{t.backend}: <b>{job.effective_backend || job.mapping_backend}</b></span>
      </div>
      <div className="resultHero">
        <div className="heroFlow">
          <div className="heroCol">
            <small>{t.sourceB}</small>
            <b>{job.source_label}</b>
          </div>
          <ArrowRight className="heroArrow" size={26} />
          <div className="heroCol">
            <small>{t.finalA}</small>
            <b className="heroFinal">{job.final_label}</b>
          </div>
        </div>
        <div className="heroConf">
          <Confidence value={job.confidence} t={t} />
          <p>{localizeMessage(job.reasons?.[0], t)}</p>
        </div>
      </div>
      {job.candidates?.length > 1 && (
        <div className="candList">
          <span className="candTitle">{t.candidates}</span>
          {job.candidates.map((c: any, i: number) => (
            <span className="candItem"><b>{i + 1}</b> {c.contig}:{c.start.toLocaleString()}–{c.end.toLocaleString()} <small>{c.evidence}</small></span>
          ))}
        </div>
      )}
      <div className="metaRow">
        <span>{t.orientation}: <b>{t[syntenyKey] || SYNTENY_LABEL[job.synteny_state] || job.synteny_state}</b></span>
        <span>{t.uniqueAnchors}: <b>{uniqA.length} / {anchors.length}</b> <span className="ori">{fwd}→ · {rev}←</span></span>
        {uniqM.length > 0 && <span>{t.markerEvidence}: <b>{uniqM.length}</b></span>}
        {job.target_contig && <span>{t.targetContig}: <b>{job.target_contig}</b></span>}
      </div>
      <EvidenceTrack job={job} anchors={uniqA} markers={uniqM} t={t} />
      <div className="plots">
        <DotPlot anchors={uniqA} t={t} />
        <TwoTrack job={job} anchors={uniqA} t={t} />
      </div>
      <div className="split">
        <div>
          <h2>{t.anchorEvidence}</h2>
          <HitTable rows={anchors} t={t} kind="anchor" />
        </div>
        <div>
          <h2>{t.markerHits}</h2>
          <HitTable rows={job.marker_hits || []} t={t} kind="marker" />
        </div>
      </div>
      <div className="split">
        <div>
          <h2>{t.warnings}</h2>
          <div className="warnings">
            {(job.warnings || []).map((x: string) => (
              <p><TriangleAlert size={15} />{localizeMessage(x, t)}</p>
            ))}
            {!job.warnings?.length && <p className="ok">{t.noWarnings}</p>}
          </div>
          <h2>{t.downloads}</h2>
          <div className="downloads">
            {job.files?.filter((x: string) => !x.startsWith("logs/")).map((f: string) => (
              <a href={`/api/jobs/${job.job_id}/files/${f}`}><Download size={14} />{f}</a>
            ))}
          </div>
        </div>
        <div>
          <h2>{t.parameters}</h2>
          <ParamSummary job={job} t={t} />
        </div>
      </div>
    </>
  );
}
function arrowPath(x: number, y: number, forward: boolean) {
  const w = 8;
  return forward
    ? `M ${x - w} ${y} L ${x + w - 3} ${y} M ${x + w - 4} ${y - 4} L ${x + w} ${y} L ${x + w - 4} ${y + 4}`
    : `M ${x + w} ${y} L ${x - w + 3} ${y} M ${x - w + 4} ${y - 4} L ${x - w} ${y} L ${x - w + 4} ${y + 4}`;
}
function EvidenceTrack({ job, anchors, markers, t }: any) {
  const final = job.final;
  const pts: number[] = [...anchors.map((h: any) => h.start), ...markers.map((h: any) => h.start)];
  if (final) pts.push(final.start, final.end);
  if (!pts.length) return null;
  let lo = Math.min(...pts), hi = Math.max(...pts);
  if (hi <= lo) hi = lo + 1;
  const pad = (hi - lo) * 0.06; lo -= pad; hi += pad;
  const X0 = 150, X1 = 930;
  const ax = (v: number) => X0 + ((v - lo) / (hi - lo)) * (X1 - X0);
  const ticks = [0, 1, 2, 3, 4].map((i) => lo + ((hi - lo) * i) / 4);
  const lift = job.evidence?.liftover;
  return (
    <div className="evCard">
      <div className="plotHead">
        <h2>{t.evidenceSummary}</h2>
        <div className="dotLegend">
          <span><i className="dot fwd" />{t.forward}</span>
          <span><i className="dot rev" />{t.reverse}</span>
          <span><i className="dot band" />{t.inferredInterval}</span>
        </div>
      </div>
      <svg viewBox="0 0 1000 190" className="evTrack" role="img">
        {final && <rect className="evBand" x={ax(final.start)} y={34} width={Math.max(2, ax(final.end) - ax(final.start))} height={140} />}
        <line className="evAxis" x1={X0} y1={44} x2={X1} y2={44} />
        {ticks.map((v) => (
          <g><line className="evAxis" x1={ax(v)} y1={41} x2={ax(v)} y2={47} /><text className="evTick" x={ax(v)} y={30}>{fmtBp(Math.round(v))}</text></g>
        ))}
        <text className="evLabel" x={14} y={80}>{t.liftover}</text>
        {lift ? <rect className="evLift" x={ax(lift.start)} y={70} width={Math.max(2, ax(lift.end) - ax(lift.start))} height={12} /> : <text className="evNA" x={X0} y={80}>{t.notRun}</text>}
        <text className="evLabel" x={14} y={118}>{t.markerEvidence}</text>
        {markers.map((h: any) => <line className="evMark" x1={ax(h.start)} y1={106} x2={ax(h.start)} y2={128} />)}
        <text className="evSupport" x={X1 + 12} y={118}>{markers.length}</text>
        <text className="evLabel" x={14} y={158}>{t.syntenyEvidence}</text>
        {anchors.map((h: any) => <path className={h.strand === "-" ? "evArrowRev" : "evArrowFwd"} d={arrowPath(ax(h.start), 156, h.strand !== "-")} />)}
        <text className="evSupport" x={X1 + 12} y={158}>{anchors.length}</text>
      </svg>
    </div>
  );
}
function DotPlot({ anchors, t }: any) {
  if (!anchors.length) return <div className="plot"><h2>{t.anchorPlot}</h2><div className="empty">{t.notAvailable}</div></div>;
  const xs = anchors.map((h: any) => h.source_start || 0), ys = anchors.map((h: any) => h.start);
  const minX = Math.min(...xs), maxX = Math.max(...xs), minY = Math.min(...ys), maxY = Math.max(...ys);
  const px = (v: number) => 48 + ((v - minX) / ((maxX - minX) || 1)) * 360;
  const py = (v: number) => 150 - ((v - minY) / ((maxY - minY) || 1)) * 132;
  return (
    <div className="plot">
      <div className="plotHead">
        <h2>{t.anchorPlot}</h2>
        <div className="dotLegend"><span><i className="dot fwd" />{t.forward}</span><span><i className="dot rev" />{t.reverse}</span></div>
      </div>
      <svg viewBox="0 0 430 180" role="img">
        <path className="axis" d="M44 12V152H424" />
        {anchors.map((h: any) => <circle className={h.strand === "-" ? "cRev" : "cFwd"} cx={px(h.source_start || 0)} cy={py(h.start)} r="4" />)}
        <text className="axisTxt" x="215" y="174">{t.sourceCoordinate}</text>
        <text className="axisTxt" x="10" y="82" transform="rotate(-90 10 82)">{t.targetCoordinate}</text>
      </svg>
    </div>
  );
}
function TwoTrack({ job, anchors, t }: any) {
  const src = job.source, fin = job.final;
  if (!src || !fin) return <div className="plot"><h2>{t.intervalPlot}</h2><div className="empty">{t.notAvailable}</div></div>;
  const X0 = 20, X1 = 410;
  const bx = (v: number) => X0 + ((v - src.start) / ((src.end - src.start) || 1)) * (X1 - X0);
  const axf = (v: number) => X0 + ((v - fin.start) / ((fin.end - fin.start) || 1)) * (X1 - X0);
  return (
    <div className="plot">
      <h2>{t.intervalPlot}</h2>
      <svg viewBox="0 0 430 150" role="img">
        <text className="trkLabel" x={X0} y={22}>{t.sourceB} · {src.contig}</text>
        <rect className="trkBg" x={X0} y={28} width={X1 - X0} height={12} />
        <rect className="trkFill" x={X0} y={28} width={X1 - X0} height={12} />
        {anchors.map((h: any) => {
          const b = bx(h.source_start || src.start), a = axf(h.start);
          if (b < X0 - 2 || b > X1 + 2 || a < X0 - 2 || a > X1 + 2) return null;
          return <line className={h.strand === "-" ? "linkRev" : "linkFwd"} x1={b} y1={40} x2={a} y2={104} />;
        })}
        <rect className="trkBg" x={X0} y={104} width={X1 - X0} height={12} />
        <rect className="trkFill" x={X0} y={104} width={X1 - X0} height={12} />
        <text className="trkLabel" x={X0} y={134}>{t.finalA} · {fin.contig}</text>
      </svg>
    </div>
  );
}
function HitTable({ rows, t, kind }: any) {
  if (!rows.length) return <div className="tableWrap"><div className="empty">{t.notAvailable}</div></div>;
  return (
    <div className="tableWrap compact">
      <table>
        <thead><tr><th>{kind === "marker" ? t.marker : t.anchor}</th><th>B → A</th><th>{t.strandCol}</th><th>{t.identityCoverage}</th><th>{t.hits}</th></tr></thead>
        <tbody>
          {rows.slice(0, 12).map((h: any) => (
            <tr>
              <td>{h.query_id}</td>
              <td className="mono">{h.source_start ? h.source_start.toLocaleString() : "—"} → {h.start.toLocaleString()}</td>
              <td><span className={h.strand === "-" ? "strand rev" : "strand fwd"}>{h.strand === "-" ? "←" : "→"}</span></td>
              <td>{h.identity}% / {h.coverage}%</td>
              <td className={h.hit_count > 1 ? "multi" : ""}>{h.hit_count}</td>
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}
function ParamSummary({ job, t }: any) {
  const p = job.params || {};
  const rows: [string, any][] = [
    [t.backend, job.effective_backend || job.mapping_backend],
    [t.anchorSpacing, p.anchor_spacing_mb != null ? p.anchor_spacing_mb + " Mb" : "—"],
    [t.minimumAnchors, p.min_anchors ?? "—"],
    [t.identityCoverage, `${p.min_identity ?? "—"}% / ${p.min_coverage ?? "—"}%`],
    [t.target, job.target_reference],
    [t.targetContig, job.target_contig || t.allContigs],
    [t.runAt, fmtDate(job.created_at)],
    [t.duration, fmtDur(job.duration_sec)],
  ];
  return <div className="paramTable">{rows.map(([k, v]) => <div><span>{k}</span><b>{v}</b></div>)}</div>;
}
function Settings({ health, t }: any) {
  const providerText: Record<string,string> = {auto:t.providerAuto,windows:t.providerWindows,wsl:t.providerWsl,ncbi:t.providerNcbi,exact:t.providerExact};
  return (
    <>
      <div className="titleRow">
        <div>
          <h1>{t.settings}</h1>
          <p>{t.settingsDesc}</p>
        </div>
      </div>
      <h2>{t.providers}</h2>
      <div className="toolList">
        {Object.entries(health?.providers || {}).map(([name, v]: any) => (
          <div>
            <code>{v.label}</code>
            <Status ok={health?.provider_status?.[name]?.available} t={t} />
            <span>{providerText[name] || v.scope}</span>
          </div>
        ))}
      </div>
      <h2>{t.tools}</h2>
      <div className="toolList">
        {Object.entries(health?.tools || {}).map(([name, v]: any) => (
          <div>
            <code>{name}</code>
            <Status ok={v.available} t={t} />
            <span>{v.path || t.notFound}</span>
          </div>
        ))}
      </div>
      <div className="notice">{t.liftover}: {health?.liftover_enabled ? t.enabled : t.disabled}</div>
      <div className="notice">{t.providerNotice}</div>
    </>
  );
}
function HowItWorks({ t }: any) {
  const ja = t.high === "高";
  return (
    <>
      <div className="titleRow">
        <div>
          <h1>{t.howto}</h1>
          <p>{ja ? "QTLiftが区間をどのように対応付けているかの仕組み。" : "How QTLift maps an interval between references."}</p>
        </div>
      </div>
      <div className="doc">
        <section>
          <h2>{ja ? "目的" : "Goal"}</h2>
          <p>{ja
            ? "ソース参照Bで報告されたQTL/候補区間が、基準参照Aのどこに当たるかを、複数の独立した根拠から推定します。塩基単位で完全なブレークポイントの主張ではなく、根拠と信頼度つきの実用的な推定区間を返します。"
            : "Given a QTL/candidate interval reported on source reference B, QTLift infers where it lies on canonical reference A using several independent lines of evidence. It reports a practical inferred interval with rationale and confidence — not a base-perfect breakpoint claim."}</p>
        </section>
        <section>
          <h2>{ja ? "全体の流れ" : "Pipeline"}</h2>
          <ol>
            <li>{ja ? "ソースB区間に重なる遺伝子をGFFから抽出する。" : "Extract genes overlapping the source-B interval from the GFF."}</li>
            <li>{ja ? "アンカー遺伝子を選び、その配列(CDS優先)を基準AへBLASTで対応付ける。" : "Select anchor genes and map their sequence (CDS preferred) onto reference A with BLAST."}</li>
            <li>{ja ? "マーカー配列(任意)も同様に対応付ける。" : "Map optional marker sequences the same way."}</li>
            <li>{ja ? "アンカーの並び順・向きからシンテニー(共線性)を評価する。" : "Evaluate synteny (collinearity) from anchor order and orientation."}</li>
            <li>{ja ? "マーカー・シンテニーの各根拠を統合し、信頼度を判定する。" : "Reconcile marker and synteny evidence and score confidence."}</li>
          </ol>
        </section>
        <section>
          <h2>{ja ? "アンカーの配置方法" : "How anchors are placed"}</h2>
          <p>{ja
            ? "アンカーは「最小アンカー間隔」を主軸に選びます。左端の遺伝子から順に走査し、直前に採用したアンカーから最小間隔以上離れた遺伝子だけを採用します。上限本数は設けないため、区間が広いほどアンカーは比例して増えます。"
            : "Anchors are chosen primarily by a minimum spacing. Walking from the left, a gene is kept only when it is at least the minimum separation from the previously kept anchor. There is no upper cap, so a wider interval yields proportionally more anchors."}</p>
          <ul>
            <li>{ja ? "区間の両端の遺伝子とピーク遺伝子は、向きの基準として常にアンカーに含めます。" : "The two interval-boundary genes and the peak gene are always included to anchor orientation."}</li>
            <li>{ja ? "区間が狭い/遺伝子が少なく本数が「最低アンカー本数」に満たない場合は、間隔を詰めて最低本数分の遺伝子を区間全体に均等配置し、結果をある程度担保します。" : "If a narrow or gene-sparse interval falls below the minimum anchor count, the spacing tightens to spread that many genes evenly, guaranteeing a usable result."}</li>
            <li>{ja ? "プリセットは最小間隔と最低本数の組み合わせです（下表）。" : "Presets are combinations of minimum spacing and minimum count (table below)."}</li>
          </ul>
          <div className="tableWrap">
            <table className="presetTable">
              <thead><tr><th>{ja ? "プリセット" : "Preset"}</th><th>{ja ? "最小アンカー間隔" : "Min spacing"}</th><th>{ja ? "最低本数" : "Min anchors"}</th><th>{ja ? "挙動" : "Behaviour"}</th></tr></thead>
              <tbody>
                {[
                  [ja ? "自動（既定）" : "Auto (default)", "0.5 Mb", "6", ja ? "区間幅に追従。遺伝子が50個以下なら全遺伝子を使用（最も正確）。" : "Follows interval width; maps all genes when ≤50 (most accurate)."],
                  [ja ? "高速" : "Fast", "0.5 Mb", "4", ja ? "本数少なめ・最速。" : "Fewer anchors, fastest."],
                  [ja ? "標準" : "Standard", "0.25 Mb", "6", ja ? "バランス型。" : "Balanced."],
                  [ja ? "精密" : "Precise", "0.1 Mb", "10", ja ? "密にアンカーを配置。" : "Dense anchoring."],
                  [ja ? "手動" : "Manual", "—", "—", ja ? "最小間隔と最低本数を自分で指定。" : "Set minimum spacing and minimum anchors yourself."],
                  [ja ? "区間内の全遺伝子" : "All genes", "—", "—", ja ? "区間に重なる全遺伝子をマップ。" : "Maps every gene overlapping the interval."],
                ].map((r: any) => (<tr><td><b>{r[0]}</b></td><td>{r[1]}</td><td>{r[2]}</td><td>{r[3]}</td></tr>))}
              </tbody>
            </table>
          </div>
          <p>{ja ? "いずれのプリセットも上限本数は設けません。区間が広いほどアンカーは比例して増え、狭い区間では最低本数を下回らないよう間隔を詰めて確保します。" : "No preset caps the count: a wider interval yields proportionally more anchors, and a narrow one tightens the spacing so it never drops below the minimum."}</p>
        </section>
        <section>
          <h2>{ja ? "シンテニー(向き)の判定" : "Synteny (orientation)"}</h2>
          <p>{ja
            ? "一意にマッピングできたアンカーを、ソース座標順に並べたときのターゲット座標の増減で向きを判定します。"
            : "Using uniquely-mapped anchors ordered by source coordinate, the direction is inferred from whether target coordinates increase or decrease."}</p>
          <ul>
            <li><b>forward</b> — {ja ? "同じ向き(共線)。" : "same orientation (collinear)."}</li>
            <li><b>reverse</b> — {ja ? "逆位様(逆向き共線)。" : "inversion-like (reverse-collinear)."}</li>
            <li><b>partial</b> — {ja ? "一部のみ共線。" : "only partially collinear."}</li>
            <li><b>split</b> — {ja ? "複数のターゲット染色体に分散。要手動確認。" : "spread across multiple target contigs; manual review."}</li>
            <li><b>failed</b> — {ja ? "一意なアンカーが不足。" : "too few unique anchors."}</li>
          </ul>
        </section>
        <section>
          <h2>{ja ? "根拠の統合と信頼度" : "Evidence & confidence"}</h2>
          <p>{ja
            ? "マーカー区間とシンテニー区間を統合し、独立した根拠がどれだけ一致するかで信頼度を決めます。高信頼には、2種類以上の根拠が同じ染色体上で重なり、4個以上の一意なアンカーが共線性を支持することが必要です。分裂・失敗・不整合は信頼度を下げるか手動確認を促します。"
            : "Marker and synteny intervals are reconciled; confidence depends on how well independent lines agree. High confidence requires at least two evidence classes overlapping on one contig with four or more unique collinear anchors. Split, failed, or inconsistent evidence lowers confidence or forces manual review."}</p>
          <div className="callout"><b>{ja ? "重要" : "Note"}:</b> {ja
            ? "マーカー（配列）を指定しない場合、独立した根拠はシンテニーの1つだけになるため、信頼度は最大でも Medium（中）に留まります。High（高）にするには左右外側マーカーやピークマーカーのDNA配列を追加してください。"
            : "Without marker sequences, synteny is the only evidence class, so confidence is capped at Medium. Add flanking or peak marker DNA sequences to reach High."}</div>
        </section>
        <section>
          <h2>{ja ? "対象染色体の絞り込み・座標規則" : "Target-chromosome focus & coordinates"}</h2>
          <ul>
            <li>{ja ? "ターゲットA側の対象染色体を選ぶと、BLASTヒットをその染色体に限定し、他染色体の偽ヒットによるノイズを除きます。" : "Selecting a target-A chromosome restricts BLAST hits to it, removing off-chromosome noise."}</li>
            <li>{ja ? "UI/APIの座標は1始まり両端含む。BEDは0始まり半開区間。" : "UI/API coordinates are 1-based inclusive; BED is 0-based half-open."}</li>
            <li>{ja ? "マーカーはDNA配列(FASTAまたは生塩基)のみ。向き・相補鎖は自動で両鎖を探索します。" : "Markers are DNA sequences only (FASTA or raw bases); both strands are searched automatically."}</li>
            <li>{ja ? "結果は推定研究区間であり、存在しないリフトオーバー/アライメント結果は決して捏造しません。" : "Results are inferred research intervals; missing liftover/alignment results are never fabricated."}</li>
          </ul>
        </section>
      </div>
    </>
  );
}
createRoot(document.getElementById("root")!).render(
  <React.StrictMode>
    <App />
  </React.StrictMode>,
);
