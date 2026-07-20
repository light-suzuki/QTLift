from __future__ import annotations
import hashlib, os, shutil, subprocess, tempfile, time
from pathlib import Path
from .models import Hit
FIELDS = "qseqid sseqid stitle sstart send pident qstart qend qlen"


def _find_blast_entry(metadata: str, contig: str) -> str | None:
    for line in metadata.splitlines():
        entry, _, title = line.partition("\t")
        if title.split()[0] == contig:
            return entry
    return None


def _wsl_contig_db(prefix: list[str], target_db: str, target_contig: str) -> str:
    """Create/reuse a BLAST DB containing only the selected target contig."""
    db_info = subprocess.run(prefix + ["blastdbcmd", "-db", target_db, "-info"], check=True, text=True,
                             capture_output=True, timeout=60).stdout
    key = hashlib.sha256(f"{target_db}\n{target_contig}\n{db_info}".encode()).hexdigest()[:20]
    user = subprocess.run(prefix + ["whoami"], check=True, text=True, capture_output=True).stdout.strip()
    cache_base = f"/home/{user}/.cache/qtlift/blast_contigs"
    root = f"{cache_base}/{key}"
    db = f"{root}/subject"
    exists = subprocess.run(prefix + ["test", "-f", f"{root}/READY"]).returncode == 0
    if exists:
        return db
    lock = f"{root}.lock"
    subprocess.run(prefix + ["mkdir", "-p", cache_base], check=True, timeout=30)
    acquired = False
    for _ in range(120):
        if subprocess.run(prefix + ["mkdir", lock], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode == 0:
            acquired = True
            break
        if subprocess.run(prefix + ["test", "-f", f"{root}/READY"]).returncode == 0:
            return db
        time.sleep(1)
    if not acquired:
        raise RuntimeError(f"Timed out waiting for target-contig BLAST cache: {target_contig}")
    temp_root = f"{root}.tmp-{os.getpid()}"
    metadata = subprocess.run(prefix + ["blastdbcmd", "-db", target_db, "-entry", "all", "-outfmt", "%i\t%t"],
                              check=True, text=True, capture_output=True, timeout=60).stdout
    entry = _find_blast_entry(metadata, target_contig)
    if not entry:
        subprocess.run(prefix + ["rmdir", lock], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        raise RuntimeError(f"Target contig is not present in the BLAST database: {target_contig}")
    try:
        subprocess.run(prefix + ["rm", "-rf", temp_root], check=True, timeout=30)
        subprocess.run(prefix + ["mkdir", "-p", temp_root], check=True, timeout=30)
        fasta = f"{temp_root}/subject.fa"; temp_db = f"{temp_root}/subject"
        helper = Path(__file__).resolve().parents[2] / "scripts" / "extract_blast_contig.py"
        helper_text = str(helper.resolve())
        helper_wsl = f"/mnt/{helper_text[0].lower()}/{helper_text[3:].replace(chr(92), '/')}"
        subprocess.run(prefix + ["python3", helper_wsl, target_db, target_contig, fasta], check=True, timeout=600)
        subprocess.run(prefix + ["makeblastdb", "-in", fasta, "-dbtype", "nucl", "-parse_seqids", "-out", temp_db],
                       check=True, timeout=1200, capture_output=True)
        subprocess.run(prefix + ["rm", "-f", fasta], check=True, timeout=30)
        subprocess.run(prefix + ["touch", f"{temp_root}/READY"], check=True, timeout=30)
        subprocess.run(prefix + ["rm", "-rf", root], check=True, timeout=30)
        subprocess.run(prefix + ["mv", temp_root, root], check=True, timeout=30)
        return db
    finally:
        subprocess.run(prefix + ["rm", "-rf", temp_root], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run(prefix + ["rmdir", lock], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
def _parse(text: str, source_start: int | None) -> list[Hit]:
    rows=[]
    for line in text.splitlines():
        c=line.split("\t")
        if len(c)!=9: continue
        # `-db` search yields a real accession in stitle but a BL_ORD_ID in sseqid; `-subject`
        # search is the opposite (sseqid holds the contig name, stitle is "N/A"). Prefer the
        # stitle accession, falling back to sseqid so the target contig name is always resolved.
        sseqid, stitle = c[1], c[2]
        contig = stitle.split()[0] if stitle and stitle != "N/A" else sseqid
        s1,s2=int(c[3]),int(c[4]); q1,q2,qlen=int(c[6]),int(c[7]),int(c[8])
        coverage=100.0*(abs(q2-q1)+1)/max(1,qlen)
        rows.append(Hit(c[0],contig,min(s1,s2),max(s1,s2),"+" if s1<=s2 else "-",float(c[5]),coverage,
                        method="blastn",source_start=source_start,query_start=min(q1,q2),query_end=max(q1,q2),query_length=qlen))
    for row in rows: row.hit_count=len(rows)
    return rows


def _interval_union_length(intervals: list[tuple[int, int]]) -> int:
    total = 0
    for start, end in sorted(intervals):
        if not total:
            current_start, current_end = start, end
            total = current_end - current_start + 1
            continue
        if start > current_end + 1:
            total += end - start + 1
            current_start, current_end = start, end
        elif end > current_end:
            total += end - current_end
            current_end = end
    return total


def _query_overlap_fraction(hit: Hit, members: list[Hit]) -> float:
    if hit.query_start is None or hit.query_end is None:
        return 0.0
    overlap = sum(max(0, min(hit.query_end, x.query_end or 0) - max(hit.query_start, x.query_start or 1) + 1)
                  for x in members if x.query_start is not None and x.query_end is not None)
    return overlap / max(1, hit.query_end - hit.query_start + 1)


def _collapse_loci(rows: list[Hit], max_gap: int = 100_000) -> list[Hit]:
    if not rows:
        return []
    clusters: list[list[Hit]] = []
    for row in sorted(rows, key=lambda h: (h.contig, h.strand, h.start, h.end)):
        members = clusters[-1] if clusters else []
        previous = members[-1] if members else None
        same_locus = (previous and previous.contig == row.contig and previous.strand == row.strand
                      and row.start - max(x.end for x in members) <= max_gap
                      and _query_overlap_fraction(row, members) <= 0.20)
        if same_locus:
            members.append(row)
        else:
            clusters.append([row])
    loci: list[Hit] = []
    for members in clusters:
        first = members[0]
        query_intervals = [(x.query_start, x.query_end) for x in members if x.query_start is not None and x.query_end is not None]
        aligned = _interval_union_length(query_intervals) if query_intervals else 0
        qlen = first.query_length or max((end for _, end in query_intervals), default=1)
        weights = [max(1, (x.query_end or 0) - (x.query_start or 1) + 1) for x in members]
        identity = sum(x.identity*w for x,w in zip(members,weights))/sum(weights)
        loci.append(Hit(first.query_id, first.contig, min(x.start for x in members), max(x.end for x in members),
                        first.strand, round(identity,3), round(100.0*aligned/max(1,qlen),3), method=first.method,
                        source_start=first.source_start, query_start=min((x.query_start or 1) for x in members),
                        query_end=max((x.query_end or 1) for x in members), query_length=qlen))
    loci.sort(key=lambda h: (h.coverage * h.identity, h.end - h.start), reverse=True)
    for locus in loci:
        locus.hit_count = len(loci)
    return loci
def blast_many(queries: list[tuple[str, str, int | None]], target_fasta: str, backend: str, distro: str | None = None, target_db: str | None = None, target_contig: str | None = None, min_identity: float = 0.0, min_coverage: float = 0.0) -> dict[str, list[Hit]]:
    if backend not in ("windows","wsl"): raise ValueError(f"Unsupported BLAST backend: {backend}")
    with tempfile.TemporaryDirectory(prefix="qtlift-") as tmp:
        query=Path(tmp)/"query.fa"; query_text="".join(f">{query_id}\n{sequence}\n" for query_id, sequence, _ in queries); query.write_text(query_text,encoding="ascii")
        if backend=="windows":
            exe=shutil.which("blastn")
            if not exe: raise RuntimeError("Windows blastn is unavailable")
            cmd=[exe,"-query",str(query),"-subject",str(Path(target_fasta).resolve()),"-outfmt",f"6 {FIELDS}","-dust","no","-max_target_seqs","20"]
        else:
            wsl=shutil.which("wsl.exe") or shutil.which("wsl")
            if not wsl: raise RuntimeError("wsl.exe is unavailable")
            prefix=[wsl]+(["-d",distro] if distro else [])
            def wslpath(p: str) -> str: return subprocess.run(prefix+["wslpath","-a","-u",p.replace("\\","/")],check=True,text=True,capture_output=True).stdout.strip()
            if target_db and target_contig:
                target_args = ["-db", _wsl_contig_db(prefix, target_db, target_contig)]
            else:
                target_args = ["-db", target_db] if target_db else ["-subject", wslpath(str(Path(target_fasta).resolve()))]
            cmd=prefix+["blastn","-query","-",*target_args,"-outfmt",f"6 {FIELDS}","-dust","no","-max_target_seqs","20"]
        proc=subprocess.run(cmd,check=True,text=True,capture_output=True,input=query_text if backend=="wsl" else None,timeout=300)
        starts = {query_id: source_start for query_id, _, source_start in queries}
        grouped = {query_id: [] for query_id, _, _ in queries}
        for line in proc.stdout.splitlines():
            query_id = line.split("\t", 1)[0]
            grouped.setdefault(query_id, []).extend(_parse(line, starts.get(query_id)))
        # WSL database searches are restricted before BLAST when target_contig is set.
        # Keep this equality check as validation and for non-database backends.
        def collapse_and_filter(rows: list[Hit]) -> list[Hit]:
            hsps = [r for r in rows if (not target_contig or r.contig == target_contig) and r.identity >= min_identity]
            return [locus for locus in _collapse_loci(hsps) if locus.coverage >= min_coverage]
        return {query_id: collapse_and_filter(rows) for query_id, rows in grouped.items()}


def blast_subject(query_id: str, sequence: str, target_fasta: str, backend: str, source_start: int | None = None, distro: str | None = None, target_db: str | None = None, target_contig: str | None = None, min_identity: float = 0.0, min_coverage: float = 0.0) -> list[Hit]:
    return blast_many([(query_id, sequence, source_start)], target_fasta, backend, distro, target_db, target_contig, min_identity, min_coverage).get(query_id, [])
