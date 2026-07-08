# QTLift

[![CI](https://github.com/light-suzuki/QTLift/actions/workflows/ci.yml/badge.svg)](https://github.com/light-suzuki/QTLift/actions/workflows/ci.yml)
[![Release](https://img.shields.io/github/v/release/light-suzuki/=semver)](https://github.com/light-suzuki/QTLift/releases)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)


**Lift a QTL / candidate interval from one reference genome onto another — for any species — with evidence and a confidence score.**

🌐 **[English](#english) · [日本語](#日本語)**

QTLift is a small, local, Windows-first web application. You give it an interval reported on a *source* assembly **B** (for example, the confidence interval of a QTL, a candidate region, or a locus of interest) and it infers **where that interval lies on a canonical *target* assembly A**. It combines several independent lines of evidence — anchor genes / synteny, optional sequence markers, and optional whole-genome liftover — and reports the inferred interval together with its orientation, rationale, warnings, and a High / Medium / Low / Manual-check confidence.

It is **species-agnostic**: anything with a FASTA genome and a GFF3 annotation works — plants, animals, fungi, microbes. Nothing about the tool is tied to a particular organism.

> Built with "vibe coding" — designed and implemented collaboratively with **Claude** and **[OpenAI Codex](https://openai.com/codex/)**. It is open source (MIT); please fork it, adapt it to your organism and workflow, and send improvements back. 🙌

---

## English

### What it is for

Reference genomes are re-assembled and re-versioned constantly. A QTL interval, marker, or candidate gene published against one assembly often needs to be re-located on a newer or different assembly before you can continue the analysis. QTLift automates that re-location and, crucially, **tells you how much to trust the answer** instead of silently emitting coordinates.

Results are **inferred research intervals**, not base-perfect breakpoint claims. QTLift never fabricates evidence: when a tool or an alignment is missing, it says so in the warnings.

### Requirements

- **Windows 10/11** (this is the supported, first-class target)
- **Python 3.13+** and **Node.js 20+** on PATH
- **BLAST+** (`blastn`) — see [External tools](#external-tools). *You provide this.*
- **minimap2** (optional) — only needed for whole-genome liftover. *You provide this.*
- **Genome data** — one FASTA + one GFF3 per reference. *You provide this* (see [Genome data requirements](#genome-data-requirements)).

The application code is self-contained; the genomes and the alignment tools are supplied by you.

### Quick start (Windows)

```powershell
git clone https://github.com/light-suzuki/QTLift.git
cd QTLift
powershell -ExecutionPolicy Bypass -File .\start-qtlift.ps1
```

The first launch creates a Python virtual environment, installs the Python and frontend dependencies, builds the UI, and starts the server. Then open **http://127.0.0.1:8765**.

To stop only QTLift:

```powershell
powershell -ExecutionPolicy Bypass -File .\stop-qtlift.ps1
```

Prefer a one-click, no-terminal launch? Double-click **`QTLift.vbs`** (opens QTLift in Chrome app mode).

Want to see it work immediately with no real data? The bundled artificial sample runs out of the box:

```powershell
py -3.13 .\scripts\create_sample_data.py
$env:PYTHONPATH="$PWD\backend"
py -3.13 .\scripts\run_sample.py
py -3.13 -m unittest discover -s tests -v
```

In the UI, point the genome root at `sample_data\genomes`, pick `RefA` as target A and `RefB` as source B, enter `Chr1:100-850`, and run.

### Genome data requirements

Choose a **genome root folder** whose child folders each hold **one reference**. Each reference folder needs exactly two things:

```text
genomes/
  MyReferenceA/
    genomeA.fa        # or .fasta / .fna, optionally .gz     ← required (sequence)
    genomeA.gff3      # or .gff,           optionally .gz     ← required (annotation)
  MyReferenceB/
    genomeB.fa
    genomeB.gff3
```

| File | Required? | What it provides |
|------|-----------|------------------|
| **FASTA** (`.fa` / `.fasta` / `.fna`, optionally `.gz`) | ✅ Yes | The genome sequence. |
| **GFF3** (`.gff3` / `.gff`, optionally `.gz`) | ✅ Yes | Gene models. QTLift reads `gene` features and their `CDS`; anchors are built from CDS when available (falls back to the whole gene span). |
| **`.fai` index** (e.g. from `samtools faidx genomeA.fa`) | ⭐ Recommended for large genomes | Random access, so QTLift never has to scan multi-GB FASTA files. |
| **BLAST+ nucleotide DB** (`makeblastdb`) | ⭐ Recommended for large genomes | Fast, indexed target search. Without one, QTLift falls back to `blastn -subject` against the FASTA (slower). |

Notes:
- **Contig / chromosome names in the GFF3 should match the FASTA** (e.g. both use `chr1`, or both use `NC_000001.1`). Mismatched names are the most common cause of "no genes found".
- Coordinates you enter in the UI/API are **1-based inclusive**. BED output is 0-based half-open.
- Keep your real genomes **outside this repository**. Point QTLift at wherever they live on disk.

### External tools

QTLift shells out to standard bioinformatics tools that **you install**. It detects them on PATH (or you can pass explicit paths), and any missing tool only disables the evidence step that needs it — the app still starts and runs.

| Tool | Needed for | Required? |
|------|-----------|-----------|
| `blastn` (BLAST+) | Anchor-gene and marker mapping | ✅ Core |
| `makeblastdb` (BLAST+) | Building an indexed target database | ⭐ Recommended |
| `minimap2` | Whole-genome liftover (independent 3rd evidence class) | Optional |

On Windows, running BLAST+ and minimap2 **inside WSL** is well supported (select the *WSL BLAST+* provider); native Windows builds work too.

#### How to install NCBI BLAST+

1. Go to the NCBI BLAST+ downloads page: **https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/**
2. Download the Windows installer (`ncbi-blast-*-win64.exe`) — or, for WSL/Linux, the `*-x64-linux.tar.gz` archive.
3. Install it, then make sure the `bin` folder (containing `blastn.exe` / `makeblastdb.exe`) is on your PATH. Verify with:
   ```powershell
   blastn -version
   ```
4. (Recommended) Build an indexed database for each large target genome:
   ```powershell
   makeblastdb -in genomeA.fa -dbtype nucl -parse_seqids -out genomeA_db
   ```
   `-parse_seqids` keeps the real sequence names available to QTLift.

minimap2 (optional, for liftover): **https://github.com/lh3/minimap2** — download a release binary or `conda install -c bioconda minimap2`, and put it on PATH.

### What if my assembly has no chromosomes (only scaffolds/contigs)?

QTLift works fine without chromosome-scale sequences — it maps anchor genes to whatever contigs they land on. Practical guidance:

- **Leave "Target A chromosome" set to "All sequences"** so the search is not restricted.
- If anchors scatter across many small scaffolds, the synteny state may come back **`split`** and confidence drops to **Manual check** — that is the honest, expected outcome for a fragmented target. Add **sequence markers** to pin the region, and inspect the anchor table.
- If you have a **chromosome-scale assembly for the target A**, prefer it — a contiguous target gives cleaner synteny and higher confidence. The *source* B can still be scaffold-level.
- For a genome with thousands of scaffolds, use the **search box** under the chromosome dropdown to find a specific scaffold by name.

### How it works (in brief)

1. Extract genes overlapping the source-B interval from the GFF3.
2. Select anchor genes by a minimum genomic spacing (no upper cap; a floor guarantees a usable count on short intervals) and BLAST their CDS onto reference A.
3. Optionally map user-supplied **marker sequences** the same way, and optionally run **whole-genome liftover** with minimap2.
4. Evaluate synteny (forward / reverse / partial / split / failed) from anchor order and orientation.
5. Reconcile the evidence and assign confidence. **High** confidence requires at least two independent evidence classes agreeing — so **without markers or liftover, confidence is capped at Medium** by design.

The running app has a **"How it works"** tab with a per-preset anchor-placement table and the full evidence/confidence rules.

### Made with

QTLift was built with **"vibe coding"** — iteratively, in natural language, pair-designing and pair-programming with **Claude** and **Codex**. The scientific contract (never fabricate evidence; always surface uncertainty) was kept front and center throughout.

### Contributing / License

**MIT licensed — please use it, fork it, and adapt it to your own species and pipelines.** Pull requests are very welcome: new evidence sources, better synteny models, non-Windows launchers, additional annotation dialects, and UI improvements are all fair game. See [`AGENTS.md`](AGENTS.md) for an implementation guide written so that a human **or an AI agent** can pick up the codebase and extend it.

---

## 日本語

**QTLift は、ある参照ゲノム B 上で報告された QTL / 候補区間が、別の基準ゲノム A のどこに当たるかを、根拠と信頼度つきで推定するローカルWebアプリです。あらゆる生物種に対応します。**

再アセンブリ・バージョン更新のたびに、旧アセンブリで報告された QTL 区間・マーカー・候補遺伝子を新しいアセンブリ座標へ載せ替える必要が生じます。QTLift はこの載せ替えを自動化し、さらに**その結果をどこまで信頼してよいか**を明示します(座標を黙って出すことはしません)。

**種非依存**です。FASTA ゲノムと GFF3 アノテーションがあれば、植物・動物・菌類・微生物、何でも動きます。特定の生物に依存する要素は一切ありません。

結果は**推定研究区間**であり、塩基単位で完全なブレークポイントの主張ではありません。ツールやアライメントが欠ければ、警告として必ず明示します(根拠の捏造はしません)。

### 動作要件

- **Windows 10/11**(第一級のサポート対象)
- PATH 上の **Python 3.13+** と **Node.js 20+**
- **BLAST+**(`blastn`) — [外部ツール](#外部ツール)参照。**利用者が用意**します。
- **minimap2**(任意) — 全ゲノムリフトオーバーにのみ必要。**利用者が用意**します。
- **ゲノムデータ** — 参照ごとに FASTA 1つ + GFF3 1つ。**利用者が用意**します([ゲノムデータ要件](#ゲノムデータ要件)参照)。

アプリ本体は自己完結しています。ゲノムとアライメントツールは利用者側でご用意ください。

### クイックスタート（Windows）

```powershell
git clone https://github.com/light-suzuki/QTLift.git
cd QTLift
powershell -ExecutionPolicy Bypass -File .\start-qtlift.ps1
```

初回起動時に Python 仮想環境の作成・依存関係のインストール・UI のビルド・サーバー起動を自動で行います。完了したら **http://127.0.0.1:8765** を開いてください。

QTLift だけを停止:

```powershell
powershell -ExecutionPolicy Bypass -File .\stop-qtlift.ps1
```

ターミナルなしのワンクリック起動は **`QTLift.vbs`** をダブルクリック(Chrome アプリモードで開きます)。

実データなしで即試したい場合、同梱の人工サンプルがそのまま動きます:

```powershell
py -3.13 .\scripts\create_sample_data.py
$env:PYTHONPATH="$PWD\backend"
py -3.13 .\scripts\run_sample.py
py -3.13 -m unittest discover -s tests -v
```

UI ではゲノムルートに `sample_data\genomes` を指定し、ターゲット A に `RefA`、ソース B に `RefB` を選び、`Chr1:100-850` を入力して実行します。

### ゲノムデータ要件

**ゲノムルートフォルダー**を選び、その直下の各フォルダーを 1 つの参照とします。各参照フォルダーに必要なのは 2 つだけです。

```text
genomes/
  MyReferenceA/
    genomeA.fa        # または .fasta / .fna（.gz 可）  ← 必須（配列）
    genomeA.gff3      # または .gff          （.gz 可）  ← 必須（アノテーション）
  MyReferenceB/
    genomeB.fa
    genomeB.gff3
```

| ファイル | 必須? | 役割 |
|------|------|------|
| **FASTA**（`.fa`/`.fasta`/`.fna`、`.gz` 可） | ✅ 必須 | ゲノム配列。 |
| **GFF3**（`.gff3`/`.gff`、`.gz` 可） | ✅ 必須 | 遺伝子モデル。`gene` と `CDS` を読み、アンカーは CDS 優先(無ければ遺伝子全長)。 |
| **`.fai` インデックス**（`samtools faidx genomeA.fa` 等） | ⭐ 大きなゲノムでは推奨 | ランダムアクセス。数GBの FASTA を走査せずに済みます。 |
| **BLAST+ 核酸DB**（`makeblastdb`） | ⭐ 大きなゲノムでは推奨 | 高速なインデックス検索。無ければ `blastn -subject`(低速)にフォールバック。 |

補足:
- **GFF3 の contig / 染色体名は FASTA と一致**させてください(両方 `chr1`、または両方 `NC_000001.1` 等)。名前の不一致は「遺伝子が見つからない」最頻の原因です。
- UI/API の座標は **1始まり・両端含む**。BED 出力は 0始まり半開区間です。
- **実ゲノムはこのリポジトリの外**に置き、ディスク上の場所を QTLift に指定してください。

### 外部ツール

QTLift は**利用者がインストールした**標準的なバイオインフォツールを呼び出します。PATH 上から検出(明示的にパス指定も可)し、欠けているツールはその根拠ステップだけを無効化します — アプリの起動・動作は妨げません。

| ツール | 用途 | 必須? |
|------|------|------|
| `blastn`（BLAST+） | アンカー遺伝子・マーカーのマッピング | ✅ 中核 |
| `makeblastdb`（BLAST+） | インデックス済みDBの構築 | ⭐ 推奨 |
| `minimap2` | 全ゲノムリフトオーバー(独立した第3の根拠) | 任意 |

Windows では BLAST+ / minimap2 を **WSL 内**で動かす構成を十分にサポートしています(「WSL版 BLAST+」プロバイダを選択)。Windows ネイティブ版でも動作します。

#### NCBI BLAST+ の入れ方

1. NCBI BLAST+ ダウンロードページを開く: **https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/**
2. Windows 用インストーラ（`ncbi-blast-*-win64.exe`）— WSL/Linux 用なら `*-x64-linux.tar.gz` — を取得。
3. インストール後、`blastn.exe` / `makeblastdb.exe` のある `bin` フォルダーを PATH に追加。確認:
   ```powershell
   blastn -version
   ```
4. （推奨）大きなターゲットゲノムごとにインデックスDBを構築:
   ```powershell
   makeblastdb -in genomeA.fa -dbtype nucl -parse_seqids -out genomeA_db
   ```
   `-parse_seqids` を付けると本来の配列名が QTLift から使えます。

minimap2（任意・リフトオーバー用）: **https://github.com/lh3/minimap2** — リリースバイナリ、または `conda install -c bioconda minimap2` を PATH に。

### 染色体が無い（scaffold/contig だけの）アセンブリの場合

染色体スケールの配列が無くても QTLift は動きます — アンカー遺伝子が乗った contig に対してマッピングします。実務的な指針:

- **「ターゲットA 対象染色体」は「全ての配列」のまま**にして検索を絞り込まないでください。
- アンカーが多数の小さな scaffold に散る場合、シンテニーは **`split`** となり信頼度は **要手動確認** に下がります — 断片化したターゲットでは正直で妥当な結果です。**マーカー配列**を加えて区間を固定し、アンカー表を確認してください。
- **ターゲット A に染色体スケールのアセンブリがある**なら、それを使う方が綺麗なシンテニー・高い信頼度になります。ソース B は scaffold レベルでも構いません。
- scaffold が数千あるゲノムでは、染色体ドロップダウン下の**検索欄**で名前から探せます。

### 仕組み（概要）

1. GFF3 からソース B 区間に重なる遺伝子を抽出。
2. 最小間隔でアンカー遺伝子を選び(上限なし・短い区間には最低本数を保証)、その CDS を参照 A へ BLAST。
3. 任意で**マーカー配列**も同様にマッピングし、任意で minimap2 による**全ゲノムリフトオーバー**を実行。
4. アンカーの並び順・向きからシンテニー(順/逆/部分/分散/失敗)を評価。
5. 根拠を統合し信頼度を判定。**High** には独立した 2 種類以上の根拠の一致が必要 — つまり**マーカーもリフトオーバーも無い場合、信頼度は設計上 Medium 止まり**です。

起動後のアプリには**「仕組み」タブ**があり、プリセット別のアンカー配置表と、根拠・信頼度の全ルールを掲載しています。

### つくりかた

QTLift は **「Vibe コーディング」** で作りました — 自然言語で反復しながら、**Claude** と **Codex** とペア設計・ペアプログラミングして実装しています。「根拠を捏造しない・不確実性を必ず示す」という科学的な約束を最後まで最優先にしました。

### コントリビュート / ライセンス

**MIT ライセンスです。ぜひ使って、フォークして、あなたの生物種やパイプラインに合わせて自由に改変してください。** プルリクエスト大歓迎です — 新しい根拠ソース、より良いシンテニーモデル、非Windowsランチャー、別のアノテーション方言対応、UI改善など、何でも歓迎します。人間でも **AI エージェント**でもコードベースを引き継いで拡張できるよう、実装ガイドを [`AGENTS.md`](AGENTS.md) に用意しています。
