![genome-toolkit](assets/banner.png)

# genome-toolkit

![license](https://img.shields.io/badge/license-MIT-blue)
![shell](https://img.shields.io/badge/bash-portable-44cc88)
![input](https://img.shields.io/badge/input-whole--genome%20sequencing-8a2be2)
![privacy](https://img.shields.io/badge/data-never%20committed-e4405f)
![tools](https://img.shields.io/badge/bcftools%20%C2%B7%20samtools-htslib-104756)

Turn **your own whole-genome sequencing** into things you can actually read — QC
numbers, a trait report, a pharmacogenomic (drug-response) card, and a ClinVar
screen — and serve it to your **AI agent as a private MCP server** you can ask
about your genome. A handful of small, composable bash scripts that **reproduce
the whole analysis from scratch** on any Mac or Linux box. Works with WGS from
**any provider** (the author's happens to be 30× from sequencing.com, but
nothing here is tied to that).

> ### 🔒 Privacy is the whole point
> This repository contains **only code**. No reads, no variant calls, no
> alignments, no personal reports — ever. A broad [`.gitignore`](.gitignore)
> *and* an independent [pre-commit guard](.githooks/pre-commit) enforce it.
> Point the scripts at your local data; every result lands in `analysis/`,
> which is ignored. The MCP server runs locally and opens no network port by
> default. **Your genome never leaves your machine.**

---

## ✅ Do I have the right kind of data?

This toolkit needs **whole-genome sequencing (WGS)** — your actual genome, read
base by base. Any one of these works (more of them = more you can run):

| You have… | File looks like | Unlocks |
|---|---|---|
| A **VCF / gVCF** (variant calls) | `*.vcf.gz`, `*.g.vcf.gz` (tens–hundreds of MB) | the **quick-start**: QC, traits, pharmacogenomics, ClinVar, lookups |
| Raw **FASTQ** reads | `*_R1*.fastq.gz` + `*_R2*` (many GB) | re-alignment **and** the second-pass specialty callers |
| An aligned **BAM / CRAM** | `*.bam` / `*.cram` (many GB) | the second-pass specialty callers directly |

From any provider or pipeline that lets you **download your data** — e.g.
Sequencing.com, Dante Labs, Nucleus, ySeq, a research/clinical core, or your own
DRAGEN / GATK / DeepVariant run. The scripts **auto-detect** common file names
(sequencing.com, DRAGEN, GATK, DeepVariant, or a plain `genome.vcf.gz`).

> ### ❌ What does *not* work: 23andMe / AncestryDNA / MyHeritage
> Those "raw data" downloads are from a **genotyping array (a SNP chip)** — they
> read ~600–700K hand-picked spots (a tiny ~15–25 MB text file), **less than 0.1%
> of your genome**. That's not sequencing, can't be converted to it, and this
> toolkit can't use it. Rule of thumb: **WGS files are gigabytes; array files are
> megabytes.**

A note on depth and build (kept honest): the **quick-start runs from any
reasonable WGS VCF**, but single-call confidence tracks the underlying coverage,
and very low-pass/imputed data is shakier on rare variants. The **second-pass
callers genuinely want ~30× reads** (treat anything below ~20× as unreliable)
and are built for **GRCh38 (hg38)** — the build of most current consumer WGS.
Trait/PGx lookups are by rsID and build-agnostic; ClinVar + the second pass are
GRCh38. The toolkit reads and reports your build during inspection.

---

## Quick start

```bash
# 0. (recommended) tell the scripts where your data is — else they look in ./
export DATA_DIR=/path/to/your/genome      # the folder with your VCF / FASTQ

# 1. install the core toolchain (bcftools / samtools / htslib via Homebrew or conda)
scripts/00_install_tools.sh

# 2. turn on the privacy pre-commit guard (so you never commit genetic data)
scripts/install_git_hooks.sh

# 3. run the whole quick-start pipeline (minutes, from your VCF)
scripts/run_all.sh
```

Everything generated goes to `analysis/` (gitignored). If auto-detection can't
find your VCF, the scripts list what they *did* see and tell you to set
`VCF_SNP=/path/to/your.vcf.gz`. That's it.

### 🤖 Or just hand this repo to your AI agent

The whole point is that an agent can set this up *for* you. Open this folder in
**Claude Code** (or Cursor, or any MCP-capable agent), and paste something like:

> *I have whole-genome sequencing data in `/path/to/my/genome` (a `.vcf.gz`, and
> maybe FASTQ). Set up this genome-toolkit: install the tools, run the
> quick-start on my data, build the insights file, and start the MCP server so I
> can ask you about my genome. Tell me if you need anything from me.*

The agent has what it needs: the scripts auto-detect your files, print the exact
MCP config when ready, and keep every output local. From then on you can ask
your agent things like *"any genetic concern with clopidogrel?"* and it queries
your genome directly.

---

## What you get

| Output | From | What it tells you |
|---|---|---|
| **Data inventory** | `01_inspect.sh` | files, sizes, read length, reference build, sample id |
| **QC summary** | `04_vcf_stats.sh` | SNP/indel counts, Ts/Tv, het:hom, inferred sex, SV/CNV tallies |
| **Trait report** | `06_trait_report.sh` | your genotype at well-known trait variants (APOE, MTHFR, lactase…) |
| **Pharmacogenomic card** | `07_pgx_report.sh` | predicted drug-metabolizer phenotypes (CPIC genes) — pharmacist-ready |
| **ClinVar report** | `20_clinvar_annotate.sh` | every ClinVar tier you carry, ranked by review confidence |
| **Genotype lookup** | `05_genotype_lookup.sh` | your genotype at any rsID / position / region, on demand |

Example — the pharmacogenomic profile card *(illustrative format, not real data)*:

```
| Gene    | Genotype | Predicted phenotype       | Key drugs              |
|---------|----------|---------------------------|------------------------|
| CYP2C19 | *1/*2    | Intermediate metabolizer  | clopidogrel, SSRIs     |
| DPYD    | *1/*1    | Normal metabolizer        | fluorouracil (safe)    |
| SLCO1B1 | *1/*1    | Normal function           | simvastatin (safe)     |
```

(QC sanity-ranges assume ~30× WGS of one person; exome or low-pass data will
read differently.)

---

## The scripts

```
scripts/
├── 00_install_tools.sh      install bcftools/samtools/htslib (+ --align / --callers)
├── 01_inspect.sh            characterize every input file
├── 02_index_vcfs.sh         bgzip + tabix index for fast queries
├── 03_extract_variants.sh   gVCF → compact variants.vcf.gz (+ PASS subset)
├── 04_vcf_stats.sh          QC: counts, Ts/Tv, het:hom, inferred sex, SV/CNV
├── 05_genotype_lookup.sh    look up your genotype by rsID / position / region
├── 06_trait_report.sh       genotype at a curated set of well-known variants
├── 07_pgx_report.sh         CPIC pharmacogenomic profile card
├── 10_align_from_fastq.sh   FROM SCRATCH: FASTQ → CRAM → coverage (hours)
│   ── second pass: specialty callers on the CRAM (see docs/SECOND_PASS.md) ──
├── 11_cyp2d6.sh             CYP2D6 diplotype + copy number  (Cyrius)
├── 13_hla.sh                HLA drug-hypersensitivity alleles  (T1K)
├── 14_repeats.sh            repeat-expansion disorders  (ExpansionHunter)
├── 15_smn.sh                SMN1/SMN2 — SMA carrier status  (SMNCopyNumberCaller)
├── 16_gba.sh                GBA — Gaucher / Parkinson's  (Gauchian)
├── 17_mtdna.sh              mtDNA haplogroup + MT-RNR1  (bcftools + Haplogrep)
├── 20_clinvar_annotate.sh   full ClinVar screen, every tier, star-ranked
├── 50_serve_mcp.sh          serve your insights to any AI agent over MCP
├── run_all.sh               the quick-start chain (01→07)
├── install_git_hooks.sh     turn on the privacy pre-commit guard
└── lib/common.sh            shared config + auto-detection of your files
```

### Explore on demand

```bash
scripts/05_genotype_lookup.sh rs429358            # APOE — Alzheimer's
scripts/05_genotype_lookup.sh rs1801133 rs671     # MTHFR + alcohol flush
scripts/05_genotype_lookup.sh 19:44905000-44910000   # a whole region
```

### Reproduce the alignment from raw reads (optional, heavy — needs FASTQ)

The quick-start works from variant calls. If you have the **raw FASTQ**, you can
re-align from scratch — to verify your vendor's calls, compute true coverage, or
unlock the second pass. One command (downloads the reference, several hours,
~80 GB scratch; auto-picks `bwa` on low-RAM machines so it won't OOM):

```bash
scripts/00_install_tools.sh --align     # bwa-mem2 / bwa, fastp, mosdepth, …
scripts/10_align_from_fastq.sh
```

### Second pass: the genes a flat VCF can't type

A few of the most important loci are structurally hard — a near-identical
pseudogene, copy-number variation, an HLA gene with thousands of alleles, a
tandem repeat — so a generic SNP caller returns a dead-end or a wrong answer.
The classic example is **CYP2D6** (codeine, tamoxifen, atomoxetine, many
antidepressants), which a flat VCF reports as *"incompletely typed."* With an
aligned **CRAM** (from your FASTQ, ~30× and GRCh38), purpose-built callers
resolve them:

```bash
scripts/00_install_tools.sh --callers   # ExpansionHunter, T1K, Gauchian, deps
scripts/11_cyp2d6.sh    # CYP2D6 diplotype + copy number  (Cyrius)
scripts/13_hla.sh       # HLA-B*57:01 (abacavir), B*15:02/A*31:01 (carbamazepine)…
scripts/14_repeats.sh   # Huntington, ALS/FTD (C9orf72), Fragile-X, SCAs, DM1…
scripts/15_smn.sh       # SMN1/SMN2 — spinal muscular atrophy carrier status
scripts/16_gba.sh       # GBA — Gaucher disease / Parkinson's risk
scripts/17_mtdna.sh     # mtDNA haplogroup + MT-RNR1 (aminoglycoside ototoxicity)
scripts/50_serve_mcp.sh --rebuild        # MCP picks up the new calls (auto-restarts the service)
```

Each writes a normalized `analysis/calls/*.json`; `build_insights.py` ingests
them and the MCP's `get_second_pass_status` shows what's been run vs is pending.
Full rationale and per-gene detail in [`docs/SECOND_PASS.md`](docs/SECOND_PASS.md).

---

## Serve it to your AI agent (MCP)

Turn your processed genome into a private knowledge source any agent can query —
[Claude Desktop](https://modelcontextprotocol.io), Claude Code, Cursor, or
anything that speaks [MCP](https://modelcontextprotocol.io):

```bash
scripts/50_serve_mcp.sh        # builds a structured insights file + prints the config to paste
```

It distills everything into `analysis/genome_insights.json` and exposes a local
**stdio** MCP server — your agent launches it as a subprocess, so **no network
port is ever opened and nothing leaves your machine**. The agent gets tools like
`check_drug("clopidogrel")`, `get_pharmacogenomics()`, `search_clinvar(gene=…)`,
and a live `lookup_variant(rsID/region)`. Full setup + per-client config in
[`mcp/README.md`](mcp/README.md). To run it always-on for a remote agent over
**Tailscale**, see [`docs/DEPLOY.md`](docs/DEPLOY.md).

```
"Any genetic concern with clopidogrel?"  →  check_drug → CYP2C19 *1/*2,
                                             intermediate metabolizer → prefer prasugrel/ticagrelor
```

## How the data flows

```
 FASTQ (raw reads) ──10──► CRAM (aligned) ──► 11 CYP2D6 · 13 HLA · 14 repeats ──┐
                                              15 SMN · 16 GBA · 17 mtDNA         │ analysis/calls/*.json
 your VCF/gVCF ──03──► variants.vcf.gz ─────┐                                   │
                                            ├─04─► QC summary                   │
                                            ├─06─► trait report                 │
                                            ├─07─► pharmacogenomic card         │
                                            ├─20─► ClinVar screen               │
                                            └─50─► genome_insights.json ◄────────┘
                                                        └─► MCP server ─► your agent
```

Most analyses (the whole right side) run in minutes from the VCF. Only the left
side (FASTQ → CRAM → second pass) needs the raw reads. See
[`docs/PIPELINE.md`](docs/PIPELINE.md) for methods and caveats.

---

## Requirements

- **macOS or Linux**, bash, and either **Homebrew** or **conda** (for the tools).
- Core analysis needs only `bcftools` / `samtools` / `htslib` — installed by
  `scripts/00_install_tools.sh`. Runs on a laptop.
- Re-alignment additionally needs the `--align` toolchain and ~80 GB free disk
  (a roomy external drive is fine — point `OUT_DIR` at it).
- The second-pass callers need the `--callers` tier (conda/bioconda), an aligned
  **CRAM at ~30× on GRCh38**, and several GB of working space.
- Public reference data (the genome reference, ClinVar) is **downloaded on
  demand**, never committed.

---

## How privacy is enforced

1. **`.gitignore`** denies every sequence / variant / alignment extension, the
   `reports/` and `analysis/` directories, the cloned `tools/`, and anything
   carrying a sample id or owner name.
2. **`.githooks/pre-commit`** independently rejects any staged file that *looks*
   like genomic data — by extension, by sample id, or by sniffing VCF/BAM/CRAM
   content. Enable it with `scripts/install_git_hooks.sh`.
3. The only data committed on purpose is **public reference knowledge** —
   `assets/traits/*.tsv`, `assets/pharmacogenes/*.tsv` (rsID → trait/drug maps,
   CYP2D6 + HLA allele functions), and `assets/repeats/*.tsv` (repeat-disorder
   thresholds), which contain **no individual's genotype**.

Sanity-check at any time:

```bash
git ls-files | grep -Ei '\.(vcf|fq|fastq|bam|cram|gvcf)(\.gz)?$' || echo "clean ✓"
```

---

## ⚠️ Disclaimer

For **personal education and curiosity only**. **Not medical advice**, not a
diagnostic, not validated for clinical use. Everything here is *predicted from
consumer-grade WGS*, not produced in a CLIA/CAP lab; single-variant and
short-read calls contain errors and reflect statistical associations, not
certainties. Talk to a physician, pharmacist, or certified genetic counselor —
and confirm anything actionable in a certified lab — before making any decision.

## License

[MIT](LICENSE) — the scripts are yours to use and adapt. Your genome is yours alone.
