# 🧬 genome-toolkit

Reusable scripts for exploring a personal **30× whole-genome sequence** from
[sequencing.com](https://sequencing.com) — from raw reads to a readable
trait & pharmacogenomic report — built so the whole analysis can be
**reproduced from scratch** on any machine.

> **Privacy by design.** This repository contains *only code*. No reads, no
> variant calls, no alignments, and no personal reports are ever committed.
> A broad [`.gitignore`](.gitignore) plus a [pre-commit guard](.githooks/pre-commit)
> enforce that. Point the scripts at your local data and the results land in
> `analysis/`, which is ignored.

---

## What the data looks like

A sequencing.com 30× WGS download contains:

| File | What it is |
|---|---|
| `*.1.fq.gz`, `*.2.fq.gz` | Raw paired-end reads, 150 bp (~60 GB total) |
| `*.snp-indel.genome.vcf.gz` | Small-variant **gVCF** (SNPs + indels, every position) |
| `*.cnv.vcf.gz` | Copy-number variants (Canvas) |
| `*.sv.vcf.gz` | Structural variants (Manta) |
| `reports/*.pdf` | Vendor health-report PDFs |

Reference build: **GRCh38.p13**, contigs named `1`…`22`, `X`, `Y` (no `chr`).
The scripts auto-detect all of this — nothing is hard-coded to one person.

## Quick start (minutes, from the VCFs)

```bash
# 1. install the core toolchain (bcftools / samtools / htslib via Homebrew)
scripts/00_install_tools.sh

# 2. activate the privacy pre-commit guard
scripts/install_git_hooks.sh

# 3. run the whole quick-start pipeline
scripts/run_all.sh
```

By default the scripts look for your data in the repo root (where the
sequencing.com download lands). To keep data elsewhere:

```bash
DATA_DIR=/path/to/your/genome scripts/run_all.sh
```

Everything generated goes to `analysis/` (gitignored).

## Scripts

| Script | Does |
|---|---|
| `00_install_tools.sh` | Install tools. `--align` adds the FASTQ re-alignment toolchain. |
| `01_inspect.sh` | Characterize every input: sizes, formats, read length, build, sample. |
| `02_index_vcfs.sh` | Ensure VCFs are BGZF + tabix-indexed for fast queries. |
| `03_extract_variants.sh` | gVCF → compact `variants.vcf.gz` + `variants.pass.vcf.gz`. |
| `04_vcf_stats.sh` | QC: SNP/indel counts, Ts/Tv, het:hom, per-chrom, inferred sex, SV/CNV counts. |
| `05_genotype_lookup.sh` | Look up your genotype by **rsID / position / region**. |
| `06_trait_report.sh` | Your genotype at a curated set of well-known variants. |
| `10_align_from_fastq.sh` | **From scratch:** FASTQ → CRAM → coverage (hours; downloads GRCh38). |
| `20_clinvar_annotate.sh` | Flag (likely) pathogenic ClinVar variants you carry (downloads ClinVar). |
| `run_all.sh` | The quick-start chain (01→04 + 06). |
| `install_git_hooks.sh` | Turn on the privacy pre-commit guard. |

### Explore individual variants

```bash
scripts/05_genotype_lookup.sh rs429358            # APOE — Alzheimer's
scripts/05_genotype_lookup.sh rs1801133 rs671     # MTHFR, alcohol flush
scripts/05_genotype_lookup.sh 19:44905000-44910000   # a whole region
```

### Reproduce the alignment from raw reads

```bash
scripts/00_install_tools.sh --align          # bwa-mem2, fastp, mosdepth, …
scripts/10_align_from_fastq.sh               # downloads GRCh38, aligns, ~hours
```

## How privacy is enforced

1. **`.gitignore`** denies every sequence/variant/alignment extension, the
   `reports/` dir, the `analysis/` output dir, and anything carrying the
   sample id or owner name.
2. **`.githooks/pre-commit`** independently rejects any staged file that looks
   like genomic data — by extension, by sample id, or by sniffing VCF/BAM/CRAM
   content. Install it with `scripts/install_git_hooks.sh`.
3. The only data committed on purpose is **public reference knowledge** in
   `assets/traits/known_variants.tsv` (rsID → trait maps — not anyone's genotype).

Verify nothing sensitive is staged at any time:

```bash
git status --porcelain        # data files should never appear here
git ls-files | grep -Ei 'vcf|fq|bam|cram' || echo "clean ✓"
```

## ⚠️ Disclaimer

This is for **personal education and curiosity only**. It is **not medical
advice** and **not a diagnostic**. Single-variant trait associations are
simplifications; real risk depends on many genes, ancestry, and environment,
and consumer-grade calls contain errors. Talk to a physician or a certified
genetic counselor before making any health decision. See `docs/PIPELINE.md`
for methods and caveats.
