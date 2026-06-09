![genome-toolkit](assets/banner.png)

# genome-toolkit

![license](https://img.shields.io/badge/license-MIT-blue)
![shell](https://img.shields.io/badge/bash-portable-44cc88)
![reference](https://img.shields.io/badge/reference-GRCh38-8a2be2)
![privacy](https://img.shields.io/badge/data-never%20committed-e4405f)
![tools](https://img.shields.io/badge/bcftools%20%C2%B7%20samtools-htslib-104756)

Turn a consumer **30× whole-genome sequence** (e.g. from
[sequencing.com](https://sequencing.com)) into things you can actually read —
QC numbers, a trait report, a pharmacogenomic profile card, and a ClinVar
screen — with a handful of small, composable scripts that **reproduce the whole
analysis from scratch** on any machine.

> ### 🔒 Privacy is the whole point
> This repository contains **only code**. No reads, no variant calls, no
> alignments, no personal reports — ever. A broad [`.gitignore`](.gitignore)
> *and* an independent [pre-commit guard](.githooks/pre-commit) enforce it.
> Point the scripts at your local data; every result lands in `analysis/`,
> which is ignored. Your genome never touches the repo.

---

## Quick start

```bash
# 1. install the toolchain (bcftools / samtools / htslib via Homebrew or conda)
scripts/00_install_tools.sh

# 2. turn on the privacy pre-commit guard
scripts/install_git_hooks.sh

# 3. run the whole quick-start pipeline (minutes, from the VCFs)
scripts/run_all.sh
```

The scripts **auto-detect** your data files by pattern — nothing is hard-coded
to one person or filename. Keep your download in the repo root, or point at it:

```bash
DATA_DIR=/path/to/your/genome scripts/run_all.sh
```

Everything generated goes to `analysis/` (gitignored). That's it.

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

---

## The scripts

```
scripts/
├── 00_install_tools.sh      install bcftools/samtools/htslib (+ --align toolchain)
├── 01_inspect.sh            characterize every input file
├── 02_index_vcfs.sh         bgzip + tabix index for fast queries
├── 03_extract_variants.sh   gVCF → compact variants.vcf.gz (+ PASS subset)
├── 04_vcf_stats.sh          QC: counts, Ts/Tv, het:hom, inferred sex, SV/CNV
├── 05_genotype_lookup.sh    look up your genotype by rsID / position / region
├── 06_trait_report.sh       genotype at a curated set of well-known variants
├── 07_pgx_report.sh         CPIC pharmacogenomic profile card
├── 10_align_from_fastq.sh   FROM SCRATCH: FASTQ → CRAM → coverage (hours)
├── 20_clinvar_annotate.sh   full ClinVar screen, every tier, star-ranked
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

### Reproduce the alignment from raw reads (optional, heavy)

The quick-start works from the variant calls. To independently re-align the raw
FASTQ — to verify the vendor's calls, compute true coverage, or run specialized
callers — it's one command (downloads GRCh38, several hours, ~80 GB scratch):

```bash
scripts/00_install_tools.sh --align     # bwa-mem2, fastp, mosdepth, …
scripts/10_align_from_fastq.sh
```

---

## How the data flows

```
 FASTQ (raw reads) ──10──► CRAM (aligned) ──┐
                                            ├─► specialized callers (DeepVariant, Cyrius…)
 vendor gVCF ──03──► variants.vcf.gz ───────┤
                                            ├─04─► QC summary
                                            ├─06─► trait report
                                            ├─07─► pharmacogenomic card
                                            └─20─► ClinVar screen
```

Most analyses (the whole right side) run in minutes from the gVCF. Only the
left side (FASTQ → CRAM) is the long, optional re-alignment. See
[`docs/PIPELINE.md`](docs/PIPELINE.md) for methods and caveats.

---

## Requirements

- **macOS or Linux**, bash, and either **Homebrew** or **conda** (for the tools).
- Core analysis needs only `bcftools` / `samtools` / `htslib` — installed by
  `scripts/00_install_tools.sh`.
- Re-alignment additionally needs the `--align` toolchain and ~80 GB free disk.
- Public reference data (GRCh38, ClinVar) is **downloaded on demand**, never committed.

---

## How privacy is enforced

1. **`.gitignore`** denies every sequence / variant / alignment extension, the
   `reports/` and `analysis/` directories, and anything carrying a sample id or
   owner name.
2. **`.githooks/pre-commit`** independently rejects any staged file that *looks*
   like genomic data — by extension, by sample id, or by sniffing VCF/BAM/CRAM
   content. Enable it with `scripts/install_git_hooks.sh`.
3. The only data committed on purpose is **public reference knowledge** —
   `assets/traits/*.tsv` and `assets/pharmacogenes/*.tsv` (rsID → trait/drug
   maps), which contain **no individual's genotype**.

Sanity-check at any time:

```bash
git ls-files | grep -Ei '\.(vcf|fq|fastq|bam|cram|gvcf)(\.gz)?$' || echo "clean ✓"
```

---

## ⚠️ Disclaimer

For **personal education and curiosity only**. **Not medical advice**, not a
diagnostic, not validated for clinical use. Single-variant and consumer-grade
calls contain errors and reflect statistical associations, not certainties.
Pharmacogenomic predictions here are *predicted from WGS*, not a clinical PGx
test. Talk to a physician, pharmacist, or certified genetic counselor — and
confirm anything actionable in a CLIA/CAP lab — before making any decision.

## License

[MIT](LICENSE) — the scripts are yours to use and adapt. Your genome is yours alone.
