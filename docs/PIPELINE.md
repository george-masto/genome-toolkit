# Pipeline & methods

How the scripts turn a sequencing.com 30× WGS download into analysis outputs,
and the caveats that come with each step.

## Inputs

sequencing.com performs secondary analysis and ships:

- **Raw reads** — 150 bp paired-end FASTQ, MGI/DNBSEQ platform.
- **Small variants** — a **gVCF** (`*.snp-indel.genome.vcf.gz`): every position
  is represented, either as a variant row or a hom-ref block (`END=…`, `MIN_DP`).
  dbSNP rsIDs are pre-annotated. FORMAT carries `GT:DP:GQ:AD`.
- **CNVs** — called by **Canvas 1.40**; `<CN0>`, `<DUP>`, etc.
- **SVs** — called by **Manta 1.6**; `DEL`, `DUP:TANDEM`, `INV`, `BND`.

Reference: **GRCh38.p13**, bare contig names (`1`, not `chr1`).

## Steps

### 01 — inspect
Reads only the first FASTQ read (stream terminated early with `head`) and the
VCF headers. Records sizes, formats (`htsfile`), read length, sample id, build,
and caller versions to `analysis/01_data_inventory.md`.

### 02 — index
sequencing.com already ships BGZF, so this just writes `.tbi` sidecars next to
the originals. A non-BGZF file would be re-compressed into `analysis/` (originals
are never modified).

### 03 — extract variants
The gVCF is mostly hom-ref blocks. This step:
1. `bcftools view -e 'ALT="."'` — drop hom-ref blocks, keep variant sites.
2. `bcftools norm -m -any` — split multiallelic records into biallelic rows.
3. (optional) `-f REF.fa` — left-align indels if a reference FASTA is available.

Outputs `variants.vcf.gz` (all) and `variants.pass.vcf.gz` (FILTER=PASS).

> **Why split multiallelics?** Downstream rsID/genotype logic is simplest when
> every row has exactly one ALT. Left-alignment is skipped unless you provide a
> reference (run `10_align_from_fastq.sh`, or set `REF_FASTA=`), which only
> affects how some indels are represented, not SNP lookups.

### 04 — QC stats
`bcftools stats` plus derived metrics. Sanity benchmarks for a human 30× WGS:

| Metric | Typical | Meaning if off |
|---|---|---|
| SNPs | 3.5–4.5 M | very high → false positives; low → under-calling |
| Indels | 0.6–1.0 M | — |
| Ts/Tv | ~2.0–2.1 | < 1.8 → noise / artifacts |
| Het:Hom | ~1.3–2.0 | extreme → contamination or ancestry mismatch |

**Genetic sex** is inferred from chrX heterozygosity and chrY call count
(heuristic, not clinical). **SV/CNV** rows are tallied by type.

### 06 — trait report
For each variant in `assets/traits/known_variants.tsv`, resolves the rsID against
your own variants (authoritative) and reports genotype + copies of the effect
allele. Sites absent from the variants file are homozygous reference (the
hom-ref blocks were stripped in step 03). The script flags any REF/ALT mismatch
between the catalog and your data, so a mis-curated row can't be silently
misread.

### 07 — pharmacogenomic profile
Screens a curated **CPIC pharmacogene** catalog
(`assets/pharmacogenes/pgx_variants.tsv`) — only star-allele-defining variants a
single-sample short-read VCF can read reliably — and predicts a
metabolizer/function phenotype per gene (activity-score model for the CYP genes,
bespoke rules for VKORC1, UGT1A1, CYP3A5, NAT2, SLCO1B1, etc.). Outputs a
pharmacist-ready profile card with CPIC drug implications, separating
non-default phenotypes from the reassuring normals.

> **CYP2D6 is reported as incomplete on purpose.** Its star alleles depend on
> copy number and the CYP2D7 pseudogene, which short reads can't resolve — a
> dedicated CNV/paralog-aware caller (Cyrius, Aldy) on the CRAM is required for a
> true CYP2D6 phenotype. This is one reason to run step 10.

### 10 — re-align from scratch (optional, heavy)
`fastp` → `bwa-mem2 mem` → `samtools fixmate|sort|markdup` → **CRAM**, then
`mosdepth` coverage. Downloads the GRCh38 no-alt analysis set. Several hours and
~50–80 GB of working disk. Use this to independently reproduce alignments rather
than trusting the vendor's BAM (which isn't shipped anyway).

### 20 — ClinVar screen (optional)
Annotates your PASS variants with ClinVar `CLNSIG`/`CLNDN` and keeps the
(likely) pathogenic ones you actually carry. **Most hits are recessive carrier
status or low-penetrance** — not a diagnosis.

## Limitations

- Single-sample, no trio/phasing → no de-novo or compound-het resolution.
- Consumer short-read calls miss repeat expansions, complex SVs, and much of the
  pharmacogenome (e.g. *CYP2D6* star-alleles need specialized callers).
- Trait associations are population-statistical and ancestry-dependent.
- Nothing here is validated for clinical use. Confirm anything actionable in a
  CLIA/CAP lab.
