# Pipeline & methods

How the scripts turn whole-genome sequencing into analysis outputs, and the
caveats that come with each step. The reference dataset this was developed
against is a 30× WGS download from sequencing.com, but the scripts read whatever
you point them at — see Inputs.

## Inputs

**Required:** one small-variant **VCF/gVCF** with genotypes (`GT`), ideally
dbSNP-annotated — rsIDs drive the trait / PGx / ClinVar matching, so a VCF with
no rsIDs yields near-empty trait/PGx reports. The scripts auto-detect common
names — sequencing.com (`*.snp-indel.genome.vcf.gz`), DRAGEN
(`*.hard-filtered.vcf.gz`), GATK/DeepVariant (`*.g.vcf.gz`, `*deepvariant*`), or
a plain `genome.vcf.gz` — and fall back to any other `*.vcf.gz`; override with
`VCF_SNP=`. A gVCF (hom-ref blocks via `END=…`/`MIN_DP`) and a plain VCF both work.

**Optional, auto-detected if present:**
- **Raw reads** — paired FASTQ (`*_R1*`/`*_R2*`, `.fq.gz` or `.fastq.gz`) — for
  re-alignment and the second-pass callers.
- **CNVs** — a `*.cnv.vcf.gz` (e.g. Canvas: `<CN0>`, `<DUP>`).
- **SVs** — a `*.sv.vcf.gz` (e.g. Manta: `DEL`, `DUP:TANDEM`, `INV`, `BND`).

**Reference build:** GRCh38 is assumed (ClinVar annotation and the second-pass
callers are GRCh38-specific); trait/PGx lookups are by rsID and build-agnostic.
The toolkit reads and records your build and contig style (`chr1` vs bare `1`)
during inspection. The reference dataset was GRCh38.p13, bare contig names.

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
