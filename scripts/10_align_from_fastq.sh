#!/usr/bin/env bash
# ============================================================================
# 10_align_from_fastq.sh — reproduce the alignment from raw reads, from scratch.
#
# FASTQ (R1+R2)  ->  fastp QC/trim  ->  bwa-mem2  ->  sort  ->  mark duplicates
#                ->  CRAM (ref-compressed)  ->  mosdepth coverage
#
# This is the heavy step: ~60 GB of reads, a ~3 GB reference download, several
# hours of CPU, and ~50–80 GB of working disk. It is NOT run by the quick-start;
# call it explicitly when you want to regenerate alignments yourself.
#
# Requires the alignment toolchain:  scripts/00_install_tools.sh --align
#
# Usage:
#   scripts/10_align_from_fastq.sh            # prompts before the long run
#   FORCE=1 THREADS=8 scripts/10_align_from_fastq.sh   # non-interactive
# ============================================================================
source "$(dirname "$0")/lib/common.sh"
for t in bwa-mem2 fastp samtools mosdepth; do require_tool "$t"; done

[[ -n "$FASTQ_R1" && -f "$FASTQ_R1" ]] || die "No R1 FASTQ detected in $DATA_DIR"
[[ -n "$FASTQ_R2" && -f "$FASTQ_R2" ]] || die "No R2 FASTQ detected in $DATA_DIR"

THREADS="${THREADS:-$(sysctl -n hw.ncpu 2>/dev/null || nproc)}"
SAMPLE="$(sample_name)"; [[ "$SAMPLE" == "unknown" ]] && SAMPLE="sample"
ALIGN_DIR="$OUT_DIR/align"; mkdir -p "$ALIGN_DIR" "$REF_DIR"

# GRCh38 (no-alt analysis set) — matches the build sequencing.com used.
REF_URL="${REF_URL:-https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz}"
REF_FA="$REF_DIR/GRCh38_no_alt.fa"

cat >&2 <<EOF

  This will:
    • download the GRCh38 reference (~900 MB compressed) into $REF_DIR
    • build a bwa-mem2 index (~10 GB, one-time)
    • align $(basename "$FASTQ_R1") + R2 with $THREADS threads
    • mark duplicates and write a CRAM to $ALIGN_DIR
    • compute coverage with mosdepth

  Expect several hours and ~50–80 GB of working disk.

EOF
if [[ "${FORCE:-0}" != "1" ]]; then
  read -r -p "  Proceed? [y/N] " ans
  [[ "$ans" =~ ^[Yy]$ ]] || { log "Aborted."; exit 0; }
fi

# --- 1. reference -----------------------------------------------------------
if [[ ! -f "$REF_FA" ]]; then
  log "Downloading GRCh38 reference…"
  curl -L --fail "$REF_URL" -o "$REF_FA.gz"
  gunzip "$REF_FA.gz"
fi
[[ -f "$REF_FA.fai" ]] || samtools faidx "$REF_FA"
if [[ ! -f "$REF_FA.bwt.2bit.64" ]]; then
  log "Building bwa-mem2 index (one-time, ~10 GB)…"
  bwa-mem2 index "$REF_FA"
fi

# --- 2. fastp QC + adapter/quality trim -------------------------------------
log "fastp QC + trimming…"
R1T="$ALIGN_DIR/${SAMPLE}.trim.R1.fq.gz"
R2T="$ALIGN_DIR/${SAMPLE}.trim.R2.fq.gz"
fastp -i "$FASTQ_R1" -I "$FASTQ_R2" -o "$R1T" -O "$R2T" \
  --thread "$THREADS" \
  --json "$ALIGN_DIR/${SAMPLE}.fastp.json" \
  --html "$ALIGN_DIR/${SAMPLE}.fastp.html"

# --- 3. align -> sort -> markdup -> CRAM ------------------------------------
RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:DNBSEQ\tLB:${SAMPLE}"
CRAM="$ALIGN_DIR/${SAMPLE}.markdup.cram"
log "Aligning with bwa-mem2 (this is the long part)…"
bwa-mem2 mem -t "$THREADS" -R "$RG" "$REF_FA" "$R1T" "$R2T" \
  | samtools fixmate -m -@ "$THREADS" - - \
  | samtools sort -@ "$THREADS" -m 1G - \
  | samtools markdup -@ "$THREADS" --reference "$REF_FA" - "$CRAM"
samtools index "$CRAM"
ok "Wrote $CRAM"

# --- 4. QC: flagstat + coverage ---------------------------------------------
samtools flagstat "$CRAM" > "$ALIGN_DIR/${SAMPLE}.flagstat"
log "Computing coverage with mosdepth…"
( cd "$ALIGN_DIR" && mosdepth -t "$THREADS" -n --fast-mode --by 500 "${SAMPLE}" "$CRAM" )
ok "Coverage summary: $ALIGN_DIR/${SAMPLE}.mosdepth.summary.txt"

cat >&2 <<EOF

  Done. Next steps (separate scripts you can add):
    • call variants from $CRAM with your caller of choice
      (DeepVariant, GATK HaplotypeCaller, or bcftools mpileup) and
    • set REF_FASTA=$REF_FA before re-running 03_extract_variants.sh to
      enable indel left-alignment.
EOF
