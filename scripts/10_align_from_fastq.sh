#!/usr/bin/env bash
# ============================================================================
# 10_align_from_fastq.sh - reproduce the alignment from raw reads, from scratch.
#
# FASTQ (R1+R2)  ->  fastp QC/trim  ->  bwa-mem2  ->  sort  ->  mark duplicates
#                ->  CRAM (ref-compressed)  ->  mosdepth coverage
#
# This is the heavy step: ~60 GB of reads, a ~3 GB reference download, several
# hours of CPU, and ~50-80 GB of working disk. It is NOT run by the quick-start;
# call it explicitly when you want to regenerate alignments yourself.
#
# Requires the alignment toolchain:  scripts/00_install_tools.sh --align
#
# Usage:
#   scripts/10_align_from_fastq.sh            # prompts before the long run
#   FORCE=1 THREADS=8 scripts/10_align_from_fastq.sh   # non-interactive
# ============================================================================
source "$(dirname "$0")/lib/common.sh"
for t in fastp samtools mosdepth; do require_tool "$t"; done

[[ -n "$FASTQ_R1" && -f "$FASTQ_R1" ]] || die "No R1 FASTQ detected in $DATA_DIR"
[[ -n "$FASTQ_R2" && -f "$FASTQ_R2" ]] || die "No R2 FASTQ detected in $DATA_DIR"

THREADS="${THREADS:-$(sysctl -n hw.ncpu 2>/dev/null || nproc)}"
SAMPLE="$(sample_name)"; [[ "$SAMPLE" == "unknown" ]] && SAMPLE="sample"
ALIGN_DIR="$OUT_DIR/align"; mkdir -p "$ALIGN_DIR" "$REF_DIR"

# --- aligner selection (RAM-aware) ------------------------------------------
# bwa-mem2 is faster, but building its human-genome index needs ~40-70 GB RAM
# and OOM-kills on smaller machines. Fall back to classic bwa there (~4 GB to
# index, a bit slower to align, identical downstream - the specialty callers are
# validated on bwa-mem alignments). Override with ALIGNER=bwa|bwa-mem2.
MEM_GB=$(( $(sysctl -n hw.memsize 2>/dev/null || echo 0) / 1073741824 ))
if [[ -z "${ALIGNER:-}" ]]; then
  if   have bwa-mem2 && (( MEM_GB >= 48 )); then ALIGNER="bwa-mem2"
  elif have bwa;                            then ALIGNER="bwa"
  elif have bwa-mem2;                        then ALIGNER="bwa-mem2"
  else die "Need 'bwa' or 'bwa-mem2'. Install: scripts/00_install_tools.sh --align"
  fi
fi
require_tool "$ALIGNER"
log "Aligner: $ALIGNER (system RAM ${MEM_GB} GB, $THREADS threads)"

# GRCh38 (no-alt analysis set) - matches the build sequencing.com used.
REF_URL="${REF_URL:-https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz}"
REF_FA="$REF_DIR/GRCh38_no_alt.fa"

cat >&2 <<EOF

  This will:
    * download the GRCh38 reference (~900 MB compressed) into $REF_DIR
    * build the aligner index (one-time; bwa ~4 GB RAM, bwa-mem2 ~40+ GB RAM)
    * align $(basename "$FASTQ_R1") + R2 with $THREADS threads
    * mark duplicates and write a CRAM to $ALIGN_DIR
    * compute coverage with mosdepth

  Expect several hours and ~50-80 GB of working disk.

EOF
if [[ "${FORCE:-0}" != "1" ]]; then
  read -r -p "  Proceed? [y/N] " ans
  [[ "$ans" =~ ^[Yy]$ ]] || { log "Aborted."; exit 0; }
fi

# --- 1. reference -----------------------------------------------------------
if [[ ! -f "$REF_FA" ]]; then
  log "Downloading GRCh38 reference..."
  curl -L --fail "$REF_URL" -o "$REF_FA.gz"
  gunzip "$REF_FA.gz"
fi
[[ -f "$REF_FA.fai" ]] || samtools faidx "$REF_FA"
# Index sentinel differs by aligner (bwa: .bwt; bwa-mem2: .bwt.2bit.64).
if [[ "$ALIGNER" == "bwa-mem2" ]]; then INDEX_SENTINEL="$REF_FA.bwt.2bit.64"
else INDEX_SENTINEL="$REF_FA.bwt"; fi
if [[ ! -f "$INDEX_SENTINEL" ]]; then
  log "Building $ALIGNER index (one-time)..."
  "$ALIGNER" index "$REF_FA"
fi

# --- 2. fastp QC + adapter/quality trim -------------------------------------
log "fastp QC + trimming..."
R1T="$ALIGN_DIR/${SAMPLE}.trim.R1.fq.gz"
R2T="$ALIGN_DIR/${SAMPLE}.trim.R2.fq.gz"
fastp -i "$FASTQ_R1" -I "$FASTQ_R2" -o "$R1T" -O "$R2T" \
  --thread "$THREADS" \
  --json "$ALIGN_DIR/${SAMPLE}.fastp.json" \
  --html "$ALIGN_DIR/${SAMPLE}.fastp.html"

# --- 3. align -> sort -> markdup -> CRAM ------------------------------------
RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:DNBSEQ\tLB:${SAMPLE}"
CRAM="$ALIGN_DIR/${SAMPLE}.markdup.cram"
# samtools sort buffers are per-thread; cap total so the aligner (~5 GB for the
# human index) + sort fit in RAM on a small machine. Sort temp -> ALIGN_DIR
# (which may be a roomy external drive), not the default /tmp.
SORT_THREADS="${SORT_THREADS:-$(( THREADS < 4 ? THREADS : 4 ))}"
SORT_MEM="${SORT_MEM:-768M}"
log "Aligning with $ALIGNER (this is the long part)..."
"$ALIGNER" mem -t "$THREADS" -R "$RG" "$REF_FA" "$R1T" "$R2T" \
  | samtools fixmate -m -@ "$THREADS" - - \
  | samtools sort -@ "$SORT_THREADS" -m "$SORT_MEM" -T "$ALIGN_DIR/sorttmp.${SAMPLE}" - \
  | samtools markdup -@ "$THREADS" --reference "$REF_FA" - "$CRAM"
samtools index "$CRAM"
ok "Wrote $CRAM"

# --- 4. QC: flagstat + coverage ---------------------------------------------
samtools flagstat "$CRAM" > "$ALIGN_DIR/${SAMPLE}.flagstat"
log "Computing coverage with mosdepth..."
# mosdepth needs the reference to decode a CRAM (-f / REF_PATH), else it errors.
( cd "$ALIGN_DIR" && mosdepth -t "$THREADS" -f "$REF_FA" -n --fast-mode --by 500 "${SAMPLE}" "$CRAM" )
ok "Coverage summary: $ALIGN_DIR/${SAMPLE}.mosdepth.summary.txt"

cat >&2 <<EOF

  Done. This CRAM unlocks the second-pass specialty callers - the genes a flat
  VCF can't type. Install them once, then run any/all:
      scripts/00_install_tools.sh --callers
      scripts/11_cyp2d6.sh    # CYP2D6 diplotype + copy number (Cyrius)
      scripts/13_hla.sh       # HLA drug-hypersensitivity alleles (T1K)
      scripts/14_repeats.sh   # repeat-expansion disorders (ExpansionHunter)
      scripts/15_smn.sh       # SMN1/2 - SMA carrier (SMNCopyNumberCaller)
      scripts/16_gba.sh       # GBA - Gaucher / Parkinson's (Gauchian)
      scripts/17_mtdna.sh     # mtDNA haplogroup + MT-RNR1 (bcftools + Haplogrep)
      scripts/50_serve_mcp.sh --rebuild     # MCP picks up the new calls
  See docs/SECOND_PASS.md. (Or call variants from $CRAM with DeepVariant/GATK and
  set REF_FASTA=$REF_FA before re-running 03_extract_variants.sh for left-alignment.)
EOF
