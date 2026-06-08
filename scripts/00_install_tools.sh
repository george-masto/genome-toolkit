#!/usr/bin/env bash
# ============================================================================
# 00_install_tools.sh — install the bioinformatics toolchain.
#
# Two tiers:
#   CORE  (required for all VCF analysis): bcftools, samtools, htslib
#   ALIGN (only for re-aligning FASTQ from scratch): bwa-mem2, fastp, mosdepth,
#         fastqc, multiqc, seqkit  — installed via bioconda if conda is present.
#
# Usage:
#   scripts/00_install_tools.sh          # core tools only
#   scripts/00_install_tools.sh --align  # core + alignment toolchain
# ============================================================================
source "$(dirname "$0")/lib/common.sh"

WANT_ALIGN=0
[[ "${1:-}" == "--align" ]] && WANT_ALIGN=1

# --- Core tools via Homebrew ------------------------------------------------
install_core() {
  log "Installing core tools (bcftools, samtools, htslib)…"
  if have brew; then
    local missing=()
    for t in bcftools samtools; do have "$t" || missing+=("$t"); done
    have bgzip || missing+=(htslib)
    if [[ ${#missing[@]} -gt 0 ]]; then
      brew install "${missing[@]}"
    else
      ok "core tools already present"
    fi
  elif have conda; then
    conda install -y -c bioconda -c conda-forge bcftools samtools htslib
  else
    die "Need Homebrew or conda to install tools. See https://brew.sh"
  fi
  for t in bcftools samtools bgzip tabix; do
    have "$t" && ok "$t -> $(command -v "$t")" || warn "$t still missing"
  done
}

# --- Alignment toolchain via bioconda ---------------------------------------
install_align() {
  log "Installing alignment toolchain (bioconda)…"
  have conda || die "conda required for the alignment toolchain. Install miniconda first."
  conda install -y -c bioconda -c conda-forge \
    bwa-mem2 fastp mosdepth fastqc multiqc seqkit
  for t in bwa-mem2 fastp mosdepth fastqc seqkit; do
    have "$t" && ok "$t -> $(command -v "$t")" || warn "$t still missing"
  done
}

install_core
[[ $WANT_ALIGN -eq 1 ]] && install_align
log "Done. Re-run with --align to add the FASTQ re-alignment toolchain."
