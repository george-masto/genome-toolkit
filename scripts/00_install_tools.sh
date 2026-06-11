#!/usr/bin/env bash
# ============================================================================
# 00_install_tools.sh - install the bioinformatics toolchain.
#
# Three tiers:
#   CORE    (required for all VCF analysis): bcftools, samtools, htslib
#   ALIGN   (only for re-aligning FASTQ from scratch): bwa-mem2, fastp, mosdepth,
#           fastqc, multiqc, seqkit  - installed via bioconda if conda is present.
#   CALLERS (second-pass specialty callers that run on the CRAM): ExpansionHunter,
#           T1K (HLA), Gauchian (GBA), plus the python deps Cyrius/SMN need. Cyrius
#           and SMNCopyNumberCaller are git-cloned on demand by their wrappers.
#
# Usage:
#   scripts/00_install_tools.sh            # core tools only
#   scripts/00_install_tools.sh --align    # core + alignment toolchain
#   scripts/00_install_tools.sh --callers  # core + second-pass specialty callers
# ============================================================================
source "$(dirname "$0")/lib/common.sh"

WANT_ALIGN=0; WANT_CALLERS=0
for a in "$@"; do
  case "$a" in
    --align)   WANT_ALIGN=1 ;;
    --callers) WANT_CALLERS=1 ;;
  esac
done

# --- Core tools via Homebrew ------------------------------------------------
install_core() {
  log "Installing core tools (bcftools, samtools, htslib)..."
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
  log "Installing alignment toolchain (bioconda)..."
  have conda || die "conda required for the alignment toolchain. Install miniconda first."
  conda install -y -c bioconda -c conda-forge \
    bwa-mem2 fastp mosdepth fastqc multiqc seqkit
  for t in bwa-mem2 fastp mosdepth fastqc seqkit; do
    have "$t" && ok "$t -> $(command -v "$t")" || warn "$t still missing"
  done
}

# --- Second-pass specialty callers ------------------------------------------
install_callers() {
  log "Installing second-pass specialty callers..."
  have conda || die "conda required for the caller toolchain (bioconda packages). Install miniconda first."
  # conda-installable callers (ExpansionHunter, T1K) + samtools for read extraction
  conda install -y -c bioconda -c conda-forge \
    expansionhunter t1k samtools \
    pysam numpy scipy statsmodels         # deps Cyrius / SMNCopyNumberCaller / Gauchian need
  # Gauchian is pip-only
  if have pip; then pip install -q gauchian || warn "pip install gauchian failed - install it manually."
  else warn "pip not found; install Gauchian with: pip install gauchian"; fi
  for t in ExpansionHunter run-t1k gauchian; do
    have "$t" && ok "$t -> $(command -v "$t")" || warn "$t still missing"
  done
  cat >&2 <<'EOF'
  Notes:
    * Cyrius (CYP2D6) and SMNCopyNumberCaller (SMA) are git-cloned automatically by
      scripts/11_cyp2d6.sh and scripts/15_smn.sh into tools/ on first run.
    * HLA (T1K) builds its IPD-IMGT/HLA index once on first run of 13_hla.sh.
    * Optional for 17_mtdna.sh haplogroup (needs Java): Haplogrep 3 on PATH, OR the
      portable Haplogrep 2 jar (macOS has no v3 build) at tools/haplogrep/haplogrep.jar:
        mkdir -p tools/haplogrep && curl -fsSL -o tools/haplogrep/haplogrep.zip \
          https://github.com/seppinho/haplogrep-cmd/releases/download/v2.4.0/haplogrep.zip \
          && (cd tools/haplogrep && unzip -o haplogrep.zip)
      (MT-RNR1 pharmacogenomics is computed regardless.) Contamination: Haplocheck
      (github.com/genepi/haplocheck), also optional.
EOF
}

install_core
[[ $WANT_ALIGN -eq 1 ]] && install_align
[[ $WANT_CALLERS -eq 1 ]] && install_callers
log "Done. Tiers:  --align (re-align FASTQ)   --callers (CYP2D6/HLA/repeats/SMN/GBA/mtDNA)."
