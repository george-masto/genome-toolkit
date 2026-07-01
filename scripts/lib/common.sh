#!/usr/bin/env bash
# ============================================================================
# common.sh - shared config, paths, and helpers for all pipeline scripts.
#
# Source this from every script:  source "$(dirname "$0")/lib/common.sh"
#
# It auto-detects the input data files by glob, so nothing is hard-coded to
# one person's filenames. Override any path with an environment variable.
# ============================================================================
set -euo pipefail

# --- Locate repo + data dirs ------------------------------------------------
# REPO_ROOT = the git repo (one level up from scripts/)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# --- User config (optional) -------------------------------------------------
# The one file most users edit: copy config.example.sh -> config.sh and point
# DATA_DIR at your data. Sourced here so its values seed the ${VAR:-default}
# settings below. config.sh uses ": ${VAR:=…}", so environment variables still
# win over it. config.sh is gitignored (it may hold personal paths).
if [[ -f "$REPO_ROOT/config.sh" ]]; then
  # shellcheck source=/dev/null
  source "$REPO_ROOT/config.sh"
fi

# DATA_DIR = where the raw .fq.gz / .vcf.gz live. Defaults to repo root
# (that's where sequencing.com's download landed), override with env var.
DATA_DIR="${DATA_DIR:-$REPO_ROOT}"

# OUT_DIR = all generated output. Gitignored. Override with env var.
OUT_DIR="${OUT_DIR:-$REPO_ROOT/analysis}"

# ASSETS_DIR = committed public reference data (rsID maps, etc.)
ASSETS_DIR="$REPO_ROOT/assets"

# REF_DIR = downloaded reference genome cache (gitignored)
REF_DIR="${REF_DIR:-$REPO_ROOT/reference}"

mkdir -p "$OUT_DIR"

# --- Pretty logging ---------------------------------------------------------
_c_blue='\033[0;34m'; _c_green='\033[0;32m'; _c_yellow='\033[0;33m'
_c_red='\033[0;31m';  _c_dim='\033[2m';      _c_off='\033[0m'
log()   { printf "${_c_blue}==>${_c_off} %s\n" "$*" >&2; }
ok()    { printf "${_c_green}  ✓${_c_off} %s\n" "$*" >&2; }
warn()  { printf "${_c_yellow}  ! ${_c_off}%s\n" "$*" >&2; }
die()   { printf "${_c_red}  ✗ %s${_c_off}\n" "$*" >&2; exit 1; }
hr()    { printf "${_c_dim}%s${_c_off}\n" "----------------------------------------------------------------------" >&2; }

# --- Tool checks ------------------------------------------------------------
have()         { command -v "$1" >/dev/null 2>&1; }
require_tool() { have "$1" || die "Required tool '$1' not found. Run scripts/00_install_tools.sh"; }

# Second-pass callers all need the CRAM (+ its reference). Fail with the exact
# next step rather than a cryptic tool error when the alignment hasn't been run.
require_cram() {
  [[ -n "${CRAM:-}" && -f "$CRAM" ]] || die "No CRAM/BAM found in $ALIGN_DIR.
     The specialty callers run on the aligned reads, not the gVCF. First:
       scripts/00_install_tools.sh --align && scripts/10_align_from_fastq.sh"
  [[ -f "$REF_FA" ]] || warn "Reference $REF_FA not found - CRAM decoding may fail.
     It is created by scripts/10_align_from_fastq.sh (set REF_FA to override)."
  mkdir -p "$CALLS_DIR"
}

# Stamp a normalized call file's generated_utc without needing python/date math
# inside each wrapper. Usage: now_utc -> ISO-8601 string.
now_utc() { date -u +%Y-%m-%dT%H:%M:%SZ; }

# Die with an actionable message when no small-variant VCF/gVCF was found — list
# what *was* in DATA_DIR so the user can point VCF_SNP= at the right file.
require_vcf_snp() {
  [[ -n "${VCF_SNP:-}" && -f "$VCF_SNP" ]] && return 0
  warn "No small-variant VCF/gVCF auto-detected in $DATA_DIR"
  shopt -s nullglob; local found=("$DATA_DIR"/*.vcf.gz); shopt -u nullglob
  if [[ ${#found[@]} -gt 0 ]]; then
    printf '     Found these .vcf.gz files — set VCF_SNP= to the small-variant one:\n' >&2
    printf '       %s\n' "${found[@]}" >&2
  else
    printf '     No .vcf.gz files here. Point at your data: DATA_DIR=/path/to/your/genome\n' >&2
    printf '     This toolkit needs whole-genome SEQUENCING (a VCF/gVCF or FASTQ),\n' >&2
    printf '     not a genotyping-array file (23andMe / AncestryDNA raw data).\n' >&2
  fi
  exit 1
}

# --- Auto-detect the input files (first match of each glob) -----------------
# Use nullglob so a missing file yields an empty string instead of the pattern.
shopt -s nullglob
_first() { local a=("$@"); [[ ${#a[@]} -gt 0 ]] && printf '%s' "${a[0]}" || printf ''; }

# FASTQ R1/R2 — match both .fq.gz and .fastq.gz, and the common R1/_1/.1 styles.
FASTQ_R1="${FASTQ_R1:-$(_first "$DATA_DIR"/*_R1*.f*q.gz "$DATA_DIR"/*_1.f*q.gz "$DATA_DIR"/*.1.f*q.gz "$DATA_DIR"/*.R1.f*q.gz)}"
FASTQ_R2="${FASTQ_R2:-$(_first "$DATA_DIR"/*_R2*.f*q.gz "$DATA_DIR"/*_2.f*q.gz "$DATA_DIR"/*.2.f*q.gz "$DATA_DIR"/*.R2.f*q.gz)}"
# CNV / SV first, so the small-variant fallback below can exclude them.
VCF_CNV="${VCF_CNV:-$(_first "$DATA_DIR"/*.cnv.vcf.gz "$DATA_DIR"/*cnv*.vcf.gz)}"
VCF_SV="${VCF_SV:-$(_first "$DATA_DIR"/*.sv.vcf.gz "$DATA_DIR"/*sv*.vcf.gz)}"
# Primary small-variant VCF/gVCF. Try common vendor/pipeline names, then fall
# back to ANY other *.vcf.gz that isn't the CNV/SV file — so it works for
# sequencing.com (*snp-indel*), DRAGEN (*.hard-filtered), GATK/DeepVariant
# (*.g.vcf.gz / *deepvariant*), or a plain genome.vcf.gz. Override with VCF_SNP=.
VCF_SNP="${VCF_SNP:-$(_first \
  "$DATA_DIR"/*snp-indel*.vcf.gz "$DATA_DIR"/*.hard-filtered.vcf.gz \
  "$DATA_DIR"/*.g.vcf.gz "$DATA_DIR"/*.gvcf.gz "$DATA_DIR"/*.genome.vcf.gz \
  "$DATA_DIR"/*deepvariant*.vcf.gz "$DATA_DIR"/*dragen*.vcf.gz "$DATA_DIR"/*haplotypecaller*.vcf.gz)}"
if [[ -z "$VCF_SNP" ]]; then
  for _v in "$DATA_DIR"/*.vcf.gz; do
    [[ -f "$_v" && "$_v" != "$VCF_CNV" && "$_v" != "$VCF_SV" ]] && { VCF_SNP="$_v"; break; }
  done
fi

# --- Second-pass (CRAM-based) caller infrastructure -------------------------
# The alignment from 10_align_from_fastq.sh lands here, and is the shared input
# for every specialty caller (Cyrius, ExpansionHunter, SMN, GBA, HLA, mtDNA).
ALIGN_DIR="${ALIGN_DIR:-$OUT_DIR/align}"
CRAM="${CRAM:-$(_first "$ALIGN_DIR"/*.markdup.cram "$ALIGN_DIR"/*.cram "$DATA_DIR"/*.cram "$DATA_DIR"/*.bam)}"
# The GRCh38 no-alt reference that 10_align downloaded; CRAM decoding needs it.
REF_FA="${REF_FA:-$REF_DIR/GRCh38_no_alt.fa}"
# Normalized caller outputs that build_insights.py ingests (analysis/calls/*.json).
CALLS_DIR="${CALLS_DIR:-$OUT_DIR/calls}"
shopt -u nullglob

# --- Derived: sample name + reference build (read from the VCF itself) ------
# Cached so we don't re-shell out repeatedly.
sample_name() {
  if [[ -n "${_SAMPLE_NAME:-}" ]]; then printf '%s' "$_SAMPLE_NAME"; return; fi
  if [[ -n "$VCF_SNP" && -f "$VCF_SNP" ]] && have bcftools; then
    # `|| true`: head closing the pipe early can SIGPIPE bcftools (exit 141) on a
    # multi-sample VCF, which would trip `set -o pipefail`. Match ref_build()'s style.
    _SAMPLE_NAME="$(bcftools query -l "$VCF_SNP" 2>/dev/null | head -1 || true)"
  fi
  printf '%s' "${_SAMPLE_NAME:-unknown}"
}

ref_build() {
  if [[ -n "${_REF_BUILD:-}" ]]; then printf '%s' "$_REF_BUILD"; return; fi
  if [[ -n "$VCF_SNP" && -f "$VCF_SNP" ]]; then
    _REF_BUILD="$(gzcat "$VCF_SNP" 2>/dev/null | grep -m1 '^##reference=' | sed 's/##reference=//' || true)"
  fi
  printf '%s' "${_REF_BUILD:-unknown}"
}

# Contig naming style: "chr1" vs "1". Needed so downstream regions match.
chrom_prefix() {
  if [[ -n "${_CHR_PREFIX+x}" ]]; then printf '%s' "$_CHR_PREFIX"; return; fi
  _CHR_PREFIX=""
  if [[ -n "$VCF_SNP" && -f "$VCF_SNP" ]] && have bcftools; then
    if bcftools view -h "$VCF_SNP" 2>/dev/null | grep -q '^##contig=<ID=chr'; then _CHR_PREFIX="chr"; fi
  fi
  printf '%s' "$_CHR_PREFIX"
}

# Print a one-line summary of what was auto-detected.
show_inputs() {
  hr
  log "Detected inputs (DATA_DIR=$DATA_DIR)"
  printf '    %-12s %s\n' "Sample:"   "$(sample_name)"       >&2
  printf '    %-12s %s\n' "Reference:" "$(ref_build)"        >&2
  printf '    %-12s %s\n' "FASTQ R1:"  "${FASTQ_R1:-<none>}" >&2
  printf '    %-12s %s\n' "FASTQ R2:"  "${FASTQ_R2:-<none>}" >&2
  printf '    %-12s %s\n' "SNP/indel:" "${VCF_SNP:-<none>}"  >&2
  printf '    %-12s %s\n' "CNV:"       "${VCF_CNV:-<none>}"  >&2
  printf '    %-12s %s\n' "SV:"        "${VCF_SV:-<none>}"   >&2
  hr
}
