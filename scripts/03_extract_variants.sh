#!/usr/bin/env bash
# ============================================================================
# 03_extract_variants.sh — turn the big gVCF (every position, hom-ref blocks
# and all) into a compact variants-only VCF that downstream steps can use.
#
# Produces (all in analysis/, gitignored):
#   variants.vcf.gz        every non-reference site, multiallelics split
#   variants.pass.vcf.gz   the above, FILTER=PASS only
#
# If a reference FASTA is available (scripts/10_align_from_fastq.sh downloads
# one, or set REF_FASTA=/path/to/GRCh38.fa) indels are also left-aligned.
# ============================================================================
source "$(dirname "$0")/lib/common.sh"
require_tool bcftools

require_vcf_snp   # dies with an actionable message (lists DATA_DIR's .vcf.gz) if none
[[ -f "$VCF_SNP.tbi" ]] || { log "Indexing first…"; tabix -p vcf "$VCF_SNP"; }

ALL="$OUT_DIR/variants.vcf.gz"
PASS="$OUT_DIR/variants.pass.vcf.gz"
THREADS="${THREADS:-4}"

# Optional reference for left-alignment.
REF_FASTA="${REF_FASTA:-}"
if [[ -z "$REF_FASTA" ]]; then
  shopt -s nullglob
  cand=("$REF_DIR"/*.fa "$REF_DIR"/*.fasta "$REF_DIR"/*.fna)
  shopt -u nullglob
  [[ ${#cand[@]} -gt 0 ]] && REF_FASTA="${cand[0]}"
fi

log "Extracting variant sites from $(basename "$VCF_SNP")…"
# -e 'ALT="."'  drops gVCF hom-ref blocks (ALT is a literal '.')
# norm -m -any  splits multiallelic records into biallelic rows
norm_cmd=(bcftools norm -m -any --threads "$THREADS")
if [[ -n "$REF_FASTA" && -f "$REF_FASTA" ]]; then
  ok "Left-aligning against $(basename "$REF_FASTA")"
  norm_cmd+=(-f "$REF_FASTA")
else
  warn "No reference FASTA — skipping left-alignment (split-only). Set REF_FASTA to enable."
fi

bcftools view -e 'ALT="."' --threads "$THREADS" "$VCF_SNP" \
  | "${norm_cmd[@]}" -Oz -o "$ALL"
tabix -f -p vcf "$ALL"
n_all=$(bcftools index -n "$ALL")
ok "variants.vcf.gz: $n_all variant records"

log "Filtering to PASS only…"
bcftools view -f PASS --threads "$THREADS" -Oz -o "$PASS" "$ALL"
tabix -f -p vcf "$PASS"
n_pass=$(bcftools index -n "$PASS")
ok "variants.pass.vcf.gz: $n_pass PASS records"

log "Done. Use these with 04_vcf_stats.sh, 05_genotype_lookup.sh, 06_trait_report.sh."
