#!/usr/bin/env bash
# ============================================================================
# 20_clinvar_annotate.sh — flag variants you carry that ClinVar classifies as
# pathogenic or likely-pathogenic.
#
# Downloads NCBI's public ClinVar VCF (GRCh38, ~100 MB), intersects it with
# your PASS variants, and writes the matches with their clinical significance
# and condition. ClinVar uses bare contig names (1,2,…X,Y), matching this data.
#
# Output: analysis/20_clinvar_hits.tsv  (gitignored — your genotype).
#
# ⚠️  A ClinVar "pathogenic" tag does NOT mean you are affected: many are
#     recessive (one copy = carrier), low-penetrance, or condition-specific.
#     Use as a starting point for discussion with a professional, nothing more.
#
# Run 03_extract_variants.sh first.
# ============================================================================
source "$(dirname "$0")/lib/common.sh"
require_tool bcftools

PASS="${PASS:-$OUT_DIR/variants.pass.vcf.gz}"
[[ -f "$PASS" ]] || die "Missing $PASS — run scripts/03_extract_variants.sh first."

CLINVAR_DIR="$REF_DIR/clinvar"; mkdir -p "$CLINVAR_DIR"
CLINVAR_VCF="$CLINVAR_DIR/clinvar.vcf.gz"
CLINVAR_URL="${CLINVAR_URL:-https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz}"

if [[ ! -f "$CLINVAR_VCF" ]]; then
  log "Downloading ClinVar (GRCh38)…"
  curl -L --fail "$CLINVAR_URL"     -o "$CLINVAR_VCF"
  curl -L --fail "$CLINVAR_URL.tbi" -o "$CLINVAR_VCF.tbi" || tabix -p vcf "$CLINVAR_VCF"
fi

OUT="$OUT_DIR/20_clinvar_hits.tsv"
log "Annotating your variants with ClinVar significance…"

# Pull clinical significance (CLNSIG) + disease (CLNDN) onto your variants,
# then keep only (likely) pathogenic calls where you actually carry an ALT.
bcftools annotate -a "$CLINVAR_VCF" -c INFO/CLNSIG,INFO/CLNDN "$PASS" 2>/dev/null \
  | bcftools view -i 'INFO/CLNSIG ~ "athogenic"' 2>/dev/null \
  | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%GT]\t%INFO/CLNSIG\t%INFO/CLNDN\n' 2>/dev/null \
  | awk -F'\t' 'BEGIN{print "chrom\tpos\trsid\tref\talt\tgenotype\tclinical_significance\tcondition"}
                $6!="0/0" && $6!="0|0" {print}' \
  > "$OUT" || true

n=$(($(wc -l < "$OUT" | tr -d ' ') - 1))
ok "Wrote $OUT  ($n (likely) pathogenic ClinVar matches you carry)"
warn "Most matches are carrier-status or low-penetrance. Do NOT self-diagnose."
[[ $n -gt 0 ]] && { hr; column -t -s$'\t' "$OUT" | head -40 >&2; }
