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
# Annotate with ClinVar significance (CLNSIG), disease (CLNDN), gene (GENEINFO),
# and review status (CLNREVSTAT — ClinVar's confidence / "star" level).
ANNOT="$OUT_DIR/.clinvar_annot.vcf.gz"
bcftools annotate -a "$CLINVAR_VCF" \
  -c INFO/CLNSIG,INFO/CLNDN,INFO/GENEINFO,INFO/CLNREVSTAT \
  -Oz -o "$ANNOT" "$PASS" 2>/dev/null
tabix -f -p vcf "$ANNOT" 2>/dev/null || true

# Real actionable tiers only: "Pathogenic", "Likely_pathogenic",
# "Pathogenic/Likely_pathogenic" (incl. *,_low_penetrance). Explicitly EXCLUDE
# "Conflicting_classifications_of_pathogenicity" — that's the disagree-among-
# submitters tier, not an actionable call. Keep only sites you carry an ALT for.
# NOTE: bcftools '&& ... !~' does NOT reliably exclude a substring here, so we
# CHAIN a second view with -e (exclude) instead of one compound expression.
# emit() reads an already-filtered VCF stream on stdin and writes a TSV to $1.
emit() {  # $1 = outfile
  { echo -e "stars\tchrom\tpos\trsid\tgene\tref\talt\tgenotype\tzygosity\tclinical_significance\treview_status\tcondition"
    bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/GENEINFO\t%REF\t%ALT\t[%GT]\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/CLNDN\n' 2>/dev/null \
      | awk -F'\t' 'BEGIN{OFS="\t"}
          $7!="0/0"&&$7!="0|0"&&$7!="0"{
            z=($7=="1/1"||$7=="1|1")?"homozygous":(($7=="1")?"hemizygous":"heterozygous")
            gsub(/:[0-9]+/,"",$4)            # GENEINFO "GENE:id" -> "GENE"
            # ClinVar review status -> star confidence (0-4)
            st=0
            if      ($9 ~ /practice_guideline/)                 st=4
            else if ($9 ~ /reviewed_by_expert_panel/)           st=3
            else if ($9 ~ /multiple_submitters,_no_conflicts/)  st=2
            else if ($9 ~ /single_submitter/ || $9 ~ /conflicting/) st=1
            else                                                st=0
            print st,$1,$2,$3,$4,$5,$6,$7,z,$8,$9,$10 }' \
      | sort -t$'\t' -k1,1nr -k9,9 -k2,2V -k3,3n   # highest confidence first
  } > "$1"
}

# Actionable tier: contains "pathogenic", but NOT the conflicting/uncertain bucket.
bcftools view -i 'INFO/CLNSIG ~ "[Pp]athogenic"' "$ANNOT" 2>/dev/null \
  | bcftools view -e 'INFO/CLNSIG ~ "Conflicting"' 2>/dev/null \
  | emit "$OUT"
n=$(($(wc -l < "$OUT" | tr -d ' ') - 1)); [[ $n -lt 0 ]] && n=0
n_hom=$(awk -F'\t' 'NR>1 && $8=="homozygous"' "$OUT" | wc -l | tr -d ' ')

# Lower-confidence "conflicting" tier, set aside for context.
CONFLICT="$OUT_DIR/20_clinvar_conflicting.tsv"
bcftools view -i 'INFO/CLNSIG ~ "Conflicting"' "$ANNOT" 2>/dev/null | emit "$CONFLICT"
n_conf=$(($(wc -l < "$CONFLICT" | tr -d ' ') - 1)); [[ $n_conf -lt 0 ]] && n_conf=0
rm -f "$ANNOT" "$ANNOT.tbi"

ok "Wrote $OUT"
ok "  $n Pathogenic / Likely-pathogenic variants you carry ($n_hom homozygous, rest heterozygous = usually recessive carrier)"
warn "  Set aside $n_conf 'conflicting classification' variants (submitters disagree) -> $CONFLICT"
warn "Heterozygous + recessive condition = healthy CARRIER, not affected. Do NOT self-diagnose; review with a genetic counselor."
[[ $n -gt 0 ]] && { hr; column -t -s$'\t' "$OUT" | cut -c1-200 | head -40 >&2; }
