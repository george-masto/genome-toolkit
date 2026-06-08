#!/usr/bin/env bash
# ============================================================================
# 20_clinvar_annotate.sh — categorize EVERY ClinVar tier you carry, not just
# pathogenic: pharmacogenomics (drug_response), risk-factor, protective,
# association, uncertain, conflicting, benign.
#
# Downloads NCBI's public ClinVar VCF (GRCh38, ~190 MB), annotates your PASS
# variants, buckets each carried variant into one tier, and rates each by
# ClinVar review confidence (0–4 "stars"). ClinVar uses bare contig names
# (1,2,…X,Y), matching this data.
#
# Outputs (all gitignored — your genotype):
#   analysis/20_clinvar_report.md            human-readable, every tier
#   analysis/20_clinvar_all_categorized.tsv  one row per carried ClinVar variant
#   analysis/20_clinvar_pathogenic.tsv       the pathogenic tier
#   analysis/20_clinvar_drug_response.tsv    the pharmacogenomic tier
#
# ⚠️  A ClinVar tag is NOT a diagnosis. ★ matters: most 0★ rows are common
#     variants with a single unreviewed submission. Heterozygous + recessive
#     condition = healthy carrier. Discuss anything actionable with a pro.
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

log "Annotating your variants with ClinVar (CLNSIG, CLNDN, GENEINFO, CLNREVSTAT)…"
ANNOT="$OUT_DIR/.clinvar_annot.vcf.gz"
bcftools annotate -a "$CLINVAR_VCF" \
  -c INFO/CLNSIG,INFO/CLNDN,INFO/GENEINFO,INFO/CLNREVSTAT \
  -Oz -o "$ANNOT" "$PASS" 2>/dev/null
tabix -f -p vcf "$ANNOT" 2>/dev/null || true

# Build ONE categorized table of every ClinVar-annotated variant you carry.
# Each row is bucketed into a single tier and given a 0-4 "star" confidence
# from CLNREVSTAT. Columns: bucket, stars, chrom, pos, gene, rsid, genotype,
# zygosity, clinical_significance, condition.
ALLCATS="$OUT_DIR/20_clinvar_all_categorized.tsv"
log "Categorizing every tier you carry…"
bcftools view -i 'INFO/CLNSIG!="."' "$ANNOT" 2>/dev/null \
 | bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/GENEINFO\t%REF\t%ALT\t[%GT]\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/CLNDN\n' 2>/dev/null \
 | awk -F'\t' 'BEGIN{OFS="\t"}
     $7!="0/0"&&$7!="0|0"&&$7!="0" {
       chrom=$1;pos=$2;rsid=$3;gene=$4;ref=$5;alt=$6;gt=$7;sig=$8;rev=$9;dn=$10
       gsub(/:[0-9]+/,"",gene); split(gene,gg,"|"); gene=gg[1]; if(gene=="")gene="-"
       n=split(gt,a,/[\/|]/); gb=""
       for(i=1;i<=n;i++){ b=(a[i]==0?ref:(a[i]==1?alt:"?")); gb=gb (i>1?"/":"") b }
       z=(gt=="1/1"||gt=="1|1")?"hom":((gt=="1")?"hemi":"het")
       st=0
       if      (rev~/practice_guideline/)                st=4
       else if (rev~/reviewed_by_expert_panel/)          st=3
       else if (rev~/multiple_submitters,_no_conflicts/) st=2
       else if (rev~/single_submitter/||rev~/conflicting/) st=1
       if      (sig~/[Pp]athogenic/ && sig!~/Conflicting/)  bk="1_Pathogenic"
       else if (sig~/drug_response|confers_sensitivity/)    bk="2_DrugResponse"
       else if (sig~/risk_allele|risk_factor|Established_risk/) bk="3_RiskFactor"
       else if (sig~/protective/)                           bk="4_Protective"
       else if (sig~/association|Affects/)                  bk="5_Association"
       else if (sig~/Conflicting/)                          bk="8_Conflicting"
       else if (sig~/Uncertain/)                            bk="6_Uncertain"
       else if (sig~/[Bb]enign/)                            bk="9_Benign"
       else                                                 bk="7_Other"
       ncd=split(dn,cc,"|"); cond=""
       for(i=1;i<=ncd;i++){ t=cc[i]
         if(t==""||t=="not_provided"||t=="not_specified"||t=="See_cases")continue
         if(bk=="2_DrugResponse"&&t~/[Rr]esponse/){cond=t;break}
         if(cond=="")cond=t }
       if(cond=="")cond="-"; gsub(/_/," ",cond); if(length(cond)>55)cond=substr(cond,1,52)"…"
       print bk,st,chrom,pos,gene,rsid,gb,z,sig,cond }' \
 > "$ALLCATS"

rm -f "$ANNOT" "$ANNOT.tbi"

# ---- markdown report ------------------------------------------------------
REPORT="$OUT_DIR/20_clinvar_report.md"
count_of(){ awk -F'\t' -v b="$1" '$1==b' "$ALLCATS" | wc -l | tr -d ' '; }
summary_row(){ local c; c=$(count_of "$1"); [[ "$c" -gt 0 ]] && echo "| $2 | $c |"; }
render(){ # $1 bucket  $2 title  $3 desc  $4 max(0=all)
  local key="$1" title="$2" desc="$3" max="${4:-0}" n; n=$(count_of "$key")
  echo "## $title  ($n)"; echo
  [[ -n "$desc" ]] && { echo "$desc"; echo; }
  if [[ "$n" -eq 0 ]]; then echo "_None._"; echo; return 0; fi
  echo "| ★ | Gene | Variant | Genotype | ClinVar significance | Condition / drug |"
  echo "|---|---|---|---|---|---|"
  awk -F'\t' -v b="$key" '$1==b{printf "%s\t%s\t%s\t%s (%s)\t%s\t%s\n",$2,$5,$6,$7,$8,$9,$10}' "$ALLCATS" \
    | sort -t$'\t' -k1,1nr -k2,2 \
    | { if [[ "$max" -gt 0 ]]; then head -n "$max"; else cat; fi; } \
    | awk -F'\t' '{printf "| %s | %s | %s | %s | %s | %s |\n",$1,$2,$3,$4,$5,$6}'
  [[ "$max" -gt 0 && "$n" -gt "$max" ]] && { echo; echo "_…and $((n-max)) more in \`$(basename "$ALLCATS")\`._"; }
  echo
}
{
  echo "# ClinVar report — every tier — $(sample_name)"
  echo
  echo "_Reference: $(ref_build). Generated by 20_clinvar_annotate.sh. Gitignored._"
  echo
  echo "> ⚠️ **Educational only — not medical advice.** A ClinVar tag is not a"
  echo "> diagnosis. ★ = ClinVar review confidence (0–4); 0★ = a single unreviewed"
  echo "> submission, often a common variant. Heterozygous + a recessive condition"
  echo "> = healthy **carrier**, not affected. Confirm anything actionable with a"
  echo "> clinician or pharmacist (pharmacogenomics: see CPIC / cpicpgx.org)."
  echo
  echo "## Summary — all ClinVar tiers among variants you carry"
  echo
  echo "| Tier | Variants you carry |"
  echo "|---|---|"
  summary_row "1_Pathogenic" "Pathogenic / Likely-pathogenic"
  summary_row "2_DrugResponse" "Drug response (pharmacogenomics)"
  summary_row "3_RiskFactor" "Risk-factor alleles"
  summary_row "4_Protective" "Protective alleles"
  summary_row "5_Association" "Disease-association alleles"
  summary_row "6_Uncertain" "Uncertain significance"
  summary_row "8_Conflicting" "Conflicting classifications"
  summary_row "7_Other" "Other"
  summary_row "9_Benign" "Benign / Likely-benign"
  echo
  render "1_Pathogenic" "Pathogenic / Likely-pathogenic" \
    "Disease variants. Sort by ★ — most low-★ rows are common variants." 0
  render "2_DrugResponse" "Drug response (pharmacogenomics)" \
    "How you may respond to specific drugs. The actionable category — bring high-★ rows for established pharmacogenes (CYP2C19, CYP2C9, VKORC1, CYP2D6, CYP2B6, UGT1A1, NAT2, CYP3A5, DPYD, TPMT) to a pharmacist." 0
  render "4_Protective" "Protective alleles" \
    "Variants associated with reduced risk." 0
  render "3_RiskFactor" "Risk-factor alleles" \
    "Common-variant risk modifiers (small effects, not Mendelian)." 40
  render "5_Association" "Disease-association alleles" \
    "GWAS-style associations; mostly small statistical effects." 40
  echo "## Lower-signal tiers (counts only)"
  echo
  echo "Not enumerated here (mostly common-variant noise) — full rows are in \`$(basename "$ALLCATS")\`:"
  echo
  echo "- Uncertain significance: **$(count_of 6_Uncertain)**"
  echo "- Conflicting classifications: **$(count_of 8_Conflicting)**"
  echo "- Benign / Likely-benign: **$(count_of 9_Benign)**"
  echo "- Other: **$(count_of 7_Other)**"
} > "$REPORT"

# focused TSVs for the actionable tiers
awk -F'\t' '$1=="1_Pathogenic"'   "$ALLCATS" > "$OUT_DIR/20_clinvar_pathogenic.tsv"
awk -F'\t' '$1=="2_DrugResponse"' "$ALLCATS" > "$OUT_DIR/20_clinvar_drug_response.tsv"

ok "Wrote $REPORT"
ok "  Pathogenic/Likely-path: $(count_of 1_Pathogenic) | drug-response: $(count_of 2_DrugResponse) | protective: $(count_of 4_Protective) | risk: $(count_of 3_RiskFactor) | association: $(count_of 5_Association)"
ok "  Full categorized table: $ALLCATS"
warn "★ matters: most 0★ rows are common variants with one unreviewed submission. Not a diagnosis."
