#!/usr/bin/env bash
# ============================================================================
# 07_pgx_report.sh — a clean, pharmacist-ready pharmacogenomic profile card.
#
# Screens the curated CPIC pharmacogene catalog (assets/pharmacogenes/) against
# your genome, predicts a metabolizer/function phenotype per gene, and prints
# the drug implications — WITHOUT the low-evidence research-panel noise that
# clutters the raw ClinVar drug_response tier.
#
# Output: analysis/07_pgx_report.md  (gitignored — your genotype).
#
# ⚠️  Predicted phenotypes from a single-sample short-read VCF. CYP2D6 cannot
#     be fully typed here. EDUCATIONAL ONLY — confirm against CPIC (cpicpgx.org)
#     with a pharmacist before changing anything.
#
# Run 03_extract_variants.sh first.
# ============================================================================
source "$(dirname "$0")/lib/common.sh"
require_tool bcftools

VARIANTS="${VARIANTS:-$OUT_DIR/variants.vcf.gz}"
CATALOG="${PGX_CATALOG:-$ASSETS_DIR/pharmacogenes/pgx_variants.tsv}"
REPORT="$OUT_DIR/07_pgx_report.md"
[[ -f "$VARIANTS" ]] || die "Missing $VARIANTS — run scripts/03_extract_variants.sh first."
[[ -f "$CATALOG"  ]] || die "Missing PGx catalog $CATALOG"

rsfile="$(mktemp)"; foundfile="$(mktemp)"
trap 'rm -f "$rsfile" "$foundfile"' EXIT
grep -v '^#' "$CATALOG" | awk -F'\t' 'NF>=7{print $2}' | sort -u > "$rsfile"
log "Screening $(wc -l < "$rsfile" | tr -d ' ') pharmacogene variants…"
bcftools view "$VARIANTS" 2>/dev/null \
  | bcftools query -i "ID=@$rsfile" -f '%ID\t[%GT]\n' 2>/dev/null > "$foundfile" || true
ok "$(wc -l < "$foundfile" | tr -d ' ') present (rest = homozygous reference = *1/normal)"

awk -F'\t' -v sample="$(sample_name)" '
  function star_diplo(gene,   k,arr,i,rs,c,nm,cnt,out) {
    k=split(geners[gene],arr," "); cnt=0; out=""
    for(i=1;i<=k;i++){ rs=arr[i]; if(rs=="")continue; c=copies[gene,rs]+0
      while(c-->0){ nm=defn[gene,rs]; sub(/^[A-Z0-9]+/,"",nm); out=out (out==""?"":"/") nm; cnt++ } }
    while(cnt<2){ out="*1" (out==""?"":"/") out; cnt++ }   # pad *1 at front (*1/*2)
    return out
  }
  function metab(act){
    if(act<=0.5)return "Poor metabolizer"
    if(act<=1.5)return "Intermediate metabolizer"
    if(act<=2.25)return "Normal metabolizer"
    if(act<=2.75)return "Rapid metabolizer"
    return "Ultrarapid metabolizer"
  }
  # ---- file 1: found genotypes (rsid GT) ----
  FNR==NR { gtv[$1]=$2; next }
  /^#/ || NF<7 { next }
  {
    gene=$1; rsid=$2; def=$3; fn=$4; model=$5; drugs=$6; action=$7
    if(!(gene in seen)){ order[++ng]=gene; seen[gene]=1; gmodel[gene]=model; gdrugs[gene]=drugs; gaction[gene]=action }
    c=0; if(rsid in gtv){ n=split(gtv[rsid],a,/[\/|]/); for(i=1;i<=n;i++) if(a[i]=="1") c++ }
    f=1.0; if(fn=="no_function")f=0.0; else if(fn=="decreased")f=0.5; else if(fn=="increased")f=1.5
    actadj[gene]+=c*(f-1.0)
    copies[gene,rsid]=c; defn[gene,rsid]=def; geners[gene]=geners[gene] rsid " "
    if(c>0){ ncar[gene]++ }
    if(model=="nat2"   && c>0) slow[gene]+=c
    if(model=="g6pd"   && c>0) g6[gene]+=c
    if(model=="vkorc1") vk[gene]=c
    if(model=="cyp4f2") c4[gene]=c
    if(model=="ugt")    ug[gene]=c
    if(model=="cyp3a5") c35[gene]=c
    if(model=="slco")   sl[gene]=c
    if(model=="ifnl3")  il[gene]=c
    if(model=="cyp2d6_partial" && c>0){ nm=def; sub(/^[A-Z0-9]+/,"",nm); d6=d6 (d6==""?"":", ") nm }
  }
  END {
    OFS="\t"
    for(gi=1; gi<=ng; gi++){
      g=order[gi]; m=gmodel[g]; act=2.0+actadj[g]
      diplo=""; ph=""; neutral=0
      if(m=="score"){ diplo=star_diplo(g); ph=metab(act); neutral=(ph=="Normal metabolizer") }
      else if(m=="lof"){ diplo=star_diplo(g); ph=metab(act); if(ph=="Rapid metabolizer"||ph=="Ultrarapid metabolizer")ph="Normal metabolizer"; neutral=(ph=="Normal metabolizer") }
      else if(m=="vkorc1"){ c=vk[g]+0; diplo=(c==0?"-1639 GG":c==1?"-1639 GA":"-1639 AA"); ph=(c==0?"Normal warfarin sensitivity":c==1?"Increased sensitivity (lower dose)":"High sensitivity (lowest dose)"); neutral=(c==0) }
      else if(m=="cyp4f2"){ c=c4[g]+0; diplo=(c==0?"*1/*1":c==1?"*1/*3":"*3/*3"); ph=(c==0?"Normal":"Carries *3 (minor: slightly higher warfarin dose)"); neutral=(c==0) }
      else if(m=="slco"){ c=sl[g]+0; diplo=(c==0?"*1/*1":c==1?"*1/*5":"*5/*5"); ph=(c==0?"Normal function":c==1?"Decreased function (higher statin myopathy risk)":"Poor function (high statin myopathy risk)"); neutral=(c==0) }
      else if(m=="ugt"){ c=ug[g]+0; diplo=(c==0?"*1/*1":c==1?"*1/*28":"*28/*28"); ph=(c==0?"Normal":c==1?"Intermediate (Gilbert\47s carrier)":"Reduced (Gilbert\47s; irinotecan caution)"); neutral=(c==0) }
      else if(m=="cyp3a5"){ c=c35[g]+0; diplo=(c==2?"*3/*3":c==1?"*1/*3":"*1/*1"); ph=(c==2?"Non-expressor (standard tacrolimus dose)":c==1?"Intermediate expressor":"Expressor (higher tacrolimus dose)"); neutral=(c==2) }
      else if(m=="nat2"){ s=slow[g]+0; diplo=sprintf("%d slow allele(s)",s); ph=(s>=2?"Slow acetylator":s==1?"Intermediate acetylator":"Rapid/normal acetylator"); neutral=(s==0) }
      else if(m=="ifnl3"){ c=il[g]+0; diplo=(c==0?"CC":c==1?"CT":"TT"); ph=(c==0?"Favorable (CC)":c==1?"Intermediate (CT)":"Unfavorable (TT)"); neutral=(c==0) }
      else if(m=="g6pd"){ c=g6[g]+0; diplo=(c==0?"none detected":"variant present"); ph=(c==0?"No common deficiency variant (normal)":"Possible deficiency — confirm by enzyme assay"); neutral=(c==0) }
      else if(m=="cyp2d6_partial"){ diplo=(d6==""?"no *4/*10/*41":"carries "d6); ph="Incompletely typed (copy-number not callable)"; neutral=2 }
      else { diplo=star_diplo(g); ph="see details"; neutral=0 }
      print g, diplo, ph, neutral, gdrugs[g], gaction[g]
    }
  }
' "$foundfile" "$CATALOG" > "$foundfile.rows"

# ---- assemble markdown ----
{
  echo "# Pharmacogenomic profile — $(sample_name)"
  echo
  echo "_Generated by 07_pgx_report.sh. Gitignored — your genotype._"
  echo
  echo "> ⚠️ **Educational, predicted from WGS — not a clinical PGx test.** Confirm"
  echo "> against CPIC (cpicpgx.org) with a pharmacist before any decision. CYP2D6"
  echo "> is only partially typeable from this data."
  echo
  echo "## Profile card"
  echo
  echo "| Gene | Genotype | Predicted phenotype | Key drugs | Implication |"
  echo "|---|---|---|---|---|"
  awk -F'\t' '{printf "| %s | %s | %s | %s | %s |\n",$1,$2,$3,$5,$6}' "$foundfile.rows"
  echo
  echo "## Bring these to a prescriber / pharmacist (non-default phenotypes)"
  echo
  awk -F'\t' '$4==0{printf "- **%s — %s** (%s). Drugs: %s\n",$1,$3,$2,$5}' "$foundfile.rows"
  echo
  echo "## Normal / no risk allele detected (reassuring)"
  echo
  awk -F'\t' '$4==1{a=a (a==""?"":", ") $1} END{if(a!="")print "Standard dosing expected for: " a "."}' "$foundfile.rows"
  echo
  echo "## Not fully typeable from this data"
  echo
  awk -F'\t' '$4==2{printf "- **%s** — %s. %s\n",$1,$3,$6}' "$foundfile.rows"
} > "$REPORT"

ok "Wrote $REPORT"
hr
cat "$REPORT" >&2
