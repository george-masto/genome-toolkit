#!/usr/bin/env bash
# ============================================================================
# 05_genotype_lookup.sh — look up your genotype at a SNP, position, or region.
#
# Usage:
#   scripts/05_genotype_lookup.sh rs429358            # by dbSNP rsID
#   scripts/05_genotype_lookup.sh 19:44908684         # single position
#   scripts/05_genotype_lookup.sh 19:44905000-44910000  # a region
#   scripts/05_genotype_lookup.sh rs429358 rs7412 rs1801133   # several at once
#
# Region/position queries use the tabix index (instant). rsID queries scan the
# variants file (a few seconds). Reports genotype, zygosity, depth and quality.
# Run 03_extract_variants.sh first.
# ============================================================================
source "$(dirname "$0")/lib/common.sh"
require_tool bcftools

VARIANTS="${VARIANTS:-$OUT_DIR/variants.vcf.gz}"
[[ -f "$VARIANTS" ]] || die "Missing $VARIANTS — run scripts/03_extract_variants.sh first."
[[ $# -ge 1 ]] || die "Give one or more rsIDs / positions / regions. See header for usage."

CP="$(chrom_prefix)"

# Normalize a region arg so the user can type "chr19:.." or "19:.." either way.
norm_region() {
  local r="$1"
  r="${r#chr}"                 # strip any chr prefix the user typed
  [[ -n "$CP" ]] && r="${CP}${r}"   # add back the build's actual prefix
  printf '%s' "$r"
}

zygosity() {
  case "$1" in
    0/0|0\|0) echo "homozygous reference" ;;
    0/1|1/0|0\|1|1\|0) echo "heterozygous" ;;
    1/1|1\|1) echo "homozygous alternate" ;;
    1/2|2/1|*\|*) echo "compound / phased — see GT" ;;
    ./.) echo "no call" ;;
    *) echo "$1" ;;
  esac
}

# bcftools query format: tab-separated fields we print for each hit.
FMT='%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%GT]\t[%DP]\t[%GQ]\t[%AD]\n'

print_hits() {  # reads VCF lines on stdin via bcftools query
  while IFS=$'\t' read -r chrom pos id ref alt gt dp gq ad; do
    [[ -z "$chrom" ]] && continue
    printf '  %-5s %-10s %-12s %s>%s\n' "$chrom" "$pos" "${id:-.}" "$ref" "$alt"
    printf '        genotype: %s  (%s)\n' "$gt" "$(zygosity "$gt")"
    printf '        depth: %s   GQ: %s   allelic depths: %s\n\n' "${dp:-.}" "${gq:-.}" "${ad:-.}"
  done
}

lookup_one() {
  local q="$1"
  hr
  if [[ "$q" == rs* ]]; then
    log "rsID $q — scanning variants…"
    bcftools view "$VARIANTS" 2>/dev/null \
      | bcftools query -i "ID==\"$q\"" -f "$FMT" 2>/dev/null | print_hits \
      || true
  else
    local region; region="$(norm_region "$q")"
    log "region $region — tabix query…"
    bcftools query -r "$region" -f "$FMT" "$VARIANTS" 2>/dev/null | print_hits || true
  fi
}

found_any=0
for q in "$@"; do
  out="$(lookup_one "$q")"
  printf '%s\n' "$out"
  [[ -n "$(printf '%s' "$out" | grep -E 'genotype:')" ]] && found_any=1
done
hr
[[ $found_any -eq 1 ]] || warn "No variant records matched. (Absent rsID often means homozygous reference — i.e. the common genotype — since gVCF hom-ref blocks were stripped.)"
