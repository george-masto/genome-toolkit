#!/usr/bin/env bash
# ============================================================================
# 02_index_vcfs.sh — ensure every VCF is BGZF-compressed and tabix-indexed so
# downstream tools (region queries, rsID lookups) are fast.
#
# sequencing.com ships BGZF already, so this usually just writes .tbi sidecars
# next to the originals. If a file is plain gzip, it is re-compressed as BGZF
# into analysis/ (originals are never modified).
# ============================================================================
source "$(dirname "$0")/lib/common.sh"
require_tool bcftools
require_tool tabix

index_one() {
  local vcf="$1"
  [[ -n "$vcf" && -f "$vcf" ]] || return 0

  if [[ "$(htsfile "$vcf" 2>/dev/null)" != *BGZF* ]]; then
    warn "$(basename "$vcf") is not BGZF — re-compressing into $OUT_DIR"
    local rebg="$OUT_DIR/$(basename "${vcf%.gz}").vcf.gz"
    gzcat "$vcf" | bgzip > "$rebg"
    vcf="$rebg"
  fi

  if [[ -f "$vcf.tbi" && "$vcf.tbi" -nt "$vcf" ]]; then
    ok "$(basename "$vcf") already indexed"
  else
    log "Indexing $(basename "$vcf")…"
    tabix -f -p vcf "$vcf"
    ok "wrote $(basename "$vcf").tbi"
  fi
}

show_inputs
index_one "$VCF_SNP"
index_one "$VCF_CNV"
index_one "$VCF_SV"
log "All VCFs indexed."
