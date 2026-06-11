#!/usr/bin/env bash
# ============================================================================
# run_all.sh - quick-start: everything that runs in minutes from the VCFs.
# (Does NOT re-align FASTQ - that's the separate, hours-long 10_align_from_fastq.sh.)
# ============================================================================
source "$(dirname "$0")/lib/common.sh"
here="$(dirname "$0")"

log "1/6  inspect inputs";        "$here/01_inspect.sh"
log "2/6  index VCFs";            "$here/02_index_vcfs.sh"
log "3/6  extract variants";      "$here/03_extract_variants.sh"
log "4/6  QC stats";              "$here/04_vcf_stats.sh"
log "5/6  trait report";          "$here/06_trait_report.sh"
log "6/6  pharmacogenomics";      "$here/07_pgx_report.sh"

hr
ok "Quick-start complete. Outputs are in: $OUT_DIR  (gitignored)"
log "Explore individual variants:  scripts/05_genotype_lookup.sh rs429358"
log "Optional ClinVar screen:      scripts/20_clinvar_annotate.sh"
log "Re-align from raw reads:      scripts/00_install_tools.sh --align && scripts/10_align_from_fastq.sh"
log "Second-pass specialty callers: docs/SECOND_PASS.md  (CYP2D6, HLA, repeats, SMN, GBA, mtDNA - need the CRAM)"
