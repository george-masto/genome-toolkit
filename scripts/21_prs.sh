#!/usr/bin/env bash
# ============================================================================
# 21_prs.sh - Polygenic risk scores (PRS) from this genome. GEO-163.
#
# Produces an ancestry-matched risk PERCENTILE for a trait (the lens single-
# variant ClinVar findings can't give) and wires it into the MCP as
# get_polygenic_scores  (analysis/calls/prs.json -> build_insights.py -> tool).
#
# THE CRITICAL LESSON (why this is more than "plink2 --score"):
#   A variant-only VCF CANNOT score PRS correctly. For ~half of a score's
#   variants the effect allele IS the reference allele, so a homozygous-reference
#   genotype = TWO effect alleles. A variant-only file omits hom-ref sites, so the
#   scorer counts them as 0 instead of 2 (this genome lost ~80% of its effect-
#   allele dosage that way: 567k vs a full genome's ~2.7M). You MUST reconstruct
#   dense genotypes at the score positions from the gVCF's hom-ref blocks (filled
#   against the reference FASTA) before scoring.
#
# VALIDATED PIPELINE (2026-06-30, PGS000018 metaGRS_CAD -> EUR pctile (example result)):
#   1. Toolchain: JDK17, nextflow (official installer -> tools/nextflow), Miniforge
#      (~/miniforge3, conda+mamba), plink2 + zstd via conda, bcftools/samtools brew.
#   2. pgsc_calc (conda profile; the docker daemon won't run headless) downloads +
#      harmonizes the PGS Catalog scorefile -> normalised_<PGS>_hmPOS_GRCh38.txt.gz.
#      Build a plink2 scorefile: ID=chr:pos, A1=effect_allele, BETA=effect_weight.
#   3. Reference: HGDP+1kGP panel
#        ftp.ebi.ac.uk/pub/databases/spot/pgs/resources/pgsc_HGDP+1kGP_v1.tar.zst
#      Extract GRCh38_HGDP+1kGP_ALL.{pgen,pvar.zst,psam} + king.cutoff. Score all
#      samples densely:  plink2 --pfile REF vzs --set-all-var-ids @:# --rm-dup
#      exclude-all --score SF 1 2 3 header cols=+scoresums  -> ref_scores.sscore.
#      plink2 ... --freq -> ref.afreq (subset to score positions to stay small).
#   4. Subject DENSE genotypes: download GRCh38 primary-assembly FASTA (Ensembl,
#      contig naming "1" to match the gVCF), then
#        bcftools view -R score_pos GVCF -Ou \
#         | bcftools convert --gvcf2vcf -f FASTA -Ou \
#         | bcftools view -T score_pos -Oz -o tmp ; bcftools sort -> subject.dense.vcf.gz
#      (restrict to score positions BEFORE sort, or whole-block expansion blows up
#      to GBs and OOMs plink2). hom-ref sites are now EXPLICIT 0/0.
#   5. Score subject + RANK: count effect-allele dosage per variant from the dense
#      VCF (python; sidesteps a plink2 single-sample --read-freq segfault), then
#      rank vs the EUR reference subset (relateds via king.cutoff excluded) ->
#      percentile + z. SANITY GATE: subject dosage_sum must be ~2.5-2.8M (a full
#      genome), NOT ~500k (the variant-only artifact). Write analysis/calls/prs.json.
#   6. scripts/50_serve_mcp.sh --rebuild   to push live.
#
# Per-stage working scripts (exactly what ran) are in scripts/_prs/. Large inputs
# (reference panel, FASTA, dense VCF) are gitignored; re-fetch per steps 3-4.
# ============================================================================
set -uo pipefail
cat <<'EOF'
21_prs.sh is the documented canonical PRS pipeline (see header).
Validated: PGS000018 (Inouye metaGRS_CAD) -> EUR percentile (example result).
Result served at analysis/calls/prs.json via the MCP tool get_polygenic_scores.
Per-stage scripts: scripts/_prs/   |   push live: scripts/50_serve_mcp.sh --rebuild
To add a trait: pick its PGS Catalog id, repeat steps 2-5 with that scorefile.
EOF
