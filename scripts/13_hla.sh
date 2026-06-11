#!/usr/bin/env bash
# ============================================================================
# 13_hla.sh - HLA typing for drug-hypersensitivity pharmacogenetics, with T1K.
#
# HLA genes have thousands of alleles; a SNP/indel VCF cannot resolve a 2-field
# type. Yet a handful of HLA alleles are among the most actionable PGx results
# there are: HLA-B*57:01 (abacavir - contraindicated), HLA-B*15:02 / HLA-A*31:01
# (carbamazepine SJS/TEN), HLA-B*58:01 (allopurinol). T1K types HLA straight
# from WGS reads; we extract the MHC + unmapped reads, type them, and flag any
# carried risk allele (assets/pharmacogenes/hla_pgx.tsv).
#
# Output:  analysis/calls/hla.json   (build_insights.py -> hla; check_drug()
#          surfaces the drug links; get_hla() lists them).
#
#   scripts/10_align_from_fastq.sh   ->   scripts/13_hla.sh   ->   --rebuild
#
# Env: HLA_IDX (T1K index dir, default tools/hlaidx - built once via t1k-build.pl),
#      THREADS, MHC_REGION (default GRCh38 chr6:28510000-33480000)
# ============================================================================
source "$(dirname "$0")/lib/common.sh"
require_tool samtools
have run-t1k || die "T1K (run-t1k) not found. Install: scripts/00_install_tools.sh --callers
     (conda install -c bioconda t1k), or set PATH to a T1K build."
require_cram
[[ -f "$REF_FA" ]] || die "Need the reference FASTA at $REF_FA for CRAM decoding."

THREADS="${THREADS:-$(sysctl -n hw.ncpu 2>/dev/null || nproc)}"
SAMPLE="$(sample_name)"; [[ "$SAMPLE" == "unknown" ]] && SAMPLE="sample"
WORK="$CALLS_DIR/hla.work"; mkdir -p "$WORK"
OUT_JSON="$CALLS_DIR/hla.json"
HLA_IDX="${HLA_IDX:-$REPO_ROOT/tools/hlaidx}"
PGX_TSV="$ASSETS_DIR/pharmacogenes/hla_pgx.tsv"

# --- 1. one-time IPD-IMGT/HLA index ----------------------------------------
if [[ ! -f "$HLA_IDX/hlaidx_dna_seq.fa" ]]; then
  have t1k-build.pl || die "t1k-build.pl not on PATH (comes with T1K)."
  log "Building the IPD-IMGT/HLA index once into $HLA_IDX (downloads, a few min)..."
  mkdir -p "$(dirname "$HLA_IDX")"
  ( cd "$(dirname "$HLA_IDX")" && t1k-build.pl -o "$(basename "$HLA_IDX")" --download IPD-IMGT/HLA ) \
    || die "t1k-build.pl failed (needs internet for the IMGT/HLA download)."
fi

# --- 2. extract MHC + unmapped reads -> paired FASTQ -----------------------
# chr-prefix aware MHC window; add unmapped reads (HLA reads often misalign).
PREFIX=""; samtools idxstats "$CRAM" 2>/dev/null | cut -f1 | grep -qiE '^chr' && PREFIX="chr"
MHC_REGION="${MHC_REGION:-${PREFIX}6:28510000-33480000}"
log "Extracting MHC ($MHC_REGION) + unmapped reads..."
samtools view -@ "$THREADS" -T "$REF_FA" -u "$CRAM" "$MHC_REGION" > "$WORK/region.bam"
samtools view -@ "$THREADS" -T "$REF_FA" -f 4 -u "$CRAM"          > "$WORK/unmapped.bam"
samtools merge -@ "$THREADS" -f -u "$WORK/merged.bam" "$WORK/region.bam" "$WORK/unmapped.bam"
samtools collate -@ "$THREADS" -Ou "$WORK/merged.bam" \
  | samtools fastq -@ "$THREADS" -1 "$WORK/hla_R1.fq.gz" -2 "$WORK/hla_R2.fq.gz" \
      -0 /dev/null -s /dev/null -n -

# --- 3. T1K with the WGS preset, 2-field resolution ------------------------
log "Typing HLA with T1K..."
run-t1k -1 "$WORK/hla_R1.fq.gz" -2 "$WORK/hla_R2.fq.gz" \
  --preset hla-wgs -f "$HLA_IDX/hlaidx_dna_seq.fa" \
  --alleleDigitUnits 2 --alleleDelimiter : \
  -t "$THREADS" --od "$WORK" -o "$SAMPLE" \
  || die "T1K failed. Check $WORK."

# --- 4. normalize -> analysis/calls/hla.json -------------------------------
GT_TSV="$(_first "$WORK/${SAMPLE}_genotype.tsv" "$WORK"/*_genotype.tsv)"
[[ -f "$GT_TSV" ]] || die "No T1K genotype TSV in $WORK."
GENERATED_UTC="$(now_utc)" GT_TSV="$GT_TSV" PGX_TSV="$PGX_TSV" \
python3 - "$OUT_JSON" <<'PY' || die "Could not parse T1K output."
import json, os, sys
out_path=sys.argv[1]
# T1K genotype tsv: gene, num_diff, allele_1, abund_1, qual_1, allele_2, abund_2, qual_2, secondary
diplo, carried = {}, set()
for line in open(os.environ["GT_TSV"]):
    p=line.rstrip("\n").split("\t")
    if len(p) < 3: continue
    gene=p[0]; calls=[]
    for ai,qi in ((2,4),(5,7)):
        if len(p)>qi:
            allele=p[ai].strip()
            try: q=float(p[qi])
            except (ValueError,IndexError): q=1.0
            if allele and allele not in (".","") and q>0:    # README: drop quality<=0
                calls.append(allele); carried.add(allele)
    if calls: diplo[gene]=calls

# flag PGx risk alleles from the public table
alleles=[]
for line in open(os.environ["PGX_TSV"]):
    if line.startswith("#") or not line.strip(): continue
    p=line.rstrip("\n").split("\t")
    if len(p)<3: continue
    allele,drugs,impl=p[0],p[1],p[2]
    # carried if any typed 2-field token starts with the risk allele (e.g. HLA-B*57:01[:01])
    is_carried=any(c.startswith(allele) or allele.startswith(c) for c in carried)
    alleles.append({"allele":allele,"carried":is_carried,
                    "drugs":[d.strip() for d in drugs.replace(";",",").split(",") if d.strip()],
                    "actionable":is_carried,
                    "implication":impl if is_carried else f"Not detected. ({impl})"})

flagged=[a["allele"] for a in alleles if a["carried"]]
call={"schema":"genome-call/1","tool":"T1K","kind":"hla",
      "status":"ok" if diplo else "no_call",
      "diplotypes":diplo,"alleles":alleles,
      "summary":(f"PGx risk allele(s) present: {', '.join(flagged)}" if flagged
                 else "no catalogued HLA drug-hypersensitivity risk allele detected"),
      "generated_utc":os.environ["GENERATED_UTC"]}
json.dump(call,open(out_path,"w"),indent=2)
print(f"HLA: typed {len(diplo)} genes; PGx flags: {', '.join(flagged) or 'none'}")
PY

ok "Wrote $OUT_JSON"
log "Rebuild insights:  scripts/50_serve_mcp.sh --rebuild"
