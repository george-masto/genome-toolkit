#!/usr/bin/env bash
# ============================================================================
# 17_mtdna.sh - mitochondrial PGx + haplogroup from the chrM reads.
#
# Three things a diploid SNP VCF handles poorly on the mitochondrion:
#   (1) MT-RNR1 m.1555A>G / m.1494C>T - CPIC-actionable: these confer risk of
#       aminoglycoside-induced hearing loss (avoid gentamicin/amikacin/etc.).
#   (2) mtDNA haplogroup (maternal lineage) - Haplogrep 3, if installed.
#   (3) heteroplasmy / contamination QC - Haplocheck, if installed.
# mtDNA is haploid and multi-copy, so we call it with --ploidy 1 and read allele
# depths to catch heteroplasmy rather than collapsing to a 0/1 call.
#
# Output:  analysis/calls/mito.json   (build_insights.py -> mito; check_drug()
#          surfaces the aminoglycoside link).
#
#   scripts/10_align_from_fastq.sh   ->   scripts/17_mtdna.sh   ->   --rebuild
#
# Env: HAPLOGREP (haplogrep3 binary; optional), HAPLOCHECK (optional)
# ============================================================================
source "$(dirname "$0")/lib/common.sh"
require_tool samtools; require_tool bcftools
have python3 || die "python3 required."
require_cram
[[ -f "$REF_FA" ]] || die "Need the reference FASTA at $REF_FA for chrM mpileup."

SAMPLE="$(sample_name)"; [[ "$SAMPLE" == "unknown" ]] && SAMPLE="sample"
WORK="$CALLS_DIR/mito.work"; mkdir -p "$WORK"
OUT_JSON="$CALLS_DIR/mito.json"

# --- detect the mito contig name (GRCh38 no-alt uses chrM) ------------------
MT="$(samtools idxstats "$CRAM" 2>/dev/null | cut -f1 | grep -iE '^(chrM|chrMT|MT|M)$' | head -1)"
[[ -n "$MT" ]] || die "No mitochondrial contig found in $CRAM."
log "Mitochondrial contig: $MT"

# --- 1. MT-RNR1 PGx sites (rCRS 1494, 1555); haploid; keep allele depths -----
PGX_TSV="$WORK/mtrnr1.tsv"
bcftools mpileup -f "$REF_FA" -r "$MT:1480-1570" -a FORMAT/AD,FORMAT/DP \
    -d 50000 -q 20 -Q 20 "$CRAM" 2>/dev/null \
  | bcftools call -m --ploidy 1 2>/dev/null \
  | bcftools query -f '%POS\t%REF\t%ALT[\t%GT\t%AD\t%DP]\n' > "$PGX_TSV" || true

# --- 2. full chrM VCF for haplogroup ---------------------------------------
CHRM_VCF="$WORK/${SAMPLE}.chrM.vcf.gz"
bcftools mpileup -f "$REF_FA" -r "$MT" -a FORMAT/AD,FORMAT/DP -d 50000 -q 20 -Q 20 "$CRAM" 2>/dev/null \
  | bcftools call -m --ploidy 1 -Oz -o "$CHRM_VCF" 2>/dev/null && bcftools index -t "$CHRM_VCF" 2>/dev/null || true

# --- 3. haplogroup (optional) ----------------------------------------------
HAPLO_OUT="$WORK/${SAMPLE}.haplogroup.txt"
HG="${HAPLOGREP:-haplogrep3}"
if have "$HG" && [[ -f "$CHRM_VCF" ]]; then
  log "Classifying haplogroup with Haplogrep 3..."
  "$HG" classify --tree phylotree-rcrs@17.0 --input "$CHRM_VCF" --output "$HAPLO_OUT" --extend-report \
    2>/dev/null || warn "haplogrep3 failed; haplogroup will be omitted."
else
  warn "haplogrep3 not found - skipping haplogroup (install: get.haplogrep.org). MT-RNR1 PGx still computed."
fi

# --- 4. contamination/heteroplasmy (optional) ------------------------------
HC_DIR="$WORK/haplocheck"
if have "${HAPLOCHECK:-haplocheck}" && [[ -f "$CHRM_VCF" ]]; then
  log "Estimating contamination with Haplocheck..."
  mkdir -p "$HC_DIR"
  "${HAPLOCHECK:-haplocheck}" --out "$HC_DIR/contamination.txt" "$CHRM_VCF" 2>/dev/null || warn "haplocheck failed; omitting."
fi

# --- 5. normalize -> analysis/calls/mito.json ------------------------------
GENERATED_UTC="$(now_utc)" PGX_TSV="$PGX_TSV" HAPLO_OUT="$HAPLO_OUT" HC_TXT="$HC_DIR/contamination.txt" \
python3 - "$OUT_JSON" <<'PY' || die "Could not assemble mito call."
import json, os, sys, csv
out_path=sys.argv[1]

# --- MT-RNR1 sites ---
SITES={1555:("A","G","m.1555A>G"), 1494:("C","T","m.1494C>T")}
found={}
pgx=os.environ.get("PGX_TSV","")
if pgx and os.path.exists(pgx):
    for line in open(pgx):
        p=line.rstrip("\n").split("\t")
        if len(p)<3: continue
        try: pos=int(p[0])
        except ValueError: continue
        if pos not in SITES: continue
        ref,alt,gt = p[1], p[2], (p[3] if len(p)>3 else ".")
        ad = (p[4] if len(p)>4 else ".")
        dp = (p[5] if len(p)>5 else ".")
        found[pos]=(ref,alt,gt,ad,dp)

def site_call(pos):
    exp_ref,exp_alt,name=SITES[pos]
    drugs=["aminoglycosides (gentamicin, tobramycin, amikacin, streptomycin, etc.)"]
    impl_base=("MT-RNR1 variant associated with aminoglycoside-induced hearing loss (CPIC). "
               "If present, avoid aminoglycosides unless no alternative (then monitor closely).")
    if pos not in found:
        return {"variant":name,"position":pos,"present":False,"vaf":None,"depth":None,
                "drugs":drugs,"actionable":False,
                "implication":f"Reference at m.{pos} (no coverage record) - no increased risk indicated. "+impl_base}
    ref,alt,gt,ad,dp = found[pos]
    alts = alt.split(",") if alt not in (".","") else []
    ai = alts.index(exp_alt) if exp_alt in alts else -1   # AD is [ref, alt1, alt2...]
    vaf=None
    try:
        a=[int(x) for x in ad.split(",")]
        if ai>=0 and len(a) >= ai+2 and sum(a)>0: vaf=round(a[ai+1]/sum(a),3)
    except (ValueError,ZeroDivisionError): vaf=None
    has_alt = (ai>=0) and gt not in ("0",".","0/0")
    try: depth=int(dp)
    except (TypeError,ValueError): depth=None
    present = bool(has_alt and (vaf is None or vaf>=0.03))
    note=""
    if present and vaf is not None and vaf<0.95: note=f" HETEROPLASMIC (VAF {vaf})."
    return {"variant":name,"position":pos,"present":present,"vaf":vaf,"depth":depth,
            "drugs":drugs,"actionable":present,
            "implication":(impl_base+note) if present else f"Reference/absent at m.{pos}. "+impl_base}

mt_rnr1=[site_call(1555), site_call(1494)]

# --- haplogroup ---
haplogroup=None; hg_quality=None
hp=os.environ.get("HAPLO_OUT","")
if hp and os.path.exists(hp):
    rows=list(csv.DictReader(open(hp),delimiter="\t"))
    if rows:
        r=rows[0]
        haplogroup=r.get("Haplogroup") or r.get("haplogroup")
        try: hg_quality=float(r.get("Quality") or r.get("quality"))
        except (TypeError,ValueError): hg_quality=None

# --- contamination ---
contamination=None
hc=os.environ.get("HC_TXT","")
if hc and os.path.exists(hc):
    rows=list(csv.DictReader(open(hc),delimiter="\t"))
    if rows:
        r=rows[0]
        contamination={"status":r.get("Contamination Status") or r.get("Contamination"),
                       "level":r.get("Contamination Level")}

actionable=any(s["actionable"] for s in mt_rnr1)
call={"schema":"genome-call/1","tool":"bcftools + Haplogrep3","kind":"mito",
      "status":"ok" if found or haplogroup else "no_call",
      "mt_rnr1":mt_rnr1,"haplogroup":haplogroup,"haplogroup_quality":hg_quality,
      "contamination":contamination,"actionable":actionable,
      "implication":("MT-RNR1 aminoglycoside-ototoxicity variant detected - see mt_rnr1." if actionable
                     else "No MT-RNR1 aminoglycoside-risk variant detected."),
      "generated_utc":os.environ["GENERATED_UTC"]}
json.dump(call,open(out_path,"w"),indent=2)
flag=[s["variant"] for s in mt_rnr1 if s["present"]]
print(f"mtDNA: haplogroup={haplogroup} MT-RNR1 flags={flag or 'none'}")
PY

ok "Wrote $OUT_JSON"
log "Rebuild insights:  scripts/50_serve_mcp.sh --rebuild"
