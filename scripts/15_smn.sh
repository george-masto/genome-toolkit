#!/usr/bin/env bash
# ============================================================================
# 15_smn.sh - SMN1/SMN2 copy number with Illumina SMNCopyNumberCaller.
#
# Spinal muscular atrophy is a copy-number disease of SMN1, which is ~identical
# to its paralog SMN2 - so it is invisible to a SNP VCF and needs a paralog/CN
# aware caller. SMN1 copy number: 0 = affected, 1 = CARRIER, 2 = typical. This
# is one of the most common serious recessive carrier states; carrier screening
# is standard of care preconception.
#
# Output:  analysis/calls/smn.json   (build_insights.py -> carrier_status,
#          surfaced by the MCP's get_carrier_status tool).
#
#   scripts/10_align_from_fastq.sh   ->   scripts/15_smn.sh   ->   --rebuild
#
# Env: SMN_DIR (checkout dir, default tools/SMNCopyNumberCaller, auto-cloned),
#      THREADS, NO_CLONE=1
# ============================================================================
source "$(dirname "$0")/lib/common.sh"
require_tool samtools; require_tool git
have python3 || die "python3 required."
require_cram

THREADS="${THREADS:-$(sysctl -n hw.ncpu 2>/dev/null || nproc)}"
SAMPLE="$(sample_name)"; [[ "$SAMPLE" == "unknown" ]] && SAMPLE="sample"
SMN_DIR="${SMN_DIR:-$REPO_ROOT/tools/SMNCopyNumberCaller}"
WORK="$CALLS_DIR/smn.work"; mkdir -p "$WORK"
OUT_JSON="$CALLS_DIR/smn.json"

if [[ ! -f "$SMN_DIR/smn_caller.py" ]]; then
  [[ "${NO_CLONE:-0}" == "1" ]] && die "SMNCopyNumberCaller not found at $SMN_DIR (NO_CLONE set)."
  log "Cloning Illumina SMNCopyNumberCaller..."
  git clone --depth 1 https://github.com/Illumina/SMNCopyNumberCaller.git "$SMN_DIR" \
    || die "Could not clone SMNCopyNumberCaller; set SMN_DIR to a checkout."
fi

[[ -f "$REF_FA.fai" ]] || { [[ -f "$REF_FA" ]] && samtools faidx "$REF_FA"; }
export REF_PATH="$REF_FA"
MANIFEST="$WORK/manifest.txt"; printf '%s\n' "$CRAM" > "$MANIFEST"

log "Running SMNCopyNumberCaller (SMN1/SMN2 copy number)..."
python3 "$SMN_DIR/smn_caller.py" \
  --manifest "$MANIFEST" --genome 38 \
  --prefix "$SAMPLE" --outDir "$WORK" --threads "$THREADS" \
  || die "SMNCopyNumberCaller failed. Check $WORK and that the CRAM is GRCh38."

CY_JSON="$(_first "$WORK/${SAMPLE}.json" "$WORK"/*.json)"
CY_TSV="$(_first "$WORK/${SAMPLE}.tsv" "$WORK"/*.tsv)"
GENERATED_UTC="$(now_utc)" SAMPLE="$SAMPLE" CY_JSON="$CY_JSON" CY_TSV="$CY_TSV" \
python3 - "$OUT_JSON" <<'PY' || die "Could not parse SMN output."
import json, os, sys, csv
out_path = sys.argv[1]; sample = os.environ.get("SAMPLE", "sample")
cyj, cyt = os.environ.get("CY_JSON",""), os.environ.get("CY_TSV","")
def pick(d,*names,default=None):
    low={k.lower():v for k,v in d.items()}
    for n in names:
        if n.lower() in low: return low[n.lower()]
    return default
def to_int(x):
    try: return int(float(x)) if x not in (None,"","None",".") else None
    except (TypeError,ValueError): return None

rec=None
if cyj and os.path.exists(cyj):
    data=json.load(open(cyj))
    rec=data.get(sample) or (next(iter(data.values())) if isinstance(data,dict) and data else None)
if rec is None and cyt and os.path.exists(cyt):
    rows=list(csv.DictReader(open(cyt),delimiter="\t"))
    if rows: rec=next((r for r in rows if r.get("Sample")==sample),rows[0])

if not isinstance(rec,dict):
    json.dump({"schema":"genome-call/1","tool":"SMNCopyNumberCaller","kind":"carrier",
               "gene":"SMN1","status":"error","note":"no parseable record",
               "generated_utc":os.environ["GENERATED_UTC"]},open(out_path,"w"),indent=2); sys.exit(0)

smn1=to_int(pick(rec,"SMN1_CN","SMN1","Total_SMN1_CN"))
smn2=to_int(pick(rec,"SMN2_CN","SMN2"))
carrier_field=pick(rec,"isCarrier","Carrier")
is_carrier = (str(carrier_field).lower() in ("true","1","yes")) if carrier_field is not None else (smn1==1)
if smn1 is None: result="SMN1 copy number indeterminate"
elif smn1==0: result="SMN1 copy number 0 - consistent with SMA (affected); confirm clinically"
elif smn1==1: result="SMN1 copy number 1 - SMA CARRIER"
else: result=f"SMN1 copy number {smn1} - typical (not a carrier by CN)"
status="ok" if smn1 is not None else "no_call"
actionable = bool(smn1 is not None and smn1<=1)
call={"schema":"genome-call/1","tool":"SMNCopyNumberCaller","kind":"carrier","gene":"SMN1",
      "status":status,"smn1_cn":smn1,"smn2_cn":smn2,"is_carrier":bool(is_carrier),
      "result":result,"condition":"Spinal muscular atrophy (recessive)","actionable":actionable,
      "drugs":[],"implication":("Carrier/affected by copy number - confirm in a clinical lab; "
              "relevant to reproductive carrier screening. Note: the 2+0 'silent carrier' "
              "configuration is not detectable by copy number alone."),
      "raw":rec,"generated_utc":os.environ["GENERATED_UTC"]}
json.dump(call,open(out_path,"w"),indent=2)
print(f"SMN1_CN={smn1} SMN2_CN={smn2} carrier={is_carrier} ({status})")
PY

ok "Wrote $OUT_JSON"
log "Rebuild insights:  scripts/50_serve_mcp.sh --rebuild"
