#!/usr/bin/env bash
# ============================================================================
# 16_gba.sh - GBA1 genotyping with Illumina Gauchian.
#
# GBA1 sits beside the GBAP1 pseudogene and carries recombinant alleles
# (RecNciI, RecTL, p.L483P=legacy L444P...) that a generic SNP caller mis-types or
# misses entirely. Gauchian is the GBA analogue of Cyrius - paralog-aware. GBA1
# variants cause Gaucher disease (recessive) AND are the single strongest common
# genetic risk factor for Parkinson's disease, so it is worth resolving even
# though it is not a PGx dosing gene.
#
# Output:  analysis/calls/gba.json   (build_insights.py -> carrier_status).
#
#   scripts/10_align_from_fastq.sh   ->   scripts/16_gba.sh   ->   --rebuild
#
# Env: GAUCHIAN (the gauchian command, default: gauchian on PATH), THREADS
# ============================================================================
source "$(dirname "$0")/lib/common.sh"
require_tool samtools
have python3 || die "python3 required."
require_cram
GAUCHIAN="${GAUCHIAN:-gauchian}"
have "$GAUCHIAN" || die "gauchian not found. Install: scripts/00_install_tools.sh --callers
     (pip install gauchian), or set GAUCHIAN=/path/to/gauchian."

THREADS="${THREADS:-$(sysctl -n hw.ncpu 2>/dev/null || nproc)}"
SAMPLE="$(sample_name)"; [[ "$SAMPLE" == "unknown" ]] && SAMPLE="sample"
WORK="$CALLS_DIR/gba.work"; mkdir -p "$WORK"
OUT_JSON="$CALLS_DIR/gba.json"

[[ -f "$REF_FA.fai" ]] || { [[ -f "$REF_FA" ]] && samtools faidx "$REF_FA"; }
export REF_PATH="$REF_FA"
MANIFEST="$WORK/manifest.txt"; printf '%s\n' "$CRAM" > "$MANIFEST"

log "Running Gauchian (GBA1 / GBAP1)..."
"$GAUCHIAN" \
  --manifest "$MANIFEST" --genome 38 \
  --reference "$REF_FA" \
  --prefix "$SAMPLE" --outDir "$WORK" --threads "$THREADS" \
  || die "Gauchian failed. Check $WORK and that the CRAM is GRCh38 + indexed."

# --- normalize: Gauchian writes <prefix>.tsv (one row/sample) ---------------
CY_TSV="$(_first "$WORK/${SAMPLE}.tsv" "$WORK"/*.tsv)"
GENERATED_UTC="$(now_utc)" SAMPLE="$SAMPLE" CY_TSV="$CY_TSV" \
python3 - "$OUT_JSON" <<'PY' || die "Could not parse Gauchian output."
import json, os, sys, csv
out_path=sys.argv[1]; sample=os.environ.get("SAMPLE","sample")
cyt=os.environ.get("CY_TSV","")
def err(note):
    json.dump({"schema":"genome-call/1","tool":"Gauchian","kind":"carrier","gene":"GBA1",
               "status":"error","note":note,"generated_utc":os.environ["GENERATED_UTC"]},
              open(out_path,"w"),indent=2); sys.exit(0)
if not (cyt and os.path.exists(cyt)): err("no Gauchian TSV")
rows=list(csv.DictReader(open(cyt),delimiter="\t"))
if not rows: err("empty Gauchian TSV")
row=next((r for r in rows if r.get("Sample")==sample), rows[0])

def by(sub, default=None):                      # match a column by substring
    for k,v in row.items():
        if sub.lower() in k.lower(): return v
    return default
def by_start(pre, default=None):                # match a column by prefix (exact-ish)
    for k,v in row.items():
        if k.lower().startswith(pre.lower()): return v
    return default
def truthy(v): return str(v).strip().lower() in ("true","1","yes")
def real(v):   return v not in (None,"","None",".","NA")

is_biallelic = by("is_biallelic")
is_carrier   = by("is_carrier")
cn           = by("CN(")              # CN(GBA+GBAP1)
# the STANDALONE recombinant-allele column starts with "GBAP1-like"; the
# is_biallelic(...)/is_carrier(...) columns also CONTAIN that token, so match by
# prefix to avoid grabbing a boolean column.
recomb       = by_start("GBAP1-like")
other        = by("other_unphased")

variants=[v for v in (recomb, other) if real(v)]
no_call = not real(cn) and not variants
if no_call:
    status, result, actionable = "no_call", "GBA1 not callable (no copy-number evidence)", False
else:
    status="ok"
    bits=[]
    if real(recomb): bits.append(f"recombinant/exon9-11: {recomb}")
    if real(other):  bits.append(f"other variant(s): {other}")
    result = ("GBA1: " + "; ".join(bits)) if bits else "GBA1: no pathogenic variant detected"
    actionable = truthy(is_biallelic) or truthy(is_carrier) or bool(variants)

call={"schema":"genome-call/1","tool":"Gauchian","kind":"carrier","gene":"GBA1",
      "status":status,
      "is_biallelic":truthy(is_biallelic),"is_carrier":truthy(is_carrier),
      "combined_cn":cn if real(cn) else None,
      "variants":variants,"result":result,
      "condition":"Gaucher disease (recessive); Parkinson's disease risk factor",
      "actionable":bool(actionable),"drugs":[],
      "implication":("GBA1 variant(s) -> Gaucher carrier/affected and elevated Parkinson's risk. "
                     "Paralog-aware call from WGS; confirm in a clinical lab. Not a PGx dosing gene. "
                     "Legacy names: p.L483P=L444P, p.N409S=N370S."),
      "raw":row,"generated_utc":os.environ["GENERATED_UTC"]}
json.dump(call,open(out_path,"w"),indent=2)
print(f"GBA1: biallelic={call['is_biallelic']} carrier={call['is_carrier']} variants={variants} ({status})")
PY

ok "Wrote $OUT_JSON"
log "Rebuild insights:  scripts/50_serve_mcp.sh --rebuild"
