#!/usr/bin/env bash
# ============================================================================
# 11_cyp2d6.sh - genotype CYP2D6 from the aligned reads with Illumina Cyrius.
#
# This is the one gene a flat short-read VCF genuinely cannot type: the CYP2D7
# pseudogene plus copy-number, hybrid, and duplication alleles defeat a generic
# SNP caller. Cyrius is purpose-built for exactly this. It reads the CRAM from
# 10_align_from_fastq.sh and emits a full diplotype + copy number.
#
# Output:  analysis/calls/cyp2d6.json   (normalized - build_insights.py ingests
#          it and translates the diplotype to a CPIC activity-score phenotype,
#          flipping the MCP's CYP2D6 entry from "incompletely typed" to a real
#          predicted metabolizer status).
#
#   FROM RAW READS:  scripts/10_align_from_fastq.sh   (produces the CRAM)
#   THEN:            scripts/11_cyp2d6.sh
#   THEN:            scripts/50_serve_mcp.sh --rebuild  (MCP picks up the call)
#
# Env:
#   CYRIUS_DIR   path to a Cyrius checkout (default: tools/Cyrius, auto-cloned)
#   THREADS      worker threads (default: all cores)
#   NO_CLONE=1   fail instead of auto-cloning Cyrius
# ============================================================================
source "$(dirname "$0")/lib/common.sh"
require_tool samtools
require_tool git
have python3 || die "python3 required for Cyrius."
require_cram

THREADS="${THREADS:-$(sysctl -n hw.ncpu 2>/dev/null || nproc)}"
SAMPLE="$(sample_name)"; [[ "$SAMPLE" == "unknown" ]] && SAMPLE="sample"
CYRIUS_DIR="${CYRIUS_DIR:-$REPO_ROOT/tools/Cyrius}"
WORK="$CALLS_DIR/cyp2d6.work"; mkdir -p "$WORK"
OUT_JSON="$CALLS_DIR/cyp2d6.json"

# --- 1. obtain Cyrius -------------------------------------------------------
if [[ ! -f "$CYRIUS_DIR/star_caller.py" ]]; then
  [[ "${NO_CLONE:-0}" == "1" ]] && die "Cyrius not found at $CYRIUS_DIR (NO_CLONE set)."
  log "Cloning Illumina Cyrius into $CYRIUS_DIR..."
  mkdir -p "$(dirname "$CYRIUS_DIR")"
  git clone --depth 1 https://github.com/Illumina/Cyrius.git "$CYRIUS_DIR" \
    || die "Could not clone Cyrius. Clone it manually and set CYRIUS_DIR."
fi
# Cyrius deps: numpy, scipy, pysam, statsmodels. Install into the user site if
# missing (non-fatal hint if pip is unavailable).
python3 - <<'PY' 2>/dev/null || warn "Could not verify Cyrius deps; if it fails, pip install numpy scipy pysam statsmodels"
import importlib, sys
missing = [m for m in ("numpy","scipy","pysam","statsmodels") if not importlib.util.find_spec(m)]
sys.exit(1 if missing else 0)
PY

# --- 2. reference for CRAM decoding ----------------------------------------
# pysam decodes CRAM against the reference; make sure it is findable. Pass
# Cyrius an explicit --reference if this version supports it, and also export
# REF_PATH as a belt-and-suspenders fallback for htslib.
[[ -f "$REF_FA" ]] || warn "Reference $REF_FA missing - CRAM decode may fail."
[[ -f "$REF_FA.fai" ]] || { [[ -f "$REF_FA" ]] && samtools faidx "$REF_FA"; }
export REF_PATH="$REF_FA"
REF_FLAG=()
if python3 "$CYRIUS_DIR/star_caller.py" --help 2>&1 | grep -qiE -- '--reference|-r,'; then
  REF_FLAG=(--reference "$REF_FA")
fi

# --- 3. run Cyrius ----------------------------------------------------------
MANIFEST="$WORK/manifest.txt"; printf '%s\n' "$CRAM" > "$MANIFEST"
log "Running Cyrius on $(basename "$CRAM") (CYP2D6; a few minutes)..."
python3 "$CYRIUS_DIR/star_caller.py" \
  --manifest "$MANIFEST" \
  --genome 38 \
  --prefix "$SAMPLE" \
  --outDir "$WORK" \
  --threads "$THREADS" \
  "${REF_FLAG[@]+${REF_FLAG[@]}}" \
  || die "Cyrius failed. Check $WORK and that the CRAM is GRCh38 + indexed."

# --- 4. normalize Cyrius output -> analysis/calls/cyp2d6.json ---------------
# Per Illumina Cyrius: the diplotype + filter live in the TSV (3 columns:
# Sample / Genotype / Filter); the JSON holds QC metrics (Total_CN, Spacer_CN,
# CNV_group, Coverage_MAD). So read Genotype/Filter from the TSV and Total_CN
# from the JSON. Genotype "None" = no-call; an ambiguous call is two diplotypes
# joined by ';'; a Not_assigned filter joins alleles with '_' (kept verbatim,
# scored as indeterminate downstream).
CY_JSON="$(_first "$WORK/${SAMPLE}.json" "$WORK"/*.json)"
CY_TSV="$(_first "$WORK/${SAMPLE}.tsv" "$WORK"/*.tsv)"

GENERATED_UTC="$(now_utc)" SAMPLE="$SAMPLE" CY_JSON="$CY_JSON" CY_TSV="$CY_TSV" \
python3 - "$OUT_JSON" <<'PY' || die "Could not parse Cyrius output."
import json, os, sys, csv

out_path = sys.argv[1]
sample = os.environ.get("SAMPLE", "sample")
cyj, cyt = os.environ.get("CY_JSON", ""), os.environ.get("CY_TSV", "")

def err(note):
    json.dump({"schema": "genome-call/1", "tool": "Cyrius", "kind": "pharmacogene",
               "gene": "CYP2D6", "status": "error", "note": note,
               "generated_utc": os.environ["GENERATED_UTC"]}, open(out_path, "w"), indent=2)
    sys.exit(0)

# diplotype + filter from the TSV (Sample / Genotype / Filter)
diplo, filt = None, "PASS"
if cyt and os.path.exists(cyt):
    rows = list(csv.DictReader(open(cyt), delimiter="\t"))
    row = next((r for r in rows if r.get("Sample") == sample), rows[0] if rows else None)
    if row:
        diplo = row.get("Genotype")
        filt = row.get("Filter") or "PASS"
if diplo is None:
    err("Cyrius TSV had no Genotype column/row")

# QC metrics from the JSON (Total_CN etc.), keyed by sample
tcn, raw = None, {}
if cyj and os.path.exists(cyj):
    data = json.load(open(cyj))
    raw = (data.get(sample) or (next(iter(data.values())) if isinstance(data, dict) and data else {})) or {}
    if isinstance(raw, dict):
        v = raw.get("Total_CN", raw.get("Total_CN_raw"))
        try:
            tcn = int(round(float(v))) if v not in (None, "", "None", ".") else None
        except (TypeError, ValueError):
            tcn = None

ok = bool(diplo) and str(diplo).strip().lower() not in ("none", "no_call", "", ".")
call = {
    "schema": "genome-call/1", "tool": "Cyrius", "kind": "pharmacogene",
    "gene": "CYP2D6",
    "status": "ok" if ok else "no_call",
    "diplotype": str(diplo).strip() if ok else None,
    "total_cn": tcn,
    "filter": str(filt).strip(),
    "raw": raw,
    "generated_utc": os.environ["GENERATED_UTC"],
}
json.dump(call, open(out_path, "w"), indent=2)
print(f"CYP2D6: {call['diplotype']}  (CN={call['total_cn']}, filter={call['filter']}, status={call['status']})")
PY

ok "Wrote $OUT_JSON"
hr
log "Now rebuild insights so the MCP upgrades CYP2D6:"
log "    scripts/50_serve_mcp.sh --rebuild"
