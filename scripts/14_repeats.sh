#!/usr/bin/env bash
# ============================================================================
# 14_repeats.sh - repeat-expansion disorders with Illumina ExpansionHunter.
#
# Tandem-repeat loci (CAG, CGG, GGGGCC...) are completely invisible to a SNP/indel
# VCF - no SNP caller even looks at them. ExpansionHunter genotypes the repeat
# COUNT at known disease loci from the WGS reads: Huntington (HTT), ALS/FTD
# (C9orf72), Fragile-X (FMR1), myotonic dystrophy (DMPK), Friedreich (FXN), the
# spinocerebellar ataxias, and more. We then classify each against standard
# thresholds (assets/repeats/repeat_loci.tsv).
#
# Output:  analysis/calls/repeats.json   (build_insights.py -> repeat_expansions,
#          surfaced by the MCP's get_repeat_expansions tool).
#
#   scripts/10_align_from_fastq.sh   ->   scripts/14_repeats.sh   ->   --rebuild
#
# Env: EH (ExpansionHunter binary, default: ExpansionHunter on PATH),
#      CATALOG (GRCh38 variant catalog json; auto-located or downloaded), THREADS
# ============================================================================
source "$(dirname "$0")/lib/common.sh"
require_tool samtools
require_cram
EH="${EH:-ExpansionHunter}"
have "$EH" || die "ExpansionHunter not found. Install: scripts/00_install_tools.sh --callers
     (conda install -c bioconda expansionhunter), or set EH=/path/to/ExpansionHunter."

SAMPLE="$(sample_name)"; [[ "$SAMPLE" == "unknown" ]] && SAMPLE="sample"
WORK="$CALLS_DIR/repeats.work"; mkdir -p "$WORK"
OUT_JSON="$CALLS_DIR/repeats.json"
LOCI_TSV="$ASSETS_DIR/repeats/repeat_loci.tsv"

# --- pick the catalog matching the CRAM's contig naming ---------------------
# ExpansionHunter ships TWO GRCh38 catalogs: hg38/ (chr-prefixed: chr1, chrM...)
# and grch38/ (bare: 1, 2, MT...). A GRCh38 no-alt CRAM from 10_align is
# chr-prefixed, so default to hg38; auto-fall-back to grch38 for bare contigs.
[[ -f "$REF_FA" ]] || die "ExpansionHunter needs the reference FASTA at $REF_FA (from 10_align_from_fastq.sh)."
if samtools idxstats "$CRAM" 2>/dev/null | cut -f1 | grep -qiE '^chr'; then SUB="hg38"; else SUB="grch38"; fi
log "CRAM contigs look like '$SUB' style - using the $SUB variant catalog."
# Search known locations for an EXISTING catalog. (_first can't be used here: it
# returns its first arg unconditionally, and with nullglob off - which it is when
# sourced - an unmatched glob/path is not a real file, so it would mask a cached
# catalog and force a re-download every run / fail offline.)
if [[ -z "${CATALOG:-}" || ! -f "$CATALOG" ]]; then
  ehdir="$(dirname "$(command -v "$EH")")"
  CATALOG=""
  shopt -s nullglob
  for c in \
    "$ehdir"/../share/ExpansionHunter*/variant_catalog/$SUB/variant_catalog.json \
    "$ehdir"/../variant_catalog/$SUB/variant_catalog.json \
    "$REPO_ROOT"/tools/ExpansionHunter*/variant_catalog/$SUB/variant_catalog.json \
    "$REPO_ROOT"/tools/eh_catalog/${SUB}_variant_catalog.json; do
    [[ -f "$c" ]] && { CATALOG="$c"; break; }
  done
  shopt -u nullglob
fi
if [[ -z "$CATALOG" || ! -f "$CATALOG" ]]; then
  CATALOG="$REPO_ROOT/tools/eh_catalog/${SUB}_variant_catalog.json"
  if [[ ! -f "$CATALOG" ]]; then           # reuse a prior download / work offline
    mkdir -p "$(dirname "$CATALOG")"
    log "Downloading ExpansionHunter $SUB variant catalog..."
    curl -fsSL -o "$CATALOG" \
      "https://raw.githubusercontent.com/Illumina/ExpansionHunter/master/variant_catalog/$SUB/variant_catalog.json" \
      || die "Could not fetch the $SUB variant catalog; set CATALOG=/path/to/variant_catalog.json"
  fi
fi
log "Variant catalog: $CATALOG"

# --- sex (drives chrX loci: FMR1, AR) ---------------------------------------
# auto-detect from chrY mapped-read count unless SEX is set.
if [[ -z "${SEX:-}" ]]; then
  ycount="$(samtools idxstats "$CRAM" 2>/dev/null | awk -F'\t' '$1 ~ /^(chr)?Y$/ {s+=$3} END{print s+0}')"
  if [[ "${ycount:-0}" -gt 1000 ]]; then SEX="male"; else SEX="female"; fi
fi
log "Using --sex $SEX (chrY reads drive the call; override with SEX=male|female)."

# --- run ExpansionHunter ----------------------------------------------------
log "Running ExpansionHunter (repeat loci; can take a while)..."
THREADS_FLAG=()
"$EH" --help 2>&1 | grep -qiE -- '--threads' && THREADS_FLAG=(--threads "${THREADS:-$(sysctl -n hw.ncpu 2>/dev/null || nproc)}")
"$EH" \
  --reads "$CRAM" \
  --reference "$REF_FA" \
  --variant-catalog "$CATALOG" \
  --output-prefix "$WORK/$SAMPLE" \
  --sex "$SEX" \
  --analysis-mode seeking \
  "${THREADS_FLAG[@]+${THREADS_FLAG[@]}}" \
  || die "ExpansionHunter failed. Check $WORK and that the CRAM is GRCh38 + indexed."

# --- classify against the threshold table -> normalized json ----------------
EH_OUT="$WORK/$SAMPLE.json"
[[ -f "$EH_OUT" ]] || die "No ExpansionHunter JSON at $EH_OUT."
GENERATED_UTC="$(now_utc)" LOCI_TSV="$LOCI_TSV" EH_OUT="$EH_OUT" \
python3 - "$OUT_JSON" <<'PY' || die "Could not parse ExpansionHunter output."
import json, os, sys
out_path = sys.argv[1]
eh = json.load(open(os.environ["EH_OUT"]))
loci_results = eh.get("LocusResults", {})
# lower-cased index so catalog LocusId case differences don't matter
lr_idx = {k.lower(): v for k, v in loci_results.items()}

def max_repeat(locus_block):
    best = None
    for var in (locus_block.get("Variants") or {}).values():
        gt = var.get("Genotype")
        if not gt:
            continue
        for a in str(gt).split("/"):
            try:
                n = int(a)
            except ValueError:
                continue
            best = n if best is None else max(best, n)
    return best

loci_out = []
for line in open(os.environ["LOCI_TSV"]):
    if line.startswith("#") or not line.strip():
        continue
    p = line.rstrip("\n").split("\t")
    if len(p) < 7:
        continue
    locus_id, gene, unit, condition, normal_max, path_min, inh = p[:7]
    blk = lr_idx.get(locus_id.lower())
    if not blk:
        continue
    mx = max_repeat(blk)
    if mx is None:
        continue
    nmax, pmin = int(normal_max), int(path_min)
    if mx >= pmin:
        status, actionable = "EXPANDED (pathogenic range)", True
    elif mx > nmax:
        status, actionable = "intermediate / premutation - review", True
    else:
        status, actionable = "normal", False
    # report the full genotype string for transparency
    gts = [v.get("Genotype") for v in (blk.get("Variants") or {}).values() if v.get("Genotype")]
    loci_out.append({
        "locus": locus_id, "gene": gene, "unit": unit, "condition": condition,
        "genotype": "/".join(gts) if gts else None, "max_repeats": mx,
        "normal_max": nmax, "pathogenic_min": pmin, "status": status,
        "actionable": actionable, "inheritance": inh,
        "note": "Short-read estimate; very large expansions are flagged but not precisely sized - confirm clinically.",
    })

flagged = [l for l in loci_out if l["actionable"]]
call = {
    "schema": "genome-call/1", "tool": "ExpansionHunter", "kind": "repeat_expansion",
    "status": "ok" if loci_out else "no_call",
    "loci": loci_out,
    "summary": (f"{len(flagged)} locus/loci above normal range" if flagged
                else f"all {len(loci_out)} loci in normal range"),
    "generated_utc": os.environ["GENERATED_UTC"],
}
json.dump(call, open(out_path, "w"), indent=2)
print(f"ExpansionHunter: {len(loci_out)} loci typed, {len(flagged)} flagged "
      f"({', '.join(l['gene'] for l in flagged) or 'none'})")
PY

ok "Wrote $OUT_JSON"
log "Rebuild insights:  scripts/50_serve_mcp.sh --rebuild"
