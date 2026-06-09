#!/usr/bin/env python3
"""
genome-toolkit MCP server — serves YOUR processed genome insights to any agent.

Local stdio server. Reads analysis/genome_insights.json (built by
build_insights.py) and answers structured questions; also does live variant
lookups against your variants.vcf.gz. Nothing leaves your machine.

Run (usually your MCP client launches this for you):
    python3 mcp/server.py

Env:
    GENOME_INSIGHTS   path to genome_insights.json (default: ../analysis/…)
"""
import json
import os
import re
import subprocess
import tempfile
from pathlib import Path

from mcp.server.fastmcp import FastMCP

# MCP clients often spawn servers with a minimal PATH — make sure bcftools
# (used for live variant lookups) is findable regardless.
for _d in ("/opt/homebrew/bin", "/usr/local/bin", "/usr/bin", "/bin"):
    if Path(_d).exists() and _d not in os.environ.get("PATH", "").split(os.pathsep):
        os.environ["PATH"] = _d + os.pathsep + os.environ.get("PATH", "")

INSIGHTS_PATH = Path(os.environ.get(
    "GENOME_INSIGHTS",
    Path(__file__).resolve().parents[1] / "analysis" / "genome_insights.json",
))

mcp = FastMCP("genome-toolkit")


def _load():
    if not INSIGHTS_PATH.exists():
        return None
    return json.loads(INSIGHTS_PATH.read_text())


DATA = _load()


def _err():
    return {"error": f"No insights found at {INSIGHTS_PATH}. "
                     f"Run: python3 mcp/build_insights.py (after scripts/run_all.sh)."}


# ---------------------------------------------------------------- tools -------
@mcp.tool()
def get_summary() -> dict:
    """High-level summary of the genome: sample, reference build, QC metrics
    (SNP/indel counts, Ts/Tv, het:hom), inferred genetic sex, and how many
    pharmacogenomic/trait/ClinVar findings are available."""
    if not DATA:
        return _err()
    cv = DATA.get("clinvar", {})
    return {
        "sample": DATA["sample"],
        "reference": DATA["reference"],
        "generated_utc": DATA["generated_utc"],
        "qc": DATA["qc"],
        "counts": {
            "traits": len(DATA["traits"]),
            "pharmacogenes": len(DATA["pharmacogenomics"]),
            "pgx_actionable": sum(1 for g in DATA["pharmacogenomics"] if g["actionable"]),
            "clinvar_tiers": cv.get("tier_counts") if cv.get("built") else "not built",
        },
        "disclaimer": DATA["disclaimer"],
    }


@mcp.tool()
def get_pharmacogenomics(actionable_only: bool = False) -> list:
    """Pharmacogenomic profile: predicted metabolizer/function phenotype per
    gene with drug implications (CPIC-style). Set actionable_only=True to get
    only non-default phenotypes worth flagging to a prescriber."""
    if not DATA:
        return [_err()]
    rows = DATA["pharmacogenomics"]
    if actionable_only:
        rows = [r for r in rows if r["actionable"]]
    return rows


@mcp.tool()
def check_drug(drug: str) -> dict:
    """Given a drug name, return the pharmacogenes that affect it and this
    person's predicted phenotype + guidance. Example: check_drug('clopidogrel')."""
    if not DATA:
        return _err()
    q = drug.strip().lower()
    hits = []
    for g in DATA["pharmacogenomics"]:
        if any(q in d.lower() or d.lower() in q for d in g["drugs"]):
            hits.append({"gene": g["gene"], "genotype": g["genotype"],
                         "phenotype": g["phenotype"], "actionable": g["actionable"],
                         "implication": g["implication"]})
    return {"drug": drug, "matches": hits,
            "note": "Predicted from WGS — confirm against CPIC (cpicpgx.org) with a pharmacist."
                    if hits else "No catalogued pharmacogene for this drug in the local panel."}


@mcp.tool()
def get_traits(category: str = "") -> list:
    """Trait-associated genotypes (APOE, MTHFR, lactase, caffeine, etc.).
    Optional category filter, e.g. 'Pharmacogenomics', 'Nutrition', 'Trait'."""
    if not DATA:
        return [_err()]
    rows = DATA["traits"]
    if category:
        c = category.lower()
        rows = [r for r in rows if c in r["category"].lower()]
    return rows


@mcp.tool()
def search_clinvar(gene: str = "", tier: str = "", min_stars: int = 0) -> dict:
    """Search ClinVar findings this person carries. tier in {pathogenic,
    drug_response, protective}. Filter by gene symbol and/or minimum star
    confidence (0-4). Returns tier counts plus matching rows."""
    if not DATA:
        return _err()
    cv = DATA.get("clinvar", {})
    if not cv.get("built"):
        return {"error": "ClinVar not built. Run scripts/20_clinvar_annotate.sh then rebuild insights."}
    tiers = [tier] if tier else ["pathogenic", "drug_response", "protective"]
    g = gene.strip().upper()
    rows = []
    for t in tiers:
        for r in cv.get(t, []):
            if g and g not in r["gene"].upper():
                continue
            if r["stars"] < min_stars:
                continue
            rows.append({**r, "tier": t})
    return {"tier_counts": cv["tier_counts"], "matches": rows}


@mcp.tool()
def lookup_variant(query: str) -> list:
    """Live genotype lookup by rsID (e.g. 'rs429358'), single position
    ('19:44908684'), or region ('19:44905000-44910000'). Reads your
    variants.vcf.gz directly, so it works for ANY variant, not just catalogued
    ones. An absent rsID usually means homozygous reference."""
    if not DATA:
        return [_err()]
    vcf = DATA["paths"]["variants_vcf"]
    fmt = "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t[%GT]\\t[%DP]\\t[%GQ]\\n"
    if query.strip().lower().startswith("rs"):
        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
            f.write(query.strip())
            tmp = f.name
        try:
            cmd = (f"bcftools view {vcf} 2>/dev/null | "
                   f"bcftools query -i 'ID==\"{query.strip()}\"' -f '{fmt}' 2>/dev/null")
            res = subprocess.run(cmd, shell=True, text=True, capture_output=True).stdout
        finally:
            os.unlink(tmp)
    else:
        region = query.strip()
        if not re.match(r"^(chr)?[\dXYMT]+:\d+(-\d+)?$", region):
            return [{"error": f"Unrecognized query '{query}'. Use an rsID, 'chr:pos', or 'chr:start-end'."}]
        cmd = f"bcftools query -r '{region}' -f '{fmt}' {vcf} 2>/dev/null"
        res = subprocess.run(cmd, shell=True, text=True, capture_output=True).stdout

    out = []
    for line in res.splitlines():
        p = line.split("\t")
        if len(p) < 6:
            continue
        gt = p[5]
        bases = "/".join(p[3] if a == "0" else p[4] if a == "1" else "?"
                         for a in gt.replace("|", "/").split("/"))
        zyg = ("homozygous-ref" if gt in ("0/0", "0|0")
               else "homozygous-alt" if gt in ("1/1", "1|1")
               else "hemizygous" if gt in ("0", "1")
               else "heterozygous")
        out.append({"chrom": p[0], "pos": int(p[1]), "rsid": p[2], "ref": p[3],
                    "alt": p[4], "genotype": bases, "zygosity": zyg,
                    "depth": p[6] if len(p) > 6 else None,
                    "gq": p[7] if len(p) > 7 else None})
    if not out:
        return [{"note": f"No variant record for '{query}'. If an rsID, you are "
                         f"likely homozygous reference (the common genotype)."}]
    return out


# ------------------------------------------------------------ resources -------
@mcp.resource("genome://summary")
def res_summary() -> str:
    return json.dumps(get_summary(), indent=2)


@mcp.resource("genome://pharmacogenomics")
def res_pgx() -> str:
    return json.dumps(get_pharmacogenomics(), indent=2)


@mcp.resource("genome://traits")
def res_traits() -> str:
    return json.dumps(get_traits(), indent=2)


if __name__ == "__main__":
    mcp.run()
