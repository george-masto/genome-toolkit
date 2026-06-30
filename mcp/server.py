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


def _variants_vcf():
    """Resolve the variants file portably: prefer the copy sitting next to the
    insights JSON (so moving the whole folder to another machine just works),
    fall back to the absolute path recorded at build time."""
    local = INSIGHTS_PATH.parent / "variants.vcf.gz"
    if local.exists():
        return str(local)
    return (DATA or {}).get("paths", {}).get("variants_vcf", str(local))


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


def _drug_match(q: str, drugs) -> bool:
    return any(q in d.lower() or d.lower() in q for d in (drugs or []))


@mcp.tool()
def get_polygenic_scores() -> dict:
    """Polygenic risk scores (PRS): thousands of common variants aggregated into a
    calibrated, ancestry-matched risk percentile — the lens single-variant ClinVar
    findings cannot give. Each entry has the score id/name, trait, percentile (within
    the matching ancestry group), z-score, interpretation, method, and caveats.
    Currently: coronary artery disease (PGS000018, the Inouye metaGRS)."""
    if not DATA:
        return _err()
    return {"polygenic_scores": DATA.get("polygenic_scores", []),
            "disclaimer": DATA.get("disclaimer", ""),
            "note": "Percentiles are ancestry-matched and population-RELATIVE, not "
                    "diagnostic; converting to absolute risk needs an age/sex baseline. "
                    "PRS requires DENSE genotypes — a variant-only VCF scores it "
                    "incorrectly (omits hom-ref sites where the effect allele is the "
                    "reference allele)."}


@mcp.tool()
def check_drug(drug: str) -> dict:
    """Given a drug name, return everything in this person's genome that affects
    it: pharmacogene metabolizer phenotype (CPIC), HLA hypersensitivity-risk
    alleles (e.g. abacavir/HLA-B*57:01, carbamazepine/HLA-B*15:02), and any
    copy-number-aware caller results (CYP2D6 via Cyrius, extra CYPs via Aldy).
    Example: check_drug('clopidogrel'), check_drug('codeine'), check_drug('abacavir')."""
    if not DATA:
        return _err()
    q = drug.strip().lower()
    hits = []
    # 1. core pharmacogenes (+ any CYP2D6 entry upgraded by Cyrius) and Aldy extras
    for g in DATA.get("pharmacogenomics", []) + DATA.get("advanced_pgx", []):
        if _drug_match(q, g.get("drugs")):
            hits.append({"gene": g.get("gene"), "genotype": g.get("genotype"),
                         "phenotype": g.get("phenotype"), "actionable": g.get("actionable"),
                         "method": g.get("method", "SNP panel"),
                         "implication": g.get("implication")})
    # 2. HLA hypersensitivity-risk alleles (only meaningful from a CRAM-based caller)
    for h in DATA.get("hla", []):
        if _drug_match(q, h.get("drugs")):
            hits.append({"hla_allele": h.get("allele"), "carried": h.get("carried", h.get("status")),
                         "actionable": h.get("actionable"), "method": "HLA caller",
                         "implication": h.get("implication")})
    # 3. mitochondrial MT-RNR1 (aminoglycoside ototoxicity)
    for m in DATA.get("mito", []):
        for s in m.get("mt_rnr1", []):
            if _drug_match(q, s.get("drugs")):
                hits.append({"variant": s.get("variant"), "present": s.get("present"),
                             "actionable": s.get("actionable"), "method": "mtDNA (MT-RNR1)",
                             "implication": s.get("implication")})
    return {"drug": drug, "matches": hits,
            "note": "Predicted from WGS — confirm against CPIC (cpicpgx.org) with a pharmacist."
                    if hits else "No catalogued pharmacogene or HLA risk allele for this drug in the local panel."}


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
def get_second_pass_status() -> dict:
    """What of the CRAM-based specialty-caller panel has actually been run vs is
    still pending. These callers genotype things a flat short-read VCF cannot —
    CYP2D6 copy number (Cyrius), HLA alleles, repeat expansions, SMN1/2 (SMA),
    GBA (Gaucher/Parkinson's), mtDNA. Use this to tell whether e.g. CYP2D6 is a
    real call yet or still 'incompletely typed'."""
    if not DATA:
        return _err()
    return DATA.get("second_pass", {"available": {}, "pending": {},
                                    "note": "schema predates second-pass; rebuild insights."})


@mcp.tool()
def get_carrier_status() -> dict:
    """Recessive carrier-screening results from copy-number/paralog-aware callers
    that a SNP VCF misses — e.g. SMN1/SMN2 (spinal muscular atrophy) and GBA
    (Gaucher). Empty until the relevant scripts/1x_*.sh have run on the CRAM."""
    if not DATA:
        return _err()
    rows = DATA.get("carrier_status", [])
    return {"carrier_status": rows,
            "note": "Carrier screening from WGS — confirm clinically before reproductive decisions."
                    if rows else "Not yet computed. Run scripts/15_smn.sh / 16_gba.sh on the CRAM, then rebuild insights."}


@mcp.tool()
def get_repeat_expansions() -> dict:
    """Repeat-expansion disorder genotypes (HTT/Huntington, C9orf72/ALS-FTD,
    FMR1/Fragile-X, DM1, SCAs, Friedreich's) from ExpansionHunter. These STR loci
    are invisible to a SNP VCF. Empty until scripts/14_repeats.sh has run."""
    if not DATA:
        return _err()
    rows = DATA.get("repeat_expansions", [])
    return {"repeat_expansions": rows,
            "note": "Repeat sizes from short-read WGS; very large expansions are flagged "
                    "but not precisely sized — confirm flagged loci clinically."
                    if rows else "Not yet computed. Run scripts/14_repeats.sh on the CRAM, then rebuild insights."}


@mcp.tool()
def get_hla() -> dict:
    """HLA alleles relevant to drug-hypersensitivity pharmacogenetics
    (HLA-B*57:01 abacavir; B*15:02 / A*31:01 carbamazepine; B*58:01 allopurinol)
    plus disease-associated HLA. Invisible to a SNP VCF. Empty until
    scripts/13_hla.sh has run. For a specific drug, prefer check_drug()."""
    if not DATA:
        return _err()
    rows = DATA.get("hla", [])
    return {"hla": rows,
            "note": "HLA typed from WGS reads — confirm before acting on a hypersensitivity flag."
                    if rows else "Not yet computed. Run scripts/13_hla.sh on the CRAM, then rebuild insights."}


@mcp.tool()
def get_mitochondrial() -> dict:
    """Mitochondrial findings: haplogroup, MT-RNR1 variants (m.1555A>G / m.1494C>T
    — CPIC: aminoglycoside-induced hearing loss), and heteroplasmy/contamination.
    Empty until scripts/17_mtdna.sh has run."""
    if not DATA:
        return _err()
    rows = DATA.get("mito", [])
    return {"mito": rows,
            "note": "Mitochondrial calls from WGS chrM reads."
                    if rows else "Not yet computed. Run scripts/17_mtdna.sh on the CRAM, then rebuild insights."}


@mcp.tool()
def lookup_variant(query: str) -> list:
    """Live genotype lookup by rsID (e.g. 'rs429358'), single position
    ('19:44908684'), or region ('19:44905000-44910000'). Reads your
    variants.vcf.gz directly, so it works for ANY variant, not just catalogued
    ones. An absent rsID usually means homozygous reference."""
    if not DATA:
        return [_err()]
    vcf = _variants_vcf()
    fmt = "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t[%GT]\\t[%DP]\\t[%GQ]\\n"
    rsid = query.strip()
    if rsid.lower().startswith("rs"):
        # An rsID is 'rs' + digits. This strict check is also the injection guard:
        # `rsid` is interpolated into the shell command below, so anything that is
        # not rs<digits> must be rejected before it can reach the shell.
        if not re.fullmatch(r"rs\d+", rsid, re.IGNORECASE):
            return [{"error": f"Invalid rsID '{query}'. Expected 'rs' followed by digits, e.g. rs429358."}]
        cmd = (f"bcftools view {vcf} 2>/dev/null | "
               f"bcftools query -i 'ID==\"{rsid}\"' -f '{fmt}' 2>/dev/null")
        res = subprocess.run(cmd, shell=True, text=True, capture_output=True).stdout
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


# ---------------------------------------------------- semantic search ---------
# Loaded lazily on first use: the embedding model + the prebuilt local index.
_EMB = {"ready": False, "model": None, "vecs": None, "metas": None}


def _ensure_embeddings():
    if _EMB["ready"]:
        return True
    npz = INSIGHTS_PATH.parent / "genome_embeddings.npz"
    jsonl = INSIGHTS_PATH.parent / "genome_embeddings.jsonl"
    if not npz.exists() or not jsonl.exists():
        return False
    try:
        import numpy as np
        from fastembed import TextEmbedding
    except ImportError:
        return False
    d = np.load(npz, allow_pickle=True)
    _EMB["vecs"] = d["vectors"]
    _EMB["metas"] = [json.loads(line) for line in jsonl.read_text().splitlines()]
    _EMB["model"] = TextEmbedding(model_name=str(d["model"]))
    _EMB["np"] = np
    _EMB["ready"] = True
    return True


@mcp.tool()
def semantic_search(query: str, top_k: int = 10, source: str = "", min_stars: int = 0) -> dict:
    """Fuzzy / meaning-based search over YOUR genome findings — finds matches by
    concept, not keywords. Use for open-ended questions where the exact words
    won't appear in the data, e.g. 'connective tissue', 'heart rhythm',
    'serotonin', 'inflammation'. Optional filters: source in
    {clinvar, pharmacogenomics, trait}, and min_stars (ClinVar confidence 0-4).
    For an exact rsID/gene/drug, prefer lookup_variant / check_drug / search_clinvar."""
    if not _ensure_embeddings():
        return {"error": "Semantic index not built. Run: scripts/50_serve_mcp.sh --embeddings"}
    np = _EMB["np"]
    q = np.asarray(next(iter(_EMB["model"].embed([query]))), dtype="float32")
    q /= (np.linalg.norm(q) + 1e-9)
    sims = _EMB["vecs"] @ q
    order = np.argsort(-sims)
    out = []
    for i in order:
        m = _EMB["metas"][int(i)]
        if source and m.get("source") != source:
            continue
        if m.get("stars", 99) < min_stars:
            continue
        out.append({"score": round(float(sims[int(i)]), 3), **m})
        if len(out) >= top_k:
            break
    return {"query": query, "matches": out,
            "note": "Semantically ranked over your findings. Educational only; "
                    "confirm anything actionable with a clinician."}


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


def _run_http():
    """Serve over streamable-HTTP. Bind to a Tailscale IP (or localhost) so the
    endpoint is reachable only inside your tailnet. Optional bearer-token auth
    via MCP_AUTH_TOKEN for defense in depth on top of the network boundary."""
    import uvicorn
    from mcp.server.transport_security import TransportSecuritySettings
    host = os.environ.get("MCP_HOST", "127.0.0.1")
    port = int(os.environ.get("MCP_PORT", "8787"))
    token = os.environ.get("MCP_AUTH_TOKEN", "").strip()

    # FastMCP's DNS-rebinding Host check is a browser-localhost defense; it 421s
    # any non-localhost Host (e.g. a Tailscale IP). Our clients are token-authed
    # agents on a private tailnet, not browsers — so either pin allowed hosts via
    # MCP_ALLOWED_HOSTS or disable the check (network scope + token are the guard).
    allowed = os.environ.get("MCP_ALLOWED_HOSTS", "").strip()
    if allowed:
        mcp.settings.transport_security = TransportSecuritySettings(
            enable_dns_rebinding_protection=True,
            allowed_hosts=[h.strip() for h in allowed.split(",") if h.strip()],
            allowed_origins=["*"])
    else:
        mcp.settings.transport_security = TransportSecuritySettings(
            enable_dns_rebinding_protection=False)
    app = mcp.streamable_http_app()  # MCP mounted at /mcp
    if token:
        from starlette.middleware.base import BaseHTTPMiddleware
        from starlette.responses import JSONResponse

        class _Bearer(BaseHTTPMiddleware):
            async def dispatch(self, request, call_next):
                if request.headers.get("authorization", "") != f"Bearer {token}":
                    return JSONResponse({"error": "unauthorized"}, status_code=401)
                return await call_next(request)

        app.add_middleware(_Bearer)
        print(f"[genome-toolkit] HTTP on {host}:{port}/mcp (bearer-token auth ON)", flush=True)
    else:
        print(f"[genome-toolkit] HTTP on {host}:{port}/mcp "
              f"(no token — relying on network scope; set MCP_AUTH_TOKEN to add one)", flush=True)
    uvicorn.run(app, host=host, port=port, log_level="warning")


if __name__ == "__main__":
    transport = os.environ.get("MCP_TRANSPORT", "stdio").lower()
    if transport in ("http", "streamable-http"):
        _run_http()
    else:
        mcp.run()
