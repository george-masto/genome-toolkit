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
