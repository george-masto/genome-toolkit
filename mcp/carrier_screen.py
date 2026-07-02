"""SNV-based recessive carrier screen.

The toolkit's carrier_status section is fed only by the copy-number / paralog
callers (SMN1 → SMA, GBA1 → Gaucher). But a person's ClinVar pathogenic /
likely-pathogenic hits in recessive-disease genes (CFTR, HBB, DYSF, …) are
equally valid carrier findings and were not being surfaced. This module
intersects those ClinVar hits (analysis/20_clinvar_all_categorized.tsv) with a
curated autosomal-recessive panel (assets/carriers/recessive_panel.tsv) and
classifies by zygosity:

    heterozygous          -> CARRIER (one copy; typically unaffected)
    homozygous            -> biallelic — possibly AFFECTED, clinical review
    2+ P/LP hits, one gene -> possible compound het, clinical review

Restricting to a curated recessive panel is deliberate — it keeps the screen
clinically meaningful and stops every "pathogenic"-labelled dominant/risk
variant from being mislabelled a carrier result. Pure stdlib. EDUCATIONAL ONLY;
confirm any hit with clinical testing + genetic counseling.
"""
from __future__ import annotations
from pathlib import Path

_ASSETS = Path(__file__).resolve().parents[1] / "assets" / "carriers"

# analysis/20_clinvar_all_categorized.tsv column order (see build_insights.build_clinvar)
_COLS = ("bucket", "stars", "chrom", "pos", "gene", "rsid", "geno", "zyg", "sig", "cond")
PATHOGENIC_BUCKET = "1_Pathogenic"   # the categorizer's P + likely-P tier


def load_panel(assets_dir=None) -> dict:
    path = Path(assets_dir or _ASSETS) / "recessive_panel.tsv"
    panel = {}
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        p = line.split("\t")
        if p[0] == "gene":       # header
            continue
        if len(p) >= 3:
            panel[p[0]] = {"condition": p[1], "inheritance": p[2],
                           "note": p[3] if len(p) > 3 else ""}
    return panel


def parse_clinvar_tsv(path) -> list:
    rows = []
    p = Path(path)
    if not p.exists():
        return rows
    for line in p.read_text().splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < len(_COLS):
            continue
        rows.append(dict(zip(_COLS, parts)))
    return rows


def _stars(row) -> int:
    try:
        return int(row.get("stars", 0))
    except (TypeError, ValueError):
        return 0


def screen(clinvar_rows, panel, min_stars: int = 1) -> list:
    """Intersect P/LP ClinVar rows with the recessive panel, classify by zygosity.
    Returns one entry per affected gene, sorted carriers-then-review, gene A→Z."""
    by_gene = {}
    for r in clinvar_rows:
        if r.get("bucket") != PATHOGENIC_BUCKET:
            continue
        if _stars(r) < min_stars:
            continue
        gene = r.get("gene")
        if gene not in panel:
            continue
        by_gene.setdefault(gene, []).append(r)

    results = []
    for gene, hits in by_gene.items():
        zygs = [h.get("zyg", "").lower() for h in hits]
        if any(z == "hom" for z in zygs):
            classification = "biallelic_possible"
            interp = ("Homozygous pathogenic variant — possible AFFECTED status "
                      "(not just carrier). Clinical confirmation strongly advised.")
        elif len(hits) >= 2:
            classification = "possible_compound_het"
            interp = ("Two pathogenic variants in one recessive gene — possible "
                      "compound heterozygote (affected) if in trans. Phase unknown "
                      "from short reads; clinical confirmation advised.")
        else:
            classification = "carrier"
            interp = "Heterozygous carrier of one pathogenic variant — typically unaffected."
        info = panel[gene]
        if info["inheritance"] == "XLR":
            interp += " X-linked: carrier impact depends on sex (males hemizygous can be affected)."
        results.append({
            "gene": gene,
            "condition": info["condition"],
            "inheritance": info["inheritance"],
            "classification": classification,
            "interpretation": interp,
            "gene_note": info["note"],
            "variants": [{"rsid": h.get("rsid"), "pos": f"{h.get('chrom')}:{h.get('pos')}",
                          "zygosity": h.get("zyg"), "significance": h.get("sig"),
                          "stars": _stars(h), "condition_label": h.get("cond")}
                         for h in hits],
        })
    order = {"biallelic_possible": 0, "possible_compound_het": 1, "carrier": 2}
    results.sort(key=lambda r: (order.get(r["classification"], 9), r["gene"]))
    return results


def summarize(clinvar_tsv_path, assets_dir=None, min_stars: int = 1) -> dict:
    panel = load_panel(assets_dir)
    rows = parse_clinvar_tsv(clinvar_tsv_path)
    results = screen(rows, panel, min_stars=min_stars)
    carriers = [r for r in results if r["classification"] == "carrier"]
    review = [r for r in results if r["classification"] != "carrier"]
    return {
        "results": results,
        "carrier_count": len(carriers),
        "review_flag_count": len(review),
        "panel_size": len(panel),
        "min_stars": min_stars,
        "summary": (f"{len(carriers)} recessive carrier finding(s)"
                    + (f", {len(review)} needing clinical review" if review else "")
                    + f" across a {len(panel)}-gene panel."
                    if results else
                    f"No pathogenic/likely-pathogenic carrier findings in the "
                    f"{len(panel)}-gene recessive panel (≥{min_stars}★)."),
        "disclaimer": "Educational only. A heterozygous carrier is typically unaffected but "
                      "may pass the variant to children; two carriers of the same gene have "
                      "recurrence risk. SNV-based — misses CNV/deletion carriers (e.g. alpha-"
                      "thalassemia, DMD dup/del) and can't phase compound hets. Confirm with "
                      "clinical testing + genetic counseling.",
    }
