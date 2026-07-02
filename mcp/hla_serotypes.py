"""Derive HLA serotypes (DR/DQ) and disease associations from T1K molecular
diplotypes.

T1K (scripts/13_hla.sh) already types every HLA gene it sees — class I
(HLA-A/B/C) and class II (HLA-DRB1/3/4/5, DQA1, DQB1, DPA1/B1) — and stores the
2-field diplotypes in analysis/calls/hla.json. The rest of the pipeline used to
keep only the drug-hypersensitivity (PGx) alleles. This module surfaces the
class-II types the way clinicians think about them — as serotypes (DR3, DR4,
DQ2, DQ8, ...) — and attaches the well-established class-II disease
associations, so an agent can answer "do I carry DR3/DR4?" or "am I
celiac-permissive (DQ2.5/DQ8)?" directly.

Pure stdlib; reference tables live in assets/hla/. EDUCATIONAL ONLY — HLA types
are susceptibility markers, not diagnoses.
"""
from __future__ import annotations
import os
from pathlib import Path

_ASSETS = Path(__file__).resolve().parents[1] / "assets" / "hla"

CLASS_I = ("HLA-A", "HLA-B", "HLA-C")
CLASS_II = ("HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5",
            "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1")


def _rows(path: Path):
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        parts = line.split("\t")
        if parts and parts[0] in ("locus", "tokens_all"):  # header
            continue
        yield parts


def load_serotype_map(assets_dir=None):
    """Return the allele-prefix -> serotype table, longest prefix first so that
    a second-field split (DQB1*03:02 -> DQ8) beats the group default
    (DQB1*03 -> DQ3)."""
    path = Path(assets_dir or _ASSETS) / "hla_serotypes.tsv"
    out = []
    for p in _rows(path):
        if len(p) < 3:
            continue
        locus, prefix, serotype = p[0], p[1], p[2]
        broad = p[3] if len(p) > 3 else ""
        out.append((prefix, serotype, broad))
    out.sort(key=lambda r: -len(r[0]))
    return out


def _bare(allele: str) -> str:
    """'HLA-DRB1*16:01' -> 'DRB1*16:01'; take the representative call if T1K gave
    a comma-joined ambiguity list."""
    a = (allele or "").split(",")[0].strip()
    return a[4:] if a.startswith("HLA-") else a


def allele_to_serotype(allele: str, smap) -> tuple[str | None, str]:
    """Longest-prefix match of a 2-field allele to its serotype. Returns
    (serotype, broad_group); (None, '') if unmapped."""
    bare = _bare(allele)
    for prefix, serotype, broad in smap:  # already longest-first
        if bare.startswith(prefix):
            return serotype, broad
    return None, ""


def derive_serotypes(diplotypes: dict, assets_dir=None) -> dict:
    """Map molecular diplotypes -> serotypes. Returns class_i / class_ii detail,
    the flat DR and DQ serotype lists, and a token set used for disease
    matching (serotypes + bare 2-field alleles)."""
    smap = load_serotype_map(assets_dir)
    diplotypes = diplotypes or {}

    class_i = {g: [_bare(a) for a in diplotypes.get(g, [])] for g in CLASS_I if g in diplotypes}

    class_ii, tokens = {}, set()
    dr_serotypes, dq_serotypes = [], []
    for gene in CLASS_II:
        calls = diplotypes.get(gene)
        if not calls:
            continue
        entries = []
        for a in calls:
            bare = _bare(a)
            serotype, broad = allele_to_serotype(a, smap)
            entries.append({"allele": bare, "serotype": serotype, "broad_group": broad})
            tokens.add(bare)                              # e.g. DRB1*04:01
            if serotype:
                tokens.add(serotype)                      # e.g. DR4 / DQ8 / DR52
            if broad:
                tokens.add(broad)
            # "DR type" / "DQ type" as a clinician means them = the DRB1 and DQB1
            # loci only. DRB3/4/5 (DR52/53/51) and DQA1 are supplementary and stay
            # in class_ii detail + tokens, but don't clutter the headline lists.
            if gene == "HLA-DRB1" and serotype and serotype not in dr_serotypes:
                dr_serotypes.append(serotype)
            elif gene == "HLA-DQB1" and serotype and serotype not in dq_serotypes:
                dq_serotypes.append(serotype)
        class_ii[gene] = entries

    return {
        "class_i": class_i,
        "class_ii": class_ii,
        "dr_serotypes": dr_serotypes,
        "dq_serotypes": dq_serotypes,
        "tokens": tokens,
    }


def disease_associations(tokens: set, assets_dir=None) -> list:
    """Fire every rule whose tokens are all present in `tokens`."""
    path = Path(assets_dir or _ASSETS) / "hla_disease.tsv"
    hits = []
    for p in _rows(path):
        if len(p) < 4:
            continue
        req = [t.strip() for t in p[0].split(",") if t.strip()]
        if req and all(t in tokens for t in req):
            hits.append({
                "disease": p[1], "direction": p[2], "effect": p[3],
                "note": p[4] if len(p) > 4 else "",
                "source": p[5] if len(p) > 5 else "",
                "matched_on": req,
            })
    return hits


def celiac_read(dq_serotypes: list) -> dict:
    """Celiac disease is near-obligate on HLA-DQ2.5 or DQ8. Absence of both is a
    strong negative (negative predictive value ~99%)."""
    permissive = [s for s in ("DQ2", "DQ8") if s in dq_serotypes]
    if permissive:
        return {"celiac_permissive": True, "haplotypes": permissive,
                "interpretation": f"Carries a celiac-permissive HLA type ({', '.join(permissive)}). "
                                  "Permissive, NOT diagnostic — most DQ2/DQ8 carriers never develop celiac."}
    return {"celiac_permissive": False, "haplotypes": [],
            "interpretation": "No HLA-DQ2.5 or DQ8. Celiac disease is essentially excluded on "
                              "genetic grounds (~99% of celiac patients carry DQ2/DQ8)."}


def summarize(diplotypes: dict, assets_dir=None) -> dict:
    """One call: serotypes + disease associations + celiac read + a short line."""
    s = derive_serotypes(diplotypes, assets_dir)
    diseases = disease_associations(s["tokens"], assets_dir)
    celiac = celiac_read(s["dq_serotypes"])
    dr = ", ".join(s["dr_serotypes"]) or "none typed"
    dq = ", ".join(s["dq_serotypes"]) or "none typed"
    has_dr3 = "DR3" in s["dr_serotypes"]
    has_dr4 = "DR4" in s["dr_serotypes"]
    summary = f"Class-II serotypes — DR: {dr}; DQ: {dq}. " \
              f"DR3: {'present' if has_dr3 else 'absent'}; DR4: {'present' if has_dr4 else 'absent'}."
    return {
        "class_i": s["class_i"],
        "class_ii": s["class_ii"],
        "dr_serotypes": s["dr_serotypes"],
        "dq_serotypes": s["dq_serotypes"],
        "disease_associations": diseases,
        "celiac": celiac,
        "summary": summary,
        "disclaimer": "Educational only. HLA serotypes are susceptibility markers, not diagnoses; "
                      "T1K types from short reads — confirm anything actionable with clinical HLA typing.",
    }
