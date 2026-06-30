#!/usr/bin/env python3
"""
build_insights.py — distill the toolkit's analysis into ONE structured file:
    analysis/genome_insights.json

This is the machine-readable knowledge source the MCP server serves. It is
written into analysis/ (gitignored) — your genome never enters the repo.

It re-derives everything from your own data + the committed public catalogs:
  - QC headline numbers + inferred sex        (from variants.pass.vcf.gz)
  - trait genotypes                            (assets/traits/known_variants.tsv)
  - pharmacogenomic phenotypes                 (assets/pharmacogenes/pgx_variants.tsv)
  - ClinVar tiers                              (analysis/20_clinvar_all_categorized.tsv, if built)

Usage:  python3 mcp/build_insights.py
Requires: bcftools on PATH, and scripts/03_extract_variants.sh already run.
"""
import json, os, re, subprocess, sys, tempfile
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get("OUT_DIR", REPO / "analysis"))
ASSETS = REPO / "assets"
VARIANTS = OUT / "variants.vcf.gz"
PASS = OUT / "variants.pass.vcf.gz"
CLINVAR_TSV = OUT / "20_clinvar_all_categorized.tsv"
# Normalized output of the second-pass specialty callers (Cyrius, ExpansionHunter,
# SMNCopyNumberCaller, Gauchian, HLA, mtDNA). Each writes analysis/calls/<name>.json;
# build_insights ingests whatever is present and degrades gracefully when absent.
CALLS_DIR = OUT / "calls"


def sh(cmd: str) -> str:
    return subprocess.run(cmd, shell=True, text=True, capture_output=True).stdout


def need(path: Path, hint: str):
    if not path.exists():
        sys.exit(f"Missing {path}\n  -> {hint}")


need(VARIANTS, "run scripts/03_extract_variants.sh first")
if not PASS.exists():
    PASS = VARIANTS


# ---------- helpers ----------------------------------------------------------
def chrom_prefix() -> str:
    hdr = sh(f"bcftools view -h {VARIANTS} 2>/dev/null")
    return "chr" if "##contig=<ID=chr" in hdr else ""


CP = chrom_prefix()


def sample_and_ref():
    sample = sh(f"bcftools query -l {VARIANTS} 2>/dev/null").strip().splitlines()
    sample = sample[0] if sample else "unknown"
    ref = ""
    for line in sh(f"bcftools view -h {VARIANTS} 2>/dev/null").splitlines():
        if line.startswith("##reference="):
            ref = line.split("=", 1)[1]
            break
    return sample, ref


def query_genotypes(rsids):
    """One scan of the variants file for a set of rsIDs -> {rsid: {...}}."""
    out = {}
    if not rsids:
        return out
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        f.write("\n".join(sorted(rsids)))
        tmp = f.name
    try:
        res = sh(f"bcftools view {VARIANTS} 2>/dev/null | "
                 f"bcftools query -i 'ID=@{tmp}' "
                 f"-f '%ID\\t%CHROM\\t%POS\\t%REF\\t%ALT\\t[%GT]\\t[%DP]\\n' 2>/dev/null")
        for line in res.splitlines():
            p = line.split("\t")
            if len(p) < 6:
                continue
            rid, chrom, pos, ref, alt, gt = p[0], p[1], p[2], p[3], p[4], p[5]
            dp = p[6] if len(p) > 6 else "."
            out.setdefault(rid, {"chrom": chrom, "pos": int(pos), "ref": ref,
                                 "alt": alt, "gt": gt, "dp": dp})
    finally:
        os.unlink(tmp)
    return out


def geno_bases(rec):
    """Return genotype as bases, e.g. 'C/T', from a found record."""
    alleles = []
    for a in rec["gt"].replace("|", "/").split("/"):
        alleles.append(rec["ref"] if a == "0" else rec["alt"] if a == "1" else "?")
    return "/".join(alleles)


# ---------- QC ---------------------------------------------------------------
def build_qc():
    stats = sh(f"bcftools stats -s - {PASS} 2>/dev/null")
    qc = {"snps": None, "indels": None, "ts_tv": None, "het_hom": None}
    for line in stats.splitlines():
        if line.startswith("SN") and "number of SNPs:" in line:
            qc["snps"] = int(line.split("\t")[-1])
        elif line.startswith("SN") and "number of indels:" in line:
            qc["indels"] = int(line.split("\t")[-1])
        elif line.startswith("TSTV"):
            qc["ts_tv"] = float(line.split("\t")[4])
        elif line.startswith("PSC"):
            c = line.split("\t")
            hom, het = int(c[4]), int(c[5])
            qc["het_hom"] = round(het / hom, 3) if hom else None
    # inferred sex from chrX ploidy + chrY calls
    try:
        xgt = sh(f"bcftools view -r {CP}X {PASS} 2>/dev/null | "
                 f"bcftools query -f '[%GT]\\n' 2>/dev/null").splitlines()
        hap = sum(1 for g in xgt if g in ("0", "1"))
        het = sum(1 for g in xgt if g in ("0/1", "1/0", "0|1", "1|0"))
        ycount = len(sh(f"bcftools view -r {CP}Y {PASS} 2>/dev/null | "
                        f"bcftools query -f '%CHROM\\n' 2>/dev/null").splitlines())
        if het > 1000:
            qc["inferred_sex"] = "XX (female)"
        elif hap > 1000 and ycount > 100:
            qc["inferred_sex"] = "XY (male)"
        else:
            qc["inferred_sex"] = "uncertain"
    except Exception:
        qc["inferred_sex"] = "unknown"
    return qc


# ---------- traits -----------------------------------------------------------
def load_catalog(path):
    rows = []
    for line in path.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        rows.append(line.split("\t"))
    return rows


def build_traits(found):
    cat = load_catalog(ASSETS / "traits" / "known_variants.tsv")
    out = []
    for r in cat:
        if len(r) < 10:
            continue
        rsid, gene, chrom, pos, ref, alt, eff, category, trait, effect = r[:10]
        if rsid in found:
            rec = found[rsid]
            g = geno_bases(rec)
            copies = g.split("/").count(eff)
        else:
            g = f"{ref}/{ref}"
            copies = 2 if eff == ref else 0  # effect allele is ALT -> 0 when hom-ref
        out.append({"rsid": rsid, "gene": gene, "category": category,
                    "trait": trait, "genotype": g, "effect_allele": eff,
                    "effect_allele_copies": copies, "effect": effect})
    return out


# ---------- second-pass specialty calls --------------------------------------
def load_calls():
    """Read every analysis/calls/*.json the specialty callers have produced.
    Returns {stem: parsed}. Malformed/partial files are skipped, never fatal —
    a half-written call must not break the whole insights rebuild."""
    calls = {}
    if not CALLS_DIR.exists():
        return calls
    for p in sorted(CALLS_DIR.glob("*.json")):
        try:
            calls[p.stem] = json.loads(p.read_text())
        except (json.JSONDecodeError, OSError):
            continue
    return calls


def load_cyp2d6_table():
    """CPIC CYP2D6 allele activity values + activity-score->phenotype thresholds,
    from the committed public reference table."""
    path = ASSETS / "pharmacogenes" / "cyp2d6_allele_functions.tsv"
    act, thresholds = {}, []
    if not path.exists():
        return act, thresholds
    for line in path.read_text().splitlines():
        if line.startswith("#PHENO\t"):
            _, mx, label = line.split("\t")
            thresholds.append((float(mx), label))
        elif line.startswith("#") or not line.strip():
            continue
        else:
            p = line.split("\t")
            if len(p) >= 3:
                act[p[0]] = float(p[2])
    thresholds.sort(key=lambda t: t[0])
    return act, thresholds


D6_ACT, D6_THRESH = load_cyp2d6_table()


def cyp2d6_phenotype(diplotype):
    """Translate a Cyrius diplotype -> (activity_score, phenotype, all_known).
    Handles 'xN' duplication multipliers and '+' structural tandems (e.g.
    '*36+*10/*1', '*2x2/*4'). Returns (None, None, False) if any allele token is
    not in the CPIC table, so we never invent a phenotype for an unknown allele."""
    if not D6_ACT or not D6_THRESH:
        return None, None, False
    # if Cyrius reported ambiguous alternates, score the first candidate
    dip = re.split(r"\s+or\s+|;", diplotype.strip())[0].strip()
    sides = dip.split("/")
    if len(sides) != 2:
        return None, None, False
    total = 0.0
    for side in sides:                       # each side is one haplotype
        for tok in side.split("+"):          # '+' joins tandem gene copies
            m = re.match(r"(\*[0-9A-Za-z.]+?)(?:x(\d+))?$", tok.strip())
            if not m:
                return None, None, False
            core = re.sub(r"\.\d+$", "", m.group(1))   # suballele *4.001 -> core *4
            if core not in D6_ACT:
                return None, None, False
            total += D6_ACT[core] * int(m.group(2) or 1)
    pheno = next((label for mx, label in D6_THRESH if total <= mx), D6_THRESH[-1][1])
    return round(total, 3), pheno, True


# ---------- pharmacogenomics (CPIC engine) -----------------------------------
FUNCVAL = {"no_function": 0.0, "decreased": 0.5, "increased": 1.5}


def metab(act):
    if act <= 0.5: return "Poor metabolizer"
    if act <= 1.5: return "Intermediate metabolizer"
    if act <= 2.25: return "Normal metabolizer"
    if act <= 2.75: return "Rapid metabolizer"
    return "Ultrarapid metabolizer"


def build_pgx(found, cyp2d6_call=None):
    cat = load_catalog(ASSETS / "pharmacogenes" / "pgx_variants.tsv")
    genes = {}
    order = []
    for r in cat:
        if len(r) < 7:
            continue
        gene, rsid, defines, fn, model, drugs, action = r[:7]
        if gene not in genes:
            genes[gene] = {"model": model, "drugs": drugs, "action": action,
                           "act": 2.0, "carried": [], "slow": 0, "g6pd": 0,
                           "special": 0}
            order.append(gene)
        g = genes[gene]
        copies = 0
        if rsid in found:
            copies = found[rsid]["gt"].replace("|", "/").split("/").count("1")
        f = FUNCVAL.get(fn, 1.0)
        g["act"] += copies * (f - 1.0)
        if copies:
            star = defines[len(gene):] if defines.startswith(gene) else defines
            g["carried"].append((star.strip(), copies))
        if model == "nat2" and copies:
            g["slow"] += copies
        if model == "g6pd" and copies:
            g["g6pd"] += copies
        if model in ("vkorc1", "cyp4f2", "ugt", "cyp3a5", "slco", "ifnl3"):
            g["special"] = copies  # single-variant genes

    out = []
    for gene in order:
        g = genes[gene]
        m = g["model"]
        act = g["act"]
        diplo, ph, neutral = "", "", False
        method = "SNP panel (short-read VCF)"
        extra = {}
        if m in ("score", "lof"):
            # build *1/* diplotype
            alleles = []
            for star, c in g["carried"]:
                alleles += [star] * c
            while len(alleles) < 2:
                alleles.insert(0, "*1")
            diplo = "/".join(alleles)
            ph = metab(act)
            if m == "lof" and ph in ("Rapid metabolizer", "Ultrarapid metabolizer"):
                ph = "Normal metabolizer"
            neutral = (ph == "Normal metabolizer")
        elif m == "vkorc1":
            c = g["special"]; diplo = ["-1639 GG", "-1639 GA", "-1639 AA"][c]
            ph = ["Normal warfarin sensitivity", "Increased sensitivity (lower dose)",
                  "High sensitivity (lowest dose)"][c]; neutral = (c == 0)
        elif m == "cyp4f2":
            c = g["special"]; diplo = ["*1/*1", "*1/*3", "*3/*3"][c]
            ph = "Normal" if c == 0 else "Carries *3 (minor: slightly higher warfarin dose)"
            neutral = (c == 0)
        elif m == "slco":
            c = g["special"]; diplo = ["*1/*1", "*1/*5", "*5/*5"][c]
            ph = ["Normal function", "Decreased function (higher statin myopathy risk)",
                  "Poor function (high statin myopathy risk)"][c]; neutral = (c == 0)
        elif m == "ugt":
            c = g["special"]; diplo = ["*1/*1", "*1/*28", "*28/*28"][c]
            ph = ["Normal", "Intermediate (Gilbert's carrier)",
                  "Reduced (Gilbert's; irinotecan caution)"][c]; neutral = (c == 0)
        elif m == "cyp3a5":
            c = g["special"]; diplo = ["*1/*1", "*1/*3", "*3/*3"][c]
            ph = ["Expressor (higher tacrolimus dose)", "Intermediate expressor",
                  "Non-expressor (standard tacrolimus dose)"][c]; neutral = (c == 2)
        elif m == "nat2":
            s = g["slow"]; diplo = f"{s} slow allele(s)"
            ph = "Slow acetylator" if s >= 2 else "Intermediate acetylator" if s == 1 else "Rapid/normal acetylator"
            neutral = (s == 0)
        elif m == "ifnl3":
            c = g["special"]; diplo = ["CC", "CT", "TT"][c]
            ph = ["Favorable (CC)", "Intermediate (CT)", "Unfavorable (TT)"][c]; neutral = (c == 0)
        elif m == "g6pd":
            c = g["g6pd"]; diplo = "none detected" if c == 0 else "variant present"
            ph = "No common deficiency variant (normal)" if c == 0 else "Possible deficiency — confirm by enzyme assay"
            neutral = (c == 0)
        elif m == "cyp2d6_partial":
            call = cyp2d6_call if isinstance(cyp2d6_call, dict) else None
            if call and call.get("status") == "ok" and call.get("diplotype"):
                # A real CYP2D6 call from Cyrius — translate it via the CPIC
                # activity-score table. This is the upgrade that flips CYP2D6
                # from "incompletely typed" to a usable predicted phenotype.
                score, pheno, known = cyp2d6_phenotype(call["diplotype"])
                # Cyrius Total_CN is the COMBINED CYP2D6+CYP2D7 copy number (≈4 is
                # normal: 2 gene + 2 pseudogene), NOT a CYP2D6 duplication. The
                # diplotype itself encodes any CYP2D6 duplication via xN, so show
                # the diplotype clean and keep the combined CN only as a QC field.
                cn = call.get("total_cn")
                diplo = call["diplotype"]
                method = f"Cyrius (CYP2D6-specific caller){'' if call.get('filter','PASS') == 'PASS' else ' - ' + str(call.get('filter'))}"
                impl = ("Typed by Cyrius from WGS reads (CYP2D6/CYP2D7-paralog and "
                        "copy-number aware) — a real upgrade over the SNP panel. "
                        "Research-grade: confirm in a CLIA lab before prescribing. "
                        "Phenotype drives CPIC dosing for codeine/tramadol (avoid in "
                        "poor & ultrarapid metabolizers), tamoxifen, atomoxetine, and "
                        "many antidepressants/antipsychotics.")
                if known:
                    ph = pheno + (f" (activity score {score})" if score is not None else "")
                    neutral = (pheno == "Normal metabolizer")
                    extra = {"activity_score": score, "cyrius_total_cn_d6d7": cn,
                             "filter": call.get("filter", "PASS"), "implication": impl}
                else:
                    ph = "Typed by Cyrius; phenotype indeterminate (uncharacterized allele)"
                    neutral = None
                    extra = {"cyrius_total_cn_d6d7": cn, "filter": call.get("filter", "PASS"),
                             "implication": impl}
            else:
                alleles = [s for s, c in g["carried"] for _ in range(c)]
                diplo = "carries " + ", ".join(alleles) if alleles else "no *4/*10/*41 detected"
                ph = "Incompletely typed (copy-number not callable)"; neutral = None
        else:
            diplo, ph = "?", "see details"
        out.append({"gene": gene, "genotype": diplo, "phenotype": ph,
                    "actionable": (neutral is False), "typed": (neutral is not None),
                    "method": method,
                    "drugs": [d.strip() for d in g["drugs"].replace(";", ",").split(",") if d.strip()],
                    "implication": g["action"], **extra})
    return out


# ---------- ClinVar ----------------------------------------------------------
def build_clinvar():
    if not CLINVAR_TSV.exists():
        return {"built": False,
                "note": "Run scripts/20_clinvar_annotate.sh to populate ClinVar tiers."}
    bucket_label = {"1_Pathogenic": "pathogenic", "2_DrugResponse": "drug_response",
                    "3_RiskFactor": "risk_factor", "4_Protective": "protective",
                    "5_Association": "association", "6_Uncertain": "uncertain",
                    "7_Other": "other", "8_Conflicting": "conflicting", "9_Benign": "benign"}
    counts, rows = {}, {"pathogenic": [], "drug_response": [], "protective": []}
    for line in CLINVAR_TSV.read_text().splitlines():
        p = line.split("\t")
        if len(p) < 10:
            continue
        bucket, stars, chrom, pos, gene, rsid, geno, zyg, sig, cond = p[:10]
        label = bucket_label.get(bucket, bucket)
        counts[label] = counts.get(label, 0) + 1
        if label in rows:
            rows[label].append({"gene": gene, "rsid": rsid, "genotype": geno,
                                "zygosity": zyg, "stars": int(stars),
                                "significance": sig, "condition": cond})
    for k in rows:
        rows[k].sort(key=lambda r: -r["stars"])
    return {"built": True, "tier_counts": counts, **rows}


# ---------- advanced (second-pass caller) sections ---------------------------
# The full panel of specialty callers that run on the re-aligned CRAM. The stem
# is the analysis/calls/<stem>.json filename each wrapper writes.
SECOND_PASS_PANEL = {
    "cyp2d6":  ("11_cyp2d6.sh",  "CYP2D6 diplotype + copy number (Cyrius) — codeine/tramadol, "
                                 "tamoxifen, atomoxetine, many antidepressants/antipsychotics"),
    "hla":     ("13_hla.sh",     "HLA pharmacogenetic alleles — HLA-B*57:01 (abacavir), "
                                 "B*15:02/A*31:01 (carbamazepine), B*58:01 (allopurinol)"),
    "repeats": ("14_repeats.sh", "Repeat-expansion disorders — HTT, C9orf72, FMR1, SCAs, DM1, FXN"),
    "smn":     ("15_smn.sh",     "SMN1/SMN2 copy number — spinal muscular atrophy carrier status"),
    "gba":     ("16_gba.sh",     "GBA — Gaucher disease / Parkinson's risk (Gauchian)"),
    "mito":    ("17_mtdna.sh",   "Mitochondrial haplogroup + MT-RNR1 (aminoglycoside ototoxicity)"),
}


def build_advanced(calls):
    """Route the specialty-caller outputs into MCP-facing sections. Permissive by
    design: each caller owns its own payload shape; if a call bundles many items
    under 'loci'/'alleles'/'items' we flatten those, else we keep the object."""
    sections = {"carrier_status": [], "repeat_expansions": [], "hla": [],
                "mito": [], "advanced_pgx": [], "polygenic_scores": []}
    routed = {"carrier": "carrier_status", "repeat_expansion": "repeat_expansions",
              "hla": "hla", "mito": "mito", "pharmacogene": "advanced_pgx",
              "polygenic": "polygenic_scores"}
    for c in calls.values():
        if not isinstance(c, dict):
            continue
        dest = routed.get(c.get("kind"))
        if not dest:
            continue
        # the CYP2D6 pharmacogene call is consumed by build_pgx — don't double-list it
        if c.get("kind") == "pharmacogene" and c.get("gene") == "CYP2D6":
            continue
        items = c.get("loci") or c.get("alleles") or c.get("items")
        sections[dest].extend(items) if isinstance(items, list) else sections[dest].append(c)
    return sections


def build_second_pass_summary(calls):
    """What of the CRAM-based specialty panel has actually been run vs is pending.
    Lets an agent see at a glance whether e.g. CYP2D6 is a real call yet."""
    available, pending = {}, {}
    for stem, (script, desc) in SECOND_PASS_PANEL.items():
        c = calls.get(stem)
        if isinstance(c, dict) and c.get("status") == "ok":
            available[stem] = desc
        else:
            pending[stem] = {"run": f"scripts/{script}", "yields": desc}
    return {
        "available": available,
        "pending": pending,
        "note": "These callers genotype things a flat short-read VCF cannot "
                "(copy number, paralog-confounded genes, HLA, repeat expansions). "
                "They consume the CRAM from scripts/10_align_from_fastq.sh; run the "
                "listed script, then rebuild insights (scripts/50_serve_mcp.sh --rebuild).",
    }


# ---------- assemble ---------------------------------------------------------
def main():
    sample, ref = sample_and_ref()
    trait_cat = load_catalog(ASSETS / "traits" / "known_variants.tsv")
    pgx_cat = load_catalog(ASSETS / "pharmacogenes" / "pgx_variants.tsv")
    rsids = {r[0] for r in trait_cat if len(r) >= 10} | {r[1] for r in pgx_cat if len(r) >= 7}
    print(f"Scanning {len(rsids)} catalog variants…", file=sys.stderr)
    found = query_genotypes(rsids)

    calls = load_calls()           # second-pass specialty-caller outputs, if any
    advanced = build_advanced(calls)

    insights = {
        "schema_version": 2,
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sample": sample,
        "reference": ref,
        "disclaimer": "Educational only. Not medical advice. Predicted from "
                      "consumer WGS; confirm anything actionable with a clinician "
                      "or pharmacist (pharmacogenomics: CPIC / cpicpgx.org).",
        "qc": build_qc(),
        "traits": build_traits(found),
        "pharmacogenomics": build_pgx(found, calls.get("cyp2d6")),
        "clinvar": build_clinvar(),
        "carrier_status": advanced["carrier_status"],
        "repeat_expansions": advanced["repeat_expansions"],
        "hla": advanced["hla"],
        "mito": advanced["mito"],
        "advanced_pgx": advanced["advanced_pgx"],
        "polygenic_scores": advanced["polygenic_scores"],
        "second_pass": build_second_pass_summary(calls),
        "paths": {"variants_vcf": str(VARIANTS)},
    }
    out = OUT / "genome_insights.json"
    out.write_text(json.dumps(insights, indent=2))
    pgx_act = sum(1 for g in insights["pharmacogenomics"] if g["actionable"])
    d6 = next((g for g in insights["pharmacogenomics"] if g["gene"] == "CYP2D6"), {})
    sp_avail = ", ".join(insights["second_pass"]["available"]) or "none yet"
    print(f"Wrote {out}", file=sys.stderr)
    print(f"  sample={sample} sex={insights['qc'].get('inferred_sex')} "
          f"traits={len(insights['traits'])} pgx={len(insights['pharmacogenomics'])} "
          f"(actionable={pgx_act}) clinvar={'yes' if insights['clinvar'].get('built') else 'no'}",
          file=sys.stderr)
    print(f"  CYP2D6={d6.get('genotype')} [{d6.get('method')}]  "
          f"second-pass available: {sp_avail}", file=sys.stderr)


if __name__ == "__main__":
    main()
