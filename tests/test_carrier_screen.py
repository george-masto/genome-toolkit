"""Unit tests for mcp/carrier_screen.py — SNV recessive carrier screen.
Pure logic on a SYNTHETIC ClinVar table (no real genotype).
Run: `python -m pytest tests/` or `python tests/test_carrier_screen.py`.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "mcp"))
import carrier_screen as cs  # noqa: E402

PANEL = cs.load_panel()
FIXTURE = ROOT / "tests" / "fixtures" / "clinvar_synthetic.tsv"


def _by_gene(results):
    return {r["gene"]: r for r in results}


def test_panel_loads_recessive_genes():
    assert "CFTR" in PANEL and "DYSF" in PANEL
    assert PANEL["DYSF"]["inheritance"] == "AR"
    assert PANEL["G6PD"]["inheritance"] == "XLR"


def test_single_het_is_carrier():
    g = _by_gene(cs.screen(cs.parse_clinvar_tsv(FIXTURE), PANEL))
    assert g["CFTR"]["classification"] == "carrier"
    assert len(g["CFTR"]["variants"]) == 1          # 0-star duplicate filtered out


def test_homozygous_flags_biallelic():
    g = _by_gene(cs.screen(cs.parse_clinvar_tsv(FIXTURE), PANEL))
    assert g["HEXA"]["classification"] == "biallelic_possible"


def test_two_hits_flag_compound_het():
    g = _by_gene(cs.screen(cs.parse_clinvar_tsv(FIXTURE), PANEL))
    assert g["DYSF"]["classification"] == "possible_compound_het"
    assert len(g["DYSF"]["variants"]) == 2


def test_benign_and_nonpanel_ignored():
    g = _by_gene(cs.screen(cs.parse_clinvar_tsv(FIXTURE), PANEL))
    assert "SAMD11" not in g          # benign
    assert "BRCA1" not in g           # pathogenic but dominant, not on the recessive panel


def test_star_filter_drops_zero_star():
    rows = cs.parse_clinvar_tsv(FIXTURE)
    # With a high threshold, the 1-star DYSF hit drops -> DYSF no longer compound
    g = _by_gene(cs.screen(rows, PANEL, min_stars=2))
    assert len(g["DYSF"]["variants"]) == 1
    assert g["DYSF"]["classification"] == "carrier"


def test_summarize_counts():
    out = cs.summarize(FIXTURE)
    assert out["carrier_count"] == 1                 # CFTR
    assert out["review_flag_count"] == 2             # HEXA + DYSF
    assert out["panel_size"] == len(PANEL)


def test_xlr_gets_sex_dependent_note():
    rows = [{"bucket": "1_Pathogenic", "stars": "2", "chrom": "X", "pos": "154532000",
             "gene": "G6PD", "rsid": "rs5030868", "geno": "A/G", "zyg": "het",
             "sig": "Pathogenic", "cond": "G6PD deficiency"}]
    r = cs.screen(rows, PANEL)[0]
    assert "X-linked" in r["interpretation"]


def test_missing_clinvar_file_is_safe():
    out = cs.summarize(ROOT / "tests" / "fixtures" / "does_not_exist.tsv")
    assert out["results"] == [] and out["carrier_count"] == 0


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn(); print(f"  ok  {fn.__name__}")
    print(f"\n{len(fns)}/{len(fns)} carrier-screen tests passed.")
