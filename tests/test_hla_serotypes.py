"""Unit tests for mcp/hla_serotypes.py — HLA class-II serotype derivation and
disease-association mapping. Pure logic on SYNTHETIC diplotypes (no real
genotype). Run: `python -m pytest tests/` or `python tests/test_hla_serotypes.py`.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "mcp"))
import hla_serotypes as hs  # noqa: E402

SMAP = hs.load_serotype_map()


def _sero(allele):
    return hs.allele_to_serotype(allele, SMAP)[0]


def test_drb1_group_to_dr_serotype():
    assert _sero("HLA-DRB1*16:01") == "DR16"
    assert _sero("HLA-DRB1*11:04") == "DR11"
    assert _sero("HLA-DRB1*03:01") == "DR3"
    assert _sero("HLA-DRB1*04:01") == "DR4"
    assert _sero("HLA-DRB1*15:01") == "DR15"


def test_dqb1_second_field_splits_beat_group_default():
    # DQB1*03 group default is DQ3, but the 2-field splits must win
    assert _sero("HLA-DQB1*03:01") == "DQ7"
    assert _sero("HLA-DQB1*03:02") == "DQ8"
    assert _sero("HLA-DQB1*03:03") == "DQ9"
    assert _sero("HLA-DQB1*02:01") == "DQ2"
    assert _sero("HLA-DQB1*05:02") == "DQ5"


def test_secondary_drb_loci_map_to_dr51_52_53():
    assert _sero("HLA-DRB3*02:02") == "DR52"
    assert _sero("HLA-DRB4*01:01") == "DR53"
    assert _sero("HLA-DRB5*02:02") == "DR51"


def test_comma_joined_ambiguity_list_uses_first_call():
    ambiguous = "HLA-DRB4*01:01,HLA-DRB4*01:02,HLA-DRB4*01:03"
    assert _sero(ambiguous) == "DR53"


def test_unmapped_allele_returns_none():
    assert _sero("HLA-DRB1*99:99") is None


def test_derive_collects_dr_and_dq_lists():
    diplo = {"HLA-DRB1": ["HLA-DRB1*16:01", "HLA-DRB1*11:04"],
             "HLA-DQB1": ["HLA-DQB1*05:02", "HLA-DQB1*03:01"]}
    d = hs.derive_serotypes(diplo)
    assert set(d["dr_serotypes"]) == {"DR16", "DR11"}
    assert set(d["dq_serotypes"]) == {"DQ5", "DQ7"}
    assert "DRB1*04:01" not in d["tokens"]


def test_supplementary_loci_do_not_pollute_headline_serotypes():
    # DRB3/4/5 (DR52/53/51) and DQA1 must stay out of the DR/DQ headline lists,
    # which are DRB1- and DQB1-only. (Regression: DQA1*01 once leaked into dq.)
    diplo = {"HLA-DRB1": ["HLA-DRB1*16:01", "HLA-DRB1*11:04"],
             "HLA-DRB3": ["HLA-DRB3*02:02"], "HLA-DRB4": ["HLA-DRB4*01:01"],
             "HLA-DRB5": ["HLA-DRB5*02:02"],
             "HLA-DQA1": ["HLA-DQA1*01:02", "HLA-DQA1*05:05"],
             "HLA-DQB1": ["HLA-DQB1*05:02", "HLA-DQB1*03:01"]}
    d = hs.derive_serotypes(diplo)
    assert d["dr_serotypes"] == ["DR16", "DR11"]
    assert d["dq_serotypes"] == ["DQ5", "DQ7"]
    # supplementary serotypes still available as tokens
    assert {"DR52", "DR53", "DR51"} <= d["tokens"]


def test_disease_dr3_dr4_fires_type1_diabetes():
    diplo = {"HLA-DRB1": ["HLA-DRB1*03:01", "HLA-DRB1*04:01"],
             "HLA-DQB1": ["HLA-DQB1*02:01", "HLA-DQB1*03:02"]}
    out = hs.summarize(diplo)
    diseases = {a["disease"] for a in out["disease_associations"]}
    assert "Type 1 diabetes" in diseases
    assert out["celiac"]["celiac_permissive"] is True          # DQ2 + DQ8 both present
    t1d = next(a for a in out["disease_associations"] if a["disease"] == "Type 1 diabetes"
               and a["matched_on"] == ["DR3", "DR4"])
    assert t1d["effect"] == "very_high"


def test_celiac_excluded_when_no_dq2_or_dq8():
    # A DR16/DR11, DQ5/DQ7 person (no DQ2, no DQ8): celiac essentially excluded
    diplo = {"HLA-DRB1": ["HLA-DRB1*16:01", "HLA-DRB1*11:04"],
             "HLA-DQB1": ["HLA-DQB1*05:02", "HLA-DQB1*03:01"]}
    out = hs.summarize(diplo)
    assert out["celiac"]["celiac_permissive"] is False
    assert "excluded" in out["celiac"]["interpretation"].lower()
    assert "DR3: absent" in out["summary"] and "DR4: absent" in out["summary"]


def test_ms_and_narcolepsy_allele_rules():
    diplo = {"HLA-DRB1": ["HLA-DRB1*15:01", "HLA-DRB1*01:01"],
             "HLA-DQB1": ["HLA-DQB1*06:02", "HLA-DQB1*05:01"]}
    out = hs.summarize(diplo)
    diseases = {a["disease"] for a in out["disease_associations"]}
    assert "Multiple sclerosis" in diseases
    assert "Narcolepsy type 1" in diseases


def test_empty_diplotypes_is_safe():
    out = hs.summarize({})
    assert out["dr_serotypes"] == [] and out["dq_serotypes"] == []
    assert out["disease_associations"] == []
    assert out["celiac"]["celiac_permissive"] is False


def test_fixture_file_round_trips():
    import json
    fx = Path(__file__).resolve().parent / "fixtures" / "hla_synthetic.json"
    diplo = json.loads(fx.read_text())["diplotypes"]
    out = hs.summarize(diplo)
    assert "DR3" in out["dr_serotypes"] and "DR4" in out["dr_serotypes"]
    assert any(a["disease"] == "Type 1 diabetes" for a in out["disease_associations"])


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    passed = 0
    for fn in fns:
        fn()
        passed += 1
        print(f"  ok  {fn.__name__}")
    print(f"\n{passed}/{len(fns)} HLA serotype tests passed.")
