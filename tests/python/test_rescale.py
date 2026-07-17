"""V/J usage rescaling — swap a model's protocol-dependent usage, keep the junction model.

Usage is protocol-dependent (5'RACE vs DNA multiplex amplify different V's); the junction model
is not. These pin that rescale_usage moves ONLY the usage.
"""
from __future__ import annotations

import polars as pl
import pytest

from vdjtools.model import load_bundled, rescale_usage


def _sample(v_calls, j_calls=None) -> pl.DataFrame:
    j_calls = j_calls or ["TRBJ2-7*01"] * len(v_calls)
    return pl.DataFrame({
        "junction_aa": [f"CASS{i}EQYF" for i in range(len(v_calls))],
        "junction_nt": ["ACG"] * len(v_calls),
        "v_call": v_calls, "j_call": j_calls,
        "duplicate_count": [1] * len(v_calls), "frequency": [1.0 / len(v_calls)] * len(v_calls),
    })


def _gene_mass(model, table, col):
    t = model.tables[table].with_columns(pl.col(col).str.split("*").list.first().alias("g"))
    return {r["g"]: r["p"] for r in t.group_by("g").agg(pl.col("p").sum()).iter_rows(named=True)}


def test_rescale_sets_v_usage_to_the_sample():
    m = load_bundled("TRB", "olga")
    s = _sample(["TRBV19*01"] * 75 + ["TRBV20-1*01"] * 25)
    r = rescale_usage(m, s)
    g = _gene_mass(r, "v_choice", "v_allele")
    assert g["TRBV19"] == pytest.approx(0.75)
    assert g["TRBV20-1"] == pytest.approx(0.25)
    assert sum(g.values()) == pytest.approx(1.0)
    assert g.get("TRBV5-1", 0.0) == 0.0          # absent from the sample -> no mass


def test_rescale_leaves_the_junction_model_untouched():
    """The whole point: usage moves, the recombination machinery does not."""
    m = load_bundled("TRB", "olga")
    r = rescale_usage(m, _sample(["TRBV19*01"] * 10 + ["TRBV20-1*01"] * 10))
    for ev in ("v_3_del", "j_5_del", "d_del", "vd_ins", "dj_ins", "vd_dinucl", "dj_dinucl",
               "d_gene", "n_d"):
        assert r.tables[ev].equals(m.tables[ev]), f"{ev} changed — rescale touched the junction model"


def test_rescale_preserves_within_gene_allele_split():
    """A gene's new mass is split across its alleles in the model's existing proportions.

    The sample cannot resolve alleles (mismapping on short reads), so the model's own split is
    the best available and must survive.
    """
    m = load_bundled("TRB", "olga")
    v = m.tables["v_choice"].with_columns(pl.col("v_allele").str.split("*").list.first().alias("g"))
    multi = [r["g"] for r in v.group_by("g").len().iter_rows(named=True) if r["len"] > 1]
    gene = next(g for g in multi if v.filter(pl.col("g") == g)["p"].sum() > 0)
    before = v.filter(pl.col("g") == gene)
    ratio_before = (before["p"] / before["p"].sum()).to_list()

    r = rescale_usage(m, _sample([f"{gene}*01"] * 10))
    after = r.tables["v_choice"].with_columns(pl.col("v_allele").str.split("*").list.first().alias("g")).filter(pl.col("g") == gene)
    assert after["p"].sum() == pytest.approx(1.0)
    assert (after["p"] / after["p"].sum()).to_list() == pytest.approx(ratio_before)


def test_rescale_j_and_selective_flags():
    m = load_bundled("TRB", "olga")
    s = _sample(["TRBV19*01"] * 10, ["TRBJ2-7*01"] * 6 + ["TRBJ1-1*01"] * 4)
    r = rescale_usage(m, s)
    gj = _gene_mass(r, "j_choice", "j_allele")
    assert gj["TRBJ2-7"] == pytest.approx(0.6)
    assert gj["TRBJ1-1"] == pytest.approx(0.4)
    # v=False leaves v_choice alone
    assert rescale_usage(m, s, v=False).tables["v_choice"].equals(m.tables["v_choice"])
    assert rescale_usage(m, s, j=False).tables["j_choice"].equals(m.tables["j_choice"])


def test_rescale_rejects_an_unusable_sample():
    m = load_bundled("TRB", "olga")
    # An all-null V column is genuinely unusable.
    s = _sample(["TRBV6-2*01"] * 3).with_columns(pl.lit(None, dtype=pl.String).alias("v_call"))
    with pytest.raises(ValueError, match="no usable"):
        rescale_usage(m, s)


def test_ambiguous_only_sample_is_usable_via_fractional_split():
    """An all-ambiguous sample is NOT unusable — the ties still carry gene information."""
    m = load_bundled("TRB", "olga")
    r = rescale_usage(m, _sample(["TRBV6-2*01,TRBV6-3*01"] * 5))    # was: raised; now: 0.5/0.5
    g = (r.tables["v_choice"].with_columns(pl.col("v_allele").str.split("*").list.first().alias("g"))
         .group_by("g").agg(pl.col("p").sum()))
    mass = {row["g"]: row["p"] for row in g.iter_rows(named=True)}
    assert mass.get("TRBV6-2", 0) == pytest.approx(0.5)
    assert mass.get("TRBV6-3", 0) == pytest.approx(0.5)


def test_ambiguous_calls_are_split_fractionally_not_dropped():
    """A comma tie is allocated 1/k per named gene, not discarded (duplicated loci tie constantly)."""
    from vdjtools.model.rescale import _empirical

    # 4 clonotypes: 2 unambiguous IGHV3-23, 2 ties naming IGHV3-23 + IGHV3-23D.
    s = pl.DataFrame({
        "junction_aa": [f"CAS{i}F" for i in range(4)],
        "junction_nt": ["ACG"] * 4,
        "v_call": ["IGHV3-23*01", "IGHV3-23*01",
                   "IGHV3-23*01,IGHV3-23D*01", "IGHV3-23*01,IGHV3-23D*01"],
        "j_call": ["IGHJ4*02"] * 4,
        "duplicate_count": [1] * 4, "frequency": [0.25] * 4,
    })
    p = _empirical(s, "v_call")
    # votes: IGHV3-23 = 2 + 0.5 + 0.5 = 3 ; IGHV3-23D = 0.5 + 0.5 = 1 ; total 4
    assert p["IGHV3-23"] == pytest.approx(0.75)
    assert p["IGHV3-23D"] == pytest.approx(0.25)
    assert sum(p.values()) == pytest.approx(1.0)
    # dropping the ties (the old behaviour) would have given IGHV3-23=1.0 and hidden IGHV3-23D
    assert "IGHV3-23D" in p, "the duplicated-locus paralog was dropped"
