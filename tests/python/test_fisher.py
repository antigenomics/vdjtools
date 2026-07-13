"""Incidence-based Fisher association (Emerson 2017) — oracle pins + option coverage."""
from __future__ import annotations

import polars as pl
import pytest
from scipy.stats import fisher_exact

from vdjtools.biomarker import fisher_association
from vdjtools.io import schema as S

_COLS = ["sample_id", "v_call", "j_call", "junction_aa", "duplicate_count"]


def _cohort():
    """20 subjects (10 CMV+, 10 CMV−) with planted enriched/depleted/background features."""
    pos = [f"pos{i}" for i in range(10)]
    neg = [f"neg{i}" for i in range(10)]
    rows = []
    # X: enriched — 8/10 CMV+, 1/10 CMV−  (allele suffix present, must be stripped)
    rows += [(s, "TRBV1*01", "TRBJ1", "CASSXF", 10) for s in pos[:8] + neg[:1]]
    # Z: depleted — 1/10 CMV+, 8/10 CMV−
    rows += [(s, "TRBV4", "TRBJ4", "CASSZF", 10) for s in pos[:1] + neg[:8]]
    # B: ubiquitous background — all 20 subjects, no association
    rows += [(s, "TRBV2", "TRBJ2", "CASSBG", 5) for s in pos + neg]
    # R: private singleton — one subject, filtered out by min_incidence
    rows += [("pos0", "TRBV3", "TRBJ3", "CASSRARE", 1)]
    df = pl.DataFrame(rows, schema=_COLS, orient="row")
    pheno = pl.DataFrame({"sample_id": pos + neg, "cmv": [True] * 10 + [False] * 10})
    return df, pheno


def test_recovers_enriched_feature():
    df, pheno = _cohort()
    res = fisher_association(df, pheno, pheno_col="cmv")  # default: full key, one-tailed greater
    x = res.filter(pl.col("junction_aa") == "CASSXF").row(0, named=True)
    assert (x["incidence"], x["n_pos_present"], x["n_neg_present"]) == (9, 8, 1)
    assert x["direction"] == "enriched"
    assert x["log2_or"] > 0 and x["p_value"] < 0.05
    # ubiquitous background is not associated; private singleton is filtered.
    bg = res.filter(pl.col("junction_aa") == "CASSBG").row(0, named=True)
    assert bg["p_value"] > 0.5
    assert res.filter(pl.col("junction_aa") == "CASSRARE").height == 0


def test_hypergeom_tail_matches_scipy_fisher_exact():
    """The vectorised hypergeometric one-tailed p == scipy.stats.fisher_exact, both tails."""
    df, pheno = _cohort()
    for alt in ("greater", "less"):
        res = fisher_association(df, pheno, pheno_col="cmv", alternative=alt)
        for cdr3 in ("CASSXF", "CASSZF", "CASSBG"):
            r = res.filter(pl.col("junction_aa") == cdr3).row(0, named=True)
            a, b = r["n_pos_present"], r["n_neg_present"]
            c, d = r["n_pos"] - a, r["n_neg"] - b
            _, p_sp = fisher_exact([[a, b], [c, d]], alternative=alt)
            assert r["p_value"] == pytest.approx(p_sp, rel=1e-9), (cdr3, alt)


def test_depleted_direction_and_less_tail():
    df, pheno = _cohort()
    res = fisher_association(df, pheno, pheno_col="cmv", alternative="less")
    z = res.filter(pl.col("junction_aa") == "CASSZF").row(0, named=True)
    assert z["direction"] == "depleted" and z["log2_or"] < 0 and z["p_value"] < 0.05


def test_two_sided_is_doubled_one_tail():
    df, pheno = _cohort()
    g = fisher_association(df, pheno, pheno_col="cmv", alternative="greater")
    two = fisher_association(df, pheno, pheno_col="cmv", alternative="two-sided")
    xg = g.filter(pl.col("junction_aa") == "CASSXF")["p_value"].item()
    xt = two.filter(pl.col("junction_aa") == "CASSXF")["p_value"].item()
    assert xt == pytest.approx(min(1.0, 2 * xg))  # enriched: greater is the smaller tail


def test_vj_match_requirement_via_key():
    """The key selects the V/J match requirement; CDR3-only drops the gene columns."""
    df, pheno = _cohort()
    full = fisher_association(df, pheno, pheno_col="cmv")
    assert {"v_call", "j_call"} <= set(full.columns)
    cdr_only = fisher_association(df, pheno, pheno_col="cmv", key=(S.JUNCTION_AA,))
    assert "v_call" not in cdr_only.columns and "j_call" not in cdr_only.columns
    assert {"junction_aa", "incidence", "p_value", "q_value"} <= set(cdr_only.columns)


def test_min_incidence_and_lazyframe_input():
    df, pheno = _cohort()
    res = fisher_association(df.lazy(), pheno, pheno_col="cmv", min_incidence=20)
    assert set(res["junction_aa"]) == {"CASSBG"}  # only the all-subject background survives


def test_unknown_phenotype_labels_excluded():
    """A null phenotype label removes that subject from both classes (n_tot shrinks)."""
    df, pheno = _cohort()
    pheno = pheno.with_columns(
        pl.when(pl.col("sample_id") == "neg9").then(None).otherwise(pl.col("cmv")).alias("cmv"))
    res = fisher_association(df, pheno, pheno_col="cmv")
    assert res["n_pos"][0] + res["n_neg"][0] == 19  # neg9 dropped


def test_empty_result_matches_nonempty_schema():
    """A too-high min_incidence yields an empty frame with the SAME columns as a non-empty
    run, so per-partition results concat without a polars ShapeError."""
    df, pheno = _cohort()
    full = fisher_association(df, pheno, pheno_col="cmv")
    empty = fisher_association(df, pheno, pheno_col="cmv", min_incidence=999)
    assert empty.height == 0 and empty.columns == full.columns
    assert pl.concat([full, empty]).height == full.height  # concat-safe


def test_bad_args():
    df, pheno = _cohort()
    with pytest.raises(ValueError):
        fisher_association(df, pheno, pheno_col="cmv", alternative="nope")
    with pytest.raises(ValueError):
        fisher_association(df, pheno, pheno_col="cmv", key=(S.V_CALL,))  # no cdr3_aa
