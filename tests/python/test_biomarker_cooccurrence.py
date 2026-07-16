"""Tests for vdjtools.biomarker.cooccurrence — α-β pairing and same-chain co-specificity."""
import numpy as np
import polars as pl
import pytest

from vdjtools.biomarker import cooccurrence


def _c(sid, v, j, aa):
    return dict(sample_id=sid, v_call=v, j_call=j, junction_aa=aa,
                junction_nt="ACG", duplicate_count=10, frequency=0.0)


def _paired_cohort(seed=1, n=50, n_pair=24):
    """n subjects with both chains; a planted α-β pair co-occurs in n_pair of them."""
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n):
        sid = f"S{i:02d}"
        if i < n_pair:                                     # the true pair
            rows.append(_c(sid, "TRAV1", "TRAJ1", "CAVRF"))
            rows.append(_c(sid, "TRBV1", "TRBJ1", "CASSPF"))
        rows.append(_c(sid, "TRAV2", "TRAJ2", "CAVBG"))    # background, both chains everywhere
        rows.append(_c(sid, "TRBV2", "TRBJ2", "CASSBG"))
        if rng.random() < 0.5:                             # independent extras
            rows.append(_c(sid, "TRAV3", "TRAJ3", "CAVXX"))
        if rng.random() < 0.5:
            rows.append(_c(sid, "TRBV3", "TRBJ3", "CASSYY"))
    return pl.DataFrame(rows)


def test_alpha_beta_recovers_planted_pair():
    cohort = _paired_cohort()
    r = cooccurrence(cohort, chain_a="TRA", chain_b="TRB", min_incidence=3,
                     min_cooccurrence=3, evalue=True)
    top = r.sort("p_value").row(0, named=True)
    assert (top["a_junction_aa"], top["b_junction_aa"]) == ("CAVRF", "CASSPF")
    assert top["n_ab"] == 24 and top["theta"] > 1.8 and top["q_value"] < 1e-6
    # an independent pair is not significant
    ind = r.filter((pl.col("a_junction_aa") == "CAVXX") & (pl.col("b_junction_aa") == "CASSYY"))
    assert ind.height == 0 or ind["q_value"][0] > 0.1
    assert {"expected", "e_value"} <= set(r.columns)               # evalue columns present


def test_universe_is_subjects_with_both_chains():
    cohort = _paired_cohort(n=50)
    # add 10 alpha-only subjects — they must not enter the paired universe (n stays 50)
    extra = pl.DataFrame([_c(f"A{i}", "TRAV1", "TRAJ1", "CAVRF") for i in range(10)])
    r = cooccurrence(pl.concat([cohort, extra]), chain_a="TRA", chain_b="TRB",
                     min_incidence=3, min_cooccurrence=3)
    assert r["n"][0] == 50


def test_same_chain_upper_triangle_no_self_pairs():
    # plant two TRB clonotypes that co-occur (co-specific), plus independent background
    rows = []
    for i in range(40):
        sid = f"S{i:02d}"
        if i < 20:
            rows.append(_c(sid, "TRBV1", "TRBJ1", "CASSA"))
            rows.append(_c(sid, "TRBV2", "TRBJ2", "CASSB"))       # co-occurs with CASSA
        rows.append(_c(sid, "TRBV9", "TRBJ9", "CASSBG"))
    cohort = pl.DataFrame(rows)
    r = cooccurrence(cohort, chain_a="TRB", chain_b=None, min_incidence=3, min_cooccurrence=3)
    # no self pairs
    assert (r["a_junction_aa"] != r["b_junction_aa"]).all()
    pair = r.filter(pl.col("a_junction_aa").is_in(["CASSA", "CASSB"])
                    & pl.col("b_junction_aa").is_in(["CASSA", "CASSB"]))
    assert pair.height == 1 and pair["n_ab"][0] == 20 and pair["p_value"][0] < 1e-6


def test_empty_result_stable_schema():
    # nothing clears min_cooccurrence -> empty frame, columns still present
    cohort = pl.DataFrame([_c("S0", "TRAV1", "TRAJ1", "CAVA"), _c("S0", "TRBV1", "TRBJ1", "CASSB"),
                           _c("S1", "TRAV2", "TRAJ2", "CAVC"), _c("S1", "TRBV2", "TRBJ2", "CASSD")])
    r = cooccurrence(cohort, chain_a="TRA", chain_b="TRB", min_incidence=2, min_cooccurrence=2)
    assert r.height == 0 and "theta" in r.columns and "a_junction_aa" in r.columns


def test_bad_test_arg():
    with pytest.raises(ValueError, match="test must be"):
        cooccurrence(_paired_cohort(), test="nope")
