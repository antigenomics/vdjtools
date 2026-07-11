"""Tests for vdjtools.stats.rarefaction — the unified iNEXT R/E entry point.

``rarefaction`` is now the single canonical entry (default ``q=0, base="size"``
= the vdjtools-original richness curve, computed by the validated iNEXT ``q=0``
estimator). It dispatches to the size/coverage engines in
``vdjtools.stats.inext`` (oracle-pinned in ``test_inext.py``); these tests cover
the wrapper, the default, and its hand-checkable q=0 values.
"""
import math

import numpy as np
import polars as pl

from vdjtools.io import schema as S
from vdjtools import stats


def _frame(counts):
    n = len(counts)
    df = pl.DataFrame({
        S.V_CALL: ["TRBV1"] * n, S.J_CALL: ["TRBJ1"] * n,
        S.CDR3_AA: [f"CASS{i}" for i in range(n)], S.COUNT: counts,
    })
    return S.add_locus(S.normalize(df, recompute_freq=True))


def test_rarefaction_default_is_size_q0_richness():
    # Default q=0, base="size" is the vdjtools-original richness R/E curve and is
    # exactly the validated iNEXT q=0 engine.
    df = _frame([50, 20, 10, 5, 3, 2, 1, 1])
    rc = stats.rarefaction(df, se=False)
    assert set(rc["order_q"].unique().to_list()) == {0}
    assert rc.columns == [
        "order_q", "m", "method", "sample_coverage", "qD", "qD_lo", "qD_hi"]
    # Identical to a direct inext(q=0) call.
    ic = stats.inext(df, q=0, se=False)
    assert rc.equals(ic)


def test_rarefaction_methods_and_observed():
    df = _frame([50, 20, 10, 5, 3, 2, 1, 1])
    n = sum([50, 20, 10, 5, 3, 2, 1, 1])
    rc = stats.rarefaction(df, q=0, se=False)
    assert set(rc["method"].unique()) == {"rarefaction", "observed", "extrapolation"}
    # observed (m=n) point equals observed richness Sobs=8.
    obs = rc.filter(pl.col("method") == "observed")
    assert obs["m"][0] == n
    assert obs["qD"][0] == 8.0


def test_rarefaction_hill_profile():
    # q accepts a tuple -> full Hill-number profile.
    df = _frame([50, 20, 10, 5, 3, 2, 1, 1])
    rc = stats.rarefaction(df, q=(0, 1, 2), se=False)
    assert set(rc["order_q"].unique().to_list()) == {0, 1, 2}


def test_rarefaction_coverage_base():
    # base="coverage" dispatches to the coverage engine (m becomes Float64).
    df = _frame([50, 20, 10, 5, 3, 2, 1, 1])
    rc = stats.rarefaction(df, q=0, base="coverage", coverages=[0.9, 0.95], se=False)
    assert rc.columns == [
        "order_q", "sample_coverage", "m", "method", "qD", "qD_lo", "qD_hi"]
    assert rc["m"].dtype == pl.Float64
    assert np.isclose(sorted(rc["sample_coverage"].to_list()), [0.9, 0.95]).all()


def test_rarefaction_monotone_nondecreasing():
    # Richness (q=0) is non-decreasing in sampling depth m.
    df = _frame([40, 20, 10, 8, 5, 3, 2, 1, 1, 1])
    rc = stats.rarefaction(df, q=0, knots=30, se=False).sort("m")
    vals = rc["qD"].to_list()
    assert all(b >= a - 1e-9 for a, b in zip(vals, vals[1:]))


def test_rarefaction_ci_brackets_estimate():
    df = _frame([10, 5, 3, 2, 1, 1])
    rc = stats.rarefaction(df, q=0, knots=15, nboot=30, seed=0)
    assert (rc["qD_lo"] <= rc["qD"]).all()
    assert (rc["qD"] <= rc["qD_hi"]).all()
    assert (rc["qD_lo"] >= 0.0).all()


def test_rarefaction_interpolation_hand_value():
    # [1,1,1,1]: N=4, Sobs=4. Exact Hurlbert E[S] at subsample m=2 is
    #   4 * (1 - C(3,2)/C(4,2)) = 4 * (1 - 3/6) = 2.0.
    # Still valid under the unified iNEXT q=0 interpolation (MVUE == Hurlbert).
    rc = stats.rarefaction(_frame([1, 1, 1, 1]), q=0, sizes=[2], se=False)
    assert rc["method"][0] == "rarefaction"
    assert math.isclose(rc["qD"][0], 2.0, rel_tol=1e-12)


def test_rarefaction_extrapolation_hand_value():
    # [1,1,2,3,5]: n=12, Sobs=5. At m=24 (2n) the unified iNEXT q=0 extrapolation
    # (beta method: qD = D_obs + (D_asy - D_obs)*(1-(1-beta)^m*)) gives the value
    # below and exceeds Sobs. (This supersedes the old interim chaoE pin.)
    rc = stats.rarefaction(_frame([1, 1, 2, 3, 5]), q=0, sizes=[24], se=False)
    pt = rc.filter(pl.col("m") == 24)
    assert pt["method"][0] == "extrapolation"
    assert math.isclose(pt["qD"][0], 6.188008015307441, rel_tol=1e-12)
    assert pt["qD"][0] > 5.0                                # > Sobs at m > n
