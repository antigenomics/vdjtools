"""Tests for vdjtools.stats.rarefaction — analytic interpolation/extrapolation."""
import math

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


def test_rarefaction_endpoints_and_kinds():
    df = _frame([50, 20, 10, 5, 3, 2, 1, 1])
    n = sum([50, 20, 10, 5, 3, 2, 1, 1])
    rc = stats.rarefaction(df, steps=20)
    kinds = set(rc["kind"].unique())
    assert kinds == {"interpolated", "observed", "extrapolated"}
    # x = 0 -> richness 0
    assert rc.filter(pl.col("x") == 0)["mean"][0] == 0.0
    # observed point equals Sobs at depth n
    obs = rc.filter(pl.col("kind") == "observed")
    assert obs["x"][0] == n
    assert obs["mean"][0] == 8.0


def test_rarefaction_monotone_nondecreasing():
    df = _frame([40, 20, 10, 8, 5, 3, 2, 1, 1, 1])
    rc = stats.rarefaction(df, steps=30).sort("x")
    means = rc["mean"].to_list()
    assert all(b >= a - 1e-9 for a, b in zip(means, means[1:]))


def test_rarefaction_ci_brackets_mean():
    df = _frame([10, 5, 3, 2, 1, 1])
    rc = stats.rarefaction(df, steps=15)
    assert (rc["ci_lo"] <= rc["mean"]).all()
    assert (rc["ci_hi"] >= rc["mean"]).all()


def test_rarefaction_interpolation_hand_value():
    # [1,1,1,1]: N=4, Sobs=4. Exact Hurlbert E[S] at subsample m=2 is
    #   4 * (1 - C(3,2)/C(4,2)) = 4 * (1 - 3/6) = 2.0
    # steps=5 samples x in {0,2,4,6,8}, so the x=2 interpolation point exists.
    rc = stats.rarefaction(_frame([1, 1, 1, 1]), steps=5)
    pt = rc.filter(pl.col("x") == 2)
    assert pt["kind"][0] == "interpolated"
    assert math.isclose(pt["mean"][0], 2.0, rel_tol=1e-12)


def test_rarefaction_extrapolation_hand_value():
    # [1,1,2,3,5]: n=12, Sobs=5, F1=2, F2=1, F0=0.5. At x=24 (2n) the extrapolated
    # richness equals the analytic Chao extrapolation 5 + 0.5*(1 - (2/3)^12), and
    # exceeds Sobs. steps=5 samples x in {0,6,12,18,24}, so x=24 exists.
    rc = stats.rarefaction(_frame([1, 1, 2, 3, 5]), steps=5, extrapolate_to=24)
    pt = rc.filter(pl.col("x") == 24)
    assert pt["kind"][0] == "extrapolated"
    assert math.isclose(pt["mean"][0], 5.496146326685371, rel_tol=1e-12)
    assert pt["mean"][0] > 5.0                              # > Sobs at x > N
