"""Tests for vdjtools.stats.rarefaction — analytic interpolation/extrapolation."""
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
