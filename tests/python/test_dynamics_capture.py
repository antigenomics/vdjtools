"""Tests for the VDJtrack capture model + metaclonotype-grouped dynamics."""
import numpy as np
import polars as pl
import pytest

from vdjtools.dynamics import (
    SIZE_CLASSES,
    capture_paired_test,
    capture_rates,
    capture_test,
    poisson_capture,
    size_class,
)
from vdjtools.dynamics import test_metaclonotypes as run_metaclonotypes  # aliased: pytest collects test_*
from vdjtools.io.schema import COUNT, J_CALL, JUNCTION_AA, V_CALL


def _frame(cdrs, counts, v="TRBV9", j="TRBJ2-3"):
    return pl.DataFrame({JUNCTION_AA: cdrs, V_CALL: [v] * len(cdrs),
                         J_CALL: [j] * len(cdrs), COUNT: list(counts)})


# --------------------------------------------------------------------------- capture model
def test_poisson_capture_formula():
    assert poisson_capture(0.0, 1e6) == 0.0
    assert abs(poisson_capture(1e-5, 1e6) - (1 - np.exp(-10))) < 1e-12
    # monotone increasing in both frequency and depth
    assert poisson_capture(1e-6, 1e6) < poisson_capture(1e-5, 1e6) < poisson_capture(1e-4, 1e6)


def test_size_class_buckets():
    df = pl.DataFrame({COUNT: [1, 2, 3, 4, 99]}).with_columns(size_class().alias("sc"))
    assert df["sc"].to_list() == ["singleton", "doubleton", "tripleton", "large", "large"]
    assert SIZE_CLASSES == ("singleton", "doubleton", "tripleton", "large")


def test_capture_rates_shape_and_bounds():
    cdrs = [f"CASS{c}F" for c in "ACDEFGHIK"]
    pre = _frame(cdrs, [1, 1, 2, 2, 3, 3, 4, 5, 6])
    post = _frame(cdrs[:5], [10] * 5)                 # first 5 recaptured
    r = capture_rates(pre, post)
    assert set(r["size_class"].unique()) <= set(SIZE_CLASSES)
    assert (r["capture_rate"] >= 0).all() and (r["capture_rate"] <= 1).all()
    assert (r["ci_lo"] <= r["capture_rate"]).all() and (r["capture_rate"] <= r["ci_hi"]).all()
    assert (r["n_captured"] <= r["n_total"]).all()


def _cohort_rates(seed=0):
    """6 donors; a 'specific' group persists far more than background at every size."""
    rng = np.random.default_rng(seed)
    parts = []
    for d in range(6):
        n = 400
        cdrs = [f"CASS{i:04d}F".replace("0", "A").replace("1", "C").replace("2", "D")
                .replace("3", "E").replace("4", "F").replace("5", "G").replace("6", "H")
                .replace("7", "K").replace("8", "L").replace("9", "M") for i in range(n)]
        counts = rng.integers(1, 6, n)
        pre = _frame(cdrs, counts).with_columns(
            pl.Series("specific", ["yes"] * 100 + ["no"] * 300))
        p = np.where(np.arange(n) < 100, 0.9, 0.3)
        keep = rng.random(n) < p
        post = _frame([c for c, k in zip(cdrs, keep) if k],
                      [int(x) for x, k in zip(counts, keep) if k])
        parts.append(capture_rates(pre, post, group_col="specific", donor=f"D{d}"))
    return pl.concat(parts)


def test_capture_test_detects_group_effect():
    rates = _cohort_rates()
    coef = capture_test(rates, group_col="group")
    g = coef.filter(pl.col("term").str.starts_with("group:"))
    assert g.height == 1
    assert abs(g["estimate"][0]) > 0.5 and g["p_value"][0] < 0.05


def test_capture_paired_test_per_bucket():
    rates = _cohort_rates()
    paired = capture_paired_test(rates)
    assert set(paired["size_class"]) <= set(SIZE_CLASSES)
    assert (paired.drop_nulls("p_value")["p_value"] < 0.05).any()


def test_capture_test_null_when_groups_equal():
    """No group effect when both groups recapture at the same rate."""
    rng = np.random.default_rng(1)
    parts = []
    for d in range(6):
        n = 400
        cdrs = [f"CASS{chr(65 + i % 20)}{chr(65 + (i // 20) % 20)}F" for i in range(n)]
        cdrs = [f"CASS{s}" for s in
                ["".join(rng.choice(list("ACDEFGHIKLMNPQRSTVWY"), 5)) for _ in range(n)]]
        counts = rng.integers(1, 6, n)
        pre = _frame(cdrs, counts).with_columns(
            pl.Series("specific", (["yes"] * 200 + ["no"] * 200)))
        keep = rng.random(n) < 0.5                       # same rate, independent of group
        post = _frame([c for c, k in zip(cdrs, keep) if k],
                      [int(x) for x, k in zip(counts, keep) if k])
        parts.append(capture_rates(pre, post, group_col="specific", donor=f"D{d}"))
    coef = capture_test(pl.concat(parts), group_col="group")
    g = coef.filter(pl.col("term").str.starts_with("group:"))
    assert g["p_value"][0] > 0.05                         # no spurious signal


# --------------------------------------------------------------------------- metaclonotype grouping
@pytest.fixture
def _vdjmatch():
    pytest.importorskip("vdjmatch")
    try:
        from vdjtools.biomarker.metaclonotype import _require_vdjmatch
        _require_vdjmatch()
    except ImportError:
        pytest.skip("vdjmatch not importable")


def test_metaclonotypes_collapses_1mm_family(_vdjmatch):
    rng = np.random.default_rng(0)
    AA = list("ACDEFGHIKLMNPQRSTVWY")
    fam = [f"CASSL{c}PGATNEKLFF" for c in "AGSTNDEQ"]     # 1 substitution apart
    bg = ["CASR" + "".join(rng.choice(AA, 6)) + "QYF" for _ in range(200)]
    a = _frame(fam + bg, [3] * len(fam) + [100] * len(bg), v="TRBV7-2", j="TRBJ1-1")
    b = _frame(fam + bg, [300] * len(fam) + [100] * len(bg), v="TRBV7-2", j="TRBJ1-1")
    res = run_metaclonotypes(a, b, neff=None)
    fam_row = res.filter(pl.col("n_variants") >= len(fam))
    assert fam_row.height == 1                            # the family is ONE metaclonotype
    assert fam_row["dynamics"][0] == "expanded"
    assert fam_row["count_b"][0] > fam_row["count_a"][0]
