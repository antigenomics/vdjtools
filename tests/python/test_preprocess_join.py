"""Tests for vdjtools.preprocess.join — geometric-mean joint table."""
import math

import polars as pl

from vdjtools.io import schema as S
from vdjtools.preprocess.join import JITTER
from vdjtools import preprocess as pp


def _sample(cdr3, counts, v=None, j=None):
    n = len(cdr3)
    df = pl.DataFrame({
        S.V_CALL: v or ["TRBV1"] * n, S.J_CALL: j or ["TRBJ1"] * n,
        S.CDR3_AA: cdr3, S.COUNT: counts,
    })
    return S.add_locus(S.normalize(df, recompute_freq=True))


def _row(df, cdr3):
    return df.filter(pl.col(S.CDR3_AA) == cdr3).to_dicts()[0]


def test_join_geometric_mean_and_normalized_count():
    # s1 freqs: A=0.1 B=0.9 ; s2 freqs: A=0.4 B=0.6
    s1 = _sample(["A", "B"], [10, 90])
    s2 = _sample(["A", "B"], [40, 60])
    j = pp.join_samples([s1, s2], key="aa", min_samples=2)

    base_a = math.sqrt((0.1 + JITTER) * (0.4 + JITTER))
    base_b = math.sqrt((0.9 + JITTER) * (0.6 + JITTER))
    total = base_a + base_b
    ra, rb = _row(j, "A"), _row(j, "B")
    assert math.isclose(ra[S.FREQ], base_a / total, rel_tol=1e-9)     # geometric mean, normalised
    assert math.isclose(rb[S.FREQ], base_b / total, rel_tol=1e-9)
    assert math.isclose(ra[S.FREQ] + rb[S.FREQ], 1.0, rel_tol=1e-9)   # joint freqs sum to 1
    # normalised count: smallest base -> 1, other floored to base/min_base
    assert ra[S.COUNT] == 1                                           # A is the smaller base
    assert rb[S.COUNT] == math.floor(base_b / base_a)
    assert ra["incidence"] == 2 and rb["incidence"] == 2
    assert math.isclose(ra["freq_0"], 0.1, rel_tol=1e-12)
    assert math.isclose(ra["freq_1"], 0.4, rel_tol=1e-12)


def test_join_min_samples_filters_private_clones():
    s1 = _sample(["A", "B", "P"], [10, 40, 50])              # P is private to s1
    s2 = _sample(["A", "B"], [30, 70])
    j = pp.join_samples([s1, s2], key="aa", min_samples=2)
    assert set(j[S.CDR3_AA].to_list()) == {"A", "B"}         # P dropped (incidence 1)


def test_join_min_samples_one_keeps_private():
    s1 = _sample(["A", "P"], [50, 50])
    s2 = _sample(["A"], [100])
    j = pp.join_samples([s1, s2], key="aa", min_samples=1)
    assert set(j[S.CDR3_AA].to_list()) == {"A", "P"}


def test_join_min_samples_above_n_empty_result():
    # Only 2 samples but min_samples=3 -> nothing can pass; the empty-result path must
    # still return the full joint schema (key + per-sample freq/count + incidence/freq/count).
    s1 = _sample(["A", "B"], [10, 90])
    s2 = _sample(["A", "B"], [40, 60])
    j = pp.join_samples([s1, s2], key="aa", min_samples=3)
    assert j.height == 0
    assert set(j.columns) == {S.CDR3_AA, "freq_0", "freq_1", "count_0", "count_1",
                              "incidence", S.FREQ, S.COUNT}


def test_join_named_columns():
    s1 = _sample(["A", "B"], [10, 90])
    s2 = _sample(["A", "B"], [40, 60])
    j = pp.join_samples([s1, s2], key="aa", names=["donor1", "donor2"])
    assert "freq_donor1" in j.columns and "count_donor2" in j.columns
