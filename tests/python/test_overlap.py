"""Tests for vdjtools.overlap.metrics — exact-match D/F/F2/R overlap."""
import math

import polars as pl

from vdjtools.io import schema as S
from vdjtools import overlap as O


def _sample(cdr3, counts, v=None, j=None):
    n = len(cdr3)
    df = pl.DataFrame({
        S.V_CALL: v or ["TRBV1"] * n, S.J_CALL: j or ["TRBJ1"] * n,
        S.CDR3_AA: cdr3, S.COUNT: counts,
    })
    return S.add_locus(S.normalize(df, recompute_freq=True))


def test_overlap_metrics_hand_values():
    a = _sample(["CASSL", "CASSF", "CASSX"], [10, 30, 60])   # freqs .1 .3 .6
    b = _sample(["CASSL", "CASSF", "CASSY"], [40, 40, 20])   # freqs .4 .4 .2
    shared, m = O.overlap_pair(a, b)
    assert m["d1"] == 3 and m["d2"] == 3 and m["d12"] == 2
    assert math.isclose(m["D"], 2 / 9, rel_tol=1e-12)
    assert math.isclose(m["F"], math.sqrt(0.4 * 0.8), rel_tol=1e-12)
    assert math.isclose(m["F2"], 0.2 + math.sqrt(0.12), rel_tol=1e-12)
    assert shared.height == 2
    assert set(shared[S.CDR3_AA].to_list()) == {"CASSL", "CASSF"}


def test_overlap_no_shared():
    a = _sample(["CASSA"], [1])
    b = _sample(["CASSB"], [1])
    m = O.overlap_metrics(a, b)
    assert m["d12"] == 0 and m["D"] == 0.0 and m["F"] == 0.0 and m["F2"] == 0.0
    assert m["R"] == 0.0


def test_overlap_correlation_positive_on_log_freq():
    # concordant frequencies over >2 shared clonotypes -> positive Pearson R
    a = _sample(["A", "B", "C", "D"], [1, 10, 100, 1000])
    b = _sample(["A", "B", "C", "D"], [2, 20, 200, 2000])
    m = O.overlap_metrics(a, b, key=(S.CDR3_AA,))
    assert m["d12"] == 4
    assert m["R"] > 0.99                                   # perfectly rank/scale concordant


def test_overlap_cdr3_only_key():
    a = _sample(["CASSL"], [1], v=["TRBV1"], j=["TRBJ1"])
    b = _sample(["CASSL"], [1], v=["TRBV9"], j=["TRBJ9"])   # same cdr3, different V/J
    assert O.overlap_metrics(a, b)["d12"] == 0             # default key includes V/J
    assert O.overlap_metrics(a, b, key=(S.CDR3_AA,))["d12"] == 1
