"""Tests for vdjtools.overlap.fuzzy — edit-distance overlap delegated to vdjmatch.

The whole point of fuzzy overlap is that it catches near-variant clonotype pairs that
exact overlap misses; these tests pin that behaviour and the summary metrics. vdjmatch
(the ``overlap`` extra) is guarded with ``importorskip``.
"""
import numpy as np
import polars as pl
import pytest

from vdjtools.io import schema as S
from vdjtools import overlap as O

pytest.importorskip("vdjmatch")


def _sample(cdr3, counts, v=None, j=None):
    n = len(cdr3)
    df = pl.DataFrame({
        S.V_CALL: v or ["TRBV1"] * n, S.J_CALL: j or ["TRBJ1"] * n,
        S.JUNCTION_AA: cdr3, S.COUNT: counts,
    })
    return S.add_locus(S.normalize(df, recompute_freq=True))


def test_fuzzy_finds_pair_exact_misses():
    # a and b share only a single-substitution variant pair (…GF vs …GY).
    a = _sample(["CASSLAPGF", "CASADEF"], [10, 5])
    b = _sample(["CASSLAPGY", "CQRSTV"], [8, 3])

    # Exact CDR3 overlap sees nothing.
    assert O.overlap_metrics(a, b, key=(S.JUNCTION_AA,))["d12"] == 0

    # Fuzzy (1 substitution) recovers the pair.
    fo = O.fuzzy_overlap(a, b, scope="1,0,0,1")
    assert fo.height == 1
    row = fo.row(0, named=True)
    assert row["a_cdr3"] == "CASSLAPGF" and row["b_cdr3"] == "CASSLAPGY"
    assert row["n_subs"] == 1
    # counts/freqs joined back from the source samples.
    assert row["count_a"] == 10 and row["count_b"] == 8


def test_fuzzy_overlap_metrics_sane():
    a = _sample(["CASSLAPGF", "CASADEF"], [10, 5])
    b = _sample(["CASSLAPGY", "CQRSTV"], [8, 3])
    m = O.fuzzy_overlap_metrics(a, b, scope="1,0,0,1")
    assert m["pairs"] == 1
    assert m["frac_a_matched"] == 0.5 and m["frac_b_matched"] == 0.5
    # fuzzy_F = sqrt(freq_a_matched * freq_b_matched) = sqrt((10/15)*(8/11)) = 0.69631.
    assert np.isclose(m["fuzzy_F"], np.sqrt((10 / 15) * (8 / 11)), atol=1e-6)


def test_fuzzy_no_match_empty_and_zero_metrics():
    a = _sample(["CASSLAPGF"], [10])
    b = _sample(["CQRSTVWKY"], [8])   # far apart, no 1-sub match
    fo = O.fuzzy_overlap(a, b, scope="1,0,0,1")
    assert fo.height == 0
    m = O.fuzzy_overlap_metrics(a, b, scope="1,0,0,1")
    assert m == {"pairs": 0, "frac_a_matched": 0.0, "frac_b_matched": 0.0, "fuzzy_F": 0.0}


def test_fuzzy_exact_pair_is_zero_subs():
    # An identical clonotype is a 0-substitution fuzzy match.
    a = _sample(["CASSLAPGF"], [10])
    b = _sample(["CASSLAPGF"], [8])
    fo = O.fuzzy_overlap(a, b, scope="1,0,0,1")
    assert fo.height == 1 and fo.row(0, named=True)["n_subs"] == 0
