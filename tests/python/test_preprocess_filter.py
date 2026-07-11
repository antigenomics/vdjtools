"""Tests for vdjtools.preprocess.filter — functional / frequency / segment / by-sample."""
import math

import polars as pl

from vdjtools.io import schema as S
from vdjtools import preprocess as pp


def _sample(cdr3, counts, v=None, j=None, d=None):
    n = len(cdr3)
    df = pl.DataFrame({
        S.V_CALL: v or ["TRBV1"] * n, S.D_CALL: d or ["TRBD1"] * n,
        S.J_CALL: j or ["TRBJ1"] * n, S.CDR3_AA: cdr3, S.COUNT: counts,
    })
    return S.add_locus(S.normalize(df, recompute_freq=True))


def test_filter_functional_drops_stop_and_oof():
    # CA*SF has a stop codon; CASaL carries a lowercase nucleotide (legacy out-of-frame
    # marker [atgc#~_?]); CASSL is clean coding.
    a = _sample(["CASSL", "CA*SF", "CASaL"], [10, 20, 30])
    coding = pp.filter_functional(a, keep="coding")
    assert coding[S.CDR3_AA].to_list() == ["CASSL"]
    assert math.isclose(coding[S.FREQ].sum(), 1.0, rel_tol=1e-12)
    noncoding = pp.filter_functional(a, keep="noncoding")
    assert set(noncoding[S.CDR3_AA].to_list()) == {"CA*SF", "CASaL"}


def test_filter_frequency_min_freq():
    a = _sample(["A", "B", "C"], [1, 9, 90])                 # freqs .01 .09 .90
    out = pp.filter_frequency(a, min_freq=0.05)
    assert set(out[S.CDR3_AA].to_list()) == {"B", "C"}       # drops the .01 clone


def test_filter_frequency_top_quantile():
    # counts 60,30,10 -> freqs .6,.3,.1. top_quantile=0.25 keeps only clones whose
    # cumulative read mass stays <= 0.25 -> just the top (0.6 > 0.25 excludes rest,
    # but the top itself: cumulative 0.6 > 0.25 -> even it fails? legacy keeps while
    # cumulative <= q; here 0.6 > 0.25 so nothing passes). Use q=0.65 -> keep top only.
    a = _sample(["T", "M", "S"], [60, 30, 10])
    q = pp.filter_frequency(a, top_quantile=0.65)
    assert q[S.CDR3_AA].to_list() == ["T"]                   # cum .6<=.65, +.3 -> .9>.65
    q2 = pp.filter_frequency(a, top_quantile=0.95)
    assert q2[S.CDR3_AA].to_list() == ["T", "M"]             # .6,.9<=.95 ; +.1 excluded


def test_filter_segment_keep_and_remove():
    a = _sample(["A", "B", "C"], [1, 1, 1],
                v=["TRBV12-3*01", "TRBV20-1*01", "TRBV12-4*01"])
    keep = pp.filter_segment(a, v=["TRBV12"])               # prefix, allele-insensitive
    assert set(keep[S.V_CALL].to_list()) == {"TRBV12-3*01", "TRBV12-4*01"}
    remove = pp.filter_segment(a, v=["TRBV12"], keep=False)
    assert remove[S.V_CALL].to_list() == ["TRBV20-1*01"]


def test_filter_segment_vj_conjunction():
    a = _sample(["A", "B"], [1, 1], v=["TRBV1", "TRBV1"], j=["TRBJ1-1", "TRBJ2-1"])
    out = pp.filter_segment(a, v=["TRBV1"], j=["TRBJ2"])    # V AND J must match
    assert out[S.J_CALL].to_list() == ["TRBJ2-1"]


def test_filter_by_sample_keep_and_remove():
    a = _sample(["CASSL", "CASSF", "CASSX"], [1, 2, 3])
    other = _sample(["CASSL", "CASSF"], [5, 5])
    keep = pp.filter_by_sample(a, other, keep=True)
    assert set(keep[S.CDR3_AA].to_list()) == {"CASSL", "CASSF"}
    remove = pp.filter_by_sample(a, other, keep=False)
    assert remove[S.CDR3_AA].to_list() == ["CASSX"]


def test_filter_by_sample_key_includes_vj():
    a = _sample(["CASSL"], [1], v=["TRBV1"], j=["TRBJ1"])
    other = _sample(["CASSL"], [1], v=["TRBV9"], j=["TRBJ1"])   # same aa, different V
    assert pp.filter_by_sample(a, other, keep=True).height == 0  # default key has V/J
    assert pp.filter_by_sample(a, other, keep=True, key=(S.CDR3_AA,)).height == 1
