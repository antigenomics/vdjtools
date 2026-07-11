"""Tests for vdjtools.preprocess.downsample — resampling and top-N selection."""
import math

import polars as pl

from vdjtools.io import schema as S
from vdjtools import preprocess as pp


def _sample(cdr3, counts, v=None, j=None):
    n = len(cdr3)
    df = pl.DataFrame({
        S.V_CALL: v or ["TRBV1"] * n, S.J_CALL: j or ["TRBJ1"] * n,
        S.CDR3_AA: cdr3, S.COUNT: counts,
    })
    return S.add_locus(S.normalize(df, recompute_freq=True))


def test_downsample_reads_total_and_reproducible():
    a = _sample(["CASSA", "CASSB", "CASSC"], [10, 30, 60])   # 100 reads
    r1 = pp.downsample(a, 50, by="reads", seed=7)
    r2 = pp.downsample(a, 50, by="reads", seed=7)
    assert r1[S.COUNT].sum() == 50                            # exact target depth
    assert r1[S.COUNT].to_list() == r2[S.COUNT].to_list()     # seeded -> reproducible
    # without replacement: no clonotype can gain reads it never had
    orig = dict(zip(a[S.CDR3_AA].to_list(), a[S.COUNT].to_list()))
    for cdr3, c in zip(r1[S.CDR3_AA].to_list(), r1[S.COUNT].to_list()):
        assert c <= orig[cdr3]
    assert math.isclose(r1[S.FREQ].sum(), 1.0, rel_tol=1e-12)


def test_downsample_reads_different_seed_differs():
    a = _sample([f"C{i}" for i in range(20)], [50] * 20)     # 1000 reads
    r1 = pp.downsample(a, 200, by="reads", seed=1)[S.COUNT].to_list()
    r2 = pp.downsample(a, 200, by="reads", seed=2)[S.COUNT].to_list()
    assert r1 != r2


def test_downsample_reads_guard_size_ge_total():
    a = _sample(["CASSA", "CASSB"], [10, 30])                 # 40 reads
    out = pp.downsample(a, 40, by="reads")                    # size >= total -> unchanged
    assert out[S.COUNT].to_list() == [10, 30]


def test_downsample_clones_uniform_without_replacement():
    a = _sample(["CASSA", "CASSB", "CASSC", "CASSD"], [10, 20, 30, 40])
    out = pp.downsample(a, 2, by="clones", seed=3)
    assert out.height == 2                                    # exactly size clonotypes
    # counts are preserved (not resampled) in clones mode
    orig = dict(zip(a[S.CDR3_AA].to_list(), a[S.COUNT].to_list()))
    for cdr3, c in zip(out[S.CDR3_AA].to_list(), out[S.COUNT].to_list()):
        assert orig[cdr3] == c
    assert pp.downsample(a, 2, by="clones", seed=3)[S.CDR3_AA].to_list() == \
        out[S.CDR3_AA].to_list()                              # reproducible


def test_select_top_by_count_renormalizes():
    a = _sample(["CASSA", "CASSB", "CASSC"], [10, 30, 60])
    top = pp.select_top(a, 2)
    assert top[S.CDR3_AA].to_list() == ["CASSC", "CASSB"]     # descending by count
    assert top[S.COUNT].to_list() == [60, 30]
    assert math.isclose(top[S.FREQ].sum(), 1.0, rel_tol=1e-12)
    assert math.isclose(top[S.FREQ].to_list()[0], 60 / 90, rel_tol=1e-12)


def test_select_top_no_renormalize_preserves_freq():
    a = _sample(["CASSA", "CASSB", "CASSC"], [10, 30, 60])
    top = pp.select_top(a, 2, renormalize=False)
    assert math.isclose(top[S.FREQ].to_list()[0], 0.6, rel_tol=1e-12)   # 60/100 preserved
