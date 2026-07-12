"""Tests for vdjtools.preprocess.decontaminate — cross-sample ratio filtering."""
import polars as pl

from vdjtools.io import schema as S
from vdjtools import preprocess as pp


def _sample(cdr3, counts, nt, v=None, j=None):
    n = len(cdr3)
    df = pl.DataFrame({
        S.V_CALL: v or ["TRBV1"] * n, S.J_CALL: j or ["TRBJ1"] * n,
        S.JUNCTION_AA: cdr3, S.JUNCTION_NT: nt, S.COUNT: counts,
    })
    return S.add_locus(S.normalize(df, recompute_freq=True))


def test_decontaminate_reads_removes_dominated_clone():
    main = _sample(["CASSF", "CASSK"], [5, 100], ["AAACCC", "GGGTTT"])
    other = _sample(["CASSF"], [500], ["AAACCC"])            # 500 >= 5*20 -> remove CASSF
    out = pp.decontaminate(main, [other], ratio=20.0, by="reads")
    assert out[S.JUNCTION_NT].to_list() == ["GGGTTT"]            # only the un-dominated clone


def test_decontaminate_reads_keeps_below_ratio():
    main = _sample(["CASSF"], [100], ["AAACCC"])
    other = _sample(["CASSF"], [500], ["AAACCC"])            # 500 < 100*20 = 2000 -> keep
    assert pp.decontaminate(main, [other], ratio=20.0, by="reads").height == 1


def test_decontaminate_no_match_keeps_all():
    main = _sample(["CASSF"], [5], ["AAACCC"])
    other = _sample(["CASSK"], [999], ["GGGTTT"])            # different key -> no effect
    assert pp.decontaminate(main, [other], ratio=20.0).height == 1


def test_decontaminate_freq_mode():
    # main CASSF freq = 5/105 ~= 0.0476; other CASSF freq = 1.0. 1.0 >= 0.0476*20=0.952
    main = _sample(["CASSF", "CASSK"], [5, 100], ["AAACCC", "GGGTTT"])
    other = _sample(["CASSF"], [500], ["AAACCC"])
    out = pp.decontaminate(main, [other], ratio=20.0, by="freq")
    assert out[S.JUNCTION_NT].to_list() == ["GGGTTT"]


def test_decontaminate_key_includes_vj():
    main = _sample(["CASSF"], [5], ["AAACCC"], v=["TRBV1"])
    other = _sample(["CASSF"], [500], ["AAACCC"], v=["TRBV9"])   # same nt/aa, diff V
    assert pp.decontaminate(main, [other], ratio=20.0, by="reads").height == 1  # not matched
