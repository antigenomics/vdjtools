"""Tests for vdjtools.preprocess.pool — cross-sample pooling with incidence/convergence."""
import math

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


def _row(df, cdr3):
    return df.filter(pl.col(S.JUNCTION_AA) == cdr3).to_dicts()[0]


def test_pool_aa_incidence_occurrences_convergence():
    s1 = _sample(["CASSF", "CASSF", "CASSX"], [10, 20, 70],
                 ["AAACCC", "AAACCG", "TTTGGG"])
    s2 = _sample(["CASSF", "CASSY"], [40, 60], ["AAACCC", "GGGTTT"])
    pool = pp.pool_samples([s1, s2], key="aa")

    f = _row(pool, "CASSF")
    assert f[S.COUNT] == 70                    # 10 + 20 + 40
    assert f["incidence"] == 2                 # present in both samples
    assert f["occurrences"] == 3               # 3 aggregated rows
    assert f["convergence"] == 2               # 2 distinct nt variants (AAACCC, AAACCG)
    assert math.isclose(f[S.FREQ], 70 / 200, rel_tol=1e-12)   # pool total 70+70+60

    x = _row(pool, "CASSX")
    assert (x["incidence"], x["occurrences"], x["convergence"]) == (1, 1, 1)


def test_pool_strict_splits_variants():
    s1 = _sample(["CASSF", "CASSF"], [10, 20], ["AAACCC", "AAACCG"])
    s2 = _sample(["CASSF"], [40], ["AAACCC"])
    pool = pp.pool_samples([s1, s2], key="strict")
    # strict key = nt+V+J -> AAACCC (10+40=50) and AAACCG (20) are separate clonotypes
    by_nt = dict(zip(pool[S.JUNCTION_NT].to_list(), pool[S.COUNT].to_list()))
    assert by_nt == {"AAACCC": 50, "AAACCG": 20}
    assert set(pool["convergence"].to_list()) == {1}         # nt-level -> convergence 1


def test_pool_representative_is_most_abundant():
    # aa pool of two nt variants; representative nt should be the larger one's.
    s1 = _sample(["CASSF"], [15], ["AAACCC"])
    s2 = _sample(["CASSF"], [40], ["AAACCG"])
    pool = pp.pool_samples([s1, s2], key="aa")
    assert pool[S.JUNCTION_NT].to_list() == ["AAACCG"]           # 40 > 15
