"""Tests for vdjtools.stats.spectratype — CDR3-length distributions."""
import math

import polars as pl

from vdjtools.io import schema as S
from vdjtools import stats


def _frame():
    df = pl.DataFrame({
        S.V_CALL: ["TRBV1", "TRBV1", "TRBV2"],
        S.J_CALL: ["TRBJ1", "TRBJ1", "TRBJ2"],
        S.JUNCTION_AA: ["CASSL", "CASSF", "CASSLGF"],       # lengths 5, 5, 7
        S.JUNCTION_NT: ["A" * 15, "A" * 15, "A" * 21],
        S.COUNT: [3, 7, 5],
    })
    return S.add_locus(S.normalize(df, recompute_freq=True))


def test_spectratype_aa_reads():
    sp = stats.spectratype(_frame(), kind="aa", weight="reads")
    got = dict(zip(sp["length"].to_list(), sp["weight"].to_list()))
    assert got == {5: 10, 7: 5}                          # 3+7 at len 5, 5 at len 7


def test_spectratype_unique():
    sp = stats.spectratype(_frame(), kind="aa", weight="unique")
    got = dict(zip(sp["length"].to_list(), sp["weight"].to_list()))
    assert got == {5: 2, 7: 1}


def test_spectratype_nt_matches_aa_times_three():
    sp = stats.spectratype(_frame(), kind="nt", weight="reads")
    got = dict(zip(sp["length"].to_list(), sp["weight"].to_list()))
    assert got == {15: 10, 21: 5}


def test_vj_spectratype_breakdown():
    vs = stats.vj_spectratype(_frame(), kind="aa", weight="reads")
    rows = {(r["v_call"], r["j_call"], r["length"]): r["weight"] for r in vs.iter_rows(named=True)}
    assert rows[("TRBV1", "TRBJ1", 5)] == 10
    assert rows[("TRBV2", "TRBJ2", 7)] == 5


def test_spectratype_freq_weight_hand_value():
    # counts [3,7,5] (total 15) at lengths 5,5,7 -> freqs .2, .4667, .3333.
    # weight='freq' sums frequency per length: len5 = 3/15 + 7/15 = 10/15,
    # len7 = 5/15. Frequencies are recomputed over the whole sample, not per length.
    sp = stats.spectratype(_frame(), kind="aa", weight="freq", by_locus=False)
    got = dict(zip(sp["length"].to_list(), sp["weight"].to_list()))
    assert math.isclose(got[5], 10 / 15, rel_tol=1e-12)
    assert math.isclose(got[7], 5 / 15, rel_tol=1e-12)
