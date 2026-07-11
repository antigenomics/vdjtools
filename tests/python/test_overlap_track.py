"""Tests for vdjtools.overlap.track — per-sample clonotype frequency trajectories."""
import polars as pl

from vdjtools.io import schema as S
from vdjtools import overlap as O


def _sample(cdr3, counts):
    n = len(cdr3)
    df = pl.DataFrame({S.V_CALL: ["TRBV1"] * n, S.J_CALL: ["TRBJ1"] * n,
                       S.CDR3_AA: cdr3, S.COUNT: counts})
    return S.add_locus(S.normalize(df, recompute_freq=True))


def test_track_pivots_frequency_across_samples():
    t0 = _sample(["CASSA", "CASSB"], [90, 10])       # freqs .9 .1
    t1 = _sample(["CASSA", "CASSC"], [50, 50])       # freqs .5 .5
    tr = O.track_clonotypes({"t0": t0, "t1": t1}, order=["t0", "t1"])

    assert tr.columns == [S.CDR3_AA, S.V_CALL, S.J_CALL,
                          "freq_t0", "freq_t1", "freq_sum"]
    # union of clonotypes present in >=1 sample.
    assert set(tr[S.CDR3_AA].to_list()) == {"CASSA", "CASSB", "CASSC"}
    row = {r[S.CDR3_AA]: r for r in tr.iter_rows(named=True)}
    # CASSA present in both; CASSB only in t0 (0.0 in t1); CASSC only in t1.
    assert row["CASSA"]["freq_t0"] == 0.9 and row["CASSA"]["freq_t1"] == 0.5
    assert row["CASSB"]["freq_t1"] == 0.0
    assert row["CASSC"]["freq_t0"] == 0.0
    # sorted by summed frequency descending: CASSA (1.4) first.
    assert tr[S.CDR3_AA][0] == "CASSA"


def test_track_top_n_and_list_input():
    samples = [_sample(["A", "B", "C"], [100, 10, 1]),
               _sample(["A", "B", "D"], [100, 10, 1])]
    tr = O.track_clonotypes(samples, top=2, key=(S.CDR3_AA,))
    assert tr.height == 2                             # top-2 by summed frequency
    assert tr[S.CDR3_AA].to_list()[0] == "A"          # most abundant across samples
    # list input names samples "0","1".
    assert "freq_0" in tr.columns and "freq_1" in tr.columns


def test_track_default_order():
    tr = O.track_clonotypes({"s1": _sample(["A"], [1]), "s2": _sample(["B"], [1])})
    assert tr.columns[-3:] == ["freq_s1", "freq_s2", "freq_sum"]


def test_track_empty_inputs():
    """No usable samples (empty dict, or an order that filters everything out) return
    an empty frame with the key schema — the out-is-None early-return branch."""
    for empty in (O.track_clonotypes({}),
                  O.track_clonotypes({"s1": _sample(["A"], [1])}, order=["missing"])):
        assert empty.height == 0
        assert empty.columns == [S.CDR3_AA, S.V_CALL, S.J_CALL]
