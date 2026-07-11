"""Tests for vdjtools.preprocess.correct — frequency-based error correction.

The abundance-ratio merge core is tested directly (no seqtree needed); the
end-to-end path with the seqtree neighbour search is guarded by ``importorskip``.
"""
import math

import numpy as np
import polars as pl
import pytest

from vdjtools.io import schema as S
from vdjtools import preprocess as pp
from vdjtools.preprocess.correct import _corrected_counts


def _sample(cdr3, counts, nt, v=None, j=None):
    n = len(cdr3)
    df = pl.DataFrame({
        S.V_CALL: v or ["TRBV1"] * n, S.J_CALL: j or ["TRBJ1"] * n,
        S.CDR3_AA: cdr3, S.CDR3_NT: nt, S.COUNT: counts,
    })
    return S.add_locus(S.normalize(df, recompute_freq=True))


def test_merge_core_parent_absorbs_child():
    # counts [1000, 10], 1 mismatch, ratio 0.05 -> child 10 < 0.05*1000 = 50 -> merge.
    counts = np.array([1000, 10])
    neighbours = [[(1, 1)], [(0, 1)]]
    lrt = -math.log10(0.05)
    out = _corrected_counts(counts, neighbours, lrt)
    assert out.tolist() == [1010, 0]           # parent gains child, child removed


def test_merge_core_ambiguous_kept():
    # child 100 is NOT below 0.05*1000 = 50 -> neither merges, both survive.
    counts = np.array([1000, 100])
    neighbours = [[(1, 1)], [(0, 1)]]
    lrt = -math.log10(0.05)
    out = _corrected_counts(counts, neighbours, lrt)
    assert out.tolist() == [1000, 100]


def test_merge_core_two_mismatch_threshold():
    # 2 mismatches: threshold ratio^2 = 0.0025 -> child must be < 2.5 to merge.
    lrt = -math.log10(0.05)
    assert _corrected_counts(np.array([1000, 2]), [[(1, 2)], [(0, 2)]], lrt).tolist() == [1002, 0]
    assert _corrected_counts(np.array([1000, 5]), [[(1, 2)], [(0, 2)]], lrt).tolist() == [1000, 5]


def test_correct_end_to_end_merges_error():
    pytest.importorskip("seqtree")
    s = _sample(["CASSP", "CASSE"], [1000, 10], ["ACGTACGT", "ACGTACGA"])  # 1 mismatch
    out = pp.correct(s, max_mismatches=2, ratio=0.05, same_vj=True)
    assert out.height == 1
    assert out[S.COUNT].to_list() == [1010]
    assert out[S.CDR3_NT].to_list() == ["ACGTACGT"]           # parent survives
    assert math.isclose(out[S.FREQ].to_list()[0], 1.0, rel_tol=1e-12)


def test_correct_same_vj_isolates_blocks():
    pytest.importorskip("seqtree")
    # same nt-distance error but different J -> not compared when same_vj=True.
    s = _sample(["CASSP", "CASSE"], [1000, 10], ["ACGTACGT", "ACGTACGA"],
                j=["TRBJ1", "TRBJ2"])
    out = pp.correct(s, same_vj=True)
    assert out.height == 2                                    # segments differ -> no merge
