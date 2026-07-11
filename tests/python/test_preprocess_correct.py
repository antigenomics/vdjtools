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


def test_merge_core_chain_read_non_conservation():
    # 3-node chain 1000 -> 40 -> 1, all decisions on the ORIGINAL counts:
    #   node 0 (1000): child 40 < 0.05*1000=50 -> absorbs 40 -> 1040.
    #   node 1 (40): sees bigger parent 1000 (40 < 50) -> removed (0), returns
    #     immediately WITHOUT absorbing its own child 1.
    #   node 2 (1): sees bigger parent 40 (1 < 2) -> removed (0).
    # The grandchild's read (1) is dropped rather than passed up: reads are not
    # conserved, and the result is independent of node order.
    lrt = -math.log10(0.05)
    r = _corrected_counts(np.array([1000, 40, 1]),
                          [[(1, 1)], [(0, 1), (2, 1)], [(1, 1)]], lrt)
    assert r.tolist() == [1040, 0, 0]


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


def test_correct_default_is_segment_agnostic():
    pytest.importorskip("seqtree")
    # Default (same_vj=False, legacy fidelity): a 1-substitution error is merged even
    # across DIFFERENT V segments (TRBV1 vs TRBV2). counts [1000, 10] -> single 1010.
    s = _sample(["CASSP", "CASSE"], [1000, 10], ["ACGTACGT", "ACGTACGA"],
                v=["TRBV1", "TRBV2"])
    out = pp.correct(s)                                       # no same_vj -> default False
    assert out.height == 1
    assert out[S.COUNT].to_list() == [1010]
    assert out[S.CDR3_NT].to_list() == ["ACGTACGT"]           # parent survives
