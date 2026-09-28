"""``infer_nt_batch``: the per-row loop's answer, faster -- never a different answer.

Filed as antigenomics/vdjtools#181. ``infer_nt`` wrapped a native DP in per-row Python, and on a
real annotation table that wrapper was the cost: 86% of vdjdb-db's build at ~200,000 distinct
junctions. Batching is therefore a speed change with a hard correctness contract -- every field of
every row must equal what the per-row loop returned -- which is what these tests pin.

The multi-allele case is not decoration. It is the only mode where a row pools more candidates than
``n_best``, and reporting the post-truncation count there was a real bug this file caught.
"""
from __future__ import annotations

import inspect

import pytest

from vdjtools.model import infer_nt, infer_nt_batch, load_bundled

#: Real human TRB junctions with gene-level V/J, as an aligner reports them.
CASES = [("CASSQDLNTEAFF", "TRBV4-3*01", "TRBJ1-1*01"),
         ("CASSLGQGAYEQYF", "TRBV5-1*01", "TRBJ2-7*01"),
         ("CASSPRTGELFF", "TRBV7-9*01", "TRBJ2-2*01")]

#: All four call modes infer_nt accepts, because they take different paths through the flatten.
MODES = [("CASSQDLNTEAFF", "TRBV4-3*01", "TRBJ1-1*01"),               # pinned
         ("CASSLGQGAYEQYF", "TRBV5-1*01,TRBV6-5*01", "TRBJ2-7*01"),   # multi-allele V
         ("CASSPRTGELFF", None, "TRBJ2-2*01"),                        # V marginalized
         ("CASSPRTGELFF", None, None)]                                # both marginalized

FIELDS = ("cdr3_nt", "v_call", "j_call", "v_end", "j_start", "d_call", "d_start", "d_end",
          "pgen", "scenario_p", "n_candidates", "runner_up_pgen")


@pytest.fixture(scope="module")
def trb():
    return load_bundled("TRB")


def test_it_is_reachable_from_vdjtools_model():
    import vdjtools.model as M

    assert "infer_nt_batch" in M.__all__ and hasattr(M, "infer_nt_batch")


def test_the_batch_is_the_per_row_loop_exactly(trb):
    """Field for field, over every call mode. Batching may not move a single returned number."""
    cases = MODES * 25
    got = infer_nt_batch(trb, [c[0] for c in cases], v=[c[1] for c in cases],
                         j=[c[2] for c in cases])
    assert got.height == len(cases)
    assert list(got.columns) == list(FIELDS)
    for i in (0, 1, 2, 3, 61, 98):
        aa, v, j = cases[i]
        ref = infer_nt(trb, aa, v=v, j=j)
        row = got.row(i, named=True)
        assert ref is not None, (aa, v, j)
        for f in FIELDS:
            assert row[f] == pytest.approx(getattr(ref, f)) if isinstance(
                getattr(ref, f), float) else row[f] == getattr(ref, f), f


def test_n_candidates_counts_every_distinct_candidate_not_just_the_rescored_ones(trb):
    """A multi-allele call pools more candidates than ``n_best``; the count is pre-truncation.

    Reporting the re-scored count instead is silent -- every other field stays right -- and it is
    exactly what this batch got wrong first time round.
    """
    aa, v, j = "CASSGRGLNTEAFF", "TRBV19*01,TRBV6-5*01", "TRBJ1-1*01"
    ref = infer_nt(trb, aa, v=v, j=j, n_best=8)
    got = infer_nt_batch(trb, [aa], v=[v], j=[j], n_best=8)
    assert ref.n_candidates > 8, "case no longer pools more than n_best; pick another"
    assert got.row(0, named=True)["n_candidates"] == ref.n_candidates


def test_the_result_does_not_depend_on_the_thread_count(trb):
    """The work is disjoint per row, so this must be identical -- not close -- at any width."""
    cases = MODES * 25
    aa = [c[0] for c in cases]
    v = [c[1] for c in cases]
    j = [c[2] for c in cases]
    one = infer_nt_batch(trb, aa, v=v, j=j, threads=1)
    many = infer_nt_batch(trb, aa, v=v, j=j, threads=8)
    auto = infer_nt_batch(trb, aa, v=v, j=j, threads=0)
    assert one.equals(many) and one.equals(auto)


def test_a_declined_row_is_a_null_row_at_its_own_index(trb):
    """``Scenario | None`` batches to a null ROW, never an absent one -- the caller joins by index.

    Deliberately unlike ``best_aa_scenarios_batch``, which signals a decline by absence because it
    returns k rows per query. Here one input row owes one output row.
    """
    aa = ["CASSQDLNTEAFF", "CASSXQDLNTEAFF", "", "CASSPRTGELFF"]
    got = infer_nt_batch(trb, aa, v=["TRBV4-3*01", "TRBV4-3*01", None, "TRBV7-9*01"],
                         j=["TRBJ1-1*01", "TRBJ1-1*01", None, "TRBJ2-2*01"])
    assert got.height == 4
    assert got["cdr3_nt"][1] is None and got["cdr3_nt"][2] is None
    assert got["pgen"][1] is None and got["n_candidates"][2] is None
    for i in (0, 3):                                   # the neighbours are untouched
        ref = infer_nt(trb, aa[i], v=["TRBV4-3*01", None, None, "TRBV7-9*01"][i],
                       j=["TRBJ1-1*01", None, None, "TRBJ2-2*01"][i])
        assert got["cdr3_nt"][i] == ref.cdr3_nt


def test_a_length_mismatch_is_rejected_rather_than_zipped_short(trb):
    with pytest.raises(ValueError, match="same length"):
        infer_nt_batch(trb, ["CASSPRTGELFF", "CASSPRTGELFF"], v=["TRBV7-9*01"])
    with pytest.raises(ValueError, match="same length"):
        infer_nt_batch(trb, ["CASSPRTGELFF", "CASSPRTGELFF"], j=["TRBJ2-2*01"])


def test_infer_nt_keeps_its_signature():
    """The single-row entry point is unchanged public API; it now delegates to the batch core.

    Pinned as an interface fact rather than trusted to review -- a parameter that survives in a
    signature while being ignored has shipped from this repo before.
    """
    p = inspect.signature(infer_nt).parameters
    assert list(p) == ["model_or_prep", "cdr3_aa", "v", "j", "n_best", "keep", "top_vj"]
    assert p["n_best"].default == 8 and p["keep"].default == 1e-2 and p["top_vj"].default == 8


def test_the_batch_does_not_carry_the_reference_paths_dead_knobs():
    """``keep``/``top_vj`` are reference-search-only; the native path ignores them."""
    p = inspect.signature(infer_nt_batch).parameters
    assert "keep" not in p and "top_vj" not in p
    assert p["threads"].default == 0 and p["n_best"].default == 8


def test_more_threads_buy_wall_time_but_cannot_buy_all_of_it(trb):
    """The GIL is released for the native stages, so threads must help -- measurably, not notionally.

    The bar is deliberately **not** "quadrupling the workers quarters the wall time", which is the
    right test for a fully native batch (``test_thread_scaling.py``) and the wrong one here. The
    codon reconstruction is pure Python and holds the GIL: measured at 35% of the per-row cost on
    human TRB and 87% on TRA, which caps a 4-thread run at ~1.9x by Amdahl and makes the TRA batch
    essentially flat in ``threads``. Measured 1.86x on 4,000 real VDJdb TRB keys, 16-core M-series.

    So this asserts the guard is there, not that it scales freely. Do not "fix" this test upward
    without first moving the reconstruction off the GIL.
    """
    import time

    aa = [c[0] for c in CASES] * 120
    v = [c[1] for c in CASES] * 120
    j = [c[2] for c in CASES] * 120
    infer_nt_batch(trb, aa[:8], v=v[:8], j=j[:8])           # warm the pack/prep caches

    t = time.perf_counter()
    one = infer_nt_batch(trb, aa, v=v, j=j, threads=1)
    t1 = time.perf_counter() - t
    t = time.perf_counter()
    many = infer_nt_batch(trb, aa, v=v, j=j, threads=4)
    t4 = time.perf_counter() - t

    assert one.equals(many), "threads changed the answer"
    assert t1 / t4 > 1.25, f"4 threads bought only {t1 / t4:.2f}x -- is the GIL still released?"
