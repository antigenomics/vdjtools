"""``columns=`` must skip the work, not merely drop the result.

The Pgen block is ~96% of ``vsig``'s cost on a seven-locus sample and the number is irreducible
(see :func:`vdjtools.signature.blocks.pgen_block` -- it is the D trim state space, 9,212 states on
IGH against 297 on TRB). Declining the block is therefore the only large saving available, and a
``columns=`` that computed it and then dropped it would look identical in the output while saving
nothing. These tests assert on *execution*, not on elapsed time: a stopwatch cannot tell a skipped
block from a fast one, and would flake on a loaded box.
"""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from vdjtools.signature import assemble, layout as L
from vdjtools.signature import blocks as B


def _frame(n, v, j, c=None, seed=0):
    r = np.random.default_rng(seed)
    aa = list("ACDEFGHIKLMNPQRSTVWY")
    return pl.DataFrame({
        "v_call": [v] * n, "j_call": [j] * n, "c_call": [c] * n,
        "junction_aa": ["C" + "".join(r.choice(aa, 12)) + "F" for _ in range(n)],
        "duplicate_count": np.ceil(r.zipf(1.5, n).clip(1, 500)).astype(int).tolist(),
    }, schema_overrides={"c_call": pl.Utf8})


@pytest.fixture
def sample():
    return {"TRB": _frame(300, "TRBV20-1", "TRBJ2-2", seed=1),
            "IGH": _frame(200, "IGHV1-2", "IGHJ4", "IGHM", seed=2)}


NO_PGEN = [c for c in L.columns("standard", "vsig") if not c.startswith("vsig:pgen:")]


def _identity(x):
    """Module-level so ``functools.partial`` over it survives pickling under ``spawn``."""
    return x


def test_declining_pgen_does_not_call_pgen_block(sample, monkeypatch):
    """The point of the parameter. A booby-trapped block proves it never ran."""
    def boom(*a, **k):
        raise AssertionError("pgen_block ran despite being excluded by columns=")

    monkeypatch.setattr(B, "pgen_block", boom)
    out = assemble.vsig(sample, tier="standard", columns=NO_PGEN)
    assert list(out) == NO_PGEN


def test_asking_for_pgen_still_calls_it(sample, monkeypatch):
    """The other half of the guard: without the exclusion the block must still run."""
    monkeypatch.setattr(B, "pgen_block", lambda *a, **k: (_ for _ in ()).throw(RuntimeError("ran")))
    with pytest.raises(RuntimeError, match="ran"):
        assemble.vsig(sample, tier="standard")


def test_skipping_a_block_moves_no_other_column(sample):
    """Blocks must be independent: dropping one cannot perturb the rest."""
    full = assemble.vsig(sample, tier="standard")
    cut = assemble.vsig(sample, tier="standard", columns=NO_PGEN)
    moved = [k for k in cut
             if not (cut[k] == full[k] or (np.isnan(cut[k]) and np.isnan(full[k])))]
    assert moved == []


def test_every_optional_block_can_be_declined(sample, monkeypatch):
    """Not just Pgen. Each gated block is booby-trapped in turn."""
    for name, chan in (("qc_block", "qc"), ("depth_block", "depth"), ("clon_block", "clon"),
                       ("len_block", "len"), ("pgen_block", "pgen"), ("aa_block", "aa"),
                       ("pchem_block", "pchem"), ("iso_block", "iso"), ("shm_block", "shm")):
        keep = [c for c in L.columns("full", "vsig") if f":{chan}:" not in c]
        with monkeypatch.context() as m:
            m.setattr(B, name, lambda *a, _n=name, **k: (_ for _ in ()).throw(
                AssertionError(f"{_n} ran despite being excluded")))
            out = assemble.vsig(sample, tier="full", columns=keep)
        assert list(out) == keep, chan


def test_mask_estimable_survives_dropping_the_div_columns(sample):
    """``mask:*:estimable`` is derived from the div block, so div stays even when div is not kept."""
    keep = [c for c in L.columns("standard", "vsig") if ":div:" not in c]
    out = assemble.vsig(sample, tier="standard", columns=keep)
    assert out["vsig:mask:TRB:estimable"] == 1.0
    assert out["vsig:mask:TRA:estimable"] == 0.0


def test_columns_are_returned_in_layout_order_not_caller_order(sample):
    """The contract is positional: a caller's ordering must not become the frame's ordering."""
    shuffled = list(reversed(NO_PGEN))
    out = assemble.vsig(sample, tier="standard", columns=shuffled)
    assert list(out) == NO_PGEN


def test_columns_outside_the_tier_are_ignored(sample):
    """An intersection, not a union: a ``full``-tier name cannot widen a ``core`` request."""
    out = assemble.vsig(sample, tier="core",
                        columns=[*L.columns("core", "vsig"), "vsig:pchem:TRB:cdr3_kf1"])
    assert list(out) == L.columns("core", "vsig")


def test_cohort_columns_skips_the_block_too(sample, monkeypatch):
    """``vsig_cohort`` must push the subset into the per-sample call, not filter afterwards."""
    monkeypatch.setattr(B, "pgen_block", lambda *a, **k: (_ for _ in ()).throw(
        AssertionError("pgen_block ran in the cohort path")))
    out = assemble.vsig_cohort({"s1": sample, "s2": sample}, tier="standard", columns=NO_PGEN)
    assert out.height == 2
    assert out.columns == ["sample_id", *NO_PGEN]


def test_pgen_n_max_reaches_the_block(sample, monkeypatch):
    """The knob is only useful if it arrives; it used to stop at ``vsig``'s signature."""
    seen = {}

    def spy(df, locus, *, q05=None, n_max=2000, threads=0):
        seen[locus] = n_max
        return {"mean_log10": np.nan, "sd_log10": np.nan, "frac_atypical": np.nan}

    monkeypatch.setattr(B, "pgen_block", spy)
    assemble.vsig(sample, tier="standard", pgen_n_max=137)
    assert seen == {"TRB": 137, "IGH": 137}


def test_columns_and_deferred_samples_survive_the_process_pool(sample):
    """``columns=`` and the deferred callable both have to cross a spawn boundary.

    The in-process path exercises neither pickling nor the worker's own resolution, so a subset
    that works at ``n_jobs=1`` can still die at ``n_jobs=2``.
    """
    import functools

    # functools.partial over a module-level function, not a lambda: under `spawn` the deferred
    # callable is pickled, and a lambda or a closure cannot be. That is the shape real callers
    # use (`functools.partial(io.read, path)`), and it is the only shape that works.
    items = [(f"s{i}", functools.partial(_identity, sample)) for i in range(4)]
    out = assemble.vsig_cohort(items, tier="core", n_jobs=2)
    assert out.height == 4 and out.columns == ["sample_id", *L.columns("core", "vsig")]

    core_no_div = [c for c in L.columns("core", "vsig") if ":div:" not in c]
    sub = assemble.vsig_cohort(items, tier="core", n_jobs=2, columns=core_no_div)
    assert sub.columns == ["sample_id", *core_no_div]
    for c in sub.columns[1:]:
        assert sub[c].equals(out[c], null_equal=True), c


def test_a_deferred_sample_is_read_inside_the_worker(sample):
    """The documented ``O(n_jobs)`` memory path: a zero-argument callable per sample.

    It was documented on ``vsig_cohort`` and ``signature_cohort`` but implemented only in
    mirpy's ``_one_rsig``, so passing one died on ``'function' object has no attribute
    'columns'`` -- the low-memory path was unusable rather than merely slow.
    """
    calls = []

    def deferred():
        calls.append(1)
        return sample

    out = assemble.vsig_cohort({"s1": deferred}, tier="core")
    assert out.height == 1 and calls == [1]
    assert out.columns == ["sample_id", *L.columns("core", "vsig")]
