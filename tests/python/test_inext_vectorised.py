"""``_rtd_moment`` sums over the abundance spectrum in one flat pass; it must equal the loop.

The loop it replaced was ~200 round trips through numpy per call on a real repertoire -- that is
the *distinct count* cardinality, not the clonotype count -- and it is called 28 times per sample
(orders q=0 and q=1, floor and ceil of a non-integer m, seven loci). Measured 2026-09-26,
``div_block`` over seven loci: 58.3 -> 3.4 ms on a 415-clonotype sample and 286.9 -> 70.8 ms on a
4,863-clonotype one.

A speedup that moves a diversity value is worthless here -- ``vsig:div:*`` is compared against a
frozen scale reference -- so these pin the equality, at the tolerance floating-point reassociation
actually allows, and not the time.
"""
from __future__ import annotations

import importlib

import numpy as np
import pytest

I = importlib.import_module("vdjtools.stats.inext")


def _loop(cs, fs, n, m, gfun):
    """The pre-3.16.0 body, verbatim."""
    total = 0.0
    log_cnm = I._logbinom(n, m)
    for cval, fcount in zip(cs, fs):
        kmax = int(min(cval, m))
        k = np.arange(1, kmax + 1)
        logterm = I._logbinom(cval, k) + I._logbinom(n - cval, m - k) - log_cnm
        total += fcount * float(np.sum(gfun(k) * np.exp(logterm)))
    return total


def _g(q, m):
    if q == 0:
        return lambda k: np.ones_like(k, dtype=np.float64)
    return lambda k: -(k / m) * np.log(k / m)


@pytest.mark.parametrize("species,a", [(60, 1.3), (200, 1.5), (2000, 1.8), (500, 2.5)])
@pytest.mark.parametrize("q", [0, 1])
def test_the_flat_pass_equals_the_loop(species, a, q):
    x = np.ceil(np.random.default_rng(species).zipf(a, species).clip(1, 5000)).astype(np.int64)
    n = int(x.sum())
    cs, fs = I._spectrum(x)
    for m in (max(1, n // 8), max(1, n // 2), n - 1):
        g = _g(q, m)
        want, got = _loop(cs, fs, n, m, g), I._rtd_moment(cs, fs, n, m, g)
        assert got == pytest.approx(want, rel=1e-11, abs=1e-12), (species, q, m)


def test_chunking_does_not_change_the_answer():
    """The flat buffer is bounded; the bound must be invisible in the result."""
    x = np.ceil(np.random.default_rng(7).zipf(1.4, 400).clip(1, 40000)).astype(np.int64)
    n = int(x.sum())
    cs, fs = I._spectrum(x)
    m, g = n // 3, _g(1, n // 3)
    want = I._rtd_moment(cs, fs, n, m, g)
    original = I._FLAT_TERMS
    try:
        for budget in (1, 7, 64, 1 << 14):          # force many chunks, including one per group
            I._FLAT_TERMS = budget
            assert I._rtd_moment(cs, fs, n, m, g) == pytest.approx(want, rel=1e-12)
    finally:
        I._FLAT_TERMS = original


def test_chunks_cover_every_group_exactly_once():
    """A dropped or repeated group would silently under- or over-count the moment."""
    original = I._FLAT_TERMS
    try:
        for budget in (1, 3, 1000):
            I._FLAT_TERMS = budget
            for kmax in ([1], [5, 5, 5], [1, 1, 900, 1], list(range(1, 40))):
                arr = np.array(kmax, dtype=np.int64)
                got = list(I._flat_chunks(arr))
                assert got[0][0] == 0 and got[-1][1] == arr.size
                assert all(b == c for (_, b), (c, _) in zip(got, got[1:]))
                assert all(a < b for a, b in got), "an empty chunk would loop forever"
    finally:
        I._FLAT_TERMS = original


def test_a_spectrum_with_nothing_to_sum_is_zero():
    """``m`` below every count leaves no k range; the loop returned 0.0 by falling through."""
    cs, fs = np.array([3.0, 4.0]), np.array([2.0, 1.0])
    assert I._rtd_moment(cs, fs, 10.0, 0, _g(0, 1)) == 0.0


def test_diversity_estimates_are_unchanged_end_to_end():
    """The number the signature actually ships, not the private helper."""
    import polars as pl

    from vdjtools.signature import blocks as B

    r = np.random.default_rng(3)
    df = pl.DataFrame({
        "v_call": ["TRBV20-1"] * 300, "j_call": ["TRBJ2-2"] * 300, "c_call": [None] * 300,
        "junction_aa": ["C" + "".join(r.choice(list("ACDEFGHIKLMNPQRSTVWY"), 12)) + "F"
                        for _ in range(300)],
        "duplicate_count": np.ceil(r.zipf(1.5, 300).clip(1, 900)).astype(int).tolist(),
    }, schema_overrides={"c_call": pl.Utf8})
    clean, _ = B.sanitise(df)

    fast = B.div_block(clean, 0.20, tier_full=True)
    vec = I._rtd_moment
    try:
        I._rtd_moment = _loop
        slow = B.div_block(clean, 0.20, tier_full=True)
    finally:
        I._rtd_moment = vec

    for k in fast:
        if np.isfinite(fast[k]) or np.isfinite(slow[k]):
            assert fast[k] == pytest.approx(slow[k], rel=1e-10, nan_ok=True), k
