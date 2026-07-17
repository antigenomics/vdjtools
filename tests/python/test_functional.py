"""Leinster-Cobbold similarity-aware (functional) diversity.

The load-bearing test is ``Z = I``: L&C's ᑫD^Z reduces to the plain Hill numbers exactly when
nothing resembles anything but itself. That is the identity the whole construction is defined
against, so it doubles as a free oracle — richness / exp(Shannon) / inverse-Simpson are already
implemented independently in ``stats.diversity`` and ``stats.inext``.
"""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from vdjtools.stats import functional_diversity, inverse_simpson, observed_richness
from vdjtools.stats.functional import _profile


def _sample(n=60, seed=0) -> pl.DataFrame:
    """n distinct clonotypes with a skewed count profile."""
    rng = np.random.default_rng(seed)
    aa = set()
    while len(aa) < n:
        aa.add("CAS" + "".join(rng.choice(list("ACDEFGHIKLMNPQRSTVWY"), 8)) + "F")
    aa = sorted(aa)
    counts = (rng.pareto(1.2, n) * 10 + 1).astype(np.int64)
    return pl.DataFrame({
        "junction_aa": aa,
        "junction_nt": ["ACG"] * n,
        "v_call": ["TRBV19*01"] * n,
        "j_call": ["TRBJ2-7*01"] * n,
        "duplicate_count": counts,
        "frequency": counts / counts.sum(),
    })


def test_identity_kernel_recovers_hill_numbers():
    """Z = I  =>  ᑫD^Z == the plain Hill numbers. The defining special case."""
    df = _sample()
    got = functional_diversity(df, q=(0, 1, 2), kernel="identity")
    d = dict(zip(got["q"].to_list(), got["diversity"].to_list()))

    p = (df["duplicate_count"].to_numpy() / df["duplicate_count"].sum())
    assert d[0.0] == pytest.approx(observed_richness(df["duplicate_count"].to_numpy()))
    assert d[1.0] == pytest.approx(float(np.exp(-np.sum(p * np.log(p)))))       # exp(Shannon)
    assert d[2.0] == pytest.approx(inverse_simpson(df["duplicate_count"].to_numpy()))


def test_identity_kernel_rao_is_one_minus_simpson():
    """Z = I  =>  Rao's Q = 1 - Σp² — the Gini-Simpson index."""
    df = _sample()
    got = functional_diversity(df, q=(1,), kernel="identity")
    p = df["duplicate_count"].to_numpy() / df["duplicate_count"].sum()
    assert got["rao"][0] == pytest.approx(1.0 - float(np.sum(p**2)))


def test_similarity_never_increases_diversity():
    """Folding in real similarity can only merge clonotypes, never split them.

    Z >= I entrywise => (Zp)_i >= p_i => every ᑫD^Z <= the plain Hill number. This is the
    property that makes the measure worth having: it is what "functional" buys over richness.
    """
    pytest.importorskip("seqtree", reason="seqtree is the similarity kernel engine")
    df = _sample()
    plain = functional_diversity(df, q=(0, 1, 2), kernel="identity")["diversity"].to_list()
    fuzzy = functional_diversity(df, q=(0, 1, 2), kernel="exp")["diversity"].to_list()
    for q, (a, b) in enumerate(zip(fuzzy, plain)):
        assert a <= b + 1e-9, f"q={q}: similarity-aware {a} exceeded plain Hill {b}"
    assert fuzzy[0] < plain[0]           # some CDR3s really are neighbours; not a no-op kernel


def test_profile_matches_the_closed_form():
    """_profile against a hand-evaluated ᑫD^Z on a 3-clonotype toy, incl. the q→1 limit."""
    p = np.array([0.5, 0.3, 0.2])
    Z = np.array([[1.0, 0.5, 0.0], [0.5, 1.0, 0.0], [0.0, 0.0, 1.0]])
    zp = Z @ p
    assert _profile(p, zp, 0) == pytest.approx(float(np.sum(p / zp)))
    assert _profile(p, zp, 1) == pytest.approx(float(np.exp(-np.sum(p * np.log(zp)))))
    assert _profile(p, zp, 2) == pytest.approx(1.0 / float(np.sum(p * zp)))
    # q→1 is a genuine limit, not a separate formula: approach it from both sides.
    lo, hi = _profile(p, zp, 0.999), _profile(p, zp, 1.001)
    assert lo == pytest.approx(_profile(p, zp, 1), rel=1e-3)
    assert hi == pytest.approx(_profile(p, zp, 1), rel=1e-3)


def test_zero_abundance_clonotypes_are_dropped():
    """A p_i = 0 clonotype must not contribute — and must not put 0*log 0 in the sum."""
    p = np.array([0.6, 0.0, 0.4])
    zp = np.array([0.6, 0.4, 0.4])
    assert np.isfinite(_profile(p, zp, 0))
    assert np.isfinite(_profile(p, zp, 1))
    assert _profile(np.zeros(3), zp, 1) == 0.0


def test_presence_weight_is_the_unweighted_profile():
    """weight='presence' ignores abundance: Z=I then gives richness at every q."""
    df = _sample()
    got = functional_diversity(df, q=(0, 1, 2), kernel="identity", weight="presence")
    for d in got["diversity"].to_list():
        assert d == pytest.approx(df.height)
