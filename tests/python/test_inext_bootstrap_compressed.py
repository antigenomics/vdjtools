"""Software oracles for the compressed iNEXT augmented-assemblage sampler."""
import importlib

import numpy as np
import pytest

from vdjtools import _core, stats

ref = importlib.import_module("vdjtools.stats.inext")


def test_small_bootstrap_keeps_seeded_expanded_draws():
    x = np.array([1, 1, 1, 2, 3, 7], dtype=float)
    p = ref._bootstrap_probs(x)
    expected = np.random.default_rng(7).multinomial(int(x.sum()), p / p.sum(), size=8)
    actual = list(stats.bootstrap_abundances(x, nboot=8, seed=7))
    assert len(actual) == 8
    assert all(np.array_equal(a, b[b > 0]) for a, b in zip(actual, expected))


def test_compressed_bootstrap_matches_exact_occupancy_moments(monkeypatch):
    x = np.array([1, 1, 1, 2, 3, 7], dtype=float)
    p = ref._bootstrap_probs(x)
    p /= p.sum()
    n = int(x.sum())
    # Force compression on a small population whose complete law is tractable.
    monkeypatch.setattr(ref, "_BOOTSTRAP_MAX_UNSEEN", 0)
    draws = stats.bootstrap_abundances(x, nboot=10000, seed=3)
    moments = np.array([(len(draw), np.count_nonzero(draw == 1)) for draw in draws])
    expected = [np.sum(1 - (1 - p) ** n), np.sum(n * p * (1 - p) ** (n - 1))]
    np.testing.assert_allclose(moments.mean(axis=0), expected, atol=0.04, rtol=0)


def test_large_singletons_keep_chao1_and_bound_population():
    x = np.ones(15000)
    p, k = ref._bootstrap_population(x)
    assert k == 112485001
    assert p.size == x.size + 1
    assert ref._asymptotic(x, 0) == 15000 + (14999 ** 2) / 2
    a = list(stats.bootstrap_abundances(x, nboot=3, seed=9))
    b = list(stats.bootstrap_abundances(x, nboot=3, seed=9))
    assert len(a) == len(b) == 3
    for first, second in zip(a, b):
        assert np.array_equal(first, second)
        assert first.dtype == np.int64
        assert np.all(first > 0) and first.sum() == 15000
        assert first.size <= 15000
    se = np.asarray(_core.inext_bootstrap(x.tolist(), [0], [7500., 15000., 30000.], 3, 9))
    assert se.shape == (1, 3)
    assert np.all(np.isfinite(se)) and np.all(se >= 0)


def test_bootstrap_count_validation():
    for x in ([1.5, 2], [1, float("inf")]):
        with pytest.raises(ValueError, match="finite integers"):
            list(stats.bootstrap_abundances(x))
    with pytest.raises(ValueError, match="nboot"):
        list(stats.bootstrap_abundances([1, 2], nboot=-1))


def test_native_compressed_bootstrap_matches_python():
    x = np.ones(100)  # ceil(f0) = 4901: compression in both implementations
    native = np.asarray(_core.inext_bootstrap(x.tolist(), [0], [100.], 4000, 7))
    python = ref._bootstrap_se(x, [100.], [0], 4000, 7)
    np.testing.assert_allclose(native, python, rtol=0.10)
