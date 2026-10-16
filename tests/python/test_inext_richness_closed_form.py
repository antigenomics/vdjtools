"""Richness identity against exact finite-population combinatorics."""
import importlib
import math
from collections import Counter
from fractions import Fraction

import numpy as np
import pytest

from vdjtools import _core

ref = importlib.import_module("vdjtools.stats.inext")


@pytest.mark.parametrize("counts,depths", [
    ([1, 1, 2, 4], list(range(9))),
    ([8], [0, 1, 2, 4, 7, 8]),
    ([3, 3, 3], [0, 1, 2, 4, 8, 9]),
    ([1] * 15000, [0, 1, 2, 7500, 14999, 15000]),
    ([2] * 7500, [0, 1, 2, 7500, 14999, 15000]),
    ([40, 80, 120], [0, 1, 35, 100, 220, 239, 240]),
])
def test_native_richness_matches_exact_combinatorics(counts, depths):
    n = sum(counts)
    actual = _core.inext_curve(counts, [0], depths).qD[0]
    expected = []
    for m in depths:
        denominator = math.comb(n, m)
        expected.append(float(sum(
            frequency * Fraction(denominator - (math.comb(n - c, m) if n - c >= m else 0),
                                 denominator)
            for c, frequency in Counter(counts).items())))
    np.testing.assert_allclose(actual, expected, rtol=1e-9, atol=1e-10)


def test_native_richness_integer_endpoints_are_exact():
    counts = [1, 1, 2, 4]
    got = _core.inext_curve(counts, [0], [0, 1, 7, 8]).qD[0]
    assert got == [0.0, 1.0, 3.75, 4.0]


def test_native_richness_preserves_fractional_count_reference():
    counts = np.array([1.5, 2.5])
    sizes = [0., 1., 2., 3., 4.]
    actual = _core.inext_curve(counts.tolist(), [0], sizes).qD[0]
    expected = [ref._diversity_at(counts, m, 0) for m in sizes]
    np.testing.assert_allclose(actual, expected, rtol=1e-9, atol=1e-10)
