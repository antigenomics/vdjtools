"""Tests for vdjtools.stats.diversity — hand-computed estimator values."""
import math

import numpy as np
import polars as pl

from vdjtools.io import schema as S
from vdjtools import stats


def _frame(counts):
    n = len(counts)
    df = pl.DataFrame({
        S.V_CALL: ["TRBV1"] * n, S.J_CALL: ["TRBJ1"] * n,
        S.CDR3_AA: [f"CASS{i}" for i in range(n)], S.COUNT: counts,
    })
    return S.add_locus(S.normalize(df, recompute_freq=True))


def test_chao1_bias_corrected():
    # counts [100,1,1,2,1]: Sobs=5, F1=3, F2=1 -> 5 + 3*2/(2*2) = 6.5
    assert stats.chao1(np.array([100, 1, 1, 2, 1])) == 6.5
    # all singletons [1,1,1]: F1=3,F2=0 -> 3 + 3*2/(2*1) = 6
    assert stats.chao1(np.array([1, 1, 1])) == 6.0
    # no singletons -> F0 = 0 -> Chao1 == Sobs
    assert stats.chao1(np.array([5, 5, 5])) == 3.0


def test_uniform_hill_numbers():
    c = np.array([1, 1, 1, 1])                # H = ln 4, uniform
    assert math.isclose(stats.shannon_wiener(c), 4.0, rel_tol=1e-12)
    assert math.isclose(stats.normalized_shannon_wiener(c), 1.0, rel_tol=1e-12)
    assert math.isclose(stats.inverse_simpson(c), 4.0, rel_tol=1e-12)


def test_inverse_simpson_value():
    c = np.array([100, 1, 1, 2, 1])
    # Sum f^2 = (100^2+1+1+4+1)/105^2 = 10007/11025
    assert math.isclose(stats.inverse_simpson(c), 11025 / 10007, rel_tol=1e-12)


def test_d50_dominance_fraction():
    # top clone alone covers 100/105 >= 0.5 -> 1 clone / 5 = 0.2
    assert stats.d50(np.array([100, 1, 1, 2, 1])) == 0.2
    # perfectly even 4 clones: need 2 of 4 to reach 0.5 -> 0.5
    assert stats.d50(np.array([1, 1, 1, 1])) == 0.5


def test_efron_thisted_at_least_observed():
    c = np.array([1, 1, 2, 3, 5])
    est = stats.efron_thisted(c)
    assert est >= stats.observed_richness(c)
    assert math.isfinite(est)


def test_chao_e_extrapolates_up_and_validates():
    c = np.array([1, 1, 2, 3, 5])
    n = int(c.sum())
    assert stats.chao_e(c, extrapolate_to=n) >= stats.observed_richness(c)
    # extrapolating farther yields >= closer target
    assert stats.chao_e(c, 2 * n) >= stats.chao_e(c, n)


def test_diversity_stats_frame_shape():
    df = _frame([100, 1, 1, 2, 1])
    out = stats.diversity_stats(df)
    assert out.height == 1
    assert out["reads"][0] == 105
    assert out["observed_diversity"][0] == 5
    assert out["chao1"][0] == 6.5
    for col in ("chaoE", "efron_thisted", "shannon_wiener",
                "normalized_shannon_wiener", "inverse_simpson", "d50"):
        assert col in out.columns
