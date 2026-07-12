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
        S.JUNCTION_AA: [f"CASS{i}" for i in range(n)], S.COUNT: counts,
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
    # legacy getDxxIndex = 1 - k/Sobs. Top clone alone covers 100/105 >= 0.5,
    # so k=1 of Sobs=5 -> 1 - 1/5 = 0.8
    assert stats.d50(np.array([100, 1, 1, 2, 1])) == 0.8
    # perfectly even 4 clones: need k=2 of 4 to reach 0.5 -> 1 - 2/4 = 0.5
    assert stats.d50(np.array([1, 1, 1, 1])) == 0.5


def _efron_reference(counts, max_depth=20, cv_threshold=0.05):
    """Independent re-derivation of the Efron-Thisted Euler/CV series (for pinning)."""
    from collections import Counter
    sobs = len(counts)
    fx = Counter(int(c) for c in counts)
    s = float(sobs)
    for depth in range(1, max_depth + 1):
        h = [0.0] * depth
        for y in range(1, depth + 1):
            for x in range(1, y + 1):
                coef = math.comb(y - 1, x - 1)
                h[x - 1] += coef if x % 2 == 1 else -coef
        nx = [fx.get(y, 0) for y in range(1, depth + 1)]
        s = sobs + sum(h[i] * nx[i] for i in range(depth))
        d = math.sqrt(sum(h[i] * h[i] * nx[i] for i in range(depth)))
        if s != 0 and d / s >= cv_threshold:
            break
    return float(s)


def test_efron_thisted_all_singletons_closed_form():
    # With the (non-legacy) max_count cap removed, only legacy's CV stopping rule
    # remains. For an all-singletons sample the CV = sqrt(Sobs)/(2*Sobs) fires at
    # depth 1 (= 0.158 >= 0.05 for Sobs=10), so the estimate is Sobs + Sobs = 2*Sobs.
    # (This is the true legacy-faithful value; an unconditional run to max_depth
    #  would instead give Sobs*(1+max_depth)=210, but the CV rule prevents that.)
    assert stats.efron_thisted(np.array([1] * 10)) == 20.0


def test_efron_thisted_matches_independent_reference():
    c = np.array([5, 4, 3, 3, 2, 2, 1])
    assert stats.efron_thisted(c) == _efron_reference(c) == 8.0


def test_chao_e_extrapolation_hand_value():
    c = np.array([1, 1, 2, 3, 5])                     # n=12, Sobs=5, F1=2, F2=1, F0=0.5
    n = int(c.sum())
    # Chao extrapolation: 5 + 0.5*(1 - (1 - 2/(12*0.5))^12) = 5 + 0.5*(1 - (2/3)^12)
    assert math.isclose(stats.chao_e(c, 24), 5.496146326685371, rel_tol=1e-9)
    # extrapolate_to == n gives m*=0 -> exactly Sobs (the old test's "tautology")
    assert stats.chao_e(c, n) == float(stats.observed_richness(c))
    assert stats.chao_e(c, 2 * n) >= stats.chao_e(c, n)


def test_diversity_degenerate_empty_and_single():
    empty = np.array([], dtype=np.int64)
    assert stats.observed_richness(empty) == 0
    assert stats.chao1(empty) == 0.0
    assert stats.chao_e(empty) == 0.0
    assert stats.efron_thisted(empty) == 0.0
    assert stats.inverse_simpson(empty) == 0.0
    assert stats.d50(empty) == 0.0
    # single clonotype: Sobs fall-throughs and freq == 1.0
    one = np.array([5], dtype=np.int64)
    assert stats.observed_richness(one) == 1
    assert stats.chao1(one) == 1.0                    # F0=0 -> Chao1 == Sobs
    assert stats.inverse_simpson(one) == 1.0
    assert stats.d50(one) == 0.0                      # 1 - 1/1
    df = _frame([5])
    assert df[S.FREQ].to_list() == [1.0]              # single clonotype -> freq 1.0


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
