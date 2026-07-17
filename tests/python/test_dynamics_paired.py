"""Tests for vdjtools.dynamics.paired — N_eff recovery, p-value calibration, 5-way classes."""
import numpy as np
import polars as pl
import pytest
from scipy.stats import kstest

from vdjtools.dynamics import estimate_neff
from vdjtools.dynamics import test_pair as run_pair   # aliased: pytest collects `test_*` names
from vdjtools.dynamics.paired import (CLASSES, DEFAULT_KEY, _bin_index, _downscale,
                                      _joined)
from vdjtools.io import schema as S


def _frame(counts: np.ndarray) -> pl.DataFrame:
    """A canonical clonotype frame from a per-clone count vector (zeros dropped)."""
    nz = np.nonzero(counts)[0]
    return pl.DataFrame({
        S.JUNCTION_AA: [f"CASS{i}F" for i in nz],
        S.V_CALL: ["TRBV7-9"] * nz.size,
        S.J_CALL: ["TRBJ1-1"] * nz.size,
        S.COUNT: counts[nz].astype(np.int64),
    })


def _two_step(f, n_s1, n_seq, rng):
    """The thesis's sampling model (Fig. 2.5): f -> N_S1 molecules -> N_seq reads."""
    mol = rng.multinomial(n_s1, f)
    return rng.multinomial(n_seq, mol / mol.sum())


def _pair(n_s1, n_seq, seed=0, n_clones=200_000, expand=None):
    """A replicate pair from one planted repertoire; `expand` boosts a slice of clones in b."""
    rng = np.random.default_rng(seed)
    f = rng.pareto(1.2, n_clones) + 1.0          # heavy-tailed, repertoire-like
    f /= f.sum()
    ca = _two_step(f, n_s1, n_seq, rng)
    fb = f.copy()
    if expand is not None:
        fb[:expand] *= 100.0
        fb /= fb.sum()
    cb = _two_step(fb, n_s1, n_seq, rng)
    truth = 1 / (1 / n_s1 + 1 / n_seq)
    return _frame(ca), _frame(cb), truth


# --- N_eff estimation (rung -1: recover a PLANTED truth) ---------------------------------

@pytest.mark.parametrize("n_s1,n_seq", [(50_000, 2_000_000), (200_000, 1_000_000),
                                        (1_000_000, 300_000)])
def test_estimate_neff_recovers_planted_truth(n_s1, n_seq):
    # The unambiguous check: simulate the thesis's own two-step process and assert the
    # estimator recovers 1/N_eff = 1/N_S1 + 1/N_seq. Everything downstream conflates the fit,
    # the pre-filter, the downscale and the Fisher test; this isolates the fit.
    # Tolerance is empirical, not aspirational: measured 0.95-1.07x across these regimes.
    a, b, truth = _pair(n_s1, n_seq)
    got = estimate_neff(a, b)
    assert 0.85 * truth < got < 1.15 * truth, f"N_eff {got:,.0f} vs planted {truth:,.0f}"


def _thesis_fit(a, b, min_count=2, bins=25, min_bin=10):
    """The thesis's own estimator (p.19): bin by f_B, mean+var of f_A, then x2."""
    j = _joined(a, b, list(DEFAULT_KEY))
    fa, fb = j["f_a"].to_numpy(), j["f_b"].to_numpy()
    sel = fb > min_count / int(j["count_b"].sum())
    fa, fb = fa[sel], fb[sel]
    idx = _bin_index(fb, bins)
    acc = [(np.log(fa[idx == g].var(ddof=1)) - np.log(fa[idx == g].mean()), int((idx == g).sum()))
           for g in range(bins)
           if (idx == g).sum() >= min_bin and fa[idx == g].var(ddof=1) > 0]
    return 2.0 * np.exp(-np.average([x for x, _ in acc], weights=[w for _, w in acc]))


def test_thesis_x2_overshoots_on_a_heavy_tail_which_is_why_we_fit_the_difference():
    # Pins the REASON dynamics.estimate_neff departs from the thesis's letter. The x2 corrects
    # for binning on the noisy f_B instead of the true f, and is exact only when the clone-size
    # distribution is flat -- then the bin contamination equals the sampling variance. Real
    # repertoires are heavy-tailed, the contamination is smaller, and a fixed x2 overshoots.
    # If this ever stops overshooting, the deviation is no longer justified: revisit.
    a, b, truth = _pair(200_000, 2_000_000)
    thesis = _thesis_fit(a, b)
    ours = estimate_neff(a, b)
    assert thesis / truth > 1.15, f"thesis form no longer overshoots ({thesis/truth:.2f}x)"
    assert abs(ours / truth - 1) < abs(thesis / truth - 1)      # ours is closer to the truth


def test_thesis_x2_is_exact_on_a_flat_clone_size_distribution():
    # The other side of the same mechanism, and the evidence that the x2 is not simply wrong:
    # with ONE true frequency there is no mis-binning to correct, so the x2 has nothing to undo
    # and overshoots by exactly 2. That is what makes it a distribution-dependent correction
    # rather than a constant -- and what makes a fixed factor unsafe on real data.
    rng = np.random.default_rng(3)
    f = np.full(200_000, 1 / 200_000)
    n_s1, n_seq = 200_000, 1_000_000
    a, b = _frame(_two_step(f, n_s1, n_seq, rng)), _frame(_two_step(f, n_s1, n_seq, rng))
    truth = 1 / (1 / n_s1 + 1 / n_seq)
    assert _thesis_fit(a, b) / truth == pytest.approx(2.0, rel=0.15)   # x2 overshoots by 2x
    assert estimate_neff(a, b) / truth == pytest.approx(1.0, rel=0.15)  # difference: unbiased


def test_estimate_neff_raises_rather_than_guessing():
    # A pair too shallow to fit must raise. A silent fallback (e.g. "use N_seq") would return
    # a number that is wrong by the entire oversampling factor, with no error.
    tiny = _frame(np.array([3, 2, 1]))
    with pytest.raises(ValueError, match="shallow|bins"):
        estimate_neff(tiny, tiny)


def test_estimate_neff_ref_n_is_free_not_hardcoded():
    # The thesis hard-codes 200,000 for the outlier pre-filter; that is a property of its
    # cohort, not a constant (it is 33x one of our datasets' whole libraries). Auto-derivation
    # must land near an explicitly-pinned value on data where the two agree.
    a, b, truth = _pair(200_000, 1_000_000)
    auto = estimate_neff(a, b, ref_n=None)
    pinned = estimate_neff(a, b, ref_n=truth)
    assert auto == pytest.approx(pinned, rel=0.25)


# --- the downscale ------------------------------------------------------------------------

def test_downscale_is_deterministic_not_a_resample():
    # preprocess.downsample draws a fresh hypergeometric sample; this must NOT. Drawing again
    # would stack a second noise process on the one N_eff already models.
    f = np.array([0.5, 0.25, 0.25])
    first = _downscale(f, 1000)
    assert np.array_equal(first, _downscale(f, 1000))     # same input -> same output, always
    assert first.tolist() == [500, 250, 250]


# --- p-value calibration (the acceptance test, thesis Fig. 2.13) ---------------------------

def test_pvalues_are_calibrated_on_a_replicate_pair_and_skewed_without_downscaling():
    # THE gate. On a replicate pair there are no true dynamics, so every "significant" clone
    # is a false positive and the p-values must not pile up near 0.
    #
    # NOT a two-sided KS against U(0,1): exact tests on small discrete counts are conservative
    # and atomic, so a two-sided KS rejects a *correct* engine outright. What must hold is the
    # one-sided claim -- no EXCESS of small p -- plus a directly-measured false-positive rate.
    #
    # The negative control is mandatory: the same assertions must FAIL when the downscale is
    # skipped on read counts. Without it the test would also pass on a no-op.
    a, b, _ = _pair(200_000, 2_000_000)           # 10x oversampled: N_eff << N_seq

    out = run_pair(a, b)
    p = out.filter(pl.col("dynamics") != "untested")["p_value"].to_numpy()
    assert p.size > 1000
    assert kstest(p, "uniform", alternative="greater").pvalue > 0.01   # no excess of small p
    for lvl in (0.01, 0.05):
        assert (p < lvl).mean() <= 1.5 * lvl, f"FPR {(p < lvl).mean():.4f} at alpha={lvl}"

    # Negative control: assume N_seq is the sample size (i.e. no downscale) and the same
    # checks must break -- that is the thesis's Fig. 2.13 top row.
    raw = run_pair(a, b, neff=None)
    praw = raw.filter(pl.col("dynamics") != "untested")["p_value"].to_numpy()
    assert (praw < 0.01).mean() > 3 * 0.01, "no-downscale control did not inflate the FPR"


def test_replicate_pair_calls_almost_nothing_significant():
    # The same claim in the units that matter: on a replicate pair, ~no clonotype should be
    # called changed at q<0.01.
    a, b, _ = _pair(200_000, 2_000_000)
    out = run_pair(a, b)
    changed = out.filter(pl.col("dynamics").is_in(["emergent", "expanded",
                                                   "contracted", "vanishing"]))
    assert changed.height / out.height < 0.01


# --- the 5-way classification -------------------------------------------------------------

def test_classes_partition_the_frame():
    a, b, _ = _pair(200_000, 1_000_000)
    out = run_pair(a, b)
    assert set(out["dynamics"].unique()) <= set(CLASSES)
    assert out.height == sum(out.filter(pl.col("dynamics") == c).height for c in CLASSES)


def test_untested_is_the_floor_not_a_verdict():
    # Clones below the testability floor must be `untested`, never `persistent`: "too rare to
    # test" and "tested, unchanged" are different claims and only one of them is evidence.
    a, b, _ = _pair(200_000, 1_000_000)
    out = run_pair(a, b, min_total=6)
    u = out.filter(pl.col("dynamics") == "untested")
    assert (u["count_a"] + u["count_b"] < 6).all()
    t = out.filter(pl.col("dynamics") != "untested")
    assert (t["count_a"] + t["count_b"] >= 6).all()


def test_emergent_and_vanishing_require_a_real_zero():
    a, b, _ = _pair(200_000, 1_000_000, expand=300)
    out = run_pair(a, b)
    assert (out.filter(pl.col("dynamics") == "emergent")["count_a"] == 0).all()
    assert (out.filter(pl.col("dynamics") == "vanishing")["count_b"] == 0).all()
    # expanded/contracted are the both-present cases
    both = out.filter(pl.col("dynamics").is_in(["expanded", "contracted"]))
    assert ((both["count_a"] > 0) & (both["count_b"] > 0)).all()


def test_planted_expansion_is_recovered_as_expanded_or_emergent():
    # Direction is b-relative-to-a: clones boosted 100x in b must come back as expanded (or
    # emergent, if they were absent from a), never contracted.
    a, b, _ = _pair(200_000, 1_000_000, expand=300)
    out = run_pair(a, b)
    boosted = out.filter(pl.col(S.JUNCTION_AA).str.extract(r"CASS(\d+)F")
                         .cast(pl.Int64) < 300)
    called = boosted.filter(pl.col("dynamics").is_in(["expanded", "emergent"]))
    assert called.height > 0.3 * boosted.height, "planted expansion not recovered"
    assert boosted.filter(pl.col("dynamics") == "contracted").height == 0
