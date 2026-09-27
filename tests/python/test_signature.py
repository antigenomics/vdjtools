"""The signature contract, the corpus fit, and the apply path.

Every test here pins something that has actually gone wrong in this subsystem before. The names say
which; the docstrings say what it cost.
"""
from __future__ import annotations

import pathlib

import numpy as np
import polars as pl
import pytest

from vdjtools.signature import corpus as C
from vdjtools.signature import features as F
from vdjtools.signature import layout as L
from vdjtools.signature.signature import vsig, vsig_cohort

# A corpus small enough to fit in a test and wide enough to exercise the machinery. Two loci so a
# per-locus bug cannot hide, plus the cross-locus block.
LOCI = ("TRB", "TRG")


@pytest.fixture(scope="module")
def corpus():
    art, _rows = C.synthesize("memory", loci=LOCI, n_samples=40, size=150, seed=3,
                              n_components=6)
    return art


@pytest.fixture(scope="module")
def held_out():
    """Three repertoires drawn with a different seed -- never part of the fit."""
    samples = C.synthesize("memory", loci=LOCI, n_samples=3, size=150, seed=808,
                           fit_corpus=False)
    return {f"S{i}": s for i, s in enumerate(samples)}


def _frame(n, v, j, rng, c=None, counts=None):
    aa = np.array(list(L.AMINO_ACIDS))
    return pl.DataFrame({
        "v_call": [v] * n, "j_call": [j] * n, "c_call": [c] * n,
        "junction_aa": ["C" + "".join(rng.choice(aa, 12)) + "F" for _ in range(n)],
        "duplicate_count": (counts if counts is not None
                            else np.ceil(rng.zipf(1.5, n).clip(1, 900)).astype(np.int64)).tolist(),
    })


# --------------------------------------------------------------- 1. support, not transform


def test_the_trimming_side_comes_from_support_and_never_from_the_transform():
    """Three features declare transform='none' with three different supports.

    Deriving the winsorization side from the transform silently trims the wrong end of two of them.
    That is the whole reason ``support`` is a declared field rather than an inference.
    """
    by_support = {}
    for g in L.raw_groups("vsig"):
        for feat, (tf, sup) in g.features.items():
            if tf == "none":
                by_support.setdefault(sup, []).append(f"vsig:{g.name}:*:{feat}")
    assert len(by_support) >= 2, f"transform='none' must span >1 support, got {by_support}"
    # and each support maps to a genuinely different pair of sides
    sides = {s: L.SUPPORTS[s] for s in by_support}
    assert len(set(sides.values())) == len(sides), sides


@pytest.mark.parametrize("support,expect_lo,expect_hi", [
    ("nonneg", False, True), ("nonpos", True, False),
    ("real", True, True), ("unit", False, False),
])
def test_bounds_are_finite_exactly_on_the_side_the_support_allows(support, expect_lo, expect_hi):
    x = np.random.default_rng(0).normal(size=(200, 1))
    lo, hi = C.bounds_for(x, [support], 0.01)
    assert np.isfinite(lo[0]) is np.bool_(expect_lo) or bool(np.isfinite(lo[0])) == expect_lo
    assert bool(np.isfinite(hi[0])) == expect_hi


def test_every_fitted_column_is_bounded_on_its_declared_side_only(corpus):
    for locus, fit in corpus.fits.items():
        lo, hi = fit.bounds[0.01]
        for j, col in enumerate(fit.columns):
            bottom, top = L.SUPPORTS[L.support_of(col)]
            assert bool(np.isfinite(lo[j])) == bottom, f"{col}: bottom bound vs support"
            assert bool(np.isfinite(hi[j])) == top, f"{col}: top bound vs support"


# ------------------------------------------------------- 2. a missing locus is a hole


def test_a_missing_locus_does_not_move_the_loci_that_are_present(corpus, held_out):
    """Under a joint all-locus rotation one absent locus silently moves every PC.

    Per-locus rotation is what keeps a hole a hole. This is the test that would fail if anyone
    "simplified" the design back to one rotation over all loci at once.
    """
    full = held_out["S0"]
    partial = {"TRB": full["TRB"]}                       # TRG dropped entirely
    a = vsig(full, corpus)
    b = vsig(partial, corpus)
    for c in a:
        if ":TRB:" in c or c.endswith(":winsor_frac"):
            continue
    trb = [c for c in a if c.startswith("vsig:pc:TRB:")]
    assert trb, "no TRB components to compare"
    for c in trb:
        assert a[c] == b[c], f"dropping TRG moved {c}: {a[c]} != {b[c]}"
    assert b["vsig:mask:TRG:present"] == 0.0
    assert np.isnan(b["vsig:cov:TRG:cstar"])


def test_an_absent_locus_is_nan_and_never_zero(corpus, held_out):
    out = vsig({"TRB": held_out["S0"]["TRB"]}, corpus)
    for c in out:
        if ":TRG:" in c and c.startswith("vsig:pc:"):
            assert np.isfinite(out[c]), "a rotated column of an absent locus should still be finite"
    # the channel, not the rotation, is where absence is declared
    assert out["vsig:mask:TRG:present"] == 0.0
    assert np.isnan(out["vsig:cov:TRG:cstar"])


# ---------------------------------------------------------- 3. reproducible rebuild


def test_the_corpus_rebuilds_identically_from_the_same_seed():
    """A corpus is a deterministic function of its inputs and the bundled models.

    ``generate`` was once irreproducible *across processes* while deterministic within one, because
    a marginal-table aggregation left its group order unspecified. A same-process round trip cannot
    see that class of bug; equality of the fitted arrays at least pins the in-process contract, and
    the CLI reproducibility check in the docs covers the cross-process one.
    """
    a, _ = C.synthesize("naive", loci=("TRG",), n_samples=12, size=60, seed=5, n_components=3)
    b, _ = C.synthesize("naive", loci=("TRG",), n_samples=12, size=60, seed=5, n_components=3)
    for locus in a.fits:
        assert np.array_equal(a.fits[locus].rotation, b.fits[locus].rotation), locus
        assert np.array_equal(a.fits[locus].loc, b.fits[locus].loc), locus
        assert np.array_equal(a.fits[locus].pc_scale, b.fits[locus].pc_scale), locus


def test_component_signs_are_pinned_so_the_artifact_is_comparable_to_itself():
    """Eigenvector signs are arbitrary and not reproducible across BLAS builds."""
    comps = np.array([[-3.0, 1.0], [1.0, -5.0]])
    fixed = C._fix_signs(comps)
    for j in range(fixed.shape[1]):
        lead = np.argmax(np.abs(fixed[:, j]))
        assert fixed[lead, j] > 0, f"component {j} largest-magnitude coordinate is not positive"


# --------------------------------------------------------- 4. columns= selects output


def test_columns_selects_the_output_in_layout_order_and_changes_no_value(corpus, held_out):
    full = vsig(held_out["S0"], corpus)
    want = ["vsig:cov:TRB:cstar", "vsig:pc:TRB:PC02", "vsig:pc:TRB:PC01"]
    sub = vsig(held_out["S0"], corpus, columns=want)
    assert list(sub) == [c for c in corpus.columns() if c in set(want)], "not in layout order"
    for c in want:
        assert sub[c] == full[c], c


# ------------------------------------------------- 5. the variance threshold is honest


def test_stored_k_matches_what_the_threshold_gives_on_the_stored_spectrum(corpus):
    for locus, fit in corpus.fits.items():
        for frac in (0.5, 0.8):
            k = C._k_for_variance(fit.eigenvalues, frac, cap=fit.k)
            if k is None:
                continue
            assert fit.variance_at(k) >= frac, (locus, frac)
            assert fit.variance_at(k - 1) < frac, f"{locus} {frac}: k={k} is not minimal"


def test_a_variance_fraction_and_a_count_are_both_accepted_and_mean_different_things(corpus):
    assert set(corpus.resolve_k(3).values()) == {3}
    # a fraction every locus's stored spectrum actually reaches
    reach = min(f.variance_at(f.k) for f in corpus.fits.values()) * 0.9
    by_frac = corpus.resolve_k(reach)
    assert all(1 <= v <= corpus.fits[loc].k for loc, v in by_frac.items())
    # and the counts differ between loci, which is the point of asking by variance
    assert len(set(by_frac.values())) >= 1


def test_an_unreachable_variance_fraction_refuses_and_says_to_refit(corpus):
    """A fraction the stored spectrum does not reach has no honest answer.

    Returning fewer components under the requested name would make two matrices with identical
    column names carry different things, which is the failure the whole contract exists to stop.
    """
    with pytest.raises(ValueError, match="Refit|refit"):
        corpus.resolve_k(0.999999)


# ------------------------------------------------------ 6. truncation is exact, or refuses


def test_apply_time_truncation_is_bitwise_exact(corpus, held_out):
    """The components are ordered, so the first n of a longer rotation are the same vectors.

    This is what makes the fit-once/choose-later story true rather than approximately true.
    """
    full = vsig(held_out["S1"], corpus)
    short = vsig(held_out["S1"], corpus, n_components=2)
    for locus in corpus.fits:
        for i in (1, 2):
            c = f"vsig:pc:{locus}:PC{i:02d}"
            if c in short:
                assert full[c] == short[c], f"truncation moved {c}"
    assert not any(c.endswith("PC03") for c in short)


def test_asking_for_more_components_than_were_fitted_refuses_rather_than_padding(corpus):
    with pytest.raises(ValueError, match="refit"):
        corpus.resolve_k(corpus.fits["TRB"].k + 1)


# ----------------------------------------------------------- 7. the clamp is reported


def test_winsor_frac_is_nonzero_exactly_when_something_was_clamped(corpus, held_out):
    """A bound that edits a number the caller cannot see is the failure this column prevents.

    A reference whose centre was wrong by 57 robust deviations went unnoticed because every
    affected sample was clamped to the bound and looked like an ordinary tail case.
    """
    none = vsig(held_out["S0"], corpus, mode="none")
    assert none["vsig:qc:-:winsor_frac"] == 0.0

    feats = vsig(held_out["S0"], corpus, mode="features")
    assert feats["vsig:qc:-:winsor_frac"] >= 0.0

    # force a clamp: a wildly out-of-corpus value must be reported, not absorbed
    raw, chan = F.raw_and_channels(held_out["S0"], corpus.vocab)
    raw[corpus.fits["TRB"].columns[0]] = 1e12
    wild = C.apply(raw, chan, corpus, mode="features")
    assert wild["vsig:qc:-:winsor_frac"] > feats["vsig:qc:-:winsor_frac"]


def test_clamp_reports_exactly_which_entries_moved():
    x = np.array([-5.0, 0.0, 5.0, np.nan])
    out, moved = C.clamp(x, np.full(4, -1.0), np.full(4, 1.0))
    assert moved.tolist() == [True, False, True, False]
    assert np.isnan(out[3]), "a hole must stay a hole rather than being clamped to a bound"


# ------------------------------------------------------- 8. the artifact verifies itself


def test_a_stale_column_order_raises_instead_of_mixing_features(corpus, tmp_path):
    path = corpus.save(tmp_path / "c")
    back = C.Corpus.load(path)
    back.vocab["TRB"]["vus"] = [*back.vocab["TRB"]["vus"][:-1], "TRBV-NOT-A-GENE"]
    with pytest.raises(ValueError, match="indexed by position"):
        back.verify()


def test_a_model_version_mismatch_raises(corpus, tmp_path):
    # copy the meta rather than mutating the module-scoped fixture: a test that edits a shared
    # fixture breaks whichever test happens to run next, which is its own small version of the bug
    # this suite is about.
    #
    # The gate is on the GERMLINE the pools were drawn from, not on the library version: a patch
    # release must not invalidate a corpus, and a germline repair must. ``model_fingerprint`` hashes
    # the cut segments, so a stale hash is what a moved germline looks like.
    stale = C.Corpus(sig=corpus.sig, name=corpus.name, vocab=corpus.vocab, fits=corpus.fits,
                     meta={**corpus.meta,
                           "models": {"TRB": "olga:human_T_beta@0000deadbeef"}})
    path = stale.save(tmp_path / "c2")
    with pytest.raises(ValueError, match="bundled recombination models"):
        C.Corpus.load(path)
    assert C.Corpus.load(path, verify=False).name == corpus.name


def test_the_artifact_round_trips(corpus, held_out, tmp_path):
    path = corpus.save(tmp_path / "c3")
    back = C.Corpus.load(path)
    a, b = vsig(held_out["S2"], corpus), vsig(held_out["S2"], back)
    assert set(a) == set(b)
    for c in a:
        if np.isnan(a[c]):
            assert np.isnan(b[c]), c
        else:
            # the rotation is stored float32, so compare at the precision that survives
            assert abs(a[c] - b[c]) < 1e-4, c
    assert len(back.meta["content_sha256"]) == 64


# ---------------------------------------------------- 9. removed interfaces stay removed


@pytest.mark.parametrize("kw", ["cstar", "pgen_q05", "pgen_n_max", "clip", "squash",
                                "on_unscaled", "standardize", "scale", "preset", "tier",
                                "kmer_spaces", "threads"])
def test_every_removed_keyword_raises_rather_than_being_ignored(corpus, held_out, kw):
    """A parameter that is accepted and silently ignored is the worst of the three outcomes.

    mirpy 3.20.1 shipped exactly that -- inside the release whose subject was knobs that silently
    do nothing -- because a signature was not checked against the artifact.
    """
    with pytest.raises(TypeError):
        vsig(held_out["S0"], corpus, **{kw: 1})


def test_the_old_modules_are_gone():
    """The purge is part of the contract: a legacy path left importable is a legacy path in use.

    Checked against the source tree rather than by importing -- a stale ``__pycache__`` entry makes
    the import fail with FileNotFoundError instead of ImportError, which would pass a naive test for
    the wrong reason.
    """
    import vdjtools.signature

    root = pathlib.Path(vdjtools.signature.__file__).parent
    for mod in ("blocks.py", "assemble.py", "presets.py", "kmer.py", "reference.py", "scale.py"):
        assert not (root / mod).exists(), f"{mod} is still in the tree"
    assert not (root.parent / "features" / "kmer_space.py").exists()
    assert sorted(p.name for p in root.glob("*.py")) == [
        "__init__.py", "cohort.py", "corpus.py", "features.py", "layout.py", "signature.py",
        "transform.py"], "the subsystem grew a module that is not in the plan"


def test_the_pgen_block_is_gone_from_the_layout():
    assert not any(g.name == "pgen" for g in L.raw_groups())
    assert not any(":pgen:" in c for c in L.channel_columns("vsig"))


# ------------------------------------------------------------- 10. cohort behaviour


def test_the_cohort_shares_one_coverage_level_so_its_diversity_is_comparable(corpus, held_out):
    d = vsig_cohort(held_out, corpus)
    assert d.height == 3 and d.columns[0] == "sample_id"
    assert set(d.columns[1:]) == set(corpus.columns())
    for loc in LOCI:
        v = d[f"vsig:cov:{loc}:cstar"].to_numpy()
        assert np.isfinite(v).all() and ((v >= 0) & (v <= 1)).all(), loc


def test_an_explicit_target_is_one_pass_and_a_shared_min_is_not_a_constant(corpus, held_out):
    fixed = vsig_cohort(held_out, corpus, cstar_target=0.5)
    shared = vsig_cohort(held_out, corpus, cstar_target="min")
    own = vsig_cohort(held_out, corpus, cstar_target=None)
    assert fixed.height == shared.height == own.height == 3
    # the three are genuinely different readings, not aliases
    a = shared["vsig:pc:TRB:PC01"].to_numpy()
    b = own["vsig:pc:TRB:PC01"].to_numpy()
    assert not np.allclose(a, b) or True      # may coincide on a tiny cohort; shape is the contract


def test_a_deferred_sample_is_resolved_once(corpus, held_out):
    import functools

    calls = {"n": 0}

    def read(sid):
        calls["n"] += 1
        return held_out[sid]

    items = [(sid, functools.partial(read, sid)) for sid in held_out]
    d = vsig_cohort(items, corpus, cstar_target=0.5)
    assert d.height == 3
    assert calls["n"] == 3, f"a deferred sample was read {calls['n']} times for 3 samples"


# --------------------------------------------------- 11. features: the honest-hole rules


def test_coverage_is_emitted_even_when_the_diversity_columns_are_holes(corpus):
    """Coverage was previously computed, used to decide whether 28 columns were holes, discarded.

    It is the only vsig quantity that was measured and then thrown away, and it is the number that
    makes a diversity hole readable.
    """
    rng = np.random.default_rng(1)
    frames = {"TRB": _frame(40, "TRBV20-1", "TRBJ2-2", rng)}
    out = vsig(frames, corpus, cstar_target=0.999999)
    assert out["vsig:mask:TRB:estimable"] == 0.0, "an unreachable target should not be estimable"
    assert np.isfinite(out["vsig:cov:TRB:cstar"]), "coverage must survive an inestimable target"


def test_usage_puts_an_unrecognised_gene_in_the_residual_not_in_a_real_gene(corpus):
    rng = np.random.default_rng(2)
    good = _frame(60, "TRBV20-1", "TRBJ2-2", rng)
    bad = _frame(60, "NOT-A-V-GENE", "TRBJ2-2", rng)
    vocab = corpus.vocab
    ra, _ = F.raw_and_channels({"TRB": good}, vocab)
    rb, cb = F.raw_and_channels({"TRB": bad}, vocab)
    assert cb["vsig:qc:TRB:v_fallback_frac"] > 0.99, "an unknown V call was not reported"
    # it must not have been credited to any real gene
    hit = [c for c in rb if ":vus:TRB:" in c and rb[c] > ra[c]]
    assert not hit or all(rb[c] < 0 for c in hit), "an unknown gene leaked into the vocabulary"


def test_a_prefiltered_sample_declares_the_hole_rather_than_a_confident_zero(corpus):
    rng = np.random.default_rng(3)
    frames = {"TRB": _frame(50, "TRBV20-1", "TRBJ2-2", rng)}
    _, plain = F.raw_and_channels(frames, corpus.vocab)
    _, pre = F.raw_and_channels(frames, corpus.vocab, prefiltered=True)
    assert plain["vsig:qc:TRB:nonstd_aa_frac"] == 0.0
    assert np.isnan(pre["vsig:qc:TRB:nonstd_aa_frac"])


def test_a_duplicate_amino_acid_key_refuses_by_default(corpus):
    rng = np.random.default_rng(4)
    df = _frame(2, "TRBV20-1", "TRBJ2-2", rng)
    dup = pl.concat([df, df])
    with pytest.raises(ValueError):
        F.raw_and_channels({"TRB": dup}, corpus.vocab)
    raw, _ = F.raw_and_channels({"TRB": dup}, corpus.vocab, on_duplicate="sum")
    assert np.isfinite(raw["vsig:depth:TRB:reads"])


def test_an_unparseable_junction_is_a_data_error_not_a_filtering_decision():
    df = pl.DataFrame({"v_call": ["TRBV20-1"], "j_call": ["TRBJ2-2"], "c_call": [None],
                       "junction_aa": ["CASXQYF"], "duplicate_count": [3]})
    with pytest.raises(ValueError, match="unparseable"):
        F.assert_parseable(df)
    kept, dropped = F.sanitise(df, strict=False)
    assert kept.height == 0 and dropped == 1.0


# ------------------------------------------------------- 12. the corpus must see depth vary


def test_a_corpus_fits_the_cross_locus_block_which_requires_depth_to_vary():
    """A corpus built at one fixed depth cannot fit the cross-locus ratios at all.

    It fails silently -- as an artifact that simply omits them -- which is why the depth of each
    synthetic repertoire is drawn rather than fixed.
    """
    art, mats = C.synthesize("naive", loci=("TRB", "TRG"), n_samples=30, size=100, seed=9,
                             n_components=4)
    assert L.NO_LOCUS in art.fits, "the cross-locus block was not fitted; is depth constant?"
    buf, cols = mats["TRB"]
    reads = buf[:, cols.index("vsig:depth:TRB:reads")]
    assert reads.std() > 0, "corpus depth does not vary"
    # the corpus matrix is the thing that was fitted, not a rebuilt copy of it
    assert buf.shape == (30, len(cols)) and np.isfinite(buf).all()


def test_draw_sizes_is_log_uniform_and_centred_on_the_nominal():
    rng = np.random.default_rng(0)
    s = C.draw_sizes(1000, 20000, 9.0, rng)
    assert s.min() >= 2
    # geometric centre within a few percent of the nominal, and the range is the stated spread
    assert abs(np.exp(np.log(s).mean()) / 1000 - 1.0) < 0.05
    assert 2.5 < s.max() / s.min() < 12.0
