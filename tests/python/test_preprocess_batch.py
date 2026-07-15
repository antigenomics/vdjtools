"""Tests for vdjtools.preprocess.batch — VJ-usage batch-effect correction."""
import math

import polars as pl
import pytest

from vdjtools import preprocess as pp


def _long():
    """2 batches x 2 samples, 2 VJ genes; batches have opposite J-usage bias."""
    rows = []

    def add(sid, batch, counts):
        for (v, j), c in counts.items():
            rows.append(dict(sample_id=sid, batch=batch, v_call=v, j_call=j,
                             cdr3_aa="CASS" + j[-1], cdr3_nt="ACGACG",
                             duplicate_count=c, frequency=0.0))

    add("s1", "A", {("TRBV1", "TRBJ1"): 90, ("TRBV1", "TRBJ2"): 10})
    add("s2", "A", {("TRBV1", "TRBJ1"): 85, ("TRBV1", "TRBJ2"): 15})
    add("s3", "B", {("TRBV1", "TRBJ1"): 12, ("TRBV1", "TRBJ2"): 88})
    add("s4", "B", {("TRBV1", "TRBJ1"): 8, ("TRBV1", "TRBJ2"): 92})
    return pl.DataFrame(rows)


def _batch_divergence(res, col):
    """Sum of |mean_batchA - mean_batchB| per gene for a usage column."""
    m = res.group_by(["batch", "v_call", "j_call"]).agg(pl.col(col).mean().alias("m"))
    piv = m.pivot(values="m", index=["v_call", "j_call"], on="batch")
    return float((piv["A"] - piv["B"]).abs().sum())


def test_batch_correction_shrinks_divergence():
    res = pp.correct_vj_usage(_long(), batch_col="batch")
    before = _batch_divergence(res, "p")
    after = _batch_divergence(res, "p_corrected")
    assert before > 1.0                            # strong batch bias in raw usage
    assert after < before * 0.05                   # correction removes most of it
    assert after < 1e-3


def test_batch_corrected_usage_is_normalized():
    res = pp.correct_vj_usage(_long(), batch_col="batch")
    tot = res.group_by(["sample_id", "locus"]).agg(pl.col("p_corrected").sum().alias("s"))
    assert all(abs(v - 1.0) < 1e-9 for v in tot["s"].to_list())   # per-sample probabilities


def test_batch_preserves_within_batch_variation():
    # s1 and s2 are both batch A but differ; correction must not collapse them.
    res = pp.correct_vj_usage(_long(), batch_col="batch")
    s1 = res.filter((pl.col("sample_id") == "s1") & (pl.col("j_call") == "TRBJ1"))["p_corrected"][0]
    s2 = res.filter((pl.col("sample_id") == "s2") & (pl.col("j_call") == "TRBJ1"))["p_corrected"][0]
    assert s1 != s2


def test_batch_corrected_value_pins_grand_mean():
    # Value pin that guards the `+ mu_grand` restore term (step 4). With the term,
    # s1/TRBJ1 corrects to 0.524859 (plain-mean default); dropping `+ mu_grand` shifts
    # it, so this assertion fails if the grand-mean restore is ever removed.
    res = pp.correct_vj_usage(_long(), batch_col="batch")
    row = res.filter((pl.col("sample_id") == "s1") & (pl.col("j_call") == "TRBJ1"))
    assert math.isclose(row["p_corrected"][0], 0.524859, rel_tol=1e-5)


def test_batch_accepts_list_of_frames():
    long = _long()
    frames = [long.filter(pl.col("sample_id") == s) for s in ("s1", "s2", "s3", "s4")]
    res = pp.correct_vj_usage(frames, batch_col="batch")
    assert res.height == 8                          # 4 samples x 2 genes


# --------------------------------------------------------------------------------------
# transform="sigmoid" (Vlasova et al. 2026: z-score + grand-mean-preserving sigmoid) and
# apply_vj_correction (rescale + roulette-wheel resample of the clonotype table).
# --------------------------------------------------------------------------------------

def _cohort(n_per_batch=6, reads=2000, seed=1):
    """A realistic cohort: n samples/batch, opposite J bias, per-sample jitter."""
    import numpy as np
    rng = np.random.default_rng(seed)
    rows = []
    for batch, base_j1 in (("A", 0.80), ("B", 0.20)):
        for k in range(n_per_batch):
            j1 = float(np.clip(base_j1 + rng.normal(0, 0.04), 0.02, 0.98))
            c1 = int(reads * j1)
            for (v, j), c in {("TRBV1", "TRBJ1"): c1, ("TRBV1", "TRBJ2"): reads - c1}.items():
                rows.append(dict(sample_id=f"{batch}{k}", batch=batch, v_call=v, j_call=j,
                                 junction_aa="CASS" + j[-1], junction_nt="ACG",
                                 duplicate_count=c, frequency=0.0))
    return pl.DataFrame(rows)


def test_sigmoid_rejects_unknown_transform():
    with pytest.raises(ValueError, match="transform must be"):
        pp.correct_vj_usage(_long(), batch_col="batch", transform="bogus")


def test_sigmoid_removes_divergence_and_preserves_grand_mean():
    long = _cohort()
    res = pp.correct_vj_usage(long, batch_col="batch", transform="sigmoid")
    # 1) batch divergence collapses
    assert _batch_divergence(res, "p_corrected") < _batch_divergence(res, "p") * 0.05
    # 2) per-(sample, locus) normalisation
    tot = res.group_by(["sample_id", "locus"]).agg(pl.col("p_corrected").sum().alias("s"))
    assert all(abs(v - 1.0) < 1e-9 for v in tot["s"].to_list())
    # 3) grand-mean preserved: mean corrected usage across samples == pooled raw P_avg
    pooled_j1 = (long.filter(pl.col("j_call") == "TRBJ1")["duplicate_count"].sum()
                 / long["duplicate_count"].sum())
    mean_corr_j1 = res.filter(pl.col("j_call") == "TRBJ1")["p_corrected"].mean()
    assert abs(mean_corr_j1 - pooled_j1) < 0.02


def test_sigmoid_value_pin():
    # Guards the 2*P_avg/(1+exp(-Z)) formula with the plain-mean/σ default; s1/TRBJ1 ->
    # 0.658608 (differs from the location path's 0.524859, so this fails if the sigmoid
    # map is ever altered).
    res = pp.correct_vj_usage(_long(), batch_col="batch", transform="sigmoid")
    row = res.filter((pl.col("sample_id") == "s1") & (pl.col("j_call") == "TRBJ1"))
    assert math.isclose(row["p_corrected"][0], 0.658608, rel_tol=1e-5)


def test_winsor_q_knob_matches_legacy_winsorized_values():
    # winsor_q=0.025 restores the legacy winsorized mean/σ (the noisy-features regime);
    # the default winsor_q=None is the paper's plain mean/σ. Guards both code paths.
    loc = pp.correct_vj_usage(_long(), batch_col="batch", winsor_q=0.025)
    sig = pp.correct_vj_usage(_long(), batch_col="batch", transform="sigmoid", winsor_q=0.025)

    def s1j1(res):
        return res.filter((pl.col("sample_id") == "s1") & (pl.col("j_call") == "TRBJ1"))["p_corrected"][0]

    assert math.isclose(s1j1(loc), 0.524768, rel_tol=1e-5)
    assert math.isclose(s1j1(sig), 0.666926, rel_tol=1e-5)


def _sample_frame(long, sid):
    """Extract one sample's canonical clonotype frame from a cohort long frame."""
    from vdjtools.io import schema as S
    df = long.filter(pl.col("sample_id") == sid).drop(["sample_id", "batch"])
    return S.normalize(df, recompute_freq=True)


def test_apply_resample_preserves_total_and_corrects_usage():
    long = _cohort()
    cu = pp.correct_vj_usage(long, batch_col="batch", transform="sigmoid")
    a0 = _sample_frame(long, "A0")
    out = pp.apply_vj_correction(a0, cu.filter(pl.col("sample_id") == "A0"), seed=0)
    # total reads preserved by the multinomial roulette wheel
    assert int(out["duplicate_count"].sum()) == int(a0["duplicate_count"].sum())
    # J1 usage moves from the batch-biased raw value toward the corrected target
    raw = a0.filter(pl.col("j_call") == "TRBJ1")["duplicate_count"].sum() / a0["duplicate_count"].sum()
    new = out.filter(pl.col("j_call") == "TRBJ1")["duplicate_count"].sum() / out["duplicate_count"].sum()
    target = cu.filter((pl.col("sample_id") == "A0") & (pl.col("j_call") == "TRBJ1"))["p_corrected"][0]
    assert abs(new - target) < abs(raw - target)          # closer to target than raw
    assert out.columns == a0.columns                       # canonical schema preserved
    # seeded → reproducible
    assert out.equals(pp.apply_vj_correction(a0, cu.filter(pl.col("sample_id") == "A0"), seed=0))


def test_apply_resample_false_is_deterministic_expected_counts():
    long = _cohort()
    cu = pp.correct_vj_usage(long, batch_col="batch", transform="sigmoid")
    a0 = _sample_frame(long, "A0")
    det = pp.apply_vj_correction(a0, cu.filter(pl.col("sample_id") == "A0"), resample=False)
    assert int(det["duplicate_count"].sum()) == int(a0["duplicate_count"].sum())   # exact total
    j1 = det.filter(pl.col("j_call") == "TRBJ1")["duplicate_count"].sum() / det["duplicate_count"].sum()
    target = cu.filter((pl.col("sample_id") == "A0") & (pl.col("j_call") == "TRBJ1"))["p_corrected"][0]
    assert abs(j1 - target) < 1e-3                          # expected counts hit the target


def test_apply_requires_sample_id_for_multi_sample_usage():
    long = _cohort()
    cu = pp.correct_vj_usage(long, batch_col="batch", transform="sigmoid")   # 12 samples
    a0 = _sample_frame(long, "A0")
    with pytest.raises(ValueError, match="multiple samples"):
        pp.apply_vj_correction(a0, cu)                      # ambiguous — needs sample_id=
    out = pp.apply_vj_correction(a0, cu, sample_id="A0")    # explicit selection works
    assert int(out["duplicate_count"].sum()) == int(a0["duplicate_count"].sum())
