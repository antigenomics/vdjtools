"""Tests for vdjtools.preprocess.batch — VJ-usage batch-effect correction."""
import math

import polars as pl

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
    # s1/TRBJ1 corrects to 0.524768; dropping `+ mu_grand` shifts it to 0.553693,
    # so this assertion fails if the grand-mean restore is ever removed.
    res = pp.correct_vj_usage(_long(), batch_col="batch")
    row = res.filter((pl.col("sample_id") == "s1") & (pl.col("j_call") == "TRBJ1"))
    assert math.isclose(row["p_corrected"][0], 0.524768, rel_tol=1e-5)


def test_batch_accepts_list_of_frames():
    long = _long()
    frames = [long.filter(pl.col("sample_id") == s) for s in ("s1", "s2", "s3", "s4")]
    res = pp.correct_vj_usage(frames, batch_col="batch")
    assert res.height == 8                          # 4 samples x 2 genes
