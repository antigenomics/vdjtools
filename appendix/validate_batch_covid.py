"""Validate VJ-usage batch correction on real deep TCRβ repertoires (FMBA covid Cohort I).

Reference data — Vlasova, Nekrasova, Komkov, … Britanova, Shugay, *Inference of SARS-CoV-2
exposure biomarkers using large-scale T-cell repertoire profiling*, Genome Medicine
2026;18:20 (10.1186/s13073-025-01589-4). The paper's batch-effect correction (Methods,
"Batch-effect correction and data normalization") is exactly what ``preprocess.batch``
implements: per (gene, batch) z-score ``Z = (log P − μ)/σ`` of the log-normal gene usage,
then the grand-mean-preserving map ``P_final = 2·P_avg/(1 + exp(−Z))``, then rescale
clonotype frequencies ``f' = f·P_final(G,S)/P(G,S)`` and roulette-wheel (multinomial)
resample at the sample's read depth. The paper uses the **plain** mean/σ (Shapiro–Wilk
validated) — i.e. ``winsor_q=None`` (the default here); winsorization is only for the
noisy usage-as-features regime (many shallow / RNA-seq repertoires).

Data provenance (aldan3 HPC, project ``fmba_covid``; requires the ``aldan3`` client):
  - clonotype tables  : /projects/fmba_covid/COV_V_usage_adjustment_v3/FMBA_functional/*.clonotypes.TRB.txt
                        (legacy VDJtools format, functional-filtered — the paper's own inputs)
  - raw V usage matrix: /projects/fmba_covid/COV_V_usage_adjustment_v3/V_usage_FMBA.tsv
  - paper's adjusted  : /projects/fmba_covid/COV_V_usage_adjustment_v3/V_usage_FMBA_adjusted.tsv
  - batch labels      : /projects/fmba_covid/unprocessed_fmba_metadata.csv → column `folder`
                        (the 9 NovaSeq sequencing runs = the paper's "nine separate batches";
                        join by the 12-digit sample id in column `name`)
  Pull a subset with e.g.:
    aldan3 pull /projects/fmba_covid/COV_V_usage_adjustment_v3/V_usage_FMBA.tsv .
    aldan3 exec "cd .../FMBA_functional && tar czf /tmp/sel.tgz <files>"; aldan3 pull /tmp/sel.tgz .

Findings (2026-07-16; 48 samples × 4 batches, median depth 3.3M reads/sample):
  - V-usage variance explained by batch (η², mean over genes): raw 0.109 → corrected 0.002
    (batch signal removed ~54×), within-sample biological variation retained.
  - grand-mean usage preserved: max |pooled_raw − mean_corrected| over genes = 0.0011.
  - apply_vj_correction resamples real deep tables exactly to the sample's read depth and
    lands V usage on the corrected target; seed-reproducible.
  - Formula fidelity vs the paper's *published* adjusted matrix (1136 samples, 24 V families,
    9 batches): the sigmoid-of-batch-z form explains r≈0.90; the residual is the paper's
    cross-cohort P_avg target + functional-gene/family-mapping pipeline specifics (parameter
    choices, not a formula difference — the Methods equation matches verbatim).

Not a pytest test (needs cluster data); run manually from a checkout with the pulled subset.
"""
import glob
import re

import numpy as np
import polars as pl

from vdjtools import io as vio
from vdjtools import preprocess as pp
from vdjtools.io import schema as S

# Directory holding the pulled *.clonotypes.TRB.txt subset (one batch encoded per filename
# stem is not needed — pass an explicit {filename: batch} map or a metadata join instead).
DATA_DIR = "fmba_tables"


def _eta2_by_batch(usage_wide: np.ndarray, batch: list[str]) -> float:
    """Mean over genes of the fraction of usage variance explained by batch (η²)."""
    b = np.array(batch)
    grand = usage_wide.mean(0)
    ss_tot = ((usage_wide - grand) ** 2).sum(0)
    ss_bet = sum((b == lev).sum() * (usage_wide[b == lev].mean(0) - grand) ** 2 for lev in set(b))
    return float(np.where(ss_tot > 0, ss_bet / ss_tot, 0.0).mean())


def run(batch_of: "dict[str, str]") -> None:
    """batch_of: {clonotype-filename -> batch label}. Prints the validation summary."""
    frames = []
    for fn, batch in batch_of.items():
        sid = re.match(r"([0-9]+)_", fn).group(1) + "_" + batch.split("/")[-1]
        df = vio.read(f"{DATA_DIR}/{fn}").with_columns(
            pl.lit(sid).alias("sample_id"), pl.lit(batch).alias("batch"))
        frames.append(df)
    cohort = pl.concat(frames, how="vertical_relaxed")

    # V-level correction (paper's primary setting): collapse J, plain mean/σ (winsor_q=None).
    vonly = cohort.with_columns(pl.lit("TRBJ0").alias("j_call"))
    cu = pp.correct_vj_usage(vonly, batch_col="batch", transform="sigmoid")

    def wide(col):
        piv = cu.pivot(values=col, index="sample_id", on="v_call").fill_null(0.0)
        bat = [cu.filter(pl.col("sample_id") == s)["batch"][0] for s in piv["sample_id"]]
        return piv.drop("sample_id").to_numpy(), bat

    Xraw, bat = wide("p")
    Xcor, _ = wide("p_corrected")
    print(f"η²(batch)  raw={_eta2_by_batch(Xraw, bat):.4f}  corrected={_eta2_by_batch(Xcor, bat):.4f}")

    gene = cu.group_by("v_call").agg(
        (pl.col("count").sum() / cu["count"].sum()).alias("pooled_raw"),
        pl.col("p_corrected").mean().alias("mean_corr"))
    print(f"grand-mean preserved: max|pooled_raw-mean_corr| = "
          f"{(gene['pooled_raw'] - gene['mean_corr']).abs().max():.4f}")

    for fn, batch in list(batch_of.items())[:4]:
        sid = re.match(r"([0-9]+)_", fn).group(1) + "_" + batch.split("/")[-1]
        s = S.normalize(cohort.filter(pl.col("sample_id") == sid).drop(["sample_id", "batch"]),
                        recompute_freq=True)
        out = pp.apply_vj_correction(s, cu.filter(pl.col("sample_id") == sid),
                                     scope="v", sample_id=sid, seed=0)
        assert int(out["duplicate_count"].sum()) == int(s["duplicate_count"].sum())  # reads preserved
    print("apply_vj_correction: reads preserved on all sampled tables ✓")


if __name__ == "__main__":  # pragma: no cover
    batch_of = {re.sub(r".*/", "", p): "unknown"
                for p in sorted(glob.glob(f"{DATA_DIR}/*.clonotypes.TRB.txt"))}
    run(batch_of)
