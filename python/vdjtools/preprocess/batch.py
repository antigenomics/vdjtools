"""VJ-usage batch-effect correction (pure polars).

Different sequencing batches carry systematic V/J gene-usage biases (primer mixes,
amplification, extraction). This module removes that batch-specific offset so that
per-sample VJ usage becomes comparable across batches.

Method
------
The routine reimplements — fresh, in polars — the log-space, batch-statistics idea
behind mirpy's ``mir.basic.gene_usage.compute_batch_corrected_gene_usage`` (which
computes per-``(locus, gene, batch)`` winsorized statistics of ``log p`` and turns
them into z-scores). Here we apply the classic **location adjustment** (the
location term of ComBat; Johnson, Li & Rabinovic, *Biostatistics* 2007) directly on
gene-usage log-probabilities:

1. Per ``(sample, locus)`` build a Laplace-smoothed VJ-usage probability over the
   *union* of genes seen in that locus (so genes absent from a sample still get a
   pseudocount): ``p = (count + pseudocount) / (total + pseudocount * n_genes)``.
2. ``log_p = ln p``.
3. For each ``(locus, gene)`` compute the winsorized (2.5/97.5%) **batch mean**
   ``mu_batch`` over samples of the same batch, and the winsorized **grand mean**
   ``mu_grand`` over all samples (mirpy's robust log-space statistic).
4. Remove the batch offset and restore the grand level:
   ``log_corrected = log_p - mu_batch + mu_grand``.
5. ``p_corrected = exp(log_corrected)``, renormalised per ``(sample, locus)`` to sum
   to 1.

Because every batch's per-gene mean is mapped onto the shared grand mean, a
systematic batch difference in a gene's usage is removed while within-batch,
between-sample variation is preserved. Only the location is adjusted; the scale
(``sigma``) term is intentionally left out — with the handful of samples per batch
typical here, a per-gene variance estimate is too noisy to divide by safely.
Alleles are stripped for the usage key (mirpy convention).
"""
from __future__ import annotations

import polars as pl

from ..io.schema import COUNT, J_CALL, V_CALL, strip_allele


def _winsorized_mean(value: str, group: "list[str]",
                     lower_q: float = 0.025, upper_q: float = 0.975) -> pl.Expr:
    """Windowed winsorized mean of ``value`` over ``group`` (linear quantiles)."""
    lo = pl.col(value).quantile(lower_q, interpolation="linear").over(group)
    hi = pl.col(value).quantile(upper_q, interpolation="linear").over(group)
    return pl.col(value).clip(lo, hi).mean().over(group)


def correct_vj_usage(samples_or_df: "pl.DataFrame | list[pl.DataFrame]", batch_col: str,
                     sample_col: str = "sample_id", weighted: bool = True,
                     pseudocount: float = 1.0) -> pl.DataFrame:
    """Batch-correct per-sample V-J gene usage by log-space location adjustment.

    Args:
        samples_or_df: A single long clonotype frame carrying ``sample_col`` and
            ``batch_col`` columns, or a list of such frames (concatenated).
        batch_col: Column naming each sample's batch.
        sample_col: Column naming the sample (default ``"sample_id"``).
        weighted: If ``True`` (default) usage counts reads (``duplicate_count``); if
            ``False`` it counts clonotypes.
        pseudocount: Laplace smoothing constant added per gene (default ``1.0``).

    Returns:
        A long frame with one row per ``(sample, locus, v_call, j_call)`` over the
        union of genes per locus, with columns: ``sample_id``, ``batch``,
        ``locus``, ``v_call``, ``j_call``, ``count``, ``p`` (smoothed raw usage
        probability) and ``p_corrected`` (batch-corrected, renormalised usage
        probability), sorted by ``sample_id, locus, v_call, j_call``.

    Raises:
        ValueError: If ``sample_col`` or ``batch_col`` is missing.
    """
    long = (pl.concat(samples_or_df, how="vertical_relaxed")
            if isinstance(samples_or_df, list) else samples_or_df)
    for c in (sample_col, batch_col):
        if c not in long.columns:
            raise ValueError(f"missing required column {c!r}")

    long = long.with_columns(
        strip_allele(pl.col(V_CALL)).alias(V_CALL),
        strip_allele(pl.col(J_CALL)).alias(J_CALL),
        pl.col(V_CALL).str.slice(0, 3).alias("locus"),
    )
    count_expr = pl.col(COUNT).sum() if weighted else pl.len()
    usage = (
        long.group_by([sample_col, batch_col, "locus", V_CALL, J_CALL])
        .agg(count_expr.alias("count"))
    )

    # Full (sample x gene) grid within each locus so absent genes get a pseudocount.
    genes = usage.select(["locus", V_CALL, J_CALL]).unique()
    n_genes = genes.group_by("locus").agg(pl.len().alias("n_genes"))
    sample_loci = usage.select([sample_col, batch_col, "locus"]).unique()
    grid = sample_loci.join(genes, on="locus", how="inner")
    full = (
        grid.join(usage, on=[sample_col, batch_col, "locus", V_CALL, J_CALL], how="left")
        .with_columns(pl.col("count").fill_null(0))
        .join(n_genes, on="locus", how="left")
    )
    full = full.with_columns(
        (pl.col("count") + pseudocount).sum().over([sample_col, "locus"]).alias("_denom_raw")
    )
    # denom = total + pseudocount * n_genes == sum(count + pseudocount) over the grid
    full = full.with_columns(
        ((pl.col("count") + pseudocount) / pl.col("_denom_raw")).alias("p")
    ).with_columns(pl.col("p").log().alias("log_p"))

    full = full.with_columns(
        _winsorized_mean("log_p", ["locus", V_CALL, J_CALL, batch_col]).alias("mu_batch"),
        _winsorized_mean("log_p", ["locus", V_CALL, J_CALL]).alias("mu_grand"),
    )
    full = full.with_columns(
        (pl.col("log_p") - pl.col("mu_batch") + pl.col("mu_grand")).exp().alias("_pc_raw")
    )
    full = full.with_columns(
        (pl.col("_pc_raw") / pl.col("_pc_raw").sum().over([sample_col, "locus"])).alias("p_corrected")
    )

    return (
        full.rename({sample_col: "sample_id", batch_col: "batch"})
        .select(["sample_id", "batch", "locus", V_CALL, J_CALL, "count", "p", "p_corrected"])
        .sort(["sample_id", "locus", V_CALL, J_CALL])
    )
