"""VJ-usage batch-effect correction + clonotype-table application (pure polars/numpy).

Different sequencing batches carry systematic V/J gene-usage biases (primer mixes,
amplification, extraction). This module removes that batch-specific offset so that
per-sample VJ usage becomes comparable across batches, and can then push the corrected
usage back onto a sample's clonotype table (rescale + resample).

Two stages
----------
1. :func:`correct_vj_usage` — batch-correct per-sample V-J gene usage. Two transforms:

   - ``transform="location"`` (default) — the classic **location adjustment** (the
     location term of ComBat; Johnson, Li & Rabinovic, *Biostatistics* 2007) on
     gene-usage log-probabilities. Per ``(locus, gene, batch)`` the winsorized batch
     mean ``mu_batch`` of ``log p`` is replaced by the winsorized grand mean
     ``mu_grand``: ``log_corrected = log_p - mu_batch + mu_grand``, then ``exp`` and
     renormalise. Location only — no scale term.
   - ``transform="sigmoid"`` — the **σ-standardised, grand-mean-preserving** correction
     of Vlasova, Nekrasova, Komkov, … Britanova, Shugay, *Genome Medicine* 2026;18:20
     (DOI 10.1186/s13073-025-01589-4). Per ``(locus, gene, batch)`` compute a winsorized
     **z-score** ``Z = (log p - mu_batch) / sigma_batch`` (capped at ``±z_cap``), then map
     it back to a probability with a sigmoid that preserves the pooled grand-mean usage
     ``P_avg(gene)``: ``P_final = 2 * P_avg / (1 + exp(-Z))`` (``Z=0`` → ``P_avg``),
     renormalised per ``(sample, locus)``. Reference: legacy mirpy v2
     ``mir.basic.gene_usage.compute_batch_corrected_gene_usage`` (winsorized ``mu``/``sigma``,
     z-cap, pooled ``P_avg``) — note that legacy's own ``pfinal`` uses ``p*exp(Z)``; the
     ``2*P_avg*sigmoid(Z)`` map here is the paper's Methods formula.

2. :func:`apply_vj_correction` — apply a sample's corrected usage back to its clonotype
   table: reweight each clonotype by ``P_final(G) / P(G)`` for its gene, then (by default)
   roulette-wheel **resample** to a new integer-count table with the corrected usage.
   Port of legacy mirpy v2 ``mir.common.sampling.resample_to_gene_usage``.

Alleles are stripped for the usage key (mirpy convention); clonotypes keep their original
allele calls.
"""
from __future__ import annotations

import numpy as np
import polars as pl

from ..io.schema import COUNT, J_CALL, V_CALL, recompute_frequency, strip_allele


def _winsorized_mean(value: str, group: "list[str]",
                     lower_q: float = 0.025, upper_q: float = 0.975) -> pl.Expr:
    """Windowed winsorized mean of ``value`` over ``group`` (linear quantiles)."""
    lo = pl.col(value).quantile(lower_q, interpolation="linear").over(group)
    hi = pl.col(value).quantile(upper_q, interpolation="linear").over(group)
    return pl.col(value).clip(lo, hi).mean().over(group)


def _winsorized_std(value: str, group: "list[str]",
                    lower_q: float = 0.025, upper_q: float = 0.975) -> pl.Expr:
    """Windowed winsorized sample SD (ddof=1) of ``value`` over ``group``; null → 0.

    Matches legacy ``_winsorized_mean_std`` (``np.std(ddof=1)``, ``0.0`` for a
    single-element group).
    """
    lo = pl.col(value).quantile(lower_q, interpolation="linear").over(group)
    hi = pl.col(value).quantile(upper_q, interpolation="linear").over(group)
    return pl.col(value).clip(lo, hi).std(ddof=1).over(group).fill_null(0.0)


def correct_vj_usage(samples_or_df: "pl.DataFrame | list[pl.DataFrame]", batch_col: str,
                     sample_col: str = "sample_id", weighted: bool = True,
                     pseudocount: float = 1.0, transform: str = "location",
                     z_cap: float = 6.0) -> pl.DataFrame:
    """Batch-correct per-sample V-J gene usage.

    Args:
        samples_or_df: A single long clonotype frame carrying ``sample_col`` and
            ``batch_col`` columns, or a list of such frames (concatenated).
        batch_col: Column naming each sample's batch.
        sample_col: Column naming the sample (default ``"sample_id"``).
        weighted: If ``True`` (default) usage counts reads (``duplicate_count``); if
            ``False`` it counts clonotypes.
        pseudocount: Laplace smoothing constant added per gene (default ``1.0``).
        transform: ``"location"`` (default) for the ComBat location adjustment
            (``log_p - mu_batch + mu_grand``), or ``"sigmoid"`` for the σ-standardised,
            grand-mean-preserving z-score map ``P_final = 2*P_avg/(1+exp(-Z))``
            (Vlasova et al. 2026). See the module docstring.
        z_cap: For ``transform="sigmoid"``, clip the z-score to ``±z_cap`` (default
            ``6.0``); ignored for ``"location"``.

    Returns:
        A long frame with one row per ``(sample, locus, v_call, j_call)`` over the
        union of genes per locus, with columns: ``sample_id``, ``batch``, ``locus``,
        ``v_call``, ``j_call``, ``count``, ``p`` (smoothed raw usage probability) and
        ``p_corrected`` (batch-corrected, renormalised usage probability), sorted by
        ``sample_id, locus, v_call, j_call``. Feed this straight into
        :func:`apply_vj_correction`.

    Raises:
        ValueError: If ``sample_col`` / ``batch_col`` is missing or ``transform`` is
            unknown.
    """
    if transform not in ("location", "sigmoid"):
        raise ValueError(f"transform must be 'location' or 'sigmoid'; got {transform!r}")
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
    sample_loci = usage.select([sample_col, batch_col, "locus"]).unique()
    grid = sample_loci.join(genes, on="locus", how="inner")
    full = (
        grid.join(usage, on=[sample_col, batch_col, "locus", V_CALL, J_CALL], how="left")
        .with_columns(pl.col("count").fill_null(0))
    )
    full = full.with_columns(
        (pl.col("count") + pseudocount).sum().over([sample_col, "locus"]).alias("_denom_raw")
    )
    # denom = total + pseudocount * n_genes == sum(count + pseudocount) over the grid
    full = full.with_columns(
        ((pl.col("count") + pseudocount) / pl.col("_denom_raw")).alias("p")
    ).with_columns(pl.col("p").log().alias("log_p"))

    if transform == "location":
        full = full.with_columns(
            _winsorized_mean("log_p", ["locus", V_CALL, J_CALL, batch_col]).alias("mu_batch"),
            _winsorized_mean("log_p", ["locus", V_CALL, J_CALL]).alias("mu_grand"),
        )
        full = full.with_columns(
            (pl.col("log_p") - pl.col("mu_batch") + pl.col("mu_grand")).exp().alias("_pc_raw")
        )
    else:  # sigmoid: z-score (with sigma) + grand-mean-preserving sigmoid map
        full = full.with_columns(
            _winsorized_mean("log_p", ["locus", V_CALL, J_CALL, batch_col]).alias("mu_batch"),
            _winsorized_std("log_p", ["locus", V_CALL, J_CALL, batch_col]).alias("sigma_batch"),
            # P_avg(gene) = pooled usage across all samples (raw counts, no pseudocount).
            (pl.col("count").sum().over(["locus", V_CALL, J_CALL])
             / pl.col("count").sum().over(["locus"])).alias("p_avg"),
        )
        z = pl.when(pl.col("sigma_batch") > 0) \
              .then((pl.col("log_p") - pl.col("mu_batch")) / pl.col("sigma_batch")) \
              .otherwise(0.0).clip(-z_cap, z_cap)
        full = full.with_columns(
            (2.0 * pl.col("p_avg") / (1.0 + (-z).exp())).alias("_pc_raw")
        )

    # Renormalise per (sample, locus); fall back to raw p if the group mass is invalid.
    grp_sum = pl.col("_pc_raw").sum().over([sample_col, "locus"])
    full = full.with_columns(
        pl.when(grp_sum > 0).then(pl.col("_pc_raw") / grp_sum)
        .otherwise(pl.col("p") / pl.col("p").sum().over([sample_col, "locus"]))
        .alias("p_corrected")
    )

    return (
        full.rename({sample_col: "sample_id", batch_col: "batch"})
        .select(["sample_id", "batch", "locus", V_CALL, J_CALL, "count", "p", "p_corrected"])
        .sort(["sample_id", "locus", V_CALL, J_CALL])
    )


def _gene_factors(corrected_usage: pl.DataFrame, scope: str) -> "tuple[list[str], pl.DataFrame]":
    """Build the per-gene correction factor ``P_corrected / P`` from a corrected-usage frame.

    Returns ``(key_cols, factor_frame)`` where ``factor_frame`` has ``key_cols`` + ``factor``.
    For ``scope="v"`` / ``"j"`` the VJ-joint usage is marginalised (summed over the other gene)
    before taking the ratio.
    """
    if scope == "vj":
        key = ["locus", V_CALL, J_CALL]
        fac = corrected_usage.group_by(key).agg(
            pl.col("p").sum().alias("_p"), pl.col("p_corrected").sum().alias("_pc"))
    elif scope in ("v", "j"):
        gene = V_CALL if scope == "v" else J_CALL
        key = ["locus", gene]
        fac = corrected_usage.group_by(key).agg(
            pl.col("p").sum().alias("_p"), pl.col("p_corrected").sum().alias("_pc"))
    else:
        raise ValueError(f"scope must be 'v', 'j', or 'vj'; got {scope!r}")
    fac = fac.with_columns(
        pl.when(pl.col("_p") > 0).then(pl.col("_pc") / pl.col("_p")).otherwise(0.0).alias("factor")
    ).select([*key, "factor"])
    return key, fac


def apply_vj_correction(sample_df: pl.DataFrame, corrected_usage: pl.DataFrame, *,
                        scope: str = "vj", weighted: bool = True, resample: bool = True,
                        sample_id: str | None = None, seed: int = 0) -> pl.DataFrame:
    """Apply batch-corrected V/J usage back to a sample's clonotype table.

    Each clonotype is reweighted by its gene's correction factor ``P_corrected(G) / P(G)``
    (``G`` = the clonotype's V/J gene at ``scope``), then either resampled or rescaled to a
    new clonotype table whose V/J usage matches the corrected usage. Port of legacy mirpy v2
    ``resample_to_gene_usage``.

    Args:
        sample_df: One sample's canonical clonotype frame (``v_call``, ``j_call``,
            ``duplicate_count``). Original allele calls are preserved on output.
        corrected_usage: Output of :func:`correct_vj_usage`. If it covers more than one
            sample, pass ``sample_id`` (or slice it) to select this sample's rows.
        scope: Gene scope for the correction key: ``"vj"`` (default, per V-J pair),
            ``"v"``, or ``"j"`` (the VJ usage is marginalised for ``"v"``/``"j"``).
        weighted: If ``True`` (default) the resampling weight is ``duplicate_count *
            factor`` (roulette-wheel over reads); if ``False`` it is ``factor`` (over
            clonotypes).
        resample: If ``True`` (default), multinomial roulette-wheel resample to a new
            integer-count table at the sample's original total read count (Vlasova et al.
            2026). If ``False``, deterministically rescale to the expected counts.
        sample_id: The sample to select from ``corrected_usage`` when it holds several.
        seed: Seed for the multinomial resample (``resample=True``).

    Returns:
        A canonical clonotype frame (same columns as ``sample_df``) with corrected
        ``duplicate_count`` and recomputed ``frequency``; zero-count clonotypes dropped.

    Raises:
        ValueError: If ``corrected_usage`` covers multiple samples and none is selected,
            if ``scope`` is unknown, or if no clonotype has a positive corrected weight.
    """
    cu = corrected_usage
    if "sample_id" in cu.columns and cu["sample_id"].n_unique() > 1:
        if sample_id is None:
            raise ValueError(
                "corrected_usage covers multiple samples; pass sample_id= or a single-sample slice")
        cu = cu.filter(pl.col("sample_id") == sample_id)
    if cu.height == 0:
        raise ValueError("corrected_usage has no rows for the requested sample")

    key, fac = _gene_factors(cu, scope)

    # Match on stripped-allele gene keys (+ locus), keeping the original rows/alleles.
    # left helper cols mirror the corrected-usage key names, so join left_on->right_on(key).
    left_map = {"locus": "_locus", V_CALL: "_vk", J_CALL: "_jk"}
    work = sample_df.with_row_index("_row").with_columns(
        strip_allele(pl.col(V_CALL)).alias("_vk"),
        strip_allele(pl.col(J_CALL)).alias("_jk"),
        pl.col(V_CALL).str.slice(0, 3).alias("_locus"),
    )
    work = work.join(fac, left_on=[left_map[k] for k in key], right_on=key,
                     how="left").sort("_row")
    factor = work["factor"].fill_null(0.0).to_numpy().astype(float)  # gene absent -> weight 0

    counts = sample_df[COUNT].to_numpy().astype(float)
    weight = factor * counts if weighted else factor
    total_w = weight.sum()
    if total_w <= 0:
        raise ValueError("no clonotype has a positive corrected weight (gene keys did not match)")
    prob = weight / total_w

    total_reads = int(sample_df[COUNT].sum())
    if resample:
        new_counts = np.random.default_rng(seed).multinomial(total_reads, prob)
    else:
        # Deterministic expected counts, largest-remainder rounded to preserve the total.
        exact = total_reads * prob
        floor = np.floor(exact).astype(np.int64)
        rem = total_reads - int(floor.sum())
        if rem > 0:
            take = np.argsort(-(exact - floor))[:rem]
            floor[take] += 1
        new_counts = floor

    out = (
        sample_df.with_columns(pl.Series(COUNT, new_counts, dtype=pl.Int64))
        .filter(pl.col(COUNT) > 0)
    )
    return recompute_frequency(out)
