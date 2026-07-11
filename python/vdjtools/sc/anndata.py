"""Bridge paired single-cell receptors into an AnnData / MuData container.

Single-cell is the shape where AnnData fits: ``obs`` is one row per observation and
here an observation is a **receptor pair** (``cell_id`` + paired ``alpha_*`` /
``beta_*`` calls), so the whole scverse ecosystem (scirpy, muon) becomes available
and a gene-expression matrix aligns naturally on ``cell_id``.

This is the opposite of the bulk-cohort path: a cohort of many *repertoires* (per-
sample clonotype tables) must NOT go in AnnData — ``obs=clonotype`` yields an
~1e9 × 100k almost-empty sparse ``X`` — use :func:`vdjtools.io.scan_cohort` (a hive-
partitioned Parquet dataset scanned as one LazyFrame) for that. Rule of thumb:
single-cell (``obs=cell``) → AnnData; bulk cohort (per-sample tables) → parquet.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

import polars as pl

if TYPE_CHECKING:  # pragma: no cover - typing only
    import anndata as ad


def to_anndata(paired: pl.DataFrame, X=None, *, index: str = "pair_id") -> "ad.AnnData":
    """Wrap :func:`vdjtools.sc.pair_chains` output as an :class:`anndata.AnnData`.

    ``obs`` is one row per receptor pair, indexed by ``pair_id`` (unique even when a
    cell yields two α/β pairs), with ``cell_id`` kept as a column so an expression
    matrix can be joined on it. With no ``X`` the result is a pure VDJ container (an
    ``n_obs × 0`` matrix); pass a cells × genes matrix aligned to ``obs`` to attach
    gene expression. For a formally multimodal object combine this with a GEX AnnData
    under ``mudata.MuData({"gex": gex, "airr": to_anndata(paired)})`` (scirpy-ready).

    Args:
        paired: Paired-receptor frame from :func:`vdjtools.sc.pair_chains`
            (``cell_id, pair_id, alpha_*, beta_*`` columns).
        X: Optional feature matrix with one row per ``obs`` (e.g. gene expression);
            defaults to an empty ``n_obs × 0`` sparse matrix.
        index: Column to use as the unique ``obs`` index (default ``"pair_id"``).

    Returns:
        An :class:`anndata.AnnData` whose ``obs`` holds the paired-chain annotation.

    Raises:
        ImportError: If ``anndata`` (the ``[sc]`` extra) is not installed.
        ValueError: If ``index`` is not a column of ``paired`` or is not unique.
    """
    try:
        import anndata as ad
        import pandas as pd  # anndata depends on pandas; avoids a pyarrow requirement
        import scipy.sparse as sp
    except ImportError as e:  # pragma: no cover - optional dep
        raise ImportError("to_anndata needs the '[sc]' extra: pip install anndata") from e

    if index not in paired.columns:
        raise ValueError(f"index column {index!r} not in paired frame: {paired.columns}")
    if paired[index].n_unique() != paired.height:
        raise ValueError(f"{index!r} is not unique; cannot index obs by it")

    # Build obs via a plain dict (not .to_pandas(), which needs pyarrow).
    obs = pd.DataFrame(paired.to_dict(as_series=False))
    obs.index = obs[index].astype(str)
    obs.index.name = index
    if X is None:
        X = sp.csr_matrix((obs.shape[0], 0), dtype="float32")
    return ad.AnnData(X=X, obs=obs)
