"""Bridge single-cell receptors into an AnnData / MuData container.

Single-cell is the shape where AnnData fits: ``obs`` is one row per observation, so the
whole scverse ecosystem (scirpy, muon) becomes available and a gene-expression matrix
aligns naturally on ``cell_id``.

Two containers, for two different jobs:

* :func:`to_scirpy` — the **scverse-native** layout. Chains live as an awkward array under
  ``adata.obsm["airr"]`` (cell → variable-length chain list → AIRR record), which is what
  scirpy ≥0.13 reads. Use this to hand data to scirpy, and :func:`from_scirpy` to get it
  back. Nothing is lost: it carries the full per-contig AIRR record.
* :func:`to_anndata` — a **flat** container, one ``obs`` row per receptor *pair*
  (``alpha_*``/``beta_*`` columns). Convenient for a quick paired-chain table or for
  attaching an expression matrix, but it is not what scirpy consumes and it keeps only the
  paired-chain fields.

Writing delegates to scirpy (``ir.io.read_airr``) because scirpy's reader is the source of
truth for its own on-disk layout and reimplementing it here would drift with their schema.
Reading is ours and needs only ``awkward``, so consuming a scirpy object costs no scirpy
install.

This is the opposite of the bulk-cohort path: a cohort of many *repertoires* (per-
sample clonotype tables) must NOT go in AnnData — ``obs=clonotype`` yields an
~1e9 × 100k almost-empty sparse ``X`` — use :func:`vdjtools.io.scan_cohort` (a hive-
partitioned Parquet dataset scanned as one LazyFrame) for that. Rule of thumb:
single-cell (``obs=cell``) → AnnData; bulk cohort (per-sample tables) → parquet.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

import polars as pl

from .airr import to_airr
from .read import CELL_ID

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


def _to_pandas(df: pl.DataFrame):
    """polars -> pandas without requiring pyarrow (see :func:`to_anndata`)."""
    import pandas as pd

    return pd.DataFrame(df.to_dict(as_series=False))


def _from_pandas(df) -> pl.DataFrame:
    """pandas -> polars without requiring pyarrow.

    pandas spells a missing value in a text column as ``NaN``, which polars rejects when the
    rest of the column is strings ("'float' object is not an instance of 'str'"), so those
    are normalised to ``None`` first. Keyed on the *numpy* kind rather than the pandas
    dtype: pandas 2 gives such a column ``object`` and pandas 3 gives it ``str``, but both
    land on an object array here.
    """
    import pandas as pd

    cols = {}
    for name in df.columns:
        arr = df[name].to_numpy()
        if arr.dtype.kind == "O":
            arr = [None if v is None or (isinstance(v, float) and pd.isna(v)) else v
                   for v in arr]
        cols[name] = arr
    return pl.DataFrame(cols)


def to_scirpy(cells: pl.DataFrame, gex=None, *, index_chains: bool = True):
    """Convert the sc long frame to a scirpy-native AnnData (or MuData with ``gex``).

    Delegates to ``scirpy.io.read_airr`` on the AIRR table from :func:`~vdjtools.sc.to_airr`,
    so the result carries scirpy's own ``obsm["airr"]`` awkward layout exactly as their
    version defines it. The whole per-contig AIRR record survives, unlike the flat
    :func:`to_anndata` view.

    Args:
        cells: Single-cell long frame (:data:`vdjtools.sc.read.SC_COLUMNS`).
        gex: Optional gene-expression :class:`anndata.AnnData`. When given, the result is a
            ``mudata.MuData({"gex": gex, "airr": ...})`` — the multimodal object scirpy
            reads with ``airr_mod="airr"``.
        index_chains: Run ``scirpy.pp.index_chains`` on the result (default ``True``), which
            populates ``obsm["chain_indices"]``. Every scirpy tool needs it, so doing it here
            is what makes the handoff seamless; pass ``False`` to index yourself.

    Returns:
        An :class:`anndata.AnnData`, or a ``mudata.MuData`` when ``gex`` is given.

    Raises:
        ImportError: If ``scirpy`` (or ``mudata``, when ``gex`` is given) is not installed.
    """
    try:
        import scirpy as ir
    except ImportError as e:  # pragma: no cover - optional dep
        raise ImportError(
            "to_scirpy needs scirpy: pip install scirpy. (from_scirpy needs only awkward.)"
        ) from e

    adata = ir.io.read_airr(_to_pandas(to_airr(cells)))
    if index_chains:
        ir.pp.index_chains(adata)
    if gex is None:
        return adata
    try:
        import mudata
    except ImportError as e:  # pragma: no cover - optional dep
        raise ImportError("passing `gex` needs mudata: pip install mudata") from e
    return mudata.MuData({"gex": gex, "airr": adata})


def from_scirpy(adata) -> pl.DataFrame:
    """Convert a scirpy AnnData / MuData back to the canonical sc long frame.

    Reads ``obsm["airr"]`` with ``awkward`` directly — **scirpy itself is not needed**, so a
    scirpy object handed to you is readable from a plain ``vdjtools[sc]`` install. The
    ``obs`` index is authoritative for ``cell_id`` (that is what AnnData keys on), so a
    stale ``cell_id`` field inside the chain records is ignored.

    Args:
        adata: An :class:`anndata.AnnData` with an ``obsm["airr"]`` awkward array, or a
            ``mudata.MuData`` carrying one in its ``"airr"`` modality.

    Returns:
        A ``pl.DataFrame`` in the canonical layout (:data:`vdjtools.sc.read.SC_COLUMNS`),
        one row per contig.

    Raises:
        ImportError: If ``awkward`` (the ``[sc]`` extra) is not installed.
        KeyError: If the object has no ``obsm["airr"]``.
    """
    try:
        import awkward as ak
    except ImportError as e:  # pragma: no cover - optional dep
        raise ImportError("from_scirpy needs the '[sc]' extra: pip install awkward") from e
    import numpy as np

    from .airr import from_airr

    if hasattr(adata, "mod"):           # MuData -> its AIRR modality
        adata = adata.mod["airr"]
    if "airr" not in adata.obsm:
        raise KeyError(
            "no obsm['airr'] on this object; scirpy >=0.13 stores chains there "
            "(older objects need scirpy.io.upgrade_schema first)"
        )

    airr = adata.obsm["airr"]
    flat = ak.flatten(airr, axis=1)
    # obs_names, repeated once per chain, is the cell key -- not any cell_id in the record.
    n_chains = np.asarray(ak.num(airr, axis=1))
    cell_ids = np.repeat(np.asarray(adata.obs_names, dtype=object), n_chains)

    data = {f: ak.to_list(flat[f]) for f in flat.fields if f != CELL_ID}
    data[CELL_ID] = list(cell_ids)
    return from_airr(pl.DataFrame(data))


def push_obs(target, df: pl.DataFrame, columns=None, *, on: str = CELL_ID):
    """Push vdjtools-computed per-cell columns into a downstream container, in place.

    The "augment" direction: take something vdjtools computed (``pgen_paired``,
    ``mispairing_flag``, a clustering) and attach it to an object the rest of an analysis is
    already built around. Works on anything with an ``obs`` (AnnData, MuData) or a
    ``metadata`` (dandelion ``Dandelion``) table; rows are matched by index, and cells the
    frame does not mention get nulls.

    Args:
        target: An :class:`anndata.AnnData` / ``mudata.MuData`` / ``Dandelion``.
        df: Per-cell frame carrying ``on`` plus the columns to attach.
        columns: Which columns to push (default: everything except ``on``).
        on: Key column in ``df`` matched against the target's index (default ``cell_id``).

    Returns:
        ``target``, mutated in place (returned for chaining).

    Raises:
        ValueError: If ``on`` is missing from ``df``, ``on`` is not unique (one row per cell
            is required -- aggregate a multi-pair frame first), or a requested column is
            absent.
        TypeError: If ``target`` has neither ``obs`` nor ``metadata``.
    """
    if on not in df.columns:
        raise ValueError(f"push_obs: {on!r} not in frame: {df.columns}")
    if df[on].n_unique() != df.height:
        raise ValueError(
            f"push_obs: {on!r} is not unique ({df.height} rows, {df[on].n_unique()} keys); "
            "aggregate to one row per cell first (a cell with two pairs yields two rows)"
        )
    cols = list(columns) if columns is not None else [c for c in df.columns if c != on]
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"push_obs: columns not in frame: {missing}")

    if hasattr(target, "obs"):
        frame = target.obs
    elif hasattr(target, "metadata"):
        frame = target.metadata
    else:
        raise TypeError(f"push_obs: {type(target).__name__} has neither .obs nor .metadata")

    pdf = _to_pandas(df.select([on, *cols]))
    pdf[on] = pdf[on].astype(str)
    aligned = pdf.set_index(on).reindex(frame.index.astype(str))
    for c in cols:
        frame[c] = aligned[c].to_numpy()
    return target
