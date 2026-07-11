"""Format auto-detection and metadata-driven batch reading.

Sits on top of the single-file readers in :mod:`vdjtools.io.read`: sniffs a file's
format, dispatches to the right reader, and loads a whole cohort described by a
metadata table into one long clonotype frame (or a per-sample dict).
"""
from __future__ import annotations

import os
from pathlib import Path

import polars as pl

from .read import read_airr, read_vdjtools
from .schema import COLUMNS, LOCUS


def _header_columns(path: str | os.PathLike) -> list[str]:
    """Return the (lower-cased) header column names of a TSV without reading rows."""
    cols = pl.read_csv(Path(path), separator="\t", n_rows=0, infer_schema_length=0).columns
    return [c.lower() for c in cols]


def sniff_format(path: str | os.PathLike) -> str:
    """Detect a clonotype file's format from its header.

    Args:
        path: Path to a clonotype table.

    Returns:
        ``"vdjtools"`` if the header looks like the native vdjtools table
        (``cdr3aa`` / ``count`` + ``cdr3nt``), ``"airr"`` if it looks like AIRR
        Rearrangement (``v_call`` / ``junction_aa`` / ``cdr3_aa``).

    Raises:
        ValueError: If neither signature is recognised.
    """
    cols = set(_header_columns(path))
    if "cdr3aa" in cols or {"count", "cdr3nt"} <= cols:
        return "vdjtools"
    if cols & {"v_call", "junction_aa", "cdr3_aa"}:
        return "airr"
    raise ValueError(f"unrecognised clonotype format; header columns: {sorted(cols)}")


def read(path: str | os.PathLike, fmt: str = "auto") -> pl.DataFrame:
    """Read a clonotype table, auto-detecting the format by default.

    Args:
        path: Path to a native vdjtools or AIRR Rearrangement table (``.gz`` ok).
        fmt: ``"auto"`` (sniff the header), ``"vdjtools"``, or ``"airr"``.

    Returns:
        Canonical clonotype frame.

    Raises:
        ValueError: If ``fmt`` is unknown or auto-detection fails.
    """
    if fmt == "auto":
        fmt = sniff_format(path)
    if fmt == "vdjtools":
        return read_vdjtools(path)
    if fmt == "airr":
        return read_airr(path)
    raise ValueError(f"fmt must be 'auto', 'vdjtools' or 'airr'; got {fmt!r}")


def read_metadata(path: str | os.PathLike) -> pl.DataFrame:
    """Read a sample metadata TSV.

    The literal string ``"nan"`` (and empty strings) are treated as null. All columns
    are read as strings so metadata joins are stable.

    Args:
        path: Path to a metadata TSV (one row per sample).

    Returns:
        A ``pl.DataFrame`` of metadata, all-Utf8, with ``"nan"`` → null.
    """
    return pl.read_csv(Path(path), separator="\t", infer_schema_length=0,
                       null_values=["nan", "NaN", ""])


def read_samples(metadata: pl.DataFrame, base_dir: str | os.PathLike,
                 sample_col: str = "sample_name", file_template: str = "{sample}.tsv.gz",
                 fmt: str = "auto", add_metadata: bool = True, as_dict: bool = False):
    """Read a batch of samples described by a metadata table.

    For each metadata row the file ``base_dir/file_template.format(sample=<sample>)``
    is read into the canonical schema, tagged with the reserved columns ``sample_id``
    (the ``sample_col`` value) and ``file_name``, and — if ``add_metadata`` — joined
    with that row's metadata columns.

    Args:
        metadata: Metadata frame (e.g. from :func:`read_metadata`).
        base_dir: Directory holding the per-sample clonotype files.
        sample_col: Metadata column carrying the sample name (default
            ``"sample_name"``).
        file_template: Filename template; ``{sample}`` is substituted with the sample
            name (default ``"{sample}.tsv.gz"``).
        fmt: Format for the readers (``"auto"`` / ``"vdjtools"`` / ``"airr"``).
        add_metadata: If ``True``, attach the metadata columns to every clonotype row.
        as_dict: If ``True``, return ``{sample_id: frame}`` instead of one long frame.

    Returns:
        One concatenated long ``pl.DataFrame`` (canonical columns + ``locus`` +
        ``sample_id`` + ``file_name`` + metadata), or a ``dict[str, pl.DataFrame]``
        if ``as_dict``.
    """
    base = Path(base_dir)
    keep = [*COLUMNS, LOCUS]
    reserved = {"sample_id", "file_name"}
    meta_cols = [c for c in metadata.columns if c not in reserved]
    frames: dict[str, pl.DataFrame] = {}
    for row in metadata.iter_rows(named=True):
        sample = row[sample_col]
        fname = file_template.format(sample=sample)
        clones = read(base / fname, fmt=fmt)
        clones = clones.select([c for c in keep if c in clones.columns])
        clones = clones.with_columns(pl.lit(sample).alias("sample_id"),
                                     pl.lit(fname).alias("file_name"))
        if add_metadata:
            clones = clones.with_columns(
                pl.lit(row[mc]).alias(mc) for mc in meta_cols if mc != sample_col
            )
        frames[str(sample)] = clones
    if as_dict:
        return frames
    return pl.concat(list(frames.values()), how="vertical_relaxed")
