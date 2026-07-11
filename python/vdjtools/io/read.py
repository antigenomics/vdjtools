"""Readers for the native vdjtools table format and AIRR Rearrangement TSV.

Both readers return the canonical clonotype frame described in :mod:`vdjtools.io.schema`.
Gzip-compressed inputs (``.gz``) are handled transparently by polars.
"""
from __future__ import annotations

import os
from pathlib import Path

import polars as pl

from . import schema
from .schema import (
    C_CALL,
    CDR3_AA,
    CDR3_NT,
    COUNT,
    D_CALL,
    J_CALL,
    V_CALL,
)

# Native vdjtools output header (com.antigenomics.vdjtools.misc.Software.VDJtools):
#   count  freq  cdr3nt  cdr3aa  v  d  j  VEnd  DStart  DEnd  JStart  [annotation...]
# Only the first seven fields carry clonotype content; VEnd..JStart are optional
# integer markup columns we do not carry into the canonical frame.
_NATIVE_MAP = {
    "count": COUNT,
    "cdr3nt": CDR3_NT,
    "cdr3aa": CDR3_AA,
    "v": V_CALL,
    "d": D_CALL,
    "j": J_CALL,
}

# AIRR Rearrangement source header -> canonical column, most specific first.
_AIRR_ALIASES: dict[str, tuple[str, ...]] = {
    V_CALL: ("v_call",),
    D_CALL: ("d_call",),
    J_CALL: ("j_call",),
    C_CALL: ("c_call",),
    # Prefer the junction (conserved anchors INCLUDED) — the canonical vdjtools /
    # AIRR ``junction_aa`` convention; fall back to the IMGT ``cdr3_aa``/``cdr3``
    # (anchors excluded) only if no junction column is present.
    CDR3_AA: ("junction_aa", "cdr3_aa"),
    CDR3_NT: ("junction", "cdr3"),
    # AIRR-hybrid exports (e.g. isalgo/airr_ankspond) carry a vdjtools-style `count`.
    COUNT: ("duplicate_count", "count", "reads"),
    # frequency is always recomputed from counts, never read from source.
}


def _read_tsv(path: str | os.PathLike, n_rows: int | None = None) -> pl.DataFrame:
    """Read a (optionally gzipped) TSV as all-Utf8; ``n_rows`` caps huge files."""
    return pl.read_csv(Path(path), separator="\t", infer_schema_length=0,
                       quote_char=None, null_values=["", "."], n_rows=n_rows)


def _first_call(expr: pl.Expr) -> pl.Expr:
    """Take the first (best) call from a comma/space-separated ambiguous call list."""
    return expr.str.split(",").list.first().str.strip_chars()


def read_vdjtools(path: str | os.PathLike, n_rows: int | None = None) -> pl.DataFrame:
    """Read a native vdjtools clonotype table into the canonical frame.

    The native header is ``count freq cdr3nt cdr3aa v d j VEnd DStart DEnd JStart``
    followed by optional annotation columns; ``.`` and empty strings mean missing.
    Ambiguous V/D/J calls (``"TRBV12-4, TRBV12-3"``) are reduced to the first
    (highest-confidence) call. ``c_call`` is always null (the format has no C gene).
    ``frequency`` is recomputed exactly from ``duplicate_count`` (the native ``freq``
    column is a rounded copy of the same ratio).

    Args:
        path: Path to a ``.txt`` or ``.txt.gz`` native vdjtools table.
        n_rows: If given, read at most this many data rows (preview huge files).

    Returns:
        Canonical clonotype frame with a derived ``locus`` column.

    Raises:
        ValueError: If the required ``count`` and ``cdr3aa`` columns are absent.
    """
    raw = _read_tsv(path, n_rows=n_rows)
    # Some legacy exports comment out the header line (``#count freq ...``); strip a
    # leading ``#`` so the first column still maps to ``count``.
    lower = {c.lower().lstrip("#"): c for c in raw.columns}
    found = {canon: lower[src] for src, canon in _NATIVE_MAP.items() if src in lower}
    if COUNT not in found or CDR3_AA not in found:
        raise ValueError(
            f"not a native vdjtools table (need 'count' and 'cdr3aa'); have {raw.columns}"
        )
    df = raw.rename({src: canon for canon, src in found.items()})
    df = df.with_columns(_first_call(pl.col(c)) for c in (V_CALL, D_CALL, J_CALL)
                         if c in df.columns)
    # The native ``freq`` column is just count/total; recompute it exactly from counts.
    df = schema.normalize(df, recompute_freq=True)
    return schema.add_locus(df)


def read_parquet(path: str | os.PathLike, n_rows: int | None = None) -> pl.DataFrame:
    """Read a Parquet clonotype table into the canonical frame.

    Parquet is the at-scale storage format for repertoire cohorts: typed, columnar,
    compressed, and — laid out as a hive-partitioned directory — scannable as one
    :class:`polars.LazyFrame` with predicate/projection pushdown (see
    :mod:`vdjtools.io.cohort`). Columns already using the canonical names (e.g. a
    file written by :meth:`polars.DataFrame.write_parquet` from a canonical frame)
    are kept as-is; AIRR source names (``junction_aa``, ``duplicate_count`` …) are
    mapped to canonical only to fill a canonical column that is otherwise absent.
    Unlike the TSV path this reads native dtypes directly (no all-Utf8 pass), so it
    is both faster and lower-peak-memory on large files.

    Args:
        path: Path to a ``.parquet`` / ``.pq`` clonotype table.
        n_rows: If given, read at most this many data rows (preview huge files).

    Returns:
        Canonical clonotype frame with a derived ``locus`` column.

    Raises:
        ValueError: If no CDR3 amino-acid column (``cdr3_aa`` / ``junction_aa``)
            is present.
    """
    df = pl.read_parquet(Path(path), n_rows=n_rows)
    have = set(df.columns)
    # Resolve each canonical column from its AIRR aliases with the SAME first-match
    # preference as read_airr (junction_aa before cdr3_aa, so the junction wins even
    # when a raw-AIRR parquet carries both). cdr3_nt is the one canonical column that
    # is not an AIRR alias, so fall back to an already-canonical column by name.
    found: dict[str, str] = {}
    for canon in schema.COLUMNS:
        src = next((s for s in _AIRR_ALIASES.get(canon, ()) if s in have), None)
        if src is None and canon in have:
            src = canon
        if src is not None:
            found[canon] = src
    if CDR3_AA not in found:
        raise ValueError(
            f"Parquet file lacks a CDR3 aa column (cdr3_aa/junction_aa); have {df.columns}"
        )
    df = df.select([pl.col(src).alias(canon) for canon, src in found.items()])
    df = df.with_columns(_first_call(pl.col(c)) for c in (V_CALL, D_CALL, J_CALL)
                         if c in df.columns)
    if COUNT not in df.columns:
        df = df.with_columns(pl.lit(1, dtype=pl.Int64).alias(COUNT))
    df = schema.normalize(df, recompute_freq=True)
    return schema.add_locus(df)


def read_airr(path: str | os.PathLike, *, collapse: bool = True,
              n_rows: int | None = None) -> pl.DataFrame:
    """Read an AIRR Rearrangement TSV into the canonical frame.

    Prefers the junction columns ``junction_aa`` / ``junction`` (conserved anchors
    INCLUDED — the canonical vdjtools / AIRR convention) over the IMGT ``cdr3_aa`` /
    ``cdr3`` columns (anchors excluded), which are used only as a fallback. Known
    limitation: a file that provides *only* the IMGT ``cdr3_aa`` yields sequences two
    residues shorter than the junction, so downstream lengths / k-mers will differ.
    The count column may be ``duplicate_count`` or a vdjtools-style ``count`` /
    ``reads`` (AIRR-hybrid exports); absent, it defaults to 1. Per-read files (one
    row per rearrangement) collapse to unique clonotypes with summed counts;
    already-aggregated files pass through unchanged. ``frequency`` is always
    recomputed after collapsing.

    Clonotype identity for collapsing is ``(v_call, j_call, cdr3_nt, cdr3_aa)`` —
    matching legacy ``Clonotype`` equality (V, J, CDR3nt); ``d_call`` and ``c_call``
    are *not* part of the identity, so a representative (first non-null) value is
    attached to each collapsed clonotype.

    Args:
        path: Path to a ``.tsv`` / ``.tsv.gz`` AIRR Rearrangement file.
        collapse: If ``True`` (default), sum the count over clonotypes identical on
            ``(v_call, j_call, cdr3_nt, cdr3_aa)``.
        n_rows: If given, read at most this many data rows (preview huge files).

    Returns:
        Canonical clonotype frame with a derived ``locus`` column.

    Raises:
        ValueError: If no CDR3 amino-acid column (``junction_aa`` or ``cdr3_aa``)
            is present.
    """
    raw = _read_tsv(path, n_rows=n_rows)
    lower = {c.lower(): c for c in raw.columns}
    found: dict[str, str] = {}
    for canon, srcs in _AIRR_ALIASES.items():
        for s in srcs:
            if s in lower:
                found[canon] = lower[s]
                break
    if CDR3_AA not in found:
        raise ValueError(
            f"AIRR file lacks a CDR3 aa column (cdr3_aa/junction_aa); have {raw.columns}"
        )
    df = raw.select([pl.col(src).alias(canon) for canon, src in found.items()])
    if COUNT not in df.columns:
        df = df.with_columns(pl.lit(1, dtype=pl.Int64).alias(COUNT))
    else:
        df = df.with_columns(pl.col(COUNT).cast(pl.Int64, strict=False).fill_null(1))

    # Legacy clonotype identity is (V, J, CDR3nt); cdr3_aa is kept in the key so
    # files with no nt column still collapse (redundant when nt is present). D and C
    # are not identity — carry a representative (first non-null) value per clonotype.
    key = [c for c in (V_CALL, J_CALL, CDR3_NT, CDR3_AA) if c in df.columns]
    if collapse:
        reps = [pl.col(c).drop_nulls().first().alias(c)
                for c in (D_CALL, C_CALL) if c in df.columns]
        df = df.group_by(key, maintain_order=True).agg(pl.col(COUNT).sum(), *reps)
    df = schema.normalize(df, recompute_freq=True)
    return schema.add_locus(df)
