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
    FREQ,
    J_CALL,
    V_CALL,
)

# Native vdjtools output header (com.antigenomics.vdjtools.misc.Software.VDJtools):
#   count  freq  cdr3nt  cdr3aa  v  d  j  VEnd  DStart  DEnd  JStart  [annotation...]
# Only the first seven fields carry clonotype content; VEnd..JStart are optional
# integer markup columns we do not carry into the canonical frame.
_NATIVE_MAP = {
    "count": COUNT,
    "freq": FREQ,
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
    # Prefer IMGT CDR3 (anchors excluded); fall back to junction (anchors included).
    CDR3_AA: ("cdr3_aa", "junction_aa"),
    CDR3_NT: ("cdr3", "junction"),
    # AIRR-hybrid exports (e.g. isalgo/airr_ankspond) carry a vdjtools-style `count`.
    COUNT: ("duplicate_count", "count", "reads"),
    FREQ: ("frequency", "freq"),
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
    lower = {c.lower(): c for c in raw.columns}
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


def read_airr(path: str | os.PathLike, *, collapse: bool = True,
              n_rows: int | None = None) -> pl.DataFrame:
    """Read an AIRR Rearrangement TSV into the canonical frame.

    Prefers the IMGT ``cdr3_aa`` / ``cdr3`` columns (anchors excluded) and falls
    back to ``junction_aa`` / ``junction`` (anchors included) when they are absent.
    The count column may be ``duplicate_count`` or a vdjtools-style ``count`` /
    ``reads`` (AIRR-hybrid exports); absent, it defaults to 1. Per-read files (one
    row per rearrangement) collapse to unique clonotypes with summed counts;
    already-aggregated files pass through unchanged. ``frequency`` is always
    recomputed after collapsing.

    Args:
        path: Path to a ``.tsv`` / ``.tsv.gz`` AIRR Rearrangement file.
        collapse: If ``True`` (default), sum the count over identical
            ``(v_call, d_call, j_call, c_call, cdr3_aa, cdr3_nt)`` clonotypes.
        n_rows: If given, read at most this many data rows (preview huge files).

    Returns:
        Canonical clonotype frame with a derived ``locus`` column.

    Raises:
        ValueError: If no CDR3 amino-acid column (``cdr3_aa`` or ``junction_aa``)
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

    key = [c for c in (V_CALL, D_CALL, J_CALL, C_CALL, CDR3_AA, CDR3_NT) if c in df.columns]
    if collapse:
        df = df.group_by(key, maintain_order=True).agg(pl.col(COUNT).sum())
    df = schema.normalize(df, recompute_freq=True)
    return schema.add_locus(df)
