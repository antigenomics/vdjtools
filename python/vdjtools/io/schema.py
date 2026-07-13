"""Canonical clonotype-frame schema (AIRR-aligned) and coercion helpers.

The basic-analytics layer speaks a single, flat clonotype frame — one row per
clonotype, AIRR Rearrangement column names, polars dtypes. Every reader emits it
and every analysis function consumes it. Kept deliberately minimal (free functions,
no classes) to mirror the vdjmatch / arda convention.

Columns:

* ``v_call, d_call, j_call, c_call`` (Utf8, nullable) — IMGT segment calls;
  ``c_call`` is frequently absent in native vdjtools data.
* ``junction_aa`` (Utf8) — the junction amino-acid sequence (conserved anchors
  Cys104 … Phe/Trp118 **INCLUDED**), per the AIRR ``junction_aa`` convention
  (equivalently the legacy vdjtools ``cdr3aa``). This is two residues longer
  than the IMGT ``cdr3_aa`` (anchors excluded); readers prefer the junction form.
* ``junction_nt`` (Utf8, nullable) — the junction nucleotide sequence (anchors
  included), matching ``junction_aa`` above. AIRR spells the nucleotide junction
  ``junction`` (no ``_nt`` suffix); readers accept that (and legacy ``cdr3_nt``)
  as input aliases.
* ``duplicate_count`` (Int64) — read/UMI count for the clonotype.
* ``frequency`` (Float64) — ``duplicate_count`` normalised within the sample.
* ``locus`` (Utf8, derived) — first three characters of ``v_call`` (``TRB``, ``IGH`` …).
"""
from __future__ import annotations

import polars as pl

V_CALL = "v_call"
D_CALL = "d_call"
J_CALL = "j_call"
C_CALL = "c_call"
JUNCTION_AA = "junction_aa"
JUNCTION_NT = "junction_nt"
COUNT = "duplicate_count"
FREQ = "frequency"
LOCUS = "locus"

#: Canonical columns in canonical order, mapped to their polars dtype.
SCHEMA: dict[str, pl.DataType] = {
    V_CALL: pl.Utf8,
    D_CALL: pl.Utf8,
    J_CALL: pl.Utf8,
    C_CALL: pl.Utf8,
    JUNCTION_AA: pl.Utf8,
    JUNCTION_NT: pl.Utf8,
    COUNT: pl.Int64,
    FREQ: pl.Float64,
}

#: Column names in canonical order.
COLUMNS: list[str] = list(SCHEMA)


def locus_of(v_call: str | None) -> str | None:
    """Return the locus (first three characters) of an IMGT V-gene call.

    Args:
        v_call: An IMGT V-gene call such as ``"TRBV12-3*01"``, or ``None``.

    Returns:
        The three-letter locus (``"TRB"``), or ``None`` if ``v_call`` is ``None``
        or shorter than three characters.

    Example:
        >>> locus_of("TRBV12-3*01")
        'TRB'
    """
    if v_call is None or len(v_call) < 3:
        return None
    return v_call[:3]


def add_locus(df: pl.DataFrame) -> pl.DataFrame:
    """Add (or overwrite) the derived ``locus`` column from ``v_call``.

    Args:
        df: A clonotype frame carrying a ``v_call`` column.

    Returns:
        The frame with a ``locus`` column (null where ``v_call`` is null).
    """
    return df.with_columns(pl.col(V_CALL).str.slice(0, 3).alias(LOCUS))


def recompute_frequency(df: pl.DataFrame) -> pl.DataFrame:
    """Recompute ``frequency`` as ``duplicate_count / sum(duplicate_count)``.

    Args:
        df: A clonotype frame with a ``duplicate_count`` column.

    Returns:
        The frame with ``frequency`` overwritten. If the total count is zero the
        frequency is set to ``0.0`` for every row.
    """
    total = df[COUNT].sum()
    if not total:
        return df.with_columns(pl.lit(0.0, dtype=pl.Float64).alias(FREQ))
    return df.with_columns((pl.col(COUNT) / pl.lit(total)).cast(pl.Float64).alias(FREQ))


def normalize(df: pl.DataFrame, *, recompute_freq: bool = False) -> pl.DataFrame:
    """Coerce an arbitrary frame to the canonical clonotype schema.

    Missing canonical columns are added as nulls, present ones are cast to their
    declared dtype (non-strict — unparseable values become null). The result is
    exactly the canonical columns in canonical order; any non-canonical columns
    (e.g. native vdjtools markup like ``VEnd``/``DStart``) are dropped.

    Args:
        df: A frame that already uses canonical column names for whatever columns
            it carries.
        recompute_freq: If ``True``, recompute ``frequency`` from ``duplicate_count``
            after coercion (use when the source lacks a trustworthy frequency).

    Returns:
        A frame with exactly the canonical columns, correctly typed and ordered.
    """
    exprs = []
    for col, dtype in SCHEMA.items():
        if col in df.columns:
            exprs.append(pl.col(col).cast(dtype, strict=False).alias(col))
        else:
            exprs.append(pl.lit(None, dtype=dtype).alias(col))
    df = df.with_columns(exprs)
    if recompute_freq:
        df = recompute_frequency(df)
    return df.select(COLUMNS)


def weight_expr(weight: str) -> pl.Expr:
    """Return the per-clonotype weight expression for an analysis mode.

    Args:
        weight: One of ``"reads"`` (weight by ``duplicate_count``), ``"unique"``
            (one per clonotype), or ``"freq"`` / ``"frequency"`` (weight by
            ``frequency``).

    Returns:
        A polars expression yielding the per-row weight.

    Raises:
        ValueError: If ``weight`` is not a recognised mode.
    """
    if weight == "reads":
        return pl.col(COUNT)
    if weight == "unique":
        return pl.lit(1, dtype=pl.Int64)
    if weight in ("freq", "frequency"):
        return pl.col(FREQ)
    raise ValueError(f"weight must be 'reads', 'unique' or 'freq'; got {weight!r}")


def strip_allele(expr: pl.Expr) -> pl.Expr:
    """Strip the IMGT allele suffix (``*01``) from a segment-call expression.

    Args:
        expr: A polars string expression over segment calls.

    Returns:
        The expression with everything from the first ``*`` onward removed
        (``TRBV12-3*01`` → ``TRBV12-3``); nulls pass through unchanged.
    """
    return expr.str.replace(r"\*.*$", "")
