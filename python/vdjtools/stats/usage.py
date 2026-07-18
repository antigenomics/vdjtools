"""V / D / J / C segment usage and V-J pairing profiles (long-format polars).

Usage is the summed weight per segment: reads (``duplicate_count``) or unique
clonotypes. Allele suffixes (``*01``) are stripped to gene level by default.
Results are long-format and normalisable to fractions by dividing within a locus.
"""
from __future__ import annotations

import polars as pl

from ..io.schema import (
    C_CALL,
    D_CALL,
    J_CALL,
    LOCUS,
    V_CALL,
    add_locus,
    column_names,
    strip_allele,
    weight_expr,
)

_SEGMENT_COL = {"v": V_CALL, "d": D_CALL, "j": J_CALL, "c": C_CALL}


def _ensure_locus(df):
    return df if LOCUS in column_names(df) else add_locus(df)


def segment_usage(df, segment: str, weight: str = "reads",
                  by_locus: bool = True, keep_allele: bool = False, by=()):
    """Segment usage profile for one of the V, D, J, or C genes.

    Args:
        df: A clonotype frame (eager ``pl.DataFrame`` or a lazy ``pl.LazyFrame`` — e.g.
            a whole cohort from :func:`vdjtools.io.scan_cohort`); the result mirrors the
            input's laziness.
        segment: One of ``"v"``, ``"d"``, ``"j"``, ``"c"``.
        weight: ``"reads"`` (sum ``duplicate_count``), ``"unique"`` (count
            clonotypes), or ``"freq"`` (sum ``frequency``).
        by_locus: If ``True``, break the profile down per locus.
        keep_allele: If ``True``, keep the IMGT allele suffix; otherwise collapse to
            gene level (default).
        by: Extra column(s) to **prepend** to the group key (e.g. ``["sample_id"]`` to
            compute the whole cohort in one grouped pass over a ``LazyFrame``). Empty
            by default — output is byte-identical to the per-sample profile.

    Returns:
        Long-format frame with the segment call column (``v_call`` …), a ``weight``
        column, the ``by`` columns, and (if ``by_locus``) a ``locus`` column; a
        ``pl.LazyFrame`` when ``df`` is lazy. Rows with a null segment call are
        dropped, so an all-null C gene yields an empty frame.

    Raises:
        ValueError: If ``segment`` is not one of v/d/j/c.
    """
    if segment not in _SEGMENT_COL:
        raise ValueError(f"segment must be one of {list(_SEGMENT_COL)}; got {segment!r}")
    col = _SEGMENT_COL[segment]
    df = _ensure_locus(df)
    gene = pl.col(col) if keep_allele else strip_allele(pl.col(col))
    df = df.with_columns(gene.alias(col), weight_expr(weight).alias("weight"))
    df = df.filter(pl.col(col).is_not_null())
    group = [*by, *([LOCUS, col] if by_locus else [col])]
    return (df.group_by(group, maintain_order=True)
              .agg(pl.col("weight").sum())
              .sort(group))


def vj_usage(df, weight: str = "reads", by_locus: bool = True,
             keep_allele: bool = False, by=()):
    """V-J pairing usage profile.

    Args:
        df: A clonotype frame (eager or lazy; see :func:`segment_usage`).
        weight: ``"reads"``, ``"unique"``, or ``"freq"`` (see :func:`segment_usage`).
        by_locus: If ``True``, break the profile down per locus.
        keep_allele: If ``True``, keep allele suffixes; otherwise collapse to gene
            level (default).
        by: Extra column(s) prepended to the group key (e.g. ``["sample_id"]`` for a
            one-pass cohort profile).

    Returns:
        Long-format frame with ``v_call``, ``j_call``, ``weight``, the ``by`` columns
        and (if ``by_locus``) ``locus`` (lazy when ``df`` is lazy). Rows with a null V
        or J call are dropped.
    """
    df = _ensure_locus(df)
    v = pl.col(V_CALL) if keep_allele else strip_allele(pl.col(V_CALL))
    j = pl.col(J_CALL) if keep_allele else strip_allele(pl.col(J_CALL))
    df = df.with_columns(v.alias(V_CALL), j.alias(J_CALL),
                         weight_expr(weight).alias("weight"))
    df = df.filter(pl.col(V_CALL).is_not_null() & pl.col(J_CALL).is_not_null())
    group = [*by, *([LOCUS, V_CALL, J_CALL] if by_locus else [V_CALL, J_CALL])]
    return (df.group_by(group, maintain_order=True)
              .agg(pl.col("weight").sum())
              .sort(group))
