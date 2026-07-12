"""CDR3-length spectratypes (long-format polars).

The spectratype is the weighted distribution of CDR3 lengths. ``vj_spectratype``
additionally breaks the distribution down by V-J pairing (legacy ``SpectratypeV``
is the V-only special case, obtained by ignoring the ``j_call`` column).
"""
from __future__ import annotations

import polars as pl

from ..io.schema import (
    JUNCTION_AA,
    JUNCTION_NT,
    J_CALL,
    LOCUS,
    V_CALL,
    add_locus,
    strip_allele,
    weight_expr,
)

_KIND_COL = {"aa": JUNCTION_AA, "nt": JUNCTION_NT}


def _length_frame(df: pl.DataFrame, kind: str, weight: str) -> pl.DataFrame:
    """Attach ``length`` (CDR3 length) and ``weight`` columns; drop null sequences."""
    if kind not in _KIND_COL:
        raise ValueError(f"kind must be 'aa' or 'nt'; got {kind!r}")
    col = _KIND_COL[kind]
    df = df if LOCUS in df.columns else add_locus(df)
    df = df.with_columns(
        pl.col(col).str.len_chars().cast(pl.Int64).alias("length"),
        weight_expr(weight).alias("weight"),
    )
    return df.filter(pl.col(col).is_not_null() & (pl.col("length") > 0))


def spectratype(df: pl.DataFrame, kind: str = "aa", weight: str = "reads",
                by_locus: bool = True) -> pl.DataFrame:
    """CDR3-length distribution.

    Args:
        df: A clonotype frame.
        kind: ``"aa"`` (length of ``junction_aa``) or ``"nt"`` (length of ``junction_nt``).
        weight: ``"reads"``, ``"unique"``, or ``"freq"``.
        by_locus: If ``True``, break the distribution down per locus.

    Returns:
        Long-format ``pl.DataFrame`` with ``length``, ``weight`` and (if ``by_locus``)
        ``locus``, sorted by length.
    """
    df = _length_frame(df, kind, weight)
    group = [LOCUS, "length"] if by_locus else ["length"]
    return (df.group_by(group, maintain_order=True)
              .agg(pl.col("weight").sum())
              .sort(group))


def vj_spectratype(df: pl.DataFrame, kind: str = "aa", weight: str = "reads",
                   by_locus: bool = True, keep_allele: bool = False) -> pl.DataFrame:
    """CDR3-length distribution broken down by V-J pairing.

    Args:
        df: A clonotype frame.
        kind: ``"aa"`` or ``"nt"`` (see :func:`spectratype`).
        weight: ``"reads"``, ``"unique"``, or ``"freq"``.
        by_locus: If ``True``, break the distribution down per locus.
        keep_allele: If ``True``, keep allele suffixes; otherwise collapse V/J to
            gene level (default).

    Returns:
        Long-format ``pl.DataFrame`` with ``v_call``, ``j_call``, ``length``,
        ``weight`` and (if ``by_locus``) ``locus``. A V-only spectratype is obtained
        by summing ``weight`` over ``j_call``.
    """
    df = _length_frame(df, kind, weight)
    v = pl.col(V_CALL) if keep_allele else strip_allele(pl.col(V_CALL))
    j = pl.col(J_CALL) if keep_allele else strip_allele(pl.col(J_CALL))
    df = df.with_columns(v.alias(V_CALL), j.alias(J_CALL))
    df = df.filter(pl.col(V_CALL).is_not_null() & pl.col(J_CALL).is_not_null())
    group = [LOCUS, V_CALL, J_CALL, "length"] if by_locus else [V_CALL, J_CALL, "length"]
    return (df.group_by(group, maintain_order=True)
              .agg(pl.col("weight").sum())
              .sort(group))
