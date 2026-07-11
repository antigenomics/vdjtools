"""CDR3 amino-acid k-mer profiles and joint V + k-mer + C feature summaries.

K-mers are overlapping sliding windows over ``cdr3_aa``. Each k-mer occurrence
carries its clonotype's weight (reads, unique, or frequency); weights are summed.
"""
from __future__ import annotations

import polars as pl

from ..io.schema import (
    C_CALL,
    CDR3_AA,
    LOCUS,
    V_CALL,
    add_locus,
    strip_allele,
    weight_expr,
)


def _explode_kmers(df: pl.DataFrame, k: int) -> pl.DataFrame:
    """Explode each clonotype's ``cdr3_aa`` into overlapping k-mers (column ``kmer``).

    Rows whose CDR3 is shorter than ``k`` (or null) contribute no k-mers.
    """
    if k < 1:
        raise ValueError(f"k must be >= 1; got {k}")
    df = df.with_columns(pl.col(CDR3_AA).str.len_chars().alias("_len"))
    df = df.filter(pl.col("_len") >= k)  # shorter/null CDR3s yield no k-mers
    df = df.with_columns(pl.int_ranges(0, pl.col("_len") - k + 1).alias("_pos"))
    df = df.explode("_pos", empty_as_null=True)
    df = df.with_columns(pl.col(CDR3_AA).str.slice(pl.col("_pos"), k).alias("kmer"))
    return df.drop("_len", "_pos")


def kmer_profile(df: pl.DataFrame, k: int = 3, weight: str = "reads",
                 by_locus: bool = True) -> pl.DataFrame:
    """CDR3 amino-acid k-mer spectrum.

    Args:
        df: A clonotype frame.
        k: K-mer length (default 3).
        weight: ``"reads"``, ``"unique"``, or ``"freq"``.
        by_locus: If ``True``, break the spectrum down per locus.

    Returns:
        Long-format ``pl.DataFrame`` with ``kmer``, ``weight`` and (if ``by_locus``)
        ``locus``, sorted by k-mer.
    """
    df = df if LOCUS in df.columns else add_locus(df)
    df = df.with_columns(weight_expr(weight).alias("weight"))
    df = _explode_kmers(df, k)
    group = [LOCUS, "kmer"] if by_locus else ["kmer"]
    return (df.group_by(group, maintain_order=True)
              .agg(pl.col("weight").sum())
              .sort(group))


def v_kmer_c_profile(df: pl.DataFrame, k: int = 3, weight: str = "reads",
                     by_locus: bool = True, keep_allele: bool = False) -> pl.DataFrame:
    """Joint (V gene, k-mer, C gene) profile — a tidy feature-matrix source.

    Produces one aggregated row per ``(v_call, kmer, c_call)`` combination (plus
    ``locus`` if requested), suitable to pivot into a feature matrix. A null
    ``c_call`` (common in native vdjtools data) is retained as its own group.

    Args:
        df: A clonotype frame.
        k: K-mer length (default 3).
        weight: ``"reads"``, ``"unique"``, or ``"freq"``.
        by_locus: If ``True``, include ``locus`` in the grouping.
        keep_allele: If ``True``, keep V allele suffixes; otherwise collapse V to
            gene level (default).

    Returns:
        Long-format ``pl.DataFrame`` with ``v_call``, ``kmer``, ``c_call``,
        ``weight`` and (if ``by_locus``) ``locus``.
    """
    df = df if LOCUS in df.columns else add_locus(df)
    v = pl.col(V_CALL) if keep_allele else strip_allele(pl.col(V_CALL))
    df = df.with_columns(v.alias(V_CALL), weight_expr(weight).alias("weight"))
    df = _explode_kmers(df, k)
    group = [LOCUS, V_CALL, "kmer", C_CALL] if by_locus else [V_CALL, "kmer", C_CALL]
    return (df.group_by(group, maintain_order=True)
              .agg(pl.col("weight").sum())
              .sort(group, nulls_last=True))
