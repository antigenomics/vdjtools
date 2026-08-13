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


def column_names(df: "pl.DataFrame | pl.LazyFrame") -> list[str]:
    """Column names of an eager **or** lazy frame.

    Uses :meth:`polars.LazyFrame.collect_schema` for a ``LazyFrame`` so the check
    does not emit polars' "resolving schema" performance warning; falls back to
    ``.columns`` for an eager ``DataFrame``.

    Args:
        df: A ``pl.DataFrame`` or ``pl.LazyFrame``.

    Returns:
        The list of column names.
    """
    return df.collect_schema().names() if isinstance(df, pl.LazyFrame) else df.columns


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


def normalize(df: pl.DataFrame, *, recompute_freq: bool = False,
              keep: tuple[str, ...] = ()) -> pl.DataFrame:
    """Coerce an arbitrary frame to the canonical clonotype schema.

    Missing canonical columns are added as nulls, present ones are cast to their
    declared dtype (non-strict — unparseable values become null). The result is
    the canonical columns in canonical order, followed by any ``keep`` columns;
    every other non-canonical column (e.g. native vdjtools markup like
    ``VEnd``/``DStart``) is dropped.

    Args:
        df: A frame that already uses canonical column names for whatever columns
            it carries.
        recompute_freq: If ``True``, recompute ``frequency`` from ``duplicate_count``
            after coercion (use when the source lacks a trustworthy frequency).
        keep: Non-canonical columns to preserve, e.g. ``("v_identity",)``. Dtypes
            are left alone — the canonical schema has nothing to say about them.

    Returns:
        A frame with the canonical columns, correctly typed and ordered, plus ``keep``.
    """
    exprs = []
    for col, dtype in SCHEMA.items():
        if col not in df.columns:
            exprs.append(pl.lit(None, dtype=dtype).alias(col))
        elif col == COUNT:
            # Route the count through Float64: read via the all-Utf8 TSV path, a count of "5000.0"
            # (pandas float-formats any integer column that once held a NaN) fails a direct
            # Utf8->Int64 cast, becomes null, and recompute_frequency then reports frequency 0.0
            # for a real clone. Float64->Int64 truncates toward zero (== io/convert.py::_to_int).
            exprs.append(pl.col(col).cast(pl.Float64, strict=False).cast(dtype, strict=False).alias(col))
        else:
            exprs.append(pl.col(col).cast(dtype, strict=False).alias(col))
    df = df.with_columns(exprs)
    if recompute_freq:
        df = recompute_frequency(df)
    # Canonical columns first, then whatever the reader was asked to keep. Without the second
    # part a `keep=` upstream is silently undone here, which is worse than never offering one:
    # the caller gets a frame with no v_identity and no error to explain it.
    extra = [c for c in keep if c in df.columns and c not in COLUMNS]
    return df.select(list(COLUMNS) + extra)


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
    """Reduce a segment-call expression to gene resolution, ambiguity-safe.

    Strips the IMGT allele suffix from **every** gene an AIRR call names, not just the first. The
    old ``\\*.*$`` regex matched from the FIRST ``*`` to end of string, so a comma-ambiguous call
    like ``IGHV3-23*01,IGHV3-23D*01`` collapsed to ``IGHV3-23`` -- silently dropping IGHV3-23D,
    which then reported zero usage across a whole cohort despite being named in tens of thousands
    of rows. Genes are de-duplicated after stripping, so an allele-level tie *within* one gene
    (``IGHV1-2*02,IGHV1-2*04``) correctly collapses to the single unambiguous gene ``IGHV1-2``,
    while a genuine cross-gene tie stays ``IGHV3-23,IGHV3-23D``.

    Args:
        expr: A polars string expression over segment calls.

    Returns:
        Each call reduced to its distinct gene(s), sorted and comma-joined
        (``TRBV12-3*01`` → ``TRBV12-3``; ``A*01,A*02`` → ``A``; ``A*01,B*01`` → ``A,B``);
        nulls pass through unchanged.
    """
    return (expr.str.split(",")
            .list.eval(pl.element().str.strip_chars().str.replace(r"\*.*$", ""))
            .list.unique().list.sort().list.join(","))


def resolve_gene(expr: pl.Expr) -> pl.Expr:
    """Reduce a segment call to exactly ONE gene: allele stripped, ambiguity resolved to the first.

    The companion to :func:`strip_allele`, and the distinction matters:

    - :func:`strip_allele` keeps a genuine cross-gene tie as ``IGHV3-23,IGHV3-23D``, because when
      you are *reporting* usage you must not invent certainty the aligner did not have.
    - :func:`resolve_gene` collapses it to ``IGHV3-23``, because when the gene is a **feature
      axis** every distinct ambiguity string otherwise becomes its own category.

    That second failure is not hypothetical. Fitting a V+k-mer vocabulary on 200 HIP samples
    produced **1,296 V "genes", 1,235 of them comma-strings** such as
    ``TRBV1,TRBV23-1,TRBV4-1,TRBV4-2,TRBV4-3``. The cost is not the 21x wider axis: it is that the
    real ``TRBV9`` bucket gets *drained*, since every TRBV9 clone that happened to be called
    ambiguously was filed elsewhere. Its features then fall below any incidence floor and vanish,
    so a cohort with clean calls is scored against columns nobody populated.

    Where the ambiguity comes from matters, because it is not a parsing artifact and will not go
    away: that cohort is Adaptive/immunoSEQ **realigned with MiXCR against IMGT from the junction
    plus short flanks**. The realignment is the better call -- MiXCR/IMGT is a sounder reference
    than Adaptive's own -- and the ambiguity is what honest calling looks like when V genes differ
    only outside the sequenced window. So first-listed is a resolution, not a correction.

    Two things it cannot fix, and which belong to the assay rather than the reference: the window
    still bounds what is resolvable, and Adaptive's multiplex V primers distort V usage
    frequencies. A V-conditioned feature axis *fitted* on such a cohort inherits both, which makes
    it a poor donor for a 5'RACE cohort whatever this function does.

    First-listed rather than dropped: an ambiguous call still carries a clonotype, and the
    aligner lists its best call first. Note this takes the first call **as written**, not
    ``strip_allele(...).list.first()`` -- ``strip_allele`` sorts, so composing the two silently
    returns the alphabetically-first gene instead (``TRBV5*01,TRBV19*03`` -> ``TRBV19``, not
    ``TRBV5``).

    Args:
        expr: A polars string expression over segment calls.

    Returns:
        One gene per row (``TRBV12-3*01`` -> ``TRBV12-3``; ``A*01,B*01`` -> ``A``); nulls pass
        through unchanged.
    """
    return (expr.str.split(",").list.first()
            .str.strip_chars().str.replace(r"\*.*$", ""))
