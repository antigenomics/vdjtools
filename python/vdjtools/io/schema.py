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

#: OPTIONAL AIRR Rearrangement annotation columns. Not in SCHEMA -- they are not required and most
#: bulk pipelines do not emit them -- but when a file DOES carry them they are authoritative and
#: :func:`vdjtools.preprocess.filter_productive` reads them in preference to re-deriving the same
#: fact from ``junction_aa``. See https://docs.airr-community.org/en/latest/datarep/rearrangements.html
#:
#: ``productive`` is the composite: an open reading frame, no defect in the start codon, splicing
#: sites or regulatory elements, no internal stop codon, and an in-frame junction. The other two
#: are components of it, useful when the composite is absent.
PRODUCTIVE = "productive"
STOP_CODON = "stop_codon"
VJ_IN_FRAME = "vj_in_frame"

#: The optional AIRR annotation columns, in the order they are preferred as evidence.
AIRR_FUNCTIONAL_COLUMNS: tuple[str, ...] = (PRODUCTIVE, STOP_CODON, VJ_IN_FRAME)

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


#: The clonotype identity key when there is no ``junction_nt`` to distinguish rows.
#: ``c_call`` is included because an isotype-resolved IGH table legitimately carries the same
#: junction under IGHM and IGHG; ``d_call`` is not, because it is an inference from the junction
#: rather than an independent observation, and aligners disagree on it for identical input.
AA_CLONOTYPE_KEY: tuple[str, ...] = (JUNCTION_AA, V_CALL, J_CALL, C_CALL)


def _aa_key(df: "pl.DataFrame | pl.LazyFrame") -> list[str]:
    """The subset of :data:`AA_CLONOTYPE_KEY` this frame actually carries."""
    have = set(column_names(df))
    return [c for c in AA_CLONOTYPE_KEY if c in have]


def has_nt_resolution(df: "pl.DataFrame | pl.LazyFrame") -> bool:
    """Whether ``junction_nt`` can tell two rows sharing an amino-acid key apart.

    True only if the column exists **and** carries at least one non-null value: a schema-conformant
    frame always has the column, and an amino-acid-collapsed export fills it entirely with nulls,
    so presence alone proves nothing.
    """
    if JUNCTION_NT not in column_names(df):
        return False
    lf = df.lazy().select(pl.col(JUNCTION_NT).is_not_null().any().alias("_any"))
    return bool(lf.collect()["_any"][0])


def duplicate_keys(df: pl.DataFrame) -> pl.DataFrame:
    """The repeated amino-acid clonotype keys in a frame, with how many rows each covers.

    Returns:
        ``<key columns>, n`` for every key occurring more than once, worst first. Empty when the
        frame resolves — which includes the case of a frame carrying ``junction_nt``, where a
        repeated amino-acid key is a real pair of nucleotide clonotypes rather than an artefact.
    """
    key = _aa_key(df)
    if df.height == 0 or not key or has_nt_resolution(df):
        return pl.DataFrame(schema={**{c: pl.Utf8 for c in key}, "n": pl.UInt32})
    return (df.group_by(key).agg(pl.len().alias("n"))
              .filter(pl.col("n") > 1).sort("n", descending=True))


def assert_resolvable(df: pl.DataFrame, *, name: str | None = None) -> None:
    """Raise if a frame without ``junction_nt`` repeats an amino-acid clonotype key.

    A **data-integrity gate, not a filter** — the sibling of
    :func:`vdjtools.signature.blocks.assert_parseable`, and it guards a question the frame cannot
    answer about itself. Without ``junction_nt`` there is no way to tell whether two rows sharing
    ``(junction_aa, v_call, j_call, c_call)`` are two nucleotide clonotypes that happen to encode
    the same peptide, or one clonotype duplicated by an export artefact. The two readings give
    different richness, different clonality, different Shannon and a different top-clone fraction,
    and until this gate existed the library counted rows and picked the first reading in silence.

    Measured cost of that silence: two exports of the same samples, one collapsed to the
    amino-acid key and one not, went through the same estimators without complaint. Richness
    differed by **1.5% overall and 5.0% in IGK**, and the affected samples then sat **3-5 robust-SD
    from the rest of the cohort** on exactly the diversity and clonality columns — read as biology
    until the two exports were diffed.

    A frame carrying ``junction_nt`` is never rejected: the duplicates are then real and the key
    that resolves them is present.

    Args:
        df: A clonotype frame.
        name: Sample identifier to name in the message, when there is one.

    Raises:
        ValueError: If the frame cannot resolve a repeated amino-acid key.
    """
    dup = duplicate_keys(df)
    if not dup.height:
        return
    where = f" in sample {name!r}" if name else ""
    ex = "; ".join(
        " ".join(str(v) for v in row[:-1]) + f" x{row[-1]}"
        for row in dup.head(3).iter_rows())
    raise ValueError(
        f"{dup.height} amino-acid clonotype key(s){where} occur more than once and the frame has "
        f"no {JUNCTION_NT} to tell the rows apart, e.g. {ex}. Two rows on one key are either two "
        "nucleotide clonotypes encoding the same peptide or one clonotype duplicated by the "
        "export, and richness, clonality, Shannon and top-clone fraction differ between those "
        "readings — so this is a question about the data, not a default to pick. Resolve it "
        f"either way: supply {JUNCTION_NT} so the rows are distinguishable, or pass "
        'on_duplicate="sum" (equivalently collapse_duplicates(df)) to add the counts together '
        "deliberately.")


def collapse_duplicates(df: pl.DataFrame) -> pl.DataFrame:
    """Sum ``duplicate_count`` over repeated amino-acid clonotype keys, keeping row order.

    The ``on_duplicate="sum"`` half of :func:`assert_resolvable` — the deliberate reading in which
    repeated rows are one clonotype the export split. Non-key columns take their value from the
    first row of each group, and ``frequency`` is recomputed because summing counts invalidates it.
    A frame that resolves (carries ``junction_nt``) or has no repeats is returned unchanged.
    """
    key = _aa_key(df)
    if df.height == 0 or not key or has_nt_resolution(df):
        return df
    out = (df.with_row_index("_i")
             .group_by(key, maintain_order=True)
             .agg(pl.col(COUNT).sum(),
                  pl.exclude([*key, COUNT]).first())
             .drop("_i").select(df.columns))
    return recompute_frequency(out) if FREQ in out.columns else out


def resolve_duplicates(df: pl.DataFrame, on_duplicate: str = "error") -> pl.DataFrame:
    """Apply an ``on_duplicate`` policy: ``"error"`` raises, ``"sum"`` collapses.

    The one entry point analysis code should call, so every path spells the policy the same way.
    """
    if on_duplicate == "error":
        assert_resolvable(df)
        return df
    if on_duplicate == "sum":
        return collapse_duplicates(df)
    raise ValueError(f'on_duplicate must be "error" or "sum", got {on_duplicate!r}')
