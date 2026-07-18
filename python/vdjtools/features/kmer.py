"""CDR3 amino-acid k-mer profiles and joint V + k-mer + C feature summaries.

K-mers are overlapping sliding windows over ``junction_aa``. Each k-mer occurrence
carries its clonotype's weight (reads, unique, or frequency); weights are summed.

**Anchors carry no information.** A junction begins with the conserved Cys104 and a few germline
V residues and ends with germline J: on a 250k human-TRB control, an N-terminal 4-mer is shared by
**31.0%** of clonotypes (``CASS`` alone by **56.5%**) while a *central* 4-mer is shared by
**0.080%** — ~386x more selective (measured in :mod:`seqtree.seeds`). So a window that starts at
index 0 spends most of its multiple-testing budget on germline. Pass ``flank`` to drop that many
residues from each end first; ``flank=4`` matches ``seqtree.seeds.core_kmers``.

``flank`` is a *fixed* trim, so it is only an approximation of the germline boundary — the real
one is per-clonotype and locus-specific (TRBJ2-3's germline is ``STDTQYF``, seven residues, so
``flank=4`` still leaves ``STD`` in the core). Where the data carries V/J markup, prefer it.
"""
from __future__ import annotations

import polars as pl

from ..io.schema import (
    C_CALL,
    JUNCTION_AA,
    LOCUS,
    V_CALL,
    add_locus,
    column_names,
    strip_allele,
    weight_expr,
)


def _explode_kmers(df: pl.DataFrame, k: int, flank: int = 0) -> pl.DataFrame:
    """Explode each clonotype's ``junction_aa`` into overlapping k-mers (column ``kmer``).

    Rows whose CDR3 is shorter than ``k`` (or, with ``flank``, whose core is) contribute no
    k-mers. ``flank`` drops that many residues from each end before windowing.
    """
    if k < 1:
        raise ValueError(f"k must be >= 1; got {k}")
    if flank < 0:
        raise ValueError(f"flank must be >= 0; got {flank}")
    # `str.len_chars()` is UInt32, so `len - 2*flank` UNDERFLOWS to ~4e9 on a junction shorter
    # than the flanks instead of going negative — cast before subtracting, then clamp at 0: such
    # a junction has NO core, and both the underflow and a bare negative make polars' str.slice
    # silently return a germline tail (`CASSFGH`, flank=4 -> "FGH"). seqtree.seeds.core_kmers
    # returns nothing there, and it is right.
    df = df.with_columns(pl.col(JUNCTION_AA).str.len_chars().cast(pl.Int64).alias("_len"))
    core = (pl.col(JUNCTION_AA) if not flank
            else pl.col(JUNCTION_AA).str.slice(
                flank, pl.max_horizontal(pl.col("_len") - 2 * flank, pl.lit(0))))
    df = df.with_columns(core.alias("_core"))
    df = df.with_columns(pl.col("_core").str.len_chars().alias("_clen"))
    df = df.filter(pl.col("_clen") >= k)  # too-short core (or CDR3) yields no k-mers
    df = df.with_columns(pl.int_ranges(0, pl.col("_clen") - k + 1).alias("_pos"))
    df = df.explode("_pos", empty_as_null=True)
    df = df.with_columns(pl.col("_core").str.slice(pl.col("_pos"), k).alias("kmer"))
    return df.drop("_len", "_core", "_clen", "_pos")


def kmer_profile(df, k: int = 3, weight: str = "reads",
                 by_locus: bool = True, flank: int = 0, by=()):
    """CDR3 amino-acid k-mer spectrum.

    Args:
        df: A clonotype frame (eager ``pl.DataFrame`` or lazy ``pl.LazyFrame`` — e.g. a
            whole cohort from :func:`vdjtools.io.scan_cohort`); result mirrors the input.
        k: K-mer length (default 3).
        weight: ``"reads"``, ``"unique"``, or ``"freq"``.
        by_locus: If ``True``, break the spectrum down per locus.
        flank: Residues to drop from each end before windowing (see the module docstring);
            ``0`` (default) keeps the whole junction, ``4`` gives the ``seqtree`` core.
        by: Extra column(s) to **prepend** to the group key (e.g. ``["sample_id"]`` for
            a one-pass cohort spectrum). Empty by default — byte-identical to per-sample.

    Returns:
        Long-format frame with ``kmer``, ``weight``, the ``by`` columns and (if
        ``by_locus``) ``locus``, sorted by the group key (lazy when ``df`` is lazy).
    """
    df = df if LOCUS in column_names(df) else add_locus(df)
    df = df.with_columns(weight_expr(weight).alias("weight"))
    df = _explode_kmers(df, k, flank)
    group = [*by, *([LOCUS, "kmer"] if by_locus else ["kmer"])]
    return (df.group_by(group, maintain_order=True)
              .agg(pl.col("weight").sum())
              .sort(group))


def v_kmer_c_profile(df, k: int = 3, weight: str = "reads",
                     by_locus: bool = True, keep_allele: bool = False,
                     flank: int = 0, by=()):
    """Joint (V gene, k-mer, C gene) profile — a tidy feature-matrix source.

    Produces one aggregated row per ``(v_call, kmer, c_call)`` combination (plus
    ``locus`` if requested), suitable to pivot into a feature matrix. A null
    ``c_call`` (common in native vdjtools data) is retained as its own group.

    Args:
        df: A clonotype frame (eager or lazy; see :func:`kmer_profile`).
        k: K-mer length (default 3).
        weight: ``"reads"``, ``"unique"``, or ``"freq"``.
        by_locus: If ``True``, include ``locus`` in the grouping.
        keep_allele: If ``True``, keep V allele suffixes; otherwise collapse V to
            gene level (default).
        flank: Residues to drop from each end before windowing (see the module docstring).
        by: Extra column(s) prepended to the group key (e.g. ``["sample_id"]`` for a
            one-pass cohort feature source).

    Returns:
        Long-format frame with ``v_call``, ``kmer``, ``c_call``, ``weight``, the ``by``
        columns and (if ``by_locus``) ``locus`` (lazy when ``df`` is lazy).
    """
    df = df if LOCUS in column_names(df) else add_locus(df)
    v = pl.col(V_CALL) if keep_allele else strip_allele(pl.col(V_CALL))
    df = df.with_columns(v.alias(V_CALL), weight_expr(weight).alias("weight"))
    df = _explode_kmers(df, k, flank)
    group = [*by, *([LOCUS, V_CALL, "kmer", C_CALL] if by_locus
                    else [V_CALL, "kmer", C_CALL])]
    return (df.group_by(group, maintain_order=True)
              .agg(pl.col("weight").sum())
              .sort(group, nulls_last=True))


def kmer_cohort(cohort: pl.DataFrame | pl.LazyFrame, k: int = 4, flank: int = 4,
                keep_allele: bool = False) -> pl.DataFrame:
    """Explode a cohort into a per-``(sample_id, v_call, kmer)`` frame for an incidence test.

    The bridge from k-mers to :func:`vdjtools.biomarker.association`: the result carries a
    ``kmer`` column alongside ``sample_id``, so a V + k-mer phenotype test is just

    .. code-block:: python

        from vdjtools.biomarker import association
        from vdjtools.features.kmer import kmer_cohort
        km = kmer_cohort(cohort, k=4, flank=4)
        res = association(km, design, key=("v_call", "kmer"), match="exact")

    ``match="exact"`` is the only valid mode here — ``fuzzy``/``1mm`` search on the CDR3, and a
    1-mismatch ball around a 4-mer is most of k-mer space.

    Pinning V matters. A short central k-mer is only selective *given* the germline context: on
    its own a 4-mer is shared by ~0.08% of a control repertoire, which over ~10⁵ features is still
    thousands of hits. ``(v_call, kmer)`` is the feature the V+k-mer search is named for.

    Args:
        cohort: A clonotype frame with ``sample_id`` (e.g. from :func:`vdjtools.io.scan_cohort`).
        k: K-mer length. Note :mod:`seqtree.seeds`' measurement that a central k-mer's median
            E-value crosses 1 at ``k=6``; below that, prune with ``seqtree.seeds.SeedIndex``
            against a real control rather than trusting a modelled background (its residual KL
            *grows* with k — D-gene germline runs correlate — so it must be counted, not fitted).
        flank: Residues dropped from each end (see the module docstring). ``4`` = seqtree's core.
        keep_allele: Keep V allele suffixes; default collapses to gene level.

    Returns:
        Unique ``(sample_id, v_call, kmer)`` rows.
    """
    from ..io.cohort import SAMPLE_ID

    df = cohort.lazy().collect() if isinstance(cohort, pl.LazyFrame) else cohort
    if SAMPLE_ID not in df.columns:
        raise ValueError(f"cohort must carry {SAMPLE_ID!r}; got {df.columns}")
    v = pl.col(V_CALL) if keep_allele else strip_allele(pl.col(V_CALL))
    df = df.with_columns(v.alias(V_CALL))
    return _explode_kmers(df, k, flank).select(SAMPLE_ID, V_CALL, "kmer").unique()
