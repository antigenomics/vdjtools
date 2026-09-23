"""Clonotype filtering (pure polars).

Reimplements the legacy vdjtools clonotype-filter family:

- :func:`filter_productive` — AIRR-productive rearrangements (supersedes ``filter_functional``,
  the legacy ``FunctionalClonotypeFilter`` / ``isCoding``).
- :func:`filter_frequency` — ``FilterByFrequency`` (``FrequencyFilter`` + ``QuantileFilter``).
- :func:`filter_segment` — ``FilterBySegment`` (``VFilter`` / ``DFilter`` / ``JFilter``).
- :func:`filter_by_sample` — ``ApplySampleAsFilter`` (``IntersectionClonotypeFilter``).

Every filter recomputes ``frequency`` within the surviving subset by default, and
:func:`filter_productive` exposes that as ``recompute_frequencies=False`` for a caller who wants
the file's own frequencies left alone.

**On the word "functional".** It is IMGT's, and IMGT applies it to a *germline gene* (F / ORF / P),
not to a rearrangement. What this module filters is AIRR **productivity** — a property of the
rearranged sequence: no stop codon, junction in frame. The two are orthogonal; a perfectly
productive rearrangement can use a pseudogene V. :func:`filter_productive` is therefore the name,
and :func:`filter_functional` is kept as a deprecated alias.
"""
from __future__ import annotations

import polars as pl

import warnings

from ..io.schema import (
    JUNCTION_AA,
    COUNT,
    D_CALL,
    FREQ,
    J_CALL,
    PRODUCTIVE,
    STOP_CODON,
    V_CALL,
    VJ_IN_FRAME,
    column_names,
    strip_allele,
    recompute_frequency,
)

#: Characters that mark a non-coding CDR3 amino-acid string. ``*`` is the stop
#: codon; the lowercase nucleotides and ``# ~ _ ?`` are the legacy out-of-frame
#: markers emitted when a junction cannot be cleanly translated (legacy
#: ``CommonUtil.OOF_SYMBOLS_POSSIBLE`` + ``STOP_CHAR``).
_NONCODING_CHARS = r"[*atgc#~_?]"


def _productive_from_annotation(df: pl.DataFrame) -> "tuple[pl.Expr, str] | None":
    """Read productivity off the file's own AIRR columns, if it carries any.

    Preference order is the composite first, then its components: ``productive`` states the whole
    of it (open reading frame, no defect in start codon / splicing / regulatory elements, no
    internal stop, junction in frame), while ``stop_codon`` and ``vj_in_frame`` each state a part.
    Re-deriving any of this from ``junction_aa`` when the file already says it is how a caller ends
    up disagreeing with their own annotation without being told.
    """
    cols = set(column_names(df))
    if PRODUCTIVE in cols:
        return pl.col(PRODUCTIVE).cast(pl.Boolean, strict=False), PRODUCTIVE
    parts, used = [], []
    if STOP_CODON in cols:
        parts.append(~pl.col(STOP_CODON).cast(pl.Boolean, strict=False))
        used.append(STOP_CODON)
    if VJ_IN_FRAME in cols:
        parts.append(pl.col(VJ_IN_FRAME).cast(pl.Boolean, strict=False))
        used.append(VJ_IN_FRAME)
    if not parts:
        return None
    expr = parts[0]
    for e in parts[1:]:
        expr = expr & e
    return expr, "+".join(used)


def productive_mask(df: pl.DataFrame) -> "tuple[pl.Expr, str]":
    """The productivity predicate for this frame, and the evidence it rests on.

    Returns ``(expr, source)`` where ``source`` is the AIRR column(s) used, or ``"junction_aa"``
    when none were present and productivity had to be derived from the amino-acid string.

    The fallback reads a stop codon as ``*`` and an out-of-frame junction as one of the legacy
    markers. It is a proxy: it cannot see a defect in a splicing site or a regulatory element,
    which AIRR's ``productive`` can.
    """
    ann = _productive_from_annotation(df)
    if ann is not None:
        expr, src = ann
        return expr.fill_null(False), src
    return (pl.col(JUNCTION_AA).is_not_null()
            & ~pl.col(JUNCTION_AA).str.contains(_NONCODING_CHARS)), JUNCTION_AA


def filter_productive(df: pl.DataFrame, keep: str = "productive", *,
                      recompute_frequencies: bool = True) -> pl.DataFrame:
    """Keep only AIRR-productive rearrangements (or only the complement).

    A rearrangement is *productive* when it can encode a receptor chain. Where the frame carries
    the AIRR annotation columns (``productive``, or ``stop_codon`` / ``vj_in_frame``) those are
    authoritative; otherwise productivity is derived from ``junction_aa``, where a stop codon is
    ``*`` and an out-of-frame junction carries one of the legacy markers ``[atgc#~_?]``. A null
    ``junction_aa`` is treated as non-productive.

    This is **not** IMGT functionality. That is a property of the germline gene (F / ORF / P) and
    is orthogonal to this one — see :func:`filter_functional_genes`.

    Args:
        df: A clonotype frame.
        keep: ``"productive"`` (default) or ``"nonproductive"`` for the complement.
        recompute_frequencies: Renormalise ``frequency`` over the survivors. **Default ``True``**,
            which is the legacy behaviour and what almost every caller wants. Pass ``False`` to
            leave the file's own frequencies untouched — useful when the frequencies are the
            quantity of interest and must stay comparable to the unfiltered file.

    Returns:
        The filtered frame.

    Raises:
        ValueError: If ``keep`` is not ``"productive"`` or ``"nonproductive"``.
    """
    if keep not in ("productive", "nonproductive"):
        raise ValueError(f"keep must be 'productive' or 'nonproductive'; got {keep!r}")
    mask, _ = productive_mask(df)
    out = df.filter(mask if keep == "productive" else ~mask)
    return recompute_frequency(out) if recompute_frequencies else out


def filter_functional(df: pl.DataFrame, keep: str = "coding") -> pl.DataFrame:
    """Deprecated alias for :func:`filter_productive`. Use that instead.

    "Functional" is IMGT's word for a germline gene; this function filters rearrangements, which
    AIRR calls *productive*. Kept for one release so existing callers do not break silently.
    """
    warnings.warn(
        "filter_functional() is deprecated; use filter_productive(). 'functional' is IMGT's term "
        "for a germline gene (F/ORF/P), while this filters rearrangements, which AIRR calls "
        "'productive'. Map keep='coding' -> keep='productive', "
        "keep='noncoding' -> keep='nonproductive'.",
        DeprecationWarning, stacklevel=2)
    if keep not in ("coding", "noncoding"):
        raise ValueError(f"keep must be 'coding' or 'noncoding'; got {keep!r}")
    return filter_productive(df, "productive" if keep == "coding" else "nonproductive")


def filter_functional_genes(df: pl.DataFrame, *, segments: tuple[str, ...] = ("V", "J"),
                            keep: tuple[str, ...] = ("F",), organism: str = "human",
                            locus: str | None = None,
                            recompute_frequencies: bool = True) -> pl.DataFrame:
    """Keep rearrangements whose germline gene calls are IMGT-functional.

    This is the **other** axis, and the one that actually deserves the word *functional*. IMGT
    classifies a germline gene as **F** (functional), **ORF** (an open reading frame, but a defect
    in splicing, regulatory elements or conserved-residue hydropathy — *not* functional), or **P**
    (pseudogene: a defect in the ORF itself). See
    https://www.imgt.org/IMGTindex/functionality.php

    It is orthogonal to :func:`filter_productive`. A rearrangement can be perfectly in frame with
    no stop codon — AIRR-productive — while using a pseudogene V; and a functional V gene can
    rearrange out of frame. Filtering one says nothing about the other.

    A call this function cannot resolve against the germline reference is **kept**, not dropped:
    an unrecognised gene name means our reference is incomplete or the caller uses a different
    nomenclature, and silently discarding those rows would be a vocabulary bug reported as biology.

    Args:
        df: A clonotype frame.
        segments: Which calls to check — any of ``"V"``, ``"D"``, ``"J"``.
        keep: IMGT functionality codes to keep. ``("F",)`` is strict; ``("F", "ORF")`` is the
            common looser choice.
        organism: Passed to the germline reference.
        locus: Locus to load the reference for. Inferred from the V calls when omitted.
        recompute_frequencies: Renormalise ``frequency`` over the survivors. Default ``True``.

    Returns:
        The filtered frame.
    """
    from ..model.reference import load_germline

    if locus is None:
        v = df[V_CALL].drop_nulls()
        if not v.len():
            return df
        locus = str(v[0])[:3].upper()
    germ = load_germline(locus, organism)
    ok = {seg: set(germ.filter((pl.col("segment") == seg)
                               & pl.col("functionality").is_in(list(keep)))["gene"].to_list())
          for seg in segments}
    known = {seg: set(germ.filter(pl.col("segment") == seg)["gene"].to_list()) for seg in segments}

    mask = pl.lit(True)
    for seg, col in (("V", V_CALL), ("D", D_CALL), ("J", J_CALL)):
        if seg not in segments or col not in column_names(df):
            continue
        gene = strip_allele(pl.col(col))
        # unknown -> kept: an unrecognised name is a vocabulary gap, not a pseudogene
        mask = mask & (gene.is_in(list(ok[seg])) | ~gene.is_in(list(known[seg]))
                       | gene.is_null())
    out = df.filter(mask)
    return recompute_frequency(out) if recompute_frequencies else out


#: Default junction_aa length bounds, INCLUSIVE, in amino acids.
#:
#: A CDR3 shorter than 5 aa cannot span the Cys104..Phe118 anchors with any diversity between them,
#: and a junction longer than 60 aa is beyond anything the germline can produce -- both are almost
#: always a misparse or a chimeric read rather than a receptor. Deliberately wide: this is a
#: sanity bound, not a biological filter, and real junctions sit far inside it.
MIN_JUNCTION_AA, MAX_JUNCTION_AA = 5, 60


def filter_length(df: pl.DataFrame, *, min_len: int = MIN_JUNCTION_AA,
                  max_len: int = MAX_JUNCTION_AA, keep: str = "within",
                  recompute_frequencies: bool = True) -> pl.DataFrame:
    """Keep clonotypes whose ``junction_aa`` length is within bounds, **inclusive**.

    Both bounds are inclusive: ``min_len=5`` keeps a 5-mer, ``max_len=60`` keeps a 60-mer.

    This is a data-sanity filter, not a receptor-biology one. It catches misparses, truncated
    reads and chimeras; it does not encode a claim about what lengths are immunologically
    interesting. Note that neither this nor :func:`filter_productive` filters on length by
    default -- nothing upstream in this package has ever imposed a length bound, so switching this
    on will change counts on any corpus that carries junk.

    **CDR3 vs junction.** These bounds are on ``junction_aa``, which *includes* the Cys104 and
    Phe118 anchors and is therefore two residues longer than the IMGT CDR3. Subtract 2 if you are
    reasoning in CDR3 lengths.

    A null ``junction_aa`` is dropped by ``keep="within"`` -- an absent junction has no length.

    Args:
        df: A clonotype frame.
        min_len: Shortest ``junction_aa`` to keep, inclusive.
        max_len: Longest ``junction_aa`` to keep, inclusive.
        keep: ``"within"`` (default) or ``"outside"`` for the complement -- useful for inspecting
            what a bound would discard before committing to it.
        recompute_frequencies: Renormalise ``frequency`` over the survivors. Default ``True``.

    Returns:
        The filtered frame.

    Raises:
        ValueError: If ``keep`` is unknown, or ``min_len`` exceeds ``max_len``.
    """
    if keep not in ("within", "outside"):
        raise ValueError(f"keep must be 'within' or 'outside'; got {keep!r}")
    if min_len > max_len:
        raise ValueError(f"min_len ({min_len}) exceeds max_len ({max_len})")
    n = pl.col(JUNCTION_AA).str.len_chars()
    within = pl.col(JUNCTION_AA).is_not_null() & (n >= min_len) & (n <= max_len)
    out = df.filter(within if keep == "within" else ~within)
    return recompute_frequency(out) if recompute_frequencies else out


def filter_frequency(df: pl.DataFrame, min_freq: float | None = None,
                     top_quantile: float | None = None) -> pl.DataFrame:
    """Keep abundant clonotypes by frequency threshold and/or top quantile.

    Reimplements ``FilterByFrequency`` (a composite of ``FrequencyFilter`` and
    ``QuantileFilter``). Both criteria, when given, are combined with AND:

    - ``min_freq``: keep clonotypes with ``frequency >= min_freq``.
    - ``top_quantile``: keep the top clonotypes (by ``duplicate_count``) whose
      *cumulative* original frequency, including the clonotype itself, is at most
      ``top_quantile`` of the full-sample total frequency. This matches the legacy
      ``QuantileFilter``: it walks the count-sorted sample accumulating frequency
      and drops the first clonotype that would push the running fraction above the
      threshold (so ``top_quantile=0.25`` keeps roughly the top 25% of the read
      mass). The denominator is the full-sample frequency total (~1.0), and only
      clonotypes that already passed ``min_freq`` contribute to the cumulative
      (legacy filters short-circuit in the order count/freq/quantile).

    Args:
        df: A clonotype frame with ``duplicate_count`` and ``frequency`` columns.
        min_freq: Minimum per-clonotype frequency (e.g. legacy default ``0.01``).
            ``None`` disables it.
        top_quantile: Top read-mass quantile to retain (e.g. legacy default
            ``0.25``). ``None`` disables it.

    Returns:
        The filtered frame, sorted by descending ``duplicate_count``, with
        ``frequency`` recomputed.
    """
    out = df.sort(COUNT, descending=True, maintain_order=True)
    total_freq = out[FREQ].sum()  # legacy parent.getFreqAsInInput(): full-sample total
    if min_freq is not None:
        out = out.filter(pl.col(FREQ) >= min_freq)
    if top_quantile is not None and total_freq:
        cutoff = top_quantile * total_freq
        out = out.filter(pl.col(FREQ).cum_sum() <= cutoff)
    return recompute_frequency(out)


def _segment_matches(col: str, names: list[str]) -> pl.Expr:
    """Boolean expression: does this segment call match any query name (prefix)?

    Incomplete query names act as wildcards (legacy ``getAtFuzzy``): a match is a
    prefix match, so ``TRBV12`` matches ``TRBV12-3*01`` (allele-insensitive) while
    ``TRBV12-3*01`` matches only itself.

    The match is per **comma-token**, not on the raw string: an AIRR call may be an
    ambiguity tie like ``IGHV3-23*01,IGHV3-23D*01``, and a ``starts_with`` on the whole
    string only ever tests the FIRST gene -- so ``filter_segment(v=["IGHV3-23D"])`` would
    match none of the tens of thousands of rows that name IGHV3-23D only in a tie.
    """
    # Regex over the whole call: `name` at the start of the string OR right after a comma
    # (+ optional whitespace). re.escape guards names containing regex metacharacters.
    import re
    pat = "(^|,)\\s*(" + "|".join(re.escape(n) for n in names) + ")"
    return pl.col(col).str.contains(pat)


def filter_segment(df: pl.DataFrame, v: list[str] | None = None,
                   d: list[str] | None = None, j: list[str] | None = None,
                   keep: bool = True) -> pl.DataFrame:
    """Keep or remove clonotypes by V/D/J segment membership.

    Reimplements ``FilterBySegment``. A clonotype *matches* when its V segment is
    in ``v`` **and** its D segment in ``d`` **and** its J segment in ``j`` (only the
    lists that are supplied constrain; unsupplied loci always pass). Matching is a
    prefix match, so incomplete names act as wildcards and are allele-insensitive
    (``TRBV12`` matches ``TRBV12-3*01``).

    Args:
        df: A clonotype frame with ``v_call`` / ``d_call`` / ``j_call`` columns.
        v: V-segment query names (prefixes). ``None`` leaves V unconstrained.
        d: D-segment query names. ``None`` leaves D unconstrained.
        j: J-segment query names. ``None`` leaves J unconstrained.
        keep: If ``True`` (default) keep matching clonotypes; if ``False`` remove
            them (legacy ``--negative``).

    Returns:
        The filtered frame with ``frequency`` recomputed.
    """
    match = pl.lit(True)
    for col, names in ((V_CALL, v), (D_CALL, d), (J_CALL, j)):
        if names:
            match = match & _segment_matches(col, names)
    out = df.filter(match if keep else ~match)
    return recompute_frequency(out)


def filter_by_sample(df: pl.DataFrame, other: pl.DataFrame, keep: bool = True,
                     key: "tuple[str, ...]" = (JUNCTION_AA, V_CALL, J_CALL)) -> pl.DataFrame:
    """Keep or remove clonotypes by exact-key presence in another sample.

    Reimplements ``ApplySampleAsFilter`` / ``IntersectionClonotypeFilter``: build
    the key set of ``other`` and keep (or, with ``keep=False``, remove) the
    clonotypes of ``df`` whose key is present in it. Matching is an exact match on
    the ``key`` columns.

    Args:
        df: The clonotype frame to filter.
        other: The filter sample; only its ``key`` columns are used.
        keep: If ``True`` (default) keep clonotypes present in ``other``; if
            ``False`` remove them (legacy ``--negative``).
        key: Columns forming the match key (default
            ``("junction_aa", "v_call", "j_call")`` — legacy "strict"-style at the aa
            level).

    Returns:
        The filtered frame with ``frequency`` recomputed.
    """
    key = list(key)
    keyset = other.select(key).unique()
    out = df.join(keyset, on=key, how="semi" if keep else "anti")
    return recompute_frequency(out)
