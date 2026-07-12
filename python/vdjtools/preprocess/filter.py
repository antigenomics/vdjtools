"""Clonotype filtering (pure polars).

Reimplements the legacy vdjtools clonotype-filter family:

- :func:`filter_functional` — ``FunctionalClonotypeFilter`` (``isCoding``).
- :func:`filter_frequency` — ``FilterByFrequency`` (``FrequencyFilter`` + ``QuantileFilter``).
- :func:`filter_segment` — ``FilterBySegment`` (``VFilter`` / ``DFilter`` / ``JFilter``).
- :func:`filter_by_sample` — ``ApplySampleAsFilter`` (``IntersectionClonotypeFilter``).

Every filter recomputes ``frequency`` within the surviving subset (legacy default;
the ``--save-freqs`` behaviour is left to the caller).
"""
from __future__ import annotations

import polars as pl

from ..io.schema import (
    JUNCTION_AA,
    COUNT,
    D_CALL,
    FREQ,
    J_CALL,
    V_CALL,
    recompute_frequency,
)

#: Characters that mark a non-coding CDR3 amino-acid string. ``*`` is the stop
#: codon; the lowercase nucleotides and ``# ~ _ ?`` are the legacy out-of-frame
#: markers emitted when a junction cannot be cleanly translated (legacy
#: ``CommonUtil.OOF_SYMBOLS_POSSIBLE`` + ``STOP_CHAR``).
_NONCODING_CHARS = r"[*atgc#~_?]"


def filter_functional(df: pl.DataFrame, keep: str = "coding") -> pl.DataFrame:
    """Keep only coding (or only non-coding) clonotypes.

    A clonotype is *coding* when it is in-frame **and** has no stop codon, derived
    from ``junction_aa``: legacy ``isCoding() = isInFrame() && isNoStop()``, where
    ``isNoStop`` means no ``*`` and ``isInFrame`` means no out-of-frame marker
    (``[atgc#~_?]``) in the amino-acid string. A null ``junction_aa`` is treated as
    non-coding.

    Args:
        df: A clonotype frame with a ``junction_aa`` column.
        keep: ``"coding"`` (default) keeps coding clonotypes; ``"noncoding"`` keeps
            the complement.

    Returns:
        The filtered frame with ``frequency`` recomputed.

    Raises:
        ValueError: If ``keep`` is not ``"coding"`` or ``"noncoding"``.
    """
    if keep not in ("coding", "noncoding"):
        raise ValueError(f"keep must be 'coding' or 'noncoding'; got {keep!r}")
    is_coding = pl.col(JUNCTION_AA).is_not_null() & ~pl.col(JUNCTION_AA).str.contains(_NONCODING_CHARS)
    out = df.filter(is_coding if keep == "coding" else ~is_coding)
    return recompute_frequency(out)


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
    prefix match on the raw call, so ``TRBV12`` matches ``TRBV12-3*01`` (allele-
    insensitive) while ``TRBV12-3*01`` matches only itself.
    """
    return pl.any_horizontal([pl.col(col).str.starts_with(name) for name in names])


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
