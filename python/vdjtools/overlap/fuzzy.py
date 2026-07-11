"""Fuzzy (edit-distance) sample overlap — delegated to vdjmatch/seqtree.

Where :mod:`vdjtools.overlap.metrics` matches clonotypes *exactly*, this module
matches CDR3s **within an edit-distance scope** (substitutions/indels), so a pair of
repertoires that share only near-variants of a clonotype still register as
overlapping. The fuzzy search itself is **not reimplemented here** — it is delegated
to :func:`vdjmatch.cluster.overlap`, which runs on the native ``seqtree`` engine.
These functions are thin polars wrappers that (a) collapse each sample to unique
CDR3s, (b) hand the CDR3 lists to vdjmatch, and (c) join the matched pairs back to
per-clonotype counts/frequencies.

Scope syntax is vdjmatch's ``"subs,ins,dels,total"`` (max substitutions, insertions,
deletions, and total edits); the default ``"1,0,0,1"`` is a single substitution.
"""
from __future__ import annotations

import numpy as np
import polars as pl

from ..io.schema import CDR3_AA, COUNT

_VDJMATCH_HINT = (
    "vdjmatch is required for vdjtools.overlap.fuzzy; install the extra with "
    "`pip install 'vdjtools[overlap]'` (or `pip install vdjmatch>=0.0.1`)."
)


def _require_vdjmatch():
    """Import and return ``vdjmatch.cluster``; raise a helpful error if missing."""
    try:
        import vdjmatch.cluster as cluster  # noqa: F401
    except ImportError as exc:  # pragma: no cover - exercised only without vdjmatch
        raise ImportError(_VDJMATCH_HINT) from exc
    return cluster


def _aggregate(df: pl.DataFrame, by: str) -> pl.DataFrame:
    """Collapse to unique ``by`` values with summed count, within-sample frequency,
    and a 0-based row index (the positional key vdjmatch returns)."""
    agg = df.group_by(by, maintain_order=True).agg(pl.col(COUNT).sum().alias("_count"))
    total = agg["_count"].sum() or 1
    return agg.with_columns(
        (pl.col("_count") / total).alias("_freq"),
        pl.int_range(pl.len(), dtype=pl.Int64).alias("_idx"),
    )


def fuzzy_overlap(a: pl.DataFrame, b: pl.DataFrame, scope: str = "1,0,0,1",
                  by: str = CDR3_AA, threads: int = 0) -> pl.DataFrame:
    """Fuzzy-matched clonotype pairs between two samples (within an edit scope).

    Both samples are collapsed to unique ``by`` values (counts summed, frequencies
    recomputed within each sample); the two CDR3 lists are handed to
    :func:`vdjmatch.cluster.overlap`, and the returned within-scope pairs are joined
    back to per-clonotype counts and frequencies.

    Args:
        a: First clonotype frame (canonical schema).
        b: Second clonotype frame.
        scope: vdjmatch edit-distance scope ``"subs,ins,dels,total"`` (default one
            substitution).
        by: The CDR3 column to match on (default ``"cdr3_aa"``).
        threads: Worker threads for the native search (``0`` = all cores).

    Returns:
        A ``pl.DataFrame``, one row per matched pair, with columns ``a_cdr3, b_cdr3,
        n_subs, score, count_a, freq_a, count_b, freq_b``. Empty (with that schema)
        when nothing matches.

    Raises:
        ImportError: If vdjmatch is not installed (see the ``overlap`` extra).
    """
    cluster = _require_vdjmatch()
    a_agg = _aggregate(a, by)
    b_agg = _aggregate(b, by)
    pairs = cluster.overlap(a_agg[by].to_list(), b_agg[by].to_list(),
                            scope=scope, threads=threads)

    out_schema = ["a_cdr3", "b_cdr3", "n_subs", "score",
                  "count_a", "freq_a", "count_b", "freq_b"]
    if pairs.height == 0:
        return pl.DataFrame(schema={
            "a_cdr3": pl.Utf8, "b_cdr3": pl.Utf8, "n_subs": pl.Int64, "score": pl.Int64,
            "count_a": pl.Int64, "freq_a": pl.Float64,
            "count_b": pl.Int64, "freq_b": pl.Float64,
        })

    a_join = a_agg.rename({"_idx": "a_idx", "_count": "count_a", "_freq": "freq_a"}).drop(by)
    b_join = b_agg.rename({"_idx": "b_idx", "_count": "count_b", "_freq": "freq_b"}).drop(by)
    return (pairs.join(a_join, on="a_idx", how="left")
                 .join(b_join, on="b_idx", how="left")
                 .select(out_schema))


def fuzzy_overlap_metrics(a: pl.DataFrame, b: pl.DataFrame, scope: str = "1,0,0,1",
                          by: str = CDR3_AA, threads: int = 0) -> dict:
    """Summary fuzzy-overlap metrics for two samples.

    Computes :func:`fuzzy_overlap` once and derives:

    - **pairs** — number of within-scope matched clonotype pairs.
    - **frac_a_matched** / **frac_b_matched** — fraction of each sample's unique
      clonotypes with at least one fuzzy match in the other.
    - **fuzzy_F** — the frequency-weighted fuzzy analogue of the exact ``F`` metric,
      ``sqrt(Σ_{a matched} freq_a · Σ_{b matched} freq_b)``, where a clonotype's
      frequency is counted once if it has any fuzzy neighbour in the other sample
      (so a clonotype matching several near-variants is not double-weighted).

    Args:
        a: First clonotype frame.
        b: Second clonotype frame.
        scope: vdjmatch edit-distance scope (see :func:`fuzzy_overlap`).
        by: The CDR3 column to match on.
        threads: Worker threads for the native search.

    Returns:
        Dict with keys ``pairs, frac_a_matched, frac_b_matched, fuzzy_F``.

    Raises:
        ImportError: If vdjmatch is not installed (see the ``overlap`` extra).
    """
    a_agg = _aggregate(a, by)
    b_agg = _aggregate(b, by)
    n_a, n_b = a_agg.height, b_agg.height
    pairs = fuzzy_overlap(a, b, scope=scope, by=by, threads=threads)

    if pairs.height == 0:
        return {"pairs": 0, "frac_a_matched": 0.0, "frac_b_matched": 0.0, "fuzzy_F": 0.0}

    # Unique matched clonotypes on each side, and their once-counted frequency mass.
    a_matched = pairs.unique(subset="a_cdr3")
    b_matched = pairs.unique(subset="b_cdr3")
    mass_a = float(a_matched["freq_a"].sum())
    mass_b = float(b_matched["freq_b"].sum())
    return {
        "pairs": pairs.height,
        "frac_a_matched": a_matched.height / n_a if n_a else 0.0,
        "frac_b_matched": b_matched.height / n_b if n_b else 0.0,
        "fuzzy_F": float(np.sqrt(mass_a * mass_b)),
    }
