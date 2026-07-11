"""Random down-sampling of clonotype frames (pure polars + numpy).

Reimplements the legacy vdjtools ``DownSampler`` / ``TopSampler`` family. Two
resampling regimes, matching the legacy ``--unweighted`` switch:

- **reads** (legacy default, weighted): draw ``size`` reads *without replacement*
  from the multiset of reads implied by ``duplicate_count``. The legacy
  ``DownSampler`` shuffles a flattened per-read array and keeps the first ``size``
  entries; the exact equivalent is the multivariate hypergeometric distribution
  (``numpy.random.Generator.multivariate_hypergeometric``). The task brief phrased
  this as "multinomial", but sampling a sequencing library to a fixed depth is a
  *without-replacement* operation — a multinomial (with replacement) could return
  more reads of a clonotype than were observed — so the hypergeometric is used to
  stay faithful to the legacy behaviour and to the biology.
- **clones** (legacy ``--unweighted``): draw ``size`` unique clonotypes *uniformly*
  at random without replacement, keeping each one's original count. Note the legacy
  clonotype-level mode is uniform, not count-weighted (weighting by count is what
  the read-level mode does).
"""
from __future__ import annotations

import numpy as np
import polars as pl

from ..io.schema import COUNT, recompute_frequency


def downsample(df: pl.DataFrame, size: int, by: str = "reads",
               seed: int = 0) -> pl.DataFrame:
    """Randomly down-sample a clonotype frame to a target size.

    Args:
        df: A clonotype frame with a ``duplicate_count`` column.
        size: Target size — number of reads (``by="reads"``) or number of unique
            clonotypes (``by="clones"``).
        by: ``"reads"`` (default) draws ``size`` reads without replacement,
            weighted by ``duplicate_count`` (multivariate hypergeometric);
            ``"clones"`` draws ``size`` unique clonotypes uniformly without
            replacement, keeping their original counts.
        seed: Seed for the numpy random generator (reproducible output).

    Returns:
        A new clonotype frame with ``frequency`` recomputed. Clonotypes that drew
        zero reads (``by="reads"``) are dropped. If ``size`` is greater than or
        equal to the available size the input is returned unchanged (legacy guard).

    Raises:
        ValueError: If ``by`` is not ``"reads"`` or ``"clones"``, or ``size`` < 0.
    """
    if size < 0:
        raise ValueError(f"size must be non-negative; got {size}")
    rng = np.random.default_rng(seed)

    if by == "reads":
        counts = df[COUNT].to_numpy()
        total = int(counts.sum())
        if size >= total:
            return df
        # numpy's default method='marginals' raises once sum(counts) >= 1e9; use the
        # slower-but-unbounded method='count' there (else keep the faster default).
        method = "count" if total >= 1_000_000_000 else "marginals"
        drawn = rng.multivariate_hypergeometric(counts.astype(np.int64), size, method=method)
        out = df.with_columns(pl.Series(COUNT, drawn, dtype=pl.Int64))
        out = out.filter(pl.col(COUNT) > 0)
        return recompute_frequency(out)

    if by == "clones":
        n = df.height
        if size >= n:
            return df
        idx = rng.choice(n, size=size, replace=False)
        idx.sort()
        out = df[idx.tolist()]
        return recompute_frequency(out)

    raise ValueError(f"by must be 'reads' or 'clones'; got {by!r}")


def select_top(df: pl.DataFrame, n: int, renormalize: bool = True) -> pl.DataFrame:
    """Select the top ``n`` clonotypes by ``duplicate_count``.

    Reimplements the legacy ``SelectTop`` / ``TopSampler`` (take the ``n`` largest
    clonotypes). Ties are broken by the frame's existing order (a stable sort).

    Args:
        df: A clonotype frame with a ``duplicate_count`` column.
        n: Number of top clonotypes to keep. If ``n`` is greater than or equal to
            the number of clonotypes, all are returned.
        renormalize: If ``True`` (legacy default), recompute ``frequency`` within
            the selected subset so it sums to 1; if ``False``, preserve the input
            frequencies (legacy ``--save-freqs``).

    Returns:
        The top-``n`` clonotype frame, sorted by descending ``duplicate_count``.
    """
    out = df.sort(COUNT, descending=True, maintain_order=True).head(n)
    if renormalize:
        out = recompute_frequency(out)
    return out
