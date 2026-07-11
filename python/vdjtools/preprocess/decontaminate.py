"""Cross-sample contamination removal (pure polars).

Reimplements the legacy ``Decontaminate`` / ``RatioFilter``: drop clonotypes that
occur in another sample at a much higher abundance, on the assumption that such a
clonotype leaked (index hopping / carry-over) into the current sample.
"""
from __future__ import annotations

import polars as pl

from ..io.schema import (
    CDR3_NT,
    COUNT,
    FREQ,
    J_CALL,
    V_CALL,
    normalize,
    recompute_frequency,
)


def decontaminate(df: pl.DataFrame, others: "list[pl.DataFrame]", ratio: float = 20.0,
                  by: str = "freq", key: "tuple[str, ...]" = (CDR3_NT, V_CALL, J_CALL)
                  ) -> pl.DataFrame:
    """Remove clonotypes dominated by an exact match in another sample.

    Reimplements ``RatioFilter``. A clonotype of ``df`` is **removed** when some
    sample in ``others`` carries the same clonotype (exact match on ``key``) at an
    abundance ``>= ratio`` times its abundance here (legacy keeps a clonotype iff
    ``max_other < here * ratio``, i.e. removes on ``>=``). Abundance is the
    within-sample ``frequency`` (``by="freq"``) or read ``duplicate_count``
    (``by="reads"``).

    Note:
        The legacy ``Decontaminate --read-based`` branch was a no-op — both CLI
        branches constructed the same frequency-based ``RatioFilter`` (a documented
        TODO). ``by="reads"`` here implements the intended read-count comparison.

    Args:
        df: The clonotype frame to decontaminate.
        others: The other samples that may be contamination sources.
        ratio: Parent-to-child abundance ratio (legacy default ``20``).
        by: ``"freq"`` (default) compares within-sample frequencies; ``"reads"``
            compares raw read counts.
        key: Exact match key (legacy default is the strict key
            ``("cdr3_nt", "v_call", "j_call")``).

    Returns:
        The decontaminated frame with ``frequency`` recomputed.

    Raises:
        ValueError: If ``by`` is not ``"freq"`` or ``"reads"``.
    """
    if by == "freq":
        col = FREQ
    elif by == "reads":
        col = COUNT
    else:
        raise ValueError(f"by must be 'freq' or 'reads'; got {by!r}")
    key = list(key)

    d = normalize(df, recompute_freq=True)
    if not others:
        return recompute_frequency(d)

    # Per-source abundance for the key, then the max across all other samples.
    per_source = [
        normalize(o, recompute_freq=True).group_by(key, maintain_order=True)
        .agg(pl.col(col).sum().alias("_ab"))
        for o in others
    ]
    max_other = (
        pl.concat(per_source, how="vertical_relaxed")
        .group_by(key, maintain_order=True)
        .agg(pl.col("_ab").max().alias("_maxo"))
    )

    joined = d.join(max_other, on=key, how="left")
    keep = pl.col("_maxo").is_null() | (pl.col("_maxo") < pl.col(col) * ratio)
    out = joined.filter(keep).drop("_maxo")
    return recompute_frequency(out)
