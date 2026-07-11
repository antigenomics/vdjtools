"""Join clonotypes across samples into a joint table (pure polars).

Reimplements the legacy ``JoinSamples`` / ``JointSample`` / ``JointClonotype``:
build the table of clonotypes present in at least ``min_samples`` samples, keep
each member's per-sample frequency, and summarise each joint clonotype by the
geometric mean of its member frequencies.
"""
from __future__ import annotations

import polars as pl

from ..io.schema import COUNT, FREQ, normalize
from .pool import resolve_key

#: Legacy ``MathUtil.JITTER`` — the epsilon added to each per-sample frequency so
#: that a geometric mean over samples where the clonotype is absent stays finite.
JITTER = 1e-9


def join_samples(samples: "list[pl.DataFrame]", key: "str | tuple[str, ...]" = "aa",
                 min_samples: int = 2, names: "list[str] | None" = None) -> pl.DataFrame:
    """Join clonotypes across samples, keyed by an overlap-type match key.

    Reimplements ``JointSample`` (legacy default key ``"aa"``, ``times-detected``
    2). A joint clonotype is kept when it is present in at least ``min_samples``
    samples (legacy ``OccurrenceJoinFilter``). Its joint frequency is the geometric
    mean of the per-sample frequencies (absent samples contribute ``JITTER``),

    ``base = (∏_i (freq_i + JITTER))**(1/n_samples)``

    then normalised so passing joint frequencies sum to 1 (legacy
    ``calcFreq = base / Σ base``). The joint count is ``floor(base / min base)`` so
    that the smallest joint clonotype has count 1 (legacy
    ``calcCount = base / min base``).

    Args:
        samples: A list of per-sample clonotype frames.
        key: Match key name or explicit column tuple (see
            :func:`vdjtools.preprocess.pool.resolve_key`).
        min_samples: Minimum number of samples a clonotype must occur in to be kept.
        names: Optional per-sample names for the ``freq_*`` / ``count_*`` columns;
            defaults to the sample indices ``0..n-1``.

    Returns:
        A joint clonotype frame: the key columns, one ``freq_<name>`` and
        ``count_<name>`` column per sample, ``incidence`` (number of samples
        present), ``frequency`` (normalised geometric-mean joint frequency) and
        ``duplicate_count`` (normalised joint count, smallest = 1), sorted by
        descending ``frequency``.

    Raises:
        ValueError: If ``names`` is given but its length differs from ``samples``.
    """
    key = resolve_key(key)
    n = len(samples)
    if names is None:
        names = [str(i) for i in range(n)]
    elif len(names) != n:
        raise ValueError(f"names has {len(names)} entries but there are {n} samples")

    joined: pl.DataFrame | None = None
    freq_cols, count_cols = [], []
    for name, s in zip(names, samples):
        fc, cc = f"freq_{name}", f"count_{name}"
        freq_cols.append(fc)
        count_cols.append(cc)
        g = (
            normalize(s, recompute_freq=True)
            .group_by(key, maintain_order=True)
            .agg(pl.col(FREQ).sum().alias(fc), pl.col(COUNT).sum().alias(cc))
        )
        joined = g if joined is None else joined.join(g, on=key, how="full", coalesce=True)

    assert joined is not None
    joined = joined.with_columns(
        [pl.col(c).fill_null(0.0) for c in freq_cols]
        + [pl.col(c).fill_null(0).cast(pl.Int64) for c in count_cols]
    )

    incidence = pl.sum_horizontal([(pl.col(c) > 0).cast(pl.Int64) for c in count_cols])
    joined = joined.with_columns(incidence.alias("incidence"))
    joined = joined.filter(pl.col("incidence") >= min_samples)
    if joined.height == 0:
        base = pl.Series("_base", [], dtype=pl.Float64)
        joined = joined.with_columns(base, pl.lit(0.0).alias(FREQ),
                                     pl.lit(0, dtype=pl.Int64).alias(COUNT))
        return joined.select([*key, *freq_cols, *count_cols, "incidence", FREQ, COUNT])

    log_sum = pl.sum_horizontal([(pl.col(c) + JITTER).log() for c in freq_cols])
    joined = joined.with_columns((log_sum / n).exp().alias("_base"))
    total = joined["_base"].sum()
    min_base = joined["_base"].min()
    joined = joined.with_columns(
        (pl.col("_base") / total).alias(FREQ),
        (pl.col("_base") / min_base).floor().cast(pl.Int64).alias(COUNT),
    )
    return (
        joined.select([*key, *freq_cols, *count_cols, "incidence", FREQ, COUNT])
        .sort(FREQ, descending=True, maintain_order=True)
    )
