"""Track clonotype frequencies across an ordered set of samples.

Reimplements the legacy ``TrackClonotypes``: given several samples (e.g. a time
course or an age series) it builds one row per clonotype and one frequency column per
sample, so a clonotype's trajectory can be read across the columns. Clonotypes
present in at least one sample are kept, sorted by their summed frequency (most
persistent/abundant first), optionally truncated to the top ``N``.
"""
from __future__ import annotations

import polars as pl

from ..io.schema import JUNCTION_AA, COUNT, J_CALL, V_CALL

#: Default clonotype match key (CDR3 aa + V + J).
DEFAULT_KEY = (JUNCTION_AA, V_CALL, J_CALL)


def _freq_by_key(df: pl.DataFrame, key: list[str], col: str) -> pl.DataFrame:
    """Collapse to unique ``key`` and return within-sample frequency as ``col``."""
    agg = df.group_by(key, maintain_order=True).agg(pl.col(COUNT).sum().alias("_c"))
    total = agg["_c"].sum() or 1
    return agg.with_columns((pl.col("_c") / total).alias(col)).drop("_c")


def track_clonotypes(samples, order=None, top: int | None = None,
                     key=DEFAULT_KEY) -> pl.DataFrame:
    """Pivot per-sample clonotype frequency into one column per sample.

    Args:
        samples: A ``dict`` mapping sample name to clonotype frame, or a ``list`` of
            frames (named ``"0".."N-1"``).
        order: Sample names giving the left-to-right column order. Defaults to the
            samples' natural order (dict insertion / list order). Names not present
            in ``samples`` are ignored.
        top: If given, keep only the ``top`` clonotypes by summed frequency.
        key: Clonotype match key (default CDR3 aa + V + J).

    Returns:
        A ``pl.DataFrame`` with the ``key`` columns, one ``freq_<sample>`` column per
        sample in ``order`` (``0.0`` where the clonotype is absent), and a
        ``freq_sum`` column, sorted by ``freq_sum`` descending. Clonotypes present in
        at least one sample are included.
    """
    key = list(key)
    if isinstance(samples, dict):
        table = dict(samples)
        default_order = list(samples.keys())
    else:
        table = {str(i): df for i, df in enumerate(samples)}
        default_order = list(table.keys())
    order = [n for n in (order if order is not None else default_order) if n in table]

    out: pl.DataFrame | None = None
    freq_cols = []
    for name in order:
        col = f"freq_{name}"
        freq_cols.append(col)
        part = _freq_by_key(table[name], key, col)
        out = part if out is None else out.join(part, on=key, how="full", coalesce=True)

    if out is None:
        return pl.DataFrame(schema={k: pl.Utf8 for k in key})

    out = out.with_columns([pl.col(c).fill_null(0.0) for c in freq_cols])
    out = out.with_columns(pl.sum_horizontal(freq_cols).alias("freq_sum"))
    out = out.sort("freq_sum", descending=True).select([*key, *freq_cols, "freq_sum"])
    return out.head(top) if top is not None else out
