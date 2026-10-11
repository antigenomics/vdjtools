"""Process samples independently, with one numerical-kernel thread per worker.

Threads share the callable and its frozen resources. Each task reads and reduces one sample,
then returns only its feature row. Ordered map preserves input order while workers take the next
sample as soon as they finish, so unequal repertoire depths do not strand idle workers.
"""
from __future__ import annotations

__all__ = ["slices", "parallel_rows", "resolve_sample"]


def resolve_sample(sample):
    """Call a deferred sample, so the read happens in the worker rather than the parent.

    Every ``*_cohort`` accepts a zero-argument callable in place of a sample, which is what keeps
    peak memory at ``O(n_jobs)`` samples instead of the whole cohort. Defined once here because
    all three cohort entry points need it and two of them used to forget: a deferred sample
    reached ``_locus_frames`` unresolved and died on ``'function' object has no attribute
    'columns'``, which made the documented low-memory path unusable rather than merely slow.

    A ``dict`` or ``DataFrame`` is never treated as deferred even if it were somehow callable --
    the sample types come first, so this cannot misfire on a real sample.

    Deferred readers may be closures or lambdas; no pickling is required.
    """
    import polars as pl

    if callable(sample) and not isinstance(sample, (dict, pl.DataFrame, pl.LazyFrame)):
        return sample()
    return sample


def slices(n: int, workers: int) -> list[tuple[int, int]]:
    """``workers`` contiguous half-open ranges covering ``0..n``, sizes differing by at most one."""
    return [(round(i * n / workers), round((i + 1) * n / workers)) for i in range(workers)]


def parallel_rows(items, fn, n_jobs: int) -> list[dict]:
    """Ordered sample map across bounded threads; ``1`` runs without a pool.

    ``0`` uses the allocated cores. Resources are shared, and inner numerical kernels
    use one thread per worker. Callbacks must not mutate shared inputs or resources.
    Worker startup and callback failures propagate; there is no serial fallback.
    """
    from concurrent.futures import ThreadPoolExecutor

    from ..cores import kernel_threads, require_serial_pools, serial_kernels

    try:
        budget = kernel_threads(n_jobs)
    except (TypeError, ValueError) as e:
        raise ValueError(f"invalid n_jobs: {e}") from e
    items = list(items)
    workers = min(len(items), budget)
    def row(item):
        with serial_kernels():
            return fn(item)

    if workers < 2:
        return [row(it) for it in items]
    require_serial_pools()

    with ThreadPoolExecutor(max_workers=workers) as executor:
        return list(executor.map(row, items))
