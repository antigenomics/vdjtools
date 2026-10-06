"""Process samples independently, with one numerical-kernel thread per worker.

Workers load the callable and its frozen resources once. Each task reads and reduces one sample,
then returns only its feature row. Ordered map preserves input order while workers take the next
sample as soon as they finish, so unequal repertoire depths do not strand idle workers.
"""
from __future__ import annotations

from contextlib import contextmanager

from ..cores import _KERNEL_ENV

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

    NOTE: at ``n_jobs != 1`` the callable is **pickled**, so it must be picklable: a
    ``functools.partial`` over a module-level function works, a lambda or a closure does not
    (``AttributeError: Can't get local object ...``). That is why the CLIs build theirs with
    ``functools.partial(_read_sample, paths)`` rather than a lambda.
    """
    import polars as pl

    if callable(sample) and not isinstance(sample, (dict, pl.DataFrame, pl.LazyFrame)):
        return sample()
    return sample


def slices(n: int, workers: int) -> list[tuple[int, int]]:
    """``workers`` contiguous half-open ranges covering ``0..n``, sizes differing by at most one."""
    return [(round(i * n / workers), round((i + 1) * n / workers)) for i in range(workers)]


#: Set in every spawned worker. A pool worker must not also ask a threaded kernel for the whole
#: machine: at ``n_jobs=cores`` that is cores x cores threads, and measured 2026-09-26 on mirpy's
#: rsig it turned a 1.47x best case into 1.03x for 4.7x the CPU. Kernels that size themselves off
#: the core count read this and take one thread instead. An env var rather than an argument
#: because the reader is several frames down inside a third-party call.
WORKER_ENV = "VDJTOOLS_POOL_WORKER"


def in_pool_worker() -> bool:
    """True inside a :func:`parallel_rows` worker process. Ask before claiming every core."""
    import os

    return bool(os.environ.get(WORKER_ENV))


@contextmanager
def single_threaded_workers():
    """Temporarily give spawned children one thread per kernel; restore the caller's env."""
    import os

    previous = {key: os.environ.get(key) for key in (*_KERNEL_ENV, WORKER_ENV)}
    try:
        os.environ.update(dict.fromkeys(previous, "1"))
        yield
    finally:
        for key, value in previous.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


_WORKER_FN = None


def _init_worker(fn):
    global _WORKER_FN
    _WORKER_FN = fn


def _row(item):
    return _WORKER_FN(item)


def parallel_rows(items, fn, n_jobs: int) -> list[dict]:
    """``[fn(item) for item in items]``, across ``n_jobs`` processes. ``1`` stays in-process.

    Args:
        items: Per-sample inputs, passed to ``fn`` one at a time.
        fn: A **picklable** callable (a module-level function, or ``functools.partial`` over one --
            not a lambda) mapping one item to one output row.
        n_jobs: Worker processes; ``0`` means every core this process may use
            (:func:`vdjtools.cores.available_cores`), ``1`` means no pool at all.

    Raises:
        RuntimeError: If the workers cannot start. It does **not** fall back to one process: a
            fallback that keeps the answer correct is how a 20x slowdown stayed invisible here.
    """
    from concurrent.futures import BrokenExecutor, ProcessPoolExecutor
    from multiprocessing import get_context

    from ..cores import available_cores

    if n_jobs < 0:
        raise ValueError("n_jobs must be non-negative (0 means all available cores)")
    items = list(items)
    workers = min(len(items), n_jobs if n_jobs > 0 else available_cores())
    if n_jobs == 1 or workers < 2:
        return [fn(it) for it in items]
    with single_threaded_workers():
        try:
            executor = ProcessPoolExecutor(max_workers=workers, mp_context=get_context("spawn"),
                                           initializer=_init_worker, initargs=(fn,))
        except (BrokenExecutor, RuntimeError) as e:
            raise _pool_error(n_jobs, e) from e
        try:
            with executor:
                return list(executor.map(_row, items, chunksize=1))
        except BrokenExecutor as e:
            raise _pool_error(n_jobs, e) from e


def _pool_error(n_jobs, error):
    return RuntimeError(
        f"worker processes failed for n_jobs={n_jobs} ({type(error).__name__}). "
        "Use an importable module guarded by if __name__ == '__main__', or pass n_jobs=1. "
        "Workers use spawn because Polars cannot safely use fork. NOT falling back to serial.")
