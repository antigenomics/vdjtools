"""Run a per-sample signature function across a cohort, in processes.

One home for the pool, because three callers need it: ``vsig_cohort`` here, ``rsig_cohort`` and
``signature_cohort`` in mirpy. Each emits its own half; none of them should carry its own copy of
the process handling.

``spawn``, not ``fork``: polars documents that it cannot be combined with ``fork`` (the child
inherits held mutexes and file locks and hangs on the first one it touches), CPython 3.12 warns
about ``fork()`` in a multi-threaded process, and 3.14 stops defaulting to it. Because ``spawn``
pays a fresh interpreter per worker, the cohort is cut into exactly as many contiguous slices as
there are workers -- one task each, so that cost is paid once per worker and not once per sample.
"""
from __future__ import annotations

__all__ = ["slices", "parallel_rows"]


def slices(n: int, workers: int) -> list[tuple[int, int]]:
    """``workers`` contiguous half-open ranges covering ``0..n``, sizes differing by at most one."""
    return [(round(i * n / workers), round((i + 1) * n / workers)) for i in range(workers)]


def _chunk(args):
    """One contiguous slice, start to finish, in one worker. ``fn`` must be picklable."""
    items, fn = args
    return [fn(it) for it in items]


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

    items = list(items)
    workers = min(len(items), n_jobs if n_jobs > 0 else available_cores())
    if n_jobs == 1 or workers < 2:
        return [fn(it) for it in items]
    chunks = [(items[a:b], fn) for a, b in slices(len(items), workers) if b > a]
    try:
        with ProcessPoolExecutor(max_workers=len(chunks), mp_context=get_context("spawn")) as ex:
            return [row for part in ex.map(_chunk, chunks) for row in part]
    except (BrokenExecutor, RuntimeError) as e:
        raise RuntimeError(
            f"could not start worker processes for n_jobs={n_jobs} ({type(e).__name__}). Workers "
            "are spawned, not forked -- polars cannot be combined with fork -- and a spawned "
            "worker re-imports the module that called this. That works from an importable module "
            "and fails from `python - <<EOF`, `python -c`, or a call at the top level of a "
            "script. Guard the call with `if __name__ == '__main__':`, or pass n_jobs=1. It is "
            "NOT falling back to serial: that is what hid a 20x slowdown here before."
        ) from e
