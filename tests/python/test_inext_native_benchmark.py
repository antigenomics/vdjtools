"""Native-vs-numpy iNEXT benchmark (guarded; off by default).

Run with::

    RUN_BENCHMARK=1 pytest tests/python/test_inext_native_benchmark.py -s

Demonstrates the "many repertoires, quickly" goal on self-contained synthetic
repertoires (deterministic; no network):

1. ``inext_bootstrap`` (native, GIL released) vs the numpy ``_bootstrap_se``
   reference on a ~50k-clone sample (nboot=50) -> per-sample speedup.
2. ``inext_batch`` (native, threaded across samples) vs a sequential numpy loop
   over ~80 samples -> cohort throughput + speedup.

Wall time via ``time.perf_counter``; peak process RSS via ``resource.getrusage``.
"""
import importlib
import os
import resource
import sys
import time

import numpy as np
import pytest

pytestmark = pytest.mark.skipif(not os.getenv("RUN_BENCHMARK"),
                                reason="set RUN_BENCHMARK=1 to run")

import vdjtools._core as core  # noqa: E402
from vdjtools import stats  # noqa: E402

ref = importlib.import_module("vdjtools.stats.inext")


def _peak_rss_mb() -> float:
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return rss / (1024 * 1024) if sys.platform == "darwin" else rss / 1024


def _timed(fn):
    t0 = time.perf_counter()
    out = fn()
    return out, time.perf_counter() - t0


def _synthetic(n_clones: int, seed: int) -> np.ndarray:
    """A heavy-tailed repertoire: singleton-rich geometric abundances (>=1)."""
    rng = np.random.default_rng(seed)
    return rng.geometric(p=0.3, size=n_clones).astype(np.float64)


def test_bootstrap_speedup():
    x = _synthetic(50_000, seed=0)
    qs = [0, 1, 2]
    n = int(x.sum())
    sizes = [float(m) for m in np.unique(
        np.floor(np.linspace(1, 2 * n, 40)).astype(np.int64))]
    nboot = 50

    _, w_py = _timed(lambda: ref._bootstrap_se(x, sizes, qs, nboot, 0))
    _, w_na = _timed(lambda: core.inext_bootstrap(x.tolist(), qs, sizes, nboot, 0))

    print(f"\nbootstrap  clones={x.size:,}  reads={n:,}  nboot={nboot}  "
          f"qs={qs}  sizes={len(sizes)}")
    print(f"  numpy   {w_py:8.3f} s")
    print(f"  native  {w_na:8.3f} s   speedup x{w_py / w_na:6.1f}   "
          f"peak_rss={_peak_rss_mb():.1f} MB")
    assert w_na > 0


def test_batch_throughput():
    n_samples = 80
    samples = [_synthetic(rng_size, seed=i)
               for i, rng_size in enumerate(
                   np.random.default_rng(1).integers(2_000, 8_000, n_samples))]
    qs = (0, 1, 2)
    nboot = 50

    # sequential numpy loop (per-sample inext with numpy bootstrap SE)
    def py_loop():
        return [ref._bootstrap_se(
            x, _sizes(x), [0, 1, 2], nboot, i) for i, x in enumerate(samples)]

    def _sizes(x):
        n = int(x.sum())
        grid = np.unique(np.concatenate([
            np.floor(np.linspace(1, 2 * n, 40)).astype(np.int64),
            np.array([n])]))
        return [float(m) for m in grid if m >= 1]

    _, w_py = _timed(py_loop)
    _, w_na = _timed(lambda: stats.inext_batch(
        samples, q=qs, se=True, nboot=nboot, seed=0, threads=0))

    total_reads = sum(int(x.sum()) for x in samples)
    print(f"\nbatch  samples={n_samples}  total_clones={sum(x.size for x in samples):,}  "
          f"total_reads={total_reads:,}  nboot={nboot}  threads={os.cpu_count()}")
    print(f"  numpy loop  {w_py:8.3f} s   {n_samples / w_py:6.1f} samples/s")
    print(f"  native batch{w_na:8.3f} s   {n_samples / w_na:6.1f} samples/s   "
          f"speedup x{w_py / w_na:6.1f}   peak_rss={_peak_rss_mb():.1f} MB")
    assert w_na > 0
