"""iNEXT throughput benchmark on a real repertoire (guarded; off by default).

Run with::

    RUN_BENCHMARK=1 pytest tests/python/test_inext_benchmark.py -s

Loads the HuggingFace ``isalgo/airr_control`` native vdjtools export (row-capped),
then times ``inext(q=(0, 1, 2))``. The point-estimate curve (``se=False``)
documents the ``O(n)``-per-size interpolation cost of the MVUE frequency-count
engine; a small bootstrap run documents the added CI cost. Skips cleanly offline.
Wall time via ``time.perf_counter``; peak process RSS via ``resource.getrusage``.
"""
import os
import resource
import sys
import time

import pytest

pytestmark = pytest.mark.skipif(not os.getenv("RUN_BENCHMARK"),
                                reason="set RUN_BENCHMARK=1 to run")

from vdjtools import io as vio, stats  # noqa: E402

CONTROL = "isalgo/airr_control"
CONTROL_CAP = 200_000


def _peak_rss_mb() -> float:
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return rss / (1024 * 1024) if sys.platform == "darwin" else rss / 1024


def _timed(fn):
    t0 = time.perf_counter()
    out = fn()
    return out, time.perf_counter() - t0


def test_inext_benchmark(hf):
    hub = pytest.importorskip("huggingface_hub")
    files = hub.list_repo_files(repo_id=CONTROL, repo_type="dataset")
    tsv = next((f for f in files if f.endswith(".vdjtools.tsv.gz")), None)
    if tsv is None:
        pytest.skip(f"no vdjtools tsv in {CONTROL}")
    path = hf(CONTROL, tsv)
    df = vio.read(path, n_rows=CONTROL_CAP)
    counts = df["duplicate_count"].to_numpy()
    n_clones = int((counts > 0).sum())
    n_reads = int(counts.sum())

    table = []  # (op, wall_s, peak_rss_mb)
    ops = [
        ("inext q=(0,1,2) knots=40 se=False", lambda: stats.inext(df, q=(0, 1, 2), se=False)),
        ("inext q=0 knots=40 se=False", lambda: stats.inext(df, q=0, se=False)),
        ("inext q=(0,1,2) nboot=20 se=True", lambda: stats.inext(df, q=(0, 1, 2), nboot=20, se=True)),
    ]
    for name, fn in ops:
        _, w = _timed(fn)
        table.append((name, w, _peak_rss_mb()))

    hdr = f"{'op':<38} {'wall_s':>9} {'peak_rss_mb':>12}"
    print(f"\ncontrol[{CONTROL_CAP}]  clonotypes={n_clones:,}  reads={n_reads:,}")
    print(hdr)
    print("-" * len(hdr))
    for op, w, rss in table:
        print(f"{op:<38} {w:>9.3f} {rss:>12.1f}")

    assert n_clones > 0
