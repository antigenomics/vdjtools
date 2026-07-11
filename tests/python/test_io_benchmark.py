"""Read + analytics throughput benchmarks (guarded; off by default).

Run with::

    RUN_BENCHMARK=1 pytest tests/python/test_io_benchmark.py -s

Data is fetched from HuggingFace (skips cleanly offline). Wall time via
``time.perf_counter``; peak process RSS via ``resource.getrusage``. Read targets
are the supported formats only: ankspond ``new/`` (AIRR-hybrid) and the native
vdjtools control export (row-capped). The ankspond ``old/`` files are a legacy
MiGEC dotted format outside the basic readers' scope.
"""
import os
import resource
import sys
import time

import pytest

pytestmark = pytest.mark.skipif(not os.getenv("RUN_BENCHMARK"),
                                reason="set RUN_BENCHMARK=1 to run")

from vdjtools import features, io as vio, stats  # noqa: E402

ANKSPOND = "isalgo/airr_ankspond"
CONTROL = "isalgo/airr_control"
CONTROL_CAP = 500_000


def _peak_rss_mb() -> float:
    """Peak resident set size of this process in MB (ru_maxrss unit differs by OS)."""
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return rss / (1024 * 1024) if sys.platform == "darwin" else rss / 1024


def _timed(fn):
    """Run ``fn`` and return ``(result, wall_seconds)``."""
    t0 = time.perf_counter()
    out = fn()
    return out, time.perf_counter() - t0


def test_io_and_analytics_benchmark(hf):
    table = []  # (dataset, n_rows, op, wall_s, rows_per_s, peak_rss_mb)

    # --- read tier 1: ankspond new/ (AIRR-hybrid, small) ---
    p1 = hf(ANKSPOND, "new/Azh_0.tsv.gz")
    df1, w1 = _timed(lambda: vio.read(p1))
    table.append(("ankspond/Azh_0", df1.height, "read", w1, df1.height / w1, _peak_rss_mb()))

    # --- read tier 2: native control export, row-capped for scale ---
    hub = pytest.importorskip("huggingface_hub")
    files = hub.list_repo_files(repo_id=CONTROL, repo_type="dataset")
    tsv = next((f for f in files if f.endswith(".vdjtools.tsv.gz")), None)
    if tsv is None:
        pytest.skip(f"no vdjtools tsv in {CONTROL}")
    p2 = hf(CONTROL, tsv)
    df2, w2 = _timed(lambda: vio.read(p2, n_rows=CONTROL_CAP))
    ds2 = f"control[{CONTROL_CAP}]"
    table.append((ds2, df2.height, "read", w2, df2.height / w2, _peak_rss_mb()))

    # --- analytics throughput on the loaded control frame ---
    ops = [
        ("diversity_stats", lambda: stats.diversity_stats(df2)),
        ("spectratype_aa", lambda: stats.spectratype(df2, kind="aa")),
        ("segment_usage_v", lambda: stats.segment_usage(df2, "v")),
        ("segment_usage_j", lambda: stats.segment_usage(df2, "j")),
        ("kmer_profile_k3", lambda: features.kmer_profile(df2, k=3)),
    ]
    for name, fn in ops:
        _, w = _timed(fn)
        table.append((ds2, df2.height, name, w, df2.height / w, _peak_rss_mb()))

    # --- structured summary table ---
    hdr = (f"{'dataset':<20} {'n_rows':>9} {'op':<18} {'wall_s':>8} "
           f"{'rows_per_s':>13} {'peak_rss_mb':>12}")
    print("\n" + hdr)
    print("-" * len(hdr))
    for ds, n, op, w, rps, rss in table:
        print(f"{ds:<20} {n:>9} {op:<18} {w:>8.3f} {rps:>13,.0f} {rss:>12.1f}")

    assert df2.height == CONTROL_CAP
    assert df1.height > 0
