"""Cohort summary-stat throughput: serial vs parallel map_samples vs fused --cohort.

Guarded; off by default. Run with::

    RUN_BENCHMARK=1 pytest tests/python/test_cohort_benchmark.py -s

Synthesises a small cohort on disk (no network), then times three cohort-diversity
strategies and reports wall time + peak process RSS:

* ``serial``     — :func:`map_samples` with one worker (the pre-change behaviour).
* ``parallel``   — :func:`map_samples` across all cores (streams + threads).
* ``cohort``     — :func:`diversity_cohort` over a hive-partitioned parquet scan
                   (one streamed pass; peak memory ~= the count spectrum).

Wall time via ``time.perf_counter``; peak RSS via ``resource.getrusage``. Small by
design so it stays inside the fast-benchmark budget on a laptop; scale ``N_SAMPLES`` /
``N_CLONES`` up only on a bigger box (heavy runs belong on the cluster, per CLAUDE.md).
"""
import gzip
import os
import resource
import sys
import time

import numpy as np
import polars as pl
import pytest

pytestmark = pytest.mark.skipif(not os.getenv("RUN_BENCHMARK"),
                                reason="set RUN_BENCHMARK=1 to run")

from vdjtools.io import ingest_cohort, map_samples, read_metadata, scan_cohort  # noqa: E402
from vdjtools.stats import diversity_cohort, diversity_stats  # noqa: E402

N_SAMPLES = 24
N_CLONES = 8_000
_VS = [f"TRBV{i}-1*01" for i in range(1, 30)]
_JS = [f"TRBJ{i}-1*01" for i in range(1, 7)]


def _peak_rss_mb() -> float:
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return rss / (1024 * 1024) if sys.platform == "darwin" else rss / 1024


def _timed(fn):
    t0 = time.perf_counter()
    out = fn()
    return out, time.perf_counter() - t0


def _write_cohort(base):
    rng = np.random.default_rng(0)
    for s in range(N_SAMPLES):
        counts = rng.integers(1, 500, N_CLONES)
        vs = rng.choice(_VS, N_CLONES)
        js = rng.choice(_JS, N_CLONES)
        lens = rng.integers(10, 18, N_CLONES)
        with gzip.open(base / f"s{s}.tsv.gz", "wt") as f:
            f.write("count\tfreq\tcdr3nt\tcdr3aa\tv\td\tj\n")
            f.writelines(f"{c}\t0\t{'A'*3}\t{'C'*L}\t{v}\t.\t{j}\n"
                         for c, v, j, L in zip(counts, vs, js, lens))
    md = base / "metadata.tsv"
    pl.DataFrame({"sample_name": [f"s{s}" for s in range(N_SAMPLES)]}).write_csv(md, separator="\t")
    return md


def test_cohort_diversity_benchmark(tmp_path):
    md = _write_cohort(tmp_path)
    items = [(f"s{s}", tmp_path / f"s{s}.tsv.gz") for s in range(N_SAMPLES)]
    cohort_dir = tmp_path / "cohort"
    ingest_cohort(read_metadata(md), tmp_path, cohort_dir,
                  sample_col="sample_name", file_template="{sample}.tsv.gz")

    _, serial_w = _timed(lambda: map_samples(diversity_stats, items, workers=1))
    par, par_w = _timed(lambda: map_samples(diversity_stats, items, workers=None))
    coh, coh_w = _timed(lambda: diversity_cohort(scan_cohort(cohort_dir, join_metadata=False)))
    rss = _peak_rss_mb()

    # all three agree on the diversity table (parallel vs serial identical; cohort exact)
    par_df = pl.concat([r.select(pl.lit(sid).alias("sample_id"), pl.all()) for sid, r in par],
                       how="vertical_relaxed").sort("sample_id")
    assert par_df.equals(coh.sort("sample_id"))

    rows = N_SAMPLES * N_CLONES
    print(f"\ncohort: {N_SAMPLES} samples x {N_CLONES} clones = {rows:,} rows; "
          f"{os.cpu_count()} cores; peak RSS {rss:.0f} MB")
    hdr = f"{'strategy':<12}{'wall_s':>9}{'speedup':>9}"
    print(hdr + "\n" + "-" * len(hdr))
    for name, w in [("serial", serial_w), ("parallel", par_w), ("cohort", coh_w)]:
        print(f"{name:<12}{w:>9.3f}{serial_w / w:>8.1f}x")
