"""Large-cohort AIRR analytics with vdjtools: Parquet + streaming polars.

The pattern for a cohort too large to hold in RAM (think 100k AIRR repertoires with
a metadata sheet). Three moves:

1. ``ingest_cohort`` streams every per-sample file through the readers ONE AT A TIME
   and writes a hive-partitioned Parquet dataset (``sample_id=<id>/part.parquet``) —
   peak memory is one sample, not the cohort. Run once.
2. ``scan_cohort`` opens the whole dataset as a single ``polars.LazyFrame``
   (``sample_id`` recovered from the path, metadata joined lazily).
3. Every analysis is a ``group_by(...).agg(...)`` collected with
   ``engine="streaming"`` — the cohort feature matrix in one pass, the cohort never
   materialised; predicate/projection pushdown means a filtered query touches only
   the partitions it needs.

Run it::

    python examples/scale_cohort.py                  # 200-sample demo (~seconds)
    python examples/scale_cohort.py --samples 5000   # watch peak RSS stay flat

Bump ``--samples`` toward 100000 for the real thing: ingest RSS stays flat (one
sample at a time), and the streamed cohort usage matrix collects without ever
holding the cohort in memory.

Single-cell is the OTHER shape: there ``obs`` is one row per CELL, so the right
container is AnnData / MuData, not a per-sample parquet cohort — see
``vdjtools.sc.to_anndata`` (the ``[sc]`` extra). Bulk per-sample cohorts must NOT go
in AnnData: obs=clonotype makes an ~1e9 × 100k almost-empty sparse X. Rule:
single-cell (obs=cell) → AnnData; bulk cohort (per-sample tables) → scan_cohort.
"""
from __future__ import annotations

import argparse
import resource
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl

from vdjtools import io as vio
from vdjtools.io import schema as S

_V = ["TRBV5-1", "TRBV7-9", "TRBV20-1", "TRBV28", "TRBV19", "TRBV6-5", "TRBV12-3"]
_J = ["TRBJ2-1", "TRBJ2-7", "TRBJ1-1", "TRBJ2-3", "TRBJ1-2"]
_AA = np.array(list("ACDEFGHIKLMNPQRSTVWY"))


def _peak_mb() -> float:
    """Peak resident set size in MB (ru_maxrss is bytes on macOS, KiB on Linux)."""
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return rss / 1e6 if sys.platform == "darwin" else rss / 1e3


def synthesize(base: Path, n_samples: int, n_clones: int, seed: int = 0) -> pl.DataFrame:
    """Write ``n_samples`` synthetic AIRR TSVs to ``base`` and return the metadata."""
    base.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    ages, groups = [], []
    for i in range(n_samples):
        k = rng.integers(n_clones // 2, n_clones + 1)
        lens = rng.integers(10, 18, size=k)
        cdr3 = ["C" + "".join(_AA[rng.integers(0, 20, size=int(le - 2))]) + "F" for le in lens]
        pl.DataFrame({
            "v_call": rng.choice(_V, size=k), "j_call": rng.choice(_J, size=k),
            "junction_aa": cdr3,
            "duplicate_count": rng.integers(1, 500, size=k),
        }).write_csv(base / f"S{i:05d}.tsv", separator="\t")
        age = int(rng.integers(1, 90))
        ages.append(str(age))
        groups.append("cord" if age < 2 else "adult" if age < 60 else "elderly")
    return pl.DataFrame({"sample_name": [f"S{i:05d}" for i in range(n_samples)],
                         "age": ages, "group": groups})


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--samples", type=int, default=200)
    ap.add_argument("--clonotypes", type=int, default=500, help="~clonotypes per sample")
    ap.add_argument("--workdir", type=Path, default=Path("examples/.data/scale"))
    args = ap.parse_args()

    raw, cohort = args.workdir / "raw", args.workdir / "cohort"
    print(f"# synthesizing {args.samples} samples × ~{args.clonotypes} clonotypes …")
    meta = synthesize(raw, args.samples, args.clonotypes)

    # 1. Streaming ingest → hive-partitioned Parquet. Peak RSS ≈ one sample.
    t0 = time.perf_counter()
    vio.ingest_cohort(meta, raw, cohort, file_template="{sample}.tsv")
    print(f"[1] ingest_cohort  {time.perf_counter()-t0:5.2f}s   peak RSS {_peak_mb():6.0f} MB "
          f"(flat vs #samples — one sample in RAM at a time)")

    # 2. One LazyFrame over the whole cohort; nothing read yet.
    lf = vio.scan_cohort(cohort)

    # 3a. Cohort V-usage matrix in ONE streamed pass (cohort never materialised).
    t0 = time.perf_counter()
    usage = (lf.group_by(["sample_id", S.V_CALL])
             .agg(pl.col(S.COUNT).sum().alias("n"))
             .collect(engine="streaming"))
    wide = usage.pivot(values="n", index="sample_id", on=S.V_CALL).fill_null(0)
    print(f"[2] streamed V-usage matrix  {time.perf_counter()-t0:5.2f}s   "
          f"shape {wide.shape}   peak RSS {_peak_mb():6.0f} MB")

    # 3b. Streamed per-sample richness + reads, metadata joined lazily.
    summ = (lf.group_by("sample_id")
            .agg(pl.len().alias("clonotypes"), pl.col(S.COUNT).sum().alias("reads"),
                 pl.col("group").first(), pl.col("age").first().cast(pl.Int32))
            .sort("age").collect(engine="streaming"))
    by_group = (summ.group_by("group").agg(pl.col("clonotypes").mean().round(1))
                .sort("group"))
    print("[3] per-sample summary joined to metadata; mean richness by age group:")
    for r in by_group.iter_rows(named=True):
        print(f"      {r['group']:8s} {r['clonotypes']:.1f}")

    # 3c. Predicate pushdown: an age-filtered query scans only matching partitions.
    elderly = (lf.filter(pl.col("group") == "elderly")
               .select(pl.col("sample_id").n_unique()).collect().item())
    print(f"[4] pushdown-filtered scan: {elderly} elderly samples "
          f"(only their partitions read)")

    # 4. Single-file Parquet read (typed, no all-Utf8 pass).
    one = next(iter(cohort.glob("sample_id=*/part.parquet")))
    df = vio.read_parquet(one)
    assert df.columns == S.COLUMNS + [S.LOCUS] and (df[S.FREQ].sum() - 1.0) < 1e-9
    print(f"[5] read_parquet single file: {df.height} clonotypes, canonical schema ✓")

    # runnable check: streamed totals equal a direct per-sample sum (no data lost).
    direct = sum(vio.read(p, fmt="airr")[S.COUNT].sum()
                 for p in sorted(raw.glob("*.tsv"))[:5])
    streamed = (lf.filter(pl.col("sample_id").is_in(
        [f"S{i:05d}" for i in range(5)])).select(pl.col(S.COUNT).sum())
        .collect().item())
    assert direct == streamed, (direct, streamed)
    print("\nOK — streamed cohort totals match per-sample sums; cohort never materialised.")


if __name__ == "__main__":
    main()
