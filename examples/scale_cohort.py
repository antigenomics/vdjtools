# vdjtools — large-cohort AIRR analytics: Parquet + streaming polars, flat memory.
# Reactive marimo app: synthesize N AIRR repertoires, ingest them ONE AT A TIME into a
# hive-partitioned Parquet dataset (peak RSS ≈ one sample, not the cohort), then run every
# analysis as a streamed group_by over a single scan_cohort LazyFrame. Slide the sample count
# up and watch peak RSS stay flat. Run with:
#     marimo edit examples/scale_cohort.py
import marimo

__generated_with = "0.23.14"
app = marimo.App(width="medium")


@app.cell
def _(mo):
    mo.md(
        r"""
        # Large-cohort analytics — Parquet + streaming polars

        The pattern for a cohort too large to hold in RAM (think 100k AIRR repertoires with a
        metadata sheet). Three moves, each demonstrated below on a synthetic cohort you size with
        the slider:

        1. **`io.ingest_cohort`** streams every per-sample file through the readers *one at a time*
           and writes a hive-partitioned Parquet dataset (`sample_id=<id>/part.parquet`) — peak
           memory is one sample, not the cohort.
        2. **`io.scan_cohort`** opens the whole dataset as a single `polars.LazyFrame`
           (`sample_id` recovered from the path, metadata joined lazily).
        3. Every analysis is a `group_by(...).agg(...)` collected with `engine="streaming"` — the
           cohort feature matrix in one pass, the cohort never materialised; predicate/projection
           pushdown means a filtered query touches only the partitions it needs.

        Slide **samples** up (toward thousands) and watch the reported peak RSS stay flat.

        > Bulk per-sample cohorts go in **Parquet**, not AnnData: `obs=clonotype` would make an
        > ~1e9 × 100k almost-empty sparse `X`. Single-cell (`obs=cell`) is the AnnData shape — see
        > `vdjtools.sc.to_anndata`. Rule: bulk cohort → `scan_cohort`; single-cell → AnnData.
        """
    )
    return


@app.cell
def _():
    # --- imports & helpers (single cell so every name is defined once) ---
    import resource
    import sys
    import time
    from pathlib import Path

    import marimo as mo
    import numpy as np
    import polars as pl

    from vdjtools import io as vio
    from vdjtools.io import schema as S

    _V = ["TRBV5-1", "TRBV7-9", "TRBV20-1", "TRBV28", "TRBV19", "TRBV6-5", "TRBV12-3"]
    _J = ["TRBJ2-1", "TRBJ2-7", "TRBJ1-1", "TRBJ2-3", "TRBJ1-2"]
    _AA = np.array(list("ACDEFGHIKLMNPQRSTVWY"))

    def peak_mb():
        """Peak resident set (MB); ru_maxrss is bytes on macOS, KiB on Linux. Process-cumulative."""
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return rss / 1e6 if sys.platform == "darwin" else rss / 1e3

    def synthesize(base, n_samples, n_clones, seed=0):
        """Write ``n_samples`` synthetic AIRR TSVs to ``base``; return the metadata frame."""
        base.mkdir(parents=True, exist_ok=True)
        rng = np.random.default_rng(seed)
        ages, groups = [], []
        for i in range(n_samples):
            k = int(rng.integers(n_clones // 2, n_clones + 1))
            lens = rng.integers(10, 18, size=k)
            cdr3 = ["C" + "".join(_AA[rng.integers(0, 20, size=int(le - 2))]) + "F" for le in lens]
            pl.DataFrame({
                "v_call": rng.choice(_V, size=k), "j_call": rng.choice(_J, size=k),
                "junction_aa": cdr3, "duplicate_count": rng.integers(1, 500, size=k),
            }).write_csv(base / f"S{i:05d}.tsv", separator="\t")
            age = int(rng.integers(1, 90))
            ages.append(str(age))
            groups.append("cord" if age < 2 else "adult" if age < 60 else "elderly")
        return pl.DataFrame({"sample_name": [f"S{i:05d}" for i in range(n_samples)],
                             "age": ages, "group": groups})

    return Path, S, mo, peak_mb, pl, synthesize, time, vio


@app.cell
def _(mo):
    n_samples = mo.ui.slider(50, 5000, value=200, step=50, label="samples (repertoires)")
    n_clones = mo.ui.slider(100, 1000, value=500, step=100, label="~clonotypes / sample")
    mo.hstack([n_samples, n_clones], justify="start", gap=1.5)
    return n_clones, n_samples


@app.cell
def _(Path, mo, n_clones, n_samples, peak_mb, synthesize, time, vio):
    _work = (mo.notebook_dir() or Path.cwd()) / ".data" / "scale"
    raw, cohort = _work / f"raw_{n_samples.value}", _work / f"cohort_{n_samples.value}"

    _t0 = time.perf_counter()
    meta = synthesize(raw, n_samples.value, n_clones.value)
    _t_syn = time.perf_counter() - _t0

    # 1. Streaming ingest → hive-partitioned Parquet. Peak RSS ≈ one sample.
    _t0 = time.perf_counter()
    if not any(cohort.glob("sample_id=*/*.parquet")):
        vio.ingest_cohort(meta, raw, cohort, file_template="{sample}.tsv")
    _t_ing = time.perf_counter() - _t0
    lf = vio.scan_cohort(cohort)                           # one LazyFrame; nothing read yet

    mo.md(
        f"Synthesised **{n_samples.value} samples × ~{n_clones.value} clonotypes** in "
        f"{_t_syn:.2f}s; `ingest_cohort` → hive Parquet in **{_t_ing:.2f}s**, "
        f"**peak RSS {peak_mb():.0f} MB** (flat vs #samples — one sample in RAM at a time). "
        f"Cohort at `{cohort.name}`."
    )
    return cohort, lf, meta


@app.cell
def _(S, lf, mo, peak_mb, pl, time):
    # Cohort V-usage matrix in ONE streamed pass (cohort never materialised).
    _t0 = time.perf_counter()
    usage = (lf.group_by(["sample_id", S.V_CALL]).agg(pl.col(S.COUNT).sum().alias("n"))
             .collect(engine="streaming"))
    wide = usage.pivot(values="n", index="sample_id", on=S.V_CALL).fill_null(0)
    mo.md(f"**Streamed V-usage matrix** in {time.perf_counter()-_t0:.2f}s — shape "
          f"`{wide.shape}`, peak RSS {peak_mb():.0f} MB. The cohort was never held in memory.")
    wide.head(6)
    return


@app.cell
def _(S, lf, mo, pl, time):
    # Streamed per-sample richness + reads, metadata joined lazily; mean richness by age group.
    _t0 = time.perf_counter()
    summ = (lf.group_by("sample_id")
            .agg(pl.len().alias("clonotypes"), pl.col(S.COUNT).sum().alias("reads"),
                 pl.col("group").first(), pl.col("age").first().cast(pl.Int32))
            .sort("age").collect(engine="streaming"))
    by_group = (summ.group_by("group").agg(pl.col("clonotypes").mean().round(1).alias("mean_clonotypes"),
                                           pl.len().alias("n_samples")).sort("group"))
    mo.md(f"**Per-sample summary** joined to metadata in {time.perf_counter()-_t0:.2f}s; "
          f"mean richness by age group:")
    by_group
    return


@app.cell
def _(lf, mo, pl):
    # Predicate pushdown: an age-filtered query scans only matching partitions.
    elderly = (lf.filter(pl.col("group") == "elderly")
               .select(pl.col("sample_id").n_unique()).collect().item())
    mo.md(f"**Pushdown-filtered scan:** {elderly} elderly samples — only their partitions were read "
          f"(predicate/projection pushdown, not a full-cohort scan).")
    return


@app.cell
def _(S, cohort, mo, vio):
    # A single Parquet partition reads back as the canonical typed schema (no all-Utf8 pass).
    _first = next(cohort.glob("sample_id=*/part.parquet"))
    _df = vio.read_parquet(_first)
    ok = _df.columns == S.COLUMNS + [S.LOCUS] and abs(_df[S.FREQ].sum() - 1.0) < 1e-6
    mo.md(f"**Single-partition `read_parquet`** → canonical schema round-trips ✓ = **{ok}** "
          f"(typed columns, `frequency` sums to 1 for that sample).")
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ---
        **Takeaway.** Ingest RSS stays flat as you slide the sample count up (one sample in RAM at
        a time), and every downstream statistic — the V-usage matrix, per-sample richness, the
        age-filtered count — collects with `engine="streaming"` without ever materialising the
        cohort. This is the shape for 100k-repertoire cohorts: persist once with `ingest_cohort`,
        then query a `scan_cohort` LazyFrame.
        """
    )
    return


if __name__ == "__main__":
    app.run()
