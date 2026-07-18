# vdjtools — how the human TCR repertoire ages, via the cohort-streaming stats.
# Reactive marimo app over the Britanova "Cord Blood to Centenarians" TRB cohort
# (`isalgo/airr_benchmark`, folder `vdjtools/`, 78 donors aged 0-103). One hive-partitioned
# Parquet scan drives every summary: diversity (diversity_cohort), the clone-size distribution
# (singleton -> hyperexpanded), CDR3 spectratype and V-gene usage, each computed as a single
# streamed group_by over the whole cohort. Run with:
#     marimo edit notebooks/aging.py
import marimo

__generated_with = "0.23.14"
app = marimo.App(width="medium")


@app.cell
def _(mo):
    mo.md(
        r"""
        # Ageing of the human TCR repertoire — a vdjtools v2 walkthrough

        As we age the naive T-cell pool contracts and a handful of clones expand, so the
        repertoire becomes **less diverse and more clonal**. This notebook reads all of that off
        the full-depth **Britanova "Cord Blood to Centenarians" TRB cohort** (`isalgo/airr_benchmark`,
        78 donors, ages **0 → 103**) using vdjtools v2's **cohort-streaming** layer: the cohort is
        ingested once into a hive-partitioned Parquet dataset, then every statistic is a *single*
        streamed `group_by` over `vdjtools.io.scan_cohort` — memory stays flat no matter how many
        samples you scan.

        - **Diversity** falls with age — `stats.diversity_cohort` computes the whole cohort's
          diversity table from the count-frequency spectrum in one streamed pass.
        - **The clone-size distribution shifts**: the singleton fraction drops and the
          *hyperexpanded* read fraction rises — the naive-to-clonal transition.
        - **Spectratype** and **V-gene usage** are one fused `group_by(["sample_id", …])` each.

        Data auto-loads from HuggingFace, preferring a local `~/hf/` or `./` copy.
        """
    )
    return


@app.cell
def _():
    # --- imports & configuration (single cell so every name is defined once) ---
    from pathlib import Path

    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import polars as pl
    from scipy.stats import spearmanr

    from vdjtools.io import ingest_cohort, read_metadata, scan_cohort
    from vdjtools.io.schema import COUNT
    from vdjtools.stats import diversity_cohort, segment_usage, spectratype

    REPO = "isalgo/airr_benchmark"
    HF_FOLDER = "vdjtools"                 # the full-depth aging cohort (not vdjtools_lite/)
    OKABE = {"blue": "#0072B2", "vermillion": "#D55E00", "green": "#009E73",
             "orange": "#E69F00", "purple": "#CC79A7", "grey": "#8C8C8C"}

    def local_base(mo_dir):
        """Return a directory holding the aging files, preferring a local mirror."""
        for root in (mo_dir, Path.cwd(), Path.home() / "hf" / "airr_benchmark"):
            if (root / HF_FOLDER / "metadata_aging.txt").exists():
                return root / HF_FOLDER
        return None

    return (COUNT, HF_FOLDER, OKABE, Path, REPO, diversity_cohort, ingest_cohort,
            local_base, mo, np, pl, read_metadata, scan_cohort, segment_usage,
            spearmanr, spectratype)


@app.cell
def _(mo):
    n_samples = mo.ui.slider(12, 78, value=40, step=2, label="samples across the age range")
    n_samples
    return (n_samples,)


@app.cell
def _(HF_FOLDER, Path, REPO, ingest_cohort, local_base, mo, n_samples, pl,
      read_metadata, scan_cohort):
    _nb = mo.notebook_dir() or Path.cwd()
    _base = local_base(_nb)
    if _base is not None:                                  # local mirror present
        meta_all = read_metadata(_base / "metadata_aging.txt")
        base = _base
    else:                                                  # fetch metadata + selected files
        import huggingface_hub as _hub
        _mp = _hub.hf_hub_download(REPO, f"{HF_FOLDER}/metadata_aging.txt", repo_type="dataset")
        meta_all = read_metadata(_mp)
        base = None
    meta_all = meta_all.with_columns(pl.col("age").cast(pl.Int64, strict=False)).sort("age")
    # even coverage of the age range
    meta = meta_all.gather_every(max(1, meta_all.height // n_samples.value)).head(n_samples.value)
    if base is None:
        import huggingface_hub as _hub
        _root = Path(_hub.snapshot_download(
            REPO, repo_type="dataset",
            allow_patterns=[f"{HF_FOLDER}/{s}.txt.gz" for s in meta["sample_id"]]))
        base = _root / HF_FOLDER

    data_dir = (mo.notebook_dir() or Path.cwd()) / ".data" / "aging_nb"
    cohort = data_dir / f"cohort_{meta.height}"
    if not any(cohort.glob("sample_id=*/*.parquet")):
        ingest_cohort(meta.select("sample_id", "age", "sex"), base, cohort,
                      sample_col="sample_id", file_template="{sample}.txt.gz", fmt="vdjtools")
    lf = scan_cohort(cohort, join_metadata=False)
    ages = meta.select(pl.col("sample_id"), pl.col("age").cast(pl.Int64))
    mo.md(f"**{meta.height} donors**, ages {meta['age'].min()}–{meta['age'].max()}. "
          f"Cohort scanned from `{cohort.name}`.")
    return ages, lf, meta


@app.cell
def _(mo):
    mo.md(r"""## 1 · Diversity falls with age (one streamed pass)""")
    return


@app.cell
def _(OKABE, ages, diversity_cohort, lf, np, pl, plt, spearmanr):
    div = diversity_cohort(lf).join(ages, on="sample_id").sort("age")
    _age = div["age"].to_numpy()
    _sh = div["shannon_wiener"].to_numpy()
    _r, _p = spearmanr(_age, _sh)
    fig1, ax1 = plt.subplots(figsize=(7.0, 4.2))
    ax1.scatter(_age, _sh, s=34, color=OKABE["vermillion"], alpha=0.8)
    _z = np.polyfit(_age, np.log(_sh), 1)
    _xs = np.linspace(_age.min(), _age.max(), 50)
    ax1.plot(_xs, np.exp(np.polyval(_z, _xs)), "--", color=OKABE["grey"])
    ax1.set(xlabel="Age (years)", ylabel="Shannon–Wiener diversity",
            title=f"TCR diversity declines with age (Spearman ρ={_r:.2f}, p={_p:.1g})")
    ax1.spines[["top", "right"]].set_visible(False)
    fig1.tight_layout()
    fig1
    return (div,)


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 2 · The clone-size distribution shifts: naive → clonal

        The single most legible aging signal. The **singleton fraction** (clonotypes seen once —
        the naive-rich tail) falls, while the **hyperexpanded fraction** (reads carried by clones
        above 1% frequency) rises. Both are one streamed `group_by(sample_id)` over the cohort.
        """
    )
    return


@app.cell
def _(COUNT, OKABE, ages, lf, pl, plt, spearmanr):
    size = (lf.group_by("sample_id").agg(
                (pl.col(COUNT) == 1).sum().alias("_s"),
                pl.len().alias("_n"),
                pl.col(COUNT).sum().alias("_reads"),
                pl.col(COUNT).filter(pl.col(COUNT) / pl.col(COUNT).sum() > 0.01)
                  .sum().alias("_hx"))
            .collect(engine="streaming")
            .with_columns((pl.col("_s") / pl.col("_n")).alias("singleton_frac"),
                          (pl.col("_hx").fill_null(0) / pl.col("_reads")).alias("hyperexpanded_frac"))
            .join(ages, on="sample_id").sort("age"))
    _age = size["age"].to_numpy()
    fig2, ax2 = plt.subplots(figsize=(7.0, 4.2))
    ax2.scatter(_age, size["singleton_frac"], s=30, color=OKABE["blue"], label="singleton fraction")
    ax2.scatter(_age, size["hyperexpanded_frac"], s=30, color=OKABE["vermillion"],
                marker="s", label="hyperexpanded read fraction (>1%)")
    _r1 = spearmanr(_age, size["singleton_frac"].to_numpy())[0]
    _r2 = spearmanr(_age, size["hyperexpanded_frac"].to_numpy())[0]
    ax2.set(xlabel="Age (years)", ylabel="Fraction",
            title=f"Naive→clonal shift  (singleton ρ={_r1:.2f}, hyperexpanded ρ={_r2:.2f})")
    ax2.legend(frameon=False, fontsize=9)
    ax2.spines[["top", "right"]].set_visible(False)
    fig2.tight_layout()
    fig2
    return (size,)


@app.cell
def _(mo):
    mo.md(r"""## 3 · CDR3-length spectratype by age group (one fused pass)""")
    return


@app.cell
def _(OKABE, ages, lf, np, pl, plt, spectratype):
    sp = (spectratype(lf, by=["sample_id"], by_locus=False, weight="freq")
          .collect(engine="streaming").join(ages, on="sample_id"))
    # bin donors into three age groups, average the length distribution within each
    _grp = sp.with_columns(pl.when(pl.col("age") < 20).then(pl.lit("0–19"))
                           .when(pl.col("age") < 60).then(pl.lit("20–59"))
                           .otherwise(pl.lit("60+")).alias("agegrp"))
    _agg = (_grp.group_by("agegrp", "length").agg(pl.col("weight").mean())
            .sort("length"))
    fig3, ax3 = plt.subplots(figsize=(7.0, 4.2))
    for _g, _c in zip(["0–19", "20–59", "60+"], [OKABE["green"], OKABE["orange"], OKABE["purple"]]):
        _s = _agg.filter(pl.col("agegrp") == _g).sort("length")
        ax3.plot(_s["length"], _s["weight"], "-o", ms=3, color=_c, label=_g)
    ax3.set(xlabel="CDR3 amino-acid length", ylabel="Mean clonotype-frequency weight",
            title="CDR3-length spectratype by age group")
    ax3.legend(title="Age", frameon=False)
    ax3.spines[["top", "right"]].set_visible(False)
    fig3.tight_layout()
    fig3
    return (sp,)


@app.cell
def _(mo):
    mo.md(
        r"""
        ---
        **Takeaway.** Every panel is one streamed `group_by` over a `scan_cohort` LazyFrame —
        `diversity_cohort` for the diversity table, a size-bucket `group_by` for the naive→clonal
        shift, and a fused `spectratype(by=["sample_id"])` for the length distribution — so peak
        memory is independent of cohort size (the point of `vdjtools.io.ingest_cohort` /
        `scan_cohort`). The coverage-standardised diversity + inter-repertoire drift (MDS) view of
        the same cohort is `examples/aging_airr_benchmark.py`.
        """
    )
    return


@app.cell
def _(mo):
    mo.md(r"""&nbsp;""")
    return


if __name__ == "__main__":
    app.run()
