"""Age-related TCR repertoire (aging) analysis with vdjtools v2.

A marimo notebook. Edit interactively with

    marimo edit examples/aging_airr_benchmark.py

or run headless with

    marimo run examples/aging_airr_benchmark.py

The data (Britanova human TRB age cohort, native vdjtools format) auto-downloads
from the HuggingFace dataset ``isalgo/airr_benchmark`` (folder ``vdjtools_lite/``)
into the gitignored ``examples/.data/aging/`` directory and is md5-cached against
the committed ``examples/aging_manifest.json`` — a second run fetches nothing.
"""
import marimo

__generated_with = "0.23.14"
app = marimo.App(width="medium")


@app.cell
def _(mo):
    mo.md(
        r"""
        # Aging of the human TCR repertoire — a vdjtools v2 walkthrough

        As we age, the naive T-cell pool contracts and a handful of clones expand,
        so the T-cell receptor (TCR) repertoire becomes **less diverse and more
        clonal**. This notebook reproduces that classic signal on the **Britanova
        human TRB age cohort** (41 donors, ages 6–90) using the basic-analytics
        layer of **vdjtools v2** (`vdjtools.io`, `vdjtools.stats`,
        `vdjtools.features`).

        The data are the legacy vdjtools *aging_lite* example in **native vdjtools
        format**, fetched from the HuggingFace dataset
        [`isalgo/airr_benchmark`](https://huggingface.co/datasets/isalgo/airr_benchmark)
        (`vdjtools_lite/`). The first cell bootstraps them into a gitignored,
        **md5-verified** local cache.
        """
    )
    return


@app.cell
def _():
    # --- imports & configuration (single cell so every name is defined once) ---
    import hashlib
    import json
    import shutil
    from pathlib import Path

    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import polars as pl
    from scipy.stats import spearmanr

    from vdjtools import io as vio
    from vdjtools.stats import estimate_d, inext_batch, rarefaction, sample_coverage

    # HuggingFace dataset coordinates (see examples/README.md and SOURCES.md).
    REPO_ID = "isalgo/airr_benchmark"
    HF_FOLDER = "vdjtools_lite"

    # Okabe–Ito colorblind-safe palette.
    OKABE = {
        "blue": "#0072B2", "vermillion": "#D55E00", "green": "#009E73",
        "orange": "#E69F00", "purple": "#CC79A7", "grey": "#8C8C8C",
    }

    def file_md5(path):
        """Streaming md5 of a file (chunked, so large gz files stay cheap)."""
        h = hashlib.md5()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()

    return (
        HF_FOLDER, OKABE, Path, REPO_ID, estimate_d, file_md5, inext_batch,
        json, mo, np, pl, plt, rarefaction, sample_coverage, shutil, spearmanr, vio,
    )


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 1 · md5-gated data bootstrap

        The committed `aging_manifest.json` maps every file
        (`metadata_aging.txt` + 41 samples) to its md5. For each file we

        1. **skip entirely** (no network) if it is already cached with the right md5, else
        2. `hf_hub_download` it, **verify** the downloaded md5 against the manifest, and copy it in.

        The cache lives in the gitignored `examples/.data/aging/`, so a re-run does
        zero downloads.
        """
    )
    return


@app.cell
def _(HF_FOLDER, Path, REPO_ID, file_md5, json, mo, shutil):
    # Resolve the notebook dir robustly (falls back to cwd outside a marimo runtime).
    _nb_dir = mo.notebook_dir() or Path.cwd()
    manifest = json.loads((_nb_dir / "aging_manifest.json").read_text())

    data_dir = _nb_dir / ".data" / "aging"
    data_dir.mkdir(parents=True, exist_ok=True)

    # huggingface_hub ships in the [examples] extra — fail gracefully if it is missing.
    try:
        import huggingface_hub as _hub
    except ImportError:
        _hub = None
    mo.stop(
        _hub is None,
        mo.md(
            """
            > **`huggingface_hub` is not installed.** This notebook fetches its data
            > from HuggingFace. Install the examples extra and re-run:
            >
            > ```
            > pip install "vdjtools[examples]"
            > ```
            """
        ),
    )

    _fetched, _skipped = 0, 0
    for _name, _want in manifest.items():
        _dest = data_dir / _name
        if _dest.exists() and file_md5(_dest) == _want:
            _skipped += 1
            continue
        _src = _hub.hf_hub_download(
            repo_id=REPO_ID, filename=f"{HF_FOLDER}/{_name}", repo_type="dataset"
        )
        _got = file_md5(_src)
        if _got != _want:
            raise ValueError(f"md5 mismatch for {_name}: want {_want}, got {_got}")
        shutil.copyfile(_src, _dest)
        _fetched += 1

    mo.md(
        f"**{_skipped + _fetched}/{len(manifest)} files present, md5 verified** — "
        f"{_fetched} fetched, {_skipped} skipped.  \n"
        f"Cache: `{data_dir}`"
    )
    return (data_dir,)


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 2 · Load the cohort

        `read_metadata` reads the sample sheet; `read_samples` reads every native
        `.txt.gz` into one **long clonotype frame** (AIRR-style canonical columns),
        tagging each row with its `sample_id` and joining the `age` / `sex`
        metadata. One call, one tidy polars frame.
        """
    )
    return


@app.cell
def _(data_dir, pl, vio):
    meta = vio.read_metadata(data_dir / "metadata_aging.txt")
    long = vio.read_samples(
        meta, data_dir, sample_col="sample_id", file_template="{sample}.txt.gz"
    ).with_columns(pl.col("age").cast(pl.Int64))
    samples = meta["sample_id"].to_list()
    return long, samples


@app.cell
def _(long, pl):
    # Per-sample summary: sequencing depth (reads) and richness (clonotypes).
    summary = (
        long.group_by("sample_id", maintain_order=True)
        .agg(
            pl.len().alias("n_clonotypes"),
            pl.col("duplicate_count").sum().alias("reads"),
            pl.col("age").first(),
            pl.col("sex").first(),
        )
        .sort("age")
    )
    summary
    return (summary,)


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 3 · Diversity vs age — the headline (coverage-standardized)

        Comparing raw diversity across samples is a trap: a deeper (or more
        completely sampled) repertoire looks more diverse just from more sequencing.
        The principled fix is **coverage standardization** — evaluate every sample's
        Hill-number diversity at a common *sample coverage* Ĉ (the fraction of the
        repertoire the reads represent), following the iNEXT framework.

        > **Why it matters here.** These *lite* samples are all downsampled to a
        > common **10,000 reads**, yet their sample coverage Ĉ ranges from ~0.06
        > (young, almost every clone seen once) to ~0.79 (old, dominated by a few
        > expanded clones) — a >10× spread. Equal depth does **not** mean equal
        > completeness, so we standardize on coverage rather than depth. We fix the
        > common coverage at **C = min Ĉ** across samples and read off Hill
        > diversity there for orders q = 0 (richness), 1 (Shannon) and 2 (Simpson).
        """
    )
    return


@app.cell
def _(estimate_d, inext_batch, long, pl, sample_coverage, samples, spearmanr):
    # Split once; measure each sample's observed coverage Ĉ(n).
    _parts = {s: long.filter(pl.col("sample_id") == s) for s in samples}
    _age = {s: int(_parts[s]["age"][0]) for s in samples}
    cov_df = pl.DataFrame(
        {
            "sample": samples,
            "age": [_age[s] for s in samples],
            "coverage": [sample_coverage(_parts[s]) for s in samples],
        }
    )
    C = float(cov_df["coverage"].min())  # common standardization coverage

    # Coverage-standardized Hill diversity (q = 0, 1, 2) with bootstrap CIs.
    _rows = []
    for s in samples:
        _d = estimate_d(_parts[s], base="coverage", level=C, q=(0, 1, 2),
                        se=True, nboot=50, seed=0)
        for _r in _d.iter_rows(named=True):
            _rows.append({
                "sample": s, "age": _age[s], "order_q": _r["order_q"],
                "m": _r["m"], "qD": _r["qD"],
                "qD_lo": _r["qD_lo"], "qD_hi": _r["qD_hi"],
            })
    div_std = pl.DataFrame(_rows)

    # Spearman(age, standardized diversity) per Hill order.
    spear = {}
    for _q in (0, 1, 2):
        _sub = div_std.filter(pl.col("order_q") == _q)
        _rr, _pp = spearmanr(_sub["age"].to_numpy(), _sub["qD"].to_numpy())
        spear[_q] = (float(_rr), float(_pp))

    # Robustness: the raw (size-based) Shannon from the fast native batch engine.
    _batch = inext_batch(long, q=(1,), knots=12, se=False, seed=0)
    _obs = (_batch.filter(pl.col("method") == "observed")
            .join(cov_df.select(["sample", "age"]), on="sample"))
    _rr, _pp = spearmanr(_obs["age"].to_numpy(), _obs["qD"].to_numpy())
    raw_shannon = (float(_rr), float(_pp))
    return C, cov_df, div_std, raw_shannon, spear


@app.cell
def _(C, cov_df, div_std, np, pl, plt, spear):
    _cov = cov_df.sort("age")
    _q1 = div_std.filter(pl.col("order_q") == 1).sort("age")
    _age = _q1["age"].to_numpy()
    _qd = _q1["qD"].to_numpy()
    _err = np.vstack([_qd - _q1["qD_lo"].to_numpy(), _q1["qD_hi"].to_numpy() - _qd])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.4))

    # (A) Motivation: sample coverage varies hugely at equal 10k-read depth.
    ax1.scatter(_cov["age"], _cov["coverage"], s=42, color="#0072B2",
                edgecolor="white", linewidth=0.6, zorder=3)
    ax1.axhline(C, color="#8C8C8C", ls="--", lw=1)
    ax1.text(0.98, C, f"  common C = {C:.3f}", va="bottom", ha="right",
             transform=ax1.get_yaxis_transform(), color="#5c5c5c", fontsize=9)
    ax1.set(xlabel="Age (years)", ylabel="Sample coverage  Ĉ(n)",
            title="Equal depth ≠ equal completeness")
    ax1.spines[["top", "right"]].set_visible(False)

    # (B) Headline: coverage-standardized Shannon diversity vs age, with CIs.
    ax2.errorbar(_age, _qd, yerr=_err, fmt="o", ms=6, color="#D55E00",
                 ecolor="#D55E0055", elinewidth=1.2, capsize=2,
                 mec="white", mew=0.6, zorder=3)
    ax2.set_yscale("log")
    ax2.set(xlabel="Age (years)",
            ylabel="Coverage-standardized Shannon  ¹D  (log)",
            title="TCR diversity declines with age")
    _r, _p = spear[1]
    ax2.text(0.04, 0.06, f"Spearman r = {_r:.2f}\np = {_p:.1e}",
             transform=ax2.transAxes, fontsize=10, va="bottom",
             bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#cccccc"))
    ax2.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    fig
    return (fig,)


@app.cell
def _(mo, raw_shannon, spear):
    _r0, _p0 = spear[0]
    _r1, _p1 = spear[1]
    _r2, _p2 = spear[2]
    _rr, _rp = raw_shannon
    mo.md(
        f"""
        **Result.** Coverage-standardized diversity falls with age across every
        Hill order — richness (q=0): r = {_r0:.2f}, p = {_p0:.1e}; Shannon (q=1):
        r = {_r1:.2f}, p = {_p1:.1e}; Simpson (q=2): r = {_r2:.2f}, p = {_p2:.1e}.
        The trend is robust to the standardization choice: the raw, size-based
        Shannon (native `inext_batch`) gives essentially the same r = {_rr:.2f}
        (p = {_rp:.1e}) — reassuring, since here the samples are already depth-matched.
        The negative correlation is the expected aging biology: older repertoires
        support far fewer *effective* clones.
        """
    )
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 4 · Rarefaction / extrapolation curves

        The classic vdjtools rarefaction plot: clonotype **richness** (Hill q=0) as
        a function of sampling depth *m*, interpolated below and extrapolated above
        the observed depth (`rarefaction(sample, q=0)`), with bootstrap confidence
        bands. Young donors climb steeply and show no sign of saturating; old donors
        plateau early — their repertoires are shallower.
        """
    )
    return


@app.cell
def _(OKABE, long, np, pl, plt, rarefaction, summary):
    _ordered = summary.sort("age")
    _young = _ordered.head(2)
    _old = _ordered.tail(2)
    _picks = (
        [(s, a, OKABE["blue"], "young") for s, a in
         zip(_young["sample_id"], _young["age"])]
        + [(s, a, OKABE["vermillion"], "old") for s, a in
           zip(_old["sample_id"], _old["age"])]
    )

    figr, axr = plt.subplots(figsize=(7.2, 4.6))
    for _s, _a, _col, _grp in _picks:
        _clones = long.filter(pl.col("sample_id") == _s)
        _rc = rarefaction(_clones, q=0, knots=25, se=True, nboot=30, seed=0)
        _n = int(_clones["duplicate_count"].sum())
        _m = _rc["m"].to_numpy()
        _qd = _rc["qD"].to_numpy()
        _interp = _rc["method"].to_numpy() != "extrapolation"
        axr.fill_between(_m, _rc["qD_lo"].to_numpy(), _rc["qD_hi"].to_numpy(),
                         color=_col, alpha=0.12, lw=0)
        axr.plot(_m[_interp], _qd[_interp], color=_col, lw=2)
        axr.plot(_m[~_interp], _qd[~_interp], color=_col, lw=2, ls="--")
        axr.plot([_n], [np.interp(_n, _m, _qd)], "o", color=_col, mec="white", ms=7)
        axr.text(_m[-1], _qd[-1], f"  {_s} ({_a}y, {_grp})",
                 color=_col, va="center", fontsize=9)

    axr.set(xlabel="Sampling depth  m  (reads)",
            ylabel="Clonotype richness  ⁰D",
            title="Rarefaction / extrapolation (q=0): young vs old donors")
    axr.margins(x=0.18)
    axr.spines[["top", "right"]].set_visible(False)
    figr.tight_layout()
    figr
    return (figr,)


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 5 · Clonal expansion vs age

        Diversity's flip side, computed straight from the canonical frame: the
        fraction of the repertoire occupied by the **top 10 clones**. Age-associated
        clonal expansions should push this up in older donors.
        """
    )
    return


@app.cell
def _(long, np, pl, plt, spearmanr):
    # frequency is per-sample normalized, so summing the 10 largest gives the
    # read share of the dominant clones.
    clonal = (
        long.group_by("sample_id", maintain_order=True)
        .agg(
            pl.col("frequency").top_k(10).sum().alias("top10_frac"),
            pl.col("age").first(),
        )
        .sort("age")
    )
    _age = clonal["age"].to_numpy()
    _frac = clonal["top10_frac"].to_numpy() * 100.0
    _r, _p = spearmanr(_age, _frac)

    figc, axc = plt.subplots(figsize=(7.0, 4.4))
    axc.scatter(_age, _frac, s=48, color="#009E73", edgecolor="white",
                linewidth=0.6, zorder=3)
    _b, _a0 = np.polyfit(_age, _frac, 1)
    _xs = np.array([_age.min(), _age.max()])
    axc.plot(_xs, _a0 + _b * _xs, color="#8C8C8C", ls="--", lw=1.4)
    axc.set(xlabel="Age (years)",
            ylabel="Top-10-clone read share (%)",
            title="Repertoires become more clonal with age")
    axc.text(0.04, 0.94, f"Spearman r = {_r:.2f}\np = {_p:.1e}",
             transform=axc.transAxes, fontsize=10, va="top",
             bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#cccccc"))
    axc.spines[["top", "right"]].set_visible(False)
    figc.tight_layout()
    figc
    return (clonal, figc)


@app.cell
def _(mo):
    mo.md(
        r"""
        ---
        **Takeaway.** Across 41 healthy donors, the TRB repertoire loses diversity
        and gains clonality with age — reproduced end-to-end from raw native
        vdjtools files with `vdjtools.io`, `vdjtools.stats` (iNEXT diversity,
        coverage standardization, rarefaction) and a direct clonal-expansion
        summary. The bootstrap is md5-cached, so re-running this notebook touches
        the network zero times.
        """
    )
    return


if __name__ == "__main__":
    app.run()
