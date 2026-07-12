"""Repertoire preprocessing pipeline with vdjtools v2.

A marimo notebook. Launch it with

    marimo edit examples/preprocess.py      # interactive editor
    marimo run  examples/preprocess.py      # read-only served app

Raw immunosequencing data is noisy: it carries non-coding rearrangements, PCR/sequencing
errors, cross-sample contamination, uneven sequencing depth, and systematic batch biases.
`vdjtools.preprocess` is a toolkit of pure-polars cleaning steps; this notebook walks a
few real Britanova aging samples through them and shows the before/after effect of each —
`filter_functional`, `correct` (error-collapse), `downsample`, `filter_frequency` /
`filter_segment`, `decontaminate`, `pool_samples` / `join_samples`, and
`correct_vj_usage` (VJ-usage batch-effect correction).

Data auto-downloads from ``isalgo/airr_benchmark`` (folder ``vdjtools/``) into the
gitignored ``examples/.data/preprocess_nb/`` cache.
"""
import marimo

__generated_with = "0.23.14"
app = marimo.App(width="medium")


@app.cell
def _(mo):
    mo.md(
        r"""
        # Cleaning a repertoire — the vdjtools v2 preprocessing pipeline

        Before any diversity, overlap, or biomarker analysis, a raw repertoire needs
        cleaning. `vdjtools.preprocess` is a set of composable, pure-polars steps over the
        canonical clonotype frame. We run a handful of real Britanova TRB samples (three
        sequencing **batches**, `A2*/A4*/A6*`) through the pipeline and show what each step
        does:

        **functional filter → error-correct → downsample → frequency/segment filter →
        decontaminate → pool / join → batch-correct VJ usage.**
        """
    )
    return


@app.cell
def _():
    from pathlib import Path

    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import polars as pl

    from vdjtools import io as vio
    from vdjtools.preprocess import (correct, correct_vj_usage, decontaminate,
                                     downsample, filter_frequency, filter_functional,
                                     filter_segment, join_samples, pool_samples)

    REPO_ID, HF_FOLDER = "isalgo/airr_benchmark", "vdjtools"
    # 3 batches × 3 samples (batch = the A2/A4/A6 prefix of the sample id).
    BATCHES = {"A2": 3, "A4": 3, "A6": 3}
    OKABE = {"blue": "#0072B2", "vermillion": "#D55E00", "green": "#009E73",
             "orange": "#E69F00", "purple": "#CC79A7", "grey": "#8C8C8C"}
    return (BATCHES, HF_FOLDER, OKABE, Path, REPO_ID, correct, correct_vj_usage,
            decontaminate, downsample, filter_frequency, filter_functional,
            filter_segment, join_samples, mo, np, pl, plt, pool_samples, vio)


@app.cell
def _(mo):
    mo.md(r"""## 1 · Load raw samples (three batches)""")
    return


@app.cell
def _(BATCHES, HF_FOLDER, Path, REPO_ID, mo, pl, vio):
    _nb_dir = mo.notebook_dir() or Path.cwd()
    data_dir = _nb_dir / ".data" / "preprocess_nb"
    data_dir.mkdir(parents=True, exist_ok=True)
    try:
        import huggingface_hub as _hub
    except ImportError:
        _hub = None
    mo.stop(_hub is None, mo.md("> `pip install \"vdjtools[examples]\"`."))

    _meta = vio.read_metadata(_hub.hf_hub_download(
        REPO_ID, f"{HF_FOLDER}/metadata_aging.txt", repo_type="dataset"))
    # Pick N samples per batch prefix.
    _picks = []
    for _pref, _k in BATCHES.items():
        _ids = _meta.filter(pl.col("sample_id").str.starts_with(_pref))["sample_id"].to_list()
        _picks += [(s, _pref) for s in _ids[:_k]]

    _root = Path(_hub.snapshot_download(
        REPO_ID, repo_type="dataset",
        allow_patterns=[f"{HF_FOLDER}/{s}.txt.gz" for s, _ in _picks]))
    raw = {}
    for _s, _b in _picks:
        raw[_s] = vio.read(_root / HF_FOLDER / f"{_s}.txt.gz", fmt="vdjtools").with_columns(
            pl.lit(_s).alias("sample_id"), pl.lit(_b).alias("batch"))
    mo.md(f"Loaded **{len(raw)} raw samples** across {len(BATCHES)} batches "
          f"({', '.join(BATCHES)}). Cache: `{data_dir}`")
    return data_dir, raw


@app.cell
def _(mo, pl, raw):
    # Raw per-sample summary: clonotypes, reads, and % non-coding (a QC red flag).
    _rows = []
    for _s, _df in raw.items():
        _nc = _df.filter(pl.col("cdr3_aa").str.contains("[*_]")).height
        _rows.append({"sample": _s, "batch": _df["batch"][0], "clonotypes": _df.height,
                      "reads": int(_df["duplicate_count"].sum()),
                      "pct_noncoding": round(100 * _nc / _df.height, 2)})
    summary = pl.DataFrame(_rows)
    mo.vstack([mo.md("**Raw samples** — note the non-coding fraction:"), summary])
    return (summary,)


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 2 · Functional filter + error correction

        `filter_functional(keep="coding")` drops out-of-frame / stop-codon rearrangements
        (the `*`/`_` CDR3s). `correct` then collapses low-count clonotypes whose `cdr3_nt`
        sit within a couple of mismatches of a much more abundant "parent" — the classic
        PCR/sequencing-error signature — merging their counts upward.
        """
    )
    return


@app.cell
def _(correct, filter_functional, mo, pl, raw):
    _s = next(iter(raw))
    _r = raw[_s]
    _coding = filter_functional(_r, keep="coding")
    _corrected = correct(_coding, max_mismatches=2, ratio=0.05)
    steps = pl.DataFrame({
        "step": ["raw", "coding only", "error-corrected"],
        "clonotypes": [_r.height, _coding.height, _corrected.height],
        "reads": [int(_r["duplicate_count"].sum()), int(_coding["duplicate_count"].sum()),
                  int(_corrected["duplicate_count"].sum())],
    })
    mo.vstack([mo.md(f"Sample **{_s}** — coding filter + error correction "
                     f"(reads conserved, clonotypes collapse):"), steps])
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 3 · Downsample + filters

        `downsample` equalises sequencing depth (multinomial resample to a common read
        count) — essential before cross-sample comparison. `filter_frequency` /
        `filter_segment` subset by abundance or by gene.
        """
    )
    return


@app.cell
def _(OKABE, downsample, filter_frequency, filter_segment, mo, np, plt, raw):
    _s, _df = next(iter(raw.items()))
    _reads = int(_df["duplicate_count"].sum())
    _depth = _reads // 4
    _ds = downsample(_df, size=_depth, by="reads", seed=0)
    _topq = filter_frequency(_df, top_quantile=0.01)          # top 1% most abundant
    _v9 = filter_segment(_df, v=["TRBV9"], keep=True)         # keep only TRBV9

    _labels = ["raw", f"downsample→{_depth:,}", "top-1% freq", "TRBV9 only"]
    _counts = [_df.height, _ds.height, _topq.height, _v9.height]
    figf, axf = plt.subplots(figsize=(7.2, 4.0))
    axf.bar(_labels, _counts, color=[OKABE["grey"], OKABE["blue"], OKABE["orange"],
                                     OKABE["green"]])
    for _i, _c in enumerate(_counts):
        axf.text(_i, _c, f"{_c:,}", ha="center", va="bottom", fontsize=9)
    axf.set(ylabel="clonotypes", title=f"Downsampling & filtering — {_s}")
    axf.set_yscale("log")
    axf.spines[["top", "right"]].set_visible(False)
    figf.tight_layout()
    figf
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 4 · Decontaminate, pool, and join

        `decontaminate` removes clonotypes that are far more abundant in a *different*
        sample (cross-sample / index-hopping contamination). `pool_samples` unions
        clonotypes summing counts; `join_samples` keeps clonotypes seen in ≥ `min_samples`
        samples as a wide sample × clonotype table (the overlap/tracking substrate).
        """
    )
    return


@app.cell
def _(decontaminate, join_samples, mo, pl, pool_samples, raw):
    _samples = list(raw.values())
    _a = _samples[0]
    _clean = decontaminate(_a, others=_samples[1:], ratio=20.0)  # 20× dominance elsewhere
    _pooled = pool_samples(_samples, key="aa")
    _joined = join_samples(_samples, key="aa", min_samples=2)

    out = pl.DataFrame({
        "operation": ["sample A (raw)", "A decontaminated",
                      "pool (all, union)", "join (present in ≥2)"],
        "clonotypes": [_a.height, _clean.height, _pooled.height, _joined.height],
    })
    mo.vstack([mo.md(f"Removed **{_a.height - _clean.height}** contaminant clonotypes from "
                     f"sample A; pooled/joined sizes across {len(_samples)} samples:"), out])
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 5 · Batch-effect correction of VJ usage

        Different sequencing batches carry systematic V/J usage biases. `correct_vj_usage`
        removes the per-`(gene, batch)` log-space offset (the ComBat location term). Below:
        PCA of samples on their V-J usage vectors, **before** vs **after** correction,
        coloured by batch. Before, samples segregate by batch; after, the batch separation
        collapses while sample-to-sample structure is preserved.
        """
    )
    return


@app.cell
def _(OKABE, correct_vj_usage, np, pl, plt, raw):
    _long = pl.concat(list(raw.values()), how="vertical_relaxed")
    _u = correct_vj_usage(_long, batch_col="batch", sample_col="sample_id")

    def _matrix(col):
        _w = _u.with_columns((pl.col("v_call") + "|" + pl.col("j_call")).alias("vj")).pivot(
            values=col, index="sample_id", on="vj", aggregate_function="first").fill_null(0.0)
        _s = _w["sample_id"].to_list()
        return _s, _w.drop("sample_id").to_numpy()

    from sklearn.decomposition import PCA
    _samps, _raw_m = _matrix("p")
    _, _cor_m = _matrix("p_corrected")
    _batch = [s.split("-")[0] for s in _samps]
    _bset = sorted(set(_batch))
    _cmap = {_b: list(OKABE.values())[_i] for _i, _b in enumerate(_bset)}
    _col = [_cmap[_b] for _b in _batch]

    figb, axb = plt.subplots(1, 2, figsize=(11, 4.6))
    for _ax, _m, _t in zip(axb, (_raw_m, _cor_m), ("raw VJ usage", "batch-corrected")):
        _xy = PCA(n_components=2, random_state=0).fit_transform(
            np.log(_m + 1e-9) - np.log(_m + 1e-9).mean(0))
        _ax.scatter(_xy[:, 0], _xy[:, 1], c=_col, s=90, edgecolor="white", linewidth=0.6)
        _ax.set(title=_t, xlabel="PC1", ylabel="PC2")
        _ax.spines[["top", "right"]].set_visible(False)
    _handles = [plt.Line2D([], [], marker="o", ls="", mfc=_cmap[_b], mec="white", label=_b)
                for _b in _bset]
    axb[1].legend(handles=_handles, title="batch", fontsize=8, frameon=False)
    figb.suptitle("VJ-usage batch effect removed by correct_vj_usage", y=1.02)
    figb.tight_layout()
    figb
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ---
        **Takeaway.** Each step is one call over the canonical clonotype frame:
        `filter_functional` drops non-coding rearrangements, `correct` collapses
        error-variants (reads conserved, clonotypes down), `downsample` equalises depth,
        `filter_frequency`/`filter_segment` subset by abundance or gene, `decontaminate`
        strips cross-sample bleed, `pool_samples`/`join_samples` combine repertoires, and
        `correct_vj_usage` removes systematic batch biases in gene usage — all pure polars,
        composable into a preprocessing pipeline ahead of `vdjtools.stats` /
        `vdjtools.overlap` / `vdjtools.biomarker`.
        """
    )
    return


if __name__ == "__main__":
    app.run()
