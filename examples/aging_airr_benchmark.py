"""Age-related TCR repertoire (aging) analysis with vdjtools v2.

A marimo notebook. Launch it with

    examples/run.sh                 # -> marimo edit (interactive)

or directly

    marimo edit examples/aging_airr_benchmark.py
    marimo run  examples/aging_airr_benchmark.py     # read-only served app

The data (full-depth Britanova human TRB "Cord Blood to Centenarians" cohort,
native vdjtools format) auto-downloads from the HuggingFace dataset
``isalgo/airr_benchmark`` (folder ``vdjtools/``) into the gitignored
``examples/.data/aging/`` directory and is md5-cached against the committed
``examples/aging_manifest.json`` — a second run fetches nothing. The full cohort
is ~0.5 GB.
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
        clonal**, and — because those expansions are largely *private* (stochastic,
        person-specific) — individual repertoires **drift apart** with age. This
        notebook reproduces all three signals on the full-depth **Britanova "Cord
        Blood to Centenarians" TRB cohort** (79 donors, ages **0 → 103**, including
        8 cord-blood samples) using the basic-analytics layer of **vdjtools v2**
        (`vdjtools.io`, `vdjtools.stats`, `vdjtools.overlap`).

        The data are the legacy vdjtools example in **native vdjtools format**,
        fetched from the HuggingFace dataset
        [`isalgo/airr_benchmark`](https://huggingface.co/datasets/isalgo/airr_benchmark)
        (folder `vdjtools/`). The first cell bootstraps them into a gitignored,
        **md5-verified** local cache (~0.5 GB total).
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
    from sklearn.manifold import MDS

    from vdjtools import io as vio
    from vdjtools.overlap import overlap_metrics
    from vdjtools.stats import estimate_d, inext_batch, rarefaction, sample_coverage

    # HuggingFace dataset coordinates (see examples/README.md and SOURCES.md).
    REPO_ID = "isalgo/airr_benchmark"
    HF_FOLDER = "vdjtools"

    # Strict exact-match overlap key: CDR3 amino acid + V + J gene.
    KEY = ("cdr3_aa", "v_call", "j_call")

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
        HF_FOLDER, KEY, MDS, OKABE, Path, REPO_ID, estimate_d, file_md5,
        inext_batch, json, mo, np, overlap_metrics, pl, plt, rarefaction,
        sample_coverage, shutil, spearmanr, vio,
    )


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 1 · md5-gated data bootstrap

        The committed `aging_manifest.json` maps every file
        (`metadata_aging.txt` + 79 samples) to its md5. For each file we

        1. **skip entirely** (no network) if it is already cached with the right md5, else
        2. `hf_hub_download` it, **verify** the downloaded md5 against the manifest, and copy it in.

        The cache lives in the gitignored `examples/.data/aging/`, so a re-run does
        zero downloads. First run pulls ~0.5 GB, so allow a few minutes.
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
        ## 2 · Load & featurize the cohort

        `read_metadata` reads the sample sheet (its header line is commented out with
        a leading `#`, which the reader strips). We then **stream** each native
        `.txt.gz` through `vio.read` once, keeping only compact per-sample artifacts —
        the clonotype **count vector** (for diversity) and a **downsampled** copy on
        the CDR3aa+V+J key (for overlap) — so the full-depth cohort (~31 M clonotypes)
        never has to sit in memory all at once.
        """
    )
    return


@app.cell
def _(data_dir, pl, vio):
    meta = vio.read_metadata(data_dir / "metadata_aging.txt").with_columns(
        pl.col("age").cast(pl.Int64)
    )
    samples = meta["sample_id"].to_list()
    ages = meta["age"].to_numpy()
    return ages, meta, samples


@app.cell
def _(KEY, data_dir, np, pl, samples, vio):
    _rng = np.random.default_rng(0)
    counts = {}          # sample_id -> full clonotype count vector (diversity)
    reads = {}           # sample_id -> total reads
    n_clones = {}        # sample_id -> #clonotypes
    _keytab = {}         # transient: aggregated CDR3aa+V+J counts (freed after cell)
    for _s in samples:
        _df = vio.read(data_dir / f"{_s}.txt.gz", fmt="vdjtools")
        _c = _df["duplicate_count"].to_numpy().astype("int64")
        counts[_s] = _c
        reads[_s] = int(_c.sum())
        n_clones[_s] = _df.height
        _keytab[_s] = (_df.group_by(list(KEY), maintain_order=True)
                       .agg(pl.col("duplicate_count").sum().alias("c")))

    # Overlap metrics are depth-sensitive, so downsample every sample to a COMMON
    # depth (the cohort minimum read count) before any cross-sample comparison.
    DEPTH = min(reads.values())
    ds = {}
    for _s in samples:
        _g = _keytab[_s]
        _p = _g["c"].to_numpy().astype("float64")
        _p = _p / _p.sum()
        _new = _rng.multinomial(DEPTH, _p)
        _keep = _new > 0
        ds[_s] = (_g.filter(pl.Series(_keep))
                  .with_columns(pl.Series("duplicate_count", _new[_keep].astype("int64")))
                  .select([*KEY, "duplicate_count"]))
    return DEPTH, counts, ds, n_clones, reads


@app.cell
def _(DEPTH, ages, mo, n_clones, pl, reads, samples):
    summary = pl.DataFrame(
        {
            "sample_id": samples,
            "age": ages,
            "reads": [reads[s] for s in samples],
            "n_clonotypes": [n_clones[s] for s in samples],
        }
    ).sort("age")
    mo.md(f"Common downsampling depth (cohort min reads) = **{DEPTH:,}**.")
    summary
    return (summary,)


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 3 · Diversity vs age — the headline (coverage-standardized)

        Comparing raw diversity across samples is a trap: a deeper or more completely
        sampled repertoire looks more diverse just from more sequencing. The
        principled fix is **coverage standardization** — evaluate every sample's
        Hill-number diversity at a common *sample coverage* Ĉ (the fraction of the
        repertoire the reads represent), following the iNEXT framework. Here the
        full-depth samples span Ĉ ≈ 0.20 → 0.97, so we fix the common coverage at
        **C = min Ĉ** and read off Hill diversity there for orders q = 0 (richness),
        1 (Shannon) and 2 (Simpson).
        """
    )
    return


@app.cell
def _(ages, counts, estimate_d, inext_batch, np, pl, sample_coverage, samples, spearmanr):
    _cov = np.array([sample_coverage(counts[s]) for s in samples])
    cov_df = pl.DataFrame({"sample": samples, "age": ages, "coverage": _cov})
    C = float(_cov.min())

    # Coverage-standardized Hill diversity (q = 0, 1, 2) with (light) bootstrap CIs.
    _rows = []
    for _i, _s in enumerate(samples):
        _d = estimate_d(counts[_s], base="coverage", level=C, q=(0, 1, 2),
                        se=True, nboot=10, seed=0)
        for _r in _d.iter_rows(named=True):
            _rows.append({
                "sample": _s, "age": int(ages[_i]), "order_q": _r["order_q"],
                "qD": _r["qD"], "qD_lo": _r["qD_lo"], "qD_hi": _r["qD_hi"],
            })
    div_std = pl.DataFrame(_rows)

    spear = {}
    for _q in (0, 1, 2):
        _sub = div_std.filter(pl.col("order_q") == _q)
        _rr, _pp = spearmanr(_sub["age"].to_numpy(), _sub["qD"].to_numpy())
        spear[_q] = (float(_rr), float(_pp))

    # Robustness: raw (size-based) Shannon from the fast native batch engine.
    _batch = inext_batch([counts[s] for s in samples], q=(1,), knots=8, se=False, seed=0)
    _obs = _batch.filter(pl.col("method") == "observed").sort("sample")
    _rr, _pp = spearmanr(ages, _obs["qD"].to_numpy())
    raw_shannon = (float(_rr), float(_pp))
    return C, cov_df, div_std, raw_shannon, spear


@app.cell
def _(C, OKABE, cov_df, div_std, np, pl, plt, spear):
    _cov = cov_df.sort("age")
    _q1 = div_std.filter(pl.col("order_q") == 1).sort("age")
    _age = _q1["age"].to_numpy()
    _qd = _q1["qD"].to_numpy()
    _err = np.vstack([_qd - _q1["qD_lo"].to_numpy(), _q1["qD_hi"].to_numpy() - _qd])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.4))

    # (A) Motivation: sample coverage varies widely across the cohort.
    ax1.scatter(_cov["age"], _cov["coverage"], s=34, color=OKABE["blue"],
                edgecolor="white", linewidth=0.5, zorder=3)
    ax1.axhline(C, color=OKABE["grey"], ls="--", lw=1)
    ax1.text(0.98, C, f"  common C = {C:.2f}", va="bottom", ha="right",
             transform=ax1.get_yaxis_transform(), color="#5c5c5c", fontsize=9)
    ax1.set(xlabel="Age (years)", ylabel="Sample coverage  Ĉ(n)",
            title="Sample completeness varies across donors")
    ax1.spines[["top", "right"]].set_visible(False)

    # (B) Headline: coverage-standardized Shannon diversity vs age, with CIs.
    ax2.errorbar(_age, _qd, yerr=_err, fmt="o", ms=5, color=OKABE["vermillion"],
                 ecolor="#D55E0055", elinewidth=1.0, capsize=2,
                 mec="white", mew=0.5, zorder=3)
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
    return


@app.cell
def _(mo, raw_shannon, spear):
    _r0, _p0 = spear[0]
    _r1, _p1 = spear[1]
    _r2, _p2 = spear[2]
    _rr, _rp = raw_shannon
    mo.md(
        f"""
        **Result.** Coverage-standardized diversity falls with age across every Hill
        order — richness (q=0): r = {_r0:.2f}, p = {_p0:.1e}; Shannon (q=1):
        r = {_r1:.2f}, p = {_p1:.1e}; Simpson (q=2): r = {_r2:.2f}, p = {_p2:.1e}.
        The trend is robust to the standardization choice: the raw, size-based
        Shannon from the native `inext_batch` gives r = {_rr:.2f} (p = {_rp:.1e}).
        On this full-depth cohort — spanning cord blood to centenarians — the decline
        is far cleaner than on shallow data: older repertoires support markedly fewer
        *effective* clones.
        """
    )
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 4 · Rarefaction / extrapolation curves

        The classic vdjtools rarefaction plot: clonotype **richness** (Hill q=0) as a
        function of sampling depth *m* (`rarefaction(sample, q=0)`), with bootstrap
        confidence bands, for the youngest (cord-blood) vs the oldest donors. Note the
        log-scaled depth axis — full-depth samples differ ~50× in sequencing depth,
        exactly the confound that the coverage standardization above removes.
        """
    )
    return


@app.cell
def _(OKABE, counts, np, plt, rarefaction, reads, summary):
    _ordered = summary.sort("age")
    _picks = (
        [(s, a, OKABE["blue"], "cord blood") for s, a in
         zip(_ordered.head(2)["sample_id"], _ordered.head(2)["age"])]
        + [(s, a, OKABE["vermillion"], "old") for s, a in
           zip(_ordered.tail(2)["sample_id"], _ordered.tail(2)["age"])]
    )

    figr, axr = plt.subplots(figsize=(7.4, 4.6))
    for _s, _a, _col, _grp in _picks:
        _rc = rarefaction(counts[_s], q=0, knots=20, se=True, nboot=10, seed=0)
        _m = _rc["m"].to_numpy()
        _qd = _rc["qD"].to_numpy()
        _interp = _rc["method"].to_numpy() != "extrapolation"
        axr.fill_between(_m, _rc["qD_lo"].to_numpy(), _rc["qD_hi"].to_numpy(),
                         color=_col, alpha=0.12, lw=0)
        axr.plot(_m[_interp], _qd[_interp], color=_col, lw=2)
        axr.plot(_m[~_interp], _qd[~_interp], color=_col, lw=2, ls="--")
        _n = reads[_s]
        axr.plot([_n], [np.interp(_n, _m, _qd)], "o", color=_col, mec="white", ms=7)
        axr.text(_m[-1], _qd[-1], f"  {_s} ({_a}y, {_grp})",
                 color=_col, va="center", fontsize=9)

    axr.set_xscale("log")
    axr.set(xlabel="Sampling depth  m  (reads, log)",
            ylabel="Clonotype richness  ⁰D",
            title="Rarefaction / extrapolation (q=0): cord blood vs old")
    axr.margins(x=0.25)
    axr.spines[["top", "right"]].set_visible(False)
    figr.tight_layout()
    figr
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 5 · Repertoire divergence with age — overlap MDS

        The headline result of this cohort: repertoires **diverge from a central
        cord-blood cluster outward with age**, driven by *private* clonal expansions.
        We measure it with pairwise **repertoire overlap** (`vdjtools.overlap`, the
        exact-match **F** frequency-overlap metric on the CDR3aa+V+J key). Because
        overlap is depth-sensitive, every sample was first **downsampled to a common
        depth** (the cohort minimum, §2) — the standard normalization before any
        cross-sample overlap comparison. We turn overlap into a distance
        (`d = -log₁₀ F`), embed all 79 samples with **metric MDS**, and colour the
        2-D map by age.
        """
    )
    return


@app.cell
def _(MDS, ages, ds, np, overlap_metrics, samples, spearmanr):
    _n = len(samples)
    F = np.zeros((_n, _n))
    for _i in range(_n):
        for _j in range(_i + 1, _n):
            _f = overlap_metrics(ds[samples[_i]], ds[samples[_j]])["F"]
            F[_i, _j] = F[_j, _i] = _f

    # overlap -> distance; floor any zero-overlap pair at the max finite distance.
    with np.errstate(divide="ignore"):
        D = -np.log10(F)
    _finite_max = D[np.isfinite(D)].max()
    D[~np.isfinite(D)] = _finite_max * 1.1
    np.fill_diagonal(D, 0.0)
    D = (D + D.T) / 2.0

    emb = MDS(n_components=2, metric="precomputed", init="random", n_init=4,
              random_state=0, normalized_stress="auto", max_iter=300).fit_transform(D)

    # Divergence: distance from the cohort centroid should grow with age; the
    # cord-blood (age 0) samples should sit centrally.
    _centroid = emb.mean(axis=0)
    dist_centroid = np.linalg.norm(emb - _centroid, axis=1)
    _r, _p = spearmanr(ages, dist_centroid)
    mds_div = (float(_r), float(_p))
    cb_dist = float(dist_centroid[ages == 0].mean())
    old_dist = float(dist_centroid[ages >= 85].mean())
    return cb_dist, emb, mds_div, old_dist


@app.cell
def _(ages, emb, mds_div, plt):
    figm, axm = plt.subplots(figsize=(7.2, 5.4))
    _sc = axm.scatter(emb[:, 0], emb[:, 1], c=ages, cmap="viridis", s=60,
                      edgecolor="white", linewidth=0.6, zorder=3)
    # ring the cord-blood (age 0) samples to show they cluster centrally
    _cb = ages == 0
    axm.scatter(emb[_cb, 0], emb[_cb, 1], s=150, facecolors="none",
                edgecolors="#d1495b", linewidths=1.6, zorder=4,
                label="cord blood (age 0)")
    _cbar = figm.colorbar(_sc, ax=axm)
    _cbar.set_label("Age (years)")
    _r, _p = mds_div
    axm.set(xlabel="MDS 1", ylabel="MDS 2",
            title="Repertoires diverge from the cord-blood centre with age")
    axm.text(0.03, 0.03,
             f"dist. from centroid vs age\nSpearman r = {_r:.2f}, p = {_p:.1e}",
             transform=axm.transAxes, fontsize=9.5, va="bottom",
             bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#cccccc"))
    axm.legend(loc="upper right", fontsize=9, frameon=False)
    axm.spines[["top", "right"]].set_visible(False)
    figm.tight_layout()
    figm
    return


@app.cell
def _(cb_dist, mds_div, mo, old_dist):
    _r, _p = mds_div
    mo.md(
        f"""
        **Result.** The expected picture holds: cord-blood and young samples cluster
        centrally (a shared naive / public repertoire), and samples scatter to the
        periphery with age. Distance from the cohort centroid correlates with age at
        **Spearman r = {_r:.2f}, p = {_p:.1e}**; mean distance-from-centroid rises from
        **{cb_dist:.2f}** for cord blood (age 0) to **{old_dist:.2f}** for the oldest
        donors (age ≥ 85). The private, stochastic nature of age-associated clonal
        expansions is exactly what makes old repertoires idiosyncratic — and pulls
        them apart in overlap space.
        """
    )
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 6 · Clonal expansion vs age

        Diversity's flip side, computed straight from the canonical frame: the read
        share of the **top 10 clones**. Age-associated clonal expansions push this up
        in older donors.
        """
    )
    return


@app.cell
def _(OKABE, ages, counts, np, plt, reads, samples, spearmanr):
    _top10 = np.array([np.sort(counts[s])[::-1][:10].sum() / reads[s] * 100.0
                       for s in samples])
    _r, _p = spearmanr(ages, _top10)

    figc, axc = plt.subplots(figsize=(7.0, 4.4))
    axc.scatter(ages, _top10, s=44, color=OKABE["green"], edgecolor="white",
                linewidth=0.6, zorder=3)
    _b, _a0 = np.polyfit(ages, _top10, 1)
    _xs = np.array([ages.min(), ages.max()])
    axc.plot(_xs, _a0 + _b * _xs, color=OKABE["grey"], ls="--", lw=1.4)
    axc.set(xlabel="Age (years)", ylabel="Top-10-clone read share (%)",
            title="Repertoires become more clonal with age")
    axc.text(0.04, 0.94, f"Spearman r = {_r:.2f}\np = {_p:.1e}",
             transform=axc.transAxes, fontsize=10, va="top",
             bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#cccccc"))
    axc.spines[["top", "right"]].set_visible(False)
    figc.tight_layout()
    figc
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ---
        **Takeaway.** Across 79 healthy donors from cord blood to centenarians, the
        TRB repertoire loses diversity, gains clonality, and drifts apart with age —
        reproduced end-to-end from raw full-depth native vdjtools files with
        `vdjtools.io`, `vdjtools.stats` (iNEXT diversity, coverage standardization,
        rarefaction) and `vdjtools.overlap` (pairwise overlap → MDS). The bootstrap
        is md5-cached, so re-running this notebook touches the network zero times.
        """
    )
    return


if __name__ == "__main__":
    app.run()
