# vdjtools — how the human TCR repertoire ages.
# A reactive marimo app over the full-depth Britanova "Cord Blood to Centenarians" TRB cohort
# (`isalgo/airr_benchmark`, folder `vdjtools/`, 78 donors aged 0-103). It loads the cohort once
# and reads off every aging signal three ways:
#   * cohort-STREAMING stats — diversity_cohort, the clone-size distribution, and the CDR3
#     spectratype, each a single streamed group_by over a `scan_cohort` LazyFrame (flat memory);
#   * coverage-standardized iNEXT diversity + rarefaction (the depth-confound-free headline);
#   * pairwise repertoire overlap -> MDS (repertoires diverge from a cord-blood centre with age).
# Data auto-loads from HuggingFace, preferring a local `~/hf/` or `./` copy. Run with:
#     marimo edit examples/aging.py
import marimo

__generated_with = "0.23.14"
app = marimo.App(width="medium")


@app.cell
def _(mo):
    mo.md(
        r"""
        # Ageing of the human TCR repertoire — a vdjtools walkthrough

        As we age the naive T-cell pool contracts and a handful of clones expand, so the T-cell
        receptor (TCR) repertoire becomes **less diverse and more clonal** — and, because those
        expansions are largely *private* (stochastic, person-specific), individual repertoires
        **drift apart** with age. This notebook reads all three signals off the full-depth
        **Britanova "Cord Blood to Centenarians" TRB cohort** (`isalgo/airr_benchmark`, 78 donors,
        ages **0 → 103**) with vdjtools, three complementary ways:

        - **Cohort-streaming stats** — the cohort is ingested once into a hive-partitioned Parquet
          dataset, then `stats.diversity_cohort`, the clone-size distribution, and the
          `spectratype` are each a *single* streamed `group_by` over `io.scan_cohort` (memory stays
          flat however many samples you scan).
        - **Coverage-standardized diversity** — the rigorous headline: Hill diversity read at a
          common *sample coverage* (iNEXT), which removes the sequencing-depth confound, plus the
          classic rarefaction/extrapolation curves.
        - **Repertoire divergence** — pairwise exact-match overlap (`vdjtools.overlap`) turned into
          a distance and embedded with **MDS**: cord-blood samples cluster centrally, older donors
          scatter to the periphery.

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
    from sklearn.manifold import MDS

    from vdjtools.io import ingest_cohort, read_metadata, scan_cohort
    from vdjtools.io.schema import COUNT, J_CALL, JUNCTION_AA, V_CALL
    from vdjtools.overlap import overlap_metrics
    from vdjtools.stats import (
        diversity_cohort,
        estimate_d,
        inext_batch,
        rarefaction,
        sample_coverage,
        spectratype,
    )

    REPO = "isalgo/airr_benchmark"
    HF_FOLDER = "vdjtools"                 # the full-depth aging cohort (not vdjtools_lite/)
    KEY = (JUNCTION_AA, V_CALL, J_CALL)    # exact overlap key: CDR3aa + V + J (== overlap DEFAULT_KEY)
    OKABE = {"blue": "#0072B2", "vermillion": "#D55E00", "green": "#009E73",
             "orange": "#E69F00", "purple": "#CC79A7", "grey": "#8C8C8C"}

    def local_base(mo_dir):
        """Return a directory holding the aging files, preferring a local mirror."""
        for root in (mo_dir, Path.cwd(), Path.home() / "hf" / "airr_benchmark"):
            if (root / HF_FOLDER / "metadata_aging.txt").exists():
                return root / HF_FOLDER
        return None

    return (COUNT, HF_FOLDER, KEY, MDS, OKABE, Path, REPO, diversity_cohort,
            estimate_d, inext_batch, ingest_cohort, local_base, mo, np,
            overlap_metrics, pl, rarefaction, read_metadata, sample_coverage,
            scan_cohort, spearmanr, spectratype)


@app.cell
def _(mo):
    n_samples = mo.ui.slider(12, 78, value=24, step=2,
                             label="samples across the age range (overlap is O(n²))")
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
def _(COUNT, KEY, lf, meta, np, pl):
    # --- eager per-sample structures for the iNEXT / rarefaction / overlap panels ---
    # One streamed collect of just the columns those panels need, partitioned per sample.
    _full = lf.select("sample_id", *KEY, COUNT).collect(engine="streaming")
    _parts = {df["sample_id"][0]: df for df in _full.partition_by("sample_id")}

    samples = meta["sample_id"].to_list()         # age-ordered (meta is sorted by age)
    age_arr = meta["age"].to_numpy()
    counts, reads, n_clones, _keytab = {}, {}, {}, {}
    for _s in samples:
        _df = _parts[_s]
        _c = _df[COUNT].to_numpy().astype("int64")
        counts[_s] = _c
        reads[_s] = int(_c.sum())
        n_clones[_s] = _df.height
        _keytab[_s] = _df.group_by(list(KEY), maintain_order=True).agg(
            pl.col(COUNT).sum().alias("c"))

    # Overlap is depth-sensitive: downsample every sample to a COMMON depth (cohort min reads)
    # on the CDR3aa+V+J key before any cross-sample comparison.
    DEPTH = min(reads.values())
    _rng = np.random.default_rng(0)
    ds = {}
    for _s in samples:
        _g = _keytab[_s]
        _p = _g["c"].to_numpy().astype("float64")
        _p = _p / _p.sum()
        _new = _rng.multinomial(DEPTH, _p)
        _keep = _new > 0
        ds[_s] = (_g.filter(pl.Series(_keep))
                  .with_columns(pl.Series(COUNT, _new[_keep].astype("int64")))
                  .select([*KEY, COUNT]))
    return DEPTH, age_arr, counts, ds, n_clones, reads, samples


@app.cell
def _(DEPTH, age_arr, mo, n_clones, pl, reads, samples):
    summary = pl.DataFrame({
        "sample_id": samples,
        "age": age_arr,
        "reads": [reads[s] for s in samples],
        "n_clonotypes": [n_clones[s] for s in samples],
    }).sort("age")
    mo.md(f"Common downsampling depth (cohort min reads) = **{DEPTH:,}**.")
    summary
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 1 · Diversity falls with age

        First the fast **streamed** view: `stats.diversity_cohort` computes the whole cohort's
        diversity table from the count-frequency spectrum in one pass over the `scan_cohort`
        LazyFrame — memory flat regardless of cohort size.
        """
    )
    return


@app.cell
def _(OKABE, ages, diversity_cohort, lf, np, plt, spearmanr):
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
            title=f"Streamed diversity_cohort: Shannon declines with age "
                  f"(Spearman ρ={_r:.2f}, p={_p:.1g})")
    ax1.spines[["top", "right"]].set_visible(False)
    fig1.tight_layout()
    fig1
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ### …now depth-standardized (iNEXT)

        Comparing raw diversity across samples is a trap: a deeper or more completely sampled
        repertoire looks more diverse just from more sequencing. The principled fix is **coverage
        standardization** — evaluate every sample's Hill diversity at a common *sample coverage* Ĉ
        (the fraction of the repertoire the reads represent), following iNEXT. Here the full-depth
        samples span a wide Ĉ, so we fix the common coverage at **C = min Ĉ** and read off Hill
        diversity there for orders q = 0 (richness), 1 (Shannon) and 2 (Simpson).
        """
    )
    return


@app.cell
def _(age_arr, counts, estimate_d, inext_batch, pl, sample_coverage, samples, spearmanr):
    import numpy as _np
    _cov = _np.array([sample_coverage(counts[s]) for s in samples])
    cov_df = pl.DataFrame({"sample": samples, "age": age_arr, "coverage": _cov})
    C = float(_cov.min())

    # Coverage-standardized Hill diversity (q = 0, 1, 2) with (light) bootstrap CIs.
    _rows = []
    for _i, _s in enumerate(samples):
        _d = estimate_d(counts[_s], base="coverage", level=C, q=(0, 1, 2),
                        se=True, nboot=10, seed=0)
        for _r in _d.iter_rows(named=True):
            _rows.append({"sample": _s, "age": int(age_arr[_i]), "order_q": _r["order_q"],
                          "qD": _r["qD"], "qD_lo": _r["qD_lo"], "qD_hi": _r["qD_hi"]})
    div_std = pl.DataFrame(_rows)

    spear = {}
    for _q in (0, 1, 2):
        _sub = div_std.filter(pl.col("order_q") == _q)
        _rr, _pp = spearmanr(_sub["age"].to_numpy(), _sub["qD"].to_numpy())
        spear[_q] = (float(_rr), float(_pp))

    # Robustness: raw (size-based) Shannon from the fast native batch engine.
    _batch = inext_batch([counts[s] for s in samples], q=(1,), knots=8, se=False, seed=0)
    _obs = _batch.filter(pl.col("method") == "observed").sort("sample")
    _rr, _pp = spearmanr(age_arr, _obs["qD"].to_numpy())
    raw_shannon = (float(_rr), float(_pp))
    return C, cov_df, div_std, raw_shannon, spear


@app.cell
def _(C, OKABE, cov_df, div_std, np, pl, plt, spear):
    _cov = cov_df.sort("age")
    _q1 = div_std.filter(pl.col("order_q") == 1).sort("age")
    _age = _q1["age"].to_numpy()
    _qd = _q1["qD"].to_numpy()
    _err = np.vstack([_qd - _q1["qD_lo"].to_numpy(), _q1["qD_hi"].to_numpy() - _qd])

    fig2, (axA, axB) = plt.subplots(1, 2, figsize=(11.5, 4.4))
    # (A) Motivation: sample coverage varies widely across the cohort.
    axA.scatter(_cov["age"], _cov["coverage"], s=34, color=OKABE["blue"],
                edgecolor="white", linewidth=0.5, zorder=3)
    axA.axhline(C, color=OKABE["grey"], ls="--", lw=1)
    axA.text(0.98, C, f"  common C = {C:.2f}", va="bottom", ha="right",
             transform=axA.get_yaxis_transform(), color="#5c5c5c", fontsize=9)
    axA.set(xlabel="Age (years)", ylabel="Sample coverage  Ĉ(n)",
            title="Sample completeness varies across donors")
    axA.spines[["top", "right"]].set_visible(False)
    # (B) Headline: coverage-standardized Shannon diversity vs age, with CIs.
    axB.errorbar(_age, _qd, yerr=_err, fmt="o", ms=5, color=OKABE["vermillion"],
                 ecolor="#D55E0055", elinewidth=1.0, capsize=2, mec="white", mew=0.5, zorder=3)
    axB.set_yscale("log")
    axB.set(xlabel="Age (years)", ylabel="Coverage-standardized Shannon  ¹D  (log)",
            title="TCR diversity declines with age")
    _r, _p = spear[1]
    axB.text(0.04, 0.06, f"Spearman r = {_r:.2f}\np = {_p:.1e}", transform=axB.transAxes,
             fontsize=10, va="bottom", bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#ccc"))
    axB.spines[["top", "right"]].set_visible(False)
    fig2.tight_layout()
    fig2
    return


@app.cell
def _(mo, raw_shannon, spear):
    _r0, _p0 = spear[0]
    _r1, _p1 = spear[1]
    _r2, _p2 = spear[2]
    _rr, _rp = raw_shannon
    mo.md(
        f"""
        **Result.** Coverage-standardized diversity falls with age across every Hill order —
        richness (q=0): r = {_r0:.2f}, p = {_p0:.1e}; Shannon (q=1): r = {_r1:.2f}, p = {_p1:.1e};
        Simpson (q=2): r = {_r2:.2f}, p = {_p2:.1e}. The trend is robust to the standardization
        choice: the raw, size-based Shannon from the native `inext_batch` gives r = {_rr:.2f}
        (p = {_rp:.1e}). Older repertoires support markedly fewer *effective* clones.
        """
    )
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 2 · Rarefaction / extrapolation curves

        The classic vdjtools rarefaction plot: clonotype **richness** (Hill q=0) as a function of
        sampling depth *m* (`rarefaction(sample, q=0)`), with bootstrap bands, for the youngest
        (cord-blood) vs the oldest donors. The depth axis is log-scaled — full-depth samples differ
        many-fold in sequencing depth, exactly the confound the coverage standardization above removes.
        """
    )
    return


@app.cell
def _(OKABE, age_arr, counts, np, plt, rarefaction, reads, samples):
    # samples is age-ordered: first two = youngest (cord blood), last two = oldest.
    _picks = ([(samples[i], int(age_arr[i]), OKABE["blue"], "cord blood") for i in (0, 1)]
              + [(samples[i], int(age_arr[i]), OKABE["vermillion"], "old") for i in (-2, -1)])
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
        axr.text(_m[-1], _qd[-1], f"  {_s} ({_a}y, {_grp})", color=_col, va="center", fontsize=9)
    axr.set_xscale("log")
    axr.set(xlabel="Sampling depth  m  (reads, log)", ylabel="Clonotype richness  ⁰D",
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
        ## 3 · The clone-size distribution shifts: naive → clonal

        Diversity's flip side. The **singleton fraction** (clonotypes seen once — the naive-rich
        tail) falls, while the **hyperexpanded fraction** (reads carried by clones above 1%
        frequency) rises — both one streamed `group_by(sample_id)` over the cohort. Alongside,
        the read share of the **top 10 clones**, straight from the per-sample count vectors.
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
    fig3, ax3 = plt.subplots(figsize=(7.0, 4.2))
    ax3.scatter(_age, size["singleton_frac"], s=30, color=OKABE["blue"], label="singleton fraction")
    ax3.scatter(_age, size["hyperexpanded_frac"], s=30, color=OKABE["vermillion"],
                marker="s", label="hyperexpanded read fraction (>1%)")
    _r1 = spearmanr(_age, size["singleton_frac"].to_numpy())[0]
    _r2 = spearmanr(_age, size["hyperexpanded_frac"].to_numpy())[0]
    ax3.set(xlabel="Age (years)", ylabel="Fraction",
            title=f"Naive→clonal shift  (singleton ρ={_r1:.2f}, hyperexpanded ρ={_r2:.2f})")
    ax3.legend(frameon=False, fontsize=9)
    ax3.spines[["top", "right"]].set_visible(False)
    fig3.tight_layout()
    fig3
    return


@app.cell
def _(OKABE, age_arr, counts, np, plt, reads, samples, spearmanr):
    _top10 = np.array([np.sort(counts[s])[::-1][:10].sum() / reads[s] * 100.0 for s in samples])
    _r, _p = spearmanr(age_arr, _top10)
    figc, axc = plt.subplots(figsize=(7.0, 4.2))
    axc.scatter(age_arr, _top10, s=44, color=OKABE["green"], edgecolor="white",
                linewidth=0.6, zorder=3)
    _b, _a0 = np.polyfit(age_arr, _top10, 1)
    _xs = np.array([age_arr.min(), age_arr.max()])
    axc.plot(_xs, _a0 + _b * _xs, color=OKABE["grey"], ls="--", lw=1.4)
    axc.set(xlabel="Age (years)", ylabel="Top-10-clone read share (%)",
            title="Repertoires become more clonal with age")
    axc.text(0.04, 0.94, f"Spearman r = {_r:.2f}\np = {_p:.1e}", transform=axc.transAxes,
             fontsize=10, va="top", bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#ccc"))
    axc.spines[["top", "right"]].set_visible(False)
    figc.tight_layout()
    figc
    return


@app.cell
def _(mo):
    mo.md(r"""## 4 · CDR3-length spectratype by age group (one fused pass)""")
    return


@app.cell
def _(OKABE, ages, lf, pl, plt, spectratype):
    sp = (spectratype(lf, by=["sample_id"], by_locus=False, weight="freq")
          .collect(engine="streaming").join(ages, on="sample_id"))
    _grp = sp.with_columns(pl.when(pl.col("age") < 20).then(pl.lit("0–19"))
                           .when(pl.col("age") < 60).then(pl.lit("20–59"))
                           .otherwise(pl.lit("60+")).alias("agegrp"))
    _agg = _grp.group_by("agegrp", "length").agg(pl.col("weight").mean()).sort("length")
    fig4, ax4 = plt.subplots(figsize=(7.0, 4.2))
    for _g, _c in zip(["0–19", "20–59", "60+"], [OKABE["green"], OKABE["orange"], OKABE["purple"]]):
        _s = _agg.filter(pl.col("agegrp") == _g).sort("length")
        ax4.plot(_s["length"], _s["weight"], "-o", ms=3, color=_c, label=_g)
    ax4.set(xlabel="CDR3 amino-acid length", ylabel="Mean clonotype-frequency weight",
            title="CDR3-length spectratype by age group")
    ax4.legend(title="Age", frameon=False)
    ax4.spines[["top", "right"]].set_visible(False)
    fig4.tight_layout()
    fig4
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 5 · Repertoires diverge with age — overlap MDS

        The headline result of this cohort: repertoires **diverge from a central cord-blood cluster
        outward with age**, driven by *private* clonal expansions. We measure it with pairwise
        exact-match **overlap** (`vdjtools.overlap`, the frequency-overlap **F** metric on the
        CDR3aa+V+J key) after equal-depth downsampling (§ setup), turn overlap into a distance
        (`d = -log₁₀ F`), embed every sample with **metric MDS**, and colour the map by age.
        """
    )
    return


@app.cell
def _(MDS, age_arr, ds, np, overlap_metrics, samples, spearmanr):
    _n = len(samples)
    F = np.zeros((_n, _n))
    for _i in range(_n):
        for _j in range(_i + 1, _n):
            _f = overlap_metrics(ds[samples[_i]], ds[samples[_j]])["F"]
            F[_i, _j] = F[_j, _i] = _f
    with np.errstate(divide="ignore"):
        D = -np.log10(F)
    _finite_max = D[np.isfinite(D)].max()
    D[~np.isfinite(D)] = _finite_max * 1.1
    np.fill_diagonal(D, 0.0)
    D = (D + D.T) / 2.0
    emb = MDS(n_components=2, metric="precomputed", init="random", n_init=4,
              random_state=0, normalized_stress="auto", max_iter=300).fit_transform(D)
    _centroid = emb.mean(axis=0)
    dist_centroid = np.linalg.norm(emb - _centroid, axis=1)
    _r, _p = spearmanr(age_arr, dist_centroid)
    mds_div = (float(_r), float(_p))
    cb_dist = float(dist_centroid[age_arr == age_arr.min()].mean())
    old_dist = float(dist_centroid[age_arr >= 85].mean()) if (age_arr >= 85).any() else float("nan")
    return cb_dist, emb, mds_div, old_dist


@app.cell
def _(age_arr, emb, mds_div, plt):
    figm, axm = plt.subplots(figsize=(7.2, 5.4))
    _sc = axm.scatter(emb[:, 0], emb[:, 1], c=age_arr, cmap="viridis", s=60,
                      edgecolor="white", linewidth=0.6, zorder=3)
    _young = age_arr == age_arr.min()
    axm.scatter(emb[_young, 0], emb[_young, 1], s=150, facecolors="none",
                edgecolors="#d1495b", linewidths=1.6, zorder=4, label="youngest")
    _cbar = figm.colorbar(_sc, ax=axm)
    _cbar.set_label("Age (years)")
    _r, _p = mds_div
    axm.set(xlabel="MDS 1", ylabel="MDS 2",
            title="Repertoires diverge from the young centre with age")
    axm.text(0.03, 0.03, f"dist. from centroid vs age\nSpearman r = {_r:.2f}, p = {_p:.1e}",
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
        **Result.** Cord-blood and young samples cluster centrally (a shared naive / public
        repertoire) and samples scatter to the periphery with age. Distance from the cohort
        centroid correlates with age at **Spearman r = {_r:.2f}, p = {_p:.1e}**; mean
        distance-from-centroid rises from **{cb_dist:.2f}** for the youngest donors to
        **{old_dist:.2f}** for the oldest (age ≥ 85). The private, stochastic nature of
        age-associated clonal expansions is exactly what makes old repertoires idiosyncratic.
        """
    )
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ---
        **Takeaway.** Across dozens of healthy donors from cord blood to centenarians, the TRB
        repertoire loses diversity, gains clonality, and drifts apart with age — reproduced
        end-to-end from raw full-depth files. The streaming panels (`diversity_cohort`, the
        clone-size `group_by`, the fused `spectratype`) keep peak memory independent of cohort size
        via `io.ingest_cohort` / `scan_cohort`; the coverage-standardized iNEXT diversity,
        rarefaction, and pairwise-overlap→MDS panels add the depth-confound-free, per-sample view.
        """
    )
    return


if __name__ == "__main__":
    app.run()
