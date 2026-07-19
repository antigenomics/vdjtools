"""Repertoire overlap — exact, fuzzy, and similarity-aware — with vdjtools v2.

A marimo notebook. Launch it with

    marimo edit examples/overlap_similarity.py     # interactive editor
    marimo run  examples/overlap_similarity.py     # read-only served app

Two T-cell repertoires can share almost no *identical* receptors yet respond to the
same antigens, because they carry **near-variants** of the same specificities. Classic
exact-match overlap is blind to this; **similarity-aware** overlap (the TINA /
Leinster-Cobbold framework, `pᵀZq` with a CDR3-similarity kernel `Z`) sees it. This
notebook contrasts three lenses on the same cohort with `vdjtools.overlap` —

- **exact** (`overlap_metrics`, `Z = I`) — shared clonotype identity (CDR3+V+J),
- **fuzzy** (`fuzzy_overlap`, `Z = 1[≤1 substitution]`) — one-mismatch neighbours,
- **similarity-weighted** (`similarity_overlap`, `Z = exp(−BLOSUM62 / τ)`) — a smooth
  sequence-similarity kernel (TINA_w cosine),

builds all-pairs distance matrices (`pairwise_distances`), embeds the cohort
(`cluster_samples`, metric MDS), and runs a convergence test on one sample (`tcrnet`).

Data: a subset of the Britanova human TRB "Cord Blood to Centenarians" cohort (native
vdjtools format) auto-downloads from the HuggingFace dataset ``isalgo/airr_benchmark``
(folder ``vdjtools/``) into the gitignored ``examples/.data/overlap_nb/`` cache.
"""
import marimo

__generated_with = "0.23.14"
app = marimo.App(width="medium")


@app.cell
def _(mo):
    mo.md(
        r"""
        # Repertoire overlap, three ways — exact, fuzzy, and similarity-aware

        How similar are two T-cell repertoires? The textbook answer counts the
        **identical** clonotypes they share. But antigen recognition is degenerate: a
        single amino-acid substitution in a CDR3 usually preserves specificity, so two
        people responding to the same pathogen carry *near-variants* of the same
        receptors that an exact-match count scores as **zero overlap**.

        `vdjtools.overlap` measures overlap under a CDR3-**similarity kernel** `Z`
        (the TINA / Leinster-Cobbold form `pᵀZq`), and the choice of `Z` slides between
        three regimes:

        | Lens | kernel `Z` | what counts as "shared" |
        |---|---|---|
        | **exact** | `I` (identity) | identical CDR3 + V + J |
        | **fuzzy** | `1[≤ 1 substitution]` | one-mismatch neighbours |
        | **similarity** | `exp(−BLOSUM62 / τ)` | smooth sequence similarity |

        We compute all three on a slice of the Britanova **"Cord Blood to
        Centenarians"** TRB cohort, embed the samples with MDS, and show where the
        similarity kernel sees structure the exact count misses.
        """
    )
    return


@app.cell
def _():
    # --- imports & configuration ---
    from pathlib import Path

    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import polars as pl

    from vdjtools import io as vio
    from vdjtools.overlap import (cluster_samples, pairwise_distances,
                                  similarity_overlap, tcrnet)

    REPO_ID, HF_FOLDER = "isalgo/airr_benchmark", "vdjtools"
    N_SAMPLES = 12          # samples spanning the age range (kept small: O(n²) pairs)
    DEPTH = 3_000           # common downsampling depth (reads) before any comparison
    # BLOSUM62-penalty neighbourhood cutoff for the exp kernel: the default retains a
    # near-dense graph (~1e-3 weight floor) and is too slow for an all-pairs demo; 40
    # keeps meaningful weights (exp(−40/14) ≈ 0.06) with a bounded neighbourhood.
    MAX_PENALTY = 40
    OKABE = {"blue": "#0072B2", "vermillion": "#D55E00", "green": "#009E73",
             "orange": "#E69F00", "purple": "#CC79A7", "grey": "#8C8C8C"}
    return (DEPTH, HF_FOLDER, MAX_PENALTY, N_SAMPLES, OKABE, Path, REPO_ID,
            cluster_samples, mo, np, pairwise_distances, pl, plt,
            similarity_overlap, tcrnet, vio)


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 1 · Load a slice of the cohort, equalised to a common depth

        Overlap is depth-sensitive, so every sample is **downsampled to a common read
        depth** before any cross-sample comparison — the standard normalization. We
        take samples evenly spaced across age (cord blood → centenarian).
        """
    )
    return


@app.cell
def _(DEPTH, HF_FOLDER, N_SAMPLES, Path, REPO_ID, mo, np, pl, vio):
    _nb_dir = mo.notebook_dir() or Path.cwd()
    data_dir = _nb_dir / ".data" / "overlap_nb"
    data_dir.mkdir(parents=True, exist_ok=True)
    try:
        import huggingface_hub as _hub
    except ImportError:
        _hub = None
    mo.stop(_hub is None, mo.md("> Install the examples extra: `pip install \"vdjtools[examples]\"`."))

    _meta_path = _hub.hf_hub_download(REPO_ID, f"{HF_FOLDER}/metadata_aging.txt",
                                      repo_type="dataset")
    _meta = vio.read_metadata(_meta_path).with_columns(pl.col("age").cast(pl.Int64)).sort("age")
    # Evenly-spaced picks across the age range.
    _idx = np.linspace(0, _meta.height - 1, N_SAMPLES).round().astype(int)
    meta = _meta[_idx]
    _ids = meta["sample_id"].to_list()

    _root = Path(_hub.snapshot_download(
        REPO_ID, repo_type="dataset",
        allow_patterns=[f"{HF_FOLDER}/{s}.txt.gz" for s in _ids]))

    reps, ages = {}, {}
    _rng = np.random.default_rng(0)
    for _s, _a in zip(_ids, meta["age"].to_list()):
        _df = vio.read(_root / HF_FOLDER / f"{_s}.txt.gz", fmt="vdjtools")
        # Keep only coding CDR3s — the fuzzy/similarity kernels run on seqtree's aa
        # alphabet, which rejects the vdjtools non-coding markers (``*`` stop, ``_``
        # frameshift). (vdjtools.preprocess.filter_functional does the same.)
        _df = _df.filter(~pl.col("cdr3_aa").str.contains("[*_]"))
        # multinomial downsample to DEPTH on the clonotype frequencies
        _p = _df["duplicate_count"].to_numpy().astype(float)
        _p = _p / _p.sum()
        _n = _rng.multinomial(min(DEPTH, int(_df["duplicate_count"].sum())), _p)
        reps[_s] = _df.filter(pl.Series(_n > 0)).with_columns(
            pl.Series("duplicate_count", _n[_n > 0].astype("int64")))
        ages[_s] = _a
    age_arr = np.array([ages[s] for s in reps])
    mo.md(f"**{len(reps)} samples** (ages {age_arr.min()}–{age_arr.max()}), "
          f"each downsampled to ≤ {DEPTH:,} reads. Cache: `{data_dir}`")
    return age_arr, ages, meta, reps


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 2 · All-pairs overlap under each kernel

        `pairwise_distances` turns each overlap similarity into a distance:
        exact/fuzzy **F** → `−log10 F`, similarity cosine → `1 − cos`. The three
        matrices below are the same 15 samples seen through progressively more
        forgiving notions of "shared".
        """
    )
    return


@app.cell
def _(MAX_PENALTY, np, pairwise_distances, pl, reps, similarity_overlap):
    # Three distance matrices (heavy cell: O(n²) pairwise similarity searches).
    d_exact = pairwise_distances(reps, metric="F")                       # Z = I
    d_fuzzy = pairwise_distances(reps, metric="F", scope="1,0,0,1")      # Z = 1[≤1 sub]

    # Similarity-weighted (exp BLOSUM62 kernel). We call similarity_overlap directly with
    # a bounded neighbourhood (MAX_PENALTY) rather than pairwise_distances' default, which
    # is near-dense (~1e-3 floor) and impractically slow all-pairs. distance = 1 − cosine.
    def _sim_dist(reps, max_penalty):
        _names = list(reps)
        _n = len(_names)
        _m = np.zeros((_n, _n))
        for _i in range(_n):
            for _j in range(_i + 1, _n):
                _s = similarity_overlap(reps[_names[_i]], reps[_names[_j]], kernel="exp",
                                        max_penalty=max_penalty, metric="cosine")["similarity"]
                _m[_i, _j] = _m[_j, _i] = 1.0 - _s
        return pl.DataFrame({"sample": _names,
                             **{_names[_j]: _m[:, _j].tolist() for _j in range(_n)}})

    d_sim = _sim_dist(reps, MAX_PENALTY)
    return d_exact, d_fuzzy, d_sim


@app.cell
def _(age_arr, d_exact, d_fuzzy, d_sim, plt, reps):
    _names = list(reps)
    _order = age_arr.argsort()                     # order rows/cols by age
    figh, axh = plt.subplots(1, 3, figsize=(13.5, 4.4))
    for _ax, _d, _t in zip(axh, (d_exact, d_fuzzy, d_sim),
                           ("exact  (−log10 F)", "fuzzy  (−log10 F, ≤1 sub)",
                            "similarity  (1 − cosine)")):
        _m = _d.select(_names).to_numpy()[_order][:, _order]
        _im = _ax.imshow(_m, cmap="viridis_r")
        _ax.set_title(_t, fontsize=10)
        _ax.set_xticks([]); _ax.set_yticks([])
        _ax.set_xlabel("samples (age →)")
        figh.colorbar(_im, ax=_ax, fraction=0.046)
    axh[0].set_ylabel("samples (age →)")
    figh.suptitle("Pairwise repertoire distance — closer = darker", y=1.02)
    figh.tight_layout()
    figh
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 3 · Similarity sees overlap that exact-match misses

        Per sample pair, the exact overlap (`F`) vs the similarity-weighted cosine. The
        points sitting on the left edge (`F ≈ 0`) but well above zero on the y-axis are
        pairs with **almost no shared identical clonotypes yet substantial
        sequence-similarity overlap** — the near-variant sharing that only the kernel
        captures.
        """
    )
    return


@app.cell
def _(OKABE, d_exact, d_sim, np, plt, reps):
    _names = list(reps)
    _F = np.power(10.0, -d_exact.select(_names).to_numpy())        # F ≈ 10^(−dist)
    _C = 1.0 - d_sim.select(_names).to_numpy()                     # cosine = 1 − dist
    _iu = np.triu_indices(len(_names), k=1)
    _x, _y = _F[_iu], _C[_iu]

    figs, axs = plt.subplots(figsize=(6.6, 5.0))
    axs.scatter(_x, _y, s=34, color=OKABE["purple"], edgecolor="white", linewidth=0.5, zorder=3)
    axs.set(xlabel="exact overlap  F  (shared identical clonotypes)",
            ylabel="similarity-weighted cosine (TINA_w)",
            title="Sparse identical sharing, substantial sequence similarity")
    axs.text(0.5, 0.9, "pairs at F ≈ 0 with cosine ≫ 0:\nnear-variant overlap exact-match can't see",
             transform=axs.transAxes, fontsize=9, va="top", ha="center", color="#5c5c5c")
    axs.spines[["top", "right"]].set_visible(False)
    figs.tight_layout()
    figs
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 4 · Embed the cohort — one operator, three kernels

        Metric MDS of each distance matrix, coloured by age — the **same** `pᵀZq` operator
        read through `Z = I` → `1[≤1 sub]` → `exp(−BLOSUM62/τ)`. The number under each panel
        is the fraction of sample **pairs with any detectable overlap**: exact aa+V+J
        overlap of TRB is sparse (public sharing is rare, so many pairs share *nothing* and
        the exact embedding is dominated by all-or-nothing distances), a one-substitution
        kernel recovers more, and the smooth BLOSUM62 kernel places **every** pair on a
        graded continuum. (At this deliberately shallow, tractable depth the age structure
        itself is faint — the point here is the *sensitivity* of each lens, not an age
        clock; the full-depth aging signal is the `aging.py` notebook.)
        """
    )
    return


@app.cell
def _(ages, cluster_samples, d_exact, d_fuzzy, d_sim, np, pl, plt, reps):
    _meta = pl.DataFrame({"sample": list(reps), "age": [ages[s] for s in reps]})
    _names = list(reps)
    _iu = np.triu_indices(len(_names), k=1)

    figm, axm = plt.subplots(1, 3, figsize=(14, 4.7))
    _sc = None
    # thresholds below which a pair counts as "overlapping": F-distance < 9 (F>0),
    # cosine-distance < 1 (cos>0).
    for _ax, _d, _t, _thr in zip(axm, (d_exact, d_fuzzy, d_sim),
                                 ("exact  Z=I", "fuzzy  Z=1[≤1 sub]",
                                  "similarity  Z=exp(−BLOSUM62/τ)"), (8.99, 8.99, 0.999)):
        _emb = cluster_samples(_d, method="mds", metadata=_meta).sort("age")
        _sc = _ax.scatter(_emb["mds1"], _emb["mds2"], c=_emb["age"], cmap="viridis", s=80,
                          edgecolor="white", linewidth=0.5, zorder=3)
        _det = 100.0 * (_d.select(_names).to_numpy()[_iu] < _thr).mean()
        _ax.set_title(f"{_t}\n{_det:.0f}% of pairs share overlap", fontsize=10)
        _ax.set_xticks([]); _ax.set_yticks([])
        _ax.spines[["top", "right"]].set_visible(False)
    figm.colorbar(_sc, ax=axm, fraction=0.02, pad=0.02).set_label("Age (years)")
    figm.suptitle("Cohort MDS — the similarity kernel connects every repertoire", y=1.04)
    figm
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 5 · Convergence within one repertoire — TCRnet

        `tcrnet` asks a related question *within* a sample: which CDR3s sit in a dense
        1-substitution neighbourhood far above what a background generation model
        predicts (`vdjmatch` / `seqtree` e-values)? These enriched hubs are the
        convergent, likely antigen-driven clusters.
        """
    )
    return


@app.cell
def _(mo, reps, tcrnet):
    _one = next(iter(reps.values()))
    try:
        net = tcrnet(_one, locus="TRB", species="human")
        _enr = net.sort("p_enrichment").head(8)
        _cols = [c for c in ("cdr3_aa", "v_call", "n_target", "n_control",
                             "E", "p_enrichment") if c in _enr.columns]
        _out = mo.vstack([
            mo.md(f"**{net.filter(net['p_enrichment'] < 0.05).height}** CDR3s enriched "
                  f"(p<0.05) of {net.height:,} tested — top convergent hubs:"),
            _enr.select(_cols)])
    except ImportError:
        _out = mo.md("Install the `[overlap]` extra (`vdjmatch`, `seqtree`) for TCRnet.")
    _out
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ---
        **Takeaway.** The same cohort, three notions of "shared": exact identity is the
        strictest and sparsest; a one-substitution kernel already recovers near-variant
        sharing; and a smooth BLOSUM62 similarity kernel (TINA_w) is the most sensitive,
        registering overlap between repertoires that share no identical clonotype at
        all. All three are one call to `vdjtools.overlap`
        (`pairwise_distances` → `cluster_samples`), differing only in the kernel `Z`;
        `tcrnet` applies the same 1-substitution neighbourhood idea *within* a sample to
        surface convergent, antigen-driven clusters.
        """
    )
    return


if __name__ == "__main__":
    app.run()
