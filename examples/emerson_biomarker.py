"""Emerson 2017 biomarker discovery — an interactive vdjtools v2 walkthrough.

A marimo notebook. Launch it with

    marimo edit examples/emerson_biomarker.py     # interactive editor
    marimo run  examples/emerson_biomarker.py     # read-only served app

It reproduces the core of Emerson et al. (*Nat Genet* 2017, doi:10.1038/ng.3822):
an incidence-based **Fisher's-exact** screen for public TCRβ chains whose presence
across subjects is associated with **CMV serostatus** or an **HLA allele**, using
`vdjtools.biomarker.fisher_association`. The two knobs the notebook makes
interactive are exactly the two options of the method — the **V/J-match
requirement** and **exact vs 1-mismatch** CDR3 matching — plus the phenotype and
the significance threshold. Hits are validated live against a local **VDJdb** dump
by CMV epitope + HLA allele.

Data: a balanced subset of the 786-subject Emerson HIP cohort auto-downloads from
the HuggingFace dataset ``isalgo/airr_hip`` into the gitignored
``examples/.data/emerson_nb/`` cache (HF verifies integrity; a re-run fetches
nothing). The full-cohort, non-interactive version is ``examples/emerson_cmv_hla.py``.
"""
import marimo

__generated_with = "0.23.14"
app = marimo.App(width="medium")


@app.cell
def _(mo):
    mo.md(
        r"""
        # Finding CMV- and HLA-associated T-cell receptors — a vdjtools v2 walkthrough

        A person's T-cell repertoire is a record of what their immune system has
        seen. Chronic infections such as **cytomegalovirus (CMV)** drive clonal
        expansions of antigen-specific T cells, and because many of the responding
        **TCRβ** chains are *public* (shared across people with the same exposure and
        HLA type), their presence/absence across a cohort carries a statistical
        signature of the phenotype.

        Emerson et al. (2017) turned this into a screen: for every public TCRβ, a
        **2×2 Fisher's-exact test** of *how many subjects carry it* vs the phenotype.
        This notebook runs that screen with **vdjtools v2**
        (`vdjtools.biomarker.fisher_association`) on the Emerson HIP cohort, and lets
        you turn the method's two knobs live —

        - **V/J-match requirement** — is a "clonotype" its CDR3 alone, or CDR3 **+ V**,
          or CDR3 **+ V + J** (Emerson's definition)?
        - **exact vs 1-mismatch** — must two subjects share the *identical* CDR3, or
          does a single amino-acid substitution still count (grouping near-variants
          into a **metaclonotype**)?

        Every significant hit is checked against **VDJdb** — do the TCRβs our
        statistics flag as CMV-associated actually appear there as CMV-specific, with
        a matching HLA restriction?
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

    from vdjtools import io as vio
    from vdjtools.io import schema as S
    from vdjtools.biomarker import fisher_association

    REPO_ID = "isalgo/airr_hip"
    # Local VDJdb slim dump (antigenomics/vdjdb-db). Validation is skipped gracefully
    # if it is absent — the screen itself needs no VDJdb.
    VDJDB = Path("/Users/mikesh/vcs/code/vdjdb-db/database/vdjdb.slim.txt")

    # Balanced subset size (half CMV+, half CMV−). Kept modest so the screen recomputes
    # interactively; the full 786-subject run is examples/emerson_cmv_hla.py.
    N_SUBJECTS = 400

    # Okabe–Ito colorblind-safe palette.
    OKABE = {"blue": "#0072B2", "vermillion": "#D55E00", "green": "#009E73",
             "orange": "#E69F00", "purple": "#CC79A7", "grey": "#8C8C8C"}
    return (N_SUBJECTS, OKABE, Path, REPO_ID, S, VDJDB, fisher_association,
            mo, np, pl, plt, vio)


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 1 · Data bootstrap

        The `metadata.txt` sheet (786 subjects) carries each subject's **CMV
        serostatus** (`+`/`−`/`NA`) and **2-digit HLA-A/B typing**. We take a balanced
        subset (half CMV+, half CMV−) and pull just those subjects' repertoires
        (`corr/HIP#####.txt.gz`, VDJtools format) from HuggingFace into the gitignored
        `examples/.data/emerson_nb/` cache. A re-run re-uses the cache and the
        already-ingested Parquet cohort.

        > **Note on the split on TAB, not comma** — the `race` metadata field itself
        > contains commas, so the sheet must be parsed tab-delimited.
        """
    )
    return


@app.cell
def _(N_SUBJECTS, Path, REPO_ID, mo, pl):
    _nb_dir = mo.notebook_dir() or Path.cwd()
    data_dir = _nb_dir / ".data" / "emerson_nb"
    data_dir.mkdir(parents=True, exist_ok=True)

    try:
        import huggingface_hub as _hub
    except ImportError:
        _hub = None
    mo.stop(
        _hub is None,
        mo.md(
            """
            > **`huggingface_hub` is not installed.** Install the examples extra and
            > re-run: `pip install "vdjtools[examples]"`.
            """
        ),
    )

    _meta_path = _hub.hf_hub_download(REPO_ID, "metadata.txt", repo_type="dataset")
    meta_all = pl.read_csv(_meta_path, separator="\t", infer_schema_length=0)  # TAB — race has commas
    _known = meta_all.filter(pl.col("cmv").is_in(["+", "-"]))
    meta = pl.concat([_known.filter(pl.col("cmv") == c).head(N_SUBJECTS // 2)
                      for c in ("+", "-")])
    _ids = meta["sample_id"].to_list()

    _root = Path(_hub.snapshot_download(
        REPO_ID, repo_type="dataset",
        allow_patterns=[f"corr/{s}.txt.gz" for s in _ids]))
    corr_dir = _root / "corr"
    mo.md(f"**{len(_ids)} subjects** — "
          f"{(meta['cmv'] == '+').sum()} CMV+, {(meta['cmv'] == '-').sum()} CMV−. "
          f"Cache: `{data_dir}`")
    return corr_dir, data_dir, meta


@app.cell
def _(corr_dir, data_dir, meta, mo, pl, vio):
    # Stream the per-sample VDJtools files into a hive-partitioned Parquet cohort once
    # (one sample in RAM at a time), then scan it as a single out-of-core LazyFrame.
    cohort = data_dir / "cohort"
    if not any(cohort.glob("sample_id=*/*.parquet")):
        vio.ingest_cohort(meta.select("sample_id", "cmv", "hla"), corr_dir, cohort,
                          sample_col="sample_id", file_template="{sample}.txt.gz",
                          fmt="vdjtools")
    lf = vio.scan_cohort(cohort, join_metadata=False)

    # Per-subject binary phenotypes: CMV +/− and HLA-A*02 present/absent (2-digit).
    cmv_ph = meta.select(
        "sample_id",
        pl.when(pl.col("cmv") == "+").then(True)
          .when(pl.col("cmv") == "-").then(False).otherwise(None).alias("pheno"))
    a02_ph = meta.select(
        "sample_id",
        pl.when(pl.col("hla").is_null() | pl.col("hla").is_in(["", "NA"])).then(None)
          .otherwise(pl.col("hla").str.contains(r"HLA-A\*02")).alias("pheno"))
    mo.md(f"Cohort ingested → `{cohort}` · phenotypes built (CMV, HLA-A*02).")
    return a02_ph, cmv_ph, lf


@app.cell
def _(VDJDB, mo, pl):
    # VDJdb CMV validation set (human TRB), if the local dump is present.
    if VDJDB.exists():
        _v = pl.read_csv(VDJDB, separator="\t", infer_schema_length=0)
        # One row per CDR3, epitopes/MHCs it was reported against joined — so the
        # validation join stays 1:1 with our hits.
        vdjdb_cmv = (_v.filter((pl.col("gene") == "TRB")
                               & (pl.col("species") == "HomoSapiens")
                               & pl.col("antigen.species").str.contains("CMV"))
                     .group_by(pl.col("cdr3")).agg(
                         pl.col("antigen.epitope").unique().sort().str.join(", ").alias("epitope"),
                         pl.col("mhc.a").unique().sort().str.join(", ").alias("mhc")))
        _msg = f"VDJdb CMV reference: **{vdjdb_cmv.height} human TRB CDR3s**."
    else:
        vdjdb_cmv = None
        _msg = f"> VDJdb not found at `{VDJDB}` — validation is skipped."
    mo.md(_msg)
    return (vdjdb_cmv,)


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 2 · The screen — turn the knobs

        `fisher_association` computes, for every feature present in **≥ `min incidence`**
        subjects, the 2×2 incidence table against the phenotype and its Fisher p-value
        (one-tailed *enrichment* for CMV, two-tailed for HLA), plus a Benjamini-Hochberg
        `q`. Set the four controls below — the screen re-runs when the phenotype,
        matching, or min-incidence change; the **threshold** only re-draws.

        > **On significance.** With millions of features, BH-`q` is extremely strict at
        > this subset size — so, like Emerson, the primary line is a **nominal
        > p-threshold** (the paper used P < 1e-4, FDR ≈ 0.14). `q` is shown as the
        > conservative secondary.
        """
    )
    return


@app.cell
def _(mo):
    phenotype = mo.ui.dropdown(
        {"CMV (enrichment, one-tailed)": "cmv",
         "HLA-A*02 (two-sided)": "hla"}, value="CMV (enrichment, one-tailed)",
        label="Phenotype")
    match = mo.ui.dropdown(
        {"exact CDR3": "exact", "1 mismatch (metaclonotype)": "1mm"},
        value="exact CDR3", label="CDR3 matching")
    key = mo.ui.dropdown(
        {"CDR3 + V + J (Emerson)": "vj", "CDR3 + V": "v", "CDR3 only": "cdr3"},
        value="CDR3 + V + J (Emerson)", label="Clonotype key (V/J match)")
    min_inc = mo.ui.slider(2, 10, value=2, label="min incidence (subjects)")
    logp = mo.ui.slider(2.0, 8.0, value=4.0, step=0.5,
                        label="−log10 p threshold")
    mo.hstack([phenotype, match, key, min_inc, logp], justify="start", gap=1.5)
    return key, logp, match, min_inc, phenotype


@app.cell
def _(S, a02_ph, cmv_ph, fisher_association, key, lf, match, min_inc, phenotype):
    _KEYS = {"vj": (S.JUNCTION_AA, S.V_CALL, S.J_CALL),
             "v": (S.JUNCTION_AA, S.V_CALL), "cdr3": (S.JUNCTION_AA,)}
    _is_cmv = phenotype.value == "cmv"
    res = fisher_association(
        lf, cmv_ph if _is_cmv else a02_ph, pheno_col="pheno",
        key=_KEYS[key.value], match=match.value,
        alternative="greater" if _is_cmv else "two-sided",
        min_incidence=min_inc.value)
    return (res,)


@app.cell
def _(pl, res):
    # Reduce for rendering: keep every interesting feature (p < 0.01) plus a background
    # sample so the volcano stays legible at millions of features.
    _sig = res.filter(pl.col("p_value") < 0.01)
    _bg = res.filter(pl.col("p_value") >= 0.01)
    if _bg.height > 15000:
        _bg = _bg.sample(15000, seed=0)
    plot_df = pl.concat([_sig, _bg])
    top = res.head(200)                       # top hits by p-value for validation
    return plot_df, top


@app.cell
def _(OKABE, logp, np, phenotype, plot_df, plt, vdjdb_cmv):
    _x = plot_df["log2_or"].to_numpy()
    _y = np.clip(-np.log10(plot_df["p_value"].to_numpy()), 0, 40)
    _hit = _y >= logp.value
    _val = (plot_df["junction_aa"].is_in(list(vdjdb_cmv["cdr3"])).to_numpy()
            if vdjdb_cmv is not None else np.zeros(len(_x), bool)) & _hit

    figv, axv = plt.subplots(figsize=(7.6, 5.2))
    axv.scatter(_x[~_hit], _y[~_hit], s=5, c="#c8ccd4", alpha=0.4, lw=0)
    axv.scatter(_x[_hit], _y[_hit], s=12, c=OKABE["vermillion"], alpha=0.75, lw=0,
                label=f"p < 1e-{logp.value:g}")
    if _val.any():
        axv.scatter(_x[_val], _y[_val], s=44, facecolors="none",
                    edgecolors=OKABE["blue"], linewidths=1.4, label="in VDJdb-CMV")
    axv.axhline(logp.value, ls="--", c="k", lw=0.6, alpha=0.5)
    axv.axvline(0, ls="-", c="k", lw=0.6, alpha=0.3)
    axv.set(xlabel="log2 odds ratio  (enriched →)", ylabel="−log10 p-value",
            title=f"{phenotype.selected_key} — associated TCRβ")
    axv.legend(fontsize=9, frameon=False, loc="upper left")
    axv.spines[["top", "right"]].set_visible(False)
    figv.tight_layout()
    figv
    return


@app.cell
def _(logp, mo, pl, res, top, vdjdb_cmv):
    _thr = 10 ** (-logp.value)
    n_hits = res.filter(pl.col("p_value") < _thr).height          # true count over all features
    n_q = res.filter(pl.col("q_value") < 0.05).height
    disp = top.filter(pl.col("p_value") < _thr)                   # top hits (≤200) to show
    _cols = (["junction_aa"] + [c for c in ("v_call", "j_call") if c in disp.columns]
             + ["incidence", "n_pos_present", "odds_ratio", "p_value", "q_value"])
    if vdjdb_cmv is not None:
        table = (disp.join(vdjdb_cmv, left_on="junction_aa", right_on="cdr3", how="left")
                 .select(_cols + ["epitope", "mhc"]).sort("p_value"))
        n_val = disp.join(vdjdb_cmv, left_on="junction_aa", right_on="cdr3",
                          how="semi")["junction_aa"].n_unique()
        _val_msg = f" · **{n_val} matched to VDJdb-CMV** (exact CDR3)"
    else:
        table = disp.select(_cols).sort("p_value")
        _val_msg = " · VDJdb validation skipped"
    _summary = mo.md(
        f"**{n_hits:,} features** at p < {_thr:.0e} "
        f"({n_q} also pass BH q<0.05) of {res.height:,} tested{_val_msg}. "
        f"Showing the top {disp.height}; VDJdb-matched rows carry the epitope + MHC "
        f"they were reported against.")
    mo.vstack([_summary, table])
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 3 · 1-mismatch to VDJdb (bonus)

        The exact-match column above undercounts biology: a CMV-driving TCRβ and a
        VDJdb entry can differ by a single residue. Switch **CDR3 matching → 1
        mismatch** to group near-variants into metaclonotypes; the cell below also asks
        the looser question directly — how many of our top hits have a VDJdb-CMV
        neighbour **within one substitution** (`vdjmatch.cluster.overlap`)?
        """
    )
    return


@app.cell
def _(mo, top, vdjdb_cmv):
    if vdjdb_cmv is None:
        _msg = "VDJdb not available — skipped."
    else:
        try:
            import vdjmatch.cluster as _vc
            _ours = top["junction_aa"].unique().to_list()
            _ref = vdjdb_cmv["cdr3"].unique().to_list()
            _pairs = _vc.overlap(_ours, _ref, scope="1,0,0,1")
            _n1 = _pairs["a_idx"].n_unique() if _pairs.height else 0
            _msg = (f"**{_n1} / {len(_ours)}** top hits have a VDJdb-CMV neighbour "
                    f"within ≤ 1 substitution (vs the smaller exact-match count above).")
        except ImportError:
            _msg = "Install the `[overlap]` extra (`vdjmatch`) for the 1-mismatch match."
    mo.md(_msg)
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ---
        **Takeaway.** A pure incidence Fisher screen — no antigen assay, no structure —
        rediscovers known CMV-specific public TCRβs (e.g. `CASSLAPGATNEKLFF` ↔ the
        pp65 epitope **NLVPMVATV**, HLA-A\*02:01) and their HLA restriction, straight
        from raw repertoires. Widening the clonotype key from CDR3-only to CDR3+V+J
        sharpens specificity; allowing a single mismatch recovers more of the VDJdb
        overlap. Built on `vdjtools.io` (streamed Parquet cohort),
        `vdjtools.biomarker.fisher_association` (vectorised hypergeometric Fisher over
        millions of features), and `vdjtools.overlap` (the native 1-mismatch matcher).
        The full 786-subject screen is `examples/emerson_cmv_hla.py`.
        """
    )
    return


if __name__ == "__main__":
    app.run()
