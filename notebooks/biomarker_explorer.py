# vdjtools — interactive biomarker-association & co-occurrence explorer.
# Reactive marimo app over the Emerson HIP cohort (isalgo/airr_hip, 786 TCRβ subjects): pick a
# condition (CMV / HLA-A*02 / CMV|HLA-A*02 CMH), a statistical test (Fisher / Chi² / Bayesian
# BF / log-odds / permutation), and a match scope (CDR3[+V[+J]], exact or 1-mismatch); see the
# volcano with a live VDJdb-CMV overlay, the top hits, and a same-chain public-TCRβ
# co-occurrence panel. Run with:  marimo edit notebooks/biomarker_explorer.py
import marimo

__generated_with = "0.23.14"
app = marimo.App(width="medium")


@app.cell
def _(mo):
    mo.md(
        r"""
        # Biomarker association & co-occurrence — a vdjtools v2 explorer

        Emerson et al. (2017) screen public **TCRβ** chains for a statistical association
        between their *presence across a cohort* and a phenotype. `vdjtools.biomarker`
        generalises that one Fisher test along four axes — this notebook makes all four
        live on the **Emerson HIP** cohort (`isalgo/airr_hip`):

        - **condition** — a binary phenotype (**CMV**), an **HLA allele** (HLA-A\*02), or a
          paired condition combined by **Cochran–Mantel–Haenszel** (CMV *stratified by*
          HLA-A\*02 carriage — does the CMV signal survive controlling for HLA?);
        - **test** — Fisher exact, χ², a Bayesian **Beta-Binomial Bayes factor**, a Bayesian
          **log-odds** posterior, or a label **permutation** null;
        - **scope** — is a clonotype its CDR3 alone, **+V**, or **+V+J** (Emerson), matched
          **exactly** or within **1 mismatch** (metaclonotypes)?
        - **co-occurrence** — which public TCRβ *pairs* recur together across subjects
          (De Witt / Howie co-specificity — subjects are the "wells").

        Every hit is checked live against a local **VDJdb** dump.
        """
    )
    return


@app.cell
def _():
    # --- imports & configuration (single cell so every name is defined once) ---
    import os
    from pathlib import Path

    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import polars as pl

    from vdjtools import io as vio
    from vdjtools.io import schema as S
    from vdjtools.biomarker import association, cooccurrence, condition  # condition = the module

    REPO_ID = "isalgo/airr_hip"
    # Local VDJdb slim dump (antigenomics/vdjdb-db); point VDJDB_PATH at your checkout.
    # Validation degrades gracefully if absent — the screen itself needs no VDJdb.
    VDJDB = Path(os.environ.get("VDJDB_PATH",
                                "~/vcs/code/vdjdb-db/database/vdjdb.slim.txt")).expanduser()
    N_SUBJECTS = 400        # balanced subset; the full 786-subject run is appendix/bench_biomarker.py

    OKABE = {"blue": "#0072B2", "vermillion": "#D55E00", "green": "#009E73",
             "orange": "#E69F00", "purple": "#CC79A7", "grey": "#8C8C8C"}
    return (N_SUBJECTS, OKABE, Path, REPO_ID, S, VDJDB, association, condition,
            cooccurrence, mo, np, pl, plt, vio)


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 1 · Data bootstrap

        `metadata.txt` (786 subjects) carries **CMV serostatus** (`+`/`−`/`NA`) and **2-digit
        HLA-A/B typing**. We take a balanced subset and pull just those repertoires
        (`corr/HIP#####.txt.gz`, VDJtools format) from HuggingFace into the gitignored
        `notebooks/.data/biomarker_nb/` cache, stream them into a hive-partitioned Parquet
        cohort once, then scan it as one out-of-core `polars` LazyFrame.
        """
    )
    return


@app.cell
def _(N_SUBJECTS, Path, REPO_ID, mo, pl):
    _nb_dir = mo.notebook_dir() or Path.cwd()
    data_dir = _nb_dir / ".data" / "biomarker_nb"
    data_dir.mkdir(parents=True, exist_ok=True)
    try:
        import huggingface_hub as _hub
    except ImportError:
        _hub = None
    mo.stop(_hub is None, mo.md('> Install the examples extra: `pip install "vdjtools[examples]"`.'))

    _meta_path = _hub.hf_hub_download(REPO_ID, "metadata.txt", repo_type="dataset")
    meta_all = pl.read_csv(_meta_path, separator="\t", infer_schema_length=0)  # TAB — race has commas
    _known = meta_all.filter(pl.col("cmv").is_in(["+", "-"]))
    meta = pl.concat([_known.filter(pl.col("cmv") == c).head(N_SUBJECTS // 2) for c in ("+", "-")])
    # Derive HLA-A*02 carriage (2-digit typing) for the HLA + CMH conditions.
    meta = meta.with_columns(
        pl.when(pl.col("hla").is_null() | pl.col("hla").is_in(["", "NA"])).then(None)
          .otherwise(pl.col("hla").str.contains(r"HLA-A\*02")).alias("a02"))
    _ids = meta["sample_id"].to_list()
    _root = Path(_hub.snapshot_download(REPO_ID, repo_type="dataset",
                                        allow_patterns=[f"corr/{s}.txt.gz" for s in _ids]))
    corr_dir = _root / "corr"
    mo.md(f"**{len(_ids)} subjects** — {(meta['cmv']=='+').sum()} CMV+, {(meta['cmv']=='-').sum()} CMV−; "
          f"{meta['a02'].sum()} HLA-A\\*02 carriers. Cache: `{data_dir}`")
    return corr_dir, data_dir, meta


@app.cell
def _(corr_dir, data_dir, meta, mo, vio):
    cohort = data_dir / "cohort"
    if not any(cohort.glob("sample_id=*/*.parquet")):
        vio.ingest_cohort(meta.select("sample_id", "cmv", "hla", "a02"), corr_dir, cohort,
                          sample_col="sample_id", file_template="{sample}.txt.gz", fmt="vdjtools")
    lf = vio.scan_cohort(cohort, join_metadata=False)
    mo.md(f"Cohort ingested → `{cohort}`.")
    return (lf,)


@app.cell
def _(VDJDB, mo, pl):
    # VDJdb CMV validation set (human TRB), if the local dump is present.
    if VDJDB.exists():
        _v = pl.read_csv(VDJDB, separator="\t", infer_schema_length=0)
        vdjdb_cmv = (_v.filter((pl.col("gene") == "TRB") & (pl.col("species") == "HomoSapiens")
                               & pl.col("antigen.species").str.contains("CMV"))
                     .group_by("cdr3").agg(
                         pl.col("antigen.epitope").unique().sort().str.join(", ").alias("epitope"),
                         pl.col("mhc.a").unique().sort().str.join(", ").alias("mhc")))
        _msg = f"VDJdb CMV reference: **{vdjdb_cmv.height} human TRB CDR3s**."
    else:
        vdjdb_cmv = None
        _msg = f"> VDJdb not found at `{VDJDB}` — validation skipped."
    mo.md(_msg)
    return (vdjdb_cmv,)


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 2 · The screen — turn the knobs

        The screen re-runs when the **condition**, **test**, **scope**, or **min incidence**
        change; the **threshold** only re-draws. With CMH the *test* knob is ignored (the
        stratified odds are always combined the Mantel–Haenszel way).

        > **On significance.** With this many features BH-`q` is strict at a 400-subject
        > subset, so — like Emerson — the primary line is a **nominal p threshold** (paper:
        > P < 1e-4). Permutation is heavier than the closed-form tests; raise *min incidence*
        > to keep it snappy.
        """
    )
    return


@app.cell
def _(mo):
    cond = mo.ui.dropdown(
        {"CMV (enrichment)": "cmv", "HLA-A*02 (two-sided)": "hla",
         "CMV | HLA-A*02 (CMH)": "cmh"}, value="CMV (enrichment)", label="Condition")
    test = mo.ui.dropdown(
        {"Fisher exact": "fisher", "Chi²": "chi2", "Bayesian BF": "bayes_bf",
         "Bayesian log-odds": "bayes_logodds", "Permutation": "permutation"},
        value="Fisher exact", label="Test")
    match = mo.ui.dropdown({"exact CDR3": "exact", "1 mismatch (metaclonotype)": "1mm"},
                           value="exact CDR3", label="CDR3 matching")
    key = mo.ui.dropdown({"CDR3 + V + J (Emerson)": "vj", "CDR3 + V": "v", "CDR3 only": "cdr3"},
                         value="CDR3 + V + J (Emerson)", label="Clonotype key")
    min_inc = mo.ui.slider(3, 15, value=5, label="min incidence (subjects)")
    logp = mo.ui.slider(2.0, 8.0, value=4.0, step=0.5, label="−log10 p threshold")
    mo.hstack([cond, test, match, key, min_inc, logp], justify="start", gap=1.2)
    return cond, key, logp, match, min_inc, test


@app.cell
def _(S, association, cond, condition, key, lf, match, meta, min_inc, test):
    _KEYS = {"vj": (S.JUNCTION_AA, S.V_CALL, S.J_CALL),
             "v": (S.JUNCTION_AA, S.V_CALL), "cdr3": (S.JUNCTION_AA,)}
    if cond.value == "cmv":
        design, alt = condition.binary(meta, "cmv"), "greater"
    elif cond.value == "hla":
        design, alt = condition.binary(meta, "a02"), "two-sided"
    else:  # CMH — CMV stratified by HLA-A*02 carriage (test knob ignored, CMH is used)
        design, alt = condition.stratified(meta, "cmv", "a02"), "greater"
    res = association(lf, design, test=test.value, key=_KEYS[key.value], match=match.value,
                      min_incidence=min_inc.value, alternative=alt, n_perm=200)
    return (res,)


@app.cell
def _(cond, np, pl, res):
    # Choose the y-axis for whichever test produced `res` (p-based, Bayes-factor, or posterior).
    _cols = res.columns
    if "p_value" in _cols:
        y = np.clip(-np.log10(res["p_value"].to_numpy()), 0, 40)
        ylab, thr_kind = "−log10 p", "p"
    elif "log_bf" in _cols:
        y = np.clip(res["log_bf"].to_numpy(), -5, 40)
        ylab, thr_kind = "log Bayes factor", "bf"
    else:  # bayes_logodds → posterior P(OR>1)
        y = np.clip(-np.log10(1 - np.clip(res["p_or_gt1"].to_numpy(), 0, 1 - 1e-12)), 0, 40)
        ylab, thr_kind = "−log10 (1 − P(OR>1))", "p"
    x = res["log2_or"].to_numpy() if "log2_or" in _cols else np.zeros(len(y))
    juncs = res["junction_aa"] if "junction_aa" in _cols else pl.Series([""] * len(y))
    is_cmv = cond.value == "cmv"
    return is_cmv, juncs, thr_kind, x, y, ylab


@app.cell
def _(OKABE, cond, is_cmv, juncs, logp, np, plt, thr_kind, vdjdb_cmv, x, y, ylab):
    _cut = logp.value if thr_kind == "p" else 3.0        # log-BF>3 ≈ strong evidence
    _hit = y >= _cut
    _val = ((juncs.is_in(list(vdjdb_cmv["cdr3"])).to_numpy() if vdjdb_cmv is not None
             else np.zeros(len(x), bool)) & _hit)
    figv, axv = plt.subplots(figsize=(7.6, 5.2))
    axv.scatter(x[~_hit], y[~_hit], s=5, c="#c8ccd4", alpha=0.4, lw=0)
    axv.scatter(x[_hit], y[_hit], s=12, c=OKABE["vermillion"], alpha=0.75, lw=0,
                label=f"{ylab.split()[0]} ≥ {_cut:g}")
    if _val.any() and is_cmv:
        axv.scatter(x[_val], y[_val], s=44, facecolors="none", edgecolors=OKABE["blue"],
                    linewidths=1.4, label="in VDJdb-CMV")
    axv.axhline(_cut, ls="--", c="k", lw=0.6, alpha=0.5)
    axv.axvline(0, ls="-", c="k", lw=0.6, alpha=0.3)
    axv.set(xlabel="log2 odds ratio  (enriched →)", ylabel=ylab,
            title=f"{cond.selected_key} — associated TCRβ")
    axv.legend(fontsize=9, frameon=False, loc="upper left")
    axv.spines[["top", "right"]].set_visible(False)
    figv.tight_layout()
    figv
    return


@app.cell
def _(is_cmv, mo, pl, res, vdjdb_cmv):
    # Top hits — the significance column depends on the test that ran.
    _score = ("p_value" if "p_value" in res.columns
              else "log_bf" if "log_bf" in res.columns else "p_or_gt1")
    _desc = _score != "p_value"
    top = res.sort(_score, descending=_desc).head(200)
    _cols = ([c for c in ("junction_aa", "v_call", "j_call") if c in top.columns]
             + [c for c in ("incidence", "n_pos_present", "odds_ratio", "or_mh", "log_bf",
                            "p_value", "q_value") if c in top.columns])
    if vdjdb_cmv is not None and is_cmv and "junction_aa" in top.columns:
        table = (top.join(vdjdb_cmv, left_on="junction_aa", right_on="cdr3", how="left")
                 .select(_cols + ["epitope", "mhc"]))
        _n_val = top.join(vdjdb_cmv, left_on="junction_aa", right_on="cdr3", how="semi").height
        _msg = f"top {top.height} of {res.height:,} tested · **{_n_val} matched to VDJdb-CMV** (exact CDR3)"
    else:
        table = top.select(_cols)
        _msg = f"top {top.height} of {res.height:,} tested"
    mo.vstack([mo.md(_msg), table])
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 3 · Public-TCRβ co-occurrence (same-chain)

        Beyond single-clone association: which public TCRβ **pairs** recur *together* across
        subjects more than chance? `cooccurrence` builds a boolean subject×feature incidence
        matrix per chain, gets every pair's co-occurrence count by one matmul, and reports the
        **lift θ = n·n_AB/(n_A·n_B)** with a Fisher p and BH-`q` (De Witt / Howie
        co-specificity; on paired TRA+TRB data the same call does in-silico α-β pairing).
        """
    )
    return


@app.cell
def _(mo):
    cc_inc = mo.ui.slider(0.05, 0.5, value=0.15, step=0.05, label="min incidence fraction")
    cc_cooc = mo.ui.slider(2, 12, value=4, label="min co-occurrences")
    mo.hstack([cc_inc, cc_cooc], justify="start", gap=1.5)
    return cc_cooc, cc_inc


@app.cell
def _(S, cc_cooc, cc_inc, cooccurrence, lf, mo, pl):
    # Same-chain TCRβ pairs (hip is TRB-only). max_features caps the matmul for interactivity.
    cc = cooccurrence(lf, chain_a="TRB", chain_b="TRB", key=(S.JUNCTION_AA, S.V_CALL, S.J_CALL),
                      min_incidence_frac=cc_inc.value, min_cooccurrence=cc_cooc.value,
                      evalue=True, max_features=1500)
    _sig = cc.filter(pl.col("q_value") < 0.05) if cc.height else cc
    _show = (_sig.select("a_junction_aa", "b_junction_aa", "n_ab", "theta", "q_value", "e_value")
             .head(15) if cc.height else cc)
    mo.vstack([mo.md(f"**{cc.height:,}** candidate TCRβ pairs · **{_sig.height}** significant "
                     f"(q<0.05); top by p:"), _show])
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ---
        **Takeaway.** One streamed incidence table drives every knob. A pure incidence screen
        rediscovers CMV-specific public TCRβs and their HLA restriction; the CMH control shows
        which of the CMV signal is *not* just HLA-A\*02 tagging along; widening the key sharpens
        specificity and 1-mismatch recovers more VDJdb overlap; and the same substrate yields
        co-occurring public-clone pairs. Built on `vdjtools.io` (streamed Parquet cohort),
        `vdjtools.biomarker.{association, cooccurrence, condition, stats}`, and the native
        1-mismatch matcher. The full 786-subject benchmark is `appendix/bench_biomarker.py`.
        """
    )
    return


if __name__ == "__main__":
    app.run()
