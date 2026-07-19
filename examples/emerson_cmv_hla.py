# vdjtools — Emerson 2017 biomarker benchmark at cohort scale: CMV/HLA-associated public TCRβ.
# Reactive marimo app over the isalgo/airr_hip cohort (up to 786 subjects, VDJtools format, with
# per-subject CMV serostatus + HLA typing). The cohort streams into a hive-partitioned Parquet
# dataset one sample at a time (ingest_cohort), is analysed as one out-of-core LazyFrame
# (scan_cohort), and millions of per-feature Fisher tests are vectorised through the hypergeometric
# tail — no per-feature Python loop, the cohort never fully in RAM. Hits are validated against
# VDJdb (fetched from the antigenomics/vdjdb-db release). Run with:
#     marimo edit examples/emerson_cmv_hla.py
import marimo

__generated_with = "0.23.14"
app = marimo.App(width="medium")


@app.cell
def _(mo):
    mo.md(
        r"""
        # Emerson CMV / HLA biomarkers at cohort scale — a vdjtools explorer

        Emerson et al. (*Nat Genet* 2017) screen public **TCRβ** chains for association with a
        phenotype across a large cohort. This notebook reproduces the core screen on the
        **`isalgo/airr_hip`** cohort (the Emerson HIP repertoires, with per-subject CMV serostatus
        and 2-digit HLA-A/B typing) with `vdjtools.biomarker.fisher_association`:

        - **CMV** — a *one-tailed* enrichment test among CMV+ subjects (Emerson's setting);
        - **HLA-A\*02** — a *two-tailed* test (positive + negative association);
        - hits are validated live against **VDJdb** (CMV epitope + HLA allele), fetched from the
          `antigenomics/vdjdb-db` release.

        **Scale is the point:** the cohort is streamed into a hive-partitioned Parquet dataset one
        sample at a time, analysed as a single out-of-core `polars` LazyFrame, and the per-feature
        Fisher tests are vectorised through the hypergeometric tail — the cohort never fully in RAM.
        Slide **subjects** up toward the full cohort; data prefers a local `./data_dump/airr_hip/`
        copy (gitignored), else fetches from HuggingFace.
        """
    )
    return


@app.cell
def _():
    # --- imports, config & helpers (single cell so every name is defined once) ---
    import resource
    import sys
    import time
    from pathlib import Path

    import marimo as mo
    import numpy as np
    import polars as pl

    from vdjtools import io as vio
    from vdjtools.io import schema as S
    from vdjtools.biomarker import fisher_association

    REPO = "isalgo/airr_hip"
    OKABE = {"ns": "#c8ccd4", "sig": "#d1495b", "val": "#00798c"}

    def peak_mb():
        """Peak resident set (MB); ru_maxrss is bytes on macOS, KiB on Linux. Process-cumulative."""
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return rss / 1e6 if sys.platform == "darwin" else rss / 1e3

    def hip_local(nb):
        """./data_dump/airr_hip/ (with metadata.txt + corr/) if present, else None."""
        for root in (Path.cwd(), nb, nb.parent):
            cand = root / "data_dump" / "airr_hip"
            if (cand / "metadata.txt").exists():
                return cand
        return None

    def phenotypes(meta):
        """Per-subject CMV (+/−) and HLA-A*02 (present/absent) phenotype tables."""
        cmv = meta.select("sample_id",
            pl.when(pl.col("cmv") == "+").then(True)
              .when(pl.col("cmv") == "-").then(False).otherwise(None).alias("cmv_pos"))
        a02 = meta.select("sample_id",
            pl.when(pl.col("hla").is_null() | pl.col("hla").is_in(["", "NA"])).then(None)
              .otherwise(pl.col("hla").str.contains(r"HLA-A\*02")).alias("a02"))
        return cmv, a02

    def vdjdb_slim():
        """Path to vdjdb.slim.txt(.gz): a local ./data_dump/ copy if present, else fetch the latest
        antigenomics/vdjdb-db release into ./data_dump/ (cached). None if unavailable."""
        nb = mo.notebook_dir() or Path.cwd()
        for root in (Path.cwd(), nb, nb.parent):
            for name in ("vdjdb.slim.txt", "vdjdb.slim.txt.gz"):
                if (root / "data_dump" / name).exists():
                    return root / "data_dump" / name
        import io as _io, json, urllib.request, zipfile
        dd = (nb.parent if nb.name == "examples" else Path.cwd()) / "data_dump"
        try:
            dd.mkdir(parents=True, exist_ok=True)
            _rel = json.load(urllib.request.urlopen(
                "https://api.github.com/repos/antigenomics/vdjdb-db/releases/latest", timeout=30))
            _url = next(a["browser_download_url"] for a in _rel["assets"] if a["name"].endswith(".zip"))
            with urllib.request.urlopen(_url, timeout=180) as _r:
                _z = zipfile.ZipFile(_io.BytesIO(_r.read()))
            _m = next(n for n in _z.namelist() if n.endswith(("vdjdb.slim.txt", "vdjdb.slim.txt.gz")))
            _dest = dd / Path(_m).name
            _dest.write_bytes(_z.read(_m))
            return _dest
        except Exception:
            return None

    def load_vdjdb_cmv(path):
        """Human TRB CMV-specific records from the VDJdb slim dump (cdr3, V, epitope, MHC)."""
        v = pl.read_csv(path, separator="\t", infer_schema_length=0)
        return (v.filter((pl.col("gene") == "TRB") & (pl.col("species") == "HomoSapiens")
                         & pl.col("antigen.species").str.contains("CMV"))
                .select(pl.col("cdr3"), S.strip_allele(pl.col("v.segm")).alias("vdjdb_v"),
                        pl.col("antigen.epitope").alias("epitope"), pl.col("mhc.a").alias("mhc"))
                .unique())

    def volcano(hits, vdjdb_cdr3, title, fdr):
        """Inline volcano: log2 odds-ratio (x) vs −log10 p (y); significant + vdjdb-validated marked."""
        import matplotlib.pyplot as plt
        h = hits.with_columns((-pl.col("p_value").log10()).alias("_y"),
                              (pl.col("q_value") < fdr).alias("_sig"))
        x = h["log2_or"].to_numpy()
        y = np.clip(h["_y"].to_numpy(), 0, 320)
        sig = h["_sig"].to_numpy()
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.scatter(x[~sig], y[~sig], s=4, c=OKABE["ns"], alpha=0.4, linewidths=0, label="ns")
        ax.scatter(x[sig], y[sig], s=8, c=OKABE["sig"], alpha=0.7, linewidths=0, label=f"q<{fdr}")
        if vdjdb_cdr3:
            in_db = h[S.JUNCTION_AA].is_in(list(vdjdb_cdr3)).to_numpy()
            val = sig & in_db
            ax.scatter(x[val], y[val], s=30, facecolors="none", edgecolors=OKABE["val"],
                       linewidths=1.3, label="sig & in vdjdb-CMV")
        ax.axhline(-np.log10(0.05), ls="--", c="k", lw=0.6, alpha=0.5)
        ax.axvline(0, ls="-", c="k", lw=0.6, alpha=0.3)
        ax.set(xlabel="log2 odds ratio  (enriched →)", ylabel="−log10 p-value", title=title)
        ax.legend(fontsize=8, frameon=False, loc="upper left")
        ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()
        return fig

    return (OKABE, Path, REPO, S, fisher_association, hip_local, load_vdjdb_cmv, mo,
            np, peak_mb, phenotypes, pl, vdjdb_slim, vio, volcano)


@app.cell
def _(mo):
    max_samples = mo.ui.slider(40, 500, value=120, step=20,
                               label="subjects (balanced CMV+/−; overlap/Fisher scale with this)")
    min_incidence = mo.ui.slider(2, 8, value=2, step=1, label="min incidence (public floor)")
    fdr = mo.ui.dropdown(["0.01", "0.05", "0.1"], value="0.05", label="FDR")
    with_1mm = mo.ui.switch(label="also run the 1-mismatch metaclonotype CMV screen")
    mo.vstack([mo.hstack([max_samples, min_incidence], justify="start", gap=1.5),
               mo.hstack([fdr, with_1mm], justify="start", gap=1.5)])
    return fdr, max_samples, min_incidence, with_1mm


@app.cell
def _(Path, REPO, hip_local, max_samples, mo, phenotypes, pl, vio):
    # Resolve the cohort (local ./data_dump/airr_hip or HuggingFace), pick a balanced subset,
    # stream it into a hive-partitioned Parquet cohort once, scan it as one LazyFrame.
    _nb = mo.notebook_dir() or Path.cwd()
    _base = hip_local(_nb)
    if _base is not None:
        meta_all = pl.read_csv(_base / "metadata.txt", separator="\t", infer_schema_length=0)
        corr_dir = _base / "corr"
    else:
        import huggingface_hub as _hub
        meta_all = pl.read_csv(_hub.hf_hub_download(REPO, "metadata.txt", repo_type="dataset"),
                               separator="\t", infer_schema_length=0)   # TAB — race field has commas
        corr_dir = None
    _known = meta_all.filter(pl.col("cmv").is_in(["+", "-"]))
    meta = pl.concat([_known.filter(pl.col("cmv") == c).head(max_samples.value // 2)
                      for c in ("+", "-")])
    if corr_dir is None:
        import huggingface_hub as _hub
        _root = Path(_hub.snapshot_download(REPO, repo_type="dataset",
                     allow_patterns=[f"corr/{s}.txt.gz" for s in meta["sample_id"]]))
        corr_dir = _root / "corr"

    cmv_ph, a02_ph = phenotypes(meta)
    cohort = (mo.notebook_dir() or Path.cwd()) / ".data" / "emerson_nb" / f"cohort_{meta.height}"
    if not any(cohort.glob("sample_id=*/*.parquet")):
        vio.ingest_cohort(meta.select("sample_id", "cmv", "hla", "age", "sex"), corr_dir, cohort,
                          sample_col="sample_id", file_template="{sample}.txt.gz", fmt="vdjtools")
    lf = vio.scan_cohort(cohort, join_metadata=False)
    mo.md(f"**{meta.height} subjects** — {(meta['cmv']=='+').sum()} CMV+, "
          f"{(meta['cmv']=='-').sum()} CMV−; {meta['hla'].str.contains(r'HLA-A\*02').sum()} "
          f"HLA-A\\*02 carriers. Cohort scanned from `{cohort.name}`.")
    return a02_ph, cmv_ph, lf, meta


@app.cell
def _(cmv_ph, fdr, fisher_association, lf, min_incidence, mo, peak_mb, pl, time):
    # CMV: one-tailed enrichment among CMV+ (Emerson's setting), full V+CDR3+J key.
    _t0 = time.perf_counter()
    cmv = fisher_association(lf, cmv_ph, pheno_col="cmv_pos", alternative="greater",
                             min_incidence=min_incidence.value)
    _fdr = float(fdr.value)
    n_sig = cmv.filter(pl.col("q_value") < _fdr).height
    mo.md(f"## CMV — one-tailed enrichment\n\n**{cmv.height} public TCRβ tested, {n_sig} "
          f"significant (q<{_fdr})** in {time.perf_counter()-_t0:.1f}s, peak RSS {peak_mb():.0f} MB.")
    cmv.sort("q_value").head(8).select("junction_aa", "v_call", "j_call", "n_pos_present",
                                       "n_neg_present", "odds_ratio", "p_value", "q_value")
    return (cmv,)


@app.cell
def _(a02_ph, fdr, fisher_association, lf, min_incidence, mo, pl, time):
    # HLA-A*02: two-tailed (positive + negative association).
    _t0 = time.perf_counter()
    hla = fisher_association(lf, a02_ph, pheno_col="a02", alternative="two-sided",
                             min_incidence=min_incidence.value)
    _n = hla.filter(pl.col("q_value") < float(fdr.value)).height
    mo.md(f"## HLA-A\\*02 — two-tailed\n\n**{hla.height} tested, {_n} significant "
          f"(q<{fdr.value})** in {time.perf_counter()-_t0:.1f}s.")
    return (hla,)


@app.cell
def _(load_vdjdb_cmv, mo, vdjdb_slim):
    # VDJdb CMV reference (human TRB): ./data_dump/ copy, else fetched from the vdjdb-db release.
    _path = vdjdb_slim()
    vdjdb = load_vdjdb_cmv(_path) if _path is not None else None
    mo.md(f"VDJdb CMV reference (`{_path.name}`): **{vdjdb.height} human TRB CDR3s**."
          if vdjdb is not None else
          "> VDJdb unavailable (no `./data_dump/` copy and the release fetch failed) — validation skipped.")
    return (vdjdb,)


@app.cell
def _(cmv, fdr, mo, pl, vdjdb):
    # Exact-CDR3 overlap of significant CMV-enriched hits with VDJdb CMV entries.
    if vdjdb is None:
        _out = mo.md("_VDJdb validation skipped._")
    else:
        _sig = cmv.filter((pl.col("q_value") < float(fdr.value)) & (pl.col("direction") == "enriched"))
        _ex = _sig.join(vdjdb, left_on="junction_aa", right_on="cdr3", how="inner")
        if _ex.height == 0:
            _out = mo.md(f"**{_sig.height} significant CMV-enriched TCRβ** — none exactly in "
                         f"VDJdb-CMV at this scale (slide subjects up).")
        else:
            _tab = (_ex.sort("q_value").head(12).select(
                "junction_aa", "v_call", "odds_ratio", "q_value", "epitope", "mhc"))
            _out = mo.vstack([
                mo.md(f"**{_ex['junction_aa'].n_unique()} / {_sig.height}** significant hits are "
                      f"exact VDJdb-CMV clonotypes — the top by q:"), _tab])
    _out
    return


@app.cell
def _(cmv, fdr, mo, vdjdb, volcano):
    _db = set(vdjdb["cdr3"].to_list()) if vdjdb is not None else None
    mo.md("## Volcano — CMV")
    volcano(cmv, _db, "CMV-associated TCRβ (Emerson HIP)", float(fdr.value))
    return


@app.cell
def _(fdr, hla, mo, volcano):
    mo.md("## Volcano — HLA-A*02")
    volcano(hla, None, "HLA-A*02-associated TCRβ (Emerson HIP)", float(fdr.value))
    return


@app.cell
def _(S, cmv_ph, fdr, fisher_association, lf, min_incidence, mo, pl, time, with_1mm):
    # Optional 1-mismatch metaclonotype CMV screen over public keys (incidence≥min).
    if not with_1mm.value:
        cmv1 = None
        _out = mo.md("_Toggle the 1-mismatch switch above to also run the metaclonotype screen._")
    else:
        _t0 = time.perf_counter()
        _pub = (lf.group_by([S.JUNCTION_AA, S.V_CALL, S.J_CALL]).agg(pl.len().alias("_n"))
                .filter(pl.col("_n") >= min_incidence.value)
                .select([S.JUNCTION_AA, S.V_CALL, S.J_CALL]))
        cmv1 = fisher_association(
            lf.join(_pub, on=[S.JUNCTION_AA, S.V_CALL, S.J_CALL], how="semi"),
            cmv_ph, pheno_col="cmv_pos", alternative="greater",
            min_incidence=min_incidence.value, match="1mm")
        _n = cmv1.filter(pl.col("q_value") < float(fdr.value)).height
        _out = mo.md(f"**1-mismatch metaclonotype screen:** {cmv1.height} metaclonotypes, "
                     f"{_n} significant (q<{fdr.value}) in {time.perf_counter()-_t0:.1f}s.")
    _out
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ---
        **Takeaway.** The incidence-based Fisher screen recovers CMV-associated public TCRβ that
        overlap VDJdb's CMV-specific clonotypes, and HLA-A\*02-associated chains — at cohort scale,
        from raw repertoires, with the cohort never held in RAM. Slide **subjects** toward the full
        cohort for the published signal; the interactive, condition×test×scope version of this
        screen is `examples/emerson_biomarker.py`.
        """
    )
    return


if __name__ == "__main__":
    app.run()
