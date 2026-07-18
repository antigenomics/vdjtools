# vdjtools — the ankylosing-spondylitis TRBV9 "AS27" motif: disease, or HLA-B27 carriage?
# Reactive marimo app over the Komech ankylosing-spondylitis TCRb cohort (`isalgo/airr_ankspond`,
# new/): reproduce the public TRBV9/TRBJ2-3 CDR3 motif (Komech et al. 2018, Rheumatology 57:1097)
# and settle its confound — B27 is ~26/27 tied to disease, so only the B27-MATCHED contrast
# (AS/B27+ vs healthy/B27+) can separate disease from carriage. Run with:
#     marimo edit notebooks/ankspond_motif.py
import marimo

__generated_with = "0.23.14"
app = marimo.App(width="medium")


@app.cell
def _(mo):
    mo.md(
        r"""
        # The ankylosing-spondylitis "AS27" motif — a vdjtools v2 explorer

        Komech et al. (2018, *Rheumatology* 57:1097) found a **public TRBV9 / TRBJ2-3 CDR3β
        motif** (`CASSVGLYSTDTQYF` and relatives) expanded in the blood and synovial fluid of
        **ankylosing-spondylitis** patients — nearly all of whom carry **HLA-B27**. That raises
        the question this notebook answers: is the motif a marker of **disease**, or merely of
        **B27 carriage**?

        Because B27 is ~26/27 confounded with AS in this cohort, "AS vs healthy" cannot tell them
        apart. The only contrast that can is **B27-matched**: AS/B27+ vs *healthy* B27+ donors —
        the paper's actual claim (the motif is absent from B27+ healthy blood). Two controls keep
        us honest: B27 carriage **among healthy donors** (should be null if it is disease), and
        the lone AS/B27− donor.

        Data auto-loads from HuggingFace (`isalgo/airr_ankspond`), preferring a local `~/hf/` copy.
        """
    )
    return


@app.cell
def _():
    # --- imports & configuration (single cell so every name is defined once) ---
    from pathlib import Path

    import marimo as mo
    import matplotlib.pyplot as plt
    import polars as pl
    from scipy.stats import fisher_exact

    from vdjtools.io.batch import read

    REPO = "isalgo/airr_ankspond"
    V_GENE, J_GENE = "TRBV9", "TRBJ2-3"
    # Komech 2018's motif IS this enumerated clonotype list (Fig.1 p.1100 + prose p.1101); the
    # paper prints no regex. All are TRBV9/TRBJ2-3.
    MOTIF = ["CASSVGLYSTDTQYF", "CASSVGLFSTDTQYF", "CASSVGVYSTDTQYF", "CASSVATYSTDTQYF",
             "CASSLGLFSTDTQYF", "CASSAGLFSTDTQYF", "CASSPGLFSTDTQYF", "CASSVGGFGDTQYF",
             "CASSAGLYSTDTQYF"]
    OKABE = {"blue": "#0072B2", "vermillion": "#D55E00", "green": "#009E73", "grey": "#8C8C8C"}

    def fetch(filename):
        """Local-first: ./<f>, ~/hf/airr_ankspond/<f>, else HuggingFace."""
        for root in (Path.cwd(), Path.home() / "hf" / "airr_ankspond"):
            if (root / filename).exists():
                return str(root / filename)
        import huggingface_hub as hub
        return hub.hf_hub_download(REPO, filename, repo_type="dataset")

    def hamming1_ball(seqs, universe):
        """Sequences in `universe` within Hamming distance 1 of any motif seq (same length)."""
        out = set(seqs) & set(universe)
        by_len = {}
        for u in universe:
            by_len.setdefault(len(u), []).append(u)
        for s in seqs:
            for u in by_len.get(len(s), ()):
                if sum(a != b for a, b in zip(s, u)) <= 1:
                    out.add(u)
        return list(out)

    return (MOTIF, OKABE, Path, REPO, V_GENE, J_GENE, fetch, hamming1_ball, mo, pl, plt,
            fisher_exact, read)


@app.cell
def _(fetch, pl, read):
    # Load the `new/` cohort: every sample tagged with its donor, plus per-donor metadata.
    meta = pl.read_csv(fetch("new/metadata.tsv"), separator="\t", infer_schema_length=0,
                       null_values=["nan"])
    _frames = []
    for _r in meta.iter_rows(named=True):
        _frames.append(read(fetch(f"new/{_r['sample_name']}.tsv.gz"))
                       .with_columns(pl.lit(_r["donor_id"]).alias("donor")))
    clones = pl.concat(_frames, how="vertical_relaxed")
    donors = meta.select("donor_id", "disease_status", "b27").unique("donor_id").rename(
        {"donor_id": "donor"})
    return clones, donors, meta


@app.cell
def _(donors, mo, pl):
    _tab = donors.group_by("disease_status", "b27").len().sort("disease_status", "b27")
    mo.vstack([
        mo.md("## 1 · The confound in one table\n\nB27 tracks disease, so only a **B27-matched** "
              "contrast is interpretable."),
        _tab,
    ])
    return


@app.cell
def _(mo):
    contrast = mo.ui.dropdown(
        {"AS/B27+ vs healthy/B27+ (B27-matched — the real test)": "matched",
         "AS/B27+ vs all healthy (confounded)": "all_hd",
         "B27+ vs B27− among HEALTHY only (carriage control)": "carriage"},
        value="AS/B27+ vs healthy/B27+ (B27-matched — the real test)", label="Contrast")
    onemm = mo.ui.switch(label="expand motif to a 1-mismatch ball")
    mo.hstack([contrast, onemm], justify="start", gap=1.5)
    return contrast, onemm


@app.cell
def _(MOTIF, V_GENE, J_GENE, clones, contrast, donors, fisher_exact, hamming1_ball,
      mo, onemm, pl):
    # V-pinned TRBV9/TRBJ2-3 CDR3 universe, optionally expanded to a 1-mismatch ball.
    _v = clones.filter((pl.col("v_call") == V_GENE) & (pl.col("j_call") == J_GENE))
    _seqs = (hamming1_ball(MOTIF, _v["junction_aa"].unique().to_list()) if onemm.value else MOTIF)
    carriers = set(_v.filter(pl.col("junction_aa").is_in(_seqs))["donor"].unique().to_list())

    def _arm(state, b27):
        d = donors.filter(pl.col("disease_status") == state)
        if b27 is not None:
            d = d.filter(pl.col("b27") == b27)
        return d["donor"].to_list()

    if contrast.value == "matched":
        pos, neg, plabel, nlabel = _arm("as", "pos"), _arm("hd", "pos"), "AS/B27+", "HD/B27+"
    elif contrast.value == "all_hd":
        pos, neg, plabel, nlabel = _arm("as", "pos"), _arm("hd", None), "AS/B27+", "all HD"
    else:                                                  # carriage among healthy
        pos, neg, plabel, nlabel = _arm("hd", "pos"), _arm("hd", "neg"), "HD/B27+", "HD/B27−"
    a, b = len(carriers & set(pos)), len(carriers & set(neg))
    orr, p = fisher_exact([[a, len(pos) - a], [b, len(neg) - b]], alternative="greater")
    _verdict = ("**disease-associated**" if contrast.value == "matched" and p < 0.05
                else "**null (not carriage)**" if contrast.value == "carriage" and p > 0.05
                else "")
    mo.md(f"## 2 · Motif carriage — {'1-mismatch ball' if onemm.value else 'exact'}\n\n"
          f"**{plabel}: {a}/{len(pos)}** carry the motif vs **{nlabel}: {b}/{len(neg)}** — "
          f"OR = **{orr:.1f}**, one-sided p = **{p:.4f}** {_verdict}")
    return carriers, neg, pos


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 3 · The motif is a convergent family (metaclonotype clustering)

        The nine motif clonotypes are near-variants of one another — a **convergent** response, not
        one public clone. `vdjtools.biomarker.metaclonotypes` clusters the AS patients' TRBV9/TRBJ2-3
        CDR3s within one substitution (single-linkage); the largest component is the AS27 family.
        """
    )
    return


@app.cell
def _(MOTIF, V_GENE, J_GENE, clones, donors, mo, pl):
    from vdjtools.biomarker.metaclonotype import metaclonotypes

    _as = donors.filter(pl.col("disease_status") == "as")["donor"].to_list()
    _v = (clones.filter((pl.col("v_call") == V_GENE) & (pl.col("j_call") == J_GENE)
                        & pl.col("donor").is_in(_as))
          .select("junction_aa", "v_call", "j_call").unique())
    meta_id = metaclonotypes(_v, scope="1,0,0,1", match_v=True, match_j=True)
    _sizes = meta_id.group_by("meta_id").len().sort("len", descending=True)
    _big = _sizes["meta_id"][0]
    fam = (meta_id.filter(pl.col("meta_id") == _big)
           .with_columns(pl.col("junction_aa").is_in(MOTIF).alias("in_komech"))
           .sort("in_komech", descending=True))
    mo.vstack([
        mo.md(f"Largest TRBV9/TRBJ2-3 metaclonotype: **{fam.height} CDR3 variants**, of which "
              f"**{int(fam['in_komech'].sum())}** are Komech-2018 clonotypes. A sample:"),
        fam.select("junction_aa", "in_komech").head(14),
    ])
    return fam, meta_id


@app.cell
def _(mo):
    mo.md(
        r"""
        ---
        **Takeaway.** The TRBV9/TRBJ2-3 AS27 motif is real and **disease-associated**: B27-matched,
        AS patients carry it far more than B27+ healthy donors (the confounded "AS vs all healthy"
        inflates the effect, and B27 carriage among the healthy is null — so it is disease, not
        carriage). It is a **convergent family** of near-variant CDR3s, recovered here by
        `biomarker.metaclonotypes`. V-pinning is load-bearing — unpinned, wrong-V convergents leak
        into the healthy arm. The full campaign (Pgen of the motif under the native model, VDJdb
        cross-check, de-novo V+k-mer discovery) is `bench/bm_ankspond.py`.
        """
    )
    return


if __name__ == "__main__":
    app.run()
