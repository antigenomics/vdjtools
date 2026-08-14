# vdjtools — the repertoire signature: transforms, V-call resolution, and the V+k-mer space.
#
# Reactive marimo app. Everything here runs on repertoires SAMPLED FROM THE BUNDLED MODELS, so
# there is no download and no cohort: the point is the feature machinery, and a generated
# repertoire exercises it exactly as a real one does.
#
# Three things worth knowing before you build features from a repertoire, each measurable in a
# few seconds here:
#   1. why the amino-acid block is arcsine-transformed and not log1p of counts;
#   2. why an ambiguous V call must be resolved, not stripped;
#   3. what the V+k-mer space is, and why you should not select its components by variance.
#
# Run with:  marimo edit examples/signature_features.py
import marimo

__generated_with = "0.23.14"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _():
    import numpy as np
    import polars as pl

    from vdjtools.model import load_bundled
    from vdjtools.model.generate import generate
    return generate, load_bundled, np, pl


@app.cell
def _(mo):
    mo.md(
        """
        # The repertoire signature — feature machinery

        `vdjtools.signature` turns one AIRR sample into a fixed, named, positional vector.
        This notebook is about the three choices inside it that are easy to get wrong.
        """
    )
    return


@app.cell
def _(mo):
    depth = mo.ui.slider(200, 50_000, value=5_000, label="deep sample size (clonotypes)")
    depth
    return (depth,)


@app.cell
def _(depth, generate, load_bundled, pl):
    # Two samples of the SAME underlying process at very different depths. Any feature that
    # separates them is measuring the sequencer, not the donor.
    model = load_bundled("TRB", source="olga")

    def clones(n, seed):
        """`generate` emits one row per rearrangement; a clonotype frame counts the duplicates."""
        return (generate(model, n, seed=seed)
                .group_by(["junction_aa", "v_call", "j_call"]).len()
                .rename({"len": "duplicate_count"})
                .with_columns((pl.col("duplicate_count")
                               / pl.col("duplicate_count").sum()).alias("frequency")))

    shallow, deep = clones(200, 1), clones(depth.value, 2)
    return clones, deep, model, shallow


@app.cell
def _(mo):
    mo.md(
        """
        ## 1. Arcsine, not log1p of counts

        The amino-acid composition block is a set of proportions with a known denominator — the
        number of junction residues actually observed. Anscombe's arcsine,

        $$\\arcsin\\sqrt{\\frac{xm + 3/8}{m + 3/4}},$$

        is depth-invariant by construction; `log1p` of the raw count is not.
        """
    )
    return


@app.cell
def _(deep, np, pl, shallow):
    from vdjtools.signature.transform import arcsine

    def aa_share(df, aa="G"):
        s = df["junction_aa"]
        m = float(s.str.len_chars().sum())
        x = float(s.str.count_matches(aa).sum()) / m
        return x, m

    rows = []
    for _name, _df in (("shallow", shallow), ("deep", deep)):
        _x, _m = aa_share(_df)
        rows.append({"sample": _name, "residues m": int(_m), "proportion p": _x,
                     "arcsine(p, m)": float(arcsine(np.array([_x]), _m)[0]),
                     "log1p(count)": float(np.log1p(_x * _m))})
    drift = pl.DataFrame(rows)
    drift
    return arcsine, drift


@app.cell
def _(drift, mo):
    _a = drift["arcsine(p, m)"].to_list()
    _l = drift["log1p(count)"].to_list()
    mo.md(
        f"""
        Same biology, two depths. `arcsine` moves **{max(_a) / min(_a):.2f}×**;
        `log1p` of the count moves **{max(_l) / min(_l):.2f}×**.

        `log1p` is not measuring composition, it is measuring how much you sequenced. And note
        `arcsine(0, m)` still depends on `m`: a residue never seen in 200 clonotypes and one never
        seen in 50,000 are different evidence, and the transform says so, where `log1p(0) = 0`
        always. A 5% winsorization would not fix this — it is a parameter fitted to some corpus's
        depth distribution, and it clips the dominant residues, which is where the biology is.
        """
    )
    return


@app.cell
def _(mo):
    mo.md(
        """
        ## 2. Resolve an ambiguous V call, do not strip it

        Amplicon data realigned from junction plus short flanks routinely calls a comma-separated
        tie. Left alone it shatters every V-keyed feature space into singleton columns.
        """
    )
    return


@app.cell
def _(pl):
    from vdjtools.io.schema import resolve_gene, strip_allele

    ambiguous = pl.DataFrame({"v_call": ["TRBV5-1*01,TRBV5-5*01", "TRBV19*01", "TRBV6-2*01,TRBV6-3*01"]})
    ambiguous.with_columns(
        resolve_gene(pl.col("v_call")).alias("resolve_gene  (use this)"),
        strip_allele(pl.col("v_call")).alias("strip_allele  (reporting)"),
    )
    return resolve_gene, strip_allele


@app.cell
def _(mo):
    mo.md(
        """
        `resolve_gene` takes the **first** gene — the aligner orders by score, so the first is the
        best-supported call. `strip_allele` **sorts** the tie so that reporting is order-insensitive,
        which is why composing the two would hand you the alphabetically first gene rather than the
        aligner's. For anything keyed on V, use `resolve_gene`.
        """
    )
    return


@app.cell
def _(mo):
    mo.md(
        """
        ## 3. The V + k-mer space

        A junction k-mer profile keyed jointly on the V gene, TF-IDF scaled, projected onto a
        frozen truncated-SVD basis. Counting is C++; the pattern string marks kept positions
        (`"xx.x"` is gapped), and `n_groups < 20` clusters residues by BLOSUM62 via classical MDS
        plus Ward linkage — no hand-picked chemistry classes.
        """
    )
    return


@app.cell
def _(mo):
    alpha = mo.ui.dropdown({"20 (plain amino acids)": 20, "8 (BLOSUM62 groups)": 8,
                            "12 (BLOSUM62 groups)": 12}, value="20 (plain amino acids)",
                           label="alphabet")
    pattern = mo.ui.dropdown(["xxxx", "xxx", "xx.x", "x.xx"], value="xxxx", label="k-mer pattern")
    mo.hstack([alpha, pattern])
    return alpha, pattern


@app.cell
def _(alpha, clones, pattern, pl):
    from vdjtools.features.kmer_space import fit_kmer_space

    # A corpus of small repertoires: the vocabulary and the IDF are set by how many SAMPLES back
    # each document frequency, not by how deep any one of them is.
    corpus = [clones(500, 100 + i) for i in range(40)]
    space = fit_kmer_space(corpus, pattern=pattern.value, n_groups=alpha.value, flank=4,
                           min_df=0.02, max_df=0.99, n_components=8, max_columns=100_000)
    pl.DataFrame([{"pattern": pattern.value,
                   "alphabet": alpha.value,
                   "code space": space.meta["code_space"],
                   "surviving columns": space.meta["surviving"],
                   "kept": space.n_columns,
                   "components": space.n_components}])
    return corpus, fit_kmer_space, space


@app.cell
def _(corpus, space):
    space.transform(corpus[0], weight="freq", residual=True)
    return


@app.cell
def _(mo):
    mo.md(
        """
        The last value is the **residual** — the norm of what the retained components threw away.
        It rides with the projection because a coordinate vector without it cannot tell you whether
        the sample was well described by the basis at all.

        ### Do not pick components by explained variance

        A truncated SVD keeps the directions of greatest variance *in the fitting corpus*, which
        for repertoires are depth, V usage and batch. A motif carried by a handful of clonotypes in
        a handful of donors is not one of those.

        Measured on ankylosing spondylitis vs healthy, both HLA-B27+, with the space fitted on a
        disjoint cohort:

        | read-out | AUC | perm. *p* |
        |---|---|---|
        | sum of the 17 columns the published motif occupies | 0.769 | **0.031** |
        | best of 64 SVD components | 0.841 | 0.20 |
        | best single vocabulary column | 0.813 | — |

        The two larger numbers are the meaningless ones: a maximum over 64 components reaches
        \\|AUC − 0.5\\| = 0.34 under a label permutation null. Any best-of-N read-out must be
        nulled or not quoted.

        So: **rare and discriminative** → keep the un-projected sparse columns and an L1 model;
        **broad compositional shift** → project, and keep only components that survive a
        study-disjoint refit, not components that reach a cumulative-variance threshold.
        """
    )
    return


@app.cell
def _(mo):
    mo.md(
        """
        ## The contract

        Column names are `vsig:<block>:<locus>:<feature>`, and the tiers are exact index subsets of
        one frozen layout — a narrower tier is a slice of a wider one, never a differently-computed
        number. `layout.columns()` reads no data at all.
        """
    )
    return


@app.cell
def _(pl):
    from vdjtools.signature import layout

    # The layout registry holds BOTH namespaces -- `rsig:` blocks are declared here too, because
    # the shared contract machinery lives in vdjtools and mir.signature registers into it. `vsig`
    # returns the `vsig:` half; the geometry half comes from `mir.signature`.
    pl.DataFrame([{"tier": t,
                   "vsig columns": sum(c.startswith("vsig:") for c in layout.columns(t)),
                   "whole contract": len(layout.columns(t))}
                  for t in ("core", "standard", "full")])
    return (layout,)


@app.cell
def _(deep):
    from vdjtools.signature import vsig

    v = vsig({"TRB": deep}, tier="core")
    {k: v[k] for k in list(v)[:12]}
    return v, vsig


if __name__ == "__main__":
    app.run()
