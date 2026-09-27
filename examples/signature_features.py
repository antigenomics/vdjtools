# vdjtools — the repertoire signature: transforms, V-call resolution, and the corpus rotation.
#
# Reactive marimo app. Everything here runs on repertoires SAMPLED FROM THE BUNDLED MODELS, so
# there is no download and no cohort: the point is the feature machinery, and a generated
# repertoire exercises it exactly as a real one does.
#
# Three things worth knowing before you build features from a repertoire, each measurable in a
# few seconds here:
#   1. why the amino-acid block is arcsine-transformed and not log1p of counts;
#   2. why an ambiguous V call must be resolved, not stripped;
#   3. k-mers as a raw feature group, and why the component count is a real choice.
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
        ## 3. k-mers, and why the component count is a real choice

        A junction k-mer profile is now just another **raw feature group** — 400 columns per locus
        at `k=2` — concatenated with usage, spectratype, composition and the count statistics, and
        rotated by the corpus's own per-locus PCA. There is no separate frozen k-mer space and no
        separate TF-IDF basis any more: one corpus, one rotation per locus, one place where
        anything is fitted.

        That makes `--components` the only knob controlling how much of the k-mer signal survives,
        and it is not a free choice.
        """
    )
    return


@app.cell
def _(mo):
    mo.md(
        """
        ### Do not pick components by explained variance alone

        A truncated SVD keeps the directions of greatest variance *in the fitting corpus*, which for
        repertoires are depth, V usage and batch. A motif carried by a handful of clonotypes in a
        handful of donors is not one of those.

        Measured on ankylosing spondylitis vs healthy, both HLA-B27+, with the basis fitted on a
        disjoint cohort:

        | read-out | AUC | perm. *p* |
        |---|---|---|
        | sum of the 17 columns the published motif occupies | 0.769 | **0.031** |
        | best of 64 SVD components | 0.841 | 0.20 |
        | best single vocabulary column | 0.813 | — |

        The two larger numbers are the meaningless ones: a maximum over 64 components reaches
        \\|AUC − 0.5\\| = 0.34 under a label permutation null. Any best-of-N read-out must be nulled
        or not quoted.

        So: **rare and discriminative** → do not rotate. Take the raw features
        (`features.raw_and_channels`) and an L1 model. **Broad compositional shift** → rotate, and
        keep only components that survive a study-disjoint refit, rather than components that reach
        a cumulative-variance threshold.

        This is why the corpus stores its full eigenvalue spectrum and why `--components` accepts a
        count as well as a fraction: the fraction is convenient, the count is what you pin once you
        have measured which components transfer.
        """
    )
    return


@app.cell
def _(mo):
    mo.md(
        """
        ## The contract

        Column names are `vsig:<block>:<locus>:<feature>` throughout, so one `parse()` reads all
        three kinds: a **raw feature** (an input to the rotation), a **rotated column**
        (`vsig:pc:<locus>:PCnn`), and a **channel** (carried through untouched).

        Two fields are declared per raw feature and neither is inferred: its `transform`, applied
        where the feature is computed while its denominator is still in scope, and its `support`,
        which decides the winsorization side. The support is separate because `transform="none"`
        spans three different supports — a log-probability is bounded above, a standard deviation
        below, a log-ratio neither — so deriving the trimming side from the transform silently
        trims the wrong end of two of them.
        """
    )
    return


@app.cell
def _(pl):
    from vdjtools.signature import layout

    # The registry holds BOTH namespaces -- `rsig:` groups are declared here too, because the
    # shared contract lives in vdjtools and mir.signature registers into it. Nothing here imports
    # mir. Widths per locus are the ROTATION INPUT: what the PCA sees, not what comes out.
    pl.DataFrame([{"group": g.name, "sig": g.sig, "loci": len(g.emitted_loci),
                   "features": len(g.features) or "germline-wide",
                   "support": g.dynamic[1] if g.dynamic else
                              ", ".join(sorted({v[1] for v in g.features.values()}))}
                  for g in layout.raw_groups()])
    return (layout,)


@app.cell
def _(mo):
    mo.md(
        """
        ## Rotating a sample: a corpus is required

        There is no default corpus. A signature is only comparable to another one rotated through
        the *same* corpus, so the artifact has to be named — and the emitted row carries which one
        it was, how much of itself got clamped by that corpus's bounds, and the coverage the sample
        actually attained.

        Build one with `vdjtools corpus --smoke -o /tmp/c.npz` (minutes), or load a bundled name.
        """
    )
    return


@app.cell
def _(deep):
    from vdjtools.signature import vsig
    from vdjtools.signature.corpus import synthesize

    # A deliberately tiny corpus so this notebook runs in seconds. The real ones are 10,000
    # repertoires; see `vdjtools corpus --help`.
    corpus, _rows = synthesize("memory", loci=("TRB",), n_samples=30, size=120, seed=1,
                               n_components=6)
    v = vsig({"TRB": deep}, corpus)
    {k: v[k] for k in list(v)[:10]}
    return corpus, v, vsig


@app.cell
def _(v):
    # The three channels a caller should read before trusting any rotated column of this row.
    {k: v[k] for k in ("vsig:cov:TRB:cstar", "vsig:mask:TRB:estimable",
                       "vsig:qc:-:winsor_frac")}
    return


if __name__ == "__main__":
    app.run()
