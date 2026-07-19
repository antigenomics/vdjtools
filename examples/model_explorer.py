# vdjtools — interactive recombination-model explorer.
# Reactive marimo app: pick a locus and model source, see its Bayes net (nodes = marginal entropy H,
# edges = mutual information I), the per-event entropy / MI tables, a marginal distribution, and an
# OLGA-vs-learned comparison. Run with:  marimo edit examples/model_explorer.py
import marimo

__generated_with = "0.23.14"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _(mo):
    import subprocess

    import matplotlib.pyplot as plt
    import polars as pl

    from vdjtools.model import analyze, list_bundled, load_bundled
    return analyze, list_bundled, load_bundled, mo, pl, plt, subprocess


@app.cell
def _(mo):
    mo.md(
        """
        # V(D)J recombination model explorer

        Explore any bundled recombination model as a **Bayes net** — nodes are recombination events
        (sized/annotated by their marginal Shannon entropy *H*), edges are the conditioning
        dependencies (annotated by mutual information *I*). Compare the **OLGA** bootstrap models to
        the **EM-learned** models fit to real out-of-frame repertoires.
        """
    )
    return


@app.cell
def _(list_bundled, mo):
    avail = list_bundled()
    source = mo.ui.dropdown(list(avail), value="olga", label="Model source")
    return avail, source


@app.cell
def _(avail, mo, source):
    loci = avail.get(source.value, [])
    locus = mo.ui.dropdown(loci, value=(loci[0] if loci else None), label="Locus")
    mo.hstack([source, locus], justify="start", gap=2)
    return (locus,)


@app.cell
def _(load_bundled, locus, mo, source):
    mo.stop(locus.value is None, mo.md("*No bundled models found — build them first.*"))
    model = load_bundled(locus.value, source.value)
    mo.md(f"**{locus.value}** · *{source.value}* · chain **{model.chain_type}** · "
          f"{model.genomic['genes_v'].height} V / {model.genomic['genes_j'].height} J genes")
    return (model,)


@app.cell
def _(analyze, mo, model, subprocess):
    # Render the Bayes net to inline SVG via the graphviz `dot` CLI (no python-graphviz dep).
    dot = analyze.bayes_net_dot(model)
    try:
        svg = subprocess.run(["dot", "-Tsvg"], input=dot, capture_output=True, text=True, check=True).stdout
        svg = svg[svg.index("<svg"):]  # drop the <?xml?> / <!DOCTYPE> preamble for clean inline embedding
        bn = mo.Html(svg)
    except (FileNotFoundError, subprocess.CalledProcessError):
        bn = mo.md("`dot` (graphviz) not found — showing the DOT source instead:").callout("warn")
        bn = mo.vstack([bn, mo.plain_text(dot)])
    mo.vstack([mo.md("### Bayes net &nbsp; *(node = event · H bits; edge = conditioning · I bits)*"), bn])
    return


@app.cell
def _(analyze, mo, model):
    ent = analyze.entropy_table(model)
    mi = analyze.mutual_information(model)
    mo.hstack([
        mo.vstack([mo.md("### Marginal entropy *H* (bits)"), ent]),
        mo.vstack([mo.md("### Mutual information *I* (bits)"), mi]),
    ], widths=[1, 1], gap=2)
    return


@app.cell
def _(mo, model):
    # Pick a marginal to plot (insertion length or a gene-usage table).
    plottable = [e for e in ("vd_ins", "dj_ins", "vj_ins", "dd_ins", "d_del", "n_d", "v_choice", "j_choice")
                 if e in model.tables]
    event = mo.ui.dropdown(plottable, value=plottable[0] if plottable else None, label="Marginal to plot")
    event
    return (event,)


@app.cell
def _(event, mo, model, plt):
    mo.stop(event.value is None)
    df = model.tables[event.value]
    fig, ax = plt.subplots(figsize=(7, 3))
    cols = df.columns
    if "length" in cols:
        ax.bar(df["length"], df["p"]); ax.set_xlabel("insertion length (nt)")
    elif "n_d" in cols:
        ax.bar(df["n_d"].cast(str), df["p"]); ax.set_xlabel("number of D segments")
    elif "ndel5" in cols:  # 2D deletion: marginal over 3'
        agg = df.group_by("ndel5").agg(__import__("polars").col("p").sum()).sort("ndel5")
        ax.bar(agg["ndel5"], agg["p"]); ax.set_xlabel("D 5' deletion (nt; neg = P-nt)")
    else:  # gene usage: top 15 by probability
        top = df.sort("p", descending=True).head(15)
        gcol = [c for c in cols if c.endswith("_allele")][-1]
        ax.barh(top[gcol].to_list()[::-1], top["p"].to_list()[::-1]); ax.set_xlabel("P")
    ax.set_ylabel("probability"); ax.set_title(f"{event.value}"); fig.tight_layout()
    ax
    return


@app.cell
def _(analyze, avail, load_bundled, locus, mo):
    # OLGA vs learned entropy side-by-side, when both are available for this locus.
    both = {s: load_bundled(locus.value, s) for s in ("olga", "learned")
            if locus.value in avail.get(s, [])}
    if len(both) == 2:
        import polars as _pl
        cmp = analyze.compare_entropy(both).with_columns(
            (_pl.col("learned") - _pl.col("olga")).round(3).alias("Δ (learned−olga)"))
        out = mo.vstack([mo.md("### OLGA vs learned — marginal entropy *H* (bits)"), cmp])
    else:
        out = mo.md("*Load a locus present in both `olga` and `learned` to see the comparison.*").callout()
    out
    return


if __name__ == "__main__":
    app.run()
