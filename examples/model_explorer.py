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


@app.cell
def _(mo):
    mo.md(
        """
        ## Annotating VDJdb — the missing nucleotides

        A VDJdb record carries `(V, J, CDR3aa)` and **no nucleotide sequence**, so none of the
        boundary markup a repertoire analysis wants is there. `model.infer_nt` reconstructs it:
        germline positions are pinned to their V/D/J segment, and every free N-region position takes
        the nucleotide the corresponding insertion model prefers
        (`P(nt₁)·∏P(nt_k | nt_{k−1})`). The result is a full annotation — `cdr3_nt`, the V/D/J
        spans, and the exact `pgen_nt` of the reconstruction.
        """
    ).callout()
    return


@app.cell
def _(mo):
    vdjdb_n = mo.ui.slider(10, 200, value=40, step=10, label="VDJdb records to annotate")
    vdjdb_n
    return (vdjdb_n,)


@app.cell
def _(locus, mo, pl):
    def vdjdb_records(locus_name, n):
        """`(cdr3_aa, v_call, j_call)` from a local ./data_dump/vdjdb.slim.txt(.gz) if present, else
        from the latest antigenomics/vdjdb-db release (cached into ./data_dump/)."""
        import io
        import json
        import urllib.request
        import zipfile
        from pathlib import Path

        d = Path("data_dump")
        src = next((d / n_ for n_ in ("vdjdb.slim.txt", "vdjdb.slim.txt.gz") if (d / n_).exists()),
                   None)
        if src is None:
            d.mkdir(exist_ok=True)
            rel = json.load(urllib.request.urlopen(
                "https://api.github.com/repos/antigenomics/vdjdb-db/releases/latest", timeout=30))
            url = next(a["browser_download_url"] for a in rel["assets"]
                       if a["name"].endswith(".zip"))
            z = zipfile.ZipFile(io.BytesIO(urllib.request.urlopen(url, timeout=180).read()))
            name = next(n_ for n_ in z.namelist()
                        if n_.endswith(("vdjdb.slim.txt", "vdjdb.slim.txt.gz")))
            src = d / Path(name).name
            src.write_bytes(z.read(name))
        gene = {"TRA": "TRA", "TRB": "TRB"}.get(locus_name)
        if gene is None:
            return None
        return (pl.read_csv(src, separator="\t", infer_schema_length=0)
                .filter((pl.col("species") == "HomoSapiens") & (pl.col("gene") == gene))
                .select(cdr3_aa=pl.col("cdr3"), v_call=pl.col("v.segm"), j_call=pl.col("j.segm"))
                .unique(subset="cdr3_aa", maintain_order=True)
                .head(n))

    try:
        recs = vdjdb_records(locus.value, 2000)
        recs_err = None
    except Exception as exc:                     # offline, or no release asset — say so, don't fail
        recs, recs_err = None, str(exc)
    mo.md(f"*VDJdb unavailable: {recs_err}*").callout("warn") if recs_err else None
    return (recs,)


@app.cell
def _(mo, model, pl, recs, vdjdb_n):
    # NOTE: VDJdb stores GENE-level calls ("TRBV9"); the model is keyed by ALLELE, and passing a gene
    # name raises by design. Resolve to the model's unique allele of that gene, skip the ambiguous.
    if recs is None:
        ann = None
        out2 = mo.md("*No VDJdb records for this locus — TRA/TRB only.*").callout()
    else:
        import time

        from vdjtools.model import infer_nt

        # NOTE: pass the Model, not a prepare()-d one — that selects the pure-Python reference
        # search, which is ~600x slower on a VDJ locus.
        by_gene = {}
        for seg in ("v", "j"):
            for a in model.genomic[f"genes_{seg}"][f"{seg}_allele"]:
                by_gene.setdefault((seg, a.split("*")[0]), []).append(a)
        rows, t0 = [], time.perf_counter()
        for r in recs.head(vdjdb_n.value).iter_rows(named=True):
            va = by_gene.get(("v", r["v_call"].split("*")[0]), [])
            ja = by_gene.get(("j", r["j_call"].split("*")[0]), [])
            if len(va) != 1 or len(ja) != 1:
                continue
            sc = infer_nt(model, r["cdr3_aa"], va[0], ja[0])
            if sc is None:
                continue
            rows.append({"cdr3_aa": r["cdr3_aa"], "cdr3_nt": sc.cdr3_nt, "v_call": sc.v_call,
                         "d_call": sc.d_call, "j_call": sc.j_call, "v_end": sc.v_end,
                         "d_start": sc.d_start, "d_end": sc.d_end, "j_start": sc.j_start,
                         "pgen": sc.pgen, "margin": round(sc.margin, 2)})
        ms = (time.perf_counter() - t0) / max(len(rows), 1) * 1000
        ann = pl.DataFrame(rows)
        out2 = mo.vstack([
            mo.md(f"### {len(ann)} VDJdb records annotated · {ms:.1f} ms/record"),
            mo.md("`margin` = Pgen of the winner over the runner-up. A margin near 1 means the "
                  "N-region is long enough that several nucleotide sequences are near-equally "
                  "likely — the reconstruction is a *most likely* one, not *the* one."),
            ann,
        ])
    out2
    return (ann,)


@app.cell
def _(ann, mo, plt):
    # Where the model puts the segment boundaries across the annotated set.
    if ann is None or ann.is_empty():
        out3 = mo.md("*Nothing annotated yet.*").callout()
    else:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.2))
        ax1.hist(ann["pgen"].log10().to_list(), bins=30, color="#0072B2")
        ax1.set_xlabel("log10 Pgen of the reconstruction"); ax1.set_ylabel("records")
        ax2.hist(ann["v_end"].to_list(), bins=range(0, 25), alpha=0.7,
                 color="#009E73", label="V end")
        ax2.hist([len(s) - t for s, t in zip(ann["cdr3_nt"], ann["j_start"])], bins=range(0, 25),
                 alpha=0.7, color="#D55E00", label="J length")
        ax2.set_xlabel("nt contributed by germline"); ax2.legend()
        fig.tight_layout()
        out3 = fig
    out3
    return


if __name__ == "__main__":
    app.run()
