"""Single-cell paired-chain TCR analysis with vdjtools v2.

A marimo notebook. Launch it with

    marimo edit examples/single_cell.py      # interactive editor
    marimo run  examples/single_cell.py      # read-only served app

10x single-cell VDJ gives paired α/β chains per cell, and — in a dCODE dextramer
experiment — an antigen label per cell. `vdjtools.sc` ingests the 10x contigs, cleans and
pairs the chains with doublet / mispairing QC, and grades a clonotype clustering against a
ground truth. Here we run the full path on the public **dCODE donor 4** dataset: ingest →
chain-multiplicity QC → resolve/pair → antigen-driven CDR3 clustering →
`cluster_eval`, and show that a 1-substitution clustering of the β CDR3s is genuinely
pure for the antigen labels (and that a random relabelling collapses the scores).

Data auto-downloads from ``isalgo/airr_benchmark`` (folder ``dcode/``) into the HuggingFace
cache. Needs the ``[sc]`` and ``[overlap]`` extras (``scikit-learn``, ``vdjmatch``, ``seqtree``).
"""
import marimo

__generated_with = "0.23.14"
app = marimo.App(width="medium")


@app.cell
def _(mo):
    mo.md(
        r"""
        # Single-cell paired-chain TCR analysis — dCODE donor 4

        A 10x single-cell VDJ run resolves both chains of a T-cell receptor **per cell**;
        a dCODE dextramer panel adds an **antigen label** per cell. This notebook takes the
        raw 10x contigs of dCODE donor 4 through `vdjtools.sc`:

        1. **ingest** the contig + consensus annotations (`read_10x`),
        2. **QC** cell composition (`chain_multiplicity` — the TRA/TRB presence quadrants),
        3. **resolve & pair** chains into receptors (`resolve_chains`, `pair_chains`),
        4. take each cell's dominant **β CDR3**, cluster by **1 substitution**, and
        5. grade that clustering against the dextramer antigen labels (`cluster_eval`).

        The test question: does an *unsupervised* CDR3 clustering recover the
        antigen-specific groups?
        """
    )
    return


@app.cell
def _():
    from pathlib import Path

    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import polars as pl

    from vdjtools import sc
    from vdjtools.sc.cluster_eval import cluster_eval

    REPO_ID = "isalgo/airr_benchmark"
    FILE = "dcode/vdj_v1_hs_aggregated_donor4_{}.csv.gz"
    # The two epitopes with enough labelled cells to cluster (EBV EBNA-3B, CMV IE-1).
    EPITOPES = ("A1101_IVTDFSVIK_EBNA-3B_EBV", "A0301_KLGGALQAK_IE-1_CMV")
    OKABE = {"blue": "#0072B2", "vermillion": "#D55E00", "green": "#009E73",
             "orange": "#E69F00", "purple": "#CC79A7", "grey": "#8C8C8C"}
    return EPITOPES, FILE, OKABE, Path, REPO_ID, cluster_eval, mo, np, pl, plt, sc


@app.cell
def _(FILE, REPO_ID, mo, pl):
    def fetch(tag):
        import huggingface_hub as hub
        return hub.hf_hub_download(REPO_ID, FILE.format(tag), repo_type="dataset")

    def antigen_labels(matrix_path):
        """Barcode → antigen: the single True ``*_binder`` column (stripped), else unassigned."""
        mat = pl.read_csv(matrix_path, infer_schema_length=0)
        binder = [c for c in mat.columns if c.endswith("_binder")]
        n_true = pl.sum_horizontal([(pl.col(c) == "True").cast(pl.Int64) for c in binder])
        picks = [pl.when(pl.col(c) == "True").then(pl.lit(c[:-len("_binder")])).alias(f"_l{i}")
                 for i, c in enumerate(binder)]
        return (mat.select("barcode", n_true.alias("_n"), *picks).with_columns(
            pl.when(pl.col("_n") == 1)
            .then(pl.coalesce([f"_l{i}" for i in range(len(binder))]))
            .otherwise(pl.lit("unassigned")).alias("antigen")).select("barcode", "antigen"))

    def components(seqs, scope="1,0,0,1"):
        """1-substitution connected-component id per input CDR3 (singletons = own id)."""
        import vdjmatch.cluster as vc
        uniq = sorted(set(seqs))
        idx = {s: i for i, s in enumerate(uniq)}
        parent = list(range(len(uniq)))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        pairs = vc.overlap(uniq, scope=scope)
        for a, b in zip(pairs["a_idx"].to_list(), pairs["b_idx"].to_list()):
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb
        comp = {s: find(idx[s]) for s in uniq}
        return [comp[s] for s in seqs]

    return antigen_labels, components, fetch


@app.cell
def _(fetch, antigen_labels, mo, sc):
    _all_contig = fetch("all_contig_annotations")
    _consensus = fetch("consensus_annotations")
    contigs = sc.read_10x(_all_contig, _consensus)
    labels = antigen_labels(fetch("binarized_matrix"))
    mo.md(f"Ingested **{contigs.height:,} contigs** across "
          f"**{contigs['cell_id'].n_unique():,} cells**; "
          f"{labels.filter(labels['antigen'] != 'unassigned').height:,} antigen-labelled barcodes.")
    return contigs, labels


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 1 · Cell-composition QC — the TRA/TRB quadrants

        `chain_multiplicity` counts each cell's light (TRA) and heavy (TRB) chains. A clean
        single T cell is **(TRA+, TRB+)** with one of each; TRB-only / TRA-only / multi-chain
        cells are dropouts or doublets. The dominant quadrant should be (TRA+, TRB+).
        """
    )
    return


@app.cell
def _(OKABE, contigs, pl, plt, sc):
    _q = sc.chain_multiplicity(contigs).with_columns(
        pl.when(pl.col("n_light") > 0).then(pl.lit("TRA+")).otherwise(pl.lit("TRA−")).alias("tra"),
        pl.when(pl.col("n_heavy") > 0).then(pl.lit("TRB+")).otherwise(pl.lit("TRB−")).alias("trb"),
    ).group_by("tra", "trb").agg(pl.col("cell_count").sum().alias("cells")).sort("cells", descending=True)

    _lab = [f"{r['tra']}, {r['trb']}" for r in _q.iter_rows(named=True)]
    _val = _q["cells"].to_list()
    figq, axq = plt.subplots(figsize=(6.8, 4.0))
    axq.bar(_lab, _val, color=[OKABE["green"] if v == max(_val) else OKABE["grey"] for v in _val])
    for _i, _v in enumerate(_val):
        axq.text(_i, _v, f"{_v:,}", ha="center", va="bottom", fontsize=9)
    axq.set(ylabel="cells", title="Chain-presence quadrants (dominant = paired T cell)")
    axq.spines[["top", "right"]].set_visible(False)
    figq.tight_layout()
    figq
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 2 · Resolve & pair chains

        `resolve_chains` picks each cell's dominant chains; `pair_chains` forms α/β receptors
        (a 2-alpha / 1-beta cell yields two Cartesian pairs, incomplete cells are dropped).
        """
    )
    return


@app.cell
def _(EPITOPES, contigs, labels, mo, pl, sc):
    # Restrict to target-antigen-labelled cells, then resolve to the per-cell top-β CDR3.
    _targets = set(labels.filter(pl.col("antigen").is_in(EPITOPES))["barcode"].to_list())
    resolved = sc.resolve_chains(contigs.filter(pl.col("cell_id").is_in(list(_targets))))
    paired = sc.pair_chains(resolved, resolve=False)

    _beta = (resolved.filter(pl.col("locus") == "TRB")
             .sort(["duplicate_count", "umi_count", "sequence_id"], descending=[True, True, False])
             .group_by("cell_id", maintain_order=True).first())
    cell_beta = (_beta.join(labels, left_on="cell_id", right_on="barcode", how="left")
                 .filter(pl.col("antigen").is_in(EPITOPES) & pl.col("cdr3_aa").is_not_null()))
    mo.md(f"Resolved **{resolved['cell_id'].n_unique():,} cells** → "
          f"**{paired.height:,} α/β pairs**; **{cell_beta.height:,} cells** carry a "
          f"target-antigen label and a β CDR3 ({cell_beta['antigen'].n_unique()} epitopes).")
    return (cell_beta,)


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 3 · Unsupervised β-CDR3 clustering vs the antigen labels

        Cluster the per-cell β CDR3s into **1-substitution connected components** (blind to
        the antigen labels), then grade that clustering against the true epitope labels with
        `cluster_eval`. A random relabelling is the null. High purity / q-measure for the
        real labels — and a collapse for the shuffled ones — means the receptor sequence
        alone recovers the antigen groups.
        """
    )
    return


@app.cell
def _(OKABE, cell_beta, cluster_eval, components, mo, np, pl, plt):
    _true = cell_beta["antigen"].to_list()
    _pred = components(cell_beta["cdr3_aa"].to_list())
    real = cluster_eval(_true, _pred)
    _shuf = _true[:]
    np.random.default_rng(0).shuffle(_shuf)
    rand = cluster_eval(_shuf, _pred)

    _keys = ["purity", "normalized_purity", "homogeneity", "parsimony", "q_measure"]
    table = pl.DataFrame({"metric": _keys,
                          "antigen_labels": [round(real[k], 3) for k in _keys],
                          "shuffled_labels": [round(rand[k], 3) for k in _keys]})

    _x = np.arange(len(_keys))
    figc, axc = plt.subplots(figsize=(8.2, 4.4))
    axc.bar(_x - 0.2, [real[k] for k in _keys], 0.4, label="antigen labels", color=OKABE["green"])
    axc.bar(_x + 0.2, [rand[k] for k in _keys], 0.4, label="shuffled (null)", color=OKABE["grey"])
    axc.set_xticks(_x); axc.set_xticklabels(_keys, rotation=20, ha="right")
    axc.set(ylabel="score", title="1-substitution CDR3 clustering vs antigen labels")
    axc.legend(frameon=False)
    axc.spines[["top", "right"]].set_visible(False)
    figc.tight_layout()
    mo.vstack([table, figc])
    return real, table


@app.cell
def _(mo, real):
    mo.md(
        f"""
        ---
        **Takeaway.** An unsupervised 1-substitution clustering of the β CDR3s is genuinely
        antigen-coherent — purity **{real['purity']:.2f}**, q-measure **{real['q_measure']:.2f}** —
        and a shuffled labelling collapses the normalized scores, confirming the signal is
        real. The whole path is `vdjtools.sc`: `read_10x` → `chain_multiplicity` →
        `resolve_chains`/`pair_chains` → `cluster_eval`, with the 1-substitution graph from
        `vdjmatch`. Paired receptors also export straight to the scverse world via
        `vdjtools.sc.to_anndata`.
        """
    )
    return


if __name__ == "__main__":
    app.run()
