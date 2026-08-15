"""Single-cell interop: CellRanger in, scirpy / dandelion / scRepertoire out.

A marimo notebook. Launch it with

    marimo edit examples/single_cell_interop.py      # interactive editor
    marimo run  examples/single_cell_interop.py      # read-only served app

The downstream single-cell immune stack — scirpy, dandelion, scRepertoire — all read the
same interchange format: a flat AIRR Rearrangement table with `sequence_id` and `cell_id`.
This notebook shows vdjtools sitting in the middle of it: ingest 10x contigs, do the QC and
paired-Pgen work, hand the result to each downstream container, read it back, and push a
computed column onto an object someone else built.

Data auto-downloads from ``isalgo/airr_benchmark`` (folder ``dcode/``) into the HuggingFace
cache, or is read from a local ``./data_dump/`` if present. Needs the ``[sc]`` extra.
scirpy and dandelion are optional — every section degrades to a note if its library is
absent, and the *reading* directions (`from_scirpy`, `read_h5ddl`) need neither.
"""
import marimo

__generated_with = "0.23.14"
app = marimo.App(width="medium")


@app.cell
def _(mo):
    mo.md(
        r"""
        # Single-cell interop — one AIRR table, four ecosystems

        vdjtools is not trying to replace scirpy, dandelion or scRepertoire. It does the
        repertoire-model work those tools do not — generation probability, mispairing QC,
        clustering evaluation — and hands the result over.

        The seam is deliberately boring: **a flat AIRR Rearrangement table carrying
        `sequence_id` and `cell_id`**. That is what all three downstream tools read, so
        there is one emitter (`sc.to_airr`) and one inverse (`sc.from_airr`), and every
        bridge is a thin adapter on top.

        1. **ingest** CellRanger contigs,
        2. **QC + score** — chain multiplicity, pairing, paired Pgen,
        3. **out and back** — scirpy, dandelion, scRepertoire,
        4. **augment** — push `pgen_paired` onto an existing AnnData.
        """
    )
    return


@app.cell
def _():
    from pathlib import Path

    import marimo as mo
    import polars as pl

    from vdjtools import sc

    REPO_ID = "isalgo/airr_benchmark"
    FILE = "dcode/vdj_v1_hs_aggregated_donor4_{}.csv.gz"
    return FILE, Path, REPO_ID, mo, pl, sc


@app.cell
def _(FILE, Path, REPO_ID, mo, sc):
    def fetch(tag):
        """Local ./data_dump/ first (symlink your copies there), else HuggingFace."""
        local = Path("data_dump") / Path(FILE.format(tag)).name
        if local.exists():
            return local
        import huggingface_hub as hub
        return hub.hf_hub_download(REPO_ID, FILE.format(tag), repo_type="dataset")

    cells = sc.read_10x(fetch("all_contig_annotations"))
    mo.md(
        f"**{cells['cell_id'].n_unique():,} cells / {cells.height:,} productive contigs** "
        f"ingested from the 10x contig annotations.\n\n"
        "For CellRanger 4.0+ prefer `sc.read_airr_cell(\"outs/airr_rearrangement.tsv\")` — "
        "it is the very same file scirpy, dandelion and scRepertoire read, so ingestion "
        "and interop agree by construction. `read_10x` is here because this public dataset "
        "predates that output."
    )
    return (cells,)


@app.cell
def _(cells, mo, sc):
    quadrants = sc.chain_multiplicity(cells)
    mo.md(
        "## 1. QC — who is a usable cell?\n\n"
        "`chain_multiplicity` counts cells by `(n_light, n_heavy)`. `(1, 1)` is a clean "
        "single T cell; `(2, 1)` is a dual-alpha cell (real biology); anything with two "
        "heavies is a doublet candidate.\n\n"
        f"{mo.as_html(quadrants)}"
    )
    return


@app.cell
def _(cells, mo, sc):
    paired = sc.flag_mispairing(sc.pair_chains(cells, locus_pair="TRA_TRB"))
    scored = sc.paired_pgen(paired)
    ok = scored.filter(~scored["mispairing_flag"])
    mo.md(
        "## 2. Score — paired generation probability\n\n"
        "`paired_pgen` gives `Pgen(alpha) * Pgen(beta)` from the native V(D)J model. A low "
        "paired Pgen means a receptor unlikely to arise by chance — the usual prior for "
        "antigen-driven selection.\n\n"
        f"**{scored.height:,} receptors**, {int(scored['mispairing_flag'].sum()):,} flagged "
        f"as suspected mispairing / ambient alpha, leaving **{ok.height:,}**.\n\n"
        "NOTE: CellRanger reports *gene*-level V/J calls (`TRBV20-1`, not `TRBV20-1*01`) "
        "while the model is keyed by allele, so `paired_pgen` resolves each gene to its "
        "representative allele first. It does **not** fall back to marginalising over all "
        "V/J — that fallback once returned a Pgen 2.38x too high with no error. Pass "
        "`resolve_genes=False` to score only exact allele matches, or `condition_vj=False` "
        "to marginalise deliberately."
        f"\n\n{mo.as_html(scored.select(['cell_id', 'alpha_junction_aa', 'beta_junction_aa', 'pgen_paired']).head())}"
    )
    return paired, scored


@app.cell
def _(cells, mo, sc):
    airr = sc.to_airr(cells)
    mo.md(
        "## 3. The interchange table\n\n"
        "One row per contig. Two spellings are reconciled here and nowhere else: AIRR says "
        "`junction` where vdjtools says `junction_nt`, and scRepertoire reads "
        "`consensus_count` where scirpy and dandelion prefer `umi_count` — so both are "
        "emitted and one file feeds all three.\n\n"
        f"`{', '.join(airr.columns)}`\n\n{mo.as_html(airr.head(3))}"
    )
    return (airr,)


@app.cell
def _(cells, mo, sc):
    try:
        adata = sc.to_scirpy(cells)
        roundtrip = sc.from_scirpy(adata)
        scirpy_msg = (
            f"`to_scirpy` built an AnnData of **{adata.n_obs:,} cells** with chains under "
            f"`obsm['airr']` (keys: `{', '.join(adata.obsm.keys())}`), which is the layout "
            "scirpy >= 0.13 reads. It delegates to `scirpy.io.read_airr`, so the object "
            "carries *their* schema rather than a copy of it that could drift.\n\n"
            f"`from_scirpy` read it back to **{roundtrip.height:,} contigs** — and needs "
            "only `awkward`, not scirpy, so consuming somebody else's AnnData costs no "
            "extra install."
        )
    except ImportError:
        adata = None
        scirpy_msg = ("_scirpy is not installed (`pip install scirpy`), so this section is "
                      "a no-op. `from_scirpy` would still work on an AnnData handed to you._")
    mo.md("## 4. scirpy / scverse\n\n" + scirpy_msg)
    return (adata,)


@app.cell
def _(cells, mo, sc):
    try:
        vdj = sc.to_dandelion(cells)
        ddl_msg = (
            f"`to_dandelion` built a `Dandelion` — **{len(vdj.data):,} contigs** in `.data`, "
            f"**{len(vdj.metadata):,} cells** in `.metadata` (dandelion derives the cell "
            "table itself).\n\n"
            "Its `.h5ddl` format is plain HDF5, so `sc.read_h5ddl(path)` reads a dandelion "
            "result back **without dandelion installed** — useful when a collaborator sends "
            "you one."
        )
    except ImportError:
        ddl_msg = ("_dandelion is not installed (`pip install sc-dandelion`), so this "
                   "section is a no-op. `sc.read_h5ddl` would still read one sent to you._")
    mo.md("## 5. dandelion\n\n" + ddl_msg)
    return


@app.cell
def _(cells, mo, sc):
    sc.write_screpertoire(cells, "screpertoire_airr.tsv")
    mo.md(
        "## 6. scRepertoire (R)\n\n"
        "Export only — no R ships with vdjtools. Wrote `screpertoire_airr.tsv`:\n\n"
        "```r\n"
        'contigs  <- loadContigs("screpertoire_airr.tsv", format = "AIRR")\n'
        "combined <- combineTCR(contigs)\n"
        "seurat   <- combineExpression(combined, seurat)\n"
        "```\n\n"
        "WARNING: `combineTCR(samples=, ID=)` rewrites barcodes to `sample_ID_barcode`. "
        "That is the usual cause of a silent barcode-join failure against a Seurat "
        "`meta.data` — pass the same `samples`/`ID` to `combineExpression`, or omit both."
    )
    return


@app.cell
def _(adata, mo, sc, scored):
    if adata is None:
        aug_msg = "_Needs scirpy; see above._"
    else:
        # One row per cell: a cell with two receptor pairs would be ambiguous, so take the
        # most likely pair per cell before pushing.
        per_cell = (scored.sort("pgen_paired", descending=True, nulls_last=True)
                    .unique(subset=["cell_id"], keep="first")
                    .select(["cell_id", "pgen_paired"]))
        sc.push_obs(adata, per_cell)
        n_pushed = int(adata.obs["pgen_paired"].notna().sum())
        aug_msg = (
            f"`push_obs` wrote `pgen_paired` into `adata.obs` for **{n_pushed:,} of "
            f"{adata.n_obs:,} cells** (cells with no complete pair stay null).\n\n"
            "This is the direction that matters in practice: an analysis is already built "
            "around somebody's AnnData or Seurat object, and you want one vdjtools-computed "
            "column in it — not to move the whole analysis. `push_obs` also takes a "
            "`Dandelion` and writes to `.metadata`."
        )
    mo.md("## 7. Augmenting an object you did not build\n\n" + aug_msg)
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## Recap

        | Direction | Function | Needs |
        |---|---|---|
        | in | `read_airr_cell` / `read_10x` / `read_arda_cells` | — |
        | interchange | `to_airr` / `from_airr` / `write_airr` | — |
        | scirpy out | `to_scirpy` | `scirpy` |
        | scirpy in | `from_scirpy` | `awkward` |
        | dandelion out | `to_dandelion` | `sc-dandelion` |
        | dandelion in | `from_dandelion` / `read_h5ddl` | — / `h5py` |
        | scRepertoire out | `write_screpertoire` | — |
        | augment | `push_obs` | — |

        The asymmetry is deliberate: **writing** a container delegates to the library that
        owns it, so vdjtools never carries a stale copy of somebody's schema — but
        **reading** one is ours and stays cheap, so a result handed to you is always
        openable.

        See `docs/singlecell.rst` for the full guide, and `examples/single_cell.py` for the
        analysis side (antigen-driven clustering graded with `cluster_eval`).
        """
    )
    return


if __name__ == "__main__":
    app.run()
