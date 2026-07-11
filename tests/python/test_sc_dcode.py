"""Real-data single-cell benchmark: TCRnet-style clustering on 10x dCODE donor4.

Fetches the dCODE donor4 10x dextramer dataset (``isalgo/airr_benchmark``) via the
``hf`` fixture and runs the full single-cell path: ``read_10x`` → ``resolve_chains`` →
per-cell top-β CDR3, antigen labels from the CITE-seq binarized matrix, 1-substitution
CDR3 connected components as the predicted clustering, then
:func:`vdjtools.sc.cluster_eval`. Skips cleanly offline / without ``huggingface_hub``
(see conftest), so the default suite stays green with no network.

Asserted directionally (the data is biology, not a fixed number): the antigen-driven
clustering has high purity, and a random relabelling collapses the *normalized* scores.
"""
from __future__ import annotations

import random

import polars as pl
import pytest

from vdjtools import sc
from vdjtools.sc.cluster_eval import cluster_eval

REPO = "isalgo/airr_benchmark"
_FILE = "dcode/vdj_v1_hs_aggregated_donor4_{}.csv.gz"

# The two epitopes with enough labelled cells to cluster (EBNA-3B and IE-1).
EPITOPES = ("A1101_IVTDFSVIK_EBNA-3B_EBV", "A0301_KLGGALQAK_IE-1_CMV")
_SORT = ["duplicate_count", "umi_count", "sequence_id"]
_DESC = [True, True, False]


def _fetch_all(hf):
    """Fetch the three donor4 files (all_contig, consensus, binarized_matrix)."""
    return (
        hf(REPO, _FILE.format("all_contig_annotations")),
        hf(REPO, _FILE.format("consensus_annotations")),
        hf(REPO, _FILE.format("binarized_matrix")),
    )


def _antigen_labels(matrix_path) -> pl.DataFrame:
    """Barcode → antigen: the single True ``*_binder`` column (stripped), else unassigned."""
    mat = pl.read_csv(matrix_path, infer_schema_length=0)
    binder = [c for c in mat.columns if c.endswith("_binder")]
    n_true = pl.sum_horizontal([(pl.col(c) == "True").cast(pl.Int64) for c in binder])
    picks = [pl.when(pl.col(c) == "True").then(pl.lit(c[:-len("_binder")])).alias(f"_l{i}")
             for i, c in enumerate(binder)]
    return (mat.select("barcode", n_true.alias("_ntrue"), *picks)
            .with_columns(
                pl.when(pl.col("_ntrue") == 1)
                .then(pl.coalesce([f"_l{i}" for i in range(len(binder))]))
                .otherwise(pl.lit("unassigned")).alias("antigen"))
            .select("barcode", "antigen"))


def _components(seqs: list[str], scope: str = "1,0,0,1") -> list[int]:
    """1-substitution connected-component id per input sequence (singletons = own id)."""
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


def test_dcode_quadrants_and_pairing_invariants(hf):
    """Dominant quadrant is (TRA+, TRB+); Cartesian pairing invariants hold on real data."""
    all_contig, consensus, _ = _fetch_all(hf)
    contigs = sc.read_10x(all_contig, consensus)
    assert contigs.height > 0

    # Dominant presence quadrant on the full donor is (TRA+, TRB+).
    quad = sc.chain_multiplicity(contigs).with_columns(
        pl.when(pl.col("n_light") > 0).then(pl.lit("TRA+")).otherwise(pl.lit("TRA-")).alias("tra"),
        pl.when(pl.col("n_heavy") > 0).then(pl.lit("TRB+")).otherwise(pl.lit("TRB-")).alias("trb"),
    ).group_by("tra", "trb").agg(pl.col("cell_count").sum()).sort("cell_count", descending=True)
    top = quad.row(0, named=True)
    assert (top["tra"], top["trb"]) == ("TRA+", "TRB+")

    # Pairing invariants on a resolved slice (keep runtime bounded).
    slice_cells = contigs["cell_id"].unique().head(4000).to_list()
    resolved = sc.resolve_chains(contigs.filter(pl.col("cell_id").is_in(slice_cells)))
    paired = sc.pair_chains(resolved, resolve=False)

    per = resolved.group_by("cell_id").agg(
        (pl.col("locus") == "TRA").sum().alias("na"),
        (pl.col("locus") == "TRB").sum().alias("nb"),
    )
    # (a) 2 alpha + 1 beta -> exactly 2 pairs for that cell.
    dual = per.filter((pl.col("na") == 2) & (pl.col("nb") == 1))
    assert dual.height > 0
    cid = dual["cell_id"][0]
    assert paired.filter(pl.col("cell_id") == cid).height == 2
    # (b) incomplete cells counted-not-emitted: emitted == cells with both sides.
    complete = set(per.filter((pl.col("na") >= 1) & (pl.col("nb") >= 1))["cell_id"].to_list())
    assert set(paired["cell_id"].unique().to_list()) == complete
    # (c) total pairs == sum of the per-cell Cartesian product na*nb.
    exp = per.filter((pl.col("na") >= 1) & (pl.col("nb") >= 1)).select(
        (pl.col("na") * pl.col("nb")).sum()).item()
    assert paired.height == exp


def test_dcode_tcrnet_clustering_purity(hf, capsys):
    """Antigen-driven CDR3 clustering: high purity; shuffled labels collapse normalized scores."""
    pytest.importorskip("vdjmatch")
    pytest.importorskip("seqtree")
    all_contig, consensus, matrix_path = _fetch_all(hf)

    contigs = sc.read_10x(all_contig, consensus)
    labels = _antigen_labels(matrix_path)

    # Restrict to target-labelled cells, then resolve to the per-cell top-beta CDR3.
    targets = set(labels.filter(pl.col("antigen").is_in(EPITOPES))["barcode"].to_list())
    resolved = sc.resolve_chains(contigs.filter(pl.col("cell_id").is_in(list(targets))))
    beta = (resolved.filter(pl.col("locus") == "TRB")
            .sort(_SORT, descending=_DESC)
            .group_by("cell_id", maintain_order=True).first())
    cell_beta = (beta.join(labels, left_on="cell_id", right_on="barcode", how="left")
                 .filter(pl.col("antigen").is_in(EPITOPES) & pl.col("cdr3_aa").is_not_null()))
    assert cell_beta.height > 500  # enough cells to be meaningful

    labels_true = cell_beta["antigen"].to_list()
    labels_pred = _components(cell_beta["cdr3_aa"].to_list())

    real = cluster_eval(labels_true, labels_pred)
    shuffled = labels_true[:]
    random.Random(0).shuffle(shuffled)
    rand = cluster_eval(shuffled, labels_pred)

    # Emit a small results table (visible with -s).
    table = pl.DataFrame({
        "metric": ["purity", "normalized_purity", "homogeneity", "parsimony", "q_measure"],
        "antigen_labels": [real[k] for k in
                           ("purity", "normalized_purity", "homogeneity", "parsimony", "q_measure")],
        "shuffled_labels": [rand[k] for k in
                            ("purity", "normalized_purity", "homogeneity", "parsimony", "q_measure")],
    })
    with capsys.disabled():
        print(f"\ndCODE donor4: {cell_beta.height} labelled beta cells, "
              f"epitopes={cell_beta['antigen'].n_unique()}")
        print(table)

    # Antigen-driven clustering is genuinely pure and coherent.
    assert real["purity"] > 0.7
    assert real["q_measure"] > 0.5
    # Random relabelling collapses the normalized / information-theoretic scores.
    assert rand["normalized_purity"] < real["normalized_purity"] - 0.3
    assert rand["homogeneity"] < real["homogeneity"] - 0.3
    assert rand["q_measure"] < real["q_measure"] - 0.3
