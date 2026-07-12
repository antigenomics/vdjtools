"""vdjtools.sc — Single-cell (AIRR Cell / 10x) paired-chain analysis.

Ingest 10x / AIRR-Cell contigs into a flat, ``cell_id``-keyed Rearrangement frame
(:mod:`~vdjtools.sc.read`), clean and pair chains with doublet / mispairing QC
(:mod:`~vdjtools.sc.pair`), and grade clonotype clusterings against a ground truth
(:mod:`~vdjtools.sc.cluster_eval`).
"""
from .cluster_eval import (
    assign_singleton_ids,
    cluster_eval,
    homogeneity,
    inverse_purity,
    normalized_inverse_purity,
    normalized_purity,
    parsimony,
    purity,
    q_measure,
)
from .pair import (
    chain_multiplicity,
    flag_mispairing,
    pair_chains,
    resolve_chains,
)
from .pgen import paired_pgen
from .anndata import to_anndata
from .read import read_10x, read_airr_cell, write_airr_cell

__all__ = [
    # ingestion
    "read_10x", "read_airr_cell", "write_airr_cell",
    # pairing / QC
    "resolve_chains", "pair_chains", "chain_multiplicity", "flag_mispairing",
    # paired-chain generation probability
    "paired_pgen",
    # cluster evaluation
    "cluster_eval", "purity", "normalized_purity", "inverse_purity",
    "normalized_inverse_purity", "homogeneity", "parsimony", "q_measure",
    "assign_singleton_ids",
    # scverse bridge
    "to_anndata",
]
