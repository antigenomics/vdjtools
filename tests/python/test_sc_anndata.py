"""Tests for the single-cell AnnData bridge (vdjtools.sc.to_anndata).

Pins that paired-receptor output lands in a valid AnnData: obs indexed by the
unique pair_id, cell_id preserved as a column, empty X by default, and that a
supplied expression matrix is attached. anndata (the [sc] extra) is importorskip'd.
"""
import polars as pl
import pytest

pytest.importorskip("anndata")

from vdjtools import sc  # noqa: E402


def _paired():
    """A minimal pair_chains-shaped frame: two cells, one of which has two pairs."""
    return pl.DataFrame({
        "cell_id": ["c1", "c2", "c2"],
        "pair_id": ["c1_1", "c2_1", "c2_2"],        # unique per receptor pair
        "alpha_v_call": ["TRAV1", "TRAV2", "TRAV3"],
        "alpha_j_call": ["TRAJ1", "TRAJ2", "TRAJ3"],
        "alpha_cdr3_aa": ["CAAF", "CABF", "CACF"],
        "alpha_umi_count": [10, 8, 3], "alpha_duplicate_count": [50, 40, 12],
        "beta_v_call": ["TRBV1", "TRBV2", "TRBV2"],
        "beta_j_call": ["TRBJ1", "TRBJ2", "TRBJ2"],
        "beta_cdr3_aa": ["CASSF", "CASTF", "CASTF"],
        "beta_umi_count": [20, 15, 15], "beta_duplicate_count": [80, 60, 60],
    })


def test_to_anndata_obs_shape_and_index():
    ad_obj = sc.to_anndata(_paired())
    assert ad_obj.n_obs == 3 and ad_obj.n_vars == 0        # pure VDJ container
    assert list(ad_obj.obs_names) == ["c1_1", "c2_1", "c2_2"]
    # cell_id preserved (the dual-pair cell c2 appears twice) for GEX alignment.
    assert ad_obj.obs["cell_id"].tolist() == ["c1", "c2", "c2"]
    assert ad_obj.obs["beta_cdr3_aa"].tolist() == ["CASSF", "CASTF", "CASTF"]


def test_to_anndata_attaches_expression_matrix():
    import numpy as np
    X = np.arange(3 * 4, dtype="float32").reshape(3, 4)   # 3 cells x 4 genes
    ad_obj = sc.to_anndata(_paired(), X=X)
    assert ad_obj.shape == (3, 4)
    assert np.allclose(ad_obj.X, X)


def test_to_anndata_rejects_nonunique_index():
    dup = _paired().with_columns(pl.lit("same").alias("pair_id"))
    with pytest.raises(ValueError, match="not unique"):
        sc.to_anndata(dup)


def test_to_anndata_index_by_cell_id_when_unique():
    single = _paired().head(1)                              # one cell, one pair
    ad_obj = sc.to_anndata(single, index="cell_id")
    assert list(ad_obj.obs_names) == ["c1"]
