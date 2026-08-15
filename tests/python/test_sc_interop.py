"""Round-trip tests against the real downstream containers (scirpy, dandelion).

Each tool is skipped if unimportable (see `_need`), so this file is optional -- the
non-skippable contract for the interchange format itself lives in `test_sc_airr.py`. What is pinned here is that the
containers those tools actually build survive a round-trip through vdjtools unchanged, and
that `push_obs` attaches computed columns to both container shapes.
"""
import polars as pl
import pytest

from vdjtools import sc
from vdjtools.sc.read import SC_COLUMNS


def _cells():
    """Two cells, one with a paired alpha/beta, one beta-only."""
    return pl.DataFrame({
        "cell_id": ["c1", "c1", "c2"],
        "sequence_id": ["c1_contig_1", "c1_contig_2", "c2_contig_1"],
        "locus": ["TRA", "TRB", "TRB"],
        "v_call": ["TRAV1-2", "TRBV20-1", "TRBV9"],
        "d_call": [None, "TRBD1", None],
        "j_call": ["TRAJ33", "TRBJ2-7", "TRBJ2-3"],
        "c_call": ["TRAC", "TRBC2", "TRBC2"],
        "junction_aa": ["CAVRDSNYQLIW", "CASSLGQAYEQYF", "CASSVDTQYF"],
        "junction_nt": ["TGTGCC", "TGTGCCAGC", "TGTGCCAGCAGC"],
        "duplicate_count": [120, 340, 88],
        "umi_count": [6, 14, 4],
        "clone_id": ["clonotype1", "clonotype1", "clonotype2"],
        "productive": [True, True, True],
    }).select(SC_COLUMNS)



def _need(name):
    """Import an optional tool, skipping if it is unusable for ANY import reason.

    Not `pytest.importorskip`: that only skips when the named module is missing, and
    re-raises when the module exists but its own dependency chain is broken -- which is the
    normal state of these packages (dandelion pulls nxviz, which breaks on matplotlib >=3.9).
    The [interop] extra is best-effort by design, so an unimportable tool must SKIP, not fail.
    """
    import importlib

    try:
        return importlib.import_module(name)
    except ImportError as e:                     # incl. a transitive ImportError
        pytest.skip(f"{name} not importable: {e}")


def _sorted(df):
    return df.sort("sequence_id")


# --------------------------------------------------------------------------- scirpy

def test_scirpy_round_trip_is_lossless():
    _need("scirpy")
    cells = _cells()
    assert _sorted(sc.from_scirpy(sc.to_scirpy(cells))).equals(_sorted(cells))


def test_to_scirpy_builds_the_obsm_airr_layout():
    _need("scirpy")
    adata = sc.to_scirpy(_cells())
    # scirpy >=0.13 keys chains off obsm["airr"], one obs row per CELL (not per contig).
    assert "airr" in adata.obsm
    assert adata.n_obs == 2
    # index_chains=True is what makes the object immediately usable by scirpy's tools.
    assert "chain_indices" in adata.obsm


def test_to_scirpy_can_skip_chain_indexing():
    _need("scirpy")
    adata = sc.to_scirpy(_cells(), index_chains=False)
    assert "chain_indices" not in adata.obsm


def test_from_scirpy_needs_no_scirpy_import(monkeypatch):
    """The read direction must work from a plain [sc] install (awkward only)."""
    _need("scirpy")
    adata = sc.to_scirpy(_cells())
    monkeypatch.setitem(__import__("sys").modules, "scirpy", None)
    assert sc.from_scirpy(adata).height == 3


def test_from_scirpy_rejects_an_object_without_airr():
    pytest.importorskip("anndata")
    import anndata as ad
    import numpy as np

    with pytest.raises(KeyError, match="obsm"):
        sc.from_scirpy(ad.AnnData(X=np.zeros((2, 1), dtype="float32")))


def test_mudata_wrapping_puts_vdj_in_the_airr_modality():
    _need("scirpy")
    _need("mudata")
    import anndata as ad
    import numpy as np

    adata = sc.to_scirpy(_cells())
    gex = ad.AnnData(X=np.zeros((adata.n_obs, 3), dtype="float32"))
    gex.obs_names = adata.obs_names
    mdata = sc.to_scirpy(_cells(), gex=gex)
    assert set(mdata.mod) == {"gex", "airr"}
    # from_scirpy unwraps the modality itself.
    assert sc.from_scirpy(mdata).height == 3


# ------------------------------------------------------------------------ dandelion

def test_dandelion_round_trip_is_lossless():
    _need("dandelion")
    cells = _cells()
    assert _sorted(sc.from_dandelion(sc.to_dandelion(cells))).equals(_sorted(cells))


def test_to_dandelion_builds_contig_and_cell_tables():
    _need("dandelion")
    vdj = sc.to_dandelion(_cells())
    assert len(vdj.data) == 3          # contig level
    assert len(vdj.metadata) == 2      # cell level


def test_read_h5ddl_needs_no_dandelion(tmp_path, monkeypatch):
    """.h5ddl is plain HDF5, so a dandelion result is readable without dandelion."""
    _need("dandelion")
    pytest.importorskip("h5py")
    out = tmp_path / "vdj.h5ddl"
    sc.to_dandelion(_cells()).write_h5ddl(out)

    monkeypatch.setitem(__import__("sys").modules, "dandelion", None)
    back = sc.read_h5ddl(out)
    assert back.columns == SC_COLUMNS
    assert _sorted(back)["junction_aa"].to_list() == _sorted(_cells())["junction_aa"].to_list()


def test_cell_id_is_recovered_from_the_contig_naming():
    """dandelion derives cell_id from `<cell>_contig_<n>`; so must we."""
    from vdjtools.sc.dandelion import _ensure_cell_id

    df = pl.DataFrame({"sequence_id": ["AAACCTGAGAAACCAT-1_contig_2"]})
    assert _ensure_cell_id(df)["cell_id"].to_list() == ["AAACCTGAGAAACCAT-1"]


# -------------------------------------------------------------------------- augment

def test_push_obs_attaches_columns_to_an_anndata():
    _need("scirpy")
    adata = sc.to_scirpy(_cells())
    scores = pl.DataFrame({"cell_id": ["c1", "c2"], "pgen_paired": [1e-9, 4e-11]})
    sc.push_obs(adata, scores)
    assert adata.obs.loc["c1", "pgen_paired"] == pytest.approx(1e-9)
    assert adata.obs.loc["c2", "pgen_paired"] == pytest.approx(4e-11)


def test_push_obs_leaves_unmentioned_cells_null():
    _need("scirpy")
    import numpy as np

    adata = sc.to_scirpy(_cells())
    sc.push_obs(adata, pl.DataFrame({"cell_id": ["c1"], "flag": [1.0]}))
    assert adata.obs.loc["c1", "flag"] == 1.0
    assert np.isnan(adata.obs.loc["c2", "flag"])


def test_push_obs_attaches_to_a_dandelion_metadata():
    _need("dandelion")
    vdj = sc.to_dandelion(_cells())
    sc.push_obs(vdj, pl.DataFrame({"cell_id": ["c1", "c2"], "pgen_paired": [1e-9, 4e-11]}))
    assert "pgen_paired" in vdj.metadata.columns


def test_push_obs_rejects_a_multi_pair_frame():
    _need("scirpy")
    adata = sc.to_scirpy(_cells())
    dup = pl.DataFrame({"cell_id": ["c1", "c1"], "x": [1.0, 2.0]})
    with pytest.raises(ValueError, match="not unique"):
        sc.push_obs(adata, dup)


def test_push_obs_rejects_a_target_with_neither_obs_nor_metadata():
    with pytest.raises(TypeError, match="neither"):
        sc.push_obs(object(), pl.DataFrame({"cell_id": ["c1"], "x": [1.0]}))
