"""Chain resolution, paired-receptor assembly, and mispairing / ambient QC."""
from __future__ import annotations

import polars as pl
import pytest

from vdjtools import sc
from vdjtools.sc.read import SC_COLUMNS

_INT = {"duplicate_count", "umi_count"}


def _sc(rows):
    """Build a canonical sc long frame from partial dict rows."""
    full = []
    for r in rows:
        base = {c: None for c in SC_COLUMNS}
        base.update(r)
        full.append(base)
    schema = {c: (pl.Int64 if c in _INT else pl.Utf8) for c in SC_COLUMNS}
    return pl.DataFrame(full, schema=schema)


def _chain(cell, seq, locus, cdr3, dup, umi, v=None, j=None):
    return dict(cell_id=cell, sequence_id=seq, locus=locus, junction_aa=cdr3,
                duplicate_count=dup, umi_count=umi,
                v_call=v or f"{locus}V1", j_call=j or f"{locus}J1")


def test_resolve_keeps_single_heavy():
    """Two beta contigs → only the top-ranked one survives."""
    df = _sc([
        _chain("c1", "b1", "TRB", "CASSHI", dup=100, umi=10),
        _chain("c1", "b2", "TRB", "CASSLO", dup=40, umi=4),
    ])
    out = sc.resolve_chains(df)
    assert out.height == 1
    assert out["sequence_id"].to_list() == ["b1"]


def test_resolve_keeps_two_light_above_threshold():
    """Dual-alpha kept when the second clears every secondary threshold."""
    df = _sc([
        _chain("c1", "b1", "TRB", "CASSB", dup=100, umi=10),
        _chain("c1", "a1", "TRA", "CAVA1", dup=80, umi=8),
        _chain("c1", "a2", "TRA", "CAVA2", dup=60, umi=6),   # 0.75 ratios, clears mins
    ])
    out = sc.resolve_chains(df)
    tra = out.filter(pl.col("locus") == "TRA")
    assert tra.height == 2
    assert set(tra["sequence_id"]) == {"a1", "a2"}


def test_resolve_drops_second_light_below_ratio():
    """Second alpha with too-low read ratio is dropped (keep top-1)."""
    df = _sc([
        _chain("c1", "a1", "TRA", "CAVA1", dup=100, umi=100),
        _chain("c1", "a2", "TRA", "CAVA2", dup=5, umi=2),    # dup ratio 0.05 < 0.1
    ])
    out = sc.resolve_chains(df)
    assert out.height == 1 and out["sequence_id"].to_list() == ["a1"]


def test_resolve_drops_second_light_below_min_umi():
    """Second alpha failing the absolute-UMI minimum is dropped even if ratios pass."""
    df = _sc([
        _chain("c1", "a1", "TRA", "CAVA1", dup=80, umi=8),
        _chain("c1", "a2", "TRA", "CAVA2", dup=40, umi=1),   # umi ratio ok but umi=1 < 2
    ])
    out = sc.resolve_chains(df)
    assert out["sequence_id"].to_list() == ["a1"]


def test_pair_chains_cartesian_two_alpha():
    """2 alpha + 1 beta → 2 pairs with _1/_2 suffixes; beta shared."""
    df = _sc([
        _chain("c1", "b1", "TRB", "CASSB", dup=100, umi=10),
        _chain("c1", "a1", "TRA", "CAVA1", dup=80, umi=8),
        _chain("c1", "a2", "TRA", "CAVA2", dup=60, umi=6),
    ])
    paired = sc.pair_chains(df)
    assert paired.height == 2
    assert paired["pair_id"].to_list() == ["c1_1", "c1_2"]
    assert set(paired["alpha_junction_aa"]) == {"CAVA1", "CAVA2"}
    assert set(paired["beta_junction_aa"]) == {"CASSB"}


def test_pair_chains_incomplete_not_emitted():
    """A cell missing one side is counted (multiplicity) but not emitted as a pair."""
    df = _sc([
        _chain("c1", "b1", "TRB", "CASSB", dup=100, umi=10),
        _chain("c1", "a1", "TRA", "CAVA1", dup=80, umi=8),
        _chain("c2", "b2", "TRB", "CASSONLY", dup=50, umi=5),   # beta only
    ])
    paired = sc.pair_chains(df)
    assert set(paired["cell_id"]) == {"c1"}
    assert paired.height == 1


def test_chain_multiplicity_quadrants():
    """Quadrant histogram reflects (n_light, n_heavy) per cell."""
    df = _sc([
        _chain("c1", "b1", "TRB", "CASSB", dup=100, umi=10),
        _chain("c1", "a1", "TRA", "CAVA1", dup=80, umi=8),   # (1 light, 1 heavy)
        _chain("c2", "b2", "TRB", "CASSC", dup=50, umi=5),   # (0 light, 1 heavy)
        _chain("c3", "a3", "TRA", "CAVA3", dup=30, umi=3),   # (1 light, 0 heavy)
    ])
    q = sc.chain_multiplicity(df)
    got = {(r["n_light"], r["n_heavy"]): r["cell_count"] for r in q.to_dicts()}
    assert got == {(1, 1): 1, (0, 1): 1, (1, 0): 1}


def _pair(cell, av, bv, adup=10, aumi=5):
    return dict(cell_id=cell, pair_id=cell,
                alpha_v_call="TRAV1", alpha_j_call="TRAJ1", alpha_junction_aa=av,
                alpha_umi_count=aumi, alpha_duplicate_count=adup,
                beta_v_call="TRBV1", beta_j_call="TRBJ1", beta_junction_aa=bv,
                beta_umi_count=5, beta_duplicate_count=10)


def test_flag_mispairing_noncanonical_alpha():
    """A beta with 3 distinct alpha → the 2 minority alphas flagged non-canonical."""
    paired = pl.from_dicts([
        _pair("c1", "A1", "B1"), _pair("c2", "A1", "B1"), _pair("c3", "A1", "B1"),
        _pair("c4", "A2", "B1"), _pair("c5", "A3", "B1"),
    ])
    res = sc.flag_mispairing(paired)
    assert res["mispairing_flag"].sum() == 2
    flagged = res.filter(pl.col("mispairing_flag"))
    assert set(flagged["alpha_junction_aa"]) == {"A2", "A3"}
    assert set(flagged["mispairing_reason"]) == {"noncanonical_alpha"}
    # canonical (majority) alpha is not flagged.
    assert res.filter(pl.col("alpha_junction_aa") == "A1")["mispairing_flag"].sum() == 0


def test_flag_mispairing_ambient_master():
    """A beta over the distinct-alpha ceiling is flagged ambient across all its rows."""
    paired = pl.from_dicts([
        _pair("c1", "A1", "B1"), _pair("c2", "A2", "B1"), _pair("c3", "A3", "B1"),
    ])
    res = sc.flag_mispairing(paired, max_slaves_per_master=2)  # B1 has 3 distinct alpha
    assert res["mispairing_flag"].all()
    assert set(res["mispairing_reason"]) == {"ambient_master"}


def test_flag_mispairing_drop_keeps_canonical():
    """drop=True removes flagged rows, keeping the canonical pairings only."""
    paired = pl.from_dicts([
        _pair("c1", "A1", "B1"), _pair("c2", "A1", "B1"), _pair("c3", "A1", "B1"),
        _pair("c4", "A2", "B1"),
    ])
    kept = sc.flag_mispairing(paired, drop=True)
    assert kept.height == 3
    assert set(kept["alpha_junction_aa"]) == {"A1"}
    assert "mispairing_flag" not in kept.columns


def test_resolve_drops_second_light_below_min_dup():
    """Second alpha passing ratio + min-UMI but failing only secondary_min_dup (absolute
    read floor) is dropped — isolates secondary_min_dup as the deciding condition."""
    df = _sc([
        _chain("c1", "a1", "TRA", "CAVA1", dup=20, umi=10),
        _chain("c1", "a2", "TRA", "CAVA2", dup=4, umi=3),   # 0.2/0.3 ratios ok, umi>=2, dup 4 < 5
    ])
    out = sc.resolve_chains(df)
    assert out.height == 1 and out["sequence_id"].to_list() == ["a1"]


def test_resolve_ranks_null_counts_last():
    """A null-count contig must not outrank a real high-count chain (nulls sort last)."""
    df = _sc([
        _chain("c1", "b_real", "TRB", "CASSREAL", dup=500, umi=10),
        _chain("c1", "b_null", "TRB", "CASSNULL", dup=None, umi=None),
    ])
    out = sc.resolve_chains(df)
    assert out.height == 1 and out["sequence_id"].to_list() == ["b_real"]


def test_resolve_and_pair_bcr_igk_igl():
    """B-cell: IGH is heavy, IGK+IGL are jointly light (dual light kept), and
    pair_chains routes IGK/IGL to alpha_* and IGH to beta_* for each BCR family."""
    df = _sc([
        _chain("b", "h", "IGH", "CARHEAVY", dup=100, umi=10),
        _chain("b", "k", "IGK", "CQKLIGHT", dup=80, umi=8),
        _chain("b", "l", "IGL", "CQLLIGHT", dup=60, umi=6),   # clears the joint 2nd-light rule
    ])
    resolved = sc.resolve_chains(df)
    assert set(resolved["locus"]) == {"IGH", "IGK", "IGL"} and resolved.height == 3
    igk = sc.pair_chains(df, locus_pair="IGH_IGK")
    assert igk.height == 1
    assert igk["alpha_junction_aa"][0] == "CQKLIGHT" and igk["beta_junction_aa"][0] == "CARHEAVY"
    igl = sc.pair_chains(df, locus_pair="IGH_IGL")
    assert igl["alpha_junction_aa"][0] == "CQLLIGHT" and igl["beta_junction_aa"][0] == "CARHEAVY"


def test_flag_mispairing_within_cell_dual_alpha_not_flagged():
    """A cell's legitimate second alpha (with the canonical alpha present in the SAME
    cell) is NOT flagged — within-cell dual-alpha is biology, not contamination."""
    dual = {**_pair("c1", "A2", "B1"), "pair_id": "c1_2"}    # c1 also carries canonical A1
    paired = pl.from_dicts([
        _pair("c1", "A1", "B1"), dual,                       # cell c1: both A1 and A2 vs B1
        _pair("c2", "A1", "B1"), _pair("c3", "A1", "B1"),    # A1 is B1's canonical slave
    ])
    res = sc.flag_mispairing(paired)
    assert res["mispairing_flag"].sum() == 0                 # nothing flagged
    # sanity: the SAME A2 in a cell lacking A1 IS flagged (cross-barcode smear).
    smear = pl.from_dicts([
        _pair("c1", "A1", "B1"), _pair("c2", "A1", "B1"), _pair("c4", "A2", "B1"),
    ])
    assert sc.flag_mispairing(smear).filter(pl.col("alpha_junction_aa") == "A2")[
        "mispairing_flag"].sum() == 1


def test_pair_chains_bad_locus_pair():
    df = _sc([_chain("c1", "b1", "TRB", "CASSB", dup=10, umi=5)])
    with pytest.raises(ValueError):
        sc.pair_chains(df, locus_pair="TRA_TRG")
