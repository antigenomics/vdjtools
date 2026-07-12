"""Metaclonotype grouping (1mm) — needs vdjmatch ([overlap] extra); skips cleanly without."""
from __future__ import annotations

import polars as pl
import pytest

pytest.importorskip("vdjmatch")  # metaclonotypes delegates to vdjmatch.cluster.overlap

from vdjtools.biomarker import fisher_association, metaclonotypes  # noqa: E402


def test_hamming1_groups_respect_v_and_j():
    """One substitution + same V/J → one metaclonotype; 2 subs or different V → separate."""
    df = pl.DataFrame({
        "junction_aa": ["CASSKF", "CASSRF", "CAXXKF", "CASSKF"],
        "v_call":  ["TRBV1",  "TRBV1",  "TRBV1",  "TRBV2"],
        "j_call":  ["TRBJ1",  "TRBJ1",  "TRBJ1",  "TRBJ1"],
    })
    m = metaclonotypes(df, scope="1,0,0,1")
    g = {(r["junction_aa"], r["v_call"]): r["meta_id"] for r in m.iter_rows(named=True)}
    assert g[("CASSKF", "TRBV1")] == g[("CASSRF", "TRBV1")]   # 1 sub, same V/J → grouped
    assert g[("CAXXKF", "TRBV1")] != g[("CASSKF", "TRBV1")]   # 2 subs → separate
    assert g[("CASSKF", "TRBV2")] != g[("CASSKF", "TRBV1")]   # different V → separate
    assert m["meta_id"].n_unique() == 3


def test_cdr3_only_grouping_ignores_v():
    """match_v/match_j False → group on CDR3 alone (V/J columns not even required)."""
    df = pl.DataFrame({"junction_aa": ["CASSKF", "CASSRF", "CAXXKF"]})
    m = metaclonotypes(df, scope="1,0,0,1", match_v=False, match_j=False)
    assert m["meta_id"].n_unique() == 2  # KF~RF merge; XXKF separate


def test_1mm_merges_incidence_vs_exact():
    """1mm collapses two 1-sub-neighbour CDR3s into one feature present in all subjects."""
    df = pl.DataFrame({
        "sample_id": ["s1", "s2", "s3", "s4"],
        "v_call": ["TRBV1"] * 4, "j_call": ["TRBJ1"] * 4,
        "junction_aa": ["CASSKF", "CASSRF", "CASSKF", "CASSRF"],  # KF ~ RF (1 sub)
        "duplicate_count": [1, 1, 1, 1],
    })
    pheno = pl.DataFrame({"sample_id": ["s1", "s2", "s3", "s4"], "y": [True, True, False, False]})

    exact = fisher_association(df, pheno, pheno_col="y", min_incidence=2)
    assert sorted(exact["incidence"].to_list()) == [2, 2]      # KF and RF separate

    mm = fisher_association(df, pheno, pheno_col="y", match="1mm", min_incidence=2)
    top = mm.row(0, named=True)
    assert top["incidence"] == 4 and top["n_members"] == 2     # merged, present in all 4
    assert "meta_id" in mm.columns and top["junction_aa"] in {"CASSKF", "CASSRF"}

    # Empty 1mm result must share the non-empty column order (concat-safe, no ShapeError).
    empty = fisher_association(df, pheno, pheno_col="y", match="1mm", min_incidence=999)
    assert empty.height == 0 and empty.columns == mm.columns
    pl.concat([mm, empty])
