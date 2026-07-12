"""Tests for vdjtools.features.physchem — property-table load and region profiles."""
import math

import polars as pl

from vdjtools.io import schema as S
from vdjtools import features as F


def _frame(cdr3, counts=None, v=None, j=None):
    n = len(cdr3)
    df = pl.DataFrame({
        S.V_CALL: v or ["TRBV1"] * n, S.J_CALL: j or ["TRBJ1"] * n,
        S.JUNCTION_AA: cdr3, S.COUNT: counts or [1] * n,
    })
    return S.add_locus(S.normalize(df, recompute_freq=True))


def test_property_table_loads():
    tbl = F.load_property_table()
    assert tbl.height == 20                                  # 20 amino acids
    assert "hydropathy" in tbl.columns and "kf10" in tbl.columns
    assert tbl.filter(pl.col("amino_acid") == "A")["hydropathy"][0] == 1.8


def test_physchem_all_region_hand_value():
    # "AC": hydropathy (1.8+2.5)/2 = 2.15, volume (67+86)/2 = 76.5
    prof = F.physchem_profile(_frame(["AC"]), group_by="locus", region="all",
                              weight="unique", properties=("hydropathy", "volume"))
    got = dict(zip(prof["property"].to_list(), prof["mean_value"].to_list()))
    assert math.isclose(got["hydropathy"], 2.15, rel_tol=1e-9)
    assert math.isclose(got["volume"], 76.5, rel_tol=1e-9)


def test_physchem_regions_on_CASSLGF():
    # CASSLGF (L=7): all=mean over 7; trimmed=[3:-3]="S"; center=[1:6]="ASSLG"
    df = _frame(["CASSLGF"])
    hp = {"C": 2.5, "A": 1.8, "S": -0.8, "L": 3.8, "G": -0.4, "F": 2.8}
    exp_all = sum(hp[a] for a in "CASSLGF") / 7
    exp_center = sum(hp[a] for a in "ASSLG") / 5
    for region, expected in (("all", exp_all), ("trimmed", -0.8), ("center", exp_center)):
        prof = F.physchem_profile(df, group_by="locus", region=region,
                                  weight="unique", properties=("hydropathy",))
        assert math.isclose(prof["mean_value"][0], expected, rel_tol=1e-9)


def test_physchem_trimmed_skips_short_cdr3():
    # both len <= 6 -> empty trimmed region -> empty result
    prof = F.physchem_profile(_frame(["CASS", "CASSL"]), group_by="locus",
                              region="trimmed", weight="unique",
                              properties=("hydropathy",))
    assert prof.height == 0


def test_physchem_weighted_by_reads():
    # two clonotypes, group by locus, region=all, hydropathy; weighted by reads
    # "A"(1.8) weight 3, "C"(2.5) weight 1 -> (1.8*3 + 2.5*1)/4 = 1.975
    prof = F.physchem_profile(_frame(["A", "C"], counts=[3, 1]), group_by="locus",
                              region="all", weight="reads", properties=("hydropathy",))
    assert math.isclose(prof["mean_value"][0], (1.8 * 3 + 2.5 * 1) / 4, rel_tol=1e-9)
