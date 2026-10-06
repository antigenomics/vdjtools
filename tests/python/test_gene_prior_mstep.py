"""Exact prior edge cases without repeated full-locus E-steps."""
from __future__ import annotations

import inspect
from types import SimpleNamespace

import polars as pl
import pytest

from vdjtools._core import make_counts
from vdjtools.model.infer import _mstep_native, infer_native
from vdjtools.model.io import from_germline
from vdjtools.model.native import pack


@pytest.fixture
def mstep_case(toy_germline_vdj):
    # Each segment has an observed functional allele, an observed nonfunctional allele,
    # and an unseen functional allele. Native Counts supplies the other array dimensions.
    extra_j = toy_germline_vdj.filter(pl.col("allele") == "TOYJ1*01").with_columns(
        allele=pl.lit("TOYJ3*01"), gene=pl.lit("TOYJ3"))
    germline = pl.concat([toy_germline_vdj, extra_j]).with_columns(
        functional=~pl.col("allele").is_in(["TOYV1*02", "TOYJ2*01"]))
    model = from_germline(germline, locus="TOY", ins_max=3)
    pm, _, _ = pack(model)
    raw = make_counts(pm)
    fields = ("v_choice", "j_choice", "v_3_del", "j_5_del", "d_gene", "d_del",
              "ins_vd", "ins_dj", "dinucl_vd", "dinucl_dj", "n_d")
    counts = SimpleNamespace(**{name: getattr(raw, name) for name in fields})
    alleles = [model.genomic[f"genes_{seg}"][f"{seg}_allele"].to_list()
               for seg in ("v", "j", "d")]
    v_counts = {"TOYV1*01": 2.0, "TOYV1*02": 1.0, "TOYV2*01": 0.0}
    j_counts = {"TOYJ1*01": 2.0, "TOYJ2*01": 1.0, "TOYJ3*01": 0.0}
    counts.v_choice = [v_counts[a] for a in alleles[0]]
    counts.j_choice = [j_counts[a] for a in alleles[1]]
    nbins = {"v": pm.nbins_v, "j": pm.nbins_j, "d5": pm.nbins_d5, "d3": pm.nbins_d3}
    return model, counts, *alleles, nbins


@pytest.mark.parametrize("prior", [0.0, 1.0, 10.0])
def test_prior_only_adds_counts_to_observed_functional_alleles(mstep_case, prior):
    tables = _mstep_native(*mstep_case, gene_prior=prior)
    for seg, observed, nonfunctional, unseen in (
        ("v", "TOYV1*01", "TOYV1*02", "TOYV2*01"),
        ("j", "TOYJ1*01", "TOYJ2*01", "TOYJ3*01"),
    ):
        masses = dict(tables[f"{seg}_choice"].iter_rows())
        assert masses[observed] == pytest.approx((2.0 + prior) / (3.0 + prior))
        assert masses[nonfunctional] == pytest.approx(1.0 / (3.0 + prior))
        assert masses[unseen] == 0.0
        assert sum(masses.values()) == pytest.approx(1.0)


def test_zero_counts_stay_zero_under_a_prior(mstep_case):
    counts = mstep_case[1]
    counts.v_choice = [0.0] * len(counts.v_choice)
    counts.j_choice = [0.0] * len(counts.j_choice)
    tables = _mstep_native(*mstep_case, gene_prior=1.0)
    for seg in ("v", "j"):
        assert tables[f"{seg}_choice"]["p"].to_list() == [0.0] * len(getattr(counts, f"{seg}_choice"))


def test_default_prior_is_zero_and_mstep_is_byte_identical(mstep_case):
    assert inspect.signature(infer_native).parameters["gene_prior"].default == 0.0
    implicit = _mstep_native(*mstep_case)
    explicit = _mstep_native(*mstep_case, gene_prior=0.0)
    assert implicit.keys() == explicit.keys()
    for name in implicit:
        assert implicit[name].equals(explicit[name]), name
