"""E-step gene masks must widen on ambiguity, never narrow.

An AIRR call like ``IGHV3-23*01,IGHV3-23D*01`` is an aligner saying "I cannot tell these apart".
Keeping only the first gene makes the true scenario unreachable, so EM misattributes the read.
IGH is where this bites: 14.5% of its non-functional clonotypes carry an ambiguous V call.
"""
from __future__ import annotations

import pytest

from vdjtools.model import load_bundled
from vdjtools.model.infer import _gene_to_alleles, call_alleles, gene_masks


@pytest.fixture(scope="module")
def trb():
    return load_bundled("TRB", "olga")


def test_allele_call_expands_to_every_allele_of_its_gene(trb):
    """A call of *03 where the truth is *01 must not exclude *01.

    The gene is looked up rather than hardcoded: OLGA's TRB index carries a single allele for
    most genes (TRBV20-1 has only *01), so a hardcoded pick tests nothing.
    """
    va = _gene_to_alleles(trb, "v")
    gene = next(g for g, alls in va.items() if len(alls) > 1)
    got = call_alleles(va, f"{gene}*99")                    # an allele the model has never seen
    assert set(got) == set(va[gene])
    assert len(got) > 1, f"{gene} has several model alleles; the mask must keep them all"


def test_comma_ambiguity_unions_both_genes_rather_than_dropping_one(trb):
    """The bug: split('*')[0] on 'A*01,B*01' yields 'A' and silently loses B."""
    va = _gene_to_alleles(trb, "v")
    a, b = "TRBV6-2", "TRBV6-3"
    assert a in va and b in va
    got = call_alleles(va, f"{a}*01,{b}*01")
    assert set(got) == set(va[a]) | set(va[b])
    assert set(va[b]) <= set(got), "the second gene of an ambiguous call was dropped"


def test_call_alleles_is_deduplicated_and_order_stable(trb):
    va = _gene_to_alleles(trb, "v")
    got = call_alleles(va, "TRBV20-1*01,TRBV20-1*03")      # same gene twice
    assert got == list(dict.fromkeys(got)), "duplicate alleles leaked into the mask"
    assert set(got) == set(va["TRBV20-1"])


def test_unknown_or_empty_call_yields_unrestricted(trb):
    va = _gene_to_alleles(trb, "v")
    assert call_alleles(va, None) == []
    assert call_alleles(va, "") == []
    assert call_alleles(va, "NOTAGENE*01") == []           # [] = unrestricted, the honest default


def test_gene_masks_threads_ambiguity_through(trb):
    masks = gene_masks(trb, ["TRBV6-2*01,TRBV6-3*01"], ["TRBJ2-7*01"])
    v_mask, j_mask, d_mask = masks[0]
    va = _gene_to_alleles(trb, "v")
    assert set(v_mask) == set(va["TRBV6-2"]) | set(va["TRBV6-3"])
    assert d_mask is None                                   # D deliberately unrestricted here
    assert j_mask


def test_whitespace_after_comma_is_tolerated(trb):
    va = _gene_to_alleles(trb, "v")
    assert call_alleles(va, "TRBV6-2*01, TRBV6-3*01") == call_alleles(va, "TRBV6-2*01,TRBV6-3*01")
