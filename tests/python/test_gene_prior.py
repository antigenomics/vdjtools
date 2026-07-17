"""Dirichlet gene prior — keeps EM from permanently deleting a real germline gene.

P(V)=0 is an ABSORBING STATE of this EM: the E-step weights every scenario by P(V), so an allele
that ever reaches zero count can never be re-attributed. One unlucky iteration kills a real gene
for good — which is how the shipped TRB model came to carry 21 of 89 V alleles, with TRBV20-1
(the most-used human TRBV) at exactly 0.
"""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from vdjtools.model import load_bundled
from vdjtools.model.generate import generate
from vdjtools.model.infer import _functional_support, infer_native


@pytest.fixture(scope="module")
def trb():
    return load_bundled("TRB", "olga")


def _v_mass(model):
    t = model.tables["v_choice"]
    return {r["v_allele"]: r["p"] for r in t.iter_rows(named=True)}


def test_gene_prior_zero_is_byte_identical_to_mle(trb):
    """The default must not perturb anything — the exact-Pgen invariant depends on it."""
    draws = generate(trb, 300, seed=0)
    seqs = [s.upper() for s in draws["junction_nt"].to_list()]
    a, _ = infer_native(trb, seqs, max_iter=2, tol=0.0, single_d=True, gene_prior=0.0)
    b, _ = infer_native(trb, seqs, max_iter=2, tol=0.0, single_d=True)      # default
    for ev in a.tables:
        assert a.tables[ev].equals(b.tables[ev]), f"{ev} differs — gene_prior=0.0 is not the default"


def test_gene_prior_keeps_unseen_functional_alleles_reachable(trb):
    """A functional allele absent from the training data must keep non-zero mass.

    Without the prior it is dead permanently; with it, the germline (not one EM run) decides what
    exists, while the data still decides usage.
    """
    draws = generate(trb, 400, seed=1)
    seqs = [s.upper() for s in draws["junction_nt"].to_list()]
    zero, _ = infer_native(trb, seqs, max_iter=3, tol=0.0, single_d=True, gene_prior=0.0)
    prior, _ = infer_native(trb, seqs, max_iter=3, tol=0.0, single_d=True, gene_prior=1.0)

    functional = _functional_support(trb, "v")
    n_zero = sum(1 for a, p in _v_mass(zero).items() if p <= 0 and a in functional)
    n_prior = sum(1 for a, p in _v_mass(prior).items() if p <= 0 and a in functional)
    assert n_prior < n_zero, f"prior did not rescue anything ({n_prior} vs {n_zero} dead alleles)"
    assert n_prior == 0, f"{n_prior} functional V alleles still at p=0 under the prior"


def test_gene_prior_gives_no_mass_to_nonfunctional_alleles(trb):
    """Pseudogenes/ORFs are deliberately excluded — the model cannot score them anyway."""
    draws = generate(trb, 200, seed=2)
    seqs = [s.upper() for s in draws["junction_nt"].to_list()]
    m, _ = infer_native(trb, seqs, max_iter=2, tol=0.0, single_d=True, gene_prior=1.0)
    functional = _functional_support(trb, "v")
    leaked = [a for a, p in _v_mass(m).items() if p > 0 and a not in functional]
    assert not leaked, f"prior leaked mass onto non-functional alleles: {leaked[:3]}"


def test_gene_prior_still_normalizes(trb):
    draws = generate(trb, 200, seed=3)
    seqs = [s.upper() for s in draws["junction_nt"].to_list()]
    m, _ = infer_native(trb, seqs, max_iter=2, tol=0.0, single_d=True, gene_prior=1.0)
    assert m.tables["v_choice"]["p"].sum() == pytest.approx(1.0)
    m.validate()
