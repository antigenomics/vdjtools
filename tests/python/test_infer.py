"""Phase 1d — EM inference, closed-loop recovery of a known model's marginals.

Generate synthetic sequences from an OLGA model, then ``infer`` must recover that model's
marginals. Recovery is asserted on the *estimable* quantities: global parameters (insertion
length, dinucleotide Markov), marginal J usage, the aggregate deletion profile, and
germline-group V usage. (Per-allele V usage is limited by germline ambiguity — alleles with
identical CDR3-region germline are indistinguishable from the CDR3 alone, an intrinsic limit
shared by OLGA/IGoR — so we validate at the germline-group level.)
"""
from __future__ import annotations

import collections
import os
from pathlib import Path

import numpy as np
import pytest

from vdjtools.model import from_olga
from vdjtools.model.generate import generate
from vdjtools.model.infer import infer

OLGA_MODELS = Path(
    os.environ.get(
        "VDJTOOLS_OLGA_MODELS",
        str(Path(__file__).resolve().parent / "fixtures" / "olga" / "default_models"),
    )
)
pytest.importorskip("olga.load_model", reason="olga (the [oracle] extra) not installed")
pytestmark = pytest.mark.skipif(not OLGA_MODELS.exists(), reason=f"OLGA models not at {OLGA_MODELS}")

TRA = OLGA_MODELS / "human_T_alpha"
TRB = OLGA_MODELS / "human_T_beta"


def _corr(a: dict, b: dict) -> float:
    keys = set(a) | set(b)
    return float(np.corrcoef([a.get(k, 0.0) for k in keys], [b.get(k, 0.0) for k in keys])[0, 1])


def _table_dict(model, event, keycols):
    return {tuple(r[:-1]): r[-1] for r in model.tables[event].select([*keycols, "p"]).iter_rows()}


def _j_marginal(model):
    pv = dict(zip(model.tables["v_choice"]["v_allele"], model.tables["v_choice"]["p"]))
    d = collections.defaultdict(float)
    for v, j, p in model.tables["j_choice"].select(["v_allele", "j_allele", "p"]).iter_rows():
        d[j] += pv.get(v, 0.0) * p
    return d


def _agg_delv(model):
    pv = dict(zip(model.tables["v_choice"]["v_allele"], model.tables["v_choice"]["p"]))
    d = collections.defaultdict(float)
    for v, n, p in model.tables["v_3_del"].select(["v_allele", "ndel", "p"]).iter_rows():
        d[n] += pv.get(v, 0.0) * p
    return d


def _grouped_vusage(model, group):
    d = collections.defaultdict(float)
    for v, p in zip(model.tables["v_choice"]["v_allele"], model.tables["v_choice"]["p"]):
        d[group[v]] += p
    return d


def test_infer_plumbing():
    """EM runs, returns a valid model, and populates the report."""
    m = from_olga(TRA, locus="TRA")
    seqs = generate(m, 120, seed=0)["junction_nt"].to_list()
    fit, rep = infer(m, seqs, max_iter=2)
    fit.validate()
    assert rep.n_iter == 2
    assert len(rep.loglik) == 2 and len(rep.gene_tv) == 2
    assert set(fit.tables) == set(m.tables)


def test_gene_masks():
    """gene_masks expands an allele call to all model alleles of that gene (handles ambiguity)."""
    from vdjtools.model.infer import gene_masks

    m = from_olga(TRB, locus="TRB")
    masks = gene_masks(m, ["TRBV20-1*03", "TRBV7-9*01"], ["TRBJ1-1*01", "TRBJ2-1*01"])
    assert all(a.startswith("TRBV20-1*") for a in masks[0][0])
    assert "TRBV20-1*01" in masks[0][0]  # true allele included even though *03 was called
    assert all(a.startswith("TRBJ1-1*") for a in masks[0][1])
    assert masks[0][2] is None  # D is left unrestricted


@pytest.mark.slow
def test_vdj_estep_accumulates():
    """Cover the VDJ EM soft-count path (_accum_vdj: D-gene / delD / VD+DJ insertion counts).

    A single E-step on one (short) TRB sequence — the full VDJ EM is validated for correctness
    but is impractically slow in pure Python (~tens of s/seq; the D×deletion×position enumeration
    is the prime C++/arda-masking target), so this just exercises the path and checks it harvests
    the D-junction event counts.
    """
    from collections import defaultdict

    from vdjtools.model.infer import _estep_seq
    from vdjtools.model.pgen import prepare

    m = from_olga(TRB, locus="TRB")
    prep = prepare(m)
    seq = min(generate(m, 60, seed=1, productive_only=True)["junction_nt"].to_list(), key=len)
    # Every VDJ scenario also harvests the D-count event, so ``n_d`` must have a bucket too
    # (the E-step emits ``("n_d", (1,))`` per single-D scenario).
    counts = {name: defaultdict(float) for name in m.tables}
    pg = _estep_seq(prep, seq, counts)
    assert pg > 0
    for ev in ("d_gene", "d_del", "vd_ins", "dj_ins", "vd_dinucl", "dj_dinucl", "n_d"):
        assert counts[ev], f"{ev} soft counts were not accumulated"


@pytest.mark.slow
def test_recovers_tra_marginals():
    m = from_olga(TRA, locus="TRA")
    seqs = generate(m, 1500, seed=11)["junction_nt"].to_list()
    fit, rep = infer(m, seqs, max_iter=8)

    # EM behaves: log-likelihood climbs after the first update; gene usage stabilizes.
    assert rep.loglik[-1] >= rep.loglik[1]
    assert rep.gene_tv[-1] < rep.gene_tv[1]

    # Global parameters recover tightly.
    assert _corr(_table_dict(fit, "vj_ins", ["length"]), _table_dict(m, "vj_ins", ["length"])) > 0.95
    assert _corr(_table_dict(fit, "vj_dinucl", ["from_nt", "to_nt"]),
                 _table_dict(m, "vj_dinucl", ["from_nt", "to_nt"])) > 0.98
    # Marginal J usage and the aggregate V-deletion profile.
    assert _corr(_j_marginal(fit), _j_marginal(m)) > 0.90
    assert _corr(_agg_delv(fit), _agg_delv(m)) > 0.90
    # V usage at the germline-group level (alleles with identical germline are indistinguishable).
    group = {v: s for v, s in zip(m.genomic["genes_v"]["v_allele"], m.genomic["genes_v"]["cut_segment"])}
    assert _corr(_grouped_vusage(fit, group), _grouped_vusage(m, group)) > 0.90
