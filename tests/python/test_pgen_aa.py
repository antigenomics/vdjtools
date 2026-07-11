"""Phase 1f (aa) — amino-acid Pgen (codon-marginalizing DP), validated against OLGA.

VJ aa Pgen is fast; VDJ aa Pgen is the reference enumeration (correct but slow — the native
port is the performance follow-up). Both match OLGA exactly.
"""
from __future__ import annotations

import itertools
import os
from collections import defaultdict
from pathlib import Path

import numpy as np
import pytest

from vdjtools.model import from_olga
from vdjtools.model.pgen import pgen_aa, pgen_nt, prepare
from vdjtools.model.reference import _CODON_TABLE

OLGA_MODELS = Path(
    os.environ.get(
        "VDJTOOLS_OLGA_MODELS",
        "/Users/mikesh/vcs/code/mirpy/mir/resources/olga/default_models",
    )
)
pytest.importorskip("olga.load_model", reason="olga (the [oracle] extra) not installed")
pytestmark = pytest.mark.skipif(not OLGA_MODELS.exists(), reason=f"OLGA models not at {OLGA_MODELS}")

import olga.generation_probability as ogp  # noqa: E402
import olga.load_model as olm  # noqa: E402

TRA = OLGA_MODELS / "human_T_alpha"
TRB = OLGA_MODELS / "human_T_beta"


def _olga_vj(sub):
    g = olm.GenomicDataVJ()
    g.load_igor_genomic_data(str(sub / "model_params.txt"), str(sub / "V_gene_CDR3_anchors.csv"), str(sub / "J_gene_CDR3_anchors.csv"))
    m = olm.GenerativeModelVJ()
    m.load_and_process_igor_model(str(sub / "model_marginals.txt"))
    return ogp.GenerationProbabilityVJ(m, g)


def _olga_vdj(sub):
    g = olm.GenomicDataVDJ()
    g.load_igor_genomic_data(str(sub / "model_params.txt"), str(sub / "V_gene_CDR3_anchors.csv"), str(sub / "J_gene_CDR3_anchors.csv"))
    m = olm.GenerativeModelVDJ()
    m.load_and_process_igor_model(str(sub / "model_marginals.txt"))
    return ogp.GenerationProbabilityVDJ(m, g)


def test_vj_aa_matches_olga():
    prep = prepare(from_olga(TRA, locus="TRA"))
    po = _olga_vj(TRA)
    for aa, V, J in [("CAVRGCNAGGTSYGKLTF", "TRAV21*01", "TRAJ52*01"),
                     ("CAASQGSNDYKLSF", "TRAV13-1*01", "TRAJ20*01")]:
        assert np.isclose(pgen_aa(prep, aa, V, J), po.compute_aa_CDR3_pgen(aa, V, J), rtol=1e-6)


def test_vdj_aa_matches_olga_short():
    prep = prepare(from_olga(TRB, locus="TRB"))
    po = _olga_vdj(TRB)
    for aa, V, J in [("CASNRAGF", "TRBV7-9*01", "TRBJ1-2*01"),
                     ("CSAIRDGVF", "TRBV20-1*01", "TRBJ2-1*01")]:
        assert np.isclose(pgen_aa(prep, aa, V, J), po.compute_aa_CDR3_pgen(aa, V, J), rtol=1e-6)


def test_aa_equals_nt_sum():
    """aa Pgen == sum of nt Pgen over synonymous codons (self-consistency), on a short VJ CDR3."""
    from vdjtools.model.generate import generate

    m = from_olga(TRA, locus="TRA")
    prep = prepare(m)
    syn = defaultdict(list)
    for cod, a in _CODON_TABLE.items():
        syn[a].append(cod)
    # a short generated CDR3 with a bounded synonymous-codon count
    for r in sorted(generate(m, 600, seed=1, productive_only=True).to_dicts(), key=lambda r: len(r["cdr3_aa"])):
        aa, V, J = r["cdr3_aa"], r["v_call"], r["j_call"]
        if int(np.prod([len(syn[a]) for a in aa])) <= 40000:
            break
    brute = sum(pgen_nt(prep, "".join(c), V, J) for c in itertools.product(*[syn[a] for a in aa]))
    assert brute > 0
    assert np.isclose(pgen_aa(prep, aa, V, J), brute, rtol=1e-9)


@pytest.mark.slow
def test_vdj_aa_beta_oracle():
    """OLGA's published aa-Pgen docstring value (14-aa CDR3; the full reference enumeration)."""
    prep = prepare(from_olga(TRB, locus="TRB"))
    p = pgen_aa(prep, "CAWSVAPDRGGYTF", "TRBV30*01", "TRBJ1-2*01")
    assert np.isclose(p, 1.203646865765782e-10, rtol=1e-9)
