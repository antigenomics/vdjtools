"""Phase 1b — reference nucleotide Pgen, validated against the OLGA oracle.

Our direct scenario-sum ``pgen_nt`` (reading from our polars tables, no OLGA at runtime)
must reproduce OLGA's ``compute_nt_CDR3_pgen`` across all 7 loci, VDJ and VJ.
"""
from __future__ import annotations

import os
import random
from pathlib import Path

import numpy as np
import pytest

from vdjtools.model import from_olga
from vdjtools.model.pgen import pgen_nt, prepare

OLGA_MODELS = Path(
    os.environ.get(
        "VDJTOOLS_OLGA_MODELS",
        "/Users/mikesh/vcs/code/mirpy/mir/resources/olga/default_models",
    )
)
pytest.importorskip("olga.load_model", reason="olga (the [oracle] extra) not installed")
pytestmark = pytest.mark.skipif(
    not OLGA_MODELS.exists(), reason=f"OLGA bootstrap models not found at {OLGA_MODELS}"
)

import olga.generation_probability as ogp  # noqa: E402
import olga.load_model as olm  # noqa: E402
import olga.sequence_generation as osg  # noqa: E402

LOCI = [
    ("TRB", "human_T_beta", "VDJ"),
    ("TRD", "human_T_delta", "VDJ"),
    ("IGH", "human_B_heavy", "VDJ"),
    ("TRA", "human_T_alpha", "VJ"),
    ("TRG", "human_T_gamma", "VJ"),
    ("IGK", "human_B_kappa", "VJ"),
    ("IGL", "human_B_lambda", "VJ"),
]


def _olga_objects(sub: str, chain: str):
    d = OLGA_MODELS / sub
    args = (str(d / "model_params.txt"), str(d / "V_gene_CDR3_anchors.csv"), str(d / "J_gene_CDR3_anchors.csv"))
    if chain == "VDJ":
        g = olm.GenomicDataVDJ(); g.load_igor_genomic_data(*args)
        m = olm.GenerativeModelVDJ(); m.load_and_process_igor_model(str(d / "model_marginals.txt"))
        return g, m, ogp.GenerationProbabilityVDJ(m, g), osg.SequenceGenerationVDJ(m, g)
    g = olm.GenomicDataVJ(); g.load_igor_genomic_data(*args)
    m = olm.GenerativeModelVJ(); m.load_and_process_igor_model(str(d / "model_marginals.txt"))
    return g, m, ogp.GenerationProbabilityVJ(m, g), osg.SequenceGenerationVJ(m, g)


def test_beta_nt_known_oracle():
    """The value published in OLGA's own docstring (no OLGA call — fast)."""
    prep = prepare(from_olga(OLGA_MODELS / "human_T_beta", locus="TRB"))
    p = pgen_nt(prep, "TGTGCCTGGAGTGTAGCTCCGGACAGGGGTGGCTACACCTTC", "TRBV30*01", "TRBJ1-2*01")
    assert np.isclose(p, 2.3986503758867323e-12, rtol=1e-9)


@pytest.mark.parametrize("locus,sub,chain", [("TRB", "human_T_beta", "VDJ"), ("TRA", "human_T_alpha", "VJ")])
def test_nt_pgen_matches_olga_fast(locus, sub, chain):
    """Fast cross-check on one VDJ + one VJ locus (masked; a couple of sequences each)."""
    g, _m, pg_olga, gen = _olga_objects(sub, chain)
    prep = prepare(from_olga(OLGA_MODELS / sub, locus=locus))
    v_names = [x[0] for x in g.genV]
    j_names = [x[0] for x in g.genJ]
    random.seed(0)
    np.random.seed(0)
    for _ in range(2):
        nt, _aa, vi, ji = gen.gen_rnd_prod_CDR3()
        V, J = v_names[vi], j_names[ji]
        mine = pgen_nt(prep, nt, V, J)
        ref = pg_olga.compute_nt_CDR3_pgen(nt, V, J)
        assert ref > 0 and np.isclose(mine, ref, rtol=1e-6), f"{locus} {nt} {V} {J}: {mine} vs {ref}"


@pytest.mark.slow
@pytest.mark.parametrize("locus,sub,chain", LOCI, ids=[x[0] for x in LOCI])
def test_nt_pgen_matches_olga(locus, sub, chain):
    g, _m, pg_olga, gen = _olga_objects(sub, chain)
    prep = prepare(from_olga(OLGA_MODELS / sub, locus=locus))
    v_names = [x[0] for x in g.genV]
    j_names = [x[0] for x in g.genJ]

    random.seed(0)
    np.random.seed(0)
    checked = 0
    for _ in range(5):
        nt, _aa, vi, ji = gen.gen_rnd_prod_CDR3()
        V, J = v_names[vi], j_names[ji]
        mine = pgen_nt(prep, nt, V, J)                 # gene-masked (fast, exercises the D/ins DP)
        ref = pg_olga.compute_nt_CDR3_pgen(nt, V, J)
        assert ref > 0 and np.isclose(mine, ref, rtol=1e-6), f"{locus} {nt} {V} {J}: {mine} vs {ref}"
        checked += 1
    assert checked == 5


@pytest.mark.slow
def test_no_mask_sum_matches_olga():
    """The unrestricted sum-over-all-genes path (no usage mask), on a VDJ and a VJ locus."""
    for sub, locus, chain in [("human_T_beta", "TRB", "VDJ"), ("human_T_alpha", "TRA", "VJ")]:
        _g, _m, pg_olga, gen = _olga_objects(sub, chain)
        prep = prepare(from_olga(OLGA_MODELS / sub, locus=locus))
        random.seed(1)
        np.random.seed(1)
        for _ in range(3):
            nt = gen.gen_rnd_prod_CDR3()[0]
            mine = pgen_nt(prep, nt)               # no V/J restriction
            ref = pg_olga.compute_nt_CDR3_pgen(nt)
            assert np.isclose(mine, ref, rtol=1e-6), f"{locus} no-mask {nt}: {mine} vs {ref}"
