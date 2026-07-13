"""Phase 1e — the ancestral generation sampler.

Checks: determinism, productive-only correctness, every generated sequence is scoreable by our
own Pgen (positive generation probability — required for EM), and gene usage tracks the model.
"""
from __future__ import annotations

import os
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

from vdjtools.model import from_olga
from vdjtools.model.generate import generate
from vdjtools.model.pgen import pgen_nt, prepare

OLGA_MODELS = Path(
    os.environ.get(
        "VDJTOOLS_OLGA_MODELS",
        str(Path(__file__).resolve().parent / "fixtures" / "olga" / "default_models"),
    )
)
pytest.importorskip("olga.load_model", reason="olga (the [oracle] extra) not installed")
pytestmark = pytest.mark.skipif(not OLGA_MODELS.exists(), reason=f"OLGA models not at {OLGA_MODELS}")


def _tra():
    return from_olga(OLGA_MODELS / "human_T_alpha", locus="TRA")


def test_deterministic():
    m = _tra()
    a = generate(m, 50, seed=1)
    assert a.equals(generate(m, 50, seed=1))
    assert not a.equals(generate(m, 50, seed=2))


def test_productive_only():
    df = generate(_tra(), 200, seed=4, productive_only=True)
    assert df["productive"].all()
    for nt, aa in zip(df["junction_nt"], df["junction_aa"]):
        assert len(nt) % 3 == 0 and "*" not in aa


def test_columns_and_calls():
    df = generate(_tra(), 20, seed=0)
    assert df.columns == ["junction_nt", "junction_aa", "v_call", "d_call", "d2_call", "j_call", "productive"]
    assert df["d_call"].is_null().all() and df["d2_call"].is_null().all()  # VJ locus has no D
    assert df["v_call"].str.starts_with("TRAV").all()


def test_generated_are_scoreable():
    """Every generated sequence must have positive Pgen under its own model (needed for EM)."""
    m = _tra()
    prep = prepare(m)
    df = generate(m, 150, seed=5)
    pg = np.array([pgen_nt(prep, r["junction_nt"], r["v_call"], r["j_call"]) for r in df.to_dicts()])
    assert (pg > 0).all(), f"{(pg == 0).sum()} generated sequences had Pgen 0"


def test_gene_usage_tracks_marginal():
    m = _tra()
    df = generate(m, 2000, seed=6)
    pv = dict(zip(m.tables["v_choice"]["v_allele"], m.tables["v_choice"]["p"]))
    counts = Counter(df["v_call"].to_list())
    total = sum(counts.values())
    genes = list(pv)
    emp = np.array([counts.get(g, 0) / total for g in genes])
    true = np.array([pv[g] for g in genes])
    assert np.corrcoef(emp, true)[0, 1] > 0.9


@pytest.mark.slow
def test_length_distribution_matches_olga():
    """CDR3-length distribution of our sampler matches OLGA's, on a VDJ locus."""
    import olga.generation_probability  # noqa: F401
    import olga.load_model as olm
    import olga.sequence_generation as osg
    from scipy import stats

    d = OLGA_MODELS / "human_T_beta"
    g = olm.GenomicDataVDJ()
    g.load_igor_genomic_data(str(d / "model_params.txt"), str(d / "V_gene_CDR3_anchors.csv"), str(d / "J_gene_CDR3_anchors.csv"))
    gm = olm.GenerativeModelVDJ()
    gm.load_and_process_igor_model(str(d / "model_marginals.txt"))
    sg = osg.SequenceGenerationVDJ(gm, g)
    np.random.seed(0)
    olga_len = [len(sg.gen_rnd_prod_CDR3()[0]) for _ in range(3000)]

    mine = generate(from_olga(d, locus="TRB"), 3000, seed=0, productive_only=True)
    mine_len = mine["junction_nt"].str.len_chars().to_list()
    # Kolmogorov–Smirnov: distributions should be statistically indistinguishable.
    assert stats.ks_2samp(olga_len, mine_len).pvalue > 0.01
