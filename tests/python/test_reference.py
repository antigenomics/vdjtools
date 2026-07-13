"""Phase 1a' — arda germline as the single source of truth + shared coordinate frame.

The load-bearing assertions: arda's ``anchor_nt`` is byte-identical to OLGA's ``anchor_index``
(so annotation/scenarios/stitching/Pgen share one frame), and ``cut_segment`` reproduces OLGA's
palindrome-extended segments. Reconciliation catalogs residual germline differences.
"""
from __future__ import annotations

import os
from pathlib import Path

import polars as pl
import pytest

from vdjtools.model import cut_segment, from_olga, load_germline, reconcile_olga, reverse_complement

OLGA_MODELS = Path(
    os.environ.get(
        "VDJTOOLS_OLGA_MODELS",
        str(Path(__file__).resolve().parent / "fixtures" / "olga" / "default_models"),
    )
)
pytest.importorskip("arda.cdr3fix", reason="arda (the [model] extra) not installed")
olm = pytest.importorskip("olga.load_model", reason="olga (the [oracle] extra) not installed")
pytestmark = pytest.mark.skipif(not OLGA_MODELS.exists(), reason=f"OLGA models not at {OLGA_MODELS}")

ALL_LOCI = ["TRA", "TRB", "TRG", "TRD", "IGH", "IGK", "IGL"]


def test_reverse_complement():
    assert reverse_complement("GTGT") == "ACAC"
    assert reverse_complement("CTAA") == "TTAG"


@pytest.mark.parametrize("locus", ALL_LOCI)
def test_load_germline_shape(locus):
    g = load_germline(locus, "human")
    assert set(g.columns) == {
        "allele", "gene", "segment", "sequence", "cdr3_anchor", "functionality", "functional", "status",
    }
    segs = set(g["segment"].unique())
    assert "V" in segs and "J" in segs
    if locus in ("TRB", "TRD", "IGH"):
        assert "D" in segs  # VDJ loci carry D germlines


def test_anchor_convention_shared_frame():
    """arda & OLGA share the anchor convention (0-based Cys104/[FW]118 offset into full germline).

    J anchors match exactly. A handful of V anchors differ — but only by whole framework codons
    (upstream IMGT-version drift), and only where the CDR3-region germline is *still identical*,
    so the CDR3 frame is unchanged. (Sources are never mixed: OLGA models use OLGA germline+anchor;
    arda-native models use arda germline+anchor.)
    """
    ref = load_germline("TRB", "human")
    arda = {(r[0], r[1]): (r[2], r[3]) for r in ref.select(["segment", "allele", "cdr3_anchor", "sequence"]).iter_rows()}
    g = olm.GenomicDataVDJ()
    d = OLGA_MODELS / "human_T_beta"
    g.load_igor_genomic_data(str(d / "model_params.txt"), str(d / "V_gene_CDR3_anchors.csv"), str(d / "J_gene_CDR3_anchors.csv"))
    olga_v = {g.genV[i][0]: g.genV[i][1] for i in range(len(g.genV))}
    olga_j = {g.genJ[i][0]: g.genJ[i][1] for i in range(len(g.genJ))}

    for seg, fname, olga_germ in [("V", "V_gene_CDR3_anchors.csv", olga_v), ("J", "J_gene_CDR3_anchors.csv", olga_j)]:
        for gene, oanchor, _func in pl.read_csv(d / fname).iter_rows():
            if (seg, gene) not in arda or not olga_germ.get(gene):
                continue
            aanchor, aseq = arda[(seg, gene)]
            if seg == "J":
                assert aanchor == oanchor, f"J anchor drift {gene}: arda {aanchor} vs olga {oanchor}"
            elif aseq == olga_germ[gene] and aanchor != oanchor:
                assert (aanchor - oanchor) % 3 == 0, f"V {gene} anchor off by non-codon: {aanchor} vs {oanchor}"


def test_cut_segment_reproduces_olga():
    """arda germline -> cut_segment equals OLGA's palindrome-extended cut segments."""
    g = olm.GenomicDataVDJ()
    d = OLGA_MODELS / "human_T_beta"
    g.load_igor_genomic_data(str(d / "model_params.txt"), str(d / "V_gene_CDR3_anchors.csv"), str(d / "J_gene_CDR3_anchors.csv"))
    olga_cutv = {g.genV[i][0]: g.cutV_genomic_CDR3_segs[i] for i in range(len(g.genV))}
    olga_cutj = {g.genJ[i][0]: g.cutJ_genomic_CDR3_segs[i] for i in range(len(g.genJ))}
    ref = load_germline("TRB", "human")
    seqs = {(r[0], r[1]): r[2] for r in ref.select(["segment", "allele", "sequence"]).iter_rows()}
    # V/J palindrome max is 4 in these models
    assert cut_segment(seqs[("V", "TRBV30*01")], "V", 4) == olga_cutv["TRBV30*01"]
    assert cut_segment(seqs[("J", "TRBJ1-2*01")], "J", 4) == olga_cutj["TRBJ1-2*01"]


def test_reconcile_olga_beta():
    """J germline is drop-in identical; the large majority of V alleles resolve in arda."""
    model = from_olga(OLGA_MODELS / "human_T_beta", locus="TRB")
    rep = reconcile_olga(model)
    # restrict to OLGA alleles that actually carry a CDR3 germline (functional in OLGA)
    used = rep.filter(pl.col("olga_len") > 0)
    j = used.filter(pl.col("segment") == "J")
    assert j["in_arda"].all() and j["germline_equal"].all(), "OLGA TRBJ germline must be arda-identical"
    v = used.filter(pl.col("segment") == "V")
    assert v["in_arda"].mean() >= 0.9, f"only {v['in_arda'].mean():.0%} of OLGA TRBV resolve in arda"


def test_load_germline_unknown_locus_raises():
    """load_germline raises a clear error on an unknown locus."""
    with pytest.raises(ValueError, match="no arda germline"):
        load_germline("XYZ", "human")
