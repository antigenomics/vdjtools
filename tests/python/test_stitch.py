"""Phase 1c — contig stitching + the arda round-trip.

``stitch_contig`` rebuilds a full nt read from (V, J, CDR3); arda must then annotate it back to
the same junction and genes, proving synthetic and real reads share one alignment frame.
"""
from __future__ import annotations

import os
from pathlib import Path

import polars as pl
import pytest

from vdjtools.model import from_olga, stitch_contig, stitch_frame
from vdjtools.model.generate import generate

OLGA_MODELS = Path(
    os.environ.get(
        "VDJTOOLS_OLGA_MODELS",
        "/Users/mikesh/vcs/code/mirpy/mir/resources/olga/default_models",
    )
)
pytest.importorskip("olga.load_model", reason="olga (the [oracle] extra) not installed")
pytestmark = pytest.mark.skipif(not OLGA_MODELS.exists(), reason=f"OLGA models not at {OLGA_MODELS}")


def test_stitch_contig_none_paths():
    """stitch_contig returns None when a gene is unknown; stitch_frame yields a null contig."""
    m = from_olga(OLGA_MODELS / "human_T_beta", locus="TRB")
    good_v = m.genomic["genes_v"]["v_allele"][0]
    assert stitch_contig(m, "NOSUCHV", "TRBJ1-1*01", "TGTGCC") is None
    assert stitch_contig(m, good_v, "NOSUCHJ", "TGTGCC") is None
    frame = pl.DataFrame({"v_call": ["NOSUCHV"], "j_call": ["TRBJ1-1*01"], "cdr3_nt": ["TGTGCC"]})
    assert stitch_frame(m, frame)["contig"][0] is None


def test_stitch_contig_structure():
    """The contig embeds the CDR3 between the V and J framework germline (no arda needed)."""
    m = from_olga(OLGA_MODELS / "human_T_beta", locus="TRB")
    gen = generate(m, 5, seed=1, productive_only=True).to_dicts()
    r = gen[0]
    contig = stitch_contig(m, r["v_call"], r["j_call"], r["cdr3_nt"])
    assert contig is not None
    assert r["cdr3_nt"] in contig
    # CDR3 starts at the V anchor; upstream is pure V framework.
    vg = {x[0]: (x[1], x[2]) for x in m.genomic["genes_v"].select(["v_allele", "full_germline", "anchor"]).iter_rows()}
    fv, av = vg[r["v_call"]]
    assert contig.startswith(fv[:av])
    assert contig[av:av + len(r["cdr3_nt"])] == r["cdr3_nt"]


@pytest.mark.slow
def test_arda_roundtrip():
    """Stitch synthetic contigs, annotate with arda, and recover the junction + genes."""
    pytest.importorskip("arda", reason="arda (the [model] extra) not installed")
    from vdjtools.model.stitch import annotate

    m = from_olga(OLGA_MODELS / "human_T_beta", locus="TRB")
    gen = generate(m, 30, seed=2, productive_only=True)
    stitched = stitch_frame(m, gen).filter(pl.col("contig").is_not_null())
    calls = annotate(stitched["contig"].to_list())

    def gene(x):
        return x.split("*")[0] if x else None

    j_recovered = sum(a == b for a, b in zip(calls["junction"], stitched["cdr3_nt"]))
    v_gene_ok = sum(gene(a) == gene(b) for a, b in zip(calls["v_call"], stitched["v_call"]))
    j_gene_ok = sum(gene(a) == gene(b) for a, b in zip(calls["j_call"], stitched["j_call"]))
    n = stitched.height
    # arda must recover the exact junction and the V/J gene for the large majority.
    assert j_recovered / n > 0.9, f"junction recovered for only {j_recovered}/{n}"
    assert v_gene_ok / n > 0.85, f"V gene recovered for only {v_gene_ok}/{n}"
    assert j_gene_ok / n > 0.85, f"J gene recovered for only {j_gene_ok}/{n}"
