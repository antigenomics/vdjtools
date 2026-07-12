"""Bundled precomputed models: the OLGA and EM-learned models ship in the package and load without
OLGA / HuggingFace at runtime, and are usable by the native Pgen."""
from __future__ import annotations

import pytest

from vdjtools.model import list_bundled, load_bundled
from vdjtools.model.bundled import LOCI, SOURCES

VDJ = {"TRB", "TRD", "IGH"}


def test_all_olga_loci_bundled():
    """Every one of the 7 human loci ships an OLGA bootstrap model."""
    assert set(list_bundled()["olga"]) == set(LOCI)


@pytest.mark.parametrize("locus", LOCI)
def test_load_olga_model(locus):
    m = load_bundled(locus, "olga")
    assert m.chain_type == ("VDJ" if locus in VDJ else "VJ")
    assert m.genomic["genes_v"].height > 0 and "v_choice" in m.tables


def test_native_pgen_on_bundled():
    """A bundled model is directly usable by the native Pgen (no external model files)."""
    native = pytest.importorskip("vdjtools.model.native")
    m = load_bundled("TRB", "olga")
    assert native.pgen_aa(m, "CASSLGF") >= 0.0


def test_bad_source_and_locus_raise():
    with pytest.raises(ValueError):
        load_bundled("TRB", "nonsense")
    with pytest.raises(FileNotFoundError):
        load_bundled("ZZZ", "olga")


@pytest.mark.parametrize("locus", LOCI)
def test_learned_models_load_and_score(locus):
    """The EM-learned models (fit to real out-of-frame reads) load and are usable by the native Pgen.

    Shipped single-D: unregularized D-D EM over-attributes tandems on real data (identifiability), so
    the learned models capture real-repertoire gene-usage/trim/insertion marginals only. Skips cleanly
    until the learned models are built and committed."""
    native = pytest.importorskip("vdjtools.model.native")
    if locus not in list_bundled().get("learned", []):
        pytest.skip(f"learned model for {locus} not built yet")
    m = load_bundled(locus, "learned")
    assert m.chain_type == ("VDJ" if locus in VDJ else "VJ")
    assert "d2_gene" not in m.tables  # shipped single-D
    df = __import__("vdjtools.model.generate", fromlist=["generate"]).generate(m, 5, seed=0)
    assert native.pgen_nt(m, df["cdr3_nt"][0].upper(), df["v_call"][0], df["j_call"][0]) >= 0.0


def test_sources_constant():
    assert SOURCES == ("olga", "learned")
