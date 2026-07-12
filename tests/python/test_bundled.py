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

    D-bearing loci carry an **arda-anchored** tandem-D event with a biologically-plausible
    ``P(n_D=2)`` (a read may be tandem only where arda called a second D — this counters the
    identifiability that inflates unregularized D-D EM). Skips until the models are built and committed."""
    native = pytest.importorskip("vdjtools.model.native")
    if locus not in list_bundled().get("learned", []):
        pytest.skip(f"learned model for {locus} not built yet")
    m = load_bundled(locus, "learned")
    is_vdj = locus in VDJ
    assert m.chain_type == ("VDJ" if is_vdj else "VJ")
    assert ("d2_gene" in m.tables) == is_vdj  # D-loci carry the tandem-D event; VJ do not
    if is_vdj:
        p2 = dict(zip(m.tables["n_d"]["n_d"].to_list(), m.tables["n_d"]["p"].to_list())).get(2, 0.0)
        assert 0.0 <= p2 < 0.1  # anchored -> plausible (not the ~0.28 of unregularized EM)
    df = __import__("vdjtools.model.generate", fromlist=["generate"]).generate(m, 5, seed=0)
    assert native.pgen_nt(m, df["cdr3_nt"][0].upper(), df["v_call"][0], df["j_call"][0]) >= 0.0


def test_sources_constant():
    assert SOURCES == ("olga", "learned")
