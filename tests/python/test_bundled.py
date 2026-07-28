"""Bundled precomputed models: the OLGA and EM-learned models ship in the package and load without
OLGA / HuggingFace at runtime, and are usable by the native Pgen."""
from __future__ import annotations

import polars as pl
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

    # V support must not have collapsed onto one gene: the shipped-broken TRB put 48.5% of all V
    # mass on TRBV19 with TRBV20-1 (the most-used human TRBV) at exactly 0. A gene-fraction bound
    # does NOT discriminate collapse from legitimate breadth -- the broken TRB was 21/59 genes
    # (36%), which overlaps a healthy but V-diverse BCR locus (IGH ~43%). So pin the two specific
    # signatures instead: (a) no single V gene dominates implausibly, and (b) the gene that was
    # famously zeroed is back. Both are locus-invariant except the named-gene check, which is
    # TRB-specific (the locus that actually shipped broken).
    vg = (m.tables["v_choice"].with_columns(pl.col("v_allele").str.split("*").list.first().alias("g"))
          .group_by("g").agg(pl.col("p").sum().alias("p")))
    top = vg.sort("p", descending=True)["p"][0]
    # TRD/TRG legitimately concentrate (TRDV1 ~0.56), so this only rules out a near-total collapse.
    assert top < 0.75, f"{locus}: top V gene holds {top:.0%} of mass (collapse onto one gene?)"
    if locus == "TRB":
        t201 = vg.filter(pl.col("g") == "TRBV20-1")["p"]
        assert t201.len() and t201[0] > 0, "TRBV20-1 is at P(V)=0 again (the shipped-model bug)"

    # Every allele that CAN be chosen must have a non-empty deletion distribution, or the
    # generative sampler draws it and cannot sample its trimming (the crash gene_prior once caused).
    dsum = m.tables["v_3_del"].group_by("v_allele").agg(pl.col("p").sum().alias("s"))
    orphan = (m.tables["v_choice"].filter(pl.col("p") > 0).join(dsum, on="v_allele", how="left")
                .filter(pl.col("s").fill_null(0.0) <= 0))
    assert orphan.height == 0, f"{locus}: {orphan.height} choosable V alleles have no deletion mass"

    if is_vdj:
        p2 = dict(zip(m.tables["n_d"]["n_d"].to_list(), m.tables["n_d"]["p"].to_list())).get(2, 0.0)
        # Lower bound too: an anchored D-D EM that regressed to a no-op would learn identically
        # zero tandem mass on every locus, and `< 0.1` alone would not catch it.
        assert 0.0 < p2 < 0.1, f"{locus}: P(n_D=2)={p2} (0 = D-D EM collapsed; >=0.1 = unregularized)"

    # Generate a real batch (not 5) so a rare choosable-but-unsamplable allele is actually hit.
    df = __import__("vdjtools.model.generate", fromlist=["generate"]).generate(m, 200, seed=0)
    assert native.pgen_nt(m, df["junction_nt"][0].upper(), df["v_call"][0], df["j_call"][0]) > 0.0


def test_sources_constant():
    assert SOURCES == ("olga", "learned", "arda")


# --- the arda-namespace set (organism-keyed; the only one with mouse) -----------------

def test_arda_set_is_reachable_and_covers_mouse():
    """The 9 shipped arda models load through the public loader, mouse included.

    They were on disk but unreachable before: ``SOURCES`` omitted ``"arda"``, and the loader
    keyed ``_bundled/<source>/<LOCUS>`` while the arda dirs are ``<organism>_<LOCUS>``.
    """
    keys = list_bundled()["arda"]
    assert set(keys) == {f"human_{loc}" for loc in LOCI} | {"mouse_TRA", "mouse_TRB"}


@pytest.mark.parametrize("locus", LOCI)
def test_load_arda_human_model(locus):
    m = load_bundled(locus, "arda")
    assert m.chain_type == ("VDJ" if locus in VDJ else "VJ")
    assert m.genomic["genes_v"].height > 0 and "v_choice" in m.tables


@pytest.mark.parametrize("locus", ["TRA", "TRB"])
def test_load_arda_mouse_model_and_generate(locus):
    """Mouse generation — impossible through the public API before, no mouse model was reachable."""
    m = load_bundled(locus, "arda", organism="mouse")
    assert m.chain_type == ("VDJ" if locus in VDJ else "VJ")
    df = __import__("vdjtools.model.generate", fromlist=["generate"]).generate(
        m, 20, seed=0, productive_only=True)
    assert df.height == 20
    assert df["v_call"].str.starts_with(f"{locus}V").all()   # arda IMGT names, mouse germline


def test_organism_is_not_silently_ignored_by_human_only_sets():
    # asking olga/learned for mouse must fail loudly, not hand back the human model
    for src in ("olga", "learned"):
        with pytest.raises(ValueError, match="human-only"):
            load_bundled("TRB", src, organism="mouse")
    with pytest.raises(FileNotFoundError, match="available"):
        load_bundled("IGH", "arda", organism="mouse")     # arda has no mouse IGH


def test_arda_and_learned_are_different_namespaces():
    """The point of the arda set: allele names come from arda, so they can differ from OLGA's."""
    a = load_bundled("TRB", "arda").tables["v_choice"]["v_allele"].to_list()
    lrn = load_bundled("TRB", "learned").tables["v_choice"]["v_allele"].to_list()
    assert set(a) != set(lrn)
