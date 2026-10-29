"""Bundled precomputed models: the OLGA and EM-learned models ship in the package and load without
OLGA / HuggingFace at runtime, and are usable by the native Pgen."""
from __future__ import annotations

import polars as pl
import pytest

from vdjtools.model import list_bundled, load_bundled
from vdjtools.model.bundled import LOCI, SOURCES

VDJ = {"TRB", "TRD", "IGH"}
MOUSE_LOCI = ('TRA', 'TRB', 'IGH', 'IGK', 'IGL')


def test_canonical_learned_catalog_covers_human_and_mouse():
    assert set(list_bundled()['isalgo']) == set(LOCI) | {f'mouse_{locus}' for locus in MOUSE_LOCI}


@pytest.mark.parametrize('organism,locus', [('human', locus) for locus in LOCI] +
                         [('mouse', locus) for locus in MOUSE_LOCI])
def test_default_resolves_canonical_learned_model(organism, locus):
    default = load_bundled(locus, organism=organism)
    explicit = load_bundled(locus, 'learned', organism=organism)
    assert default.manifest == explicit.manifest
    assert default.tables.keys() == explicit.tables.keys()
    assert all(default.tables[name].equals(explicit.tables[name]) for name in default.tables)


@pytest.mark.parametrize('locus', ['TRG', 'TRD'])
def test_unbundled_mouse_loci_fail_without_human_fallback(locus):
    # Human models exist, but cannot stand in for the deferred mouse contexts.
    with pytest.raises(FileNotFoundError, match="organism 'mouse'"):
        load_bundled(locus, organism='mouse')


def test_all_olga_loci_bundled():
    """Every one of the 7 human loci ships an OLGA bootstrap model."""
    assert set(list_bundled()["olga-legacy"]) == set(LOCI)


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

    D-bearing loci carry the joint single/tandem-D event. Preliminary ARDA D calls do not
    constrain its histories."""
    native = pytest.importorskip("vdjtools.model.native")
    assert locus in list_bundled()['isalgo']
    m = load_bundled(locus, "learned")
    is_vdj = locus in VDJ
    assert m.chain_type == ("VDJ" if is_vdj else "VJ")
    assert ("d_history" in m.tables) == is_vdj  # fresh D-bearing fits use joint histories

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
        from vdjtools.model.dd import d_count_probabilities
        probabilities = d_count_probabilities(m)
        assert sum(probabilities.values()) == pytest.approx(1.0)
        p2 = probabilities.get(2, 0.0)
        # Both histories must remain supported. A universal 10% ceiling came from
        # greedy-D2 anchoring and is biologically invalid for TRD; the joint
        # history is normalized within J, then marginalized with P(J).
        assert 0.0 < p2 < 1.0, f"{locus}: single/tandem branch collapsed: P(n_D=2)={p2}"

    # Generate a real batch (not 5) so a rare choosable-but-unsamplable allele is actually hit.
    df = __import__("vdjtools.model.generate", fromlist=["generate"]).generate(m, 200, seed=0)
    assert native.pgen_nt(m, df["junction_nt"][0].upper(), df["v_call"][0], df["j_call"][0]) > 0.0


def test_sources_constant():
    assert SOURCES == ("isalgo", "olga-legacy", "olga", "learned", "arda")


# --- the arda-namespace set (organism-keyed; the only one with mouse) -----------------

def test_arda_alias_reaches_canonical_human_and_mouse_family():
    """The old source spelling resolves to fresh isalgo fits, including all five mouse loci."""
    assert set(list_bundled()["isalgo"]) == set(LOCI) | {f"mouse_{loc}" for loc in MOUSE_LOCI}
    assert "arda" not in list_bundled()


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


def test_organism_is_not_silently_ignored_by_a_set_vendored_for_human_only():
    """Asking olga for mouse must fail loudly, not hand back the human model -- and the
    message must say the gap is in what is VENDORED, not in the upstream set.

    ⚠ OLGA is **not** a human-only project: it ships `mouse_T_alpha`, `mouse_T_beta`,
    `mouse_B_heavy`, `mouse_B_kappa` and `mouse_B_lambda`, and all five import cleanly through
    `from_olga`. The old message asserted the set itself was human-only, which is false and sends a
    reader looking for a mouse model somewhere it exists.
    """
    with pytest.raises(ValueError, match="VENDORED"):
        load_bundled("TRB", 'olga', organism="mouse")
    with pytest.raises(ValueError, match="mouse_T_beta"):
        load_bundled("TRB", 'olga', organism="mouse")
    assert load_bundled("IGH", "arda", organism="mouse").organism == "mouse"


def test_canonical_learned_uses_the_arda_namespace():
    """Working fits use ARDA names; OLGA remains an explicit historical source."""
    a = load_bundled("TRB", "arda").tables["v_choice"]["v_allele"].to_list()
    lrn = load_bundled("TRB", "learned").tables["v_choice"]["v_allele"].to_list()
    assert set(a) == set(lrn)
