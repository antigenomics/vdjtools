"""Scaffold smoke tests: the package imports and the C++ core is wired through."""
import pytest

import vdjtools
from vdjtools import _core


def test_version():
    assert vdjtools.__version__ == "2.7.0"
    assert _core.version() == "2.7.0"


def test_no_duplicated_hamming():
    # Generic sequence primitives (Hamming/edit distance) are used from the seqtree
    # dependency (`seqtree.distance`, ≥0.5.0), not duplicated in vdjtools — so
    # `import vdjtools` needs no compiled ext.
    from seqtree import distance

    assert not hasattr(vdjtools, "hamming")
    assert not hasattr(_core, "hamming")
    assert distance.hamming("CASSL", "CASSF") == 1        # the shipped primitive vdjtools relies on
    assert distance.levenshtein("kitten", "sitting") == 3


@pytest.mark.parametrize("engine", ["arda-mapper", "seqtree", "vdjmatch"])
def test_antigenomics_engines_are_base_dependencies(engine):
    """A plain ``pip install vdjtools`` must ship all three antigenomics engines.

    Every capability vdjtools advertises is delegated to one of them and none has a fallback:
    arda → germline reference + model engine (mirpy depends on vdjtools *for* the germline);
    seqtree → fuzzy search / e-values (preprocess.correct, similarity overlap, TCRnet);
    vdjmatch → overlap + TCRnet + metaclonotypes. If any regresses to an optional extra, a
    plain install silently loses documented functionality with an ImportError.
    """
    from importlib.metadata import requires

    reqs = [r for r in (requires("vdjtools") or [])
            if r.split(";")[0].strip().lower().startswith(engine)]
    assert reqs, f"{engine} is not a vdjtools requirement at all"
    # At least one requirement must be unconditional (no `extra == "..."` marker).
    assert any("extra ==" not in r for r in reqs), \
        f"{engine} is gated behind an extra, so a plain install loses its features: {reqs}"


def test_plain_install_exposes_the_germline_reference():
    """The contract mirpy relies on: germline V/D/J + CDR3 anchors resolve with no extras."""
    from vdjtools.model.reference import load_germline

    g = load_germline("TRB")
    assert set(g["segment"].unique()) >= {"V", "D", "J"}
    assert g.filter(g["segment"] == "V").height > 0
    assert (g.filter(g["segment"] == "V")["cdr3_anchor"] >= 0).any()


def test_plain_install_exposes_the_seqtree_engine():
    """Error correction (seqtree near-neighbour search) must work with no extras."""
    import polars as pl

    from vdjtools.io import schema as S
    from vdjtools.preprocess import correct

    df = S.add_locus(S.normalize(pl.DataFrame({
        S.V_CALL: ["TRBV1"] * 2, S.J_CALL: ["TRBJ1"] * 2,
        S.JUNCTION_AA: ["CASSLF", "CASSLF"], S.JUNCTION_NT: ["TGTGCCAGCAGCCTTTTT", "TGTGCCAGCAGCCTTTTA"],
        S.COUNT: [1000, 1],   # the 1-read variant is a 1-mismatch child of the parent
    }), recompute_freq=True))
    out = correct(df)                      # would ImportError if seqtree were an extra
    assert out.height <= df.height
