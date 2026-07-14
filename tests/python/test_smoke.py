"""Scaffold smoke tests: the package imports and the C++ core is wired through."""
import vdjtools
from vdjtools import _core


def test_version():
    assert vdjtools.__version__ == "2.4.0"
    assert _core.version() == "2.4.0"


def test_hamming():
    assert _core.hamming("CASSL", "CASSL") == 0
    assert _core.hamming("CASSL", "CASSF") == 1
    assert _core.hamming("CASS", "CASSL") == -1  # length mismatch
    assert vdjtools.hamming("AAA", "ABA") == 1   # re-exported at top level


def test_arda_is_a_base_dependency():
    """A plain ``pip install vdjtools`` must ship arda — i.e. the model engine AND the germline
    reference. Downstream libraries (mirpy) depend on vdjtools *for* the germline reference and
    cannot be expected to opt into an extra, so arda must never regress to an optional extra.
    """
    from importlib.metadata import requires

    arda_reqs = [r for r in (requires("vdjtools") or [])
                 if r.split(";")[0].strip().lower().startswith("arda-mapper")]
    assert arda_reqs, "arda-mapper is not a vdjtools requirement at all"
    # At least one arda requirement must be unconditional (no `extra == "..."` marker).
    assert any("extra ==" not in r for r in arda_reqs), \
        f"arda-mapper is gated behind an extra, so a plain install has no germline: {arda_reqs}"


def test_plain_install_exposes_the_germline_reference():
    """The contract mirpy relies on: germline V/D/J + CDR3 anchors resolve with no extras."""
    from vdjtools.model.reference import load_germline

    g = load_germline("TRB")
    assert set(g["segment"].unique()) >= {"V", "D", "J"}
    assert g.filter(g["segment"] == "V").height > 0
    assert (g.filter(g["segment"] == "V")["cdr3_anchor"] >= 0).any()
