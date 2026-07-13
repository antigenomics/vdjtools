"""Scaffold smoke tests: the package imports and the C++ core is wired through."""
import vdjtools
from vdjtools import _core


def test_version():
    assert vdjtools.__version__ == "2.2.0"
    assert _core.version() == "2.2.0"


def test_hamming():
    assert _core.hamming("CASSL", "CASSL") == 0
    assert _core.hamming("CASSL", "CASSF") == 1
    assert _core.hamming("CASS", "CASSL") == -1  # length mismatch
    assert vdjtools.hamming("AAA", "ABA") == 1   # re-exported at top level
