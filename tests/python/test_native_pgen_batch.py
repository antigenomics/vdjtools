"""Native batch aa-Pgen: exact vs per-sequence + thread-count invariant (Phase 13 item 5)."""
import pytest

from vdjtools.model import load_bundled, native
from vdjtools.model.generate import generate


@pytest.fixture(scope="module")
def data():
    m = load_bundled("TRB", "olga")
    g = generate(m, 128, seed=5, productive_only=True)  # >64 so the threaded path actually threads
    return m, g["junction_aa"].to_list(), g["v_call"].to_list(), g["j_call"].to_list()


@pytest.mark.parametrize("mm", [0, 1])
def test_batch_equals_serial_and_is_thread_invariant(data, mm):
    m, seqs, vs, js = data
    serial = [native.pgen_aa(m, s, v, j, mismatches=mm) for s, v, j in zip(seqs, vs, js)]
    # bitwise-identical to per-sequence, and independent of thread count.
    assert native.pgen_aa_batch(m, seqs, vs, js, mismatches=mm, threads=1) == serial
    assert native.pgen_aa_batch(m, seqs, vs, js, mismatches=mm, threads=8) == serial
    assert native.pgen_aa_batch(m, seqs, vs, js, mismatches=mm, threads=0) == serial


def test_batch_gene_agnostic(data):
    m, seqs, _, _ = data
    serial = [native.pgen_aa(m, s) for s in seqs]
    assert native.pgen_aa_batch(m, seqs, mismatches=0, threads=8) == serial


def test_batch_validates_args(data):
    m, seqs, vs, _ = data
    with pytest.raises(ValueError):
        native.pgen_aa_batch(m, seqs, mismatches=2)
    with pytest.raises(ValueError):
        native.pgen_aa_batch(m, seqs, v=vs[:3], mismatches=0)  # v length != seqs
