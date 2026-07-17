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


# --- an unknown V/J must RAISE, never silently marginalize -------------------------------------
# Regression: `vi.get(v, -1)` mapped any unrecognised name to -1 = "marginalize over all V/J", so a
# gene-level call returned the V/J-agnostic Pgen with no error. Real repertoires carry gene-level
# v_call (`TRBV9`), so this silently overstated Pgen 2.38x on the AS/B27 motif.

def test_gene_level_name_raises_not_silently_marginalizes():
    m = load_bundled("TRB", "olga")
    s = "CASSVGLYSTDTQYF"
    agnostic = native.pgen_aa(m, s, None, None)
    allele = native.pgen_aa(m, s, "TRBV9*01", "TRBJ2-3*01")
    assert allele < agnostic  # conditioning on one allele must lose mass
    with pytest.raises(KeyError, match="gene name"):
        native.pgen_aa(m, s, "TRBV9", "TRBJ2-3")
    with pytest.raises(KeyError, match="gene name"):
        native.pgen_aa_batch(m, [s], ["TRBV9"], ["TRBJ2-3*01"])
    with pytest.raises(KeyError, match="gene name"):
        native.pgen_nt(m, "TGTGCCAGCAGCGTAGGGCTTTATTCGACAGATACGCAGTATTTT", "TRBV9", "TRBJ2-3*01")


@pytest.mark.parametrize("bad", ["TRBV999*01", "", "not-a-gene"])
def test_unknown_allele_raises(bad):
    m = load_bundled("TRB", "olga")
    if bad == "":  # falsy => the documented "marginalize" path, not an error
        assert native.pgen_aa(m, "CASSVGLYSTDTQYF", bad, None) == native.pgen_aa(
            m, "CASSVGLYSTDTQYF", None, None)
        return
    with pytest.raises(KeyError):
        native.pgen_aa(m, "CASSVGLYSTDTQYF", bad, None)
