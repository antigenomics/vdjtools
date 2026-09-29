"""``pgen_nt_batch``: the batch entry point nucleotide Pgen never had.

``pgen_nt`` was the one native Pgen with no batch, so every caller that had more than one sequence
either looped in Python or wrapped that loop in a ``ThreadPoolExecutor`` handing out **one task per
sequence** -- the shape this repo rejects everywhere else. Three call sites did
(``model.score._pgen_nt_many``, the ``pgen`` CLI command, and ``sc.paired_pgen`` for its amino-acid
twin); all three are one batched call now.
"""
from __future__ import annotations

import pytest

from vdjtools.model import native
from vdjtools.model.bundled import load_bundled
from vdjtools.model.generate import generate


@pytest.fixture(scope="module")
def seqs():
    m = load_bundled("TRG", "olga")            # VJ chain: cheap per call, so the test is quick
    got = generate(m, 80, seed=4, productive_only=True)
    return m, got["junction_nt"].to_list(), got["v_call"].to_list(), got["j_call"].to_list()


def test_the_batch_is_the_per_sequence_loop_exactly(seqs):
    """Batching may not move a number. Bitwise, not approximately."""
    m, nt, v, j = seqs
    one = [native.pgen_nt(m, s, a, b) for s, a, b in zip(nt, v, j)]
    many = native.pgen_nt_batch(m, nt, v, j, threads=0)
    assert many == one


def test_marginalizing_over_v_and_j_matches_too(seqs):
    m, nt, _v, _j = seqs
    one = [native.pgen_nt(m, s) for s in nt]
    assert native.pgen_nt_batch(m, nt) == one
    assert native.pgen_nt_batch(m, nt, [None] * len(nt), [None] * len(nt)) == one


def test_the_thread_count_does_not_change_the_answer(seqs):
    m, nt, v, j = seqs
    base = native.pgen_nt_batch(m, nt, v, j, threads=1)
    for threads in (2, 4, 8, 0):
        assert native.pgen_nt_batch(m, nt, v, j, threads=threads) == base, threads


def test_a_non_acgt_base_raises_rather_than_scoring_zero(seqs):
    """A silent 0 is indistinguishable from a real answer, so it has to be an error.

    ``pgen_nt`` raises ``KeyError`` from its own encoding table; the batch encodes in C++ and
    raises ``ValueError`` naming the offending row. Both refuse; neither invents a number.
    """
    m, nt, v, j = seqs
    bad = list(nt)
    bad[3] = bad[3][:-1] + "N"
    with pytest.raises(ValueError, match="non-ACGT"):
        native.pgen_nt_batch(m, bad, v, j)


def test_a_wrong_length_call_list_raises(seqs):
    m, nt, v, j = seqs
    with pytest.raises(ValueError, match="same length"):
        native.pgen_nt_batch(m, nt, v[:-1], j)
    with pytest.raises(ValueError, match="same length"):
        native.pgen_nt_batch(m, nt, v, j[:-1])


def test_a_gene_level_call_still_refuses(seqs):
    """The model is keyed by allele and must not degrade a gene name to the agnostic value."""
    m, nt, v, j = seqs
    with pytest.raises(KeyError):
        native.pgen_nt_batch(m, nt, [x.split("*")[0] for x in v], j)


def test_an_empty_batch_is_an_empty_list(seqs):
    m, *_ = seqs
    assert native.pgen_nt_batch(m, []) == []
