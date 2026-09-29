"""The batched junction translation is public, and it agrees with the scalar pair exactly.

``_core.translate_junctions`` shipped in 4.5.0 and made the format readers 5.32x, but it was
reachable only through ``vdjtools._core`` or a private helper -- so a caller holding a nucleotide
column had the scalar ``translate`` and nothing else. That is the shape this repo keeps finding:
a fast path that exists and cannot be reached.

What is pinned here is the equality, not the speed. The scalar ``translate`` + ``to_unified_cdr3aa``
pair is the reference; the C++ has to reproduce it on in-frame sequences, on out-of-frame ones
(where the walk runs inward from both ends), and on non-ACGT input.
"""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from vdjtools.io import translate_junctions
from vdjtools.io.convert import to_unified_cdr3aa, translate


def _reference(seqs):
    return [to_unified_cdr3aa(translate(s)) if s else None for s in seqs]


@pytest.fixture(scope="module")
def cohort():
    """Every frame offset, every base, plus the degenerate inputs, in one list."""
    rng = np.random.default_rng(20260929)
    out = [None, "", "TGT", "TG", "T",
           "TGTGCCAGCAGCTTAGGACAGGCCTACGAGCAGTACTTC"]          # a real in-frame TRB junction
    for n in range(3, 60):                                      # all three offsets at every length
        for _ in range(12):
            out.append("".join(rng.choice(list("ACGT"), n)))
    for _ in range(300):                                        # non-ACGT: N and lower case
        out.append("".join(rng.choice(list("ACGTNacgt"), int(rng.integers(3, 45)))))
    return out


def test_the_batch_equals_the_scalar_pair(cohort):
    assert translate_junctions(cohort).to_list() == _reference(cohort)


def test_a_null_or_empty_nucleotide_stays_null(cohort):
    """A blank junction is absent evidence, not an empty peptide."""
    got = translate_junctions([None, "", "TGTGCC"]).to_list()
    assert got[0] is None and got[1] is None and got[2] == "CA"


def test_it_takes_a_series_or_a_plain_iterable(cohort):
    """A caller holding a column should not have to convert it, and vice versa."""
    a = translate_junctions(pl.Series(cohort))
    b = translate_junctions(cohort)
    c = translate_junctions(tuple(cohort))
    assert isinstance(a, pl.Series) and a.to_list() == b.to_list() == c.to_list()


def test_threads_do_not_change_the_answer(cohort):
    """Parallelism that moves a value is a bug with a speedup, not parallelism."""
    one = translate_junctions(cohort, threads=1).to_list()
    assert one == translate_junctions(cohort, threads=8).to_list()
    assert one == translate_junctions(cohort, threads=0).to_list()


def test_an_out_of_frame_junction_translates_inward_from_both_ends(cohort):
    """The bidirectional walk is the part polars cannot express, so it is pinned by name.

    A one-base insertion leaves both ends in frame and the middle untranslatable; the middle
    collapses to a single ``_`` rather than vanishing, which is what distinguishes this from
    simply dropping the trailing partial codon.
    """
    inframe = "TGTGCCAGCAGCTTAGGACAGGCCTACGAGCAGTACTTC"     # 39 nt -> CASSLGQAYEQYF
    oof = inframe + "A"                                        # 40 nt: one base, both ends framed
    got = translate_junctions([oof]).item()
    assert translate_junctions([inframe]).item() == "CASSLGQAYEQYF"
    assert got == to_unified_cdr3aa(translate(oof)) == "CASSLG_PTSSTS", got

    # The left walk keeps the original frame, so it is still a PREFIX of the in-frame peptide.
    assert got.split("_")[0] == "CASSLG" == "CASSLGQAYEQYF"[:6]
    # The right walk is reframed by the inserted base, so it is NOT the in-frame suffix -- which
    # is the whole reason this is a bidirectional walk and not a trailing-codon drop.
    assert got.split("_")[1] == "PTSSTS" and not "CASSLGQAYEQYF".endswith("PTSSTS")
    # And the untranslatable middle collapses to exactly one marker, never vanishes.
    assert got.count("_") == 1 and translate(oof) == "CASSLGca??ggPTSSTS"
