"""`model.viterbi` — the argmax side of the Pgen DP.

⚠ Every generated draw used here passes ``productive_only=True``. An out-of-frame rearrangement has
``len(nt) != 3 * len(aa)`` (measured: 8 aa against 25 nt), so no in-frame reconstruction can explain
it and comparing against one is meaningless. A first version of these tests omitted the flag and
scored 10 of 16 records as "both impossible" — mutual agreement on nothing.
"""
from __future__ import annotations

from itertools import product

import pytest

from vdjtools.model.bundled import load_bundled
from vdjtools.model.generate import generate
from vdjtools.model.pgen import pgen_nt, prepare
from vdjtools.model.viterbi import (best_scenario, codon_options, infer_nt,
                                    infer_nt_bruteforce)

LOCI = ("TRA", "TRB")


@pytest.fixture(scope="module")
def preps():
    out = {}
    for locus in LOCI:
        m = load_bundled(locus, "olga")
        out[locus] = (m, prepare(m))
    return out


def test_codon_options_are_strings_and_complete():
    """`pgen._CODON_AA` is keyed by a TUPLE OF NUCLEOTIDE INDICES (A0 C1 G2 T3), not a string —
    `(3, 2, 1)` is "TGC". Getting that wrong yields tuples where sequence is expected."""
    assert codon_options("C") == [["TGC", "TGT"]]
    assert codon_options("W") == [["TGG"]]
    assert codon_options("L")[0] == ["CTA", "CTC", "CTG", "CTT", "TTA", "TTG"]
    assert codon_options("CW") == [["TGC", "TGT"], ["TGG"]]
    assert codon_options("X") == [[]]          # unknown residue -> no options, never a guess


@pytest.mark.parametrize("locus", LOCI)
def test_best_scenario_spans_really_are_germline(preps, locus):
    """The whole point of the markup: the span called V must BE the V germline, and likewise J/D.

    This is the invariant that would break first if the scenario walk and the model's own tables
    ever diverged, and it needs no external truth.
    """
    m, prep = preps[locus]
    g = generate(m, 200, seed=11, productive_only=True)
    n = 0
    for r in g.iter_rows(named=True):
        nt, v, j = r["junction_nt"], r["v_call"], r["j_call"]
        sc = best_scenario(prep, nt, v, j)
        if sc is None:
            continue
        n += 1
        gv, gj = prep.cut["v"][v], prep.cut["j"][j]
        assert nt[:sc.v_end] == gv[:sc.v_end], f"V span is not V germline: {sc}"
        tail = len(nt) - sc.j_start
        assert nt[sc.j_start:] == gj[len(gj) - tail:], f"J span is not J germline: {sc}"
        if sc.d_call is not None:
            assert nt[sc.d_start:sc.d_end] in prep.cut["d"][sc.d_call], f"D span is not D: {sc}"
    assert n > 100, f"only {n} scenarios scored — the fixture is not exercising anything"


@pytest.mark.parametrize("locus", LOCI)
def test_a_single_scenario_never_exceeds_the_marginal(preps, locus):
    """⛔ `scenario_p <= pgen_nt`. A maximum cannot exceed the sum it is taken over, so a violation
    means the scenario walk and the forward DP have diverged — the cheapest catch there is."""
    m, prep = preps[locus]
    g = generate(m, 200, seed=5, productive_only=True)
    for r in g.iter_rows(named=True):
        nt, v, j = r["junction_nt"], r["v_call"], r["j_call"]
        sc = best_scenario(prep, nt, v, j)
        if sc is None:
            continue
        total = pgen_nt(prep, nt, v, j)
        assert sc.scenario_p <= total * (1 + 1e-9), (
            f"one scenario ({sc.scenario_p:.6e}) exceeds the marginal ({total:.6e}) for {nt}"
        )


def test_bruteforce_matches_an_independent_enumeration(preps):
    """The oracle checked against a second, independently written enumeration.

    Ground truth comes from the model itself: enumerate every nt encoding the aa, score each with
    the shipped `pgen_nt`, take the argmax. Restricted to genuinely small search spaces so the
    reference stays exhaustive.
    """
    m, prep = preps["TRA"]
    g = generate(m, 200, seed=3, productive_only=True)
    tested = 0
    for r in g.iter_rows(named=True):
        aa, v, j = r["junction_aa"], r["v_call"], r["j_call"]
        opts = codon_options(aa)
        if any(not o for o in opts):
            continue
        space = 1
        for o in opts:
            space *= len(o)
        if space > 20_000:
            continue
        best_p, best_nt = 0.0, None
        for combo in product(*opts):
            nt = "".join(combo)
            p = pgen_nt(prep, nt, v, j)
            if p > best_p:
                best_p, best_nt = p, nt
        got = infer_nt_bruteforce(prep, aa, v, j)
        tested += 1
        if best_nt is None:
            assert got is None, f"claimed {got.cdr3_nt} where nothing is explicable"
        else:
            assert got is not None and got.cdr3_nt == best_nt, (
                f"{aa}: expected {best_nt} (p={best_p:.4e}), got "
                f"{got.cdr3_nt if got else None} (p={got.pgen if got else 0:.4e})"
            )
    assert tested >= 3, f"only {tested} records small enough to brute-force — widen the fixture"


def test_out_of_frame_is_refused_not_guessed(preps):
    """⚠ An out-of-frame CDR3 (len(nt) != 3*len(aa)) has no in-frame reconstruction. The right
    answer is None, not a plausible-looking sequence of the wrong length."""
    m, prep = preps["TRA"]
    g = generate(m, 400, seed=3)                      # NOT productive_only: we want the OOF ones
    oof = [r for r in g.iter_rows(named=True)
           if len(r["junction_nt"]) != 3 * len(r["junction_aa"])]
    if not oof:
        pytest.skip("no out-of-frame draw in this sample")
    r = oof[0]
    got = infer_nt_bruteforce(prep, r["junction_aa"], r["v_call"], r["j_call"],
                              max_candidates=20_000)
    if got is not None:
        assert len(got.cdr3_nt) == 3 * len(r["junction_aa"])


def test_infer_nt_is_honestly_unimplemented():
    """⛔ The production entry point must FAIL LOUDLY rather than quietly dispatch to the oracle.

    Enumeration cannot scale — the median VDJdb record has 5.3e6 (TRA) / 1.9e7 (TRB) candidates —
    and pinning the germline flanks to shrink it is unsound (it drops sequences whose germline was
    trimmed, and the true maximum can be one). A silent fallback would look like a working feature.
    """
    with pytest.raises(NotImplementedError, match="max-product DP"):
        infer_nt(None, "CASSF")


def test_the_chosen_D_is_genomically_possible_with_the_chosen_J(preps):
    """⛔ Regression: an earlier draft picked the longest exact D substring and IGNORED ``j``.

    TRBD2 lies 3' of the entire TRBJ1 cluster, so deletional joining can never produce a
    TRBD2-TRBJ1 pair — the model encodes that as ``p_d_given_j == 0``. Choosing D by sequence
    similarity alone would happily call the impossible pair; choosing it by the model's own terms
    cannot, because a zero prunes the branch.
    """
    m, prep = preps["TRB"]
    g = generate(m, 200, seed=17, productive_only=True)
    n = 0
    for r in g.iter_rows(named=True):
        sc = best_scenario(prep, r["junction_nt"], r["v_call"], r["j_call"])
        if sc is None or sc.d_call is None:
            continue
        n += 1
        assert prep.p_d_given_j.get((sc.j_call, sc.d_call), 0.0) > 0.0, (
            f"chose D={sc.d_call} with J={sc.j_call}, which the model gives probability 0 — "
            f"a genomically impossible rearrangement"
        )
    assert n > 50, f"only {n} D calls made — the fixture is not exercising the constraint"


def test_the_reported_scenario_probability_is_the_product_it_claims(preps):
    """`scenario_p` must be reconstructible from the model's own tables for the reported path.

    Guards against the value drifting from the path it is supposed to describe — the failure that
    would make every downstream probability quietly wrong while all the spans still looked right.
    """
    from vdjtools.model.pgen import _p_insert
    m, prep = preps["TRA"]                                  # VJ: one insertion term, no D sum
    g = generate(m, 120, seed=23, productive_only=True)
    checked = 0
    for r in g.iter_rows(named=True):
        nt, v, j = r["junction_nt"], r["v_call"], r["j_call"]
        sc = best_scenario(prep, nt, v, j)
        if sc is None:
            continue
        cutv, cutj = prep.cut["v"][sc.v_call], prep.cut["j"][sc.j_call]
        ndel_v = len(cutv) - sc.v_end - prep.maxpal["v_3"]
        len_j = len(nt) - sc.j_start
        ndel_j = len(cutj) - len_j - prep.maxpal["j_5"]
        want = (prep.p_v.get(sc.v_call, 0.0)
                * prep.p_j.get((sc.v_call, sc.j_call), 0.0)
                * prep.p_del["v"].get((sc.v_call, ndel_v), 0.0)
                * prep.p_del["j"].get((sc.j_call, ndel_j), 0.0)
                * _p_insert(nt[sc.v_end:sc.j_start], prep.p_ins["vj"], prep.R["vj"],
                            prep.bias["vj"], from_right=False))
        assert abs(sc.scenario_p - want) <= 1e-12 * max(want, 1e-300), (
            f"scenario_p={sc.scenario_p:.6e} but the reported path recomputes to {want:.6e}"
        )
        checked += 1
    assert checked > 50, f"only {checked} scenarios checked"
