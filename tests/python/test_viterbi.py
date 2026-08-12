"""`model.viterbi` — the argmax side of the Pgen DP.

NOTE: Every generated draw used here passes ``productive_only=True``. An out-of-frame rearrangement has
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
    """WARNING: `scenario_p <= pgen_nt`. A maximum cannot exceed the sum it is taken over, so a violation
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
    """NOTE: An out-of-frame CDR3 (len(nt) != 3*len(aa)) has no in-frame reconstruction. The right
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


@pytest.mark.parametrize("locus", LOCI)
def test_infer_nt_reproduces_the_exact_oracle(preps, locus):
    """WARNING: The headline contract: on every record the exponential oracle can actually resolve,
    ``infer_nt`` must return the same sequence.

    This is what separates the search from the cheap alternative. Fixing the germline to its
    most likely trim and then picking the best codon per residue — the obvious shortcut — agrees
    with the oracle on only 9/25 (TRG) and 4/19 (TRA) of these, because a trim chosen before the
    codons pins a codon the true optimum would have trimmed away. Same unsoundness as the
    germline-pinning trick the module header rejects.
    """
    m, prep = preps[locus]
    g = generate(m, 300, seed=5, productive_only=True)
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
        want = infer_nt_bruteforce(prep, aa, v, j, max_candidates=20_000)
        if want is None or want.truncated:
            continue
        got = infer_nt(m, aa, v, j)
        tested += 1
        assert got is not None and got.cdr3_nt == want.cdr3_nt, (
            f"{aa}: oracle {want.cdr3_nt} (pgen={want.pgen:.4e}), got "
            f"{got.cdr3_nt if got else None} (pgen={got.pgen if got else 0:.4e})"
        )
        assert got.pgen == pytest.approx(want.pgen, rel=1e-9)
    assert tested >= 3, f"only {tested} records small enough to brute-force — widen the fixture"


@pytest.mark.parametrize("locus", LOCI)
def test_infer_nt_returns_a_sequence_that_translates_back(preps, locus):
    """Whatever else it does, the answer must encode the amino acids it was asked about, be
    explicable under the model, and carry a markup consistent with its own sequence."""
    from vdjtools.model.reference import translate
    m, prep = preps[locus]
    g = generate(m, 60, seed=29, productive_only=True)
    n = 0
    for r in g.iter_rows(named=True):
        aa, v, j = r["junction_aa"], r["v_call"], r["j_call"]
        sc = infer_nt(m, aa, v, j)
        if sc is None:
            continue
        n += 1
        assert translate(sc.cdr3_nt) == aa
        assert sc.v_call == v and sc.j_call == j
        assert prep.cut["v"][sc.v_call].startswith(sc.cdr3_nt[:sc.v_end])
        assert prep.cut["j"][sc.j_call].endswith(sc.cdr3_nt[sc.j_start:])
        if sc.d_call is not None and sc.d_start is not None:
            assert sc.cdr3_nt[sc.d_start:sc.d_end] in prep.cut["d"][sc.d_call]
            assert prep.p_d_given_j.get((sc.j_call, sc.d_call), 0.0) > 0.0
        # a single path can never out-score the sum it belongs to
        assert sc.scenario_p <= sc.pgen * (1 + 1e-9)
        assert sc.pgen == pytest.approx(pgen_nt(prep, sc.cdr3_nt, v, j), rel=1e-9)
    assert n > 30, f"only {n} records resolved — the fixture is not exercising infer_nt"


def test_infer_nt_refuses_what_it_cannot_encode(preps):
    """An unknown residue has no codons, so there is no sequence to return. ``None``, not a guess."""
    m, prep = preps["TRA"]
    assert infer_nt(m, "CAXVF") is None
    assert infer_nt(m, "") is None


def test_the_three_call_input_modes(preps):
    """``v=``/``j=`` accept one allele, several (comma string or list), or nothing.

    The multi case is what an AIRR ``v_call`` carries when the aligner could not choose, and it is
    NOT a pin: a list containing the true allele must reach the same answer as pinning it, since
    the extra alleles only widen the search.
    """
    m, prep = preps["TRB"]
    g = generate(m, 30, seed=41, productive_only=True)
    checked = 0
    for r in g.iter_rows(named=True):
        aa, v, j = r["junction_aa"], r["v_call"], r["j_call"]
        pinned = infer_nt(m, aa, v, j)
        if pinned is None:
            continue
        checked += 1
        other = next(a for a in prep.functional_v if a != v)
        for spec in (f"{v},{other}", [v, other]):
            multi = infer_nt(m, aa, spec, j)
            assert multi is not None
            assert multi.pgen >= pinned.pgen * (1 - 1e-9), (
                f"widening V to {spec} lost probability: {multi.pgen:.4e} < {pinned.pgen:.4e}"
            )
        assert infer_nt(m, aa, [], j) is None      # an empty candidate set explains nothing
        if checked >= 3:
            break
    assert checked >= 3


def test_unknown_calls_still_return_something_plausible(preps):
    """With no V/J the search starts from a ``top_vj`` shortlist. That is a heuristic and may miss
    the true optimum, but it must still return a sequence that translates back and is explicable."""
    from vdjtools.model.reference import translate
    m, prep = preps["TRA"]                              # VJ: no D enumeration, so this stays quick
    g = generate(m, 12, seed=43, productive_only=True)
    n = 0
    for r in g.iter_rows(named=True):
        sc = infer_nt(m, r["junction_aa"])
        if sc is None:
            continue
        n += 1
        assert translate(sc.cdr3_nt) == r["junction_aa"]
        assert sc.pgen > 0.0
    assert n >= 5, f"only {n} resolved without calls"


def test_native_search_agrees_with_the_python_reference(preps):
    """WARNING: The native argmax DP and the pure-Python scenario search must return the same
    sequence — they are independent implementations of the same maximum.

    The native one is a Pi_L*Pi_R transfer matrix over the whole scenario space; the Python one
    enumerates ``(V, delV) x (J, delJ) x (D, delD, pos)`` explicitly with branch-and-bound. They
    share no code beyond the codon DP, so agreement is real evidence rather than a tautology.
    Kept on TRA (VJ) because the Python path is ~600x slower on a VDJ locus.
    """
    m, prep = preps["TRA"]
    g = generate(m, 15, seed=47, productive_only=True)
    n = 0
    for r in g.iter_rows(named=True):
        aa, v, j = r["junction_aa"], r["v_call"], r["j_call"]
        fast, ref = infer_nt(m, aa, v, j), infer_nt(prep, aa, v, j)
        if fast is None or ref is None:
            assert fast is None and ref is None
            continue
        n += 1
        assert fast.cdr3_nt == ref.cdr3_nt, f"{aa}: native {fast.cdr3_nt} vs python {ref.cdr3_nt}"
        assert fast.pgen == pytest.approx(ref.pgen, rel=1e-9)
        assert (fast.v_call, fast.j_call, fast.v_end, fast.j_start) == (
            ref.v_call, ref.j_call, ref.v_end, ref.j_start)
    assert n >= 8, f"only {n} compared"


def test_keep_widens_the_pool_without_changing_the_answer(preps):
    """``keep`` is a speed/pool knob, not a correctness one: relaxing it may only ADD candidates.

    The argmax-only setting (``keep=1``) is allowed to differ — it skips the marginal re-score,
    which is exactly what the wider pool is for — but a wider pool must never return a sequence the
    model likes *less* than a narrower one did.
    """
    m, prep = preps["TRA"]
    g = generate(m, 40, seed=31, productive_only=True)
    for r in g.iter_rows(named=True):
        aa, v, j = r["junction_aa"], r["v_call"], r["j_call"]
        narrow = infer_nt(prep, aa, v, j, keep=1.0, n_best=1)
        wide = infer_nt(prep, aa, v, j, keep=1e-4, n_best=32)
        if narrow is None or wide is None:
            assert narrow is None and wide is None
            continue
        assert wide.n_candidates >= narrow.n_candidates
        assert wide.pgen >= narrow.pgen * (1 - 1e-9)


def test_the_chosen_D_is_genomically_possible_with_the_chosen_J(preps):
    """WARNING: Regression: an earlier draft picked the longest exact D substring and IGNORED ``j``.

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
