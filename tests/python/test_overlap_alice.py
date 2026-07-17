"""Tests for vdjtools.overlap.alice — neighbourhood enrichment against the Pgen null.

A planted convergent cluster must light up; sequences drawn from the generation model itself
must not (they ARE the null). The gene-level V/J path is pinned because every real cohort ships
gene-level calls and the native Pgen is keyed by allele.
"""
import polars as pl
import pytest

from vdjtools import overlap as O
from vdjtools.io import schema as S

pytest.importorskip("vdjmatch")

from vdjtools.model import load_bundled                       # noqa: E402
from vdjtools.model.native import pack, pgen_aa               # noqa: E402
from vdjtools.overlap.alice import _alleles_for, _ball_pgen   # noqa: E402

_AA = "ACDEFGHIKLMNPQRSTVWY"


def _sample(cdr3, v="TRBV19*01", j="TRBJ2-1*01", nt=None, count=5):
    n = len(cdr3)
    d = {S.V_CALL: [v] * n, S.J_CALL: [j] * n, S.JUNCTION_AA: list(cdr3),
         S.COUNT: [count] * n}
    if nt is not None:
        d[S.JUNCTION_NT] = list(nt)
    return S.add_locus(S.normalize(pl.DataFrame(d), recompute_freq=True))


# --- the Pgen ball ------------------------------------------------------------------------

def test_ball_pgen_is_the_closed_hamming1_ball():
    # ALICE's lambda sums Pgen over the ball INCLUDING sigma, so the degree must include sigma
    # too. Verify the native identity against brute force rather than trusting the docstring:
    # an OPEN ball would be a systematic under-count of lambda, i.e. anti-conservative.
    m = load_bundled("TRB", "olga", collapse=False)
    seq, v, j = "CASSYVGSEQFF", "TRBV19*01", "TRBJ2-1*01"
    centre = pgen_aa(m, seq, v, j)
    nbrs = [seq[:k] + r + seq[k + 1:] for k in range(len(seq)) for r in _AA if r != seq[k]]
    brute_open = sum(pgen_aa(m, s, v, j) for s in nbrs)
    got = float(_ball_pgen(m, [seq], v, j, 0)[0])
    assert got == pytest.approx(brute_open + centre, rel=1e-9)   # closed
    assert abs(got / brute_open - 1) > 1e-4                      # and NOT the open ball


def test_gene_level_call_sums_its_alleles_exactly():
    # Real cohorts ship gene-level v_call; the model is keyed by allele and the native Pgen now
    # RAISES on a gene (it used to marginalise over every allele, 183x too high). Summing a
    # gene's alleles is exact, not an approximation -- pin that.
    m = load_bundled("TRB", "olga", collapse=False)
    _, vi, _ = pack(m)
    gene = "TRBV7-9"
    alleles = _alleles_for(vi, gene)
    assert len(alleles) > 1 and all(a.startswith(gene + "*") for a in alleles)

    seq, j = "CASSYVGSEQFF", "TRBJ2-1*01"
    got = float(_ball_pgen(m, [seq], gene, j, 0)[0])
    want = sum(float(_ball_pgen(m, [seq], a, j, 0)[0]) for a in alleles)
    assert got == pytest.approx(want, rel=1e-12)


def test_alleles_for_rejects_an_unknown_call():
    _, vi, _ = pack(load_bundled("TRB", "olga", collapse=False))
    assert _alleles_for(vi, "TRBV19*01") == ["TRBV19*01"]        # an allele is itself
    with pytest.raises(KeyError, match="not in the model"):
        _alleles_for(vi, "TRBV999")


# --- the test itself ----------------------------------------------------------------------

def test_alice_flags_a_planted_convergent_cluster():
    # A clique of mutually-1-substitution CDR3s is what antigen-driven convergence looks like.
    # It must clear the Pgen null; sequences drawn FROM the generation model must not, because
    # they are literally the null this test is against.
    m = load_bundled("TRB", "olga", collapse=False)
    from vdjtools.model.generate import generate

    gen = generate(m, 4000, seed=0)
    bg = [s for s in gen[S.JUNCTION_AA].to_list()
          if s and "*" not in s and "_" not in s and len(s) == 13][:250]
    tmpl = list("CASSXPGELFFYT")
    clique = ["".join(tmpl[:4] + [r] + tmpl[5:]) for r in "ACDEFGHIK"]

    res = O.alice(_sample(clique + bg), locus="TRB")
    got = dict(zip(res[S.JUNCTION_AA].to_list(), res["p_enrichment"].to_list()))
    planted = [p for s, p in got.items() if s in clique]
    assert len(planted) >= len(clique) - 1                  # nearly all clique members tested
    assert max(planted) < 1e-6                              # and all wildly enriched
    other = [p for s, p in got.items() if s not in clique]
    assert not other or min(other) > max(planted)           # the null sits below the signal


def test_alice_counts_nucleotide_clonotypes_not_amino_acids():
    # ALICE's units are NUCLEOTIDE clonotypes: distinct nt variants of one aa sequence are
    # genuine neighbours (they are different T-cell clones), and convergent recombination is
    # exactly the signal. Two nt variants of each aa must double the degree.
    aa = ["".join(list("CASSXPGELFFYT")[:4] + [r] + list("CASSXPGELFFYT")[5:]) for r in "ACDEF"]
    one = _sample(aa, nt=[f"TGT{i:09d}" for i in range(len(aa))])
    two = _sample(aa + aa, nt=[f"TGT{i:09d}" for i in range(2 * len(aa))])
    d1 = O.alice(one, locus="TRB")["n_neighbors"].max()
    d2 = O.alice(two, locus="TRB")["n_neighbors"].max()
    assert d2 == 2 * d1


def test_alice_min_count_drops_error_variants():
    # A 1-read variant of an abundant clonotype is usually sequencing error; counting it as a
    # neighbour inflates every degree in its neighbourhood.
    aa = ["".join(list("CASSXPGELFFYT")[:4] + [r] + list("CASSXPGELFFYT")[5:]) for r in "ACDEFGH"]
    big = _sample(aa, count=10)
    small = _sample(["CASSZPGELFFYT"], count=1)
    mixed = pl.concat([big, small], how="vertical_relaxed")
    res = O.alice(mixed, locus="TRB", min_count=2)
    assert "CASSZPGELFFYT" not in res[S.JUNCTION_AA].to_list()


def test_alice_returns_q_and_picks_no_threshold():
    aa = ["".join(list("CASSXPGELFFYT")[:4] + [r] + list("CASSXPGELFFYT")[5:]) for r in "ACDEFGH"]
    res = O.alice(_sample(aa), locus="TRB")
    assert "q_value" in res.columns
    assert (res["q_value"].to_numpy() >= res["p_enrichment"].to_numpy() - 1e-12).all()
    # E = n * Q * pgen_ball, reconstructible from the returned columns
    r = res.row(0, named=True)
    assert r["E"] == pytest.approx(r["n_group"] * 9.41 * r["pgen_ball"], rel=1e-9)


def test_alice_min_degree_gates_singletons():
    # ALICE tests only d(sigma) > 2: a clonotype with no neighbours carries no evidence either
    # way, and testing it only spends multiple-testing budget.
    aa = ["".join(list("CASSXPGELFFYT")[:4] + [r] + list("CASSXPGELFFYT")[5:]) for r in "ACDEFGH"]
    lone = ["CWWWWWWWWWWWF"]
    res = O.alice(_sample(aa + lone), locus="TRB", min_degree=3)
    assert "CWWWWWWWWWWWF" not in res[S.JUNCTION_AA].to_list()
    assert (res["n_neighbors"] >= 3).all()
