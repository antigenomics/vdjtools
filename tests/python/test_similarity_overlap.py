"""Tests for vdjtools.overlap.similarity — sequence-similarity-weighted (TINA) overlap.

The design's correctness anchor is that the continuous ``exp`` kernel has two *exact*
special cases on the same code path: ``Z = I`` (``kernel="identity"``) must reproduce the
classical cosine / Morisita-Horn of the shared-frequency vectors to machine precision, and
``Z = 1[P ≤ θ]`` (``kernel="step"``) must reproduce vdjmatch's fuzzy edit-distance overlap.
Both are pinned here, alongside a hand-computed 2×2, the metric axioms, and the cluster
integration. seqtree / vdjmatch / scipy (the ``overlap`` extra) are ``importorskip``-guarded.
"""
import numpy as np
import polars as pl
import pytest

from vdjtools.io import schema as S
from vdjtools import overlap as O

pytest.importorskip("seqtree")


def _sample(cdr3, counts, v=None, j=None):
    n = len(cdr3)
    df = pl.DataFrame({
        S.V_CALL: v or ["TRBV1"] * n, S.J_CALL: j or ["TRBJ1"] * n,
        S.CDR3_AA: cdr3, S.COUNT: counts,
    })
    return S.add_locus(S.normalize(df, recompute_freq=True))


def _classic(a_cdr3, a_cnt, b_cdr3, b_cnt):
    """Classical cosine and Morisita-Horn of the shared-frequency vectors (Z=I)."""
    pa = np.array(a_cnt, float) / np.sum(a_cnt)
    qb = np.array(b_cnt, float) / np.sum(b_cnt)
    keys = list(dict.fromkeys(a_cdr3 + b_cdr3))
    va = np.array([pa[a_cdr3.index(k)] if k in a_cdr3 else 0.0 for k in keys])
    vb = np.array([qb[b_cdr3.index(k)] if k in b_cdr3 else 0.0 for k in keys])
    cos = float(va @ vb / np.sqrt((va @ va) * (vb @ vb)))
    mh = float(2 * (va @ vb) / ((va @ va) + (vb @ vb)))
    return cos, mh


# ---- Fixture shared by the identity-exactness tests ----
_A_CDR3 = ["CASSLGQAYEQYF", "CASSPGTEAFF", "CASSLAPGELFF"]
_A_CNT = [50, 30, 20]
_B_CDR3 = ["CASSLGQAYEQYF", "CASSLGTEAFF", "CASSWRDNEQFF"]
_B_CNT = [40, 35, 25]


def test_identity_kernel_recovers_exact_cosine():
    """Z=I cosine == classical cosine of the shared-frequency vectors, to 1e-12.

    This exact reduction is the correctness anchor of the whole design.
    """
    a = _sample(_A_CDR3, _A_CNT)
    b = _sample(_B_CDR3, _B_CNT)
    cos_classic, _ = _classic(_A_CDR3, _A_CNT, _B_CDR3, _B_CNT)

    res = O.similarity_overlap(a, b, kernel="identity", metric="cosine")
    # Pinned exact value (verified in Python before committing).
    assert cos_classic == pytest.approx(0.5523681766661075, abs=1e-15)
    # similarity_overlap reproduces the classical metric to machine precision.
    assert res["similarity"] == pytest.approx(cos_classic, rel=1e-12)


def test_identity_kernel_recovers_exact_morisita():
    """Z=I Morisita == 2Σpq/(Σp²+Σq²) of the shared-frequency vectors, to 1e-12."""
    a = _sample(_A_CDR3, _A_CNT)
    b = _sample(_B_CDR3, _B_CNT)
    _, mh_classic = _classic(_A_CDR3, _A_CNT, _B_CDR3, _B_CNT)

    res = O.similarity_overlap(a, b, kernel="identity", metric="morisita")
    assert mh_classic == pytest.approx(0.5517241379310345, abs=1e-15)
    assert res["similarity"] == pytest.approx(mh_classic, rel=1e-12)


def test_step_kernel_matches_vdjmatch_fuzzy_mass():
    """kernel='step', max_penalty=1 → pᵀZq equals vdjmatch's matched-pair frequency mass.

    Reproduces vdjmatch scope '1,0,0,1' (single substitution, no indels): the fixture
    includes a length-differing near pair (CASSAAAAF / CASSAAAF) that BOTH sides must
    exclude (indels prohibited), so this also guards the substitution-only reduction.
    """
    vc = pytest.importorskip("vdjmatch.cluster")
    a_cdr3 = ["CASSLGQAYEQYF", "CASSPGTEAFF", "CASSAAAAF"]
    a_cnt = [50, 30, 20]
    b_cdr3 = ["CASSLGQAYEQYF", "CASSPGTEAFY", "CASSAAAF"]
    b_cnt = [40, 35, 25]
    a = _sample(a_cdr3, a_cnt)
    b = _sample(b_cdr3, b_cnt)

    res = O.similarity_overlap(a, b, kernel="step", max_penalty=1, metric="cosine")

    pa = np.array(a_cnt, float) / np.sum(a_cnt)
    qb = np.array(b_cnt, float) / np.sum(b_cnt)
    pairs = vc.overlap(a_cdr3, b_cdr3, scope="1,0,0,1")
    mass = sum(pa[r["a_idx"]] * qb[r["b_idx"]] for r in pairs.iter_rows(named=True))

    assert mass == pytest.approx(0.305, abs=1e-12)          # verified pinned value
    assert res["pTZq"] == pytest.approx(mass, abs=1e-9)


def test_hand_computed_2x2_cosine():
    """Two 2-clonotype repertoires, off-diagonal penalty 14, τ=14 → Z_off = e⁻¹.

    The CDR3 pair CASSAEQYF / CASSDEQYF differs only by A↔D, whose BLOSUM62 Gram penalty
    is exactly 14, so the 2×2 kernel is [[1, e⁻¹], [e⁻¹, 1]] and the cosine is analytic.
    """
    s1, s2 = "CASSAEQYF", "CASSDEQYF"        # A↔D → penalty 14
    a = _sample([s1, s2], [7, 3])             # p = (0.7, 0.3)
    b = _sample([s1, s2], [4, 6])             # q = (0.4, 0.6)

    res = O.similarity_overlap(a, b, kernel="exp", tau=14, metric="cosine")

    off = np.exp(-1.0)
    Z = np.array([[1.0, off], [off, 1.0]])
    p, q = np.array([0.7, 0.3]), np.array([0.4, 0.6])
    hand = (p @ Z @ q) / np.sqrt((p @ Z @ p) * (q @ Z @ q))
    assert hand == pytest.approx(0.9208164985527538, abs=1e-12)   # matches TINA literature
    assert res["similarity"] == pytest.approx(hand, abs=1e-9)


def test_axioms_bounds_selfsim_symmetry():
    """0 ≤ S ≤ 1, S(a,a)==1, S(a,b)==S(b,a), distance == 1−S ≥ 0."""
    a = _sample(_A_CDR3, _A_CNT)
    b = _sample(_B_CDR3, _B_CNT)
    for metric in ("cosine", "morisita"):
        ab = O.similarity_overlap(a, b, kernel="exp", metric=metric)
        ba = O.similarity_overlap(b, a, kernel="exp", metric=metric)
        aa = O.similarity_overlap(a, a, kernel="exp", metric=metric)
        assert 0.0 <= ab["similarity"] <= 1.0
        assert aa["similarity"] == pytest.approx(1.0, rel=1e-12)
        assert ab["similarity"] == pytest.approx(ba["similarity"], rel=1e-12)
        assert ab["distance"] >= 0.0


def test_near_neighbour_raises_similarity_above_exact():
    """When two samples share only near-variants (no exact clonotype), the exp kernel's
    similarity is strictly above the exact-overlap (identity) value of 0."""
    a = _sample(["CASSLGQAYEQYF", "CASSPGTEAFF"], [60, 40])
    b = _sample(["CASSLGQAYEQYY", "CASSPGTEAFY"], [55, 45])   # each a 1-sub variant

    exact = O.similarity_overlap(a, b, kernel="identity", metric="cosine")["similarity"]
    near = O.similarity_overlap(a, b, kernel="exp", metric="cosine")["similarity"]
    assert exact == pytest.approx(0.0, abs=1e-15)             # no shared clonotype
    assert near > exact


def test_vj_gating_blocks_cross_v_pair():
    """The exp/step kernels are block-diagonalised by the non-cdr3 key fields.

    Same CDR3, different V: with key=(cdr3_aa,v_call,j_call) the exp kernel must NOT
    connect them — it is 0, matching both the identity kernel and the exact frequency
    overlap F (identity is the exact τ→0 limit of exp on the same spine). Dropping V/J
    from the key (cdr3-only) removes the gate and the same CDR3 gives similarity 1.
    """
    key = (S.CDR3_AA, S.V_CALL, S.J_CALL)
    a = _sample(["CASSLGYEQYF"], [10], v=["TRBV7-9"], j=["TRBJ2-1"])
    b = _sample(["CASSLGYEQYF"], [10], v=["TRBV20-1"], j=["TRBJ2-1"])

    gated = O.similarity_overlap(a, b, key=key, kernel="exp", metric="cosine")["similarity"]
    ident = O.similarity_overlap(a, b, key=key, kernel="identity",
                                 metric="cosine")["similarity"]
    _, m = O.overlap_pair(a, b, key=key)
    assert gated == 0.0                       # exp does not cross the V block
    assert ident == 0.0                       # identity matches the full key
    assert m["F"] == 0.0                       # exact frequency overlap

    # step is gated the same way; and without V/J in the key the same CDR3 is a full match.
    step = O.similarity_overlap(a, b, key=key, kernel="step", max_penalty=1,
                                metric="cosine")["similarity"]
    ungated = O.similarity_overlap(a, b, key=(S.CDR3_AA,), kernel="exp",
                                   metric="cosine")["similarity"]
    assert step == 0.0
    assert ungated == pytest.approx(1.0, abs=1e-12)


@pytest.mark.parametrize("kernel", ["identity", "step", "exp"])
def test_dense_sparse_agree_on_length_differing_fixture(kernel):
    """The sparse path is a genuine sparsification of the SAME gap-block model as dense.

    Includes a length-differing pair (CASSLGYEQYF / CASSLGGYEQYF, a single central
    insertion). Because the sparse path now re-scores every candidate with the identical
    gap-block model (same matrix / gap_open / central gap_prior) instead of using the
    prior-free trie score, dense and sparse agree to machine precision on all three
    kernels. This also exercises the otherwise-untested sparse path (coo/csr assembly,
    setdiag, sparse ``_quad``).
    """
    a = _sample(["CASSLGYEQYF", "CASSLGGYEQYF"], [3, 2])
    b = _sample(["CASSLGYEQYF", "CASSLGGYEQYF"], [1, 4])
    kw = dict(kernel=kernel, metric="cosine")
    if kernel == "step":
        kw["max_penalty"] = 1
    dense = O.similarity_overlap(a, b, dense=True, **kw)["similarity"]
    sparse = O.similarity_overlap(a, b, dense=False, **kw)["similarity"]
    assert sparse == pytest.approx(dense, abs=1e-9)


def test_presence_weight_distinct_value():
    """weight='presence' (TINA-unweighted) uses uniform 1/N weights, NOT frequency.

    For the identity kernel the two samples share exactly one of three clonotypes, so the
    uniform-weight cosine is exactly 1/3 — a value the frequency weighting (0.5524) never
    produces. Also pins self-sim==1 and the [0,1] bound (which held for 'freq' too).
    """
    a = _sample(_A_CDR3, _A_CNT)
    b = _sample(_B_CDR3, _B_CNT)
    pres = O.similarity_overlap(a, b, kernel="identity", weight="presence",
                                metric="cosine")["similarity"]
    freq = O.similarity_overlap(a, b, kernel="identity", weight="freq",
                                metric="cosine")["similarity"]
    assert pres == pytest.approx(1 / 3, abs=1e-12)          # distinguishing pinned value
    assert freq == pytest.approx(0.5523681766661075, abs=1e-12)
    assert abs(pres - freq) > 1e-6                          # presence != freq

    aa = O.similarity_overlap(a, a, kernel="exp", weight="presence", metric="cosine")
    ab = O.similarity_overlap(a, b, kernel="exp", weight="presence", metric="cosine")
    assert 0.0 <= ab["similarity"] <= 1.0
    assert aa["similarity"] == pytest.approx(1.0, rel=1e-12)


def test_nonstandard_aa_clonotype_dropped_and_renormalized():
    """A CDR3 with a non-standard symbol is dropped, and the survivors are renormalized.

    Adding one extra clonotype whose CDR3 carries an ``X``/digit must leave the overlap
    byte-identical to the sample without it: the row is filtered out *and* the remaining
    frequencies are renormalized to sum to 1 (Morisita is scale-sensitive, so a missing
    renormalization would change the value).
    """
    a = _sample(_A_CDR3, _A_CNT)
    b = _sample(_B_CDR3, _B_CNT)
    a_extra = _sample(_A_CDR3 + ["CASX1F"], _A_CNT + [99])  # extra non-standard clonotype

    base = O.similarity_overlap(a, b, kernel="identity", weight="freq", metric="morisita")
    withx = O.similarity_overlap(a_extra, b, kernel="identity", weight="freq",
                                 metric="morisita")
    assert withx["similarity"] == base["similarity"]        # byte-identical (dropped)
    assert withx["pTZq"] == base["pTZq"]                     # and survivors renormalized


def test_all_nonstandard_sample_is_zero_overlap():
    """An all-non-standard sample overlapped with a normal one is empty, not an error.

    Every CDR3 is filtered out, so both weight vectors on that side are size 0; the
    bilinear forms and denominators must degrade to similarity 0 / distance 1 without
    raising (covers the ``_quad`` size-0, ``denom > 0`` and ``total or 1`` guards).
    """
    bad = _sample(["XXX111", "1234"], [5, 5])
    good = _sample(["CASSLGQAYEQYF", "CASSPGTEAFF"], [10, 10])
    for kernel in ("identity", "exp"):
        res = O.similarity_overlap(bad, good, kernel=kernel, metric="cosine")
        assert res["similarity"] == 0.0
        assert res["distance"] == 1.0


def test_cluster_integration_similarity_cosine():
    """pairwise_distances(metric='similarity_cosine') is symmetric, zero-diagonal, and
    orders the toy samples like exact F (a closer to b than to c)."""
    pytest.importorskip("sklearn")
    a = _sample(["CASSA", "CASSB", "CASSC", "CASSD"], [10, 10, 10, 10])
    b = _sample(["CASSA", "CASSB", "CASSC", "CASSE"], [10, 10, 10, 10])
    c = _sample(["CASSX", "CASSY", "CASSZ", "CASSA"], [10, 10, 10, 10])
    toy = {"a": a, "b": b, "c": c}

    dm = O.pairwise_distances(toy, metric="similarity_cosine")
    names = dm["sample"].to_list()
    M = dm.select(names).to_numpy()
    assert np.allclose(M, M.T)
    assert np.allclose(np.diag(M), 0.0)
    ai, bi, ci = names.index("a"), names.index("b"), names.index("c")
    assert M[ai, bi] < M[ai, ci] and M[ai, bi] < M[bi, ci]

    # Morisita variant dispatches too — and orders the samples the same way.
    dm2 = O.pairwise_distances(toy, metric="similarity_morisita")
    names2 = dm2["sample"].to_list()
    M2 = dm2.select(names2).to_numpy()
    assert np.allclose(M2, M2.T) and np.allclose(np.diag(M2), 0.0)
    a2, b2, c2 = names2.index("a"), names2.index("b"), names2.index("c")
    assert M2[a2, b2] < M2[a2, c2] and M2[a2, b2] < M2[b2, c2]


@pytest.mark.slow
def test_tina_reference_arithmetic():
    """TINA reference identity aᵀCb / √(aᵀCa·bᵀCb) on a tiny hand-built C, reproduced by
    the module's bilinear machinery with an identity-key set and a supplied 2×2 kernel."""
    # Hand-built similarity C and abundance vectors.
    C = np.array([[1.0, 0.5], [0.5, 1.0]])
    a_vec = np.array([2.0, 1.0])
    b_vec = np.array([1.0, 3.0])
    tina = (a_vec @ C @ b_vec) / np.sqrt((a_vec @ C @ a_vec) * (b_vec @ C @ b_vec))
    assert tina == pytest.approx(0.8910421112136306, abs=1e-12)

    # Same arithmetic through similarity_overlap: a CDR3 pair whose BLOSUM62 penalty gives
    # Z_off = 0.5 = exp(-P/τ) ⇒ P/τ = ln 2. A↔D penalty is 14; τ = 14/ln2 ⇒ Z_off = 0.5.
    s1, s2 = "CASSAEQYF", "CASSDEQYF"
    tau = 14.0 / np.log(2.0)
    a = _sample([s1, s2], [2, 1])            # relative abundances ∝ a_vec
    b = _sample([s1, s2], [1, 3])            # ∝ b_vec
    res = O.similarity_overlap(a, b, kernel="exp", tau=tau, metric="cosine")
    # cosine is scale-invariant, so normalised abundances give the same value as a_vec/b_vec.
    assert res["similarity"] == pytest.approx(tina, abs=1e-9)
