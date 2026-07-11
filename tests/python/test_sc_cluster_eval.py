"""Clustering-evaluation metric oracle pins (from clustereval's test_scores.py)."""
from __future__ import annotations


import pytest

pytest.importorskip("sklearn")  # cluster_eval builds its contingency via sklearn ([sc] extra)

from vdjtools.sc.cluster_eval import (  # noqa: E402
    assign_singleton_ids,
    cluster_eval,
    homogeneity,
    parsimony,
    purity,
    q_measure,
)

TRUE = [1, 1, 2, 2, 3, 3]


def test_identity_all_ones():
    """pred == true → every metric is 1.0."""
    m = cluster_eval(TRUE, TRUE)
    for key in ("purity", "normalized_purity", "inverse_purity",
                "normalized_inverse_purity", "homogeneity", "parsimony", "q_measure"):
        assert m[key] == pytest.approx(1.0), key


def test_single_cluster():
    """One cluster for all: homogeneity=0, parsimony=1, norm_purity=0, norm_inv_purity=1."""
    m = cluster_eval(TRUE, [1, 1, 1, 1, 1, 1])
    assert m["homogeneity"] == pytest.approx(0.0)
    assert m["parsimony"] == pytest.approx(1.0)
    assert m["normalized_purity"] == pytest.approx(0.0)
    assert m["normalized_inverse_purity"] == pytest.approx(1.0)
    # purity floor = majority-class fraction = 2/6.
    assert m["purity"] == pytest.approx(1.0 / 3.0)
    assert m["q_measure"] == pytest.approx(0.0)  # h <= 0 → q = 0


def test_all_singletons():
    """Every item its own cluster: homogeneity=1, parsimony=0, norm_purity=1, norm_inv_purity=0."""
    m = cluster_eval(TRUE, [1, 2, 3, 4, 5, 6])
    assert m["homogeneity"] == pytest.approx(1.0)
    assert m["parsimony"] == pytest.approx(0.0, abs=1e-12)
    assert m["normalized_purity"] == pytest.approx(1.0)
    assert m["normalized_inverse_purity"] == pytest.approx(0.0, abs=1e-12)
    assert m["purity"] == pytest.approx(1.0)
    assert m["inverse_purity"] == pytest.approx(0.5)  # |unique(true)|/N = 3/6
    assert m["q_measure"] == pytest.approx(0.0, abs=1e-12)  # p <= 0 → q = 0


def test_q_measure_harmonic_mean():
    """q_measure(beta=1) is the harmonic mean of homogeneity and parsimony."""
    pred = [1, 1, 2, 2, 2, 3]
    h = homogeneity(TRUE, pred)
    p = parsimony(TRUE, pred)
    assert 0.0 < h < 1.0 and 0.0 < p < 1.0
    assert q_measure(TRUE, pred) == pytest.approx(2 * h * p / (h + p))


def test_purity_partial():
    """A partially-correct clustering has purity strictly between the trivial bounds."""
    # clusters: {1,1}, {2,2}, {3,3}→ but merge one 3 with a 2 cluster.
    pred = [10, 10, 20, 20, 20, 30]
    pur = purity(TRUE, pred)
    # cluster 20 has {2,2,3} → best=2; others pure. purity = (2+2+1)/6? no:
    # max per cluster: 10→2, 20→2 (two 2s), 30→1 → (2+2+1)/6 = 5/6.
    assert pur == pytest.approx(5.0 / 6.0)


def test_single_true_class_edges():
    """H(C)==0 (single true class) → homogeneity defined as 1.0."""
    assert homogeneity([1, 1, 1, 1], [1, 2, 1, 2]) == pytest.approx(1.0)


def test_assign_singleton_ids():
    """Sentinels become distinct negative ids; real ids pass through."""
    assert assign_singleton_ids([5, None, 5, None]) == [5, -1, 5, -2]
    assert assign_singleton_ids([1, -1, 2, -1], sentinel=-1) == [1, -1, 2, -2]
    # After expansion, an all-sentinel prediction behaves as all-singletons.
    pred = assign_singleton_ids([None, None, None, None, None, None])
    assert parsimony(TRUE, pred) == pytest.approx(0.0, abs=1e-12)


def test_parsimony_pins_lnN_denominator():
    """parsimony's denominator is ``ln N - H(C)``: a partial fragmentation yields an
    independently hand-computed value, so the ln-N basis is pinned through the library
    (this breaks if the denominator formula drifts)."""
    true = [1, 1, 2, 2, 3, 3]
    pred = [10, 10, 20, 20, 30, 40]      # classes 1,2 intact; class 3 split in two
    # n rows [2],[2],[1,1]: H(C)=ln3, H(K|C)=(1/3)ln2, denom=ln6-ln3=ln2
    #   -> parsimony = 1 - (1/3 ln2)/ln2 = 2/3
    assert parsimony(true, pred) == pytest.approx(2.0 / 3.0)
