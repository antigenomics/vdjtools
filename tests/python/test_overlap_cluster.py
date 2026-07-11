"""Tests for vdjtools.overlap.cluster — pairwise distances + MDS/hclust.

Toy samples pin the distance-matrix invariants (symmetry, zero diagonal) and the MDS
embedding shape. A guarded real-data test on the cached aging cohort reproduces the
notebook's "repertoires drift apart with age" divergence via the reusable API.
sklearn (the ``overlap`` extra) is guarded with ``importorskip``.
"""
import warnings
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from vdjtools.io import schema as S
from vdjtools import overlap as O

pytest.importorskip("sklearn")

_AGING_DIR = Path("examples/.data/aging")


def _sample(cdr3, counts):
    n = len(cdr3)
    df = pl.DataFrame({S.V_CALL: ["TRBV1"] * n, S.J_CALL: ["TRBJ1"] * n,
                       S.CDR3_AA: cdr3, S.COUNT: counts})
    return S.add_locus(S.normalize(df, recompute_freq=True))


def _toy_samples():
    a = _sample(["CASSA", "CASSB", "CASSC", "CASSD"], [10, 10, 10, 10])
    b = _sample(["CASSA", "CASSB", "CASSC", "CASSE"], [10, 10, 10, 10])
    c = _sample(["CASSX", "CASSY", "CASSZ", "CASSA"], [10, 10, 10, 10])
    return {"a": a, "b": b, "c": c}


def test_pairwise_distances_symmetric_zero_diag():
    dm = O.pairwise_distances(_toy_samples(), metric="F")
    names = dm["sample"].to_list()
    M = dm.select(names).to_numpy()
    assert np.allclose(M, M.T)                       # symmetric
    assert np.allclose(np.diag(M), 0.0)              # zero diagonal
    # a and b share 3/4, a and c share 1/4 -> a is closer to b than to c.
    ai, bi, ci = names.index("a"), names.index("b"), names.index("c")
    assert M[ai, bi] < M[ai, ci]


def test_pairwise_distances_long_form():
    dm = O.pairwise_distances(_toy_samples(), metric="F", form="long")
    assert dm.columns == ["sample_a", "sample_b", "distance"]
    assert dm.height == 9                             # 3x3
    self_rows = dm.filter(pl.col("sample_a") == pl.col("sample_b"))
    assert (self_rows["distance"] == 0.0).all()


def test_pairwise_distances_jaccard_transform():
    # jaccard distance = 1 - d12/(d1+d2-d12); a,b share 3 of 5 union -> 1 - 3/5 = 0.4
    dm = O.pairwise_distances(_toy_samples(), metric="jaccard")
    names = dm["sample"].to_list()
    M = dm.select(names).to_numpy()
    ai, bi = names.index("a"), names.index("b")
    assert np.isclose(M[ai, bi], 0.4)


def test_cluster_samples_mds_returns_2d():
    dm = O.pairwise_distances(_toy_samples(), metric="F")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")              # sklearn MDS FutureWarnings
        coords = O.cluster_samples(dm, method="mds", n_components=2)
    assert coords.columns == ["sample", "mds1", "mds2"]
    assert coords.height == 3
    assert set(coords["sample"].to_list()) == {"a", "b", "c"}


def test_cluster_samples_hclust_labels():
    dm = O.pairwise_distances(_toy_samples(), metric="F")
    out = O.cluster_samples(dm, method="hclust", n_components=2)
    assert set(out.columns) == {"sample", "cluster", "leaf_order"}
    # a and b (share 3/4) should land in the same flat cluster; c in the other.
    lab = dict(zip(out["sample"].to_list(), out["cluster"].to_list()))
    assert lab["a"] == lab["b"] and lab["c"] != lab["a"]


def test_cluster_samples_metadata_join():
    dm = O.pairwise_distances(_toy_samples(), metric="F")
    meta = pl.DataFrame({"sample": ["a", "b", "c"], "group": ["x", "x", "y"]})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        coords = O.cluster_samples(dm, method="mds", metadata=meta)
    assert "group" in coords.columns


def _load_aging(n_samples=20, depth_cap=50000, seed=0):
    """Load n age-spanning cached aging samples, downsampled to a common depth."""
    from vdjtools import io as vio

    meta = (vio.read_metadata(_AGING_DIR / "metadata_aging.txt")
            .with_columns(pl.col("age").cast(pl.Int64)))
    meta = meta.filter(pl.col("sample_id").map_elements(
        lambda s: (_AGING_DIR / f"{s}.txt.gz").exists(), return_dtype=pl.Boolean))
    if meta.height < n_samples:
        return None, None
    meta = meta.sort("age")
    idx = np.linspace(0, meta.height - 1, n_samples).round().astype(int)
    sel = meta[idx.tolist()]

    key = (S.CDR3_AA, S.V_CALL, S.J_CALL)
    rng = np.random.default_rng(seed)
    agg, reads = {}, {}
    for sid in sel["sample_id"].to_list():
        g = (vio.read(_AGING_DIR / f"{sid}.txt.gz", fmt="vdjtools")
             .group_by(list(key), maintain_order=True)
             .agg(pl.col(S.COUNT).sum().alias(S.COUNT)))
        agg[sid], reads[sid] = g, int(g[S.COUNT].sum())
    depth = min(depth_cap, min(reads.values()))
    ds = {}
    for sid, g in agg.items():
        p = g[S.COUNT].to_numpy().astype(float)
        p = p / p.sum()
        drawn = rng.multinomial(depth, p)
        keep = drawn > 0
        ds[sid] = (g.filter(pl.Series(keep))
                   .with_columns(pl.Series(S.COUNT, drawn[keep]))
                   .with_columns((pl.col(S.COUNT) / depth).alias(S.FREQ)))
    return ds, sel


@pytest.mark.skipif(not _AGING_DIR.exists(), reason="aging cohort not cached")
def test_aging_divergence_correlates_with_age():
    """Reproduce the aging notebook's divergence result through the public API:
    MDS distance-from-centroid of the pairwise-F distance matrix increases with age."""
    pytest.importorskip("scipy")
    from scipy.stats import spearmanr

    ds, sel = _load_aging(n_samples=20)
    if ds is None:
        pytest.skip("insufficient cached aging samples")

    key = (S.CDR3_AA, S.V_CALL, S.J_CALL)
    dm = O.pairwise_distances(ds, metric="F", key=key)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        coords = O.cluster_samples(dm, method="mds", n_components=2)
    xy = coords.select("mds1", "mds2").to_numpy()
    dist_centroid = np.linalg.norm(xy - xy.mean(axis=0), axis=1)

    # Align ages to the coords' sample order.
    age_by_sample = dict(zip(sel["sample_id"].to_list(), sel["age"].to_list()))
    ages = np.array([age_by_sample[s] for s in coords["sample"].to_list()])
    r, p = spearmanr(ages, dist_centroid)
    assert r > 0.4 and p < 0.05
