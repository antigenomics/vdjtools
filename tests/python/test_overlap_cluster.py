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
                       S.JUNCTION_AA: cdr3, S.COUNT: counts})
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
    # The embedding must preserve the distance structure (a,b share 3/4; a,c share
    # 1/4), else a degenerate MDS that ignored the matrix would still pass the shape.
    xy = dict(zip(coords["sample"].to_list(), coords.select("mds1", "mds2").to_numpy()))
    d = lambda x, y: float(np.linalg.norm(xy[x] - xy[y]))
    assert d("a", "b") < d("a", "c") and d("a", "b") < d("b", "c")


def test_pairwise_distances_fuzzy_scope_path():
    """The scope-driven fuzzy path (delegating to fuzzy_overlap_metrics) is exercised
    and rejects metrics other than F."""
    pytest.importorskip("vdjmatch")
    dm = O.pairwise_distances(_toy_samples(), metric="F", scope="1,0,0,1")
    names = dm["sample"].to_list()
    M = dm.select(names).to_numpy()
    assert np.allclose(M, M.T) and np.allclose(np.diag(M), 0.0)
    assert np.isfinite(M[np.triu_indices(len(names), 1)]).all()
    # Only F is supported with a fuzzy scope.
    with pytest.raises(ValueError):
        O.pairwise_distances(_toy_samples(), metric="jaccard", scope="1,0,0,1")


def test_pairwise_distances_single_sample():
    """A single sample yields a 1x1 zero-diagonal matrix (edge case)."""
    dm = O.pairwise_distances({"only": _toy_samples()["a"]}, metric="F")
    assert dm.height == 1 and dm.select("only").to_numpy()[0, 0] == 0.0


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

    key = (S.JUNCTION_AA, S.V_CALL, S.J_CALL)
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

    key = (S.JUNCTION_AA, S.V_CALL, S.J_CALL)
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


@pytest.mark.parametrize('metric', ['F', 'F2', 'D', 'jaccard'])
def test_exact_matrix_batches_keys_instead_of_joining_every_pair(monkeypatch, metric):
    from vdjtools.overlap import cluster as C
    samples = _toy_samples()
    samples['duplicates'] = pl.concat([samples['a'], samples['a'].head(1)])
    samples['empty'] = samples['a'].clear()
    samples['null'] = samples['a'].head(1).with_columns(pl.lit(None, pl.String).alias(S.V_CALL))
    samples['zero'] = samples['a'].with_columns(pl.lit(0, pl.Int64).alias(S.COUNT))
    names = list(samples)
    aggs = [C._aggregate(df, list(C.DEFAULT_KEY)) for df in samples.values()]
    expected = np.zeros((len(names), len(names)))
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            _, values = C._overlap_from_agg(aggs[i], aggs[j], list(C.DEFAULT_KEY))
            expected[i, j] = expected[j, i] = C._TRANSFORM[metric](C._exact_similarity(values, metric))

    def no_pair_join(*args, **kwargs):
        raise AssertionError('exact matrix must batch the keys, not join every sample pair')

    monkeypatch.setattr(C, '_overlap_from_agg', no_pair_join)
    got = C.pairwise_distances(samples, metric=metric)
    assert got['sample'].to_list() == names
    np.testing.assert_allclose(got.select(names).to_numpy(), expected, rtol=1e-14, atol=1e-14)
    assert got.equals(C.pairwise_distances(samples, metric=metric))


@pytest.mark.parametrize('metric', ['F', 'F2', 'D', 'jaccard', 'R'])
@pytest.mark.parametrize('case', ['empty', 'null', 'signed', 'custom-key', 'large-counts', 'fractional-counts', 'int32-counts'])
def test_exact_matrix_boundary_cases_preserve_pairwise_values(metric, case):
    from vdjtools.overlap import cluster as C
    a, b = _toy_samples()['a'], _toy_samples()['b']
    key = list(C.DEFAULT_KEY)
    if case == 'empty':
        a, b = a.clear(), b.clear()
    elif case == 'null':
        a = a.with_columns(pl.lit(None, pl.String).alias(S.V_CALL))
        b = b.with_columns(pl.lit(None, pl.String).alias(S.V_CALL))
    elif case == 'signed':
        a = a.with_columns(pl.Series(S.COUNT, [-1, 2, 3, 4]))
    elif case == 'large-counts':
        a = a.with_columns(pl.Series(S.COUNT, [2**53 + 1, 2, 3, 4]))
    elif case == 'fractional-counts':
        a = a.with_columns((pl.col(S.COUNT).cast(pl.Float64) / 2).alias(S.COUNT))
    elif case == 'int32-counts':
        a = a.with_columns(pl.col(S.COUNT).cast(pl.Int32))
    else:
        key = [S.JUNCTION_AA]
        b = b.with_columns(pl.lit(None, pl.String).alias(S.V_CALL))
    with np.errstate(invalid='ignore'):
        values = O.overlap_metrics(a, b, key=tuple(key))
        expected = C._TRANSFORM[metric](C._exact_similarity(values, metric))
        samples = {'a': a, 'b': b, 'c': a.clear(), 'd': b.clear()}
        got = C.pairwise_distances(samples, metric=metric, key=tuple(key), form='long')
    assert got.select('sample_a', 'sample_b').rows() == [(a, b) for a in samples for b in samples]
    selected = got.filter(pl.col('sample_a').is_in(['a', 'b']) & pl.col('sample_b').is_in(['a', 'b']))
    np.testing.assert_allclose(selected['distance'].to_numpy(), [0, expected, expected, 0],
                               rtol=1e-14, atol=1e-14, equal_nan=True)


def test_pairwise_matrix_rejects_output_larger_than_allocation(monkeypatch):
    monkeypatch.setattr('vdjtools.cores.available_memory', lambda: 64)
    with pytest.raises(MemoryError, match='reduce sample count'):
        O.pairwise_distances(_toy_samples())


@pytest.mark.parametrize('n', [0, 1])
def test_empty_or_single_sample_output_forms(n):
    samples = dict(list(_toy_samples().items())[:n])
    assert O.pairwise_distances(samples).shape == (n, n + 1)
    long = O.pairwise_distances(samples, form='long')
    assert long.shape == (n * n, 3)
    if n:
        assert long['distance'].to_list() == [0.0]


@pytest.mark.parametrize('metric', ['similarity_cosine', 'similarity_morisita'])
def test_delegated_distance_output_preserves_similarity(metric):
    from vdjtools.overlap import cluster as C
    a = _sample(['CASSAF', 'CASSCF', 'CASSDF'], [1, 2, 3])
    b = _sample(['CASSAF', 'CASSCF', 'CASSEF'], [3, 2, 1])
    similarity = C._delegated_similarity(a, b, metric, C.DEFAULT_KEY, None)
    got = C.pairwise_distances({'a': a, 'b': b}, metric=metric, form='long')
    assert got['distance'].to_list() == [0.0, 1 - similarity, 1 - similarity, 0.0]


@pytest.mark.parametrize('metric', ['F', 'F2'])
def test_many_shared_keys_do_not_accumulate_normalized_frequency_roundoff(metric):
    from vdjtools.overlap import cluster as C
    depth = 20000
    samples = {str(i): _sample([f'CASS{k}F' for k in range(i * 10000, i * 10000 + depth)],
                               np.ones(depth, dtype=np.int64)) for i in range(4)}
    values = O.overlap_metrics(samples['0'], samples['1'])
    expected = C._TRANSFORM[metric](values[metric])
    actual = C.pairwise_distances(samples, metric=metric)['1'][0]
    np.testing.assert_allclose(actual, expected, rtol=1e-14, atol=0)
