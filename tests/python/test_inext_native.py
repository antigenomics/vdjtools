"""Tests for the native (C++/pybind11) iNEXT kernel in ``vdjtools._core``.

The native size-based engine must reproduce the validated numpy reference in
``vdjtools.stats.inext`` (which is itself pinned to iNEXT R 3.0.2 in
``test_inext.py``). Point curves are deterministic and checked at ``rtol=1e-9``;
the bootstrap is stochastic so only structure + a seeded cross-check against the
numpy reference are asserted. The batch API is checked for one block per sample
and point-estimate parity with per-sample :func:`inext`.
"""
import importlib
from pathlib import Path

import numpy as np
import polars as pl
import pytest

import vdjtools._core as core
from vdjtools import stats

# The private numpy reference module (the package attribute ``inext`` is the
# re-exported function, so reach the module via importlib).
ref = importlib.import_module("vdjtools.stats.inext")

ASSETS = Path(__file__).resolve().parent.parent / "assets"


def _load(name: str) -> np.ndarray:
    return np.array([int(t) for t in (ASSETS / name).read_text().split()],
                    dtype=np.float64)


@pytest.fixture(scope="module")
def girdled() -> np.ndarray:
    return _load("spider_girdled.txt")


@pytest.fixture(scope="module")
def logged() -> np.ndarray:
    return _load("spider_logged.txt")


def _default_sizes(x: np.ndarray, knots: int = 40) -> list[float]:
    n = int(x.sum())
    grid = np.floor(np.linspace(1, 2 * n, knots)).astype(np.int64)
    grid = np.unique(np.concatenate([grid, np.array([n], dtype=np.int64)]))
    return [float(m) for m in grid if m >= 1]


# --------------------------------------------------------------------------- #
# native symbols are present
# --------------------------------------------------------------------------- #
def test_core_exposes_inext_symbols():
    for name in ("inext_curve", "inext_bootstrap", "inext_batch",
                 "inext_digamma", "InextCurve", "InextSample"):
        assert hasattr(core, name), name


# --------------------------------------------------------------------------- #
# digamma == scipy
# --------------------------------------------------------------------------- #
def test_digamma_matches_scipy():
    from scipy.special import digamma
    for x in [1.0, 2.0, 3.5, 6.0, 26.0, 168.0, 252.0, 5000.0]:
        assert np.isclose(core.inext_digamma(x), float(digamma(x)), rtol=1e-10, atol=1e-11)


# --------------------------------------------------------------------------- #
# native point curve == numpy reference (deterministic, rtol 1e-9)
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("name", ["spider_girdled.txt", "spider_logged.txt"])
def test_native_curve_matches_python(name):
    x = _load(name)
    qs = [0, 1, 2]
    sizes = _default_sizes(x)
    cur = core.inext_curve(x.tolist(), qs, sizes)
    qD_native = np.asarray(cur.qD)
    cov_native = np.asarray(cur.coverage)

    qD_ref = np.array([[ref._diversity_at(x, m, q) for m in sizes] for q in qs])
    cov_ref = np.array([min(max(ref._chat(x, m), 0.0), 1.0) for m in sizes])

    assert np.allclose(qD_native, qD_ref, rtol=1e-9, atol=1e-9)
    assert np.allclose(cov_native, cov_ref, rtol=1e-9, atol=1e-12)


def test_native_curve_matches_inext_public(girdled):
    # native curve must also equal the public size-based inext() point estimates
    sizes = [100.0, 168.0, 200.0, 336.0]
    cur = core.inext_curve(girdled.tolist(), [0, 1, 2], sizes)
    df = stats.inext(girdled, q=(0, 1, 2), sizes=[int(m) for m in sizes], se=False)
    got = {(r["order_q"], r["m"]): r["qD"] for r in df.iter_rows(named=True)}
    for i, q in enumerate([0, 1, 2]):
        for j, m in enumerate(sizes):
            assert np.isclose(cur.qD[i][j], got[(q, int(m))], rtol=1e-9)


# --------------------------------------------------------------------------- #
# native bootstrap: structure + seeded cross-check vs numpy reference
# --------------------------------------------------------------------------- #
def test_native_bootstrap_structure(girdled):
    qs = [0, 1, 2]
    sizes = [50.0, 100.0, 168.0, 250.0]
    se = np.asarray(core.inext_bootstrap(girdled.tolist(), qs, sizes, 50, 0))
    assert se.shape == (len(qs), len(sizes))
    assert np.all(se >= 0.0)
    assert np.all(np.isfinite(se))


def test_native_bootstrap_deterministic_for_seed(girdled):
    # same seed -> identical SE; different seed -> (almost surely) different
    a = np.asarray(core.inext_bootstrap(girdled.tolist(), [0], [100.0], 40, 7))
    b = np.asarray(core.inext_bootstrap(girdled.tolist(), [0], [100.0], 40, 7))
    c = np.asarray(core.inext_bootstrap(girdled.tolist(), [0], [100.0], 40, 8))
    assert np.array_equal(a, b)
    assert not np.array_equal(a, c)


def test_native_bootstrap_se_close_to_python(girdled):
    # Both bootstraps are stochastic with independent RNGs; with the same seed and
    # a decent nboot they agree in expectation. Deterministic at seed=0/nboot=400.
    qs = [0, 1, 2]
    sizes = [50.0, 100.0, 168.0, 250.0]
    se_native = np.asarray(core.inext_bootstrap(girdled.tolist(), qs, sizes, 400, 0))
    se_python = ref._bootstrap_se(girdled, sizes, qs, 400, 0)
    rel = np.abs(se_native - se_python) / (np.abs(se_python) + 1e-30)
    assert rel.max() < 0.15, f"max rel SE diff {rel.max():.3f}"


# --------------------------------------------------------------------------- #
# batch API
# --------------------------------------------------------------------------- #
def test_inext_batch_blocks_and_parity(girdled, logged):
    samples = [girdled, logged, girdled]
    df = stats.inext_batch(samples, q=(0, 1, 2), knots=12, se=True, nboot=20, seed=0)
    assert df.columns == ["sample", "order_q", "m", "method", "sample_coverage",
                          "qD", "qD_lo", "qD_hi"]
    # one block per sample
    assert sorted(df["sample"].unique().to_list()) == [0, 1, 2]

    # point estimates + coverage + method match per-sample inext()
    for idx, x in enumerate(samples):
        block = df.filter(pl.col("sample") == idx).sort(["order_q", "m"])
        single = stats.inext(x, q=(0, 1, 2), knots=12, se=False).sort(["order_q", "m"])
        assert np.allclose(block["qD"].to_numpy(), single["qD"].to_numpy(), rtol=1e-9)
        assert np.allclose(block["sample_coverage"].to_numpy(),
                           single["sample_coverage"].to_numpy(), rtol=1e-9, atol=1e-12)
        assert block["method"].to_list() == single["method"].to_list()


def test_inext_batch_ci_structure(girdled, logged):
    df = stats.inext_batch([girdled, logged], q=(0, 1, 2), knots=10,
                           se=True, nboot=30, seed=0)
    assert (df["qD_lo"] <= df["qD"]).all()
    assert (df["qD"] <= df["qD_hi"]).all()
    assert (df["qD_lo"] >= 0.0).all()
    assert ((df["sample_coverage"] >= 0.0) & (df["sample_coverage"] <= 1.0)).all()


def test_inext_batch_se_false_nulls(girdled, logged):
    df = stats.inext_batch([girdled, logged], q=0, knots=8, se=False)
    assert df["qD_lo"].null_count() == df.height
    assert df["qD_hi"].null_count() == df.height


def test_inext_batch_dataframe_input(girdled, logged):
    frame = pl.concat([
        pl.DataFrame({"duplicate_count": girdled.astype(np.int64),
                      "sample_id": ["A"] * girdled.size}),
        pl.DataFrame({"duplicate_count": logged.astype(np.int64),
                      "sample_id": ["B"] * logged.size}),
    ])
    df = stats.inext_batch(frame, q=(0, 1, 2), knots=10, se=False)
    assert set(df["sample"].unique().to_list()) == {"A", "B"}
    # sample "A" (girdled) matches the vector input
    a = df.filter((pl.col("sample") == "A") & (pl.col("order_q") == 0)).sort("m")
    ref_df = stats.inext(girdled, q=0, knots=10, se=False).sort("m")
    assert np.allclose(a["qD"].to_numpy(), ref_df["qD"].to_numpy(), rtol=1e-9)


def test_inext_batch_requires_sample_id():
    frame = pl.DataFrame({"duplicate_count": [1, 2, 3]})
    with pytest.raises(ValueError, match="sample_id"):
        stats.inext_batch(frame)


def test_rarefaction_batch_alias():
    assert stats.rarefaction_batch is stats.inext_batch
