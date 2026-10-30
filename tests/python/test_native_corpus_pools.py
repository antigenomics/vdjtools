"""Learned corpus pools use one native batch and preserve CPU-invariant draws."""
import pytest
from polars.testing import assert_frame_equal

from vdjtools.signature import corpus as C


def test_learned_pool_is_one_model_batch_and_identical_across_cpu_counts(tmp_path, monkeypatch):
    import vdjtools.model as models
    from vdjtools import cores

    original = models.load_bundled
    fixture = original("TRG", "olga")
    calls = []

    def load(locus, source):
        calls.append((locus, source))
        return fixture

    monkeypatch.setattr(models, "load_bundled", load)
    monkeypatch.setattr(cores, "available_cores", lambda: 2)
    monkeypatch.setattr(C, "_map_unordered", lambda *args: pytest.fail("learned generation must not spawn chunk workers"))
    got = []
    for jobs in (1, 2):
        scratch = tmp_path / str(jobs)
        scratch.mkdir()
        paths = C.build_pools(["TRG"], size=30, n_samples=4, seed=7,
                              source="learned", tmp=scratch, n_jobs=jobs)
        assert len(paths["TRG"]) == 1
        frame = C.read_pools(paths)["TRG"]
        assert len(frame) == C.pool_target("TRG", 30, 4)
        assert frame["productive"].all()
        got.append(frame)
    assert calls == [("TRG", "isalgo"), ("TRG", "isalgo")]
    assert_frame_equal(*got)
    with pytest.raises(ValueError, match="invalid n_jobs"):
        C.build_pools(["TRG"], size=30, n_samples=4, seed=7,
                      source="learned", tmp=tmp_path, n_jobs=3)
    assert len(calls) == 2, "reject CPU oversubscription before loading a model"


def test_explicit_olga_pool_retains_reference_chunk_seeds(tmp_path):
    import polars as pl
    from vdjtools.model import load_bundled
    from vdjtools.model.generate import generate

    total = C.pool_target("TRG", 20, 2)
    expected = [generate(load_bundled("TRG", "olga"), b-a, seed=7+k,
                         productive_only=True)
                for k, (a, b) in enumerate(C._ranges(total, C.POOL_CHUNKS))]
    paths = C.build_pools(["TRG"], size=20, n_samples=2, seed=7,
                          source="olga", tmp=tmp_path, n_jobs=1)
    assert_frame_equal(C.read_pools(paths)["TRG"], pl.concat(expected, rechunk=False))
