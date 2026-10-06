"""Resource budgets and annotation fields are observable API contracts."""
import polars as pl
import pytest


def test_annotation_preserves_stage_one_fields():
    from arda.cdr3fix import markup_batch
    from vdjtools.model import annotate_junctions
    aa, v, j = ["CAVTDDKIIF"], ["TRAV12-2*01"], ["TRAJ30*01"]
    want = markup_batch(pl.DataFrame({"cdr3": aa, "v": v, "j": j, "species": ["human"]}),
                        cdr3="cdr3", v="v", j="j", species="species")
    got = annotate_junctions(aa, v, j, threads=1)
    for column in ("v_fix", "j_fix", "v_canonical", "j_canonical"):
        assert got[column].to_list() == want[column].to_list()
    with pytest.raises(ValueError, match="model_source"):
        annotate_junctions(aa, v, j, model_source="typo")


def test_infer_native_honours_thread_budget(monkeypatch):
    from vdjtools import _core
    from vdjtools.model import load_bundled
    from vdjtools.model.infer import infer_native
    real = _core.estep_batch
    seen = []
    def record(*args):
        seen.append(args[6])
        return real(*args)
    monkeypatch.setattr(_core, "estep_batch", record)
    m = load_bundled("TRA")
    infer_native(m, ["TGTGCTGTGACAGATGATAAAATCATCTTT"], max_iter=1, threads=2)
    assert seen == [2]
    with pytest.raises(ValueError, match="threads"):
        infer_native(m, [], threads=-1)


def test_build_all_uses_one_kernel_per_build(monkeypatch):
    from vdjtools.model import data
    seen = []
    def fake(chain, **kw):
        seen.append(kw["threads"])
        return None, None, {}
    monkeypatch.setattr(data, "build_model", fake)
    assert len(data.build_all(["TRA", "TRB"], workers=2)) == 2
    assert seen == [1, 1]
    seen.clear()
    data.build_all(["TRA"], workers=1, threads=3)
    assert seen == [3]
    with pytest.raises(ValueError, match="workers"):
        data.build_all(workers=-1)


def test_corpus_rejects_negative_budget():
    from vdjtools.signature.corpus import _workers
    with pytest.raises(ValueError, match="n_jobs"):
        _workers(-1, 4)


def test_qc_does_not_hide_broken_reference(monkeypatch):
    from vdjtools.model import reference
    from vdjtools.signature.features import qc_channel
    def broken(*args):
        raise RuntimeError("broken artifact")
    monkeypatch.setattr(reference, "load_germline", broken)
    with pytest.raises(RuntimeError, match="broken artifact"):
        qc_channel(pl.DataFrame(), pl.DataFrame(), "TRB", 0.0)


def test_build_model_forwards_threads_to_annotation_and_em(monkeypatch):
    from vdjtools.model import data, load_bundled
    import importlib
    infer = importlib.import_module("vdjtools.model.infer")
    from vdjtools.model.infer import InferenceReport
    model = load_bundled("TRA")
    seen = []
    frame = pl.DataFrame({"junction": ["TGTGCTGTGACAGATGATAAAATCATCTTT"]})
    def prepare(*args, **kw):
        seen.append(("annotation", kw["threads"]))
        return frame
    def fit(*args, **kw):
        seen.append(("em", kw["threads"]))
        return model, InferenceReport()
    monkeypatch.setattr(data, "prepare", prepare)
    monkeypatch.setattr(data, "_filter_for_em", lambda *a: frame)
    monkeypatch.setattr(data, "_build_masks", lambda *a: None)
    monkeypatch.setattr(infer, "infer_native", fit)
    data.build_model("TRA", template=model, threads=3)
    assert seen == [("annotation", 3), ("em", 3)]


def test_native_pack_does_not_retain_an_unbounded_model_history():
    import gc
    import weakref
    from vdjtools.model import load_bundled
    from vdjtools.model.model import Model
    from vdjtools.model import native
    native._pack_cache.clear()
    base = load_bundled("TRA")
    class WatchedModel(Model):
        pass
    refs = []
    for _ in range(20):
        model = WatchedModel(manifest=base.manifest, tables=base.tables, genomic=base.genomic)
        refs.append(weakref.ref(model))
        native.pack(model)
        assert len(native._pack_cache) <= 8
    del model
    gc.collect()
    assert sum(r() is not None for r in refs) <= 8


def test_model_check_uses_explicit_reference_organism(monkeypatch):
    import importlib
    from typer.testing import CliRunner
    from vdjtools.cli import app
    from vdjtools.model import reference
    check = importlib.import_module("vdjtools.model.check")
    seen = []
    marker = pl.DataFrame({"marker": [1]})
    def germline(locus, organism):
        seen.append((locus, organism))
        return marker
    def audit(model, *, germline):
        assert germline is marker
        return pl.DataFrame(schema={"severity": pl.String})
    from types import SimpleNamespace
    monkeypatch.setattr("vdjtools.cli._model_arg", lambda spec: SimpleNamespace(locus="TRA"))
    monkeypatch.setattr(reference, "load_germline", germline)
    monkeypatch.setattr(check, "check_model", audit)
    result = CliRunner().invoke(app, ["model", "check", "TRA", "--organism", "mouse"])
    assert result.exit_code == 0, result.exception
    assert seen == [("TRA", "mouse")]


def test_unused_comparison_labels_are_deprecated():
    from vdjtools.model import load_bundled
    from vdjtools.model.analyze import compare_models
    model = load_bundled("TRA")
    with pytest.warns(DeprecationWarning, match="labels is unused"):
        out = compare_models(model, model, labels=("first", "second"))
    assert out["tv"].max() == 0
