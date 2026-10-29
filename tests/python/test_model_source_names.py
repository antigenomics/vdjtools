"""Canonical model names never relabel historical assets or silently replace missing fits."""
import importlib
import inspect
import runpy
from pathlib import Path

import pytest

from vdjtools.model import bundled
from vdjtools.signature import corpus


def test_new_default_is_separate_and_missing_fits_fail(tmp_path, monkeypatch, toy_model):
    toy_model.save(tmp_path / "_bundled" / "learned" / "TRB")
    monkeypatch.setattr(bundled, "files", lambda _: tmp_path)
    assert inspect.signature(bundled.load_bundled).parameters["source"].default == "isalgo"
    with pytest.raises(FileNotFoundError, match="isalgo"):
        bundled.load_bundled("TRB")
    for alias in ("learned", "arda"):
        with pytest.raises(FileNotFoundError, match="isalgo"):
            bundled.load_bundled("TRB", alias)
    from dataclasses import replace
    fresh = replace(toy_model, manifest=replace(toy_model.manifest, source="new-isalgo-toy"))
    fresh.save(tmp_path / "_bundled" / "isalgo" / "TRB")
    for source in ("isalgo", "learned", "arda"):
        assert bundled.load_bundled("TRB", source, collapse=False).manifest.source == "new-isalgo-toy"
    assert corpus.model_fingerprint(["TRB"], "isalgo") == corpus.model_fingerprint(["TRB"], "learned")
    assert corpus.model_fingerprint(["TRB"], "isalgo") == corpus.model_fingerprint(["TRB"], "arda")
    toy_model.save(tmp_path / "_bundled" / "isalgo" / "mouse_TRB")
    assert bundled.load_bundled("TRB", organism="mouse", collapse=False).manifest == toy_model.manifest
    assert bundled.list_bundled() == {"isalgo": ["TRB", "mouse_TRB"], "olga-legacy": []}
    generate = importlib.import_module("vdjtools.model.generate")
    routed = []
    monkeypatch.setattr(generate, "generate", lambda model, n, **kw: routed.append(
        (model.manifest.source, kw["engine"])))
    for source in ("isalgo", "learned", "arda"):
        corpus.draw_pool("TRB", 1, seed=3, source=source)
    assert routed == [("new-isalgo-toy", "native")] * 3



def test_olga_alias_preserves_frozen_fingerprints_and_sampler(monkeypatch):
    assert corpus.model_fingerprint(bundled.LOCI, "olga") == corpus.model_fingerprint(bundled.LOCI, "olga-legacy")
    generate = importlib.import_module("vdjtools.model.generate")
    calls = []
    monkeypatch.setattr(generate, "generate", lambda model, n, **kw: calls.append(kw))
    for source in ("olga", "olga-legacy"):
        corpus.draw_pool("TRB", 3, seed=42, source=source)
    assert calls[0] == calls[1]
    assert calls[0]["engine"] == "reference"


def test_new_defaults_agree_and_explicit_legacy_stays_loadable():
    from vdjtools import cli
    pgen = importlib.import_module("vdjtools.sc.pgen")
    alice = importlib.import_module("vdjtools.overlap.alice")
    assert inspect.signature(pgen.paired_pgen).parameters["source"].default == "isalgo"
    assert inspect.signature(alice.alice).parameters["source"].default == "isalgo"
    for fn in (corpus.draw_pool, corpus.synthesize, corpus.build_matrices, corpus.sample_stream, corpus.model_fingerprint):
        assert inspect.signature(fn).parameters["source"].default == "isalgo"
    assert all(source in bundled.SOURCES for source in ("olga", "learned", "arda"))
    import typer
    command = typer.main.get_command(cli.app)
    for name in ("pgen", "generate"):
        assert next(p for p in command.commands[name].params if p.name == "source").default == "isalgo"


def test_production_producer_rejects_precomputed_inputs(monkeypatch):
    monkeypatch.delenv("BUNDLE_SOURCE", raising=False)
    monkeypatch.setenv("EM_DATA_DIR", "/unused/precomputed")
    script = Path(__file__).resolve().parents[2] / "appendix" / "build_bundled_models.py"
    namespace = runpy.run_path(str(script))
    assert namespace["SOURCE"] == "isalgo"
    assert namespace["DEST"].as_posix().endswith("_bundled/isalgo")
    with pytest.raises(ValueError, match="fresh raw inputs"):
        namespace["training_clones"]("TRB")


def test_isalgo_release_catalog_has_twelve_models():
    """Release gate: absent fresh assets must remain a visible failure until fitted."""
    expected = set(bundled.LOCI) | {f"mouse_{loc}" for loc in ("TRA", "TRB", "IGH", "IGK", "IGL")}
    assert set(bundled.list_bundled()["isalgo"]) == expected
    for key in expected:
        organism, locus = ("mouse", key[6:]) if key.startswith("mouse_") else ("human", key)
        model = bundled.load_bundled(locus, organism=organism, collapse=False)
        assert (model.organism, model.locus) == (organism, locus)
