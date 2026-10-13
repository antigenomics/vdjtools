"""The existing bundle builder targets only the requested source, organisms and loci."""
import importlib.util
from pathlib import Path
from types import SimpleNamespace


def test_arda_bundle_builder_targets_human_and_mouse_tra(tmp_path, monkeypatch):
    monkeypatch.setenv('BUNDLE_SOURCE', 'arda')
    monkeypatch.setenv('BUILD_GROUPS', 'human,mouse')
    monkeypatch.setenv('LOCI', 'TRA')
    path = Path(__file__).resolve().parents[2] / 'appendix' / 'build_bundled_models.py'
    spec = importlib.util.spec_from_file_location('bundle_builder', path)
    builder = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(builder)
    monkeypatch.setattr(builder, 'DEST', tmp_path)
    monkeypatch.setattr(builder, 'training_clones', lambda locus, group: None)
    def reject_olga(*args):
        raise AssertionError('An arda rebuild must not use an OLGA template')
    monkeypatch.setattr(builder, 'seed_model', reject_olga)
    calls, saved = [], []
    def build_model(locus, **kw):
        calls.append((locus, kw['group'], kw['template']))
        return SimpleNamespace(save=saved.append), None, {
            'n_used': 1, 'n_clonotypes': 1, 'iters': 1, 'seconds': 1,
            'p_nd2': None, 'loglik_first': -2, 'loglik_last': -1,
        }
    monkeypatch.setattr(builder.data, 'build_model', build_model)
    builder.main()
    assert calls == [('TRA', 'human', None), ('TRA', 'mouse', None)]
    assert saved == [tmp_path / 'human_TRA', tmp_path / 'mouse_TRA']
