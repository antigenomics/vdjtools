"""The existing bundle builder targets only the requested source, organisms and loci."""
import importlib.util
from pathlib import Path
from types import SimpleNamespace
import pytest


@pytest.mark.parametrize('source', ['arda', 'learned'])
@pytest.mark.parametrize('want', ['TRA', None])
def test_bundle_builder_targets_human_and_mouse(tmp_path, monkeypatch, source, want):
    monkeypatch.setenv('BUNDLE_SOURCE', source)
    monkeypatch.setenv('BUILD_GROUPS', 'human,mouse')
    if want:
        monkeypatch.setenv('LOCI', want)
    else:
        monkeypatch.delenv('LOCI', raising=False)
    path = Path(__file__).resolve().parents[2] / 'appendix' / 'build_bundled_models.py'
    spec = importlib.util.spec_from_file_location('bundle_builder', path)
    builder = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(builder)
    monkeypatch.setattr(builder, 'DEST', tmp_path)
    monkeypatch.setattr(builder, 'training_clones', lambda locus, group: None)
    calls, saved = [], []
    def build_model(locus, **kw):
        calls.append((locus, kw['group'], kw['template']))
        return SimpleNamespace(save=saved.append), None, {
            'n_used': 1, 'n_clonotypes': 1, 'iters': 1, 'seconds': 1,
            'p_nd2': None, 'loglik_first': -2, 'loglik_last': -1,
        }
    monkeypatch.setattr(builder.data, 'build_model', build_model)
    builder.main()
    targets = [(group,locus) for group,loci in [('human',builder.LOCI),('mouse',builder.MOUSE_LOCI)]
               for locus in loci if want is None or locus == want]
    assert calls == [(locus,group,None) for group,locus in targets]
    assert saved == [tmp_path / (f'{group}_{locus}' if source == 'arda' or group != 'human' else locus)
                     for group,locus in targets]
