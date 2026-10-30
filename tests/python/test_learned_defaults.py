"""Working model routing without scientific generation, fitting or new assets."""
import inspect

import pytest


@pytest.mark.parametrize('organism,loci', [('human', ('TRA','TRB','TRG','TRD','IGH','IGK','IGL')),
                                         ('mouse', ('TRA','TRB','IGH','IGK','IGL'))])
def test_default_loader_routes_catalog(tmp_path, monkeypatch, organism, loci):
    from vdjtools.model import bundled
    paths = []
    monkeypatch.setattr(bundled, 'files', lambda name: tmp_path)
    monkeypatch.setattr(bundled, 'load_model', lambda path: paths.append(path) or path)
    for locus in loci:
        key = locus if organism == 'human' else 'mouse_' + locus
        path = tmp_path / '_bundled' / 'isalgo' / key
        path.mkdir(parents=True)
        assert bundled.load_bundled(locus, organism=organism, collapse=False) == path
    assert len(paths) == len(loci)
    assert bundled._bundled_key('olga', 'TRB', 'human') == 'TRB'
    assert bundled._bundled_key('arda', 'TRA', 'mouse') == 'mouse_TRA'
    with pytest.raises(ValueError):
        bundled._bundled_key('olga', 'TRB', 'mouse')


def test_working_consumer_defaults_and_explicit_cli_oracle(monkeypatch):
    from vdjtools import cli
    from vdjtools.model import bundled
    from vdjtools.overlap.alice import alice
    from vdjtools.sc.pgen import paired_pgen
    from vdjtools.signature import corpus

    for function in (bundled.load_bundled, alice, paired_pgen, corpus.model_fingerprint,
                     corpus.draw_pool, corpus.build_matrices, corpus.sample_stream, corpus.synthesize):
        assert inspect.signature(function).parameters['source'].default == 'isalgo'
    for name, function in inspect.getmembers(cli, inspect.isfunction):
        parameter = inspect.signature(function).parameters.get('source')
        if parameter is not None and hasattr(parameter.default, 'default'):
            assert parameter.default.default == 'isalgo', name
    calls = []
    monkeypatch.setattr(bundled, 'load_bundled', lambda *args, **kwargs: calls.append((args, kwargs)))
    import vdjtools.model as model
    monkeypatch.setattr(model, 'load_bundled', bundled.load_bundled)
    cli._model_arg('TRB')
    cli._model_arg('TRB:olga')
    assert [kwargs['source'] for args, kwargs in calls] == ['isalgo', 'olga']
