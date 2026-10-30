"""The accepted independent functional validator stays fresh, exact and ambiguity-aware."""
import importlib.util
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import polars as pl
import pytest

pytest.importorskip('olga')
_SPEC = importlib.util.spec_from_file_location(
    'appendix_compare_models', Path(__file__).resolve().parents[2] / 'appendix/compare_models.py')
compare = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(compare)


@pytest.mark.parametrize('chain,channels', [('VJ', ('PinsVJ',)), ('VDJ', ('PinsVD', 'PinsDJ'))])
def test_oracle_workspace_padding_preserves_every_probability(chain, channels, tmp_path, monkeypatch):
    class Loader:
        def load_igor_genomic_data(self, *args):
            pass
        def load_and_process_igor_model(self, *args):
            pass
    original = {channel: np.array([.2, .3, .5]) for channel in channels}
    monkeypatch.setattr(compare.olm, 'GenomicData' + chain, Loader)
    monkeypatch.setattr(compare.olm, 'GenerativeModel' + chain, Loader)
    monkeypatch.setattr(compare.ogp, 'GenerationProbability' + chain,
                        lambda *args: SimpleNamespace(**{k: v.copy() for k, v in original.items()}))
    oracle = compare.olga_pgen_model(tmp_path, chain)
    for channel, probabilities in original.items():
        padded = getattr(oracle, channel)
        np.testing.assert_array_equal(padded, np.r_[probabilities, np.zeros(3)])


def test_engine_oracle_import_keeps_raw_germlines_and_restores_repair(toy_model, monkeypatch, tmp_path):
    from vdjtools.model import io as model_io
    original = model_io.repair_anchors
    genomic = object()
    def load(path, **kwargs):
        assert path == tmp_path and kwargs == {'locus': 'TRB', 'organism': 'human'}
        assert model_io.repair_anchors(genomic) is genomic
        return toy_model
    monkeypatch.setattr(compare, 'from_olga', load)
    assert compare.raw_olga_model(tmp_path, 'TRB', 'human') is toy_model
    assert model_io.repair_anchors is original


def test_run_uses_raw_engine_oracle_but_keeps_loaded_performance_baseline(toy_model, monkeypatch):
    raw = object()
    monkeypatch.setattr(compare.Model, 'load', lambda *a: toy_model)
    monkeypatch.setattr(compare, 'raw_olga_model', lambda *a: raw)
    monkeypatch.setattr(compare, 'functional_queries', lambda *a: (['AAA'], ['V'], ['J'], {}))
    class Oracle:
        def compute_nt_CDR3_pgen(self, seq):
            return .2
    monkeypatch.setattr(compare, 'olga_pgen_model', lambda *a: Oracle())
    def score(model, *a, **kw):
        assert model is raw
        return [.2]
    monkeypatch.setattr(compare.native, 'pgen_nt_batch', score)
    baselines = []
    def conditioned(model, *a):
        baselines.append(model)
        return [.1]
    monkeypatch.setattr(compare, '_gene_conditioned_pgens', conditioned)
    assert compare.run('TRB')['oracle']['all_values_agree']
    assert len(baselines) == 2 and all(model is toy_model for model in baselines)


def test_functional_queries_always_reannotate_in_fresh_removed_directory(tmp_path, monkeypatch):
    monkeypatch.setenv('RAW_FASTQ_DIR', str(tmp_path / 'raw'))
    monkeypatch.setattr(compare, 'WORK', tmp_path / 'work')
    monkeypatch.setattr(compare, 'N', 2)
    directories = []
    def annotate(raw, *, out_dir, prefix, organism, cap, threads):
        assert raw == tmp_path / 'raw/mouse/TRB.functional.fq.gz'
        assert organism == 'mouse' and cap == 20000
        assert Path(out_dir).is_dir()
        assert not list(Path(out_dir).iterdir())
        directories.append(Path(out_dir))
        return pl.DataFrame({'v_call': ['V*01'] * 3, 'j_call': ['J*01'] * 3,
                             'junction': ['TGTAAATTT', 'TGTCCCTTT', 'TGTGGGTTT']})
    monkeypatch.setattr(compare.data, 'annotate_reads', annotate)
    monkeypatch.setattr(pl, 'read_parquet', lambda *a, **k: pytest.fail('cached parquet read'))
    for _ in range(2):
        sequences, vs, js, counts = compare.functional_queries('mouse', 'TRB')
        assert len(sequences) == len(vs) == len(js) == counts['selected_clonotypes'] == 2
        assert counts['eligible_clonotypes'] == 3
    assert len(set(directories)) == 2
    assert not any(path.exists() for path in directories)


def test_csv_genes_and_suffix_aliases_sum_disjoint_states_in_one_native_batch(toy_model, monkeypatch):
    from vdjtools.model import native
    sequence = 'TGTGCCAGCAACTATGGCTATACCTTT'
    vs = toy_model.genomic['genes_v']['v_allele'].to_list()
    js = toy_model.genomic['genes_j']['j_allele'].to_list()
    expected = sum(native.pgen_nt(toy_model, sequence, v, j) for v in vs for j in js)
    assert expected > 0
    original = native.pgen_nt_batch
    calls = []
    def batch(model, seqs, v, j, *, threads):
        calls.append((seqs, v, j))
        return original(model, seqs, v, j, threads=1)
    monkeypatch.setattr(native, 'pgen_nt_batch', batch)
    scores = compare._gene_conditioned_pgens(
        toy_model, [sequence] * 3,
        ['TOYV1*99, TOYV2,TOYV1*01,UNKNOWN', 'UNKNOWN', 'TOYV1'],
        ['TOYJ2*42,TOYJ1,TOYJ1*01', 'TOYJ1', 'UNKNOWN'], threads=1)
    np.testing.assert_allclose(scores, [expected, 0, 0], rtol=1e-12, atol=0)
    assert len(calls) == 1
    assert len(calls[0][0]) == len(vs) * len(js)
    assert len(set(zip(calls[0][1], calls[0][2]))) == len(vs) * len(js)
    assert None not in calls[0][1] and None not in calls[0][2]
    assert compare._pgen_summary(scores)['zero'] == 2
    summary = compare._pgen_summary(scores)
    assert summary['positive'] == 1 and summary['total'] == 3
    assert summary['mean_log10_all_is_negative_infinity']
    assert summary['mean_pgen'] == pytest.approx(expected / 3)


def test_unknown_only_never_invokes_unrestricted_native_scoring(toy_model, monkeypatch):
    monkeypatch.setattr(compare.native, 'pgen_nt_batch', lambda *a, **k: pytest.fail('unknown scored unrestricted'))
    assert compare.held_out_loglik(toy_model, ['TGTAAATTT'], ['UNKNOWN'], ['UNKNOWN']) == (None, 1, 1)


@pytest.mark.parametrize('group', ['human', 'mouse'])
@pytest.mark.parametrize('ours,theirs', [([0.1, 0], [0, 0]), ([0, 0], [0.1, 0])])
def test_oracle_assertion_includes_zero_disagreements(toy_model, monkeypatch, ours, theirs, group):
    loaded, imported = [], []
    def load(path):
        loaded.append(path)
        return toy_model
    def import_olga(path, **kw):
        imported.append((path, kw))
        return toy_model
    monkeypatch.setattr(compare, 'from_olga', import_olga)
    monkeypatch.setattr(compare.Model, 'load', load)
    oracle = iter(theirs)
    class Oracle:
        def compute_nt_CDR3_pgen(self, seq):
            return next(oracle)
    monkeypatch.setattr(compare, 'olga_pgen_model', lambda *a, **k: Oracle())
    monkeypatch.setattr(compare, 'functional_queries', lambda *a: (['AAA', 'CCC'], ['V'] * 2, ['J'] * 2, {}))
    monkeypatch.setattr(compare.native, 'pgen_nt_batch', lambda *a, **k: ours)
    with pytest.raises(AssertionError, match='exact OLGA oracle mismatch'):
        compare.run('TRB', group)
    key = 'TRB' if group == 'human' else 'mouse_TRB'
    assert compare.BUNDLED / 'isalgo' / key in loaded
    if group == 'human':
        assert compare.BUNDLED / 'olga/TRB' in loaded
        assert imported == [(compare.OLGA / 'human_T_beta', {'locus': 'TRB', 'organism': 'human'})]
    else:
        assert loaded == [compare.BUNDLED / 'isalgo/mouse_TRB']
        assert imported == [(compare.MOUSE_OLGA / 'mouse_T_beta', {'locus': 'TRB', 'organism': 'mouse'})] * 2


def test_default_context_selection_covers_exactly_twelve(monkeypatch, capsys):
    monkeypatch.delenv('BUILD_GROUPS', raising=False)
    monkeypatch.delenv('LOCI', raising=False)
    contexts = []
    def run(locus, group):
        contexts.append((group, locus))
        return {'group': group, 'locus': locus}
    monkeypatch.setattr(compare, 'run', run)
    compare.main()
    assert len(contexts) == len(set(contexts)) == 12
    assert ('human', 'TRD') in contexts and ('mouse', 'IGL') in contexts
    assert ('mouse', 'TRD') not in contexts


@pytest.mark.parametrize('values,positive,correlation', [([0, .5, .25], 2, 1.), ([0, 0, 0], 0, None)])
def test_oracle_diagnostics_name_their_positive_pair_denominator(toy_model, monkeypatch, values, positive, correlation):
    toy_model.tables['vj_ins'] = pl.DataFrame({'length': [0, 1, 2, 3], 'p': [.6, .1, .2, .1]})
    monkeypatch.setattr(compare.Model, 'load', lambda *a, **k: toy_model)
    monkeypatch.setattr(compare, 'from_olga', lambda *a, **k: toy_model)
    results = iter(values)
    class Oracle:
        def compute_nt_CDR3_pgen(self, seq):
            return next(results)
    monkeypatch.setattr(compare, 'olga_pgen_model', lambda *a, **k: Oracle())
    monkeypatch.setattr(compare, 'functional_queries', lambda *a: (
        ['AAA', 'CCC', 'GGG'], ['UNKNOWN'] * 3, ['UNKNOWN'] * 3, {}))
    monkeypatch.setattr(compare.native, 'pgen_nt_batch', lambda *a, **k: values)
    result = compare.run('TRB', 'human')
    assert result['insertion_zero'] == {'olga-legacy': {'vj_ins': .6}, 'isalgo': {'vj_ins': .6}}
    report = result['oracle']
    assert report['total'] == 3 and report['both_positive'] == positive
    assert report['log10_correlation_both_positive'] == correlation
    assert report['max_relative_error_both_positive'] == (0 if positive else None)
    assert report['ours_zero'] == report['oracle_zero'] == 3 - positive
