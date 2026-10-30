"""Sequencing DNA IUPAC symbols are always replaced with A before model inference."""
import inspect

import numpy as np
import polars as pl
import pytest

from test_insertion_ambiguity import V, J, _ambiguous_model
from vdjtools import _core
from vdjtools.model import data, native
from vdjtools.model.infer import (
    InferenceReport, _candidate_identity, candidate_masks, infer, infer_frame, infer_native,
    sanitize_junctions,
)


@pytest.mark.parametrize('symbol', list('RYSWKMBDHVN'))
@pytest.mark.parametrize('rate', [None, .03])
def test_each_iupac_base_is_identical_to_observed_a(symbol, rate):
    model = _ambiguous_model()
    sequence = V[:-2] + symbol + 'A' + J
    actual, report = infer_native(model, [sequence], init='template', max_iter=1, error_rate=rate)
    expected, reference = infer_native(model, [V[:-2] + 'AA' + J], init='template', max_iter=1, error_rate=rate)
    assert report.loglik == reference.loglik
    assert report.n_scoreable == reference.n_scoreable == [1]
    assert report.n_ambiguous == 1
    assert report.error_rate == reference.error_rate
    for name in actual.tables:
        assert actual.tables[name].equals(expected.tables[name]), name
    assert 'max_observation_states' not in report.to_dict()


@pytest.mark.parametrize('engine', [infer, infer_native])
def test_direct_reference_and_native_inference_normalize_before_initialization(engine):
    model = _ambiguous_model()
    actual, report = engine(model, [V[:-2] + 'NN' + J], init='align', max_iter=1)
    expected, reference = engine(model, [V[:-2] + 'AA' + J], init='align', max_iter=1)
    assert report.loglik == reference.loglik
    for name in actual.tables:
        assert actual.tables[name].equals(expected.tables[name]), name


@pytest.mark.parametrize('policy', [None, 'C', 'G', 'T', 'marginalize'])
def test_contrary_ambiguity_options_are_rejected_even_for_acgt(policy):
    model = _ambiguous_model()
    frame = pl.DataFrame({'junction': [V + J]})
    with pytest.raises(ValueError, match="ambiguous must be 'A'"):
        sanitize_junctions(frame, 'junction', ambiguous=policy)
    with pytest.raises(ValueError, match="ambiguous must be 'A'"):
        infer_frame(model, frame, use_calls=False, ambiguous=policy, max_iter=1)
    with pytest.raises(ValueError, match="ambiguous must be 'A'"):
        data.build_model('TOY', template=model, clones=frame, ambiguous=policy, iters=1)


def test_frame_default_and_gene_grouping_use_the_normalized_junction():
    frame = pl.DataFrame({'junction': [V[:-2] + 'NN' + J, V[:-2] + 'AA' + J],
                          'v_call': ['TOYV1*01', 'TOYV1*99'], 'j_call': ['TOYJ1*01'] * 2,
                          'count': [3, 7]})
    clones = data.unique_clonotypes(frame, gene_level=True)
    assert clones.height == 1 and clones['count'].item() == 10
    assert clones['junction'].item() == V[:-2] + 'AA' + J
    _, report = infer_frame(_ambiguous_model(), frame, use_calls=False, init='template', max_iter=1)
    _, expected = infer_native(_ambiguous_model(), [V[:-2] + 'AA' + J] * 2, init='template', max_iter=1)
    assert report.loglik == expected.loglik and report.n_scoreable == [2]


def test_full_raw_reads_are_normalized_before_reverse_complement():
    frame = pl.DataFrame({'sequence': ['NRYT'], 'rev_comp': ['T'], 'junction': ['NRY'],
                          'candidate_sequences': [['nry', 'GG']]})
    normalized = data._with_candidate_sequences(frame)
    assert normalized['sequence'].item() == 'AAAT'
    assert normalized['junction'].item() == 'AAA'
    assert normalized['candidate_sequences'].to_list() == [['AAA', 'ATTT', 'GG']]


def test_candidate_queries_and_digests_use_the_same_hard_a_policy(monkeypatch, toy_model):
    import arda
    seen = []
    def candidates(records, **kw):
        seen.extend(records)
        return {'alignments': pl.DataFrame(schema={'sequence_id': pl.String,
                'segment': pl.String, 'alleles': pl.List(pl.String)}), 'limits': {}}
    monkeypatch.setattr(arda, 'alignment_candidates', candidates, raising=False)
    masks, evidence = candidate_masks(['NR'], toy_model, query_sequences=[['NRYT', 'AAAT']], threads=1)
    assert seen == [('0', 'AAAT')]
    assert evidence['candidate_ambiguity'] == 'DNA IUPAC symbols replaced with A before search and inference'
    assert evidence['input_masks_sha256'] == _candidate_identity(['AA'], masks)
    assert _candidate_identity(['NR'], masks) == _candidate_identity(['AA'], masks)


def test_nt_pgen_and_scenarios_normalize_observations_but_not_germlines():
    from vdjtools.model.pgen import pgen_nt, prepare
    from vdjtools.model.viterbi import best_scenario
    model = _ambiguous_model()
    raw, normalized = V[:-2] + 'NA' + J, V[:-2] + 'AA' + J
    assert native._encode_observed('nryacgt') == native._encode('AAAACGT')
    with pytest.raises(KeyError):
        native._encode('N')
    assert native.pgen_nt(model, raw) == native.pgen_nt(model, normalized)
    assert pgen_nt(prepare(model), raw) == pgen_nt(prepare(model), normalized)
    np.testing.assert_array_equal(native.pgen_nt_batch(model, [raw, normalized]),
                                   [native.pgen_nt(model, normalized)] * 2)
    assert best_scenario(model, raw) == best_scenario(model, normalized)


@pytest.mark.parametrize('encoded', [[-1], [-15], [4], [-16], [127], [-128]])
def test_native_observations_reject_all_non_acgt_indices(encoded):
    packed, _, _ = native.pack(_ambiguous_model())
    with pytest.raises(ValueError, match='ACGT'):
        _core.estep_batch(packed, [encoded], [], [], [], _core.make_counts(packed))


def test_removed_completion_budget_is_not_a_model_option():
    assert 'max_observation_states' not in inspect.signature(infer_native).parameters
    assert 'max_observation_states' not in inspect.signature(data.build_model).parameters
    assert 'max_observation_states' not in InferenceReport.__dataclass_fields__
    with pytest.raises(ValueError, match='DNA IUPAC'):
        native._encode_observed('ACXG')


def test_resume_does_not_restore_the_removed_completion_budget():
    from vdjtools.model.infer import resume
    model = _ambiguous_model()
    model.training = {'runs': [{'max_observation_states': 4096}]}
    _, report = resume(model, [V[:-2] + 'NN' + J], max_iter=1)
    assert report.n_scoreable == [1]
    assert 'max_observation_states' not in report.to_dict()


def test_list_only_candidate_inputs_normalize_without_losing_null_or_empty_lists():
    frame = pl.DataFrame({'candidate_sequences': [['nRy', None, 'AC'], None, []]})
    normalized = data._with_candidate_sequences(frame)
    assert normalized['candidate_sequences'].to_list() == [['AAA', None, 'AC'], None, []]
    all_null = pl.DataFrame({'candidate_sequences': [None]})
    assert data._with_candidate_sequences(all_null)['candidate_sequences'].to_list() == [None]
    with pytest.raises(ValueError, match='DNA IUPAC'):
        data._with_candidate_sequences(pl.DataFrame({'candidate_sequences': [['ACX']]}))
    with_raw = pl.DataFrame({'candidate_sequences': [None], 'sequence': ['NRYT']},
                           schema={'candidate_sequences': pl.List(pl.String), 'sequence': pl.String})
    assert data._with_candidate_sequences(with_raw)['candidate_sequences'].to_list() == [['AAAT']]
