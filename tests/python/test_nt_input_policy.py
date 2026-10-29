"""Sequencing DNA is normalized before grouping/scoring; AA and metadata are untouched."""
import importlib

import polars as pl
import pytest

from vdjtools.io import schema


def test_scalar_and_vector_dna_policy_preserves_other_columns():
    text = 'acgtryswkmbdhvn'
    assert schema.normalize_nt(text) == 'ACGT' + 'A' * 11
    assert schema.normalize_nt(None) is None
    assert schema.normalize_nt('') == ''
    frame = pl.DataFrame({'junction_nt': [text, None, ''], 'junction_aa': ['CARN', 'RYS', 'N'],
                          'sequence_id': ['N1', 'R2', 'Y3'], 'quality': ['NN!', '?R', '']})
    result = schema.normalize_nt_frame(frame, ['junction_nt', 'absent'])
    assert result['junction_nt'].to_list() == ['ACGT' + 'A' * 11, None, '']
    assert result.drop('junction_nt').equals(frame.drop('junction_nt'))
    assert schema.normalize(frame)['junction_nt'].to_list() == result['junction_nt'].to_list()


@pytest.mark.parametrize('invalid', ['ACXG', 'AC-G', 'ACUG'])
def test_invalid_dna_is_rejected_in_scalar_and_frame(invalid):
    with pytest.raises(ValueError, match='DNA IUPAC'):
        schema.normalize_nt(invalid)
    with pytest.raises(ValueError, match='DNA IUPAC'):
        schema.normalize_nt_frame(pl.DataFrame({'junction_nt': [invalid]}), ['junction_nt'])


def test_airr_normalizes_before_grouping_and_preserves_count_and_aa(tmp_path):
    from vdjtools.io.read import read_airr
    path = tmp_path / 'synthetic.airr.tsv'
    pl.DataFrame({'v_call': ['TRBV1*01'] * 2, 'j_call': ['TRBJ1-1*01'] * 2,
                  'junction': ['ACNT', 'ACAT'], 'junction_aa': ['CARN'] * 2,
                  'duplicate_count': [2, 3], 'sequence': ['NRY', 'AAA']}).write_csv(path, separator='\t')
    result = read_airr(path, keep=('sequence',))
    assert result.height == 1
    assert result['junction_nt'].to_list() == ['ACAT']
    assert result['junction_aa'].to_list() == ['CARN']
    assert result['duplicate_count'].to_list() == [5]
    assert result['frequency'].to_list() == [1.0]
    assert result['sequence'].to_list() == ['AAA']


def test_format_finalizer_normalizes_before_collapse():
    from vdjtools.io.convert import _finalize
    frame = pl.DataFrame({'v_call': ['TRBV1'] * 2, 'j_call': ['TRBJ1-1'] * 2,
                          'junction_nt': ['ACNT', 'ACAT'], 'junction_aa': ['CARN'] * 2,
                          'duplicate_count': [7, 11]})
    result = _finalize(frame)
    assert result.height == 1
    assert result['duplicate_count'].to_list() == [18]
    assert result['junction_nt'].to_list() == ['ACAT']
    assert result['junction_aa'].to_list() == ['CARN']
    assert _finalize(frame, collapse=False).height == 2


def test_sequencing_translation_normalizes_dna_before_deriving_aa():
    from vdjtools.io.convert import translate_junctions
    assert translate_junctions(['TGN', None, ''], threads=1).to_list() == ['*', None, None]


def test_scoring_types_nt_from_column_and_preserves_explicit_aa(monkeypatch):
    score = importlib.import_module('vdjtools.model.score')
    seen = []
    monkeypatch.setattr(score, '_pgen_nt_many', lambda model, seqs, *a: seen.append(('nt', seqs)) or [0.5] * len(seqs))
    monkeypatch.setattr(score.native, 'pgen_aa_batch', lambda model, seqs, **kw: seen.append(('aa', seqs)) or [0.5] * len(seqs))
    result = score.pgen_frame(None, pl.DataFrame({'junction_nt': ['ACNR']}), use_calls=False)
    assert result['sequence'].to_list() == ['ACAA']
    score.pgen_frame(None, pl.DataFrame({'junction_aa': ['CARN']}), use_calls=False)
    score.pgen_frame(None, ['CARN'], kind='aa', use_calls=False)
    assert seen == [('nt', ['ACAA']), ('aa', ['CARN']), ('aa', ['CARN'])]
