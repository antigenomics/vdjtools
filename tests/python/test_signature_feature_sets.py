"""Feature definitions preserve fixed bases and the legacy signature contract."""
import itertools

import pytest

from vdjtools.signature import layout as L


def test_feature_definitions_and_fixed_order():
    definitions = []
    for configuration in ('basic', 'joint'):
        definitions.append(L.feature_definition(configuration))
        for motifs, context in itertools.product(('imgt11', 'blosum62'), ('none', 'gene', 'cluster')):
            motif_map = None
            if motifs == 'blosum62':
                motif_map = {'mapping': {''.join(t): 'm000' for t in itertools.product(L.AMINO_ACIDS, repeat=3)},
                             'version': 1}
            definitions.append(L.feature_definition(configuration, motifs=motifs, v_context=context,
                                motif_map=motif_map, v_maps={'TRB': {'version': 1, 'mapping': {'TRBV1': 'v000'}}}
                                if context == 'cluster' else None))
    assert len({d['identity'] for d in definitions}) == 14
    d = L.feature_definition('basic', motifs='imgt11')
    assert d['classes'] == ['AILV', 'CM', 'ST', 'HKR', 'DE', 'NQ', 'F', 'W', 'Y', 'P', 'G']
    vocab = {'vj': ['TRBV1|TRBJ1', '_unknown'], 'length': ['6', '7'],
             'motif_XXX': [''.join(t) for t in itertools.product('0123456789a', repeat=3)],
             'motif_XX_X': [''.join(t) for t in itertools.product('0123456789a', repeat=3)],
             'motif_X_XX': [''.join(t) for t in itertools.product('0123456789a', repeat=3)]}
    cols = L.raw_columns('vsig', 'TRB', vocab, feature_definition=d)
    assert sum(':motif_' in c for c in cols) == 3993
    assert cols == L.raw_columns('vsig', 'TRB', dict(reversed(list(vocab.items()))), feature_definition=d)
    assert not L.raw_columns('vsig', '-', {}, feature_definition=d)
    assert L.raw_columns('vsig', 'TRB', {'vus': ['V'], 'jus': ['J'], 'spec': ['V_6']}) == L.raw_columns(
        'vsig', 'TRB', {'vus': ['V'], 'jus': ['J'], 'spec': ['V_6']}, feature_definition=L.feature_definition())


def test_definition_rejects_incompatible_or_tampered_parameters():
    with pytest.raises(ValueError):
        L.feature_definition('basic', v_context='gene')
    with pytest.raises(ValueError):
        L.feature_definition('joint', motifs='blosum62', motif_map={'version': 2, 'mapping': {'AAA': 'm0'}})
    d = L.feature_definition('basic', motifs='imgt11')
    changed = {**d, 'shapes': ['XXXX']}
    with pytest.raises(ValueError):
        L.validate_feature_definition(changed)
