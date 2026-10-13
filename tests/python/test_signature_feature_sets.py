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


def test_three_shapes_and_per_clone_normalization():
    import polars as pl
    from vdjtools.features import kmer as K
    frame = pl.DataFrame({'junction_aa': ['CADGKF'], 'duplicate_count': [1], 'v_call': ['TRBV1*01']})
    got = K.motif_profile(frame, alphabet='aa20', weight='distinct')
    by_shape = {shape: dict(got.filter(pl.col('shape') == shape).select('motif', 'frequency').iter_rows())
                for shape in L.MOTIF_SHAPES}
    assert by_shape == {'XXX': {'ADG': 0.5, 'DGK': 0.5}, 'XX_X': {'ADK': 1.0}, 'X_XX': {'AGK': 1.0}}
    frame = pl.DataFrame({'junction_aa': ['CADGF', 'CKKKKKF'], 'duplicate_count': [1, 1],
                          'v_call': ['TRBV1', None]})
    got = K.motif_profile(frame, alphabet='aa20', weight='distinct', v_context='gene', v_map={'TRBV1': 'TRBV1'})
    triples = got.filter(pl.col('shape') == 'XXX')
    assert triples.filter(pl.col('motif') == 'ADG')['frequency'].item() == 0.5
    assert triples.filter(pl.col('motif') == 'KKK')['frequency'].item() == 0.5
    assert triples.filter(pl.col('motif') == 'KKK')['context'].item() == '_unknown'


def test_motif_eligibility_and_no_cross_clonotype_windows():
    import polars as pl
    from vdjtools.features import kmer as K
    frame = pl.DataFrame({'junction_aa': ['CADF', 'CKGF', 'CAD*KF', 'AADGKF'],
                          'duplicate_count': [1, 1, 1, 1]})
    got = K.motif_profile(frame, alphabet='aa20', weight='distinct')
    assert got.is_empty()
    valid = pl.DataFrame({'junction_aa': ['CADGKF'], 'duplicate_count': [1]})
    got = K.motif_profile(valid, alphabet='imgt11', weight='distinct')
    assert got.filter(pl.col('shape') == 'XX_X')['motif'].item() == '043'


def test_new_composition_blocks_have_residuals_and_uncalled_isotype():
    import polars as pl
    from vdjtools.signature import features as F
    d = L.feature_definition('basic')
    vocab = {'vus': ['TRBV1'], 'jus': ['TRBJ1'], 'vj': ['TRBV1|TRBJ1', '_unknown'],
             'length': [str(i) for i in L.SPECTRATYPE_LENGTHS]}
    frame = pl.DataFrame({'junction_aa': ['CADGKF', 'CADGQF', 'CADGRF'],
                          'v_call': ['TRBV1*01', None, 'TRBV1,TRBV2'],
                          'j_call': ['TRBJ1*01'] * 3, 'duplicate_count': [1] * 3})
    groups = F.composition_groups(F.work_frame(frame, 'distinct'), vocab, feature_definition=d)
    assert groups['vj']['TRBV1|TRBJ1'] == pytest.approx(1/3)
    assert groups['vj']['_unknown'] == pytest.approx(2/3)
    vocab['isotype'] = [*F.ISOTYPES, '_uncalled']
    groups = F.composition_groups(F.work_frame(frame, 'distinct'), vocab, feature_definition=d)
    assert groups['isotype']['_uncalled'] == 1.0
    assert sum(groups['length'].values()) == 1.0


def test_cluster_builder_small_matrix_and_validation(tmp_path):
    import importlib.util
    from pathlib import Path
    import numpy as np
    path = Path(__file__).parents[2] / 'appendix' / 'build_motif_vocab.py'
    assert path.exists(), 'frozen motif builder is missing'
    spec = importlib.util.spec_from_file_location('motif_builder', path)
    builder = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(builder)
    dist = np.array([[0, 1, 9, 9], [1, 0, 9, 9], [9, 9, 0, 2], [9, 9, 2, 0]], dtype=float)
    got = builder.cluster_map(['AAA', 'AAC', 'CCA', 'CCC'], dist, n_clusters=2, prefix='m')
    assert got['AAA'] == got['AAC']
    assert got['CCA'] == got['CCC']
    assert got['AAA'] != got['CCA']
    assert got == builder.cluster_map(['AAA', 'AAC', 'CCA', 'CCC'], dist, n_clusters=2, prefix='m')
    with pytest.raises(ValueError):
        builder.cluster_map(['AAA', 'AAC', 'CCA', 'CCC'], dist + 1, n_clusters=2, prefix='m')
    matrix_path = tmp_path / 'human_TRB.npz'
    np.savez(matrix_path, V_alleles=['TRBV1*01', 'TRBV2*01'], V_dist=[[0, 3], [3, 0]])
    got = builder.build_v_map(matrix_path, ['TRBV1', 'TRBV2'], organism='human', locus='TRB', n_clusters=2)
    assert got['representatives'] == {'TRBV1': 'TRBV1*01', 'TRBV2': 'TRBV2*01'}
    with pytest.raises(ValueError):
        builder.build_v_map(matrix_path, ['TRBV9'], organism='human', locus='TRB')
    unclustered = builder.build_v_map(matrix_path, ['TRBV1', 'TRBV9'], organism='human', locus='TRB',
                                     unclustered_genes=['TRBV9'])
    assert unclustered['mapping']['TRBV9'] == '_unknown'
    assert unclustered['unclustered_genes'] == ['TRBV9']
    penalties = builder.motif_distances(['AAA', 'AAC'], threads=1)
    assert penalties[0, 0] == 0
    assert penalties[0, 1] == penalties[1, 0] > 0


def test_sparse_block_fit_retains_zero_mad_bins_and_roundtrips(tmp_path):
    import numpy as np
    from vdjtools.signature import corpus as C
    d = L.feature_definition('basic')
    vocab = {'TRB': {'vj': ['V|J', '_unknown'], 'length': ['6', '7']}}
    columns = L.raw_columns('vsig', 'TRB', vocab['TRB'], feature_definition=d)
    rows = np.zeros((10, len(columns)))
    rows[:, :6] = np.arange(10)[:, None]
    rows[:, 6] = 1
    rows[-1, 6:8] = [0, 1]
    rows[:, 8] = 1
    art = C.fit_matrices({'TRB': (rows, columns)}, vocab, sig='vsig', name='sparse',
                         n_components=3, feature_definition=d)
    fit = art.fits['TRB']
    assert fit.scale[6] == fit.scale[7] > 0
    assert fit.rotation[7].any(), 'rare but variable coordinate was erased'
    assert fit.scale[8] == fit.scale[9] == 0
    raw = dict(zip(columns, rows[-1]))
    expected = C.apply(raw, {}, art)
    saved = art.save(tmp_path / 'sparse')
    loaded = C.Corpus.load(saved)
    assert loaded.meta['feature_definition'] == d
    assert loaded.columns(named=True) == art.columns(named=True)
    actual = C.apply(raw, {}, loaded)
    for col in L.pc_columns('vsig', art.k):
        # Same float32 rotation contract as test_signature.py::test_the_artifact_round_trips.
        assert actual[col] == pytest.approx(expected[col], rel=0, abs=1e-4)
    absent = C.apply({}, {}, loaded)
    assert all(np.isnan(absent[col]) for col in L.pc_columns('vsig', art.k))
    loaded.meta['feature_definition']['classes'][0] = 'A'
    with pytest.raises(ValueError):
        loaded.verify()


def test_randomized_fit_records_total_variance_and_refuses_unreachable_fraction():
    import numpy as np
    from vdjtools.signature import corpus as C
    d = L.feature_definition('basic')
    columns = ['vsig:vj:TRB:' + str(i) for i in range(20)]
    rng = np.random.default_rng(4)
    x = rng.dirichlet(np.ones(20), size=30)
    fit = C.fit_locus(x, columns, feature_definition=d, n_components=2, solver='randomized', seed=7)
    assert fit.total_variance > fit.eigenvalues.sum()
    assert fit.variance_at(2) < 1
    with pytest.raises(ValueError):
        C.fit_locus(x, columns, feature_definition=d, n_components=.99, solver='randomized', max_components=2)
    with pytest.raises(ValueError):
        C.fit_locus(x, columns, feature_definition=d, mode='pcs')


def test_configured_raw_fit_apply_uses_frozen_vocabulary_and_policy(tmp_path):
    import numpy as np
    import polars as pl
    from vdjtools.signature import features as F, corpus as C
    from vdjtools.signature.signature import raw_features_cohort, vsig_cohort
    d = L.feature_definition('basic', motifs='imgt11')
    vocab = {'TRB': F.gene_vocab('TRB', feature_definition=d)}
    samples = {}
    for i in range(5):
        frame = pl.DataFrame({'junction_aa': ['CADGKF', 'CAGGQF', 'CKKGRF'],
                              'v_call': ['TRBV1*01', 'TRBV2*01', None],
                              'j_call': ['TRBJ1-1*01'] * 3, 'duplicate_count': [1+i, 2, 4-i]})
        samples[str(i)] = {'TRB': frame}
    table = raw_features_cohort(samples, vocab=vocab, feature_definition=d, cstar_target=.5, n_jobs=1)
    columns = L.raw_columns('vsig', 'TRB', vocab['TRB'], feature_definition=d)
    assert set(columns).issubset(table.columns)
    art = C.fit_cohort(samples, vocab=vocab, feature_definition=d, cstar_target=.5, n_components=2, n_jobs=1)
    result = vsig_cohort(samples, art, n_jobs=1)
    assert result.height == 5
    assert 'rotation_space' in result.columns
    assert result['rotation_space'].n_unique() == 1
    assert result['vsig:cov:TRB:cstar'].is_finite().all()
    with pytest.raises(ValueError):
        vsig_cohort(samples, art, weight='distinct')
    missing = vsig_cohort({'absent': {}}, art)
    assert np.isnan(missing['vsig:pc:TRB:PC01'][0])


def test_cli_feature_selection_and_description(tmp_path):
    import numpy as np
    import polars as pl
    from typer.testing import CliRunner
    from vdjtools.cli import app
    from vdjtools.signature import corpus as C
    runner = CliRunner()
    frame = pl.DataFrame({'junction_aa': ['CADGKF', 'CKKGQF'], 'v_call': ['TRBV1', 'TRBV2'],
                          'j_call': ['TRBJ1-1']*2, 'duplicate_count': [2, 3]})
    source = tmp_path / 'input.tsv'
    frame.write_csv(source, separator='\t')
    output = tmp_path / 'raw.tsv'
    result = runner.invoke(app, ['features', str(source), '--feature-set', 'basic', '--coverage', '.5',
                                 '--out', str(output)])
    assert result.exit_code == 0, (result.stdout, result.exception)
    table = pl.read_csv(output, separator='\t')
    assert any(':vj:' in c for c in table.columns)
    assert not any(':vjl:' in c for c in table.columns)
    bad = runner.invoke(app, ['features', str(source), '--feature-set', 'basic', '--v-context', 'gene'])
    assert bad.exit_code != 0
    d = L.feature_definition('joint')
    vocab = {'TRB': {'vjl': ['V|J|6', '_unknown']}}
    columns = L.raw_columns('vsig', 'TRB', vocab['TRB'], feature_definition=d)
    x = np.zeros((5, len(columns))); x[:,:6] = np.arange(5)[:,None]; x[:,6] = [0,.2,.4,.6,1]; x[:,7] = 1-x[:,6]
    art = C.fit_matrices({'TRB':(x,columns)},vocab,sig='vsig',name='tiny',feature_definition=d,n_components=2)
    path = art.save(tmp_path / 'joint')
    out = tmp_path / 'describe.tsv'
    result = runner.invoke(app,['signature','--corpus',str(path),'--describe','--named','all','--out',str(out)])
    assert result.exit_code == 0, (result.stdout, result.exception)
    result = runner.invoke(app,['signature',str(source),'--corpus',str(path),'--named','all',
                               '--out',str(tmp_path / 'projected.tsv')])
    assert result.exit_code == 0, (result.stdout, result.exception)
    projected = pl.read_csv(tmp_path / 'projected.tsv', separator='\t')
    assert 'rotation_space' in projected.columns
    assert pl.read_csv(out,separator='\t')['column'].to_list() == projected.columns[1:]


def test_configured_partial_holes_use_identical_fit_and_apply_center(tmp_path):
    import numpy as np
    from vdjtools.signature import corpus as C
    definition = L.feature_definition('basic')
    vocab = {'TRB': {'vj': ['V|J', '_unknown'], 'length': ['6', '7']}}
    columns = L.raw_columns('vsig', 'TRB', vocab['TRB'], feature_definition=definition)
    x = np.zeros((9, len(columns)))
    x[:, 0] = [0, 1, 2, 3, 4, 5, 20, np.nan, np.nan]
    x[:, 1] = np.arange(9)
    art = C.fit_matrices({'TRB': (x, columns)}, vocab, sig='vsig', name='holes',
                         feature_definition=definition, n_components=2)
    def scores(corpus):
        return np.array([[C.apply(dict(zip(columns, row)), {}, corpus)[f'vsig:pc:TRB:PC{i+1:02d}']
                          for i in range(2)] for row in x])
    projected = scores(art)
    assert np.allclose(np.median(projected, axis=0), 0, atol=1e-10)
    assert np.allclose(np.median(np.abs(projected), axis=0) * C.MAD_TO_SD, 1, atol=1e-10)
    assert np.allclose(projected, scores(C.Corpus.load(art.save(tmp_path/'holes'))), atol=1e-4)


def test_configured_synthetic_corpus_freezes_coverage():
    from vdjtools.signature import corpus as C
    from vdjtools.signature.signature import vsig, vsig_cohort
    d = L.feature_definition('basic')
    art, mats = C.synthesize('naive', loci=('TRG',), n_samples=8, size=40, seed=12,
                             n_components=2, feature_definition=d)
    assert set(art.meta['cstar_target']) == {'TRG'}
    assert art.meta['weight'] == 'log2p1'
    sample = C.synthesize('naive', loci=('TRG',), n_samples=1, size=40, seed=17, fit_corpus=False)[0]
    single = vsig(sample, art)
    batch = vsig_cohort({'sample': sample}, art)
    assert single['vsig:pc:TRG:PC01'] == batch['vsig:pc:TRG:PC01'][0]


def test_configured_coverage_own_is_refused_instead_of_silently_using_corpus():
    import numpy as np
    from vdjtools.signature import corpus as C
    from vdjtools.signature.signature import vsig_cohort
    d = L.feature_definition('basic')
    vocab = {'TRB': {'vj': ['V|J', '_unknown'], 'length': ['6']}}
    columns = L.raw_columns('vsig', 'TRB', vocab['TRB'], feature_definition=d)
    x = np.zeros((4, len(columns))); x[:, 0] = np.arange(4)
    art = C.fit_matrices({'TRB':(x, columns)}, vocab, sig='vsig', name='tiny', feature_definition=d,
                         meta={'cstar_target':{'TRB':.5}}, n_components=1)
    with pytest.raises(ValueError, match='shared coverage target'):
        vsig_cohort({'sample':{}}, art, cstar_target=None)


def test_cli_saved_vocabulary_matches_reference_feature_identity(tmp_path):
    import json
    import polars as pl
    from typer.testing import CliRunner
    from vdjtools.cli import app
    from vdjtools.signature import features as F
    from vdjtools.signature.signature import raw_features_cohort
    d = L.feature_definition('basic', motifs='imgt11')
    vocab = {'TRB': F.gene_vocab('TRB', feature_definition=d)}
    frame = pl.DataFrame({'junction_aa':['CADGKF','CKKGRF'], 'v_call':['TRBV1','TRBV2'],
                          'j_call':['TRBJ1-1']*2,'duplicate_count':[2,3]})
    source = tmp_path/'sample.tsv'; frame.write_csv(source,separator='\t')
    vocabulary = tmp_path/'vocabulary.json'; vocabulary.write_text(json.dumps(vocab))
    target = tmp_path/'raw.tsv'
    result = CliRunner().invoke(app,['features',str(source),'--feature-set','basic','--motifs','imgt11',
                                     '--coverage','.5','--feature-vocabulary',str(vocabulary),'--out',str(target)])
    assert result.exit_code == 0, (result.stdout,result.exception)
    expected = raw_features_cohort({'sample':{'TRB':frame}},vocab=vocab,feature_definition=d,cstar_target=.5)
    assert pl.read_csv(target,separator='\t')['embedding_space'][0] == expected['embedding_space'][0]


def test_wide_table_pca_memory_guard_uses_sample_rank(monkeypatch):
    import polars as pl
    from vdjtools.signature.corpus import TablePCA
    requests=[]
    monkeypatch.setattr('vdjtools.cores.require_memory',lambda size,reason:requests.append(size))
    table=pl.DataFrame({'embedding_space':['test']*3,
                       **{f'phi:{i}':[0.,float(i),1.] for i in range(100)}})
    rotation=TablePCA.fit(table,n_components=2)
    assert requests[0] == 32*3*100+16*3**2
    assert rotation.components.shape == (2,100)
