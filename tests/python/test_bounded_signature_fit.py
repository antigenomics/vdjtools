"""Bounded fitting must preserve the dense coordinate system without dense workspaces."""
import numpy as np
import pytest
from vdjtools.signature import corpus as C, layout as L


def matrix_fixture(n=48):
    rng = np.random.RandomState(19)
    names = ['vsig:div:TRB:1D_c', 'vsig:div:TRB:clonality']
    parts = [rng.normal(size=(n, 2))]
    for group, width in [('vj', 7), ('length', 26), ('motif_XXX', 51), ('motif_XX_X', 51), ('motif_X_XX', 51)]:
        names += [f'vsig:{group}:TRB:{j}' for j in range(width)]
        x = rng.exponential(size=(n, width))
        x[rng.uniform(size=x.shape) < .8] = 0
        x /= np.maximum(x.sum(1, keepdims=True), 1)
        parts.append(x)
    x = np.column_stack(parts)
    x[0] = np.nan
    x[1, 9:35] = np.nan
    x[2, 0] = np.nan
    return x, names


def test_sparse_raw_storage_preserves_zero_missing_and_chunk_boundary(monkeypatch):
    from vdjtools.signature import bounded as B
    x, columns = matrix_fixture()
    monkeypatch.setattr(B, 'CHUNK_BYTES', 4096)
    raw = B.BlockMatrix.from_dense(x, columns)
    restored = np.column_stack([part for _, _, part in raw.chunks()])
    np.testing.assert_array_equal(restored, x)
    assert raw.nbytes < x.nbytes
    assert max(part.nbytes for _, _, part in raw.chunks()) <= 4096


@pytest.mark.parametrize('configuration', ['basic', 'joint'])
def test_bounded_corpus_fit_matches_dense_reference(configuration, monkeypatch, tmp_path):
    from vdjtools.signature import bounded as B
    x, columns = matrix_fixture()
    if configuration == 'joint':
        keep = [i for i, c in enumerate(columns) if ':length:' not in c]
        x = x[:, keep]
        columns = [columns[i].replace(':vj:', ':vjl:') for i in keep]
    d = L.feature_definition(configuration, motifs='imgt11')
    monkeypatch.setattr(B, 'CHUNK_BYTES', 4096)
    kw = dict(feature_definition=d, n_components=5, solver='randomized', seed=27)
    reference = C.fit_locus(x[1:], columns, bounded=False, **kw)
    actual = C.fit_locus(B.BlockMatrix.from_dense(x[1:], columns), columns, **kw)
    for name in ['loc', 'scale', 'pca_offset', 'total_variance', 'eigenvalues']:
        np.testing.assert_allclose(getattr(actual, name), getattr(reference, name), atol=1e-12, rtol=1e-12)
    for q in C.WINSOR_PS:
        for a, b in zip(actual.bounds[q], reference.bounds[q]):
            np.testing.assert_allclose(a, b, atol=1e-12, rtol=0)
    np.testing.assert_array_equal(actual.n_obs, reference.n_obs)
    np.testing.assert_allclose(actual.rotation, reference.rotation, atol=1e-5, rtol=0)
    np.testing.assert_allclose(actual.pc_loc, reference.pc_loc, atol=1e-5, rtol=0)
    np.testing.assert_allclose(actual.pc_scale, reference.pc_scale, atol=1e-5, rtol=0)
    # The oracle is also independent of the requested chunk size.
    monkeypatch.setattr(B, 'CHUNK_BYTES', 8192)
    again = C.fit_locus(x[1:], columns, **kw)
    np.testing.assert_allclose(again.rotation, actual.rotation, atol=1e-5, rtol=0)


def test_operator_randomized_svd_matches_existing_solver():
    from vdjtools.signature.bounded import randomized_svd
    from scipy.sparse.linalg import aslinearoperator
    from sklearn.utils.extmath import randomized_svd as oracle
    rng = np.random.RandomState(31)
    for shape in [(28, 73), (73, 28)]:
        x = rng.normal(size=shape)
        u, s, vt = randomized_svd(aslinearoperator(x), 7, seed=81)
        a, b, c = oracle(x, n_components=7, n_iter=4, random_state=81)
        np.testing.assert_allclose(s, b, atol=1e-12, rtol=1e-12)
        np.testing.assert_allclose(u @ np.diag(s) @ vt, a @ np.diag(b) @ c, atol=1e-12, rtol=1e-12)


def test_configured_cohort_builds_sparse_without_a_dense_sample_matrix(monkeypatch):
    import polars as pl
    from vdjtools.signature import features as F
    from vdjtools.signature.signature import vsig_cohort
    d = L.feature_definition('basic', motifs='imgt11')
    vocab = {'TRB': F.gene_vocab('TRB', feature_definition=d)}
    samples = {str(i): {'TRB': pl.DataFrame({'junction_aa': ['CADGKF', 'CAGGQF', 'CKKGRF'],
               'v_call': ['TRBV1*01', 'TRBV2*01', None], 'j_call': ['TRBJ1-1']*3,
               'duplicate_count': [i+1, 2, 8-i]})} for i in range(8)}
    kw = dict(vocab=vocab, feature_definition=d, cstar_target=.5, n_components=2,
              solver='randomized', seed=3, threads=1)
    oracle = C.fit_cohort(samples, bounded=False, **kw)
    real_full = np.full
    def no_dense(shape, *args, **kwargs):
        if isinstance(shape, tuple) and len(shape) == 2:
            assert shape != (8, len(oracle.fits['TRB'].columns)), 'dense cohort allocation'
        return real_full(shape, *args, **kwargs)
    monkeypatch.setattr(np, 'full', no_dense)
    got = C.fit_cohort(samples, **kw)
    a, b = vsig_cohort(samples, oracle), vsig_cohort(samples, got)
    signs = np.sign(np.sum(oracle.fits['TRB'].rotation * got.fits['TRB'].rotation, axis=0))
    for c in got.columns():
        sign = signs[int(c.rsplit('PC', 1)[1])-1] if c.startswith('vsig:pc:TRB:PC') else 1
        np.testing.assert_allclose(a[c], sign * b[c], rtol=0, atol=1e-4, equal_nan=True)
    assert oracle.meta['trimmed_frac'] == pytest.approx(got.meta['trimmed_frac'])
    # Sparse assembly also preserves partial holes and absent loci without zero-imputing them.
    from vdjtools.cores import available_cores
    parallel = C.fit_cohort(samples, n_jobs=min(2, available_cores()), **kw)
    np.testing.assert_array_equal(got.fits['TRB'].rotation, parallel.fits['TRB'].rotation)


def test_sparse_table_pca_matches_dense_standard_scaler():
    import polars as pl
    from vdjtools.signature import features as F
    from vdjtools.signature.signature import raw_features_cohort
    d = L.feature_definition('joint', motifs='imgt11')
    vocab = {'TRB': F.gene_vocab('TRB', feature_definition=d)}
    samples = [(str(i), lambda i=i: {'TRB': pl.DataFrame({
        'junction_aa': ['CADGKF', 'CKKGQF', 'CAGGRF'], 'v_call': ['TRBV1', 'TRBV2', None],
        'j_call': ['TRBJ1-1']*3, 'duplicate_count': [i+1, 3, 9-i]})}) for i in range(9)]
    table = raw_features_cohort(samples, vocab=vocab, feature_definition=d, cstar_target=.5)
    columns = L.raw_columns('vsig', 'TRB', vocab['TRB'], feature_definition=d)
    reference = C.TablePCA.fit(table, columns=columns, n_components=2, seed=9)
    actual = C.TablePCA.fit_features(samples, vocab=vocab, feature_definition=d,
                                    cstar_target=.5, n_components=2, seed=9)
    assert actual.columns == reference.columns and actual.space == reference.space
    for name in ['fill', 'mean', 'scale', 'center', 'explained_variance_ratio']:
        np.testing.assert_allclose(getattr(actual, name), getattr(reference, name), atol=1e-12, rtol=1e-12)
    a, b = reference.transform(table), actual.transform(table)
    signs = np.sign(np.sum(actual.components * reference.components, axis=1))
    np.testing.assert_allclose(b.select(['pc0', 'pc1']).to_numpy() * signs,
                               a.select(['pc0', 'pc1']).to_numpy(), atol=1e-5, rtol=0)


def test_sparse_quantiles_include_implicit_zeros_and_exclude_holes():
    from vdjtools.signature import bounded as B
    rng = np.random.RandomState(93)
    x = rng.normal(size=(53, 27))
    x[rng.uniform(size=x.shape) < .8] = 0
    x[0] = np.nan
    x[1, :4] = np.nan
    x[:, -1] = np.nan
    raw = B.BlockMatrix.from_dense(x, ['vsig:vj:TRB:'+str(i) for i in range(27)])
    for block in raw.blocks:
        count, quantiles = B.sparse_quantiles(block, [.01, .5, .99])
        with pytest.warns(RuntimeWarning, match='All-NaN'):
            expected = np.nanquantile(x, [.01, .5, .99], axis=0)
        np.testing.assert_array_equal(count, np.isfinite(x).sum(0))
        np.testing.assert_allclose(quantiles, expected, atol=1e-15, rtol=0, equal_nan=True)


@pytest.mark.parametrize('configuration,motifs,context', [
    (config, motifs, ctx) for config in ['basic', 'joint']
    for motifs, ctx in [('none', 'none')] + [(m, v) for m in ['imgt11', 'blosum62']
                                           for v in ['none', 'gene', 'cluster']]])
def test_all_fourteen_definitions_sparse_dense_equivalence(configuration, motifs, context):
    import itertools
    import json
    from pathlib import Path
    import polars as pl
    from vdjtools.signature.signature import vsig_cohort
    maps = json.loads((Path(C.__file__).parent.parent/'resources/signature/motif_definition_v1.json').read_text())
    vmap = {'TRB': {'version': 1, 'mapping': {'TRBV1': 'v000', 'TRBV2': 'v001'}}}
    d = L.feature_definition(configuration, motifs=motifs, v_context=context,
            motif_map=maps['motif_map'] if motifs == 'blosum62' else None,
            v_maps=vmap if context == 'cluster' else None)
    v = {'vus': ['TRBV1', 'TRBV2'], 'jus': ['TRBJ1-1']}
    lengths = [str(n) for n in L.SPECTRATYPE_LENGTHS]
    if configuration == 'basic':
        v.update(vj=['TRBV1|TRBJ1-1', 'TRBV2|TRBJ1-1', '_unknown'], length=lengths)
    else:
        v['vjl'] = [f'{gene}|TRBJ1-1|{n}' for gene in v['vus'] for n in lengths] + ['_unknown']
    if motifs != 'none':
        tokens = ([''.join(t) for t in itertools.product('0123456789a', repeat=3)] if motifs == 'imgt11'
                  else sorted(set(maps['motif_map']['mapping'].values())))
        contexts = [''] if context == 'none' else (v['vus'] if context == 'gene' else ['v000', 'v001'])+['_unknown']
        names = [f'{ctx}|{t}' if ctx else t for ctx in contexts for t in tokens]
        v.update({f'motif_{shape}': names for shape in L.MOTIF_SHAPES})
    rng = np.random.RandomState(13)
    samples = [(str(i), {'TRB': pl.DataFrame({'junction_aa': ['CADGKF', 'CKKGQF', 'CAGGRF', 'CASSIRSSYEQYF'],
               'v_call': ['TRBV1', 'TRBV2', None, 'TRBV2'], 'j_call': ['TRBJ1-1']*4,
               'duplicate_count': rng.randint(1, 40, 4)})} if i else {}) for i in range(13)]
    kw = dict(vocab={'TRB': v}, feature_definition=d, cstar_target=.5, n_components=3,
              n_jobs=1, threads=1, solver='randomized', seed=81)
    dense = C.fit_cohort(samples, bounded=False, **kw)
    bounded = C.fit_cohort(samples, **kw)
    a, b = dense.fits['TRB'], bounded.fits['TRB']
    for name in ['loc', 'scale', 'pca_offset', 'eigenvalues', 'total_variance']:
        np.testing.assert_allclose(getattr(a, name), getattr(b, name), atol=1e-12, rtol=1e-12)
    signs = np.sign(np.sum(a.rotation*b.rotation, axis=0))
    np.testing.assert_allclose(a.rotation, b.rotation*signs, atol=1e-5, rtol=0)
    expected, actual = vsig_cohort(samples, dense), vsig_cohort(samples, bounded)
    for c in bounded.columns():
        sign = signs[int(c.rsplit('PC', 1)[1])-1] if c.startswith('vsig:pc:TRB:PC') else 1
        np.testing.assert_allclose(expected[c], actual[c]*sign, atol=1e-4, rtol=0, equal_nan=True)


@pytest.mark.parametrize('dtype', [int, np.float32, np.float64])
@pytest.mark.parametrize('rows', [2, 4])
def test_diversity_clipped_implicit_zero_and_input_dtype(dtype, rows):
    columns = ['vsig:div:TRB:clonality', 'vsig:vj:TRB:a', 'vsig:vj:TRB:b']
    x = np.array([[0, 0, 1], [1, 1, 0], [0, 1, 0], [1, 0, 1]], dtype=dtype)[:rows]
    kw = dict(feature_definition=L.feature_definition('basic'), n_components=2,
              solver='randomized', seed=41)
    expected = C.fit_locus(x.astype(float), columns, bounded=False, **kw)
    actual = C.fit_locus(x, columns, **kw)
    for name in ['loc', 'scale', 'pca_offset', 'eigenvalues']:
        np.testing.assert_allclose(getattr(actual, name), getattr(expected, name), atol=1e-12, rtol=1e-12)


def test_table_pca_column_checks_are_linear(monkeypatch):
    import polars as pl
    rng = np.random.RandomState(78)
    names = [f'x{i}' for i in range(100)]
    table = pl.DataFrame(rng.normal(size=(5, 100)), schema=names).with_columns(
        pl.Series('sample_id', [str(i) for i in range(5)]), pl.lit('test').alias('embedding_space'))
    calls = []
    original = pl.DataFrame.columns.fget
    monkeypatch.setattr(pl.DataFrame, 'columns', property(lambda frame: (calls.append(1), original(frame))[1]))
    fitted = C.TablePCA.fit(table, columns=names, n_components=2)
    assert len(calls) <= 5, 'column validation rereads the whole schema per feature'
    class TrackedColumns(list):
        def __contains__(self, key):
            raise AssertionError('metadata filtering must not scan selected columns per feature')
    fitted.columns = TrackedColumns(fitted.columns)
    result = fitted.transform(table)
    assert result.columns == ['sample_id', 'embedding_space', 'pc0', 'pc1']


@pytest.mark.parametrize('shape,k', [((200, 700), 6), ((45, 80), 2), ((1050, 32), 5), ((80, 550), 70)])
def test_raw_table_operator_matches_pca_auto_solver_regimes(shape, k):
    from vdjtools.signature import bounded as B
    x = np.random.RandomState(46).normal(size=shape)
    x[np.random.RandomState(47).uniform(size=shape) < .6] = 0
    x[0, :7] = np.nan
    fill = np.nanmedian(x, axis=0)
    reference = C.fit_pca(np.where(np.isfinite(x), x, fill), k, seed=23)
    scaler, pca = reference.steps[0][1], reference.steps[1][1]
    raw = B.BlockMatrix.from_dense(x, [f'vsig:vj:TRB:{i}' for i in range(shape[1])])
    got = B.table_pca(raw, k, seed=23)
    for actual, expected in zip(got[:4], [fill, scaler.mean_, scaler.scale_, pca.mean_]):
        np.testing.assert_allclose(actual, expected, atol=1e-12, rtol=1e-12)
    np.testing.assert_allclose(got[5], pca.explained_variance_ratio_, atol=1e-12, rtol=1e-12)
    np.testing.assert_allclose(got[4], pca.components_, atol=1e-5, rtol=0)


@pytest.mark.parametrize('shape,k', [((45, 80), 2), ((200, 700), 6), ((1050, 32), 5)])
def test_raw_table_preserves_older_sklearn_sign_rules(monkeypatch, shape, k):
    import sklearn
    from sklearn.utils.extmath import randomized_svd, svd_flip
    from sklearn.preprocessing import StandardScaler
    from vdjtools.signature import bounded as B
    monkeypatch.setattr(sklearn, '__version__', '1.4.2')
    x = np.random.RandomState(91).normal(size=shape)
    z = StandardScaler().fit_transform(x)
    z -= z.mean(0)
    if max(shape) <= 500:
        u, _, vt = np.linalg.svd(z, full_matrices=False)
        _, expected = svd_flip(u, vt)
    else:
        _, _, expected = randomized_svd(z, k, random_state=23, n_iter='auto')
    raw = B.BlockMatrix.from_dense(x, [f'vsig:vj:TRB:{i}' for i in range(shape[1])])
    got = B.table_pca(raw, k, seed=23)
    np.testing.assert_allclose(got[4], expected[:k], atol=1e-5, rtol=0)


def test_tall_covariance_uses_bounded_column_products(monkeypatch):
    from vdjtools.signature import bounded as B
    x = np.random.RandomState(92).normal(size=(2000, 40))
    raw = B.BlockMatrix.from_dense(x, [f'vsig:vj:TRB:{i}' for i in range(40)])
    monkeypatch.setattr(B, 'CHUNK_BYTES', 32768)
    original, widths = B.CenteredOperator._matmat, []
    def checked(self, q):
        widths.append(q.shape[1])
        assert q.shape[1] <= B._width(x.shape[0])
        return original(self, q)
    monkeypatch.setattr(B.CenteredOperator, '_matmat', checked)
    got = B.table_pca(raw, 3, seed=23)
    expected = C.fit_pca(x, 3, seed=23).steps[1][1]
    np.testing.assert_allclose(got[4], expected.components_, atol=1e-5, rtol=0)
    assert len(widths) == 20


def test_remove_unobserved_rows_shares_sparse_values_without_mutation():
    from vdjtools.signature import bounded as B
    x, columns = matrix_fixture()
    x[[9, 19]] = np.nan
    raw = B.BlockMatrix.from_dense(x, columns)
    observed = raw.observed()
    selected = raw.take_rows(observed)
    np.testing.assert_array_equal(np.column_stack([v for _, _, v in selected.chunks()]), x[observed])
    for before, after in zip(raw.blocks, selected.blocks):
        if before.values.nnz:
            assert np.shares_memory(before.values.data, after.values.data)
        np.testing.assert_array_equal(before.values.toarray()[observed], after.values.toarray())
    np.testing.assert_array_equal(np.column_stack([v for _, _, v in raw.chunks()]), x)
    # Selecting observed rows must still remove their actual entries correctly.
    keep = np.arange(len(x)) % 2 == 0
    subset = raw.take_rows(keep)
    np.testing.assert_array_equal(np.column_stack([v for _, _, v in subset.chunks()]), x[keep])


def test_missing_row_remap_bounds_numpy_index_cast(monkeypatch):
    from vdjtools.signature import bounded as B
    x, columns = matrix_fixture()
    raw = B.BlockMatrix.from_dense(x, columns)
    take = np.take
    sizes = []
    def bounded_take(a, indices, **kw):
        sizes.append(indices.size)
        return take(a, indices, **kw)
    monkeypatch.setattr(B, 'CHUNK_BYTES', 64)
    monkeypatch.setattr(np, 'take', bounded_take)
    raw.take_rows(raw.observed())
    assert sizes and max(sizes) <= 4


@pytest.mark.parametrize('has_trim', [False, True])
def test_sparse_assembly_releases_freed_pages_at_phase_boundary(monkeypatch, has_trim):
    from types import SimpleNamespace
    from vdjtools.signature import bounded as B, cohort
    from vdjtools.signature.features import gene_vocab
    events = []
    finish, parallel = B.MatrixBuilder.finish, cohort.parallel_rows
    def joined(*a, **kw):
        result = parallel(*a, **kw)
        events.append('joined')
        return result
    def assembled(self):
        result = finish(self)
        events.append('assembled')
        return result
    def trim(pad):
        assert pad == 0 and events[-1] == 'assembled'
        events.append('released')
    monkeypatch.setattr(cohort, 'parallel_rows', joined)
    monkeypatch.setattr(B.MatrixBuilder, 'finish', assembled)
    monkeypatch.setattr(B, 'sys', SimpleNamespace(platform='linux'), raising=False)
    monkeypatch.setattr(B, 'CDLL', lambda path: SimpleNamespace(**({'malloc_trim': trim} if has_trim else {})), raising=False)
    definition = L.feature_definition('basic')
    vocab = {'TRB': gene_vocab('TRB', feature_definition=definition)}
    result = B.feature_matrices([], vocab, definition, {}, 1)
    assert result['TRB'][0].shape[0] == 0
    assert events == ['joined', 'assembled'] + (['released'] if has_trim else [])


def test_cohort_releases_original_sparse_indices_before_fitting(monkeypatch):
    import weakref
    from vdjtools.signature import bounded as B
    original_indices = []
    x, columns = matrix_fixture(8)
    def features(*a, **kw):
        raw = B.BlockMatrix.from_dense(x, columns)
        original_indices.extend(weakref.ref(b.values.indices) for b in raw.blocks)
        return {'TRB': (raw, columns)}
    def fit(mats, *a, **kw):
        selected = mats['TRB'][0]
        assert selected.shape[0] == 7
        assert all(ref() is None for ref in original_indices)
        np.testing.assert_array_equal(np.column_stack([v for _, _, v in selected.chunks()]), x[1:])
        return 'fitted'
    monkeypatch.setattr(B, 'feature_matrices', features)
    monkeypatch.setattr(C, 'fit_matrices', fit)
    assert C.fit_cohort([('sample', {})], vocab={'TRB': {}}, feature_definition=L.feature_definition('basic'),
                        cstar_target=.5, solver='randomized') == 'fitted'


def test_centered_operator_vector_hooks_accept_column_vectors():
    from scipy import sparse
    from vdjtools.signature.bounded import CenteredOperator
    x = np.array([[1., 0., 2.], [0., 3., 0.]])
    baseline = np.array([.1, .2, .3])
    offset = np.array([.4, .5, .6])
    op = CenteredOperator([(sparse.csr_matrix(x), np.ones(2), baseline)], offset)
    dense = x + baseline - offset
    for hook, matrix in [(op._matvec, dense), (op._rmatvec, dense.T)]:
        q = np.arange(matrix.shape[1], dtype=float)
        for shape in [(len(q),), (len(q), 1)]:
            np.testing.assert_allclose(hook(q.reshape(shape)), matrix @ q)
