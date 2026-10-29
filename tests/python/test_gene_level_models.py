"""Working models fit genes once and ignore allele suffixes at prediction."""
import importlib

import polars as pl
import pytest

from vdjtools.model import data, native, reference
from vdjtools.model.io import from_arda, from_germline
from vdjtools.model.generate import generate


@pytest.mark.parametrize('status', ['truncated', 'partial'])
def test_arda_scaffold_excludes_incomplete_germline_templates(monkeypatch, toy_germline, status):
    complete = toy_germline.with_columns(pl.lit('ok').alias('status'))
    partial = complete.filter(pl.col('allele') == 'TOYV1*01').with_columns(
        pl.lit('TOYV1*99').alias('allele'), pl.lit('TGT').alias('sequence'),
        pl.lit(status).alias('status'))
    monkeypatch.setattr(reference, 'load_germline', lambda *a: pl.concat([complete, partial]))
    monkeypatch.setattr(reference, 'arda_full_germline', lambda *a: {})
    for gene_level in (False, True):
        model = from_arda('TOY', gene_level=gene_level)
        assert 'TOYV1*99' not in model.genomic['genes_v']['v_allele'].to_list()


def test_gene_level_scaffold_prefers_01_then_lowest_usable(monkeypatch, toy_germline):
    gl = toy_germline.with_columns(pl.when(pl.col('allele') == 'TOYV2*01')
                                  .then(pl.lit('TOYV2*03')).otherwise(pl.col('allele')).alias('allele'))
    extra = gl.filter(pl.col('allele') == 'TOYV2*03').with_columns(pl.lit('TOYV2*02').alias('allele'))
    invalid = gl.filter(pl.col('allele') == 'TOYV2*03').with_columns(
        pl.lit('TOYV2*01').alias('allele'), pl.lit('NNN').alias('sequence'))
    monkeypatch.setattr(reference, 'load_germline', lambda *a: pl.concat([gl, extra, invalid]))
    monkeypatch.setattr(reference, 'arda_full_germline', lambda *a: {})
    m = from_arda('TOY', gene_level=True)
    assert set(m.genomic['genes_v']['v_allele']) == {'TOYV1*01', 'TOYV2*02'}
    for frame in m.genomic.values():
        assert frame['gene'].n_unique() == frame.height


def test_gene_level_clonotypes_merge_suffixes_and_preserve_counts():
    frame = pl.DataFrame({'v_call': ['TOYV1*01,TOYV2*01', 'TOYV2*99,TOYV1*03'],
                          'j_call': ['TOYJ1*01', 'TOYJ1*99'], 'junction': ['TGTGCATTT']*2,
                          'junction_aa': ['CAF']*2, 'locus': ['TOY']*2,
                          'd_call': ['TOYD1*01', 'TOYD1*99'], 'd2_call': [None]*2,
                          'count': [3, 5]})
    out = data.unique_clonotypes(frame, gene_level=True)
    assert out.height == 1 and out['count'][0] == 8
    assert out['v_call'][0] == 'TOYV1,TOYV2' and out['j_call'][0] == 'TOYJ1'
    assert out['d_call'][0] == 'TOYD1'


def test_default_raw_builder_fits_gene_level_once(monkeypatch, toy_germline):
    m = from_germline(toy_germline.unique('gene', keep='first'), locus='TOY', ins_max=3)
    io = importlib.import_module('vdjtools.model.io')
    infer = importlib.import_module('vdjtools.model.infer')
    seen = []
    def scaffold(*a, **kw):
        assert kw['gene_level'] is True
        return m
    def fit(model, seqs, **kw):
        seen.extend(seqs)
        assert all(len(mask[0]) == len(mask[1]) == 1 for mask in kw['masks'])
        return model, infer.InferenceReport()
    monkeypatch.setattr(io, 'from_arda', scaffold)
    monkeypatch.setattr(infer, 'infer_native', fit)
    monkeypatch.setattr(infer, 'candidate_masks', lambda *a, **kw:
                        ([(['TOYV1*01'], ['TOYJ1*01'], None)], {}))
    frame = pl.DataFrame({'v_call': ['TOYV1*01', 'TOYV1*99'], 'j_call': ['TOYJ1*01', 'TOYJ1*77'],
                          'junction': ['TGTGCATTT']*2, 'junction_aa': ['CAF']*2,
                          'locus': ['TOY']*2, 'd_call': [None]*2, 'd2_call': [None]*2})
    data.build_model('TOY', clones=frame, threads=1)
    assert seen == ['TGTGCATTT']


@pytest.mark.parametrize('kind', ['nt', 'aa', 'degenerate'])
def test_gene_level_predictions_ignore_suffixes(toy_germline, kind):
    m = from_germline(toy_germline.unique('gene', keep='first'), locus='TOY', ins_max=3)
    draw = generate(m, 1, seed=7, productive_only=True, engine='native', threads=1).row(0, named=True)
    v, j = draw['v_call'].split('*')[0], draw['j_call'].split('*')[0]
    sequence = draw['junction_nt'] if kind == 'nt' else draw['junction_aa']
    if kind == 'degenerate':
        sequence = list(sequence)
    scalar = getattr(native, 'pgen_'+kind if kind != 'degenerate' else 'pgen_aa_degenerate')
    batch = getattr(native, 'pgen_'+kind+'_batch' if kind != 'degenerate' else 'pgen_aa_degenerate_batch')
    expected = scalar(m, sequence, v+'*01', j+'*01')
    assert expected > 0
    assert scalar(m, sequence, v, j) == expected
    assert scalar(m, sequence, v+'*99', j+'*42') == expected
    assert batch(m, [sequence]*2, [v, v+'*99'], [j, j+'*42'], threads=1) == [expected]*2
    if kind != 'degenerate':
        pgen = importlib.import_module('vdjtools.model.pgen')
        prep = pgen.prepare(m)
        fn = getattr(pgen, 'pgen_'+kind)
        assert fn(prep, sequence, v+'*99', j+'*42') == pytest.approx(expected)


def test_gene_level_filter_retains_ambiguous_known_genes(toy_germline):
    m = from_germline(toy_germline.unique('gene', keep='first'), locus='TOY', ins_max=3)
    frame = pl.DataFrame({'v_call': ['TOYV1,TOYV2', 'UNKNOWN,TOYV1'],
                          'j_call': ['TOYJ1', 'UNKNOWN,TOYJ1'],
                          'junction': ['TGTGCATTT']*2})
    assert data._filter_for_em(frame, m).height == 2


def test_prepared_gene_clonotypes_keep_counts_and_d_ambiguity():
    frame = pl.DataFrame({'v_call': ['TOYV1*01', 'TOYV1*99'],
                          'j_call': ['TOYJ1*01', 'TOYJ1*99'], 'junction': ['TGTGCATTT']*2,
                          'd_call': ['TOYD1*01', 'TOYD2*01'],
                          'd2_call': [None, 'TOYD3*01'], 'count': [3, 5]})
    out = data.unique_clonotypes(frame, gene_level=True)
    assert out.height == 1 and out['count'][0] == 8
    assert out['d_call'][0] == 'TOYD1,TOYD2'
    assert out['d2_call'][0] == 'TOYD3'


def test_gene_aliases_in_search_boundary_and_stitch(toy_germline):
    from vdjtools.model.viterbi import best_scenario, infer_nt, infer_nt_batch
    from vdjtools.model.pgen import prepare
    from vdjtools.model.boundary import germline_boundary
    from vdjtools.model.stitch import stitch_contig

    gl = toy_germline.unique('gene', keep='first').with_columns(
        full_germline=pl.concat_str(pl.lit('AAA'), 'sequence', pl.lit('AAA')),
        cdr3_anchor=pl.when(pl.col('segment') == 'V').then(3)
        .otherwise(pl.col('sequence').str.len_bytes()))
    m = from_germline(gl, locus='TOY', ins_max=3)
    draw = generate(m, 1, seed=7, productive_only=True, engine='native', threads=1).row(0, named=True)
    v, j, aa = draw['v_call'], draw['j_call'], draw['junction_aa']
    av, aj = v.split('*')[0]+'*99', j.split('*')[0]+'*99'
    assert infer_nt(m, aa, av, aj) == infer_nt(m, aa, v, j)
    assert infer_nt(prepare(m), aa, av, aj) == infer_nt(prepare(m), aa, v, j)
    assert infer_nt_batch(m, [aa], [av], [aj], threads=1).equals(
        infer_nt_batch(m, [aa], [v], [j], threads=1))
    assert germline_boundary(m, [aa], [av], [aj]).equals(germline_boundary(m, [aa], [v], [j]))
    expected = stitch_contig(m, v, j, draw['junction_nt'])
    assert expected is not None
    assert stitch_contig(m, av, aj, draw['junction_nt']) == expected
    assert best_scenario(m, draw['junction_nt'], av, aj) == best_scenario(m, draw['junction_nt'], v, j)


def test_multiallele_models_keep_explicit_alleles_and_reject_unknowns(toy_model):
    for call in ('TOYV1', 'TOYV1*99', 'UNKNOWN*99'):
        with pytest.raises(KeyError):
            native.pgen_aa(toy_model, 'CAF', call, 'TOYJ1*01')


@pytest.mark.parametrize('organism,locus', [('human', l) for l in ('TRA','TRB','TRG','TRD','IGH','IGK','IGL')]
                         + [('mouse', l) for l in ('TRA','TRB','IGH','IGK','IGL')])
def test_supported_gene_scaffolds(organism, locus):
    from vdjtools.model.check import check_model

    m = from_arda(locus, organism, gene_level=True)
    check_model(m, germline='none', raise_on='error')
    gl = reference.load_germline(locus, organism)
    invalid = reference.validate_germline(gl).filter(pl.col('check').is_in(
        ['germline_empty_sequence', 'germline_ambiguous', 'germline_anchor_frame']))['allele']
    usable = gl.filter(~pl.col('allele').is_in(invalid.implode()) &
                       (pl.col('status').fill_null('') != 'truncated'))
    for seg, frame in m.genomic.items():
        assert frame['gene'].n_unique() == frame.height
        for gene, allele in frame.select('gene', seg.replace('genes_', '')+'_allele').iter_rows():
            options = usable.filter(pl.col('gene') == gene)['allele'].to_list()
            assert allele == min(options, key=lambda a: (int(a.split('*')[1]), a))


def test_gene_level_bundle_preserves_selected_reference_label(tmp_path, monkeypatch, toy_germline):
    from vdjtools.model import bundled

    gl = toy_germline.unique('gene', keep='first').with_columns(
        pl.col('allele').str.replace(r'\*01$', '*02'))
    m = from_germline(gl, locus='TOY', ins_max=3)
    m.save(tmp_path / '_bundled' / 'learned' / 'TOY')
    monkeypatch.setattr(bundled, 'files', lambda _: tmp_path)
    loaded = bundled.load_bundled('TOY')
    assert loaded.genomic['genes_v'].equals(m.genomic['genes_v'])
    for name, table in m.tables.items():
        assert loaded.tables[name].equals(table), name


def test_raw_prepare_merges_gene_calls_without_losing_d_alternatives(monkeypatch, tmp_path):
    frame = pl.DataFrame({'v_call': ['TOYV1*01']*2, 'j_call': ['TOYJ1*01']*2,
                          'junction': ['TGTGCATTT']*2, 'd_call': ['TOYD1*01', 'TOYD2*01'],
                          'd2_call': [None, 'TOYD3*01'], 'duplicate_count': [3, 5]})
    monkeypatch.setattr(data, 'fetch_fastq', lambda *a: tmp_path / 'raw.fq.gz')
    monkeypatch.setattr(data, 'annotate_reads', lambda *a, **kw: frame)
    out = data.prepare('human', 'TRB', 'nonfunctional', out_dir=tmp_path, gene_level=True)
    assert out['count'].to_list() == [8]
    assert out['d_call'].to_list() == ['TOYD1,TOYD2']
    assert out['d2_call'].to_list() == ['TOYD3']
