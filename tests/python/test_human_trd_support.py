"""Human TRDD order follows mapped forward genomic positions, not numeric guessing."""
from dataclasses import replace

import polars as pl
import pytest

from vdjtools.model import Model, native, pgen
from vdjtools.model.dd import to_dd, to_history
from vdjtools.model.io import from_germline
from vdjtools.model.reference import forbidden_tandem_pairs


FORWARD = {('TRDD1*01', 'TRDD2*01'), ('TRDD1*01', 'TRDD3*01'),
           ('TRDD2*01', 'TRDD3*01')}


def trd_template():
    entries = [('V', 'TRDV1*01', 'TGT'), ('D', 'TRDD1*01', 'AAA'),
               ('D', 'TRDD2*01', 'CCC'), ('D', 'TRDD3*01', 'GGG'),
               ('J', 'TRDJ1*01', 'TTT')]
    gl = pl.DataFrame([dict(segment=s, allele=a, gene=a.split('*')[0], sequence=q,
        functional=True, cdr3_anchor=-1, full_germline='') for s,a,q in entries])
    model = from_germline(gl, locus='TRD', organism='human', ins_max=0,
                         palindrome_max={'v_3':0, 'j_5':0, 'd_5':0, 'd_3':0})
    tables = dict(model.tables)
    for name in ('v_3_del', 'j_5_del'):
        tables[name] = tables[name].with_columns(p=(pl.col('ndel') == 0).cast(pl.Float64))
    for name in ('d_del',):
        tables[name] = tables[name].with_columns(
            p=((pl.col('ndel5') == 0) & (pl.col('ndel3') == 0)).cast(pl.Float64))
    return Model(manifest=model.manifest, tables=tables, genomic=model.genomic).validate()


def test_human_map_preserves_unknown_nonhuman_and_custom_support():
    alleles = ['TRDD1*02', 'TRDD2*01', 'TRDD3*01', 'TRDD4*01']
    bad = forbidden_tandem_pairs(alleles, 'TRD', 'human')
    assert ('TRDD2*01', 'TRDD1*02') in bad
    assert ('TRDD3*01', 'TRDD2*01') in bad
    assert ('TRDD1*02', 'TRDD3*01') not in bad
    assert ('TRDD4*01', 'TRDD1*02') not in bad
    for organism in (None, 'mouse', 'custom'):
        assert ('TRDD2*01', 'TRDD1*02') not in forbidden_tandem_pairs(alleles, 'TRD', organism)
    assert ('TRDD2*01', 'TRDD1*02') not in forbidden_tandem_pairs(alleles, 'CUSTOM', 'human')


def test_scaffold_and_uniform_initialization_normalize_on_human_support():
    from vdjtools.model.infer import _uniform_init
    model = to_dd(trd_template())
    for table in (model.tables['d2_gene'], _uniform_init(model)['d2_gene']):
        assert set(table.filter(pl.col('p') > 0).select('d_allele', 'd2_allele').iter_rows()) == FORWARD
        assert table.filter(pl.col('d_allele') == 'TRDD3*01')['p'].sum() == 0
    joint = to_history(trd_template(), p_nd2=.5)
    h = joint.tables['d_history']
    assert set(h.filter(pl.col('n_d') == 2).select('d_allele', 'd2_allele').iter_rows()) == FORWARD
    assert h['p'].sum() == pytest.approx(1)
    assert h.filter(pl.col('n_d') == 2)['p'].sum() == pytest.approx(.5)


@pytest.mark.parametrize('engine', ['reference', 'native'])
def test_generation_pgen_and_em_share_human_forward_support(engine):
    from vdjtools.model.generate import generate
    from vdjtools.model.infer import infer, infer_native
    model = to_history(trd_template(), p_nd2=.5)
    module = native if engine == 'native' else pgen
    prepared = model if engine == 'native' else pgen.prepare(model)
    assert module.pgen_nt(prepared, 'TGTCCCAAATTT') == 0
    assert module.pgen_nt(prepared, 'TGTAAACCCTTT') > 0
    draws = generate(model, 64, seed=9, engine=engine, productive_only=False, threads=1)
    pairs = draws.filter(pl.col('d2_call').is_not_null()).select('d_call', 'd2_call').rows()
    assert pairs and set(pairs) <= FORWARD
    fit, _ = (infer_native if engine == 'native' else infer)(model, ['TGTAAACCCTTT'],
        init='template', max_iter=1, **({'threads':1} if engine == 'native' else {}))
    learned = fit.tables['d_history']
    assert set(learned.filter((pl.col('n_d') == 2) & (pl.col('p') > 0))
               .select('d_allele', 'd2_allele').iter_rows()) <= FORWARD
    assert learned['p'].sum() == pytest.approx(1)


def test_invalid_joint_history_rejected_and_legacy_checker_sees_reverse_pair():
    from vdjtools.model.check import check_model
    joint = to_history(trd_template())
    h = joint.tables['d_history'].with_columns(
        d_allele=pl.when(pl.col('n_d') == 2).then(pl.col('d2_allele')).otherwise(pl.col('d_allele')),
        d2_allele=pl.when(pl.col('n_d') == 2).then(pl.col('d_allele')).otherwise(pl.col('d2_allele')))
    bad = Model(manifest=joint.manifest, tables={**joint.tables, 'd_history':h}, genomic=joint.genomic)
    with pytest.raises(ValueError, match='physically forbidden'):
        bad.validate()
    # Construct the same formerly accepted support using a nonhuman scaffold, then declare human.
    base = trd_template()
    custom = Model(manifest=replace(base.manifest, organism='custom'), tables=base.tables, genomic=base.genomic)
    old = to_dd(custom)
    old = Model(manifest=replace(old.manifest, organism='human'), tables=old.tables, genomic=old.genomic)
    issues = check_model(old, germline='none').filter(pl.col('check') == 'tandem_support')
    assert issues.filter(pl.col('severity') == 'error').height == 3
    assert native.pgen_nt(old, 'TGTCCCAAATTT') == pgen.pgen_nt(pgen.prepare(old), 'TGTCCCAAATTT') == 0
