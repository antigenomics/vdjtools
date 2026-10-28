"""Distinct tandem identities and known physical exclusions on finite toy models."""
import numpy as np
import polars as pl
import pytest

from vdjtools.model import Model
from vdjtools.model.data import _build_masks
from vdjtools.model.dd import to_dd
from vdjtools.model.io import from_germline


def tandem(locus='TRD', same_gene=False):
    ds = (['TRDD1*01','TRDD1*02'] if same_gene else
          ['TRBD1*01','TRBD2*01'] if locus == 'TRB' else ['TRDD1*01','TRDD2*01'])
    js = ['TRBJ1-1*01','TRBJ2-1*01'] if locus == 'TRB' else ['TRDJ1*01']
    records = [('TRDV1*01','V','TGT'), *[(d,'D',s) for d,s in zip(ds,['AAA','CCC'])],
               *[(j,'J','TTT') for j in js]]
    germline = pl.DataFrame([dict(allele=a,gene=a.split('*')[0],segment=seg,sequence=s,
        cdr3_anchor=-1,functional=True,full_germline='') for a,seg,s in records])
    model = to_dd(from_germline(germline,locus=locus,ins_max=0,
        palindrome_max={'v_3':0,'j_5':0,'d_5':0,'d_3':0}),p_nd2=.5)
    tables = dict(model.tables)
    tables['n_d'] = tables['n_d'].with_columns(p=(pl.col('n_d') == 2).cast(pl.Float64))
    tables['d_gene'] = tables['d_gene'].with_columns(p=(pl.col('d_allele') == ds[0]).cast(pl.Float64))
    tables['d2_gene'] = tables['d2_gene'].with_columns(p=(pl.col('d2_allele') == ds[1]).cast(pl.Float64))
    for name in ('v_3_del','j_5_del'):
        tables[name] = tables[name].with_columns(p=(pl.col('ndel') == 0).cast(pl.Float64))
    for name in ('d_del','d2_del'):
        tables[name] = tables[name].with_columns(p=((pl.col('ndel5') == 0)&(pl.col('ndel3') == 0)).cast(pl.Float64))
    return Model(manifest=model.manifest,tables=tables,genomic=model.genomic).validate()


def test_prepared_masks_keep_second_d_call_and_invalid_calls_distinct():
    model = tandem()
    rows = pl.DataFrame({'v_call':['TRDV1']*3,'j_call':['TRDJ1']*3,
        'd_call':['TRDD1']*3,'d2_call':['TRDD2*01','UNKNOWN',None]})
    masks = _build_masks(rows,model)
    assert [mask[2] for mask in masks] == [['TRDD1*01']]*3
    assert [mask[3] for mask in masks] == [['TRDD2*01'],[],[]]


def test_separate_masks_match_native_reference_and_legacy_d1_only():
    from collections import defaultdict
    from vdjtools import _core
    from vdjtools.model.infer import _estep_seq
    from vdjtools.model.native import pack, _encode
    from vdjtools.model.pgen import prepare
    model = tandem()
    v,j,d1,d2 = ['TRDV1*01'],['TRDJ1*01'],['TRDD1*01'],['TRDD2*01']
    for mask in [(v,j,d1),(v,j,d1,d2)]:
        counts = defaultdict(lambda:defaultdict(float))
        assert _estep_seq(prepare(model),'TGTAAACCCTTT',counts,mask) == pytest.approx(1.)
        assert counts['d2_gene'][('TRDD1*01','TRDD2*01')] == pytest.approx(1.)
    counts = defaultdict(lambda:defaultdict(float))
    assert _estep_seq(prepare(model),'TGTAAACCCTTT',counts,(v,j,d1,[])) == 0
    packed,_,_ = pack(model)
    counts = _core.make_counts(packed)
    for d2masks in ((), ([[1]],)):
        counts = _core.make_counts(packed)
        ll = _core.estep_batch(packed,[_encode('TGTAAACCCTTT')],[[0]],[[0]],[[0]],counts,
                              1,[1],*d2masks)
        assert ll == pytest.approx(0.)
        assert np.asarray(counts.d2_gene).tolist() == pytest.approx([0,1,0,0])
    from vdjtools.cores import available_cores
    counts = _core.make_counts(packed)
    _core.estep_batch(packed,[_encode('TGTAAACCCTTT')]*128,[[0]]*128,[[0]]*128,[[0]]*128,
                     counts,min(2,available_cores()),[1]*128,[[1]]*128)
    assert np.asarray(counts.d2_gene).tolist() == pytest.approx([0,128,0,0])


@pytest.mark.parametrize('engine',['reference','native'])
def test_trb_second_d_cannot_reach_j1_but_j2_remains(engine):
    from vdjtools.model.generate import generate
    from vdjtools.model import native,pgen
    model = tandem('TRB')
    for module in (native,pgen):
        prepared = pgen.prepare(model) if module is pgen else model
        assert module.pgen_nt(prepared,'TGTAAACCCTTT',j='TRBJ1-1*01') == 0
        assert module.pgen_nt(prepared,'TGTAAACCCTTT',j='TRBJ2-1*01') == pytest.approx(.5)
        assert module.pgen_aa(prepared,'CKPF',j='TRBJ1-1*01') == 0
    draws = generate(model,32,seed=9,engine=engine,productive_only=False,threads=1)
    assert draws['j_call'].to_list() == ['TRBJ2-1*01']*32


@pytest.mark.parametrize('engine',['reference','native'])
def test_same_physical_gene_across_alleles_is_not_tandem(engine):
    from vdjtools.model.generate import generate
    from vdjtools.model import native,pgen
    model = tandem(same_gene=True)
    assert native.pgen_nt(model,'TGTAAACCCTTT') == pgen.pgen_nt(pgen.prepare(model),'TGTAAACCCTTT') == 0
    with pytest.raises(RuntimeError,match='attempt budget'):
        generate(model,1,seed=9,engine=engine,productive_only=False,threads=1)


def test_checker_only_errors_on_reachable_impossible_pairs():
    from vdjtools.model.check import check_model
    model = tandem(same_gene=True)
    issues = check_model(model,germline='none').filter(pl.col('check')=='tandem_support')
    assert issues['severity'].to_list() == ['error']
    inactive = Model(manifest=model.manifest,genomic=model.genomic,tables={**model.tables,
        'n_d':model.tables['n_d'].with_columns(p=(pl.col('n_d') == 1).cast(pl.Float64))})
    assert check_model(inactive,germline='none').filter(pl.col('check')=='tandem_support').is_empty()
    inactive = Model(manifest=model.manifest,genomic=model.genomic,tables={**model.tables,
        'd_gene':model.tables['d_gene'].with_columns(p=(pl.col('d_allele') == 'TRDD1*02').cast(pl.Float64)),
        'd2_gene':model.tables['d2_gene'].with_columns(p=pl.when(pl.col('d_allele') == 'TRDD1*02')
            .then(0.).otherwise(pl.col('p')))})
    assert check_model(inactive,germline='none').filter(pl.col('check')=='tandem_support').is_empty()
    legal = tandem('TRB')
    issues = check_model(legal,germline='none').filter(pl.col('check')=='tandem_support')
    assert issues['severity'].to_list() == ['warn']
    assert issues['allele'].to_list() == ['TRBD1*01|TRBD2*01|TRBJ1-1*01']


def test_scaffold_tandem_and_uniform_initializers_keep_physical_zeros():
    from vdjtools.model.infer import _uniform_init
    from vdjtools.model.reference import forbidden_tandem_pairs
    model = tandem('TRB')
    for table in (to_dd(from_germline(pl.DataFrame([
            dict(allele=a,gene=a.split('*')[0],segment=seg,sequence=s,cdr3_anchor=-1,functional=True,full_germline='')
            for a,seg,s in [('TRDV1*01','V','TGT'),('TRBD1*01','D','AAA'),
                            ('TRBD2*01','D','CCC'),('TRBJ2-1*01','J','TTT')]]),locus='TRB',ins_max=0),
            p_nd2=.1).tables['d2_gene'],_uniform_init(model)['d2_gene']):
        bad = forbidden_tandem_pairs(table['d_allele'].unique(),'TRB')
        for d1,d2,p in table.iter_rows():
            if (d1,d2) in bad:
                assert p == 0
        assert table['p'].is_finite().all()


def test_invalid_fourth_mask_disables_tandem_in_both_infer_paths(monkeypatch):
    from collections import defaultdict
    from vdjtools import _core
    from vdjtools.model.infer import _estep_seq,infer_native
    from vdjtools.model.pgen import prepare
    model = tandem()
    mask = (['TRDV1*01'],['TRDJ1*01'],['TRDD1*01'],['UNKNOWN'])
    assert _estep_seq(prepare(model),'TGTAAACCCTTT',defaultdict(lambda:defaultdict(float)),mask) == 0
    original = _core.estep_batch
    def observed(pm,sequences,vs,js,ds,counts,threads,ddflags,d2masks):
        assert ddflags == [0] and d2masks == [[]]
        return original(pm,sequences,vs,js,ds,counts,threads,ddflags,d2masks)
    monkeypatch.setattr(_core,'estep_batch',observed)
    infer_native(model,['TGTAAACCCTTT'],masks=[mask],max_iter=1,init='template',threads=1)
