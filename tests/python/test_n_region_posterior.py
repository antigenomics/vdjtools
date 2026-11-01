"""Exact pure-N posterior through the owning native batch API; no inferred-NT truth."""
import json
from pathlib import Path
import numpy as np
import polars as pl
import pytest
from vdjtools.model import Model,load_bundled,n_region_posterior_batch
from vdjtools.model.io import from_germline


@pytest.fixture(scope='module')
def toy():
    germline=pl.DataFrame([dict(allele=a,gene=a.split('*')[0],segment=k,
        sequence=s,cdr3_anchor=-1,functional=True,full_germline='')
        for a,k,s in [('TOYV1*01','V','TGT'),('TOYJ1*01','J','TTT')]])
    model=from_germline(germline,locus='TOY',ins_max=3,
                        palindrome_max={'v_3':1,'j_5':0})
    tables=dict(model.tables)
    tables['v_3_del']=tables['v_3_del'].with_columns(p=pl.when(pl.col('ndel').is_in([-1,0])).then(.5).otherwise(0.))
    tables['j_5_del']=tables['j_5_del'].with_columns(p=(pl.col('ndel')==0).cast(pl.Float64))
    tables['vj_ins']=tables['vj_ins'].with_columns(p=pl.when(pl.col('length').is_in([2,3])).then(.5).otherwise(0.))
    for name in ['vj_dinucl','vj_first_nt']:
        tables[name]=tables[name].with_columns(p=pl.lit(.25))
    return Model(manifest=model.manifest,tables=tables,genomic=model.genomic).validate()


def test_fractional_origin_and_input_order(toy):
    result=n_region_posterior_batch(toy,['CKF','AAA'],positions=[[1,0,1,2],[1]],threads=1)
    assert result['row'].to_list()==[0,0,0,0,1]
    assert result['position'].to_list()==[1,0,1,2,1]
    assert np.allclose(result['n_probability'][:4],[.2,0,.2,0])
    assert result['n_probability'][-1] is None and not result['supported'][-1]
    all_positions=n_region_posterior_batch(toy,['CKF'])
    assert np.allclose(all_positions['n_probability'],[0,.2,0])
    assert n_region_posterior_batch(toy,[]).is_empty()


def test_vj_threaded_path_and_validation(toy):
    seqs=['CKF']*64;positions=[[0,1,2]]*64
    assert n_region_posterior_batch(toy,seqs,positions,threads=1).equals(
        n_region_posterior_batch(toy,seqs,positions,threads=4))
    with pytest.raises(ValueError):
        n_region_posterior_batch(toy,['CKF'],[[-1]])
    with pytest.raises(ValueError):
        n_region_posterior_batch(toy,['CKF'],[[3]])
    with pytest.raises(ValueError):
        n_region_posterior_batch(toy,['CKF'],[])
    with pytest.raises(ValueError):
        n_region_posterior_batch(toy,['CKF'],v=[])
    with pytest.raises(KeyError):
        n_region_posterior_batch(toy,['CKF'],v=['missing'])


@pytest.mark.parametrize('locus',['TRB','TRD','IGH'])
def test_current_tandem_models_thread_equality_and_immutable_anchors(locus):
    golden=json.loads((Path(__file__).parent/'fixtures/tandem_pgen_golden.json').read_text())[locus]
    model=load_bundled(locus,source='isalgo')
    seqs=golden['data']['junction_aa'][:4]
    positions=[[0,len(s)//2,len(s)-1] for s in seqs]
    vs=golden['data']['v_call'][:4];js=golden['data']['j_call'][:4]
    serial=n_region_posterior_batch(model,seqs,positions,vs,js,threads=1)
    parallel=n_region_posterior_batch(model,seqs,positions,vs,js,threads=4)
    assert serial.equals(parallel) and serial['supported'].all()
    values=serial['n_probability'].to_numpy().reshape(-1,3)
    assert np.all(values>=0) and np.all(values<=1)
    assert np.array_equal(values[:,[0,2]],np.zeros((len(seqs),2)))
