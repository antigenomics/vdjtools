"""Exact pure-N posterior through the owning native batch API; no inferred-NT truth."""
import json
from pathlib import Path
import numpy as np
import polars as pl
import pytest
from vdjtools.model import Model,load_bundled,n_region_posterior_batch
from vdjtools.model.io import from_germline
from vdjtools.cores import available_cores


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
        n_region_posterior_batch(toy,seqs,positions,threads=min(4,available_cores())))
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
    parallel=n_region_posterior_batch(model,seqs,positions,vs,js,threads=min(4,available_cores()))
    assert serial.equals(parallel) and serial['supported'].all()
    values=serial['n_probability'].to_numpy().reshape(-1,3)
    assert np.all(values>=0) and np.all(values<=1)
    assert np.array_equal(values[:,[0,2]],np.zeros((len(seqs),2)))
    from vdjtools.model.native import pgen_aa_batch,pgen_aa_degenerate_batch
    motifs=[list(s) for s in seqs]
    intervals=[(p[1],p[1]+1) for p in positions]
    mass=pgen_aa_degenerate_batch(model,motifs,vs,js,n_regions=intervals,threads=1)
    assert mass==pgen_aa_degenerate_batch(model,motifs,vs,js,n_regions=intervals,threads=min(4,available_cores()))
    assert np.allclose(mass,np.array(pgen_aa_batch(model,seqs,vs,js,threads=1))*values[:,1],rtol=1e-12,atol=0)


def test_masked_N_mass_equals_singleton_sum_and_keeps_default(toy):
    from vdjtools.model.native import pgen_aa_batch,pgen_aa_degenerate,pgen_aa_degenerate_batch
    aa='ACDEFGHIKLMNPQRSTVWY'
    variants=['C'+a+'F' for a in aa]
    total=np.array(pgen_aa_batch(toy,variants,threads=1))
    omega=n_region_posterior_batch(toy,variants,positions=[[1]]*20,threads=1)['n_probability'].fill_null(0).to_numpy()
    expected=float(total@omega)
    motif=['C',aa,'F']
    assert pgen_aa_degenerate(toy,motif,n_region=(1,2))==pytest.approx(expected,rel=1e-12)
    assert pgen_aa_degenerate(toy,motif,n_region=(0,1))==0
    assert pgen_aa_degenerate_batch(toy,[motif]*64,n_regions=[(1,2)]*64,threads=1)==pgen_aa_degenerate_batch(
        toy,[motif]*64,n_regions=[(1,2)]*64,threads=min(4,available_cores()))
    assert pgen_aa_degenerate_batch(toy,[motif,motif],n_regions=[None,(1,2)])==[
        pgen_aa_degenerate(toy,motif),pgen_aa_degenerate(toy,motif,n_region=(1,2))]
    assert pgen_aa_degenerate_batch(toy,[])==[]
    with pytest.raises(ValueError):
        pgen_aa_degenerate_batch(toy,[motif],n_regions=[])
    with pytest.raises(ValueError):
        pgen_aa_degenerate_batch(toy,[motif],v=[])
    for bad in [(-1,2),(2,2),(1,4),(1.5,2),(True,2)]:
        with pytest.raises(ValueError):
            pgen_aa_degenerate(toy,motif,n_region=bad)


def test_joint_N_interval_is_not_product_of_site_posteriors():
    from vdjtools.model.native import pgen_aa,pgen_aa_degenerate
    germline=pl.DataFrame([dict(allele=a,gene=a.split('*')[0],segment=k,
        sequence=s,cdr3_anchor=-1,functional=True,full_germline='')
        for a,k,s in [('TOYV1*01','V','TGT'),('TOYJ1*01','J','TTT')]])
    model=from_germline(germline,locus='TOY',ins_max=6,palindrome_max={'v_3':1,'j_5':1})
    tables=dict(model.tables)
    for name in ['v_3_del','j_5_del']:
        tables[name]=tables[name].with_columns(p=pl.when(pl.col('ndel').is_in([-1,0])).then(.5).otherwise(0.))
    tables['vj_ins']=tables['vj_ins'].with_columns(p=pl.when(pl.col('length').is_in([4,6])).then(.5).otherwise(0.))
    for name in ['vj_dinucl','vj_first_nt']:
        tables[name]=tables[name].with_columns(p=pl.lit(.25))
    model=Model(manifest=model.manifest,tables=tables,genomic=model.genomic).validate()
    sequence='CKKF'
    one=n_region_posterior_batch(model,[sequence],positions=[[1,2]])['n_probability'].to_numpy()
    joint=pgen_aa_degenerate(model,list(sequence),n_region=(1,3))/pgen_aa(model,sequence)
    assert np.isclose(pgen_aa(model,sequence),9/8192,rtol=1e-12)
    assert np.isclose(joint,1/9,rtol=1e-12)
    assert 0<joint<1 and np.allclose(one,joint,rtol=1e-12)
    assert not np.isclose(joint,float(one.prod()))
