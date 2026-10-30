"""An insertion's first nucleotide is an event even when no transition follows it."""
from dataclasses import replace
from itertools import product

import polars as pl
import pytest

from test_insertion_ambiguity import V, J, _ambiguous_model
from vdjtools.model import Model
from vdjtools.model.infer import infer, infer_native
from vdjtools.model.native import pgen_nt


def _one_insertion(sequence):
    base = _ambiguous_model()
    tables = dict(base.tables)
    tables['v_3_del'] = tables['v_3_del'].with_columns(p=(pl.col('ndel') == 0).cast(pl.Float64))
    tables['vj_ins'] = tables['vj_ins'].with_columns(p=(pl.col('length') == len(sequence)).cast(pl.Float64))
    return Model(manifest=base.manifest, tables=tables, genomic=base.genomic)


@pytest.mark.parametrize('engine', [infer, infer_native])
@pytest.mark.parametrize('insert', ['T', 'TG'])
def test_fit_learns_the_first_base_without_losing_one_base_insertions(engine, insert):
    model = _one_insertion(insert)
    sequence = V + insert + J
    assert 'vj_first_nt' in model.tables
    fitted, report = engine(model, [sequence], init='template', max_iter=2)
    assert report.n_scoreable == [1, 1]
    assert fitted.tables['vj_first_nt'].sort('nt')['p'].to_list() == pytest.approx([0, 0, 0, 1])
    assert pgen_nt(fitted, sequence) == pytest.approx(1)
    if insert == 'T':
        from vdjtools.model.generate import generate
        for generator in ('reference', 'native'):
            draws = generate(fitted, 8, seed=7, engine=generator, threads=1)
            assert draws['junction_nt'].to_list() == [sequence] * 8


def test_missing_first_base_event_keeps_legacy_stationary_convention():
    model = _one_insertion('T')
    tables = {name: table for name, table in model.tables.items() if name != 'vj_first_nt'}
    tables['vj_dinucl'] = tables['vj_dinucl'].with_columns(p=(pl.col('to_nt') + 1) / 10)
    manifest = replace(model.manifest, events={name: event for name, event in model.manifest.events.items()
                                             if name != 'vj_first_nt'})
    legacy = Model(manifest=manifest, tables=tables, genomic=model.genomic)
    assert pgen_nt(legacy, V + 'T' + J) == pytest.approx(0.4)


@pytest.mark.parametrize('engine', [infer, infer_native])
def test_unobserved_transition_rows_do_not_remove_generation_mass(engine):
    model = _one_insertion('TG')
    tables = dict(model.tables)
    tables['vj_ins'] = tables['vj_ins'].with_columns(p=pl.when(pl.col('length') > 0).then(0.5).otherwise(0.0))
    model = Model(manifest=model.manifest, tables=tables, genomic=model.genomic)
    fitted, _ = engine(model, [V + 'C' + J, V + 'TG' + J], init='template', max_iter=2)
    mass = sum(pgen_nt(fitted, V + ''.join(word) + J)
               for length in (1, 2) for word in product('ACGT', repeat=length))
    assert mass == pytest.approx(1)


def test_first_base_is_used_by_aa_dp_and_model_round_trip(tmp_path):
    from vdjtools.model.native import pgen_aa
    from vdjtools.model import pgen
    base = _one_insertion('T')
    tables = dict(base.tables)
    tables['v_3_del'] = tables['v_3_del'].with_columns(p=(pl.col('ndel') == 1).cast(pl.Float64))
    model = Model(manifest=base.manifest, tables=tables, genomic=base.genomic)
    fitted, _ = infer_native(model, [V[:-1] + 'T' + J], init='template', max_iter=1)
    fitted.save(tmp_path / 'model')
    restored = Model.load(tmp_path / 'model')
    assert restored.manifest.model_version == '2.1.0'
    assert pgen_aa(restored, 'CKPF') == pytest.approx(1)
    assert pgen.pgen_aa(pgen.prepare(restored), 'CKPF') == pytest.approx(1)


def test_entropy_uses_first_base_and_nonstationary_transition_distribution():
    from vdjtools.model.analyze import entropy_table, total_entropy
    model = _one_insertion('TG')
    tables = dict(model.tables)
    tables['vj_first_nt'] = tables['vj_first_nt'].with_columns(p=(pl.col('nt') == 3).cast(pl.Float64))
    tables['vj_dinucl'] = tables['vj_dinucl'].with_columns(
        p=pl.when((pl.col('from_nt') == 3) & pl.col('to_nt').is_in([0, 1])).then(0.5).otherwise(0.0))
    model = Model(manifest=model.manifest, tables=tables, genomic=model.genomic)
    assert 'vj_first_nt' in entropy_table(model)['event']
    contributions = dict(total_entropy(model).select('event', 'contribution_bits').iter_rows())
    assert contributions['vj_first_nt'] == 0
    assert contributions['vj_dinucl'] == pytest.approx(1)
