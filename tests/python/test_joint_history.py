"""Complete TRB histories normalize before inference, including D1-D2-J2."""
import polars as pl
import pytest
from collections import defaultdict

from vdjtools.model import Model
from vdjtools.model.dd import to_history
from vdjtools.model.io import from_germline


def history_model(p_nd2=0.2, *, joint=True, duplicate_alleles=False):
    entries = [('V', 'TRBV1*01', 'TGT'), ('D', 'TRBD1*01', 'AA'), ('D', 'TRBD2*01', 'CC'),
               ('J', 'TRBJ1-1*01', 'TTT'), ('J', 'TRBJ2-1*01', 'TTC')]
    if duplicate_alleles:
        entries += [('V', 'TRBV1*02', 'TGT'), ('D', 'TRBD1*02', 'AA'), ('J', 'TRBJ2-1*02', 'TTC')]
    gl = pl.DataFrame([dict(segment=seg, allele=allele, gene=allele.split('*')[0], sequence=seq,
                            functional=True, cdr3_anchor=-1, full_germline='')
                       for seg, allele, seq in entries])
    base = from_germline(gl, locus='TRB', ins_max=0,
                          palindrome_max={'v_3': 0, 'j_5': 0, 'd_5': 0, 'd_3': 0})
    tables = dict(base.tables)
    for name in ('v_3_del', 'j_5_del'):
        tables[name] = tables[name].with_columns(p=(pl.col('ndel') == 0).cast(pl.Float64))
    base = Model(manifest=base.manifest, tables=tables, genomic=base.genomic)
    return to_history(base, p_nd2=p_nd2) if joint else base


def test_joint_history_prior_normalizes_without_discarded_j1_tandems(tmp_path):
    model = history_model()
    history = model.tables['d_history']
    assert history.group_by('j_allele').agg(pl.col('p').sum())['p'].to_list() == pytest.approx([1, 1])
    tandem = history.filter(pl.col('n_d') == 2)
    assert tandem.select('d_allele', 'd2_allele', 'j_allele').rows() == [
        ('TRBD1*01', 'TRBD2*01', 'TRBJ2-1*01')]
    assert tandem['p'][0] == pytest.approx(0.2)
    assert history.filter(pl.col('j_allele') == 'TRBJ1-1*01')['p'].sum() == 1
    model.save(tmp_path / 'model')
    restored = Model.load(tmp_path / 'model')
    assert restored.tables['d_history'].equals(history)


def test_tandem_deletions_cannot_silently_lose_fully_trimmed_mass():
    model = history_model()
    for name in ('d_tandem_del', 'd2_del'):
        table = model.tables[name]
        assert table.filter(pl.col('ndel5') + pl.col('ndel3') >= 2).height == 0
        allele = 'd2_allele' if name == 'd2_del' else 'd_allele'
        assert table.group_by(allele).agg(pl.col('p').sum())['p'].to_list() == pytest.approx([1, 1])
    assert model.tables['d_del'].filter(pl.col('ndel5') + pl.col('ndel3') == 2)['p'].sum() > 0


def exact_outcomes(model):
    ds = dict(model.genomic['genes_d'].select('d_allele', 'cut_segment').iter_rows())
    js = dict(model.genomic['genes_j'].select('j_allele', 'cut_segment').iter_rows())
    pj = dict(model.tables['j_choice'].iter_rows())
    result = defaultdict(float)
    for h in model.tables['d_history'].iter_rows(named=True):
        first = model.tables['d_tandem_del' if h['n_d'] == 2 else 'd_del']
        for d1 in first.filter(pl.col('d_allele') == h['d_allele']).iter_rows(named=True):
            seq1 = ds[h['d_allele']][d1['ndel5']:len(ds[h['d_allele']])-d1['ndel3']]
            if h['n_d'] == 1:
                result['TGT' + seq1 + js[h['j_allele']]] += pj[h['j_allele']] * h['p'] * d1['p']
            else:
                for d2 in model.tables['d2_del'].filter(pl.col('d2_allele') == h['d2_allele']).iter_rows(named=True):
                    seq2 = ds[h['d2_allele']][d2['ndel5']:len(ds[h['d2_allele']])-d2['ndel3']]
                    result['TGT' + seq1 + seq2 + js[h['j_allele']]] += pj[h['j_allele']] * h['p'] * d1['p'] * d2['p']
    return result


def test_nt_and_aa_pgen_equal_complete_history_enumeration():
    from vdjtools.model import native, pgen
    from vdjtools.model.reference import translate
    model = history_model()
    expected = exact_outcomes(model)
    assert sum(expected.values()) == pytest.approx(1)
    prep = pgen.prepare(model)
    amino = defaultdict(float)
    for sequence, probability in expected.items():
        assert native.pgen_nt(model, sequence) == pytest.approx(probability)
        assert pgen.pgen_nt(prep, sequence) == pytest.approx(probability)
        if len(sequence) % 3 == 0:
            amino[translate(sequence)] += probability
    for sequence, probability in amino.items():
        assert native.pgen_aa(model, sequence) == pytest.approx(probability)
        assert pgen.pgen_aa(prep, sequence) == pytest.approx(probability)


@pytest.mark.parametrize('engine', ['reference', 'native'])
def test_joint_history_sampling_uses_the_same_supported_law(engine):
    from vdjtools.model.generate import generate
    model = history_model()
    expected = exact_outcomes(model)
    draws = generate(model, 256, seed=9, engine=engine, threads=1)
    assert all(expected[sequence] > 0 for sequence in draws['junction_nt'])
    tandem = draws.filter(pl.col('d2_call').is_not_null())
    assert tandem.height > 0
    assert set(tandem.select('d_call', 'd2_call', 'j_call').iter_rows()) == {
        ('TRBD1*01', 'TRBD2*01', 'TRBJ2-1*01')}


def test_joint_history_e_step_and_branch_deletion_counts_conserve_mass():
    from vdjtools import _core
    from vdjtools.model import native
    model = history_model()
    packed, _, _ = native.pack(model)
    for sequence, probability in exact_outcomes(model).items():
        counts = _core.make_counts(packed)
        ll = _core.estep_batch(packed, [native._encode(sequence)], [], [], [], counts)
        assert __import__('math').exp(ll) == pytest.approx(probability)
        single, tandem = sum(counts.history_single), sum(counts.history_tandem)
        assert single + tandem == pytest.approx(1)
        assert sum(counts.d_del) == pytest.approx(single)
        assert sum(counts.d_tandem_del) == pytest.approx(tandem)
        assert sum(counts.d2_del) == pytest.approx(tandem)


def test_reference_and_native_joint_history_em_agree():
    from vdjtools.model.infer import infer, infer_native
    model = history_model()
    sequences = list(exact_outcomes(model))
    reference, rr = infer(model, sequences, init='template', max_iter=2)
    fast, nr = infer_native(model, sequences, init='template', max_iter=2, threads=1)
    assert nr.loglik == pytest.approx(rr.loglik)
    assert nr.n_scoreable == rr.n_scoreable == [len(sequences)] * 2
    for name, frame in reference.tables.items():
        keys = [c for c in frame.columns if c != 'p']
        assert fast.tables[name].sort(keys)['p'].to_list() == pytest.approx(frame.sort(keys)['p'].to_list())
    assert sum(exact_outcomes(fast).values()) == pytest.approx(1)


def test_history_audit_rejects_bad_paths_before_scoring():
    from vdjtools.model.check import check_model
    from vdjtools.model.infer import enforce_dj_order
    from vdjtools.model import native
    model = history_model()
    check_model(model, germline='none', raise_on='error')
    tables = dict(model.tables)
    history = tables['d_history'].with_columns(
        d2_allele=pl.when(pl.col('n_d') == 2).then(pl.lit('TRBD1*01')).otherwise(pl.col('d2_allele')))
    invalid = Model(manifest=model.manifest, tables={**tables, 'd_history': history}, genomic=model.genomic)
    with pytest.raises(ValueError, match='forbidden'):
        invalid.validate()
    with pytest.raises(ValueError, match='forbidden'):
        enforce_dj_order(invalid)
    with pytest.raises(ValueError, match='forbidden'):
        native.pgen_nt(invalid, 'TGTAATTC')


def test_history_diagnostics_weight_active_branches():
    from vdjtools.model.analyze import bayes_net_dot, compare_models, entropy_table, gene_marginal, total_entropy
    from vdjtools.model.dd import d_count_probabilities
    model = history_model()
    assert d_count_probabilities(model) == pytest.approx({1: 0.9, 2: 0.1})
    assert gene_marginal(model, 'd') == pytest.approx({'TRBD1*01': 0.8, 'TRBD2*01': 0.2})
    assert gene_marginal(model, 'd2') == pytest.approx({'TRBD2*01': 1})
    conditional = dict(entropy_table(model).select('event', 'H_cond_bits').iter_rows())
    contribution = dict(total_entropy(model).select('event', 'contribution_bits').iter_rows())
    for name, probability in [('d_del', 0.9), ('d_tandem_del', 0.1), ('d2_del', 0.1)]:
        assert contribution[name] == pytest.approx(probability * conditional[name], abs=1e-6)
    assert 'd_history' in bayes_net_dot(model)
    assert compare_models(model, model)['tv'].to_list() == pytest.approx([0] * len(model.tables))


def test_history_conversion_preserves_distinct_second_d_deletion_seed():
    from vdjtools.model.dd import to_dd
    legacy = to_dd(history_model(joint=False))
    tables = dict(legacy.tables)
    tables['d2_del'] = tables['d2_del'].with_columns(
        p=((pl.col('ndel5') == 0) & (pl.col('ndel3') == 0)).cast(pl.Float64))
    model = to_history(Model(manifest=legacy.manifest, tables=tables, genomic=legacy.genomic))
    assert model.tables['d2_del'].filter((pl.col('ndel5') != 0) | (pl.col('ndel3') != 0))['p'].sum() == 0


def test_joint_history_collapse_marginalizes_identical_allele_templates():
    from vdjtools.model.collapse import collapse_alleles
    from vdjtools.model import native
    model = history_model(duplicate_alleles=True)
    collapsed = collapse_alleles(model).validate()
    assert collapsed.tables['d_history'].group_by('j_allele').agg(pl.col('p').sum())['p'].to_list() == pytest.approx([1, 1])
    assert collapsed.tables['d_history'].filter(pl.col('n_d') == 1)['d2_allele'].unique().to_list() == ['']
    for sequence, probability in exact_outcomes(model).items():
        assert native.pgen_nt(model, sequence) == pytest.approx(probability)
        assert native.pgen_nt(collapsed, sequence) == pytest.approx(probability)


def tandem_only_model():
    model = history_model()
    tables = dict(model.tables)
    tables['j_choice'] = tables['j_choice'].with_columns(p=(pl.col('j_allele') == 'TRBJ2-1*01').cast(pl.Float64))
    tables['d_history'] = tables['d_history'].with_columns(p=pl.when(pl.col('j_allele') == 'TRBJ2-1*01')
        .then((pl.col('n_d') == 2).cast(pl.Float64)).otherwise(pl.col('p')))
    for name, allele in [('d_tandem_del', 'd_allele'), ('d2_del', 'd2_allele')]:
        tables[name] = tables[name].with_columns(p=1.0 / pl.len().over(allele))
    return Model(manifest=model.manifest, tables=tables, genomic=model.genomic).validate()


def test_tandem_only_scenario_search_reconstructs_both_d_spans():
    from vdjtools.model.viterbi import best_scenario, infer_nt, infer_nt_batch
    from vdjtools.model import native
    model = tandem_only_model()
    for result in [best_scenario(model, 'TGTAACCTTC'), infer_nt(model, 'CNF')]:
        assert result is not None
        assert result.d_call == 'TRBD1*01'
        assert result.d2_call == 'TRBD2*01'
        assert result.cdr3_nt[result.d_start:result.d_end] in ('A', 'AA')
        assert result.cdr3_nt[result.d2_start:result.d2_end] in ('C', 'CC')
        assert 0 < result.scenario_p <= native.pgen_nt(model, result.cdr3_nt)
    batch = infer_nt_batch(model, ['CNF', 'CTF'], threads=1)
    assert batch['d2_call'].to_list() == ['TRBD2*01', 'TRBD2*01']


@pytest.mark.parametrize('insert_dd', [False, True])
def test_tandem_codon_max_matches_complete_nt_oracle(insert_dd):
    from itertools import product
    from vdjtools.model import native, pgen
    from vdjtools.model.viterbi import best_scenario, codon_options, infer_nt, infer_nt_bruteforce
    model = tandem_only_model()
    aa = 'CNF'
    if insert_dd:
        tables = dict(model.tables)
        tables['dd_ins'] = pl.DataFrame({'length': [1], 'p': [1.0]}, schema=tables['dd_ins'].schema)
        tables['dd_first_nt'] = tables['dd_first_nt'].with_columns(p=(pl.col('nt') == 2).cast(pl.Float64))
        model = Model(manifest=model.manifest, tables=tables, genomic=model.genomic).validate()
        aa = 'CSF'
    paths = [best_scenario(model, ''.join(c)) for c in product(*codon_options(aa))]
    expected = max((path.scenario_p for path in paths if path), default=0)
    assert expected == pytest.approx(1 / 9)  # two nonempty trim profiles, three equiprobable states each
    if insert_dd:
        assert native.pgen_nt(model, 'TGTAGCTTC') == pytest.approx(4 / 9)
    scenarios = native.best_aa_scenarios(model, aa)
    assert scenarios[0][0] == pytest.approx(expected)
    assert scenarios[0][9] == 'TRBD2*01'
    brute = infer_nt_bruteforce(model, aa)
    assert brute.d2_call == 'TRBD2*01'
    for result in (infer_nt(model, aa), infer_nt(pgen.prepare(model), aa, top_vj=0)):
        assert result.pgen == pytest.approx(brute.pgen)
        assert result.d2_call == 'TRBD2*01'


def test_empty_candidate_lists_have_the_same_reference_and_native_fallback():
    from vdjtools.model.infer import infer, infer_native
    model = history_model()
    sequences = list(exact_outcomes(model))
    masks = [([], [], None)] * len(sequences)
    _, reference = infer(model, sequences, masks=masks, init='template', max_iter=1)
    _, fast = infer_native(model, sequences, masks=masks, init='template', max_iter=1, threads=1)
    assert reference.loglik == pytest.approx(fast.loglik)
    assert reference.n_scoreable == fast.n_scoreable == [len(sequences)]


def test_history_conversion_and_single_d_selection_preserve_training():
    from vdjtools.model.infer import _maybe_promote_dd
    model = history_model(joint=False)
    model.training = {'runs': [{'run_id': 'prior-fit'}]}
    joint = to_history(model)
    assert joint.training == model.training
    assert _maybe_promote_dd(joint, True, 0.2).training == model.training
