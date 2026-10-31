"""Explicit rank-one tandem approximation; ordinary Pgen remains exact."""
import numpy as np
import pytest
from vdjtools.model import load_bundled, native
from vdjtools.model.generate import generate


def test_rank_one_batch_is_explicit_immutable_and_thread_invariant():
    model = load_bundled('IGH')
    seqs = generate(model, 32, seed=226, productive_only=True, engine='native', threads=1)['junction_aa'].to_list() + ['']
    packed, _, _ = native.pack(model)
    before = list(packed.q_tandem)
    exact = native.pgen_aa_batch(model, seqs, threads=1)
    actual = native.pgen_aa_batch(model, seqs, threads=1, approximate_tandem=True)
    assert actual == native.pgen_aa_batch(model, seqs, threads=2, approximate_tandem=True)
    assert packed.q_tandem == before
    assert native.pgen_aa_batch(model, seqs, threads=1) == exact
    assert actual[-1] == exact[-1] == 0
    assert np.all(np.asarray(actual[:-1]) > 0)
    logp = np.log(exact[:-1]); loga = np.log(actual[:-1])
    assert np.sqrt(np.sum((loga - logp)**2) / np.sum(logp**2)) <= .05


@pytest.mark.parametrize('locus,source', [('TRA', 'isalgo'), ('TRB', 'isalgo'), ('IGH', 'olga-legacy')])
def test_rank_one_rejects_unsupported_models(locus, source):
    with pytest.raises(ValueError, match='rank-one'):
        native.pgen_aa_batch(load_bundled(locus, source), ['CF'], approximate_tandem=True)


def test_rank_one_rejects_hamming_ball():
    with pytest.raises(ValueError, match='mismatches'):
        native.pgen_aa_batch(load_bundled('IGH'), ['CF'], mismatches=1, approximate_tandem=True)


def test_rank_one_cli_reports_unsupported_model(tmp_path):
    from typer.testing import CliRunner
    from vdjtools.cli import app
    path = tmp_path / 'input.tsv'
    path.write_text('junction_aa\nCF\n')
    result = CliRunner().invoke(app, ['pgen', str(path), '--model', 'TRA', '--type', 'aa', '--approximate-tandem'])
    assert result.exit_code != 0
    assert 'rank-one' in result.output


def test_rank_one_frame_and_cli_label_approximate_output(tmp_path):
    import polars as pl
    from typer.testing import CliRunner
    from vdjtools.cli import app
    from vdjtools.model.score import pgen_frame
    model = load_bundled('IGH')
    seqs = generate(model, 1, seed=226, productive_only=True, engine='native', threads=1)['junction_aa'].to_list()
    expected = native.pgen_aa_batch(model, seqs, approximate_tandem=True, threads=1)
    frame = pgen_frame(model, seqs, kind='aa', use_calls=False, approximate_tandem=True, threads=1)
    assert frame['pgen_method'].to_list() == ['rank1']
    assert frame['pgen'].to_list() == expected
    assert 'pgen_method' not in pgen_frame(model, seqs, kind='aa', use_calls=False).columns
    source = tmp_path / 'input.tsv'; output = tmp_path / 'output.tsv'
    pl.DataFrame({'junction_aa': seqs}).write_csv(source, separator='\t')
    result = CliRunner().invoke(app, ['pgen', str(source), '--model', 'IGH', '--type', 'aa',
                                     '--approximate-tandem', '--out', str(output)])
    assert result.exit_code == 0, result.stdout
    actual = pl.read_csv(output, separator='\t')
    assert actual['pgen_method'].to_list() == ['rank1']
    np.testing.assert_allclose(actual['pgen'], expected, rtol=1e-12, atol=0)
    with pytest.raises(ValueError, match='amino-acid'):
        pgen_frame(model, ['TGT'], kind='nt', approximate_tandem=True)


def test_rank_one_memory_admission_precedes_native_scoring(monkeypatch):
    model = load_bundled('IGH')
    monkeypatch.setattr('vdjtools.cores.available_memory', lambda: 0)
    with pytest.raises(MemoryError, match='rank-one tandem'):
        native.pgen_aa_batch(model, ['CF'], approximate_tandem=True)


def test_cli_fresh_exact_scoring_replaces_old_method_and_does_not_deduplicate(tmp_path, monkeypatch):
    import polars as pl
    from typer.testing import CliRunner
    from vdjtools.cli import app
    calls = []
    def fresh_batch(model, seqs, *args, **kwargs):
        calls.append(list(seqs))
        return [.25] * len(seqs)
    monkeypatch.setattr(native, 'pgen_aa_batch', fresh_batch)
    source = tmp_path / 'old.tsv'; output = tmp_path / 'fresh.tsv'
    pl.DataFrame({'junction_aa': ['CF', 'CF'], 'pgen': [1., 1.],
                  'pgen_method': ['rank1', 'rank1']}).write_csv(source, separator='\t')
    result = CliRunner().invoke(app, ['pgen', str(source), '--model', 'IGH', '--type', 'aa', '--out', str(output)])
    assert result.exit_code == 0, result.output
    frame = pl.read_csv(output, separator='\t')
    assert calls == [['CF', 'CF']]
    assert frame['pgen'].to_list() == [.25, .25]
    assert frame['pgen_method'].to_list() == ['exact', 'exact']
