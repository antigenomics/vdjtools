import numpy as np
import pytest
from vdjtools.model.selection import ALICE_TRB_SELECTION_Q, selection_mass

def test_scaling_and_validation():
    np.testing.assert_array_equal(selection_mass([0.,.1],4),[0.,.4])
    assert selection_mass(.1)==.1
    assert ALICE_TRB_SELECTION_Q==9.41
    for mass,q in [(-1,1),(float('nan'),1),(.1,0),(.1,float('inf'))]:
        with pytest.raises(ValueError):
            selection_mass(mass,q)


def test_selection_fit_conditions_on_observation():
    from vdjtools.model.selection import fit_selection_factor
    q = fit_selection_factor([1, 2, 3], [1/6]*3, total_events=6)
    assert q == pytest.approx(1.59362426004004)
    # Naive observed/expected would give two, ignoring the zero truncation.
    assert q != pytest.approx(2.)
    with pytest.raises(ValueError, match="singletons"):
        fit_selection_factor([1, 1], [.1, .2], total_events=2)
    with pytest.raises(ValueError, match="unsupported"):
        fit_selection_factor([1, 2], [0., .2], total_events=3)


def test_nt_control_events_count_variants_and_reject_ambiguous_units(tmp_path):
    from vdjtools.model.selection import control_event_counts
    path = tmp_path / "control.tsv"
    path.write_text("cdr3aa\tincidence\toccurrences\tconvergence\nCASSF\t2\t2\t1\nCASSF\t3\t3\t1\nC*FF\t4\t4\t1\n")
    counts, audit = control_event_counts(path)
    assert counts.to_dicts() == [dict(junction_aa="CASSF", events=5)]
    assert audit["source_events"] == 9 and audit["invalid_events"] == 4
    path.write_text(path.read_text().replace("CASSF\t2\t2", "CASSF\t2\t3"))
    with pytest.raises(ValueError, match="occurrence/incidence"):
        control_event_counts(path)


def test_calibration_cli_records_model_and_event_units(tmp_path):
    import json
    from typer.testing import CliRunner
    from vdjtools.cli import app
    control = tmp_path / "control.tsv"
    control.write_text("cdr3aa\tincidence\toccurrences\tconvergence\nCASSLGQAYEQYF\t2\t2\t1\nCASSLGRAYEQYF\t3\t3\t1\n")
    out = tmp_path / "result.json"
    run = CliRunner().invoke(app, ["calibrate-selection", str(control), "--out", str(out), "--locus", "TRB", "--threads", "1", "--max-events", "4", "--seed", "17"])
    assert run.exit_code == 0, run.output
    report = json.loads(out.read_text())
    assert report["selection_q"] > 0 and report["audit"]["retained_events"] == 4
    assert report["audit"]["population_events"] == 5 and report["audit"]["sampling_seed"] == 17
    assert report["source"] == "isalgo" and report["model_sha256"]
    assert report["zero_mass_junctions"] == 0


def test_selection_factor_requires_matching_bundled_model(monkeypatch):
    import hashlib
    from vdjtools.model.selection import selection_factor
    assert selection_factor('TRA',organism='human') == pytest.approx(3.46977501091145)
    for locus,source in [('IGK','isalgo'),('TRA','olga-legacy')]:
        with pytest.raises(ValueError,match='calibrated'):
            selection_factor(locus,source=source)
    class Changed:
        def hexdigest(self):return '0'*64
    monkeypatch.setattr(hashlib,'sha256',lambda data:Changed())
    with pytest.raises(ValueError,match='changed'):
        selection_factor('TRA')


def test_control_subsample_preserves_events_and_seed(tmp_path):
    from vdjtools.model.selection import control_event_counts
    path = tmp_path / "control.tsv"
    path.write_text("cdr3aa\tincidence\nCASSF\t100\nCASRF\t200\nCASQF\t300\n")
    a, audit = control_event_counts(path, max_events=100, seed=42)
    b, _ = control_event_counts(path, max_events=100, seed=42)
    assert a.equals(b)
    assert a["events"].sum() == audit["retained_events"] == 100
    assert audit["population_events"] == 600 and audit["sampling_seed"] == 42
    assert all(x <= y for x, y in zip(a["events"], [100, 300, 200]))
    full, audit = control_event_counts(path, max_events=1000, seed=42)
    assert full["events"].sum() == audit["retained_events"] == 600
    with pytest.raises(ValueError, match="positive"):
        control_event_counts(path, max_events=0)


def test_four_context_selection_factors_are_model_specific():
    from vdjtools.model import selection_factor
    expected = {('human', 'TRA'): 3.46977501091145,
                ('human', 'TRB'): 6.5945522588448275,
                ('mouse', 'TRA'): 4.249595300676978,
                ('mouse', 'TRB'): 5.356167549507866}
    for (species, locus), q in expected.items():
        assert selection_factor(locus, organism=species) == pytest.approx(q, rel=1e-12)


def test_tcrnet_cli_alice_uses_new_default(tmp_path):
    import polars as pl
    from typer.testing import CliRunner
    from vdjtools.cli import app
    from vdjtools.model import selection_factor
    sample = tmp_path / 'sample.tsv'
    pl.DataFrame({'junction_aa': ['CASSF', 'CATSF', 'CATTF'],
                  'v_call': ['TRBV19'] * 3, 'j_call': ['TRBJ2-1'] * 3,
                  'duplicate_count': [5] * 3}).write_csv(sample, separator='\t')
    out = tmp_path / 'result.tsv'
    result = CliRunner().invoke(app, ['tcrnet', str(sample), '--mode', 'alice',
        '--locus', 'TRB', '--threads', '1', '--out', str(out)])
    assert result.exit_code == 0, result.output
    scored = pl.read_csv(out, separator='\t')
    assert scored.height > 0
    assert scored['selection_q'].to_list() == [selection_factor('TRB')] * scored.height


def test_streamed_profile_events_keep_gene_strata_and_sampling_units(tmp_path):
    from vdjtools.model.selection import control_event_counts
    path = tmp_path / 'control.tsv'
    path.write_text('cdr3aa\tv\tj\tincidence\nCASSF\tV1\tJ1\t100\nCASSF\tV2\tJ1\t200\nCATSF\tV1\tJ2\t300\n')
    first, audit = control_event_counts(path, max_events=100, seed=42, preserve_vj=True)
    second, _ = control_event_counts(path, max_events=100, seed=42, preserve_vj=True)
    assert first.equals(second)
    assert first['events'].sum() == 100
    assert audit['population_events'] == 600 and audit['retained_events'] == 100
    assert audit['observed_junctions'] == first['junction_aa'].n_unique()
    assert audit['observed_strata'] == first.height
    assert set(first.select('v_call', 'j_call').iter_rows()) == {('V1', 'J1'), ('V2', 'J1'), ('V1', 'J2')}
    full, _ = control_event_counts(path, max_events=1000, seed=42, preserve_vj=True)
    assert full['events'].to_list() == [100, 200, 300]
    missing = tmp_path / 'missing.tsv'
    missing.write_text('cdr3aa\tincidence\nCASSF\t3\n')
    retained, _ = control_event_counts(missing, preserve_vj=True)
    assert retained['events'].sum() == 3 and retained['v_call'].null_count() == 1


@pytest.mark.parametrize('preserve_vj', [False, True])
def test_control_sampler_reports_large_population_limit(tmp_path, preserve_vj):
    from vdjtools.model.selection import control_event_counts
    path=tmp_path/'large.tsv'
    path.write_text('cdr3aa\tincidence\nCASSF\t2000000000\nCATSF\t2000000000\n')
    with pytest.raises(ValueError, match='subsampling requires population below 1000000000'):
        control_event_counts(path,max_events=100,preserve_vj=preserve_vj)
    full,_=control_event_counts(path,preserve_vj=preserve_vj)
    assert full['events'].sum()==4000000000


def test_full_profile_pool_bypasses_sampling_across_large_batches(tmp_path):
    from vdjtools.model.selection import control_event_counts
    path=tmp_path/'large-full.tsv'
    path.write_text('cdr3aa\tincidence\n'+'CASSF\t40000\n'*65536+'CATSF\t40000\n')
    counts,audit=control_event_counts(path,preserve_vj=True)
    assert counts['events'].sum()==65537*40000
    assert audit['retained_events']==audit['population_events']==65537*40000
