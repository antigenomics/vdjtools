"""Unknown observations are marginalized once per read, never replaced by a nucleotide."""
from math import exp

import polars as pl
import pytest

from test_insertion_ambiguity import V, J, _ambiguous_model
from vdjtools._core import estep_batch, make_counts
from vdjtools.model import native
from vdjtools.model.infer import infer_frame, infer_native


@pytest.mark.parametrize("rate", [0.0, 0.03])
def test_unknown_v_or_insert_bases_preserve_both_histories(rate):
    model = _ambiguous_model()
    packed, _, _ = native.pack(model)
    observed = native._encode_observed(V[:-2] + "NN" + J)
    for threads in (1, 2):
        counts = make_counts(packed)
        ll = estep_batch(packed, [observed] * 80, [], [], [], counts, threads, error_rate=rate)
        # Nine V bases + zero N and seven V bases + two unknown N have equal weight.
        assert exp(ll / 80) == pytest.approx(0.5 * (1 - rate) ** 10)
        assert counts.n_scoreable == 80
        assert counts.ins_vj == pytest.approx([40, 0, 40])
        assert counts.dinucl_vj == pytest.approx([2.5] * 16)
        assert counts.first_vj == pytest.approx([10] * 4)
        assert counts.error_mismatches == pytest.approx(80 * rate)
        assert counts.error_nt == pytest.approx(80 * 11)


def test_frame_default_preserves_iupac_and_reports_observation_count():
    frame = pl.DataFrame({"junction": [V[:-2] + "NN" + J]})
    _, report = infer_frame(_ambiguous_model(), frame, use_calls=False, init="template", max_iter=1)
    assert report.n_sequences == report.n_ambiguous == 1
    assert report.n_scoreable == [1]
    assert exp(report.loglik[0]) == pytest.approx(0.5)


def test_observation_expansion_has_explicit_resource_limit():
    with pytest.raises(ValueError, match="max_observation_states"):
        infer_native(_ambiguous_model(), ["N" * 7], init="template", max_iter=1)
    with pytest.raises(ValueError, match="IUPAC"):
        native._encode_observed("ACXG")


def test_python_reference_rejects_ambiguous_observations_explicitly():
    from vdjtools.model.infer import infer
    with pytest.raises(ValueError, match="infer_native"):
        infer(_ambiguous_model(), ["N" + (V + J)[1:]], init="template", max_iter=1)


@pytest.mark.parametrize("symbol,allowed", list(zip(
    "RYSWKMBDHVN", ("AG", "CT", "CG", "AT", "GT", "AC", "CGT", "AGT", "ACT", "ACG", "ACGT"))))
def test_every_iupac_code_sums_compatible_insertion_bases(symbol, allowed):
    packed, _, _ = native.pack(_ambiguous_model())
    counts = make_counts(packed)
    ll = estep_batch(packed, [native._encode_observed(V[:-2] + symbol + "A" + J)],
                     [], [], [], counts)
    assert exp(ll) == pytest.approx(len(allowed) / 64)
    assert counts.ins_vj == pytest.approx([0, 0, 1])
    assert counts.dinucl_vj[:4] == pytest.approx([
        1 / len(allowed) if base in allowed else 0 for base in "ACGT"])


def test_ambiguous_tandem_counts_are_likelihood_weighted_not_equal_completion_votes():
    from test_dd import _tiny_dd_model
    packed, _, _ = native.pack(_tiny_dd_model(pdel_full_only=False, dd_ins1=0.5))
    completions = []
    for base in "ACGT":
        counts = make_counts(packed)
        ll = estep_batch(packed, [native._encode("GG" + base + "CTCC")], [], [], [], counts,
                         error_rate=0.09)
        completions.append((exp(ll) if counts.n_scoreable else 0, counts))
    total = sum(weight for weight, _ in completions)
    actual = make_counts(packed)
    ll = estep_batch(packed, [native._encode_observed("GGNCTCC")], [], [], [], actual, error_rate=0.09)
    assert exp(ll) == pytest.approx(total)
    for field in ("d_del", "d2_del", "d2_gene", "n_d", "ins_dd", "dinucl_dd"):
        expected = [sum(w * value for (w, _), value in zip(completions, values)) / total
                    for values in zip(*(getattr(counts, field) for _, counts in completions))]
        assert getattr(actual, field) == pytest.approx(expected)
    for field in ("error_mismatches", "error_nt"):
        assert getattr(actual, field) == pytest.approx(
            sum(w * getattr(counts, field) for w, counts in completions) / total)


@pytest.mark.parametrize("encoded", [[4], [-16], [127], [-128]])
def test_invalid_native_observation_cannot_index_markov_arrays(encoded):
    packed, _, _ = native.pack(_ambiguous_model())
    with pytest.raises(ValueError, match="IUPAC"):
        estep_batch(packed, [encoded], [], [], [], make_counts(packed))
