"""Posterior trimming must not turn coincidentally matching N bases into V bases."""
from itertools import product
from math import exp, log

import polars as pl
import pytest

from vdjtools._core import estep_batch, make_counts
from vdjtools.model import Model
from vdjtools.model.io import from_germline
from vdjtools.model.native import _encode, pack


V, J = "TGTAAACCC", "TTT"


def _ambiguous_model():
    germline = pl.DataFrame([
        dict(allele=f"TOY{segment}1*01", gene=f"TOY{segment}1", segment=segment,
             sequence=sequence, cdr3_anchor=-1, functional=True, full_germline="")
        for segment, sequence in (("V", V), ("J", J))
    ])
    model = from_germline(germline, locus="TOY", ins_max=2,
                          palindrome_max={"v_3": 0, "j_5": 0})
    tables = dict(model.tables)
    tables["v_3_del"] = tables["v_3_del"].with_columns(
        p=pl.when(pl.col("ndel").is_in([0, 2])).then(0.5).otherwise(0.0))
    tables["j_5_del"] = tables["j_5_del"].with_columns(p=(pl.col("ndel") == 0).cast(pl.Float64))
    tables["vj_ins"] = tables["vj_ins"].with_columns(
        p=pl.when(pl.col("length").is_in([0, 2])).then(0.5).otherwise(0.0))
    tables["vj_dinucl"] = tables["vj_dinucl"].with_columns(p=pl.lit(0.25))
    return Model(manifest=model.manifest, tables=tables, genomic=model.genomic).validate()


def _counts(sequences):
    packed, _, _ = pack(_ambiguous_model())
    counts = make_counts(packed)
    loglik = estep_batch(packed, [_encode(s) for s in sequences], [], [], [], counts, 1)
    return loglik, counts


def test_nine_matching_v_bases_also_allow_seven_v_plus_two_n():
    loglik, counts = _counts([V + J])
    # Two histories: .5(del0)*.5(ins0), and .5(del2)*.5(ins2)/16.
    assert exp(loglik) == pytest.approx(17 / 64)
    assert counts.ins_vj == pytest.approx([16 / 17, 0.0, 1 / 17])
    assert counts.v_3_del[0] == pytest.approx(16 / 17)
    assert counts.v_3_del[2] == pytest.approx(1 / 17)


def test_exact_repertoire_preserves_zero_insertion_mass_while_greedy_inflates_it():
    words = ["".join(bases) for bases in product("ACGT", repeat=2)]
    # Exact frequencies of every outcome, not a random draw: 64 equally weighted events.
    sequences = [V + J] * 16 + [V[:-2] + J] * 16
    sequences += [prefix + word + J for prefix in (V, V[:-2]) for word in words]
    _, counts = _counts(sequences)
    assert len(sequences) == 64
    assert counts.ins_vj == pytest.approx([32.0, 0.0, 32.0])
    greedy_zero = sum(s in (V + J, V[:-2] + J) for s in sequences)
    assert greedy_zero == 33  # posterior zero-insertion mass is 32/64, not greedy 33/64.


def test_native_reports_only_positive_likelihood_reads():
    from vdjtools.model.infer import infer, infer_native

    sequences = [V + J, "A" * (len(V) + len(J))] * 40
    expected = infer(_ambiguous_model(), sequences, init="template", max_iter=1)[1]
    for threads in (1, 2):
        _, report = infer_native(_ambiguous_model(), sequences, init="template",
                                 max_iter=1, threads=threads)
        assert report.n_sequences == 80
        assert report.n_scoreable == expected.n_scoreable == [40]
        assert report.loglik == pytest.approx(expected.loglik)

    _, report = infer_native(_ambiguous_model(), ["A" * 12], init="template", max_iter=1)
    assert report.n_scoreable == [0]
    assert report.loglik == [float("-inf")]


def test_training_substitution_likelihood_does_not_change_generation_pgen():
    from vdjtools.model.infer import infer_native, resume
    from vdjtools.model.native import pgen_nt

    template = _ambiguous_model()
    tables = dict(template.tables)
    for event, column in (("v_3_del", "ndel"), ("vj_ins", "length")):
        tables[event] = tables[event].with_columns(p=(pl.col(column) == 0).cast(pl.Float64))
    model = Model(manifest=template.manifest, tables=tables, genomic=template.genomic)
    sequence = "A" + (V + J)[1:]
    packed, _, _ = pack(model)
    counts = make_counts(packed)
    ll = estep_batch(packed, [_encode(sequence)], [], [], [], counts, 1, error_rate=0.03)
    assert exp(ll) == pytest.approx(0.01 * 0.97 ** 11)
    assert counts.n_scoreable == 1
    assert counts.error_mismatches == pytest.approx(1)
    assert counts.error_nt == pytest.approx(12)
    assert pgen_nt(model, sequence) == 0
    fitted, report = infer_native(model, [sequence], init="template", max_iter=1, error_rate=0.03)
    assert report.n_scoreable == [1]
    assert report.error_rate == pytest.approx([1 / 12])
    assert fitted.manifest.error_rate == pytest.approx(1 / 12)
    assert fitted.training["runs"][-1]["error_model"] == "retained-germline-substitution"
    _, resumed = resume(fitted, [sequence], max_iter=1)
    assert resumed.error_rate == pytest.approx([1 / 12])
    assert resumed.loglik == pytest.approx([log((1 / 36) * (11 / 12) ** 11)])
    _, default_report = infer_native(model, [sequence], max_iter=1, error_rate=0.03)
    assert default_report.n_scoreable == [1]
    assert default_report.init == "uniform"


@pytest.mark.parametrize("tandem,sequence,errors", [(False, "GGCTCC", 1), (True, "GGCTCTCC", 2)])
def test_single_and_tandem_d_error_statistics(tandem, sequence, errors):
    from test_dd import _tiny_dd_model

    packed, _, _ = pack(_tiny_dd_model(p_nd2=float(tandem)))
    for threads in (1, 2):
        counts = make_counts(packed)
        ll = estep_batch(packed, [_encode(sequence)] * 80, [], [], [], counts, threads,
                         error_rate=0.03)
        probability = 0.01 ** errors * 0.97 ** (len(sequence) - errors)
        assert ll / 80 == pytest.approx(log(probability))
        assert counts.n_scoreable == 80
        assert counts.error_mismatches == pytest.approx(80 * errors)
        assert counts.error_nt == pytest.approx(80 * len(sequence))
        assert counts.n_d[2 if tandem else 1] == pytest.approx(80)


@pytest.mark.parametrize("rate", [-0.01, 1.01, float("nan"), float("inf")])
def test_invalid_training_error_rate_is_rejected(rate):
    packed, _, _ = pack(_ambiguous_model())
    with pytest.raises(ValueError, match="error_rate"):
        estep_batch(packed, [_encode(V + J)], [], [], [], make_counts(packed), error_rate=rate)


def test_tandem_error_messages_match_independent_scenario_enumeration():
    from test_dd import _tiny_dd_model

    sequence, rate = "GGACTCC", 0.09
    scenarios = []
    pieces = [("AT", 0.5), ("A", 0.25), ("T", 0.25)]
    for (left, pleft), (right, pright), insert in product(pieces, pieces, ["", *"ACGT"]):
        emitted = "GG" + left + insert + right + "CC"
        if len(emitted) != len(sequence):
            continue
        templated = [True] * (2 + len(left)) + [False] * len(insert) + [True] * (len(right) + 2)
        if any(a != b for a, b, germline in zip(emitted, sequence, templated) if not germline):
            continue
        errors = sum(a != b for a, b, germline in zip(emitted, sequence, templated) if germline)
        bases = sum(templated)
        prior = 0.5 * pleft * pright * 0.5 * (0.25 if insert else 1)
        weight = prior * (rate / 3) ** errors * (1 - rate) ** (bases - errors)
        scenarios.append((weight, errors, bases, len(insert)))
    total = sum(w for w, _, _, _ in scenarios)
    packed, _, _ = pack(_tiny_dd_model(pdel_full_only=False, dd_ins1=0.5))
    counts = make_counts(packed)
    ll = estep_batch(packed, [_encode(sequence)], [], [], [], counts, error_rate=rate)
    assert exp(ll) == pytest.approx(total)
    assert counts.error_mismatches == pytest.approx(sum(w * e for w, e, _, _ in scenarios) / total)
    assert counts.error_nt == pytest.approx(sum(w * n for w, _, n, _ in scenarios) / total)
    assert counts.ins_dd == pytest.approx([
        sum(w for w, _, _, length in scenarios if length == k) / total for k in (0, 1)])
