"""Posterior trimming must not turn coincidentally matching N bases into V bases."""
from itertools import product
from math import exp

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
