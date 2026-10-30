"""Joint span messages must equal an explicit mixture of legacy history E-steps."""
from math import exp, log

import numpy as np
import polars as pl
import pytest

from test_joint_history import history_model
from vdjtools import _core
from vdjtools.model import Model, native


def _model():
    base = history_model(duplicate_alleles=True)
    tables = dict(base.tables)
    for name in ("v_3_del", "j_5_del"):
        tables[name] = tables[name].with_columns(
            p=pl.when(pl.col("ndel").is_in([0, 1])).then(0.5).otherwise(0.0))
    for channel in ("vd", "dd", "dj"):
        name = channel + "_ins"
        tables[name] = pl.DataFrame({"length": [0, 1, 2], "p": [0.5, 0.3, 0.2]},
                                   schema=tables[name].schema)
        tables[channel + "_first_nt"] = tables[channel + "_first_nt"].with_columns(
            p=(pl.col("nt") + 1) / 10)
        tables[channel + "_dinucl"] = tables[channel + "_dinucl"].with_columns(
            p=(pl.col("to_nt") + 1) / 10)
    return Model(manifest=base.manifest, tables=tables, genomic=base.genomic).validate()


def _clone(packed):
    result = _core.PackedModel()
    for name in dir(packed):
        if not name.startswith("_"):
            setattr(result, name, getattr(packed, name))
    return result


def _enumerated_histories(model, sequence, masks, rate, allow_dd):
    """Condition the unchanged legacy kernel on each complete history, then mix once."""
    packed, _, ji = native.pack(model)
    di = {a: i for i, a in enumerate(model.genomic["genes_d"]["d_allele"])}
    shape = _core.make_counts(packed)
    fields = [name for name in dir(shape) if not name.startswith("_") and name != "n_scoreable"]
    raw = {name: np.zeros_like(np.asarray(getattr(shape, name)), dtype=float) for name in fields}
    total = 0.0
    nd = len(di)
    for h in model.tables["d_history"].iter_rows(named=True):
        tandem = h["n_d"] == 2
        if not h["p"] or (tandem and not allow_dd):
            continue
        one = _clone(packed)
        one.joint_history = False
        one.dd = tandem
        one.p_nd1, one.p_nd2 = float(not tandem), float(tandem)
        j, d1 = ji[h["j_allele"]], di[h["d_allele"]]
        one.pd_given_j = [float(k == j * nd + d1) for k in range(len(ji) * nd)]
        if tandem:
            d2 = di[h["d2_allele"]]
            one.pd2_given_d1 = [float(k == d1 * nd + d2) for k in range(nd * nd)]
            one.del_d = packed.del_d_tandem
        counts = _core.make_counts(one)
        ll = _core.estep_batch(one, [native._encode_observed(sequence)], *masks[:3], counts,
                               threads=1, d2masks=masks[3], error_rate=rate)
        if not counts.n_scoreable:
            continue
        mass = h["p"] * exp(ll)
        total += mass
        for name in fields:
            if name in ("history_single", "history_tandem", "d_tandem_del") or (tandem and name == "d_del"):
                continue
            values = np.asarray(getattr(counts, name))
            if values.size:
                raw[name] += mass * values
        if tandem:
            raw["history_tandem"][(j * nd + d1) * nd + d2] += mass
            raw["d_tandem_del"] += mass * np.asarray(counts.d_del)
        else:
            raw["history_single"][j * nd + d1] += mass
    assert total > 0
    return log(total), {name: values / total for name, values in raw.items()}


@pytest.mark.parametrize("sequence,rate,masks,allow_dd", [
    ("TGTAACCTTC", 0.0, ([], [], [], []), True),
    ("TGTAGCCTTC", 0.05, ([], [], [], []), True),
    ("TGTAANCCTTC", 0.03, ([], [[1]], [[0]], [[1]]), True),
    ("TGTAATTC", 0.04, ([[0]], [], [], []), False),
    ("TGTAACCTTC", 0.0, ([[0, 0]], [[1, 1]], [[0, 0]], [[1, 1]]), True),
])
def test_joint_messages_match_complete_history_mixture(sequence, rate, masks, allow_dd):
    model = _model()
    expected_ll, expected = _enumerated_histories(model, sequence, masks, rate, allow_dd)
    packed, _, _ = native.pack(model)
    for threads, repeats in ((1, 1), (2, 64)):
        counts = _core.make_counts(packed)
        expanded = [[row[0]] * repeats if row else [] for row in masks]
        actual = _core.estep_batch(packed, [native._encode_observed(sequence)] * repeats,
                                  *expanded[:3], counts, threads=threads,
                                  dd_allowed=[int(allow_dd)] * repeats,
                                  d2masks=expanded[3], error_rate=rate)
        assert actual / repeats == pytest.approx(expected_ll, abs=1e-11)
        assert counts.n_scoreable == repeats
        for name, values in expected.items():
            assert np.asarray(getattr(counts, name)) / repeats == pytest.approx(
                values, rel=1e-10, abs=1e-11), name
