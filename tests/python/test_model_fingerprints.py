"""New fitted-model identities cover probabilities without changing legacy OLGA IDs."""
import numpy as np
import polars as pl

from vdjtools.signature import corpus as C


def test_learned_identity_detects_parameter_changes_and_ignores_table_order(monkeypatch):
    import vdjtools.model as models
    from vdjtools.model.collapse import collapse_alleles

    raw = models.load_bundled("TRG", "olga", collapse=False)
    calls = []

    def load(locus, source, *, collapse=True):
        calls.append((locus, source, collapse))
        return collapse_alleles(raw) if collapse else raw

    monkeypatch.setattr(models, "load_bundled", load)
    before = C.model_fingerprint(["TRG"], "learned")
    table = raw.tables["v_choice"]
    probabilities = table["p"].to_numpy().copy()
    i, j = np.flatnonzero(probabilities > 0)[:2]
    delta = min(probabilities[i], probabilities[j]) / 2
    probabilities[i] += delta
    probabilities[j] -= delta
    raw.tables["v_choice"] = table.with_columns(pl.Series("p", probabilities))
    after = C.model_fingerprint(["TRG"], "learned")
    assert before != after, "changed probabilities must invalidate a fitted-model corpus identity"
    assert len(calls) == 2, "load each model only once per identity calculation"
    insertions = raw.tables["vj_ins"]
    probabilities = insertions["p"].to_numpy().copy()
    i, j = np.flatnonzero(probabilities > 0)[:2]
    probabilities[j] += probabilities[i]
    probabilities[i] = 0.0
    raw.tables["vj_ins"] = insertions.with_columns(pl.Series("p", probabilities))
    zero = C.model_fingerprint(["TRG"], "learned")
    probabilities[i] = 1e-200
    raw.tables["vj_ins"] = insertions.with_columns(pl.Series("p", probabilities))
    assert C.model_fingerprint(["TRG"], "learned") != zero
    raw.tables["vj_ins"] = insertions
    raw.tables = {name: frame.reverse().select(frame.columns[::-1])
                  for name, frame in reversed(list(raw.tables.items()))}
    assert C.model_fingerprint(["TRG"], "learned") == after


def test_explicit_olga_identity_retains_published_value():
    assert C.model_fingerprint(["TRG"], "olga") == {
        "TRG": "olga:human_T_gamma@d25719222561"}
