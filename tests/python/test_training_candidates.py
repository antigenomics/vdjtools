"""Training uses batched proposals, not preliminary assignments or missing-D2 vetoes."""
import importlib

import polars as pl
import pytest


def test_candidate_masks_keep_alternatives_and_leave_ds_unrestricted(monkeypatch, toy_model):
    import arda
    infer = importlib.import_module("vdjtools.model.infer")
    calls = []

    def candidates(records, **kw):
        records = list(records)
        calls.append((records, kw))
        rows = [dict(sequence_id=records[0][0], segment="V", alleles=["TOYV1*99", "TOYV2*01"]),
                dict(sequence_id=records[0][0], segment="J", alleles=["TOYJ1*99"]),
                dict(sequence_id=records[0][0], segment="D", alleles=["D1*01"])]
        return {"alignments": pl.DataFrame(rows), "limits": {"max_seqs": 300}}

    monkeypatch.setattr(arda, "alignment_candidates", candidates, raising=False)
    masks, evidence = infer.candidate_masks(["ACGT"] * 3, toy_model, threads=2, batch_size=2)
    assert [len(records) for records, _ in calls] == [2, 1]
    assert all(kw["threads"] == 2 for _, kw in calls)
    assert set(masks[0][0]) == {"TOYV1*01", "TOYV1*02", "TOYV2*01"}
    assert masks[0][1] == ["TOYJ1*01"]
    assert masks[1][:2] == ([], [])  # no proposal means unrestricted, not excluded
    assert all(mask[2] is None and len(mask) == 3 for mask in masks)
    assert evidence["unrestricted_v"] == evidence["unrestricted_j"] == 1
    assert evidence["d_unrestricted"] is True


def test_builder_does_not_treat_greedy_calls_as_truth(monkeypatch):
    from test_dd import _tiny_dd_model
    from vdjtools.model import data

    infer = importlib.import_module("vdjtools.model.infer")
    model = _tiny_dd_model()
    frame = pl.DataFrame({"junction": ["GGATATCC"], "v_call": ["WRONG"],
                          "j_call": [None], "d_call": ["WRONG"], "d2_call": [None]})
    evidence = {"d_unrestricted": True}
    monkeypatch.setattr(infer, "candidate_masks", lambda *a, **kw: ([(["V"], ["J"], None)], evidence))

    def fit(base, sequences, **kw):
        assert sequences == ["GGATATCC"]
        assert kw["masks"] == [(["V"], ["J"], None)]
        assert "dd_allowed" not in kw
        assert kw["error_rate"] == 0.01
        return base, infer.InferenceReport(n_scoreable=[1], n_sequences=1,
                                            candidate_search=kw["candidate_search"])

    monkeypatch.setattr(infer, "infer_native", fit)
    _, report, stats = data.build_model("TEST", template=model, clones=frame, threads=1)
    assert stats["n_used"] == stats["n_scoreable"] == 1
    assert report.candidate_search == evidence
    _, frame_report = infer.infer_frame(model, frame, error_rate=0.01, threads=1)
    assert frame_report.candidate_search == evidence


def test_resume_cannot_silently_change_candidate_support():
    from test_insertion_ambiguity import _ambiguous_model, V, J
    infer = importlib.import_module("vdjtools.model.infer")
    model, sequences = _ambiguous_model(), [V + J]
    masks = [(["TOYV1*01"], ["TOYJ1*01"], None)]
    evidence = {"input_masks_sha256": infer._candidate_identity(sequences, masks)}
    fitted, _ = infer.infer_native(model, sequences, masks=masks, candidate_search=evidence,
                                   init="template", max_iter=1)
    with pytest.raises(ValueError, match="same sequences and masks"):
        infer.resume(fitted, sequences, max_iter=1)
    with pytest.raises(ValueError, match="same sequences and masks"):
        infer.resume(fitted, ["A" + sequences[0][1:]], masks=masks, max_iter=1)
    with pytest.raises(ValueError, match="same sequences and masks"):
        infer.resume(fitted, ["A" + sequences[0][1:]], masks=masks,
                     candidate_search=evidence, max_iter=1)
    _, report = infer.resume(fitted, sequences, masks=masks, max_iter=1)
    assert report.candidate_search == evidence
