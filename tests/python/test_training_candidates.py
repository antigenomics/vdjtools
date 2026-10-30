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
                          "j_call": [None], "d_call": ["WRONG"], "d2_call": [None],
                          "candidate_sequences": [["AAAAGGATATCC", "GGATATCCTTTT"]]})
    evidence = {"d_unrestricted": True}
    def candidates(*a, **kw):
        assert kw["query_sequences"] == [["AAAAGGATATCC", "GGATATCCTTTT"]]
        return [(["V"], ["J"], None)], evidence

    monkeypatch.setattr(infer, "candidate_masks", candidates)

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


def test_full_read_alternatives_survive_repeated_clonotype_collapse():
    from vdjtools.model.data import unique_clonotypes

    frame = pl.DataFrame({"v_call": ["V*01", "V*02", "V*01"],
                          "j_call": ["J*01"] * 3, "junction": ["ACGT"] * 3,
                          "sequence": ["AACCGT", "ATGC", "AACCGT"],
                          "rev_comp": ["T", "F", "T"]})
    once = unique_clonotypes(frame, gene_level=True)
    twice = unique_clonotypes(once, gene_level=True)
    assert once["candidate_sequences"].to_list() == [["ACGGTT", "ATGC"]]
    assert twice.equals(once)
    assert once["junction"].to_list() == ["ACGT"]
    assert once["count"].to_list() == [3]


def test_full_read_candidates_union_across_flat_batches(monkeypatch, toy_model):
    import arda
    infer = importlib.import_module("vdjtools.model.infer")
    batches = []

    def candidates(records, **kw):
        batches.append(records)
        rows = [dict(sequence_id=q, segment="V", alleles=["TOYV1*99" if s == "AAAA" else "TOYV2*01"])
                for q, s in records if s != "ACGT"]
        return {"alignments": pl.DataFrame(rows, schema={"sequence_id": pl.String,
                "segment": pl.String, "alleles": pl.List(pl.String)}), "limits": {"max_seqs": 300}}

    monkeypatch.setattr(arda, "alignment_candidates", candidates)
    masks, evidence = infer.candidate_masks(["ACGT", "ACGT"], toy_model, threads=1,
        batch_size=1, query_sequences=[["AAAA", "TTTT", "AAAA"], []])
    assert [record for batch in batches for record in batch] == [("0", "AAAA"), ("1", "TTTT"), ("2", "ACGT")]
    assert masks[0][0] == ["TOYV1*01", "TOYV1*02", "TOYV2*01"]
    assert masks[1] == ([], [], None)
    assert evidence["query_count"] == 3
    import hashlib
    import json
    flat_queries = [(0, "AAAA"), (0, "TTTT"), (1, "ACGT")]
    assert evidence["query_inputs_sha256"] == hashlib.sha256(
        json.dumps(flat_queries, separators=(",", ":")).encode()).hexdigest()
    assert evidence["junction_fallback_observations"] == evidence["full_read_observations"] == 1
    assert evidence["input_masks_sha256"] == infer._candidate_identity(["ACGT", "ACGT"], masks)
    with pytest.raises(ValueError, match="one list per observation"):
        infer.candidate_masks(["ACGT"], toy_model, query_sequences=[])


def test_infer_frame_uses_oriented_full_read_only_for_proposals(monkeypatch):
    from test_dd import _tiny_dd_model
    infer = importlib.import_module("vdjtools.model.infer")
    frame = pl.DataFrame({"junction": ["GGATATCC"], "sequence": ["AACCGT"], "rev_comp": [True]})

    def candidates(sequences, model, **kw):
        assert sequences == ["GGATATCC"]
        assert kw["query_sequences"] == [["ACGGTT"]]
        return [([], [], None)], {}

    def fit(model, sequences, **kw):
        assert sequences == ["GGATATCC"]
        assert kw["error_rate"] == 0.01
        return model, None

    monkeypatch.setattr(infer, "candidate_masks", candidates)
    monkeypatch.setattr(infer, "infer_native", fit)
    infer.infer_frame(_tiny_dd_model(), frame, error_rate=0.01, threads=1)


def test_clonotype_counts_keep_valid_support_and_default_missing_reads():
    from vdjtools.model.data import unique_clonotypes
    frame = pl.DataFrame({"v_call": ["V"] * 4, "j_call": ["J"] * 4,
                          "junction": ["ACGT"] * 4,
                          "duplicate_count": ["3", "", None, "0"],
                          "count": [9, 4, None, 9]})
    once = unique_clonotypes(frame)
    assert once["count"].to_list() == [8]  # 3 + fallback count 4 + one uncounted read + explicit zero
    assert unique_clonotypes(once).equals(once)


def test_arda_masks_preserves_annotation_orientation_and_checks_legacy_rows(monkeypatch):
    import arda
    from test_dd import _tiny_dd_model
    infer = importlib.import_module("vdjtools.model.infer")

    def candidates(sequences, model, **kw):
        assert sequences == ["GGATATCC"]
        assert kw["query_sequences"] == [["ACGGTT"]]
        return [([], [], None)], {}

    monkeypatch.setattr(infer, "candidate_masks", candidates)
    monkeypatch.setattr(arda, "annotate_sequences", lambda *a, **kw:
        [{"junction": "GGATATCC", "sequence": "AACCGT", "rev_comp": "T"}])
    assert infer.arda_masks(["AACCGT"], _tiny_dd_model())[0] == ["GGATATCC"]
    monkeypatch.setattr(arda, "annotate_sequences", lambda *a, **kw:
        [{"junction": "GGATATCC", "rev_comp": "T"}])
    infer.arda_masks(["AACCGT"], _tiny_dd_model())
    with pytest.raises(ValueError, match="one row per input contig"):
        infer.arda_masks(["AACCGT", "OTHER"], _tiny_dd_model())


def test_candidate_queries_reject_string_rows(toy_model):
    infer = importlib.import_module("vdjtools.model.infer")
    with pytest.raises(TypeError, match="rows must be lists"):
        infer.candidate_masks(["ACGT"], toy_model, query_sequences=["AACCGT"])


def test_clonotype_candidate_order_is_independent_of_raw_row_order():
    from vdjtools.model.data import unique_clonotypes
    frame = pl.DataFrame({'v_call': ['V2', 'V1'], 'j_call': ['J', 'J'],
                          'junction': ['ACGT', 'TGCA'], 'sequence': ['AACGT', 'ATGCA']})
    canonical = unique_clonotypes(frame, gene_level=True)
    assert canonical['v_call'].to_list() == ['V1', 'V2']
    assert unique_clonotypes(frame.reverse(), gene_level=True).equals(canonical)


def test_unmapped_arda_junction_does_not_replace_observation_with_full_read(monkeypatch):
    import arda
    from test_dd import _tiny_dd_model
    infer = importlib.import_module('vdjtools.model.infer')
    monkeypatch.setattr(arda, 'annotate_sequences', lambda *a, **kw:
                        [{'junction': None, 'sequence': 'ACGT', 'rev_comp': False}])
    def candidates(sequences, model, **kw):
        assert sequences == ['']
        assert kw['query_sequences'] == [['ACGT']]
        return [([], [], None)], {}
    monkeypatch.setattr(infer, 'candidate_masks', candidates)
    assert infer.arda_masks(['ACGT'], _tiny_dd_model()) == ([''], [([], [], None)])
