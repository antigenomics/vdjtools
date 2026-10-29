"""Complete emission metadata and raw allele lookup survive major-only annotation defaults."""
from types import SimpleNamespace
import sys

from vdjtools.model import reference


def test_d_completeness_comes_from_exact_template_metadata(tmp_path, monkeypatch):
    import arda.cdr3fix
    import arda.paths

    (tmp_path / "d_germlines.fasta").write_text(">TRB|TRBD1*01\nACGT\n")
    monkeypatch.setattr(arda.paths, "vdj_dir", lambda organism: tmp_path)
    monkeypatch.setattr(arda.cdr3fix, "load_anchors", lambda organism: {})
    monkeypatch.setitem(sys.modules, "arda.alleles", SimpleNamespace(
        d_germline_metadata=lambda organism: {("TRBD1*01", "ACGT"): {"status": "partial"}}))
    frame = reference.load_germline.__wrapped__("TRB")
    assert frame.select("allele", "status").rows() == [("TRBD1*01", "partial")]


def test_full_germline_requests_all_actual_alleles(tmp_path, monkeypatch):
    import arda.annotate.reference

    fasta = tmp_path / "alleles.fasta"
    fasta.write_text(">s\nACGTNNNTTT\n")
    seen = []
    def load(organism, seqtype, *, major_only):
        seen.append((organism, seqtype, major_only))
        return SimpleNamespace(target_fasta=fasta, entries={"s": SimpleNamespace(
            v_sequence_end=4, v_call="TRBV1*01,TRBV1*02", j_sequence_start=8,
            j_call="TRBJ1*01", vj_end=10)})
    monkeypatch.setattr(arda.annotate.reference, "load_reference", load)
    full = reference.load_full_vj_germline.__wrapped__()
    assert seen == [("human", "nt", False)]
    assert full[("V", "TRBV1*02")] == full[("V", "TRBV1*01")] == "ACGT"
