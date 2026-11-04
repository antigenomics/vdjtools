"""ARDA independent D controls keep the legacy default command compatible."""
from types import SimpleNamespace
import pytest
from vdjtools.model import data

@pytest.mark.parametrize("map_d,map_dd,expected", [(True, True, None), (True, False, "--no-map-dd"), (False, True, "--no-map-d"), (False, False, "--no-map-d")])
def test_read_adapter_d_controls(tmp_path, monkeypatch, map_d, map_dd, expected):
    calls = []
    def run(cmd, **kwargs):
        calls.append(cmd)
        (tmp_path / "out.airr.tsv").write_text("sequence_id\tv_call\nq1\tTRBV1*01\n")
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(data.subprocess, "run", run)
    data.annotate_reads("raw.fq", out_dir=tmp_path, prefix="out", map_d=map_d, map_dd=map_dd)
    flags = [x for x in calls[0] if x in ("--no-map-d", "--no-map-dd")]
    assert flags == ([] if expected is None else [expected])


def test_prepare_forwards_d_controls(tmp_path, monkeypatch):
    captured = {}
    monkeypatch.setattr(data, "fetch_fastq", lambda *a: "raw.fq")
    monkeypatch.setattr(data, "annotate_reads", lambda *a, **kw: captured.update(kw) or "rows")
    monkeypatch.setattr(data, "unique_clonotypes", lambda rows, **kw: rows)
    assert data.prepare("human", "TRD", "nonfunctional", out_dir=tmp_path,
                        map_d=True, map_dd=False) == "rows"
    assert captured["map_d"] is True and captured["map_dd"] is False
