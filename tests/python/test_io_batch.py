"""Tests for vdjtools.io.batch — format auto-detection and metadata-driven batch reads."""
import gzip

import polars as pl

from vdjtools import io as vio
from vdjtools.io import schema as S

NATIVE_HEADER = "count\tfreq\tcdr3nt\tcdr3aa\tv\td\tj\tVEnd\tDStart\tDEnd\tJStart\n"


def _write_native(path, rows):
    """Write a tiny gzipped native table; rows are (count, cdr3nt, cdr3aa, v, d, j)."""
    with gzip.open(path, "wt") as f:
        f.write(NATIVE_HEADER)
        for count, nt, aa, v, d, j in rows:
            f.write(f"{count}\t0\t{nt}\t{aa}\t{v}\t{d}\t{j}\t3\t-1\t-1\t5\n")


def test_sniff_and_read_dispatch(tmp_path):
    native = tmp_path / "n.tsv.gz"
    _write_native(native, [(1, "TGT", "CASSL", "TRBV1", ".", "TRBJ1")])
    airr = tmp_path / "a.tsv"
    airr.write_text("v_call\tj_call\tjunction_aa\tduplicate_count\nTRBV1\tTRBJ1\tCASSLF\t3\n")
    assert vio.sniff_format(native) == "vdjtools"
    assert vio.sniff_format(airr) == "airr"
    assert vio.read(native).height == 1                   # auto dispatch
    assert vio.read(airr)[S.JUNCTION_AA].to_list() == ["CASSLF"]


def test_read_metadata_nan_to_null(tmp_path):
    p = tmp_path / "meta.tsv"
    p.write_text("sample_name\tstatus\tfraction\nS1\tas\tnan\nS2\thd\tcd8\n")
    m = vio.read_metadata(p)
    assert m["fraction"].to_list() == [None, "cd8"]


def test_read_metadata_commented_header(tmp_path):
    # Some metadata sheets comment out the header line (``#file_name\t...``); the
    # leading ``#`` on the first column name must be stripped.
    p = tmp_path / "meta.tsv"
    p.write_text("#file_name\tsample_id\tage\nS1.txt\tS1\t30\nS2.txt\tS2\t70\n")
    m = vio.read_metadata(p)
    assert m.columns[0] == "file_name"
    assert m["sample_id"].to_list() == ["S1", "S2"]


def test_read_samples_batch(tmp_path):
    _write_native(tmp_path / "S1.tsv.gz", [(4, "TGT", "CASSL", "TRBV1", ".", "TRBJ1")])
    _write_native(tmp_path / "S2.tsv.gz", [(6, "TGT", "CASSF", "TRBV2", ".", "TRBJ2"),
                                           (2, "TGA", "CASSY", "TRBV2", ".", "TRBJ2")])
    meta = pl.DataFrame({"sample_name": ["S1", "S2"], "disease_status": ["as", "hd"],
                         "b27": ["pos", "neg"]})
    long = vio.read_samples(meta, tmp_path)
    assert long.height == 3
    assert set(long["sample_id"].unique()) == {"S1", "S2"}
    assert "file_name" in long.columns
    # metadata joined correctly per sample
    assert long.filter(pl.col("sample_id") == "S1")["disease_status"].unique().to_list() == ["as"]
    assert long.filter(pl.col("sample_id") == "S2")["b27"].unique().to_list() == ["neg"]

    d = vio.read_samples(meta, tmp_path, as_dict=True)
    assert set(d) == {"S1", "S2"} and d["S2"].height == 2
