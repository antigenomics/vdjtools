"""Raw streaming must preserve order, identifiers and missing values across batches."""
import gzip
import polars as pl
import pytest
from vdjtools.io import iter_table_batches


def test_raw_batches(tmp_path):
    raw = "sequence_id\tjunction_aa\tcount\n001\tCASSF\t2\n002\tCATSF\t3\n003\t\t0\n004\tCASSF\t5\n"
    path = tmp_path / "sample.tsv.gz"
    with gzip.open(path, "wt") as stream:
        stream.write(raw)
    batches = list(iter_table_batches(path, batch_size=2))
    assert all(frame.height <= 2 for frame in batches)
    expected = pl.read_csv(path, separator="\t", infer_schema_length=0)
    assert pl.concat(batches).equals(expected)
    parquet = tmp_path / "sample.parquet"
    expected.write_parquet(parquet, row_group_size=1)
    assert pl.concat(iter_table_batches(parquet, batch_size=2, columns=["sequence_id"])).equals(expected.select("sequence_id"))
    for size in (0, -1, True):
        with pytest.raises(ValueError, match="positive integer"):
            list(iter_table_batches(path, batch_size=size))
