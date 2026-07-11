"""Conformance test for the model-reads bootstrap (private HF dataset isalgo/airr_model_read).

Skips cleanly until the owner populates the dataset (and a HF token is present): the private
repo is empty / unauthorized offline, so ``list_reads`` / ``fetch`` raise and we skip. Once
populated this asserts the fetched reads map to the canonical clonotype schema and yield CDR3
nt strings ready for EM — the same graceful-offline contract as ``test_io_hf.py``.
"""
import pytest

from vdjtools.io import schema as S
from vdjtools.model import data


def _reads():
    """List the dataset's files; skip if huggingface_hub missing or the repo is empty/unreachable."""
    pytest.importorskip("huggingface_hub")
    try:
        files = data.list_reads()
    except Exception as e:  # offline / auth / repo not found
        pytest.skip(f"model-reads dataset unavailable ({data.MODEL_READS_REPO}): {e}")
    if not files:
        pytest.skip(f"model-reads dataset {data.MODEL_READS_REPO} not populated yet")
    return files


def test_model_reads_map_to_canonical():
    files = _reads()
    df = data.load_reads(files[0], n_rows=10_000)
    assert df.height > 0
    for col, dtype in S.SCHEMA.items():
        assert col in df.columns and df[col].dtype == dtype


def test_model_reads_yield_training_cdr3():
    files = _reads()
    df = data.load_reads(files[0], n_rows=10_000)
    seqs = data.training_cdr3(df, out_of_frame=False)
    assert seqs and all(set(s.upper()) <= set("ACGTN") for s in seqs[:100])
