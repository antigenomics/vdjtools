"""Bootstrap the recombination model from real AIRR reads on HuggingFace.

The engine is trained (EM; :func:`vdjtools.model.infer.infer` / :func:`infer_native`) on
out-of-frame nucleotide CDR3 reads. Until owner data lands these are OLGA-synthetic draws (see
``SOURCES.md``); the real reads live in the **private** HF dataset ``isalgo/airr_model_read``
(owner-populated). This module fetches them into the canonical clonotype frame with the same
graceful-offline behaviour as the rest of the HF-backed IO — ``huggingface_hub`` is a lazy,
test/fetch-only import, never a runtime dependency of the model math.

Once the dataset is populated::

    from vdjtools.model import from_olga, data
    from vdjtools.model.infer import infer_native
    template = from_olga(olga_dir, locus="TRB")          # gene set + germline + event graph
    df = data.load_reads(data.list_reads()[0], locus="TRB")
    fit, report = infer_native(template, data.training_cdr3(df))
"""
from __future__ import annotations

import polars as pl

from .. import io as vio
from ..io import schema as S

#: The private HuggingFace dataset of real AIRR reads used to train the models (owner: isalgo,
#: cc-by-nc; populated by the owner). Not redistributed here — only fetched into the HF cache.
MODEL_READS_REPO = "isalgo/airr_model_read"

_DATA_SUFFIXES = (".parquet", ".tsv", ".tsv.gz", ".txt.gz")


def list_reads(repo: str = MODEL_READS_REPO) -> list[str]:
    """List the data files in the model-reads dataset (empty until the owner populates it).

    Requires network + (for a private repo) a HuggingFace token in the environment.
    """
    from huggingface_hub import list_repo_files

    return [f for f in list_repo_files(repo_id=repo, repo_type="dataset") if f.endswith(_DATA_SUFFIXES)]


def fetch(filename: str, repo: str = MODEL_READS_REPO) -> str:
    """Download one file from the model-reads dataset into the HF cache; return its local path."""
    from huggingface_hub import hf_hub_download

    return hf_hub_download(repo_id=repo, filename=filename, repo_type="dataset")


def load_reads(
    filename: str,
    *,
    repo: str = MODEL_READS_REPO,
    locus: str | None = None,
    n_rows: int | None = None,
) -> pl.DataFrame:
    """Fetch and read one model-reads file into the canonical clonotype frame.

    Args:
        filename: Path of the file within the dataset (see :func:`list_reads`).
        repo: HuggingFace dataset id (default the private :data:`MODEL_READS_REPO`).
        locus: If given, keep only rows of this locus (e.g. ``"TRB"``); the locus is derived
            from ``v_call`` by the reader.
        n_rows: Optional row cap for previewing large files.

    Returns:
        The canonical clonotype :class:`polars.DataFrame` (see :mod:`vdjtools.io.schema`).
    """
    df = vio.read(fetch(filename, repo), n_rows=n_rows)
    if locus is not None:
        df = df.filter(pl.col(S.LOCUS) == locus)
    return df


def training_cdr3(df: pl.DataFrame, *, out_of_frame: bool = True) -> list[str]:
    """CDR3 nucleotide strings for EM (:func:`vdjtools.model.infer.infer`).

    Args:
        df: A canonical clonotype frame (e.g. from :func:`load_reads`).
        out_of_frame: Keep only frameshifted reads (CDR3 nt length not a multiple of 3) — the
            unbiased training set the E-step assumes. Set ``False`` if the dataset is already
            curated to non-productive reads, or to use every read as-is.

    Returns:
        List of CDR3 nucleotide strings (nulls dropped).
    """
    nt = df[S.CDR3_NT].drop_nulls()
    if out_of_frame:
        nt = nt.filter(nt.str.len_bytes() % 3 != 0)
    return nt.to_list()
