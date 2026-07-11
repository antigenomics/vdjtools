"""Shared pytest fixtures.

The ``hf`` fixture fetches a file from a HuggingFace dataset into the local HF cache.
It skips cleanly (never fails) when ``huggingface_hub`` is not installed or the fetch
fails (offline / network / auth), so the default ``pytest tests/python`` stays green
with no network and ``huggingface_hub`` never becomes a runtime dependency.
"""
import pytest


def _hf(repo, filename):
    """Download ``filename`` from dataset ``repo``; skip the test if unavailable."""
    hub = pytest.importorskip("huggingface_hub")
    try:
        return hub.hf_hub_download(repo_id=repo, filename=filename, repo_type="dataset")
    except Exception as e:  # offline / network / auth
        pytest.skip(f"HF fetch failed ({repo}/{filename}): {e}")


@pytest.fixture
def hf():
    """Return the HuggingFace dataset-file fetch helper."""
    return _hf
