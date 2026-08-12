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


# --- toy germline / models -------------------------------------------------------------------
# A hand-built 3-gene locus. Small enough that Pgen, the free-parameter count and EM all run in
# milliseconds and can be reasoned about by hand, which is what makes it the default fixture for
# the model-workshop tests: no OLGA, no arda, no network, no bundled-model assumptions.

#: V germline must start on the conserved Cys codon (TGT/TGC) and J must end on Phe/Trp.
TOY_V = {"TOYV1*01": "TGTGCCAGC", "TOYV1*02": "TGTGCCAGT", "TOYV2*01": "TGTGCTTCC"}
TOY_J = {"TOYJ1*01": "AACTATGGCTATACCTTT", "TOYJ2*01": "AACGAGCAGTTT"}
TOY_D = {"TOYD1*01": "GGGACAGGGGGC"}


def _germline_rows(mapping, segment):
    return [{"allele": a, "gene": a.split("*")[0], "segment": segment, "sequence": s,
             "cdr3_anchor": -1, "functional": True, "full_germline": ""}
            for a, s in mapping.items()]


@pytest.fixture
def toy_germline():
    """A minimal valid VJ germline frame (3 V alleles over 2 genes, 2 J alleles)."""
    import polars as pl

    return pl.DataFrame(_germline_rows(TOY_V, "V") + _germline_rows(TOY_J, "J"))


@pytest.fixture
def toy_germline_vdj(toy_germline):
    """The toy germline plus one D allele, which makes any model built from it VDJ."""
    import polars as pl

    return pl.concat([toy_germline, pl.DataFrame(_germline_rows(TOY_D, "D"))])


@pytest.fixture
def toy_germline_more_alleles(toy_germline):
    """The toy germline plus one extra allele of an EXISTING gene (no new gene)."""
    import polars as pl

    return pl.concat([toy_germline,
                      pl.DataFrame(_germline_rows({"TOYV1*03": "TGTGCCAGA"}, "V"))])


@pytest.fixture
def toy_germline_extended(toy_germline):
    """The toy germline plus a new allele of a known gene and a whole new gene."""
    import polars as pl

    extra = {"TOYV1*03": "TGTGCCAGA", "TOYV3*01": "TGTTGGGGA"}
    return pl.concat([toy_germline, pl.DataFrame(_germline_rows(extra, "V"))])


@pytest.fixture
def toy_model(toy_germline):
    """A validated VJ model over the toy germline, with placeholder marginals."""
    from vdjtools.model.io import from_germline

    return from_germline(toy_germline, locus="TOY", ins_max=3)


@pytest.fixture
def toy_model_vdj(toy_germline_vdj):
    """A validated VDJ model over the toy germline."""
    from vdjtools.model.io import from_germline

    return from_germline(toy_germline_vdj, locus="TOY", ins_max=3)


@pytest.fixture(scope="module")
def small_model():
    """The smallest real bundled model (human TRG) — for tests that need a realistic one."""
    from vdjtools.model import load_bundled

    return load_bundled("TRG", "olga")
