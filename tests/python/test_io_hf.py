"""Real-data conformance tests bootstrapped from HuggingFace datasets.

All tests fetch via the ``hf`` fixture and skip cleanly offline / without
``huggingface_hub`` (see conftest). Data (owner isalgo, cc-by-nc-nd) is NOT
redistributed here — only fetched at test time into the HF cache.
"""
from pathlib import Path

import polars as pl
import pytest

from vdjtools import io as vio
from vdjtools.io import schema as S

ANKSPOND = "isalgo/airr_ankspond"
CONTROL = "isalgo/airr_control"

# The three smallest ankspond samples (native path new/<sample_name>.tsv.gz).
SMALL_SAMPLES = ["as_Dv_SFCD8", "as_Mart_PB_F", "as_Chaad_PB_F"]


def test_ankspond_sample_maps_to_canonical(hf):
    # NOTE: ankspond `new/` files are an AIRR-hybrid (count junction_aa v_call j_call
    # locus), not the classic native vdjtools table — auto-detect routes them to AIRR.
    path = hf(ANKSPOND, f"new/{SMALL_SAMPLES[0]}.tsv.gz")
    assert vio.sniff_format(path) == "airr"
    df = vio.read(path)
    assert df.height > 0
    for col, dtype in S.SCHEMA.items():
        assert col in df.columns and df[col].dtype == dtype
    assert df[S.LOCUS].unique().to_list() == ["TRB"]      # locus derived from v_call
    assert df[S.COUNT].min() >= 1


def test_ankspond_read_samples_with_metadata(hf, tmp_path):
    meta_path = hf(ANKSPOND, "new/metadata.tsv")
    # fetch the three small samples into the (shared) HF cache dir
    sample_paths = [hf(ANKSPOND, f"new/{s}.tsv.gz") for s in SMALL_SAMPLES]
    base_dir = Path(sample_paths[0]).parent            # the cached new/ directory

    full = vio.read_metadata(meta_path)
    subset = full.filter(pl.col("sample_name").is_in(SMALL_SAMPLES))
    assert subset.height == len(SMALL_SAMPLES)

    long = vio.read_samples(subset, base_dir, sample_col="sample_name",
                            file_template="{sample}.tsv.gz")
    # canonical schema + reserved + joined metadata columns
    for col, dtype in S.SCHEMA.items():
        assert col in long.columns and long[col].dtype == dtype
    assert "sample_id" in long.columns and "file_name" in long.columns
    assert set(long["sample_id"].unique()) == set(SMALL_SAMPLES)
    assert set(long["disease_status"].unique()) == {"as"}
    assert set(long["b27"].unique()) == {"pos"}

    # per-sample row counts > 0
    per = long.group_by("sample_id").len()
    assert (per["len"] > 0).all()

    # "nan" -> null works upstream; here the small samples carry real fraction values
    frac = {r["sample_name"]: r["fraction"] for r in subset.iter_rows(named=True)}
    assert frac["as_Dv_SFCD8"] == "cd8" and frac["as_Mart_PB_F"] == "bulk"


def test_ankspond_as_dict(hf):
    sample_paths = [hf(ANKSPOND, f"new/{s}.tsv.gz") for s in SMALL_SAMPLES]
    base_dir = Path(sample_paths[0]).parent
    meta = pl.DataFrame({"sample_name": SMALL_SAMPLES})
    d = vio.read_samples(meta, base_dir, add_metadata=False, as_dict=True)
    assert set(d) == set(SMALL_SAMPLES)
    assert all(frame.height > 0 for frame in d.values())


def test_control_native_schema_capped(hf):
    # airr_control is exactly the native vdjtools schema + extra annotation columns;
    # 56M rows total, so read a capped preview and assert clean canonical mapping.
    hub = pytest.importorskip("huggingface_hub")
    files = hub.list_repo_files(repo_id=CONTROL, repo_type="dataset")
    tsv = next((f for f in files if f.endswith(".vdjtools.tsv.gz")), None)
    if tsv is None:
        pytest.skip(f"no vdjtools tsv in {CONTROL}: {files}")
    path = hf(CONTROL, tsv)
    assert vio.sniff_format(path) == "vdjtools"
    df = vio.read(path, n_rows=50_000)
    assert df.height == 50_000
    for col, dtype in S.SCHEMA.items():
        assert col in df.columns and df[col].dtype == dtype
    assert df[S.JUNCTION_AA].null_count() == 0
    assert df[S.COUNT].min() >= 1
    assert df[S.LOCUS].str.starts_with("TR").all() or df[S.LOCUS].str.starts_with("IG").all()
