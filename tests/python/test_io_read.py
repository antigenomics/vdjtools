

# ------------------------------------------- frequencies as in the file (3.10.0, 2026-08-20)


def _airr(tmp_path, freq=True):
    import polars as pl
    d = {"v_call": ["TRBV20-1"] * 2, "j_call": ["TRBJ2-2"] * 2,
         "junction_aa": ["CASSIRSSYEQYF", "CASSLRSSYEQYF"], "duplicate_count": [90, 10]}
    if freq:
        d["frequency"] = [0.5, 0.5]          # deliberately NOT count/total
    p = tmp_path / "s.tsv"
    pl.DataFrame(d).write_csv(p, separator="\t")
    return p


def test_frequency_is_recomputed_by_default(tmp_path):
    from vdjtools import io
    assert io.read(_airr(tmp_path), fmt="airr")["frequency"].to_list() == [0.9, 0.1]


def test_the_files_frequency_survives_when_asked(tmp_path):
    """A UMI-corrected frequency is not count/total and must not be silently replaced by it."""
    from vdjtools import io
    got = io.read(_airr(tmp_path), fmt="airr", recompute_frequencies=False)["frequency"].to_list()
    assert got == [0.5, 0.5]


def test_a_file_without_a_frequency_column_still_gets_one(tmp_path):
    from vdjtools import io
    got = io.read(_airr(tmp_path, freq=False), fmt="airr",
                  recompute_frequencies=False)["frequency"].to_list()
    assert got == [0.9, 0.1]
