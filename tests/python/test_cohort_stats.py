"""Cohort summary stats: streaming + parallel + fused single-pass parity.

Pins the three properties the lazy/parallel cohort path must hold:

* :func:`vdjtools.io.map_samples` reads + reduces in parallel but returns results in
  **input (metadata) order**, thread-count-invariant, and equal to a serial loop.
* the CLI ``--threads`` knob never changes output; the ``overlap`` command's
  pre-aggregation is bit-identical to the per-pair :func:`overlap_metrics`.
* the fused ``--cohort`` / ``scan_cohort`` pass equals the per-sample path — **exact**
  for integer-weighted stats and :func:`diversity_cohort` (count-order canonicalised),
  and to ``rtol`` for float (``freq``) weights.
"""
import gzip

import numpy as np
import polars as pl
import pytest
from typer.testing import CliRunner

from vdjtools import stats
from vdjtools.cli import app
from vdjtools.io import ingest_cohort, map_samples, read_metadata, scan_cohort
from vdjtools.io.batch import read

# Deliberately non-lexicographic sample order (s10 between s2 and s1) to stress ordering.
_SAMPLES = {
    "s2": [(10, "TGTGCT", "CASSF", "TRBV12-3*01", ".", "TRBJ2-1*01"),
           (1, "TGA", "CASF", "TRBV12-3*01", ".", "TRBJ2-1*01"),
           (1, "TGC", "CAF", "TRBV7-2*01", ".", "TRBJ2-1*01"),
           (2, "TGT", "CAK", "TRBV7-2*01", ".", "TRBJ1-1*01")],
    "s10": [(5, "TGTGCT", "CASSF", "TRBV12-3*01", ".", "TRBJ2-1*01"),
            (5, "TGGGGT", "CASSY", "TRBV20-1*01", ".", "TRBJ2-1*01"),
            (3, "TGT", "CAT", "TRBV20-1*01", ".", "TRBJ2-1*01")],
    "s1": [(3, "TGT", "CASSF", "TRBV7-2*01", ".", "TRBJ1-1*01"),
           (3, "TGT", "CATG", "TRBV7-2*01", ".", "TRBJ1-1*01"),
           (4, "TGT", "CATGG", "TRBV7-2*01", ".", "TRBJ1-1*01"),
           (1, "TGT", "CA", "TRBV7-2*01", ".", "TRBJ1-1*01")],
}


@pytest.fixture
def cohort(tmp_path):
    """Write the sample .tsv.gz files + a metadata sheet + an ingested parquet cohort."""
    for name, rows in _SAMPLES.items():
        body = "".join(f"{c}\t0\t{nt}\t{aa}\t{v}\t{dd}\t{j}\n" for c, nt, aa, v, dd, j in rows)
        with gzip.open(tmp_path / f"{name}.tsv.gz", "wt") as f:
            f.write("count\tfreq\tcdr3nt\tcdr3aa\tv\td\tj\n" + body)
    md_path = tmp_path / "metadata.tsv"
    pl.DataFrame({"sample_name": list(_SAMPLES)}).write_csv(md_path, separator="\t")
    cohort_dir = tmp_path / "cohort"
    ingest_cohort(read_metadata(md_path), tmp_path, cohort_dir,
                  sample_col="sample_name", file_template="{sample}.tsv.gz")
    frames = {n: read(tmp_path / f"{n}.tsv.gz") for n in _SAMPLES}
    return tmp_path, md_path, cohort_dir, frames


# --------------------------------------------------------------------------- map_samples
def test_map_samples_preserves_input_order_and_matches_serial(cohort):
    base, _, _, frames = cohort
    items = [(n, base / f"{n}.tsv.gz") for n in _SAMPLES]  # input order s2, s10, s1
    serial = [(n, stats.diversity_stats(frames[n])) for n in _SAMPLES]
    for workers in (1, 2, 8):
        got = map_samples(stats.diversity_stats, items, workers=workers)
        assert [sid for sid, _ in got] == list(_SAMPLES)            # input order kept
        for (sid, res), (esid, eres) in zip(got, serial):
            assert sid == esid and res.equals(eres)


# --------------------------------------------------------------------------- CLI --threads
@pytest.mark.parametrize("cmd,extra", [
    ("diversity", []), ("spectratype", []),
    ("segment-usage", ["--segment", "v"]), ("overlap", []),
])
def test_cli_threads_invariant(cohort, cmd, extra):
    base, _, _, _ = cohort
    files = [str(base / f"{n}.tsv.gz") for n in _SAMPLES]
    runner = CliRunner()

    def run(args):
        r = runner.invoke(app, [cmd, *files, *extra, *args])
        assert r.exit_code == 0, (r.output, r.exception)
        return r.stdout

    ref = run([])
    for t in ("1", "2", "8"):
        assert run(["--threads", t]) == ref


def test_cli_overlap_preagg_matches_overlap_metrics(cohort):
    base, _, _, frames = cohort
    from vdjtools.overlap.metrics import overlap_metrics
    ids = list(_SAMPLES)
    expect = [{"sample_a": ids[i], "sample_b": ids[k], **overlap_metrics(frames[ids[i]], frames[ids[k]])}
              for i in range(len(ids)) for k in range(i + 1, len(ids))]
    expect_df = pl.DataFrame(expect)
    files = [str(base / f"{n}.tsv.gz") for n in _SAMPLES]
    r = CliRunner().invoke(app, ["overlap", *files])
    assert r.exit_code == 0, (r.output, r.exception)
    got = pl.read_csv(r.stdout.encode(), separator="\t")
    cols = ["sample_a", "sample_b", "D", "F", "F2", "d1", "d2", "d12"]
    assert got.select(cols).equals(expect_df.select(cols))


# --------------------------------------------------------------------------- fused cohort
def _serial(fn, frames, **kw):
    return (pl.concat([fn(f, **kw).select(pl.lit(n).alias("sample_id"), pl.all())
                       for n, f in frames.items()], how="vertical_relaxed")
            .sort(pl.all()))


def test_diversity_cohort_exact(cohort):
    _, _, cohort_dir, frames = cohort
    serial = _serial(stats.diversity_stats, frames).sort("sample_id")
    fused = stats.diversity_cohort(scan_cohort(cohort_dir, join_metadata=False)).sort("sample_id")
    assert serial.equals(fused)                        # bit-exact (count order canonicalised)


@pytest.mark.parametrize("fn_name,kw", [
    ("spectratype", dict(kind="aa", weight="reads")),
    ("spectratype", dict(kind="aa", weight="unique")),
    ("segment_usage", dict(segment="v", weight="reads")),
    ("segment_usage", dict(segment="j", weight="unique")),
    ("vj_usage", dict(weight="reads")),
])
def test_fused_group_stats_exact_int_weights(cohort, fn_name, kw):
    _, _, cohort_dir, frames = cohort
    fn = getattr(stats, fn_name)
    serial = _serial(fn, frames, **kw)
    fused = (fn(scan_cohort(cohort_dir, join_metadata=False), by=["sample_id"], **kw)
             .collect(engine="streaming").sort(pl.all()))
    assert serial.equals(fused)


def test_fused_freq_weight_matches_to_rtol(cohort):
    _, _, cohort_dir, frames = cohort
    kw = dict(kind="aa", weight="freq")
    key = ["sample_id", "locus", "length"]
    serial = _serial(stats.spectratype, frames, **kw).sort(key)
    fused = (stats.spectratype(scan_cohort(cohort_dir, join_metadata=False), by=["sample_id"], **kw)
             .collect(engine="streaming").sort(key))
    assert serial.select(key).equals(fused.select(key))
    assert np.allclose(serial["weight"].to_numpy(), fused["weight"].to_numpy(), rtol=1e-12, atol=0)


@pytest.mark.parametrize("cmd,extra", [
    ("diversity", []), ("spectratype", []), ("segment-usage", ["--segment", "v"]),
])
def test_cli_cohort_matches_serial(cohort, cmd, extra):
    base, _, cohort_dir, _ = cohort
    files = [str(base / f"{n}.tsv.gz") for n in _SAMPLES]
    runner = CliRunner()

    def run(args):
        r = runner.invoke(app, [cmd, *args, *extra])
        assert r.exit_code == 0, (r.output, r.exception)
        return pl.read_csv(r.stdout.encode(), separator="\t").sort(pl.all())

    assert run(files).equals(run(["--cohort", str(cohort_dir)]))


def test_by_default_is_byte_identical_eager(cohort):
    """The additive ``by=()`` default must not change the per-sample eager output."""
    _, _, _, frames = cohort
    f = next(iter(frames.values()))
    assert stats.spectratype(f).equals(stats.spectratype(f, by=[]))
    assert stats.segment_usage(f, "v").equals(stats.segment_usage(f, "v", by=[]))
