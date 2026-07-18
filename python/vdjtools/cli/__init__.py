"""vdjtools command-line interface (typer) — one ``vdjtools`` entry point.

Two families of commands:

* **Model engine** (OLGA/IGoR-style, on the native recombination core + built-in models):
  ``pgen`` (generation probability), ``generate`` (sample sequences), ``models`` (list built-ins).
* **Repertoire analytics** (vanilla-vdjtools-style, over sample files or a metadata table):
  ``diversity``, ``spectratype``, ``segment-usage``, ``overlap``.

Analytics commands take either a list of sample files or ``-m/--metadata <table>`` (+ ``--base-dir``),
mirroring the metadata-driven workflow of the legacy tool. Every command writes a TSV to ``-o`` or,
by default, to stdout (progress/errors go to stderr), so commands pipe cleanly.
"""
from __future__ import annotations

import functools
import sys
from pathlib import Path
from typing import Optional

import polars as pl
import typer

app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help="vdjtools — TCR/BCR repertoire analysis: Pgen, sequence generation, diversity, overlap.",
)

_SEQ_COLS = ("junction_aa", "junction_nt", "cdr3_aa", "cdr3_nt", "cdr3aa", "cdr3nt", "cdr3", "junction")


# ---------------------------------------------------------------------------- helpers
def _err(msg: str) -> None:
    typer.secho(f"error: {msg}", fg=typer.colors.RED, err=True)
    raise typer.Exit(1)


def _info(msg: str) -> None:
    typer.secho(msg, fg=typer.colors.GREEN, err=True)


def _load_model(model: Optional[str], source: str, model_path: Optional[Path]):
    from vdjtools.model import load_bundled, load_model

    if model_path is not None:
        return load_model(model_path)
    if model is None:
        _err("give a model: --model <LOCUS> (built-in) or --model-path <DIR>")
    try:
        return load_bundled(model, source=source)
    except (ValueError, FileNotFoundError) as e:  # unknown source / locus not bundled
        _err(str(e))


def _write(df: pl.DataFrame, out: Optional[Path]) -> None:
    if out is None:
        sys.stdout.write(df.write_csv(separator="\t"))
    else:
        df.write_csv(out, separator="\t")
        _info(f"wrote {df.height} rows → {out}")


def _is_nt(seq: str) -> bool:
    return bool(seq) and set(seq.upper()) <= set("ACGT")


def _read_seq_table(path: Path, column: Optional[str], no_header: bool) -> tuple[pl.DataFrame, str]:
    """Read a sequence table; return (frame, sequence-column-name). Tab-separated, strings only."""
    if no_header:
        df = pl.read_csv(path, separator="\t", has_header=False, infer_schema_length=0)
        return df, df.columns[0]
    df = pl.read_csv(path, separator="\t", infer_schema_length=0)
    if column:
        if column not in df.columns:
            _err(f"column {column!r} not found; columns: {df.columns}")
        return df, column
    low = {c.lower(): c for c in df.columns}
    col = next((low[c] for c in _SEQ_COLS if c in low), df.columns[0])
    return df, col


def _sample_items(samples, metadata, base_dir, sample_col, file_template) -> list:
    """[(sample_id, path), ...] — from a metadata table if given, else positional files.

    Returns paths, not loaded frames: the caller streams them through
    :func:`vdjtools.io.map_samples` so only ``O(workers)`` samples are ever resident.
    """
    from vdjtools.io.batch import read_metadata

    if metadata is not None:
        md = read_metadata(metadata)
        if sample_col not in md.columns:
            _err(f"--sample-col {sample_col!r} not in metadata columns: {md.columns}")
        base = base_dir or Path(metadata).parent
        return [(str(r[sample_col]), base / file_template.format(sample=r[sample_col]))
                for r in md.iter_rows(named=True)]
    if not samples:
        _err("give sample files as arguments, or -m/--metadata <table> (with --base-dir)")
    return [(Path(s).name.split(".")[0], s) for s in samples]


def _tag(df: pl.DataFrame, sample_id: str) -> pl.DataFrame:
    return df.with_columns(pl.lit(sample_id).alias("sample_id")).select(
        ["sample_id", *[c for c in df.columns if c != "sample_id"]]
    )


# ---------------------------------------------------------------------------- model commands
@app.command()
def models() -> None:
    """List the precomputed recombination models shipped with the package."""
    from vdjtools.model import list_bundled

    for src, loci in list_bundled().items():
        typer.echo(f"{src:8s} {' '.join(loci) if loci else '(none)'}")


@app.command()
def pgen(
    input: Path = typer.Argument(..., help="Table (TSV) or list of CDR3 sequences."),
    model: Optional[str] = typer.Option(None, "--model", "-m", help="Built-in locus: TRA TRB TRG TRD IGH IGK IGL."),
    source: str = typer.Option("olga", help="Built-in model set: olga | learned."),
    model_path: Optional[Path] = typer.Option(None, help="Load a custom model directory instead of a built-in."),
    column: Optional[str] = typer.Option(None, "--column", "-c", help="Sequence column (default: auto-detect / first)."),
    v_col: Optional[str] = typer.Option(None, "--v-col", help="V-allele column to condition on (default: marginalize)."),
    j_col: Optional[str] = typer.Option(None, "--j-col", help="J-allele column to condition on (default: marginalize)."),
    seq_type: str = typer.Option("auto", "--type", help="auto | aa | nt."),
    mismatches: int = typer.Option(0, help="Amino-acid Hamming ball: 0 exact, 1 sums all single-substitution neighbours."),
    no_header: bool = typer.Option(False, "--no-header", help="Input is a bare sequence list (no header row)."),
    out: Optional[Path] = typer.Option(None, "--out", "-o", help="Output TSV (default: stdout)."),
) -> None:
    """Compute generation probability (Pgen) for CDR3 sequences — like ``olga-compute_pgen``.

    Appends a ``pgen`` column. V/J are marginalized unless ``--v-col``/``--j-col`` are given.
    Nucleotide vs amino-acid is auto-detected per sequence; amino-acid input can also sum the
    Hamming-distance-1 ball with ``--mismatches 1`` (fast, native).
    """
    from vdjtools.model import native

    if seq_type not in ("auto", "aa", "nt"):
        _err("--type must be auto, aa or nt")
    m = _load_model(model, source, model_path)
    df, seqcol = _read_seq_table(input, column, no_header)
    seqs = df[seqcol].to_list()
    vs = df[v_col].to_list() if v_col and v_col in df.columns else [None] * len(seqs)
    js = df[j_col].to_list() if j_col and j_col in df.columns else [None] * len(seqs)

    pg: list[float] = []
    for s, v, j in zip(seqs, vs, js):
        if s is None or s == "":
            pg.append(0.0)
            continue
        nt = _is_nt(s) if seq_type == "auto" else (seq_type == "nt")
        if nt:
            pg.append(native.pgen_nt(m, s, v, j))
        else:
            pg.append(native.pgen_aa(m, s, v, j, mismatches=mismatches))
    _write(df.with_columns(pl.Series("pgen", pg)), out)


@app.command()
def generate(
    model: Optional[str] = typer.Option(None, "--model", "-m", help="Built-in locus: TRA TRB TRG TRD IGH IGK IGL."),
    source: str = typer.Option("olga", help="Built-in model set: olga | learned."),
    model_path: Optional[Path] = typer.Option(None, help="Load a custom model directory instead of a built-in."),
    n: int = typer.Option(100, "--number", "-n", help="Number of sequences to generate."),
    seed: Optional[int] = typer.Option(None, help="Random seed for reproducible draws."),
    productive: bool = typer.Option(False, help="Only keep in-frame, stop-free (productive) rearrangements."),
    out: Optional[Path] = typer.Option(None, "--out", "-o", help="Output TSV (default: stdout)."),
) -> None:
    """Sample recombined sequences from a model — like ``olga-generate_sequences``.

    Emits ``junction_nt, junction_aa, v_call, d_call, d2_call, j_call, productive`` (``d2_call`` is the
    tandem D on the learned D-bearing loci; null otherwise).
    """
    from vdjtools.model.generate import generate as _generate

    m = _load_model(model, source, model_path)
    _write(_generate(m, n, seed=seed, productive_only=productive), out)


# ---------------------------------------------------------------------------- analytics commands
_SAMPLES = typer.Argument(None, help="Clonotype sample files (native vdjtools or AIRR).")
_META = typer.Option(None, "--metadata", "-m", help="Metadata table (one row per sample).")
_BASE = typer.Option(None, "--base-dir", help="Directory holding the sample files (with -m).")
_SCOL = typer.Option("sample_name", "--sample-col", help="Metadata column with the sample name.")
_TMPL = typer.Option("{sample}.tsv.gz", "--file-template", help="Sample filename template (with -m).")
_FMT = typer.Option("auto", "--format", help="auto | vdjtools | airr.")
_OUT = typer.Option(None, "--out", "-o", help="Output TSV (default: stdout).")
_THREADS = typer.Option(
    0, "--threads", "-t",
    help="Worker threads over samples (0 = all cores). Lower to core count if compute-bound.",
)
_COHORT = typer.Option(
    None, "--cohort",
    help="Pre-ingested parquet cohort dir (vdjtools.io.ingest_cohort): one streamed "
         "out-of-core pass over the whole cohort instead of per-sample files.",
)


@app.command()
def diversity(
    samples: Optional[list[Path]] = _SAMPLES, metadata: Optional[Path] = _META,
    base_dir: Optional[Path] = _BASE, sample_col: str = _SCOL, file_template: str = _TMPL,
    fmt: str = _FMT, threads: int = _THREADS, cohort: Optional[Path] = _COHORT,
    out: Optional[Path] = _OUT,
) -> None:
    """Per-sample diversity (observed richness, Chao, Efron-Thisted, Shannon, Simpson, d50)."""
    from vdjtools.io.batch import map_samples
    from vdjtools.stats.diversity import diversity_cohort, diversity_stats

    if cohort is not None:
        from vdjtools.io.cohort import scan_cohort
        _write(diversity_cohort(scan_cohort(cohort, join_metadata=False)), out)
        return
    items = _sample_items(samples, metadata, base_dir, sample_col, file_template)
    rows = [_tag(res, sid) for sid, res in
            map_samples(diversity_stats, items, fmt=fmt, workers=threads or None)]
    _write(pl.concat(rows, how="vertical_relaxed"), out)


@app.command()
def spectratype(
    samples: Optional[list[Path]] = _SAMPLES, metadata: Optional[Path] = _META,
    base_dir: Optional[Path] = _BASE, sample_col: str = _SCOL, file_template: str = _TMPL,
    fmt: str = _FMT, kind: str = typer.Option("aa", help="aa | nt (length unit)."),
    weight: str = typer.Option("reads", help="reads | unique | freq."),
    threads: int = _THREADS, cohort: Optional[Path] = _COHORT, out: Optional[Path] = _OUT,
) -> None:
    """Per-sample CDR3 length spectratype."""
    from vdjtools.io.batch import map_samples
    from vdjtools.stats.spectratype import spectratype as _spec

    if cohort is not None:
        from vdjtools.io.cohort import scan_cohort
        _write(_spec(scan_cohort(cohort, join_metadata=False), kind=kind, weight=weight,
                     by=["sample_id"]).collect(engine="streaming"), out)
        return
    items = _sample_items(samples, metadata, base_dir, sample_col, file_template)
    fn = functools.partial(_spec, kind=kind, weight=weight)
    rows = [_tag(res, sid) for sid, res in
            map_samples(fn, items, fmt=fmt, workers=threads or None)]
    _write(pl.concat(rows, how="vertical_relaxed"), out)


@app.command(name="segment-usage")
def segment_usage(
    samples: Optional[list[Path]] = _SAMPLES, metadata: Optional[Path] = _META,
    base_dir: Optional[Path] = _BASE, sample_col: str = _SCOL, file_template: str = _TMPL,
    fmt: str = _FMT, segment: str = typer.Option("v", help="v | d | j | c."),
    weight: str = typer.Option("reads", help="reads | unique | freq."),
    threads: int = _THREADS, cohort: Optional[Path] = _COHORT, out: Optional[Path] = _OUT,
) -> None:
    """Per-sample V/D/J/C segment usage."""
    from vdjtools.io.batch import map_samples
    from vdjtools.stats.usage import segment_usage as _usage

    if cohort is not None:
        from vdjtools.io.cohort import scan_cohort
        _write(_usage(scan_cohort(cohort, join_metadata=False), segment=segment,
                      weight=weight, by=["sample_id"]).collect(engine="streaming"), out)
        return
    items = _sample_items(samples, metadata, base_dir, sample_col, file_template)
    fn = functools.partial(_usage, segment=segment, weight=weight)
    rows = [_tag(res, sid) for sid, res in
            map_samples(fn, items, fmt=fmt, workers=threads or None)]
    _write(pl.concat(rows, how="vertical_relaxed"), out)


@app.command()
def overlap(
    samples: Optional[list[Path]] = _SAMPLES, metadata: Optional[Path] = _META,
    base_dir: Optional[Path] = _BASE, sample_col: str = _SCOL, file_template: str = _TMPL,
    fmt: str = _FMT, threads: int = _THREADS, out: Optional[Path] = _OUT,
) -> None:
    """Exact-match pairwise repertoire overlap (D, F, F2, R) for every sample pair."""
    from vdjtools.io.batch import map_samples
    from vdjtools.overlap.metrics import DEFAULT_KEY, _aggregate, _overlap_from_agg

    items = _sample_items(samples, metadata, base_dir, sample_col, file_template)
    if len(items) < 2:
        _err("overlap needs at least two samples")
    key = list(DEFAULT_KEY)
    # Aggregate each sample ONCE (streamed + parallel), then every O(n^2) pair is a join
    # over the pre-aggregated frames — not a re-aggregation of both raw frames per pair
    # (bitwise-identical to overlap_metrics, which cluster.pairwise_distances also reuses).
    aggs = map_samples(lambda df: _aggregate(df, key), items, fmt=fmt, workers=threads or None)
    rows = []
    for i in range(len(aggs)):
        for k in range(i + 1, len(aggs)):
            (a_id, a), (b_id, b) = aggs[i], aggs[k]
            rows.append({"sample_a": a_id, "sample_b": b_id,
                         **_overlap_from_agg(a, b, key)[1]})
    _write(pl.DataFrame(rows), out)


@app.command()
def dynamics(
    pre: Path = typer.Argument(..., help="Earlier sample (e.g. pre-vaccination)."),
    post: Path = typer.Argument(..., help="Later sample (e.g. post-vaccination)."),
    fmt: str = _FMT, out: Optional[Path] = _OUT,
    neff: Optional[float] = typer.Option(
        None, "--neff", help="Pin the pair's effective sample size (default: estimate it)."),
    umi: bool = typer.Option(
        False, "--umi", help="Counts are UMI/molecule counts, not reads: skip the downscale "
                             "(there is no oversampling to undo)."),
    min_total: int = typer.Option(
        6, "--min-total", help="Testability floor: combined downscaled count below this is "
                               "reported as `untested`, not as unchanged."),
    alpha: float = typer.Option(0.01, "--alpha", help="BH FDR threshold for calling a change."),
) -> None:
    """Paired within-donor test: which clonotypes changed between two timepoints.

    Classifies every clonotype as emergent / expanded / persistent / contracted / vanishing
    (or `untested`). Depth is handled PER PAIR via the effective sample size — never by
    normalising a cohort to a common depth, which is not a defined operation here.
    """
    from vdjtools.dynamics import test_pair
    from vdjtools.io.batch import read

    a, b = read(pre, fmt=fmt), read(post, fmt=fmt)
    try:
        res = test_pair(a, b, neff=None if umi else (neff if neff is not None else "auto"),
                        min_total=min_total, alpha=alpha)
    except ValueError as e:                      # pair too shallow / too few shared clonotypes
        _err(str(e))
    _info("  ".join(f"{k}={v}" for k, v in
                    sorted(res["dynamics"].value_counts().iter_rows())))
    _write(res, out)


@app.command()
def tcrnet(
    sample: Path = typer.Argument(..., help="Clonotype sample file."),
    fmt: str = _FMT, out: Optional[Path] = _OUT,
    locus: Optional[str] = typer.Option(None, "--locus", help="Force one locus (else per-locus)."),
    species: str = typer.Option("human", "--species"),
    scope: str = typer.Option("1,0,0,1", "--scope", help="Edit scope subs,ins,dels,total."),
    threads: int = typer.Option(0, "--threads", help="0 = all cores."),
) -> None:
    """Neighbourhood enrichment against a CONTROL REPERTOIRE (TCRnet).

    The control absorbs thymic selection and endemic-pathogen expansions, which a generation
    model cannot — at the cost of needing a large, HLA-matched cohort. See `alice` for the
    generative null. Neither can see a monoclonal expansion: enrichment measures breadth.
    """
    from vdjtools.io.batch import read
    from vdjtools.overlap import tcrnet as _tcrnet

    try:
        res = _tcrnet(read(sample, fmt=fmt), locus=locus, species=species, scope=scope,
                      threads=threads)
    except (ImportError, ValueError) as e:
        _err(str(e))
    _write(res, out)


@app.command()
def alice(
    sample: Path = typer.Argument(..., help="Clonotype sample file."),
    fmt: str = _FMT, out: Optional[Path] = _OUT,
    locus: Optional[str] = typer.Option(None, "--locus", help="Force one locus (else per-locus)."),
    source: str = typer.Option("olga", "--source", help="Bundled model source; leave on olga."),
    scope: str = typer.Option("1,0,0,1", "--scope", help="Edit scope subs,ins,dels,total."),
    selection_q: float = typer.Option(9.41, "--q", help="Thymic-selection factor Q."),
    min_degree: int = typer.Option(3, "--min-degree", help="Only test clonotypes with >= this "
                                                           "many neighbours (self included)."),
    min_count: int = typer.Option(2, "--min-count", help="Ignore clonotypes below this count."),
    threads: int = typer.Option(0, "--threads", help="0 = all cores."),
) -> None:
    """Neighbourhood enrichment against a V(D)J GENERATION MODEL (ALICE).

    Controls for the intrinsic biases of recombination, but knows nothing about selection or
    about which clonotypes are already common in people — the complement of `tcrnet`. Returns
    q_value and picks no threshold: the published ones differ 100-fold and were never
    reconciled.
    """
    from vdjtools.io.batch import read
    from vdjtools.overlap import alice as _alice

    try:
        res = _alice(read(sample, fmt=fmt), locus=locus, source=source, scope=scope,
                     selection_q=selection_q, min_degree=min_degree, min_count=min_count,
                     threads=threads)
    except (ImportError, KeyError, ValueError) as e:
        _err(str(e))
    _write(res, out)
