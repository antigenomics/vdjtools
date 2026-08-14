"""vdjtools command-line interface (typer) — one ``vdjtools`` entry point.

Command families:

* **Model engine** (OLGA/IGoR-style, on the native recombination core + built-in models):
  ``pgen`` (generation probability), ``generate`` (sample sequences), ``models`` (list built-ins).
* **Data** (format conversion + preprocessing): ``convert`` (any format → canonical TSV/Parquet),
  ``downsample``, ``filter`` (coding / frequency / segment), ``pool`` (pool or incidence-join samples).
* **Repertoire analytics** (vanilla-vdjtools-style, over sample files or a metadata table):
  ``diversity``, ``spectratype``, ``segment-usage``, ``overlap``.
* **Longitudinal & enrichment**: ``dynamics`` (paired within-donor expansion test), ``tcrnet`` /
  ``alice`` (neighbourhood enrichment vs a control cohort / a generation model).

Analytics commands take either a list of sample files or ``-m/--metadata <table>`` (+ ``--base-dir``),
mirroring the metadata-driven workflow of the legacy tool, and run in parallel over samples
(``-t/--threads``) or in one streamed pass over a pre-ingested Parquet cohort (``--cohort``). Every
command writes to ``-o`` — TSV, or Parquet when the path ends in ``.parquet`` / ``.pq`` — or, by
default, to stdout (progress/errors go to stderr), so commands pipe cleanly.
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
    """Write a result frame. ``.parquet`` / ``.pq`` → Parquet; anything else (or stdout) → TSV."""
    if out is None:
        sys.stdout.write(df.write_csv(separator="\t"))
        return
    if out.suffix.lower() in (".parquet", ".pq"):
        df.write_parquet(out)
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


# ---------------------------------------------------------------------------- data commands
@app.command()
def convert(
    input: Path = typer.Argument(..., help="Clonotype file in any supported format."),
    fmt: str = _FMT, out: Optional[Path] = _OUT,
) -> None:
    """Read any supported format and write the canonical AIRR-junction table.

    Auto-detects native vdjtools / AIRR / Parquet and the third-party exports (MiXcr, MiGec,
    MiTCR, immunoSEQ, IMGT/HighV-QUEST, Vidjil, RTCR, TRUST4, arda). Output is TSV, or Parquet
    when ``-o`` ends in ``.parquet`` / ``.pq`` — the typed, columnar, at-scale store.
    """
    from vdjtools.io.batch import read

    try:
        _write(read(input, fmt=fmt), out)
    except ValueError as e:                          # unknown fmt / unrecognised header
        _err(str(e))


@app.command()
def downsample(
    input: Path = typer.Argument(..., help="Clonotype sample file."),
    size: int = typer.Argument(..., help="Target size (reads, or unique clonotypes with --clones)."),
    fmt: str = _FMT, out: Optional[Path] = _OUT,
    clones: bool = typer.Option(False, "--clones", help="Draw unique clonotypes, not reads."),
    seed: int = typer.Option(0, "--seed", help="Random seed for reproducible draws."),
) -> None:
    """Randomly down-sample a repertoire to a common depth (reads, or unique clonotypes)."""
    from vdjtools.io.batch import read
    from vdjtools.preprocess import downsample as _ds

    _write(_ds(read(input, fmt=fmt), size, by="clones" if clones else "reads", seed=seed), out)


@app.command(name="filter")
def filter_(
    input: Path = typer.Argument(..., help="Clonotype sample file."),
    fmt: str = _FMT, out: Optional[Path] = _OUT,
    coding: bool = typer.Option(False, "--coding", help="Keep only coding (in-frame, stop-free) clonotypes."),
    noncoding: bool = typer.Option(False, "--noncoding", help="Keep only NON-coding clonotypes (the complement)."),
    min_freq: Optional[float] = typer.Option(None, "--min-freq", help="Keep clonotypes with frequency >= this."),
    v: Optional[str] = typer.Option(None, "--v", help="Comma-separated V segments (prefix ok)."),
    j: Optional[str] = typer.Option(None, "--j", help="Comma-separated J segments (prefix ok)."),
    remove: bool = typer.Option(False, "--remove", help="With --v/--j: remove the listed segments instead of keeping them."),
) -> None:
    """Filter clonotypes: coding / non-coding, by frequency, and/or by V/J segment."""
    from vdjtools import preprocess
    from vdjtools.io.batch import read

    if coding and noncoding:
        _err("--coding and --noncoding are mutually exclusive")
    df = read(input, fmt=fmt)
    if coding or noncoding:
        df = preprocess.filter_functional(df, keep="coding" if coding else "noncoding")
    if min_freq is not None:
        df = preprocess.filter_frequency(df, min_freq=min_freq)
    if v or j:
        df = preprocess.filter_segment(df, v=v.split(",") if v else None,
                                       j=j.split(",") if j else None, keep=not remove)
    _write(df, out)


@app.command()
def pool(
    samples: list[Path] = typer.Argument(..., help="Two or more clonotype sample files."),
    fmt: str = _FMT, out: Optional[Path] = _OUT,
    key: str = typer.Option("aa", "--key", help="Match key: strict | nt | ntV | ntVJ | aa | aaV | aaVJ."),
    join: bool = typer.Option(False, "--join", help="Incidence join (clonotypes shared across samples) instead of a flat pool."),
    min_samples: int = typer.Option(2, "--min-samples", help="With --join: keep clonotypes seen in >= this many samples."),
) -> None:
    """Pool (sum counts) or join (incidence) clonotypes across several samples."""
    from vdjtools import preprocess
    from vdjtools.io.batch import read

    if len(samples) < 2:
        _err("pool needs at least two samples")
    frames = [read(s, fmt=fmt) for s in samples]
    try:
        res = (preprocess.join_samples(frames, key=key, min_samples=min_samples) if join
               else preprocess.pool_samples(frames, key=key))
    except ValueError as e:                          # unknown key
        _err(str(e))
    _write(res, out)


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
def signature(
    samples: Optional[list[Path]] = _SAMPLES, metadata: Optional[Path] = _META,
    base_dir: Optional[Path] = _BASE, sample_col: str = _SCOL, file_template: str = _TMPL,
    fmt: str = _FMT,
    tier: str = typer.Option("standard", help="core | standard | full — nested column sets."),
    weight: str = typer.Option("log2p1", help="Clone-size weight g: log2p1 | duplicate_count | "
                                              "distinct | log1p | anscombe."),
    preset: Optional[str] = typer.Option(None, "--preset",
                                         help="Named feature set (see `vdjtools presets`). "
                                              "Overrides --tier and selects the columns."),
    describe: bool = typer.Option(False, "--describe",
                                  help="Print the column dictionary for --tier and exit."),
    threads: int = _THREADS, out: Optional[Path] = _OUT,
) -> None:
    """One repertoire in, one row of named features out — ready for a classifier.

    Emits the `vsig` (statistics) half of the portable repertoire
    signature: a fixed, named, positional feature vector, so your matrix
    and a collaborator's are the same coordinate system. Reads AIRR
    Rearrangement, native vdjtools, Parquet and the usual third-party
    exports, auto-detected. Writes TSV, or Parquet if -o ends in .parquet.

    \b
    START HERE
      # a metadata sheet plus a directory of samples
      vdjtools signature --preset classify \\
          -m metadata.txt --base-dir samples/ -o sig.tsv
    \b
      # or just pass files
      vdjtools signature --preset compact a.tsv b.tsv.gz -o sig.tsv
    \b
      # the exact columns you will get, reading no input at all
      vdjtools signature --preset classify --describe

    \b
    PICK A PRESET rather than columns by hand (`vdjtools presets` lists all):
      compact    smallest vector that still describes a repertoire (n >= 50)
      classify   general-purpose; the usual random-forest / boosting input
      transfer   for a model that must work on ANOTHER LAB's samples
    \b
    --preset overrides --tier; with neither you get all of --tier (standard).

    THE OTHER HALF: this command is statistics only. The geometry half
    (`rsig`) needs the prototype embedding and ships in mirpy --
    `mir signature --preset classify ...` emits both halves as one vector,
    which is what you usually want for a classifier. A preset spanning both
    halves keeps only its `vsig:` columns here and says so on stderr,
    because silently returning half of what was asked for is worse than
    saying it.

    \b
    GOTCHAS
      * CDR3 vs junction. The reader prefers AIRR `junction_aa` (anchors
        INCLUDED) and falls back to IMGT `cdr3_aa` (anchors excluded), so a
        file carrying only `cdr3_aa` is two residues short everywhere --
        shifting length, k-mer and Pgen features. Check your headers first.
      * Do not PCA-project the result. Plain scaling beat projection at
        every rank tested.
      * -t/--threads defaults to all cores. Inside your own process pool,
        pass -t 1 per worker.
    """
    from vdjtools.signature import layout as L
    from vdjtools.signature import presets as P
    from vdjtools.signature import vsig

    keep = None
    if preset is not None:
        try:
            spec = P.get(preset)
        except KeyError as e:
            _err(str(e))
        tier = spec.tier
        keep = [c for c in spec.columns() if c.startswith("vsig:")]
        if not keep:
            _err(f"preset {preset!r} selects no vsig columns (it is {'+'.join(spec.sig)}); "
                 f"use `mir signature --preset {preset}` for the geometry half")
        dropped = spec.n_columns - len(keep)
        if dropped:
            typer.echo(f"preset {preset!r}: {len(keep)} vsig columns "
                       f"({dropped} rsig columns need `mir signature`)", err=True)

    if tier not in L.TIERS:
        _err(f"--tier must be one of {L.TIERS}; got {tier!r}")
    if describe:
        # The column dictionary for what will actually be emitted, preset or tier.
        d = L.describe(tier)
        _write(d.filter(pl.col('column').is_in(keep)) if keep else d, out)
        return

    from vdjtools.io.batch import map_samples

    items = _sample_items(samples, metadata, base_dir, sample_col, file_template)
    fn = functools.partial(vsig, tier=tier, weight=weight, threads=1)
    rows = [{"sample_id": sid, **res} for sid, res in
            map_samples(fn, items, fmt=fmt, workers=threads or None)]
    if not rows:
        _err("no samples produced a signature")
    cols = keep if keep is not None else L.columns(tier, "vsig")
    _write(pl.DataFrame(rows).select(["sample_id", *cols]), out)


@app.command()
def presets(
    name: Optional[str] = typer.Argument(None, help="Show one preset in full."),
    out: Optional[Path] = _OUT,
) -> None:
    """List the named feature sets for `signature`, with their rankings.

    \b
      vdjtools presets            # the table: name, rank, width, halves, scaling, summary
      vdjtools presets classify   # one preset in full — what is in it, how, and when to use it

    \b
    Ranks tell you how much to trust a choice:
      recommended  use one of these unless you have a reason not to
      specific     correct for a stated purpose and wrong outside it
      avoid        a control or a measured dead end, named so that picking it is deliberate

    The `halves` column says whether a preset needs `vsig` (this package), `rsig` (the geometry
    half, in mirpy), or both. `vdjtools signature` emits only the `vsig:` columns; use
    `mir signature` for a preset spanning both.
    """
    from vdjtools.signature import presets as P

    if name is None:
        _write(P.table().select("preset", "rank", "columns", "halves", "scaling", "summary"), out)
        return
    try:
        spec = P.get(name)
    except KeyError as e:
        _err(str(e))
    typer.echo(f"{spec.name}  [{spec.rank}]  {spec.n_columns} columns  "
               f"tier={spec.tier}  halves={'+'.join(spec.sig)}  scaling={spec.scaling}\n")
    for label, text in (("summary", spec.summary), ("features", spec.features),
                        ("how it is computed", spec.how), ("use cases", spec.use_cases),
                        ("notes", spec.notes)):
        if text:
            typer.echo(f"{label}:\n  {text}\n")


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


# ---------------------------------------------------------------------------- model workshop
model_app = typer.Typer(
    no_args_is_help=True,
    help="Recombination-model workshop: build, learn, check, compare, score, extend.\n\n"
         "A model is named either as a directory, or as LOCUS[:source[:organism]] for a built-in "
         "(TRB, TRB:learned, TRA:arda:mouse). The single-model commands also accept the "
         "-m/--model + --source + --model-path flags the top-level pgen/generate use.",
)
app.add_typer(model_app, name="model")

_MODEL_OUT = typer.Option(None, "--out", "-o", help="Output model directory.")


def _model_arg(spec: str):
    """Resolve ``LOCUS[:source[:organism]]`` or a model directory into a Model."""
    from vdjtools.model import load_bundled, load_model

    path = Path(spec)
    if path.exists() and path.is_dir():
        return load_model(path)
    parts = spec.split(":")
    locus = parts[0]
    source = parts[1] if len(parts) > 1 else "olga"
    organism = parts[2] if len(parts) > 2 else "human"
    try:
        return load_bundled(locus, source=source, organism=organism)
    except (ValueError, FileNotFoundError) as e:
        _err(f"{spec!r} is neither a model directory nor a built-in ({e})")


def _germline_arg(v, j, d, anchors, locus, organism):
    """A germline frame from user FASTA if given, else arda's for the locus/organism."""
    from vdjtools.model import load_germline, read_germline_fasta

    if v and j:
        return read_germline_fasta(v, j, d, anchors=anchors)
    if v or j:
        _err("--germline-v and --germline-j must be given together")
    if not locus:
        _err("give --germline-v/--germline-j, or --locus to take the germline from arda")
    try:
        return load_germline(locus, organism)
    except (ImportError, ValueError) as e:
        _err(str(e))


@model_app.command("list")
def model_list() -> None:
    """List the recombination models shipped with the package (same as ``vdjtools models``)."""
    models()


@model_app.command("check")
def model_check(
    spec: str = typer.Argument(..., help="Model: a directory, or LOCUS[:source[:organism]]."),
    germline: str = typer.Option("auto", help="Reference germline: auto | none | a FASTA path."),
    organism: str = typer.Option("human", help="Organism, when --germline is auto."),
    out: Optional[Path] = _OUT,
) -> None:
    """Audit a model against its manifest, its germline, and a reference library.

    Writes a tidy issue frame (severity, check, event, segment, allele, detail, value) and
    **exits 1 if any issue has severity "error"**, so it works as a gate in a build script.
    """
    from vdjtools.model.check import check_model

    m = _model_arg(spec)
    gl = germline
    if germline not in ("auto", "none"):
        from vdjtools.model import read_germline_fasta

        gl = read_germline_fasta(germline, germline)
    issues = check_model(m, germline=gl)
    _write(issues, out)
    n_err = issues.filter(pl.col("severity") == "error").height
    n_warn = issues.filter(pl.col("severity") == "warn").height
    _info(f"{n_err} error(s), {n_warn} warning(s)")
    if n_err:
        raise typer.Exit(1)


@model_app.command("template")
def model_template(
    locus: Optional[str] = typer.Option(None, "--locus", help="Locus, e.g. TRB."),
    organism: str = typer.Option("human", help="Organism (arda germline), when no FASTA is given."),
    germline_v: Optional[Path] = typer.Option(None, "--germline-v", help="V-allele FASTA."),
    germline_j: Optional[Path] = typer.Option(None, "--germline-j", help="J-allele FASTA."),
    germline_d: Optional[Path] = typer.Option(None, "--germline-d", help="D-allele FASTA (makes it VDJ)."),
    anchors: Optional[Path] = typer.Option(None, help="CDR3-anchor CSV, if the FASTAs are full-length."),
    ins_max: int = typer.Option(40, help="Largest N-region insertion in the placeholder tables."),
    out: Optional[Path] = _MODEL_OUT,
) -> None:
    """Build a model scaffold from a germline library — your own FASTA, or arda's.

    The marginals are placeholders meant to be refit with ``model learn``; their support ranges
    bound what EM can then learn, which is why ``--ins-max`` is here.
    """
    from vdjtools.model.io import from_germline

    if out is None:
        _err("give an output model directory with -o")
    gl = _germline_arg(germline_v, germline_j, germline_d, anchors, locus, organism)
    try:
        m = from_germline(gl, locus=locus or "CUSTOM", organism=organism, ins_max=ins_max,
                          strict=germline_v is not None)
    except ValueError as e:
        _err(str(e))
    m.save(out)
    _info(f"{m.chain_type} template for {m.locus}: "
          f"{m.genomic['genes_v'].height} V, {m.genomic['genes_j'].height} J -> {out}")


@model_app.command("learn")
def model_learn(
    input: Path = typer.Argument(..., help="Clonotype table (TSV/Parquet) with junction + V/J calls."),
    template: Optional[str] = typer.Option(None, "--template", "-t",
                                           help="Template model, or LOCUS[:source[:organism]]."),
    locus: Optional[str] = typer.Option(None, "--locus", help="Build an arda template for this locus."),
    organism: str = typer.Option("human", help="Organism, with --locus."),
    column: Optional[str] = typer.Option(None, "--column", "-c", help="Junction column (auto-detected)."),
    max_iter: int = typer.Option(15, help="EM iteration cap."),
    tol: float = typer.Option(1e-4, help="Stop below this relative log-likelihood improvement."),
    init: str = typer.Option("align", help="align | uniform | template (template = warm start / fine-tune)."),
    gene_prior: float = typer.Option(1.0, help="Dirichlet pseudocount over functional V/J alleles."),
    nd_prior: float = typer.Option(0.0, help="Pseudocount pushing P(n_D=2) toward 0."),
    single_d: bool = typer.Option(False, "--single-d", help="Force a strict single-D model."),
    no_calls: bool = typer.Option(False, "--no-calls", help="Ignore V/J calls (much slower)."),
    verbose: bool = typer.Option(False, "--verbose", "-v",
                                 help="Report each EM iteration as it happens (to stderr)."),
    checkpoint: Optional[Path] = typer.Option(
        None, "--checkpoint", help="Save the model after each iteration, so a long fit survives "
                                   "being interrupted (resume with --resume)."),
    checkpoint_every: int = typer.Option(1, help="Checkpoint every N iterations."),
    resume_from: Optional[Path] = typer.Option(
        None, "--resume", help="Continue a fit from a checkpoint instead of starting over."),
    out: Optional[Path] = _MODEL_OUT,
) -> None:
    """Fit a model's marginals from your own sequences by EM, writing the training log alongside.

    ``--init template`` warm-starts from the template instead of realigning, which is how you
    fine-tune an existing model on a new sample rather than fitting from scratch. ``-v`` prints the
    log-likelihood and its relative change per iteration, so a long fit is visibly converging
    rather than merely running.

    For a fit that will not finish in one sitting, ``--checkpoint DIR`` saves the model after every
    iteration and ``--resume DIR`` picks it back up — resuming reaches the same log-likelihood as an
    uninterrupted run, and the training log spans every attempt.
    """
    from vdjtools.io.batch import read as _read
    from vdjtools.model.infer import infer_frame, print_progress, training_frame

    if out is None:
        _err("give an output model directory with -o")
    if resume_from is not None:
        base = _model_arg(str(resume_from))
        init = "template"
        if checkpoint is None:
            checkpoint = resume_from
        _info(f"resuming from {resume_from}")
    elif template is None and locus is None:
        _err("give a template with --template, a locus with --locus, or --resume a checkpoint")
    elif template is not None:
        base = _model_arg(template)
    else:
        from vdjtools.model.io import from_arda

        base = from_arda(locus, organism)
    try:
        clones = _read(input, fmt="auto")
    except Exception:  # noqa: BLE001 - fall back to a plain table read
        clones = pl.read_parquet(input) if input.suffix in (".parquet", ".pq") else \
            pl.read_csv(input, separator="\t", infer_schema_length=20000)
    try:
        m, rep = infer_frame(base, clones, seq_col=column, use_calls=not no_calls,
                             max_iter=max_iter, tol=tol, init=init,
                             gene_prior=gene_prior, nd_prior=nd_prior, single_d=single_d,
                             progress=print_progress() if verbose else None,
                             checkpoint=checkpoint, checkpoint_every=checkpoint_every)
    except (ValueError, KeyError) as e:
        _err(str(e))
    m.save(out)
    log = training_frame(m)
    _info(f"{rep.n_iter} iterations, converged={rep.converged}, "
          f"loglik {log['loglik'][0]:.3f} -> {log['loglik'][-1]:.3f} -> {out}")


@model_app.command("build")
def model_build(
    chains: str = typer.Option(",".join(("TRA", "TRB")), "--chains", help="Comma-separated chains."),
    groups: str = typer.Option("human", "--groups", help="Comma-separated read groups."),
    workers: Optional[int] = typer.Option(None, help="Concurrent chain builds (default: cores/2)."),
    work_dir: Path = typer.Option(Path("/tmp/vdjtools_build"), help="Scratch dir for arda output."),
    cap: Optional[int] = typer.Option(None, help="Cap reads per chain (default: use every read)."),
    max_iter: int = typer.Option(15, help="EM iteration cap."),
    verbose: bool = typer.Option(False, "--verbose", "-v",
                                 help="Stream arda's mapping output and each EM iteration."),
    out: Optional[Path] = _MODEL_OUT,
) -> None:
    """Build models from the full AIRR read corpus: fetch, arda-map, then EM — several chains at once.

    This is the real training path (raw FASTQ from the ``isalgo/airr_model_read`` dataset), so it
    needs HuggingFace access and arda's mmseqs2. Mapping is minutes per chain and EM on a D-bearing
    locus can be far longer, so **use ``-v``** — without it the whole run is silent until a chain
    finishes, and a slow fit is indistinguishable from a stuck one.
    """
    from vdjtools.model.data import build_all

    if out is None:
        _err("give an output directory with -o")
    res = build_all([c.strip() for c in chains.split(",") if c.strip()],
                    groups=tuple(g.strip() for g in groups.split(",") if g.strip()),
                    workers=workers, out_dir=out, work_dir=work_dir, cap=cap, iters=max_iter,
                    verbose=verbose)
    ok = [r["stats"] for r in res.values() if "stats" in r]
    for key, r in res.items():
        if "error" in r:
            typer.secho(f"{key}: FAILED {r['error']}", fg=typer.colors.RED, err=True)
    if not ok:
        _err("every build failed")
    _write(pl.DataFrame(ok), None)
    _info(f"{len(ok)}/{len(res)} chain(s) built -> {out}")


@model_app.command("extend")
def model_extend(
    spec: str = typer.Argument(..., help="Model: a directory, or LOCUS[:source[:organism]]."),
    locus: Optional[str] = typer.Option(None, "--locus", help="Locus for the arda germline (default: the model's)."),
    organism: str = typer.Option("human", help="Organism for the arda germline."),
    germline_v: Optional[Path] = typer.Option(None, "--germline-v", help="V-allele FASTA."),
    germline_j: Optional[Path] = typer.Option(None, "--germline-j", help="J-allele FASTA."),
    germline_d: Optional[Path] = typer.Option(None, "--germline-d", help="D-allele FASTA."),
    anchors: Optional[Path] = typer.Option(None, help="CDR3-anchor CSV, if the FASTAs are full-length."),
    weight: float = typer.Option(1.0, help="Mass for a new allele of a known gene, relative to its gene-mate."),
    out: Optional[Path] = _MODEL_OUT,
) -> None:
    """Add alleles from a larger germline library, seeded from what the model already knows.

    Each pre-existing gene keeps its total usage — a richer library splits a gene's mass more
    finely rather than multiplying it. This seeds; follow with ``model learn --init template``.
    """
    from vdjtools.model.infer import extend_alleles

    if out is None:
        _err("give an output model directory with -o")
    m = _model_arg(spec)
    gl = _germline_arg(germline_v, germline_j, germline_d, anchors, locus or m.locus, organism)
    before = {k: v.height for k, v in m.genomic.items()}
    try:
        e = extend_alleles(m, gl, weight=weight)
    except ValueError as ex:
        _err(str(ex))
    e.save(out)
    _info(", ".join(f"{k}: {before[k]}->{v.height}" for k, v in e.genomic.items()) + f" -> {out}")


@model_app.command("rescale")
def model_rescale(
    spec: str = typer.Argument(..., help="Model: a directory, or LOCUS[:source[:organism]]."),
    samples: list[Path] = typer.Argument(..., help="Clonotype sample file(s) supplying the target usage."),
    fmt: str = _FMT,
    no_v: bool = typer.Option(False, "--no-v", help="Leave P(V) alone."),
    no_j: bool = typer.Option(False, "--no-j", help="Leave P(J) alone."),
    aggregate: str = typer.Option("pool", help="Combine several samples: pool | mean."),
    out: Optional[Path] = _MODEL_OUT,
) -> None:
    """Replace a model's V/J usage with your own sample's, keeping its junction model.

    V/J usage is protocol-dependent (5'RACE and DNA-multiplex amplify different V genes at very
    different rates); the recombination machinery underneath is not. Pass the repertoire you are
    actually going to score.
    """
    from vdjtools.io.batch import read as _read
    from vdjtools.model import rescale_usage

    if out is None:
        _err("give an output model directory with -o")
    m = _model_arg(spec)
    frames = [_read(s, fmt=fmt) for s in samples]
    try:
        r = rescale_usage(m, frames if len(frames) > 1 else frames[0],
                          v=not no_v, j=not no_j, aggregate=aggregate)
    except ValueError as e:
        _err(str(e))
    r.save(out)
    _info(f"usage rescaled from {len(frames)} sample(s) -> {out}")


@model_app.command("export")
def model_export(
    spec: str = typer.Argument(..., help="Model: a directory, or LOCUS[:source[:organism]]."),
    long: bool = typer.Option(False, "--long", help="One flat long table of every marginal."),
    format: str = typer.Option("tsv", "--format", help="Model-directory format: tsv | csv | parquet."),
    out: Optional[Path] = typer.Option(None, "--out", "-o", help="Output file (--long) or directory."),
) -> None:
    """Export a model's probabilities as tables — a hand-editable directory, or one long frame.

    Both round-trip: a TSV model directory loads straight back with ``--model-path``, and the long
    frame goes back through ``vdjtools.model.set_marginals``.
    """
    from vdjtools.model.io import marginals_frame

    m = _model_arg(spec)
    if long:
        _write(marginals_frame(m), out)
        return
    if out is None:
        _err("give an output directory with -o (or use --long to write one table)")
    try:
        m.save(out, fmt=format)
    except ValueError as e:
        _err(str(e))
    _info(f"{len(m.tables)} marginal + {len(m.genomic)} germline table(s) as {format} -> {out}")


@model_app.command("net")
def model_net(
    spec: str = typer.Argument(..., help="Model: a directory, or LOCUS[:source[:organism]]."),
    format: str = typer.Option("dot", "--format", help="dot | svg | pdf | png (non-dot needs graphviz)."),
    out: Optional[Path] = typer.Option(None, "--out", "-o", help="Output file (default: stdout, dot only)."),
) -> None:
    """Render the model's recombination Bayes net, nodes annotated with entropy and edges with MI."""
    from vdjtools.model.analyze import bayes_net_dot, render_dot

    m = _model_arg(spec)
    dot = bayes_net_dot(m)
    if format == "dot":
        if out is None:
            typer.echo(dot)
        else:
            Path(out).write_text(dot)
            _info(f"wrote {out}")
        return
    if out is None:
        _err("give an output path with -o for a rendered format")
    try:
        _info(f"wrote {render_dot(dot, out, fmt=format)}")
    except RuntimeError as e:
        _err(str(e))


@model_app.command("entropy")
def model_entropy(
    spec: str = typer.Argument(..., help="Model: a directory, or LOCUS[:source[:organism]]."),
    table: str = typer.Option("entropy", help="entropy | mi | total."),
    out: Optional[Path] = _OUT,
) -> None:
    """Information content per recombination event: entropy, mutual information, or the total."""
    from vdjtools.model.analyze import entropy_table, mutual_information, total_entropy

    m = _model_arg(spec)
    if table == "entropy":
        _write(entropy_table(m), out)
    elif table == "mi":
        _write(mutual_information(m), out)
    elif table == "total":
        t = total_entropy(m)
        _write(t, out)
        _info(f"scenario entropy {t['contribution_bits'].sum():.2f} bits "
              f"({2 ** t['contribution_bits'].sum():.3g} distinct rearrangements)")
    else:
        _err("--table must be entropy, mi or total")


@model_app.command("diversity")
def model_diversity(
    spec: str = typer.Argument(..., help="Model: a directory, or LOCUS[:source[:organism]]."),
    n: int = typer.Option(5000, "--number", "-n", help="Sequences to generate for the estimate."),
    seed: int = typer.Option(0, help="Generation seed."),
    productive: bool = typer.Option(False, "--productive", help="Estimate over productive rearrangements only."),
    out: Optional[Path] = _OUT,
) -> None:
    """Estimate total diversity: scenario entropy, sequence entropy, and effective diversity.

    Reports both Hill numbers — ``2^H`` (the usual "~10^x distinct sequences" figure) and
    ``1/E[Pgen]`` (how many draws before two coincide). Monte Carlo, so give it a seed.
    """
    from vdjtools.model.score import diversity as _div

    m = _model_arg(spec)
    try:
        d = _div(m, n=n, seed=seed, productive_only=productive)
    except ValueError as e:
        _err(str(e))
    _write(d, out)
    r = d.to_dicts()[0]
    _info(f"scenario {r['scenario_entropy_bits']:.1f} bits, sequence "
          f"{r['sequence_entropy_bits']:.2f}+-{r['sequence_entropy_se_bits']:.2f} bits -> "
          f"Shannon {r['diversity_shannon']:.3g}, Simpson {r['diversity_simpson']:.3g}")


@model_app.command("compare")
def model_compare(
    a: str = typer.Argument(..., help="First model: a directory, or LOCUS[:source[:organism]]."),
    b: str = typer.Argument(..., help="Second model."),
    by: str = typer.Option("allele", help="Align on allele | gene. Use gene across germline sources."),
    usage: Optional[str] = typer.Option(None, "--usage", help="Instead emit V/J usage side by side: v | j | d."),
    dot: Optional[Path] = typer.Option(None, "--dot", help="Also write the comparison graph here."),
    dot_format: str = typer.Option("dot", "--dot-format", help="dot | svg | pdf | png."),
    out: Optional[Path] = _OUT,
) -> None:
    """Compare two models parameter by parameter: per-event divergence and support differences.

    Jensen-Shannon is the headline (symmetric, bounded, finite when the supports differ);
    ``tv_max`` finds the one broken gene an average hides.
    """
    from vdjtools.model.analyze import compare_models, compare_net_dot, compare_usage, render_dot

    ma, mb = _model_arg(a), _model_arg(b)
    labels = (a, b)
    if usage:
        _write(compare_usage(ma, mb, usage), out)
    else:
        if by not in ("allele", "gene"):
            _err("--by must be allele or gene")
        _write(compare_models(ma, mb, labels=labels, by=by), out)
    if dot is not None:
        src = compare_net_dot(ma, mb, labels=labels)
        if dot_format == "dot":
            Path(dot).write_text(src)
            _info(f"wrote {dot}")
        else:
            try:
                _info(f"wrote {render_dot(src, dot, fmt=dot_format)}")
            except RuntimeError as e:
                _err(str(e))


@model_app.command("compare-pgen")
def model_compare_pgen(
    a: str = typer.Argument(..., help="First model: a directory, or LOCUS[:source[:organism]]."),
    b: str = typer.Argument(..., help="Second model."),
    input: Path = typer.Argument(..., help="Table (TSV) or list of CDR3 sequences."),
    column: Optional[str] = typer.Option(None, "--column", "-c", help="Sequence column (auto-detected)."),
    v_col: Optional[str] = typer.Option(None, "--v-col", help="V-call column to condition on."),
    j_col: Optional[str] = typer.Option(None, "--j-col", help="J-call column to condition on."),
    seq_type: str = typer.Option("auto", "--type", help="auto | aa | nt."),
    summary: bool = typer.Option(False, "--summary", help="One row of distribution statistics instead."),
    no_header: bool = typer.Option(False, "--no-header", help="Input is a bare sequence list."),
    out: Optional[Path] = _OUT,
) -> None:
    """Score one sequence set under two models and compare the Pgen distributions.

    ``--summary`` gives correlations, the KS statistic, and — the number that usually matters —
    how many sequences each model can score that the other assigns Pgen 0.
    """
    from vdjtools.model.score import compare_pgen, pgen_summary

    ma, mb = _model_arg(a), _model_arg(b)
    df, seqcol = _read_seq_table(input, column, no_header)
    use_calls = bool(v_col or j_col)
    try:
        cmp = compare_pgen(ma, mb, df, labels=("a", "b"), kind=seq_type, use_calls=use_calls,
                           on_unknown="marginalize", seq_col=seqcol, v_col=v_col, j_col=j_col)
    except (ValueError, KeyError) as e:
        _err(str(e))
    _write(pgen_summary(cmp) if summary else cmp, out)


@model_app.command("loglik")
def model_loglik(
    input: Path = typer.Argument(..., help="Table (TSV) or list of CDR3 sequences."),
    spec: str = typer.Argument(..., help="Model: a directory, or LOCUS[:source[:organism]]."),
    column: Optional[str] = typer.Option(None, "--column", "-c", help="Sequence column (auto-detected)."),
    v_col: Optional[str] = typer.Option(None, "--v-col", help="V-call column to condition on."),
    j_col: Optional[str] = typer.Option(None, "--j-col", help="J-call column to condition on."),
    weights_col: Optional[str] = typer.Option(None, "--weights-col", help="Per-clonotype weight column."),
    seq_type: str = typer.Option("auto", "--type", help="auto | aa | nt."),
    per_sequence: bool = typer.Option(False, "--per-sequence", help="Emit per-sequence Pgen instead."),
    no_header: bool = typer.Option(False, "--no-header", help="Input is a bare sequence list."),
    out: Optional[Path] = _OUT,
) -> None:
    """How well a model explains a sequence set: log-likelihood, free parameters, AIC and BIC.

    Nucleotide input gives a properly normalized likelihood, so BIC is meaningful; amino-acid input
    is a relative score on one fixed sequence set only. Sequences the model cannot generate are
    counted in ``n_scoreable``, never turned into -inf.
    """
    from vdjtools.model.score import model_fit, pgen_frame

    m = _model_arg(spec)
    df, seqcol = _read_seq_table(input, column, no_header)
    use_calls = bool(v_col or j_col)
    kw = dict(kind=seq_type, use_calls=use_calls, on_unknown="marginalize",
              seq_col=seqcol, v_col=v_col, j_col=j_col)
    try:
        if per_sequence:
            _write(pgen_frame(m, df, **kw), out)
            return
        fit = model_fit(m, df, weights=weights_col, **kw)
    except (ValueError, KeyError) as e:
        _err(str(e))
    _write(fit, out)
    r = fit.to_dicts()[0]
    _info(f"loglik {r['loglik_sum']:.1f} over {r['n_scoreable']:.0f}/{r['n']:.0f} scoreable, "
          f"k={r['k']}, AIC={r['aic']:.1f}, BIC={r['bic']:.1f}")


@model_app.command("log")
def model_log(
    spec: str = typer.Argument(..., help="Model: a directory, or LOCUS[:source[:organism]]."),
    out: Optional[Path] = _OUT,
) -> None:
    """Show a model's EM training log — log-likelihood per iteration, one block per run."""
    from vdjtools.model.infer import training_frame

    m = _model_arg(spec)
    log = training_frame(m)
    if not log.height:
        _err(f"{spec!r} carries no training log (it was not fitted by this tool)")
    _write(log, out)
