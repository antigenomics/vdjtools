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


def _sample_frames(samples, metadata, base_dir, sample_col, file_template, fmt) -> dict:
    """{sample_id: frame} — from a metadata table if given, else the positional sample files."""
    from vdjtools.io.batch import read, read_metadata, read_samples

    if metadata is not None:
        md = read_metadata(metadata)
        if sample_col not in md.columns:
            _err(f"--sample-col {sample_col!r} not in metadata columns: {md.columns}")
        return read_samples(md, base_dir or Path(metadata).parent, sample_col=sample_col,
                            file_template=file_template, fmt=fmt, add_metadata=False, as_dict=True)
    if not samples:
        _err("give sample files as arguments, or -m/--metadata <table> (with --base-dir)")
    return {Path(s).name.split(".")[0]: read(s, fmt=fmt) for s in samples}


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


@app.command()
def diversity(
    samples: Optional[list[Path]] = _SAMPLES, metadata: Optional[Path] = _META,
    base_dir: Optional[Path] = _BASE, sample_col: str = _SCOL, file_template: str = _TMPL,
    fmt: str = _FMT, out: Optional[Path] = _OUT,
) -> None:
    """Per-sample diversity (observed richness, Chao, Efron-Thisted, Shannon, Simpson, d50)."""
    from vdjtools.stats.diversity import diversity_stats

    frames = _sample_frames(samples, metadata, base_dir, sample_col, file_template, fmt)
    rows = [_tag(diversity_stats(df), sid) for sid, df in frames.items()]
    _write(pl.concat(rows, how="vertical_relaxed"), out)


@app.command()
def spectratype(
    samples: Optional[list[Path]] = _SAMPLES, metadata: Optional[Path] = _META,
    base_dir: Optional[Path] = _BASE, sample_col: str = _SCOL, file_template: str = _TMPL,
    fmt: str = _FMT, kind: str = typer.Option("aa", help="aa | nt (length unit)."),
    weight: str = typer.Option("reads", help="reads | unique | freq."), out: Optional[Path] = _OUT,
) -> None:
    """Per-sample CDR3 length spectratype."""
    from vdjtools.stats.spectratype import spectratype as _spec

    frames = _sample_frames(samples, metadata, base_dir, sample_col, file_template, fmt)
    rows = [_tag(_spec(df, kind=kind, weight=weight), sid) for sid, df in frames.items()]
    _write(pl.concat(rows, how="vertical_relaxed"), out)


@app.command(name="segment-usage")
def segment_usage(
    samples: Optional[list[Path]] = _SAMPLES, metadata: Optional[Path] = _META,
    base_dir: Optional[Path] = _BASE, sample_col: str = _SCOL, file_template: str = _TMPL,
    fmt: str = _FMT, segment: str = typer.Option("v", help="v | d | j | c."),
    weight: str = typer.Option("reads", help="reads | unique | freq."), out: Optional[Path] = _OUT,
) -> None:
    """Per-sample V/D/J/C segment usage."""
    from vdjtools.stats.usage import segment_usage as _usage

    frames = _sample_frames(samples, metadata, base_dir, sample_col, file_template, fmt)
    rows = [_tag(_usage(df, segment=segment, weight=weight), sid) for sid, df in frames.items()]
    _write(pl.concat(rows, how="vertical_relaxed"), out)


@app.command()
def overlap(
    samples: Optional[list[Path]] = _SAMPLES, metadata: Optional[Path] = _META,
    base_dir: Optional[Path] = _BASE, sample_col: str = _SCOL, file_template: str = _TMPL,
    fmt: str = _FMT, out: Optional[Path] = _OUT,
) -> None:
    """Exact-match pairwise repertoire overlap (D, F, F2, R) for every sample pair."""
    from vdjtools.overlap.metrics import overlap_metrics

    frames = _sample_frames(samples, metadata, base_dir, sample_col, file_template, fmt)
    items = list(frames.items())
    if len(items) < 2:
        _err("overlap needs at least two samples")
    rows = []
    for i in range(len(items)):
        for k in range(i + 1, len(items)):
            (a_id, a), (b_id, b) = items[i], items[k]
            rows.append({"sample_a": a_id, "sample_b": b_id, **overlap_metrics(a, b)})
    _write(pl.DataFrame(rows), out)
