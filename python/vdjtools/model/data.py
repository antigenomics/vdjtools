"""Bootstrap the recombination model from real AIRR reads (private HF ``isalgo/airr_model_read``).

The dataset ships raw 5'-RACE FASTQ reads per organism-group (``human``, ``human_fetal`` — TdT-low,
``mouse``) and chain, split into **functional** (productive) and **non-functional** (out-of-frame or
stop) buckets. Training / benchmarking the engine means:

1. :func:`fetch_fastq` — pull one ``{group}/{CHAIN}.{label}.fq.gz`` from the hub.
2. :func:`annotate_reads` — map the reads with **arda** (``arda rnaseq map``) to a per-read AIRR
   table with V/D/J (and **D2** for D-D joins) calls, junction, CIGARs, isotype, and productivity.
3. :func:`unique_clonotypes` — collapse to the model's clonotype identity: same V allele, J allele,
   and junction nucleotides. Reads differing only in alignment (CIGAR) or isotype (IGH ``c_call``)
   are one clonotype; isotype is dropped from the key (isotype switching is the same clonotype) and
   the collapse can be restricted to naive IgM.

The non-functional clonotypes are the unbiased EM training set; the functional ones are the
selection-shaped test set. ``huggingface_hub`` (fetch) and ``arda[rnaseq]`` (annotate; needs arda +
mmseqs2 + seqtree) are lazy, tool-only imports — never runtime dependencies of the model math.

Example::

    from vdjtools.model import from_olga, data
    from vdjtools.model.infer import infer_native
    clones = data.prepare("human", "TRB", "nonfunctional", out_dir="/tmp/arda", cap=200_000)
    template = from_olga(olga_dir, locus="TRB")
    fit, rep = infer_native(template, clones["junction"].to_list())
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import polars as pl


def _arda_exe() -> str:
    """The arda CLI next to the running interpreter (the venv's, with the rnaseq extra), or PATH."""
    cand = Path(sys.executable).with_name("arda")
    return str(cand) if cand.exists() else "arda"

#: Private HuggingFace dataset of real AIRR reads (owner isalgo, cc-by-nc-nd); fetched, not vendored.
MODEL_READS_REPO = "isalgo/airr_model_read"

GROUPS = ("human", "human_fetal", "mouse")
CHAINS = ("IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD")
LABELS = ("functional", "nonfunctional")
#: arda reference organism per group (fetal T cells use the human reference).
ORGANISM = {"human": "human", "human_fetal": "human", "mouse": "mouse"}


def fetch_fastq(group: str, chain: str, label: str, *, repo: str = MODEL_READS_REPO) -> str:
    """Download one ``{group}/{CHAIN}.{label}.fq.gz`` from the dataset; return its local path."""
    from huggingface_hub import hf_hub_download

    return hf_hub_download(repo_id=repo, filename=f"{group}/{chain}.{label}.fq.gz", repo_type="dataset")


def _subsample_fastq(src: str | Path, dst: str | Path, n: int) -> None:
    """Write the first ``n`` FASTQ records of ``src`` (optionally gzipped) to ``dst``."""
    import gzip

    opener = gzip.open if str(src).endswith(".gz") else open
    with opener(src, "rt") as fh, open(dst, "wt") as out:
        for i, line in enumerate(fh):
            if i >= 4 * n:
                break
            out.write(line)


def annotate_reads(
    fq_path: str | Path,
    *,
    out_dir: str | Path,
    prefix: str,
    organism: str = "human",
    cap: int | None = None,
    threads: int = 0,
    reconstruct: bool = False,
) -> pl.DataFrame:
    """Map a FASTQ to a per-read AIRR table with arda (``arda rnaseq map`` → ``<prefix>.airr.tsv``).

    Returns the raw per-read annotation (V/D/J and **D2** calls, junction, CIGARs, isotype,
    productivity); collapse it to clonotypes with :func:`unique_clonotypes`. We deliberately stop at
    ``map`` — the clonotype identity "same V/J allele + junction, up to CIGAR" is exactly that
    dedup, so arda's heavier error-model ``correct`` stage is not needed (and, on these single-end
    5'-RACE reads, it discards reads that map but whose mate-spanned junction it cannot reassemble).

    Args:
        fq_path: Input FASTQ (``.fq``/``.fq.gz``), e.g. from :func:`fetch_fastq`.
        out_dir: Directory for arda's output.
        prefix: Output basename.
        organism: arda reference organism (``"human"`` / ``"mouse"``).
        cap: If set, annotate only the first ``cap`` reads (bounded benchmark scale).
        threads: mmseqs threads (0 = all cores).
        reconstruct: Merge overlapping paired mates into one fragment before mapping (recovers
            longer junctions single reads don't span; needs the FASTQ to carry both mates).

    Returns:
        The per-read AIRR :class:`polars.DataFrame`. Reads that do not span a full junction have an
        empty ``junction`` and are dropped by :func:`unique_clonotypes`.

    Requires the ``arda`` CLI with the ``rnaseq`` extra (arda + mmseqs2) on ``PATH``.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    r1 = Path(fq_path)
    if cap is not None:
        r1 = out_dir / f"{prefix}.sub.fq"
        _subsample_fastq(fq_path, r1, cap)
    airr = out_dir / f"{prefix}.airr.tsv"
    cmd = [
        _arda_exe(), "rnaseq", "map", "-o", str(airr), "--r1", str(r1),
        "--organism", organism, "--threads", str(threads),
    ]
    if reconstruct:
        cmd.append("--reconstruct")
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    return pl.read_csv(airr, separator="\t", infer_schema_length=20000)


def unique_clonotypes(clones: pl.DataFrame, *, naive_igm_only: bool = False) -> pl.DataFrame:
    """Collapse an arda clonotype table to unique clonotypes = ``(v_call, j_call, junction)``.

    Isotype (``c_call``) is dropped from the key — isotype switching is the same clonotype. Rows are
    keyed to allele resolution (V/J allele) and full junction nt, so reads that differ only in
    alignment (CIGAR) or isotype collapse together; read support sums into ``count``.

    Args:
        clones: An arda clonotype frame (from :func:`annotate_reads`).
        naive_igm_only: Keep only IgM (naive B) clonotypes before collapsing (IGH only; no-op
            elsewhere) — for a naive-repertoire (pre-selection-independent) subset.

    Returns:
        Deduplicated clonotype frame: ``v_call, j_call, junction, junction_aa, locus, d_call,
        d2_call, count`` (one row per unique clonotype), sorted by descending ``count``.
    """
    df = clones
    if naive_igm_only and "c_call" in df.columns:
        df = df.filter(pl.col("c_call").fill_null("").str.starts_with("IGHM"))
    df = df.filter(pl.col("junction").is_not_null() & (pl.col("junction").str.len_bytes() > 0))
    # arda writes an empty string (not null) when a D / second D is absent; normalize so downstream
    # ``is_not_null`` / D-count logic is correct. duplicate_count can arrive typed as str.
    empty_to_null = [
        pl.when(pl.col(c).cast(pl.Utf8).str.len_bytes() > 0).then(pl.col(c)).otherwise(None).alias(c)
        for c in ("d_call", "d2_call") if c in df.columns
    ]
    if empty_to_null:
        df = df.with_columns(empty_to_null)
    cnt = pl.col("duplicate_count").cast(pl.Int64, strict=False) if "duplicate_count" in df.columns else pl.lit(1)
    out = df.group_by(["v_call", "j_call", "junction"]).agg(
        pl.col("junction_aa").first(),
        pl.col("locus").first(),
        pl.col("d_call").first(),
        pl.col("d2_call").first(),
        cnt.sum().alias("count"),
    )
    return out.sort("count", descending=True)


def prepare(
    group: str,
    chain: str,
    label: str,
    *,
    out_dir: str | Path,
    cap: int | None = None,
    reconstruct: bool = False,
    naive_igm_only: bool = False,
) -> pl.DataFrame:
    """Fetch → map (arda) → unique clonotypes for one ``(group, chain, label)`` bucket."""
    fq = fetch_fastq(group, chain, label)
    reads = annotate_reads(
        fq, out_dir=out_dir, prefix=f"{group}_{chain}_{label}",
        organism=ORGANISM[group], cap=cap, reconstruct=reconstruct,
    )
    return unique_clonotypes(reads, naive_igm_only=naive_igm_only)
