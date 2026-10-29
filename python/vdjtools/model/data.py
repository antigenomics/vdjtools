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

from vdjtools.cores import kernel_threads, require_serial_pools

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
    """Download one ``{group}/{CHAIN}.{label}.fq.gz`` from the dataset; return its local path.

    Raises:
        PermissionError: If the dataset is not readable with the current credentials. It is
            **private**, so that is the expected outcome for anyone outside the lab, and the bare
            ``huggingface_hub`` 401 does not say which repo it means or what to do instead.
    """
    from huggingface_hub import hf_hub_download

    try:
        return hf_hub_download(repo_id=repo, filename=f"{group}/{chain}.{label}.fq.gz",
                               repo_type="dataset")
    except Exception as e:                        # 401/403 from a private repo, or no token at all
        if "401" not in str(e) and "403" not in str(e) and "gated" not in str(e).lower():
            raise
        raise PermissionError(
            f"{repo!r} is a PRIVATE HuggingFace dataset and your credentials cannot read it. "
            "This is only needed to RETRAIN the recombination models from raw reads -- the "
            "bundled models ship with vdjtools, are the versioned reference, and need no rebuild "
            "(see `vdjtools models`). This is expected outside the lab and is not a bug."
        ) from e


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
    verbose: bool = False,
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
    # `arda map` is stage 1 on its own, and is what we want -- see the note above on stopping
    # before `correct`. This used to be spelled `arda rnaseq map -o ...`; arda 2.19 turned
    # `rnaseq` into the full map->assemble->correct preset that takes NO stage positional and
    # writes via `-p/--out-prefix` + `-d/--out-dir` rather than `-o`, so the old argv exits 2.
    # (`arda rnaseq map --help` still exits 0 -- --help short-circuits before argument parsing --
    # so a --help smoke test would not have caught this.)
    cmd = [
        _arda_exe(), "map", "-o", str(airr), "--r1", str(r1),
        "--organism", organism, "--threads", str(kernel_threads(threads)),
    ]
    if reconstruct:
        cmd.append("--reconstruct")
    # verbose lets arda's own mapping progress through: it is the slowest stage (minutes on a full
    # chain) and, captured, it looks identical to a hang.
    proc = subprocess.run(cmd, capture_output=not verbose, text=True)
    if proc.returncode:
        # capture_output swallows arda's own message, so a failed run used to surface as a bare
        # CalledProcessError with an exit code and nothing else -- which is exactly how the
        # `arda rnaseq map` -> `arda map` rename stayed invisible. Say what arda said, and name
        # the version, since a usage error here is almost always a CLI-version mismatch.
        try:
            from importlib.metadata import version

            installed = version("arda-mapper")
        except Exception:  # noqa: BLE001
            installed = "unknown"
        detail = (proc.stderr or proc.stdout or "").strip()
        raise RuntimeError(
            f"arda exited {proc.returncode} (arda-mapper {installed}; vdjtools needs >= 2.19, "
            f"where stage-1 mapping is `arda map`):\n"
            f"  {' '.join(cmd)}\n{detail}"
        )
    return pl.read_csv(airr, separator="\t", infer_schema_length=20000)


def _with_candidate_sequences(clones: pl.DataFrame) -> pl.DataFrame:
    """Retain coding-strand full reads as candidate proposals, separate from junctions."""
    from ..io.schema import normalize_nt_frame

    clones = normalize_nt_frame(clones, ["sequence", "junction"])
    if "candidate_sequences" in clones.columns:
        candidates = pl.col("candidate_sequences").cast(pl.List(pl.String))
        bases = pl.element().str.to_uppercase()
        invalid = candidates.list.eval(
            bases.str.contains(r"[^ACGTRYSWKMBDHVN]").fill_null(False)).list.any().fill_null(False)
        if clones.select(invalid.any()).item():
            raise ValueError("invalid DNA IUPAC nucleotide input")
        clones = clones.with_columns(candidates.list.eval(
            bases.str.replace_all(r"[RYSWKMBDHVN]", "A")))
    if "sequence" not in clones.columns:
        return clones
    seq = pl.col("sequence").cast(pl.String)
    if "rev_comp" in clones.columns:
        reverse = pl.col("rev_comp").cast(pl.String).str.to_lowercase().is_in(["t", "true", "1"])
        complement = dict(zip("ACGT", "TGCA"))
        seq = pl.when(reverse).then(seq.str.reverse().str.replace_many(complement)).otherwise(seq)
    parts = [pl.col("candidate_sequences").fill_null([])] if "candidate_sequences" in clones.columns else []
    return clones.with_columns(pl.concat_list([*parts, seq]).list.eval(
        pl.element().filter(pl.element().is_not_null() & (pl.element().str.len_bytes() > 0)))
        .list.unique().list.sort().alias("candidate_sequences"))


def unique_clonotypes(clones: pl.DataFrame, *, naive_igm_only: bool = False,
                      gene_level: bool = False) -> pl.DataFrame:
    """Collapse an arda clonotype table to unique clonotypes = ``(v_call, j_call, junction)``.

    Isotype (``c_call``) is dropped from the key — isotype switching is the same clonotype. Rows are
    keyed to allele resolution (V/J allele) and full junction nt, so reads that differ only in
    alignment (CIGAR) or isotype collapse together; read support sums into ``count``.

    Args:
        clones: An arda clonotype frame (from :func:`annotate_reads`).
        naive_igm_only: Keep only IgM (naive B) clonotypes before collapsing (IGH only; no-op
            elsewhere) — for a naive-repertoire (pre-selection-independent) subset.
        gene_level: Strip allele suffixes, canonicalize ambiguous gene sets and merge their counts.

    Returns:
        Deduplicated clonotype frame: ``v_call, j_call, junction, junction_aa, locus, d_call,
        d2_call, count`` (one row per unique clonotype), sorted by descending ``count``.
        When full reads are available, ``candidate_sequences`` retains their coding-strand
        union for preliminary gene proposals; junctions remain the EM observations.
    """
    df = _with_candidate_sequences(clones)
    if gene_level:
        df = df.with_columns([
            pl.col(c).cast(pl.Utf8).str.split(",")
            .list.eval(pl.element().str.strip_chars().str.split("*").list.first())
            .list.unique().list.sort().list.join(",").alias(c)
            for c in ("v_call", "j_call", "d_call", "d2_call") if c in df.columns])
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
    count_cols = [pl.col(c).cast(pl.Int64, strict=False)
                  for c in ("duplicate_count", "count") if c in df.columns]
    cnt = pl.coalesce([*count_cols, pl.lit(1, dtype=pl.Int64)]).sum() if count_cols else pl.len().cast(pl.Int64)
    annotations = []
    for c in ("junction_aa", "locus", "d_call", "d2_call"):
        if c not in df.columns:
            expr = pl.lit(None, dtype=pl.Utf8)
        elif gene_level and c in ("d_call", "d2_call"):
            expr = (pl.col(c).drop_nulls().str.split(",").explode(empty_as_null=True).drop_nulls()
                    .unique().sort().str.join(",").replace("", None))
        else:
            expr = pl.col(c).first()
        annotations.append(expr.alias(c))
    if "candidate_sequences" in df.columns:
        annotations.append(pl.col("candidate_sequences").explode(empty_as_null=True).drop_nulls()
                           .unique().sort().alias("candidate_sequences"))
    out = df.group_by(["v_call", "j_call", "junction"]).agg(
        *annotations, cnt.alias("count"))
    return out.sort(["count", "v_call", "j_call", "junction"],
                    descending=[True, False, False, False], nulls_last=True)


def prepare(
    group: str,
    chain: str,
    label: str,
    *,
    out_dir: str | Path,
    cap: int | None = None,
    reconstruct: bool = False,
    naive_igm_only: bool = False,
    gene_level: bool = False,
    threads: int = 0,
    verbose: bool = False,
) -> pl.DataFrame:
    """Fetch → map (arda) → unique clonotypes for one ``(group, chain, label)`` bucket."""
    fq = fetch_fastq(group, chain, label)
    reads = annotate_reads(
        fq, out_dir=out_dir, prefix=f"{group}_{chain}_{label}",
        organism=ORGANISM[group], cap=cap, reconstruct=reconstruct, threads=threads, verbose=verbose,
    )
    return unique_clonotypes(reads, naive_igm_only=naive_igm_only, gene_level=gene_level)


#: Directory of the arda-mapped clonotype examples shipped with the source tree.
PREPARED_DIR = Path(__file__).resolve().parents[3] / "tests" / "python" / "fixtures" / "model_reads"

#: FASTA header fields, in order, after the record id. See :func:`write_prepared`.
_PREPARED_FIELDS = ("v_call", "j_call", "d_call", "d2_call", "count")


def write_prepared(clones: pl.DataFrame, path: str | Path) -> Path:
    """Write a clonotype frame as the gzipped FASTA :func:`load_prepared` reads.

    FASTA rather than a table because the payload *is* a sequence set — it stays greppable,
    readable by any bioinformatics tool, and about half the size of the equivalent TSV. The V/J/D
    calls EM needs for its per-read masks ride in the header, pipe-separated:

    ``>{id}|{v_call}|{j_call}|{d_call}|{d2_call}|{count}`` then the junction nucleotides.

    Empty fields are written as empty strings and read back as nulls.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = [c for c in _PREPARED_FIELDS if c in clones.columns]
    # Built as one column, not a Python string per clonotype: a real repertoire is 10^5-10^6 rows
    # and this is the hand-off to arda, so it runs on every EM build. Every field here is a string
    # or an integer, so a polars cast and `str()` agree exactly -- with a float they would not.
    if clones.height == 0:
        text = "\n"
    else:
        rec = (pl.lit(">") + pl.int_range(pl.len()).cast(pl.Utf8) + pl.lit("|")
               + pl.concat_str([pl.col(c).cast(pl.Utf8).fill_null("") for c in cols],
                               separator="|")
               + pl.lit("\n") + pl.col("junction").cast(pl.Utf8))
        text = clones.select(rec.str.join("\n").alias("r")).item() + "\n"
    if path.suffix == ".gz":
        import gzip

        path.write_bytes(gzip.compress(text.encode()))
    else:
        path.write_text(text)
    return path


def load_prepared(group: str = "human", chain: str = "TRB", label: str = "nonfunctional", *,
                  path: str | Path | None = None) -> pl.DataFrame:
    """Load a small **arda-mapped** clonotype example — no arda, no mmseqs2, no network, seconds.

    :func:`prepare` runs the full pipeline on the raw FASTQ, which needs arda + mmseqs2 and takes
    minutes per chain. For examples, notebooks, tests and a quick look at real junctions that is
    overkill, so a few already-mapped subsets ship with the source tree as gzipped FASTA (see
    :func:`write_prepared` for the header format).

    These live under ``tests/`` and are **not** packaged into the wheel — they exist for working in
    a checkout. From an installed vdjtools, use :func:`prepare` or pass ``path=``.

    Args:
        group: ``human``, ``human_fetal`` or ``mouse``.
        chain: e.g. ``"TRB"``.
        label: ``"nonfunctional"`` (the EM training bucket) or ``"functional"``.
        path: Read this file instead of looking one up.

    Returns:
        ``junction, v_call, j_call, d_call, d2_call, count`` — the columns
        :func:`unique_clonotypes` produces, minus the ones a FASTA cannot carry.

    Raises:
        FileNotFoundError: If no example ships for that combination. Only a couple do; build the
            rest with :func:`prepare` (or the whole corpus with :func:`build_all`).

    Example:
        >>> clones = load_prepared("human", "TRB")
        >>> model, report = infer_frame("TRB", clones, max_iter=10)
    """
    from .reference import read_fasta   # arda's parser, plus gzip handling

    src = Path(path) if path else PREPARED_DIR / f"{group}.{chain}.{label}.fa.gz"
    if not src.exists():
        raise FileNotFoundError(
            f"no prepared example at {src}. Only a couple ship; build this one with "
            f"data.prepare({group!r}, {chain!r}, {label!r}, out_dir=...), or the whole corpus "
            f"with data.build_all()."
        )
    rows = []
    for header, seq in read_fasta(src):
        parts = header.split("|")[1:]
        rec = {c: (parts[i] or None) if i < len(parts) else None
               for i, c in enumerate(_PREPARED_FIELDS)}
        rec["junction"] = seq
        rows.append(rec)
    if not rows:
        raise FileNotFoundError(f"{src} holds no records")
    # Explicit schema, not inference: d_call is null for every clonotype until the first read that
    # actually carried a D, which is well past any inference window.
    schema = {"junction": pl.Utf8, **{c: pl.Utf8 for c in _PREPARED_FIELDS}}
    return (pl.DataFrame(rows, schema=schema)
            .select(["junction", *_PREPARED_FIELDS])
            .with_columns(pl.col("count").cast(pl.Int64, strict=False)))


# --- full corpus build ------------------------------------------------------------------------

#: EM defaults for a corpus build. ``gene_prior`` is not optional in spirit: ``P(V) = 0`` is an
#: absorbing state of this EM, so without a Dirichlet pseudocount over the germline's functional
#: alleles one unlucky iteration deletes a real gene permanently (human TRB kept 30 of 57 V genes
#: unregularized, having seen 54 in the data).
BUILD_DEFAULTS = {"iters": 15, "tol": 1e-4, "gene_prior": 1.0, "nd_prior": 0.0, "prune_ratio": 1e-5}


def build_model(chain: str, *, group: str = "human", template=None, clones: pl.DataFrame | None = None,
                work_dir: str | Path = "/tmp/vdjtools_build", cap: int | None = None,
                iters: int = 15, tol: float = 1e-4, single_d: bool = False,
                nd_prior: float = 0.0, gene_prior: float = 1.0, threads: int = 0,
                ambiguous: str = "A", verbose: bool = False,
                checkpoint=None, checkpoint_every: int = 1,
                error_rate: float | None = 0.01, prune_ratio: float = 1e-5):
    """Fetch, annotate and EM-fit a model for one chain — the whole corpus pipeline, end to end.

    Args:
        chain: e.g. ``"TRB"``, ``"IGH"``.
        group: ``human``, ``human_fetal`` or ``mouse``. Sets the arda reference organism.
        template: Model supplying the gene set, germline and event graph. Defaults to
            :func:`~vdjtools.model.io.from_arda` for the chain and the group's organism.
        clones: Skip fetch+annotate and fit these clonotypes instead (an ``EM_DATA_DIR`` parquet,
            a :func:`load_prepared` subset, or your own reads). Preliminary V/J/D calls are
            not constraints; fresh batched ARDA proposals supply the V/J candidate sets.
        work_dir: Scratch directory for freshly generated arda per-read output.
        cap: Read cap for annotation. ``None`` uses every read — a subsample is a silent claim that
            the tail does not matter, and the tail is where the rare V genes are.
        iters: EM iteration cap. Convergence is on relative log-likelihood, so this is a safety net.
        tol: Stop when the relative log-likelihood improvement falls below this.
        single_d: Force a strict single-D model on a D-bearing locus.
        nd_prior: Dirichlet pseudocount pushing ``P(n_D=2)`` toward 0.
        gene_prior: Dirichlet pseudocount over the germline's functional V/J alleles — see
            :data:`BUILD_DEFAULTS`.
        threads: Annotation and E-step kernel threads (``0`` = auto).
        error_rate: Initial retained-germline substitution probability, fitted by EM.
            ``None`` retains exact matching. D/D2 paths remain unrestricted by annotation.
        ambiguous: Must be ``"A"``: DNA IUPAC sequencing symbols are replaced with A.
        prune_ratio: Relative best-scenario cutoff; 0 is exact. Actual omitted mass
            is recorded separately in the training report.

    Returns:
        ``(model, report, stats)`` — ``stats`` records ``chain, group, chain_type, n_clonotypes,
        n_used, p_nd2, loglik_first, loglik_last, iters, seconds``.

    Example:
        >>> model, rep, stats = build_model("TRB", work_dir="/tmp/em")
        >>> model.save("models/TRB")
    """
    import time
    from math import isfinite

    if not isfinite(prune_ratio) or not 0 <= prune_ratio < 1:
        raise ValueError("prune_ratio must be finite and between 0 (inclusive) and 1 (exclusive)")

    from .infer import candidate_masks, infer_native, print_progress
    from .io import from_arda

    started = time.perf_counter()
    organism = ORGANISM[group]
    base = template if template is not None else from_arda(chain, organism, gene_level=True)
    if clones is None:
        clones = prepare(group, chain, "nonfunctional", out_dir=work_dir, cap=cap,
                         threads=threads, verbose=verbose, gene_level=template is None)

    if template is None:
        clones = clones.with_columns([pl.lit(None, dtype=pl.String).alias(c)
                                      for c in ("v_call", "j_call") if c not in clones.columns])
        clones = unique_clonotypes(clones, gene_level=True)
    uniq = _filter_for_em(clones, base, ambiguous, use_calls=False)
    n_all = uniq.height
    seqs = [s.upper() for s in uniq["junction"].to_list()]
    uniq = _with_candidate_sequences(uniq)
    queries = uniq["candidate_sequences"].to_list() if "candidate_sequences" in uniq.columns else None
    masks, evidence = candidate_masks(seqs, base, organism=organism, threads=threads,
                                      query_sequences=queries)

    model, rep = infer_native(base, seqs, masks=masks, max_iter=iters, tol=tol,
                              single_d=single_d, nd_prior=nd_prior, error_rate=error_rate,
                              candidate_search=evidence,
                              prune_ratio=prune_ratio,
                              gene_prior=gene_prior, threads=threads,
                              progress=print_progress(prefix=f"[{chain}] ") if verbose else None,
                              checkpoint=checkpoint, checkpoint_every=checkpoint_every)
    # Stamp the builder that produced this model. A germline defect lives in the BUILDER, so
    # without this the only way to date a shipped model is to compare its posterior against an
    # unaffected reference fit — which is how the 3.9.1 J-anchor bug had to be sized.
    from dataclasses import replace

    from .. import __version__

    model.manifest = replace(model.manifest, builder_version=__version__)
    from .dd import d_count_probabilities
    nd = d_count_probabilities(model)
    stats = {
        "chain": chain, "group": group, "chain_type": base.chain_type,
        "n_clonotypes": clones.height, "n_used": n_all, "p_nd2": nd.get(2),
        "loglik_first": rep.loglik[0] if rep.loglik else None,
        "loglik_last": rep.loglik[-1] if rep.loglik else None,
        "iters": rep.n_iter, "seconds": round(time.perf_counter() - started, 1),
        "n_scoreable": rep.n_scoreable[-1] if rep.n_scoreable else 0,
        "candidate_search": evidence,
        "prune_ratio": prune_ratio,
    }
    return model, rep, stats


def _filter_for_em(clones: pl.DataFrame, base, ambiguous: str = "A", *,
                   use_calls: bool = True) -> pl.DataFrame:
    """Keep clonotypes the template can actually score: known V/J **gene**, encodable junction.

    ``use_calls=False`` (the corpus builder) retains junctions regardless of preliminary
    V/J assignments; fresh candidate search supplies the enumeration proposals.
    With explicit call filtering, filter on GENE, not allele. arda and a bootstrap model resolve alleles differently -- arda
    calls TRBV20-1*07, which OLGA's 89-allele index does not contain -- so an allele-level
    membership test deletes the WHOLE gene before `gene_masks` ever sees it. Measured on human TRB:
    TRBV20-1, the most-used human TRBV, went to zero training clonotypes and hence P(V)=0 in the
    shipped model. `gene_masks` already maps a call to every model allele of its gene, so the read
    is perfectly usable. Gene-level keeps 32,562 clonotypes and 54 V genes vs 24,980 and 51.

    DNA IUPAC sequencing symbols are always replaced with A; other policies are rejected.
    """
    from .infer import sanitize_junctions

    if not use_calls:
        kept = clones.filter(pl.col("junction").is_not_null() &
                             (pl.col("junction").str.len_bytes() > 0))
        return sanitize_junctions(kept, "junction", ambiguous=ambiguous, where="build_model")
    vgenes = {a.split("*")[0] for a in base.genomic["genes_v"]["v_allele"].to_list()}
    jgenes = {a.split("*")[0] for a in base.genomic["genes_j"]["j_allele"].to_list()}
    vg = pl.col("v_call").str.split(",").list.eval(pl.element().str.strip_chars().str.split("*").list.first())
    jg = pl.col("j_call").str.split(",").list.eval(pl.element().str.strip_chars().str.split("*").list.first())
    kept = clones.filter((vg.list.set_intersection(list(vgenes)).list.len() > 0) &
                         (jg.list.set_intersection(list(jgenes)).list.len() > 0))
    return sanitize_junctions(kept, "junction", ambiguous=ambiguous, where="build_model")


def _build_masks(uniq: pl.DataFrame, base) -> list[tuple]:
    """Per-read E-step masks with separate D1/D2 allele sets.

    AIRR writes an aligner tie comma-separated (``IGHD2-2*01,IGHD2-2*02``), so testing the whole
    joined string against a set of single alleles never matches and falls through to an empty mask
    = every D enumerated. On IGH that was 64% of reads: not just slow (35 D genes vs 3 for TRB) but
    a worse model, because the D marginal stops being anchored to what arda actually saw. A
    genuinely absent D1 call yields ``[]`` (unrestricted). A supplied D2 column
    adds a fourth mask: absent/unknown D2 yields ``[]`` and disables tandem paths.
    """
    from .infer import gene_masks, call_alleles, _gene_to_alleles

    masks = gene_masks(base, uniq["v_call"].to_list(), uniq["j_call"].to_list())
    if base.chain_type != "VDJ":
        return masks
    index = _gene_to_alleles(base,"d")
    calls = uniq["d_call"].to_list() if "d_call" in uniq.columns else [None]*uniq.height
    masks = [(mk[0],mk[1],call_alleles(index,call)) for mk,call in zip(masks,calls)]
    if "d2_call" in uniq.columns:
        masks = [(*mk,call_alleles(index,call)) for mk,call in zip(masks,uniq["d2_call"].to_list())]
    return masks


def build_all(chains=CHAINS, *, groups=("human",), workers: int | None = None, out_dir=None,
              **kw) -> dict:
    """Build models for several chains concurrently — the full-corpus entry point.

    Each ``(group, chain)`` runs the whole :func:`build_model` pipeline. Concurrency is threads,
    not processes: arda annotation is a subprocess (so it blocks on I/O, not the GIL) and the
    native E-step releases the GIL, so a thread pool gets real parallelism without the memory cost
    of forking the germline and count arrays per chain.

    Args:
        chains: Chains to build. Defaults to all seven.
        groups: Organism groups to build each chain for.
        workers: Concurrent builds; ``None`` uses the available cores, capped at the job count.
            Each build uses one kernel thread unless ``threads`` is explicitly supplied.
        out_dir: If given, each model is saved to ``out_dir/{group}_{chain}/``.
        **kw: Passed to :func:`build_model` (``iters``, ``tol``, ``cap``, ``gene_prior``, ...).

    Returns:
        ``{f"{group}_{chain}": {"model": ..., "report": ..., "stats": ...}}`` for the builds that
        succeeded, plus ``{"error": "..."}`` entries for those that did not — one failing chain
        never aborts the rest of a 30-minute corpus run.

    Example:
        >>> results = build_all(["TRB", "TRA"], workers=2, out_dir="models")
        >>> pl.DataFrame([r["stats"] for r in results.values() if "stats" in r])
    """
    from ..signature.cohort import parallel_rows

    if workers is not None and workers < 1:
        raise ValueError("workers must be positive")
    try:
        budget = kernel_threads(0 if workers is None else workers)
    except (TypeError, ValueError) as e:
        raise ValueError(f"invalid workers: {e}") from e
    jobs = [(g, c) for g in groups for c in chains]
    if not jobs:
        return {}
    n_workers = min(len(jobs), budget)
    kw = dict(kw)
    # Parallel builds own the CPU budget; each inner kernel is serial.
    kw.setdefault("threads", 1)
    if n_workers > 1 and kw["threads"] != 1:
        raise ValueError("parallel builds require threads=1 per worker")
    kw["threads"] = kernel_threads(kw["threads"])
    if n_workers > 1:
        require_serial_pools()

    def one(job):
        group, chain = job
        key = f"{group}_{chain}"
        try:
            model, rep, stats = build_model(chain, group=group, **kw)
            if out_dir is not None:
                model.save(Path(out_dir) / key)
            return key, {"model": model, "report": rep, "stats": stats}
        except Exception as e:  # noqa: BLE001 - one chain must not abort the corpus
            return key, {"error": f"{type(e).__name__}: {e}"}

    if n_workers == 1:
        return dict(map(one, jobs))
    return dict(parallel_rows(jobs, one, n_workers))
