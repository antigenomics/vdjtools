"""Germline reference — **arda is the single source of germline truth**.

Every V/D/J germline sequence and CDR3 anchor used anywhere in vdjtools (annotation,
scenario enumeration, contig stitching, generation, and arda-native models) resolves from
arda's germline library by allele name, so the whole pipeline speaks one coordinate frame.

arda's anchor convention is **byte-identical to OLGA's**: ``anchor_nt`` is a 0-based offset
into the full germline marking the conserved Cys104 codon (V) or [FW]118 codon (J); the
CDR3-region germline is ``full[anchor:]`` for V and ``full[:anchor+3]`` for J. So no coordinate
conversion is needed between the two — this module documents and enforces that shared frame
(:func:`reconcile_olga` catalogs any residual sequence differences).

Full-length V/J germline for contig **stitching** (Phase 1c) is recovered from arda's bundled
scaffold reference by :func:`load_full_vj_germline` / :func:`arda_full_germline`: arda's
``alleles.fasta`` scaffolds are ``full_V + N-pad + full_J``, sliced per allele at the
``v_sequence_end`` / ``j_sequence_start`` boundaries from ``arda.annotate.reference``. The
sliced full germline is anchor-consistent with :func:`load_germline` (verified across all
functional V/J alleles: ``full_V[anchor:]`` == the CDR3-region germline, and the CDR3-region J
germline is a prefix of ``full_J``). Only functional/ORF alleles with complete markup have a
full germline in arda; pseudogenes (present in the CDR3 anchors) may be absent.
"""
from __future__ import annotations

from functools import lru_cache

import polars as pl

_COMP = str.maketrans("ACGT", "TGCA")

# Standard genetic code (DNA codons -> amino acid; '*' = stop).
_BASES = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_CODON_TABLE = {a + b + c: _AA[i * 16 + j * 4 + k]
                for i, a in enumerate(_BASES) for j, b in enumerate(_BASES) for k, c in enumerate(_BASES)}


def reverse_complement(seq: str) -> str:
    """Reverse complement of a nucleotide string (ACGT)."""
    return seq.translate(_COMP)[::-1]


def translate(seq: str) -> str:
    """Translate a nucleotide string to amino acids (standard code; trailing partial codon dropped)."""
    return "".join(_CODON_TABLE[seq[i:i + 3]] for i in range(0, len(seq) - len(seq) % 3, 3))


def _append_3p_palindrome(seq: str, k: int) -> str:
    n = min(len(seq), k)
    return seq + reverse_complement(seq[len(seq) - n:]) if n else seq


def _prepend_5p_palindrome(seq: str, k: int) -> str:
    n = min(len(seq), k)
    return reverse_complement(seq[:n]) + seq if n else seq


def cut_segment(seq: str, segment: str, max_pal: int) -> str:
    """Palindrome-extend a CDR3-region germline for the Pgen DP (mirrors OLGA's cutR/cutL_seq).

    Appends up to ``max_pal`` reverse-complement (P-nucleotide) bases on the trimmable end so
    a deletion *index* into the result directly counts nt removed:

    - V: append at the 3' end.  - J: prepend at the 5' end.

    Args:
        seq: CDR3-region germline (V: Cys104→3' end; J: 5'→[FW]118 codon end).
        segment: ``"V"`` or ``"J"``.
        max_pal: Maximum palindromic nt for this end.
    """
    if segment == "V":
        return _append_3p_palindrome(seq, max_pal)
    if segment == "J":
        return _prepend_5p_palindrome(seq, max_pal)
    raise ValueError(f"cut_segment: segment must be 'V' or 'J', got {segment!r} (use cut_segment_d for D)")


def cut_segment_d(seq: str, max_pal5: int, max_pal3: int) -> str:
    """Palindrome-extend a D germline on both ends (5' then 3')."""
    return _append_3p_palindrome(_prepend_5p_palindrome(seq, max_pal5), max_pal3)


@lru_cache(maxsize=16)
def load_germline(locus: str, organism: str = "human") -> pl.DataFrame:
    """V/D/J germline + CDR3 anchors for a locus, from arda (the source of truth).

    Args:
        locus: e.g. ``"TRB"``, ``"TRA"``, ``"IGH"``.
        organism: e.g. ``"human"``, ``"mouse"``.

    Returns:
        One row per allele with columns ``allele, gene, segment, sequence, cdr3_anchor,
        functionality, functional, status``. For V/J, ``sequence`` is the CDR3-region germline
        and ``cdr3_anchor`` the 0-based anchor codon offset in the *full* germline; for D,
        ``sequence`` is the full D germline, ``cdr3_anchor = -1``.

    Raises:
        ImportError: If arda is not importable (it is a base dependency; a plain
            ``pip install vdjtools`` ships it).
        ValueError: If no germline is found for ``locus`` / ``organism``.
    """
    from arda.cdr3fix import load_anchors  # base dep; imported lazily to keep `import vdjtools` light
    from arda.paths import vdj_dir
    from arda.refbuild.imgt import read_fasta

    rows = []
    for (segment, allele), an in load_anchors(organism).items():
        if an.locus != locus:
            continue
        rows.append(
            {
                "allele": allele,
                "gene": allele.split("*")[0],
                "segment": segment,
                "sequence": an.germline_nt,
                "cdr3_anchor": an.anchor_nt,
                "functionality": an.functionality,
                "functional": an.functionality == "F",
                "status": an.status,
            }
        )
    d_path = vdj_dir(organism) / "d_germlines.fasta"
    if d_path.exists():
        for header, seq in read_fasta(d_path):
            d_locus, allele = header.split("|", 1)
            if d_locus != locus:
                continue
            rows.append(
                {
                    "allele": allele,
                    "gene": allele.split("*")[0],
                    "segment": "D",
                    "sequence": seq,
                    "cdr3_anchor": -1,
                    "functionality": "",
                    "functional": True,
                    "status": "ok",
                }
            )
    if not rows:
        raise ValueError(f"no arda germline for locus {locus!r} / organism {organism!r}")
    return pl.DataFrame(rows)


@lru_cache(maxsize=4)
def load_full_vj_germline(organism: str = "human") -> dict[tuple[str, str], str]:
    """Full-length V and J germline nucleotide sequences from arda, by ``(segment, allele)``.

    arda ships full V/J germline only inside deduplicated V–J **scaffolds**
    (``database/vdj/<organism>/alleles.fasta``, keyed by opaque scaffold id); this slices each
    scaffold at the ``v_sequence_end`` / ``j_sequence_start`` boundaries from
    ``arda.annotate.reference.load_reference`` to recover the per-allele full germline (V: FR1 →
    3' end of V-REGION; J: 5' J → end of FR4). The first scaffold carrying an allele wins (all
    are byte-identical for that allele's segment).

    Args:
        organism: e.g. ``"human"``, ``"mouse"``.

    Returns:
        ``{("V"|"J", allele): full_germline_nt}``. Only functional/ORF alleles with complete
        arda markup are present; pseudogenes (in the CDR3 anchors) may be missing.

    Raises:
        ImportError: If arda is not importable (it is a base dependency; a plain
            ``pip install vdjtools`` ships it).
    """
    from arda.annotate.reference import load_reference
    from arda.refbuild.imgt import read_fasta

    ref = load_reference(organism, "nt")
    scaffolds = dict(read_fasta(ref.target_fasta))
    out: dict[tuple[str, str], str] = {}
    for sid, e in ref.entries.items():
        scaf = scaffolds.get(sid)
        if not scaf:
            continue
        v_end = getattr(e, "v_sequence_end", None)
        if v_end:
            for allele in e.v_call.split(","):
                out.setdefault(("V", allele.strip()), scaf[:v_end])
        j_start = getattr(e, "j_sequence_start", None)
        if j_start:
            end = getattr(e, "vj_end", None) or len(scaf)
            for allele in e.j_call.split(","):
                out.setdefault(("J", allele.strip()), scaf[j_start - 1:end])
    return out


@lru_cache(maxsize=8)
def arda_full_germline(locus: str, organism: str = "human") -> dict[tuple[str, str], tuple[str, int]]:
    """Stitch-ready full V/J germline + anchor for a locus, entirely from arda.

    Combines :func:`load_full_vj_germline` (full-length germline) with :func:`load_germline`
    (the CDR3-region germline) so the anchor is derived self-consistently by length — no reliance
    on a cross-source coordinate assumption:

    * **V**: ``anchor = len(full) - len(cdr3_region)`` — ``full[:anchor]`` is the framework 5'
      of the conserved Cys104.
    * **J**: ``anchor = len(cdr3_region) - 3`` — ``full[anchor + 3:]`` is the framework 3' of the
      conserved Phe/Trp118 codon.

    This is exactly the ``(full_germline, anchor)`` pair :func:`vdjtools.model.stitch.stitch_contig`
    consumes, so an arda-native model (no OLGA germline) can stitch full contigs. Alleles whose
    full germline is absent from arda (pseudogenes / incomplete markup) are skipped.

    Args:
        locus: e.g. ``"TRB"``, ``"IGH"``.
        organism: e.g. ``"human"``.

    Returns:
        ``{("V"|"J", allele): (full_germline_nt, anchor)}``.
    """
    full = load_full_vj_germline(organism)
    gl = load_germline(locus, organism)
    out: dict[tuple[str, str], tuple[str, int]] = {}
    for r in gl.filter(pl.col("segment").is_in(["V", "J"])).iter_rows(named=True):
        seg, allele, cdr3_region = r["segment"], r["allele"], r["sequence"]
        fg = full.get((seg, allele))
        if not fg or not cdr3_region:
            continue
        if seg == "V":
            if not fg.endswith(cdr3_region):
                continue  # defensive: full V must end with the CDR3-region germline
            anchor = len(fg) - len(cdr3_region)
        else:
            if not fg.startswith(cdr3_region):
                continue  # defensive: CDR3-region J must prefix the full J
            anchor = len(cdr3_region) - 3
        out[(seg, allele)] = (fg, anchor)
    return out


# --- custom germline libraries ---------------------------------------------------------------

#: Columns a germline frame must carry, and the defaults filled in for the optional ones.
GERMLINE_REQUIRED = ("allele", "segment", "sequence")
GERMLINE_OPTIONAL = {"gene": None, "functional": True, "cdr3_anchor": -1, "full_germline": ""}

#: Column schema of the tidy issue frame shared by :func:`validate_germline` and ``check_model``.
_ISSUE_SCHEMA = {"severity": pl.Utf8, "check": pl.Utf8, "event": pl.Utf8, "segment": pl.Utf8,
                 "allele": pl.Utf8, "detail": pl.Utf8, "value": pl.Float64}


def _issue(severity: str, check: str, detail: str, *, event=None, segment=None, allele=None,
           value=None) -> dict:
    """One row of the shared issue frame (same schema as :func:`vdjtools.model.check.check_model`)."""
    return {"severity": severity, "check": check, "event": event, "segment": segment,
            "allele": allele, "detail": detail, "value": value}


def normalize_germline(germline: pl.DataFrame) -> pl.DataFrame:
    """Fill a germline frame's optional columns with their defaults; returns a new frame.

    ``gene`` defaults to ``allele.split("*")[0]``. See :data:`GERMLINE_OPTIONAL` for the rest.
    """
    df = germline
    if "gene" not in df.columns:
        df = df.with_columns(pl.col("allele").str.split("*").list.first().alias("gene"))
    for col, default in GERMLINE_OPTIONAL.items():
        if col not in df.columns:
            df = df.with_columns(pl.lit(default).alias(col))
    return df


def validate_germline(germline: pl.DataFrame) -> pl.DataFrame:
    """Audit a germline frame destined for :func:`vdjtools.model.io.from_germline`.

    Catches the mistakes that otherwise produce a model that builds cleanly and scores wrongly —
    above all a **misplaced CDR3 anchor**, which shifts every deletion profile by a constant and is
    invisible downstream.

    Args:
        germline: A frame with at least ``allele, segment, sequence`` (see
            :data:`GERMLINE_REQUIRED`); optional columns are described by :data:`GERMLINE_OPTIONAL`.
            ``sequence`` is the **CDR3-region** germline (V: Cys104 codon → 3' end; J: 5' end →
            through the [FW]118 codon), or the full germline for D.

    Returns:
        Tidy issue frame ``severity, check, event, segment, allele, detail, value``; empty when the
        library is clean. ``severity == "error"`` means the frame cannot build a model.
    """
    rows: list[dict] = []
    missing = [c for c in GERMLINE_REQUIRED if c not in germline.columns]
    if missing:
        rows.append(_issue("error", "germline_columns",
                           f"missing required column(s): {', '.join(missing)}"))
        return pl.DataFrame(rows, schema=_ISSUE_SCHEMA)

    df = normalize_germline(germline)
    bad_seg = sorted(set(df["segment"].to_list()) - {"V", "D", "J"})
    if bad_seg:
        rows.append(_issue("error", "germline_segment",
                           f"segment must be V, D or J; got {bad_seg}"))
    counts = {s: df.filter(pl.col("segment") == s).height for s in ("V", "D", "J")}
    for seg in ("V", "J"):
        if not counts[seg]:
            rows.append(_issue("error", "germline_segment_missing",
                               f"no {seg} alleles — a model needs at least one", segment=seg))

    dup = (df.group_by("allele").len().filter(pl.col("len") > 1))
    for r in dup.iter_rows(named=True):
        rows.append(_issue("error", "germline_duplicate_allele",
                           f"allele {r['allele']!r} appears {r['len']} times",
                           allele=r["allele"], value=float(r["len"])))

    for r in df.iter_rows(named=True):
        allele, seg, seq = r["allele"], r["segment"], r["sequence"] or ""
        if not allele:
            rows.append(_issue("error", "germline_empty_allele", "empty allele name", segment=seg))
            continue
        if "*" not in allele:
            # The model is keyed by ALLELE and `call_alleles`/`native._gene_idx` split on "*".
            # A gene-level name here is the documented "2.38x too high" trap, surfaced at build time.
            rows.append(_issue("warn", "germline_gene_level_name",
                               f"{allele!r} has no '*NN' allele suffix; the model is allele-keyed",
                               segment=seg, allele=allele))
        if not seq:
            rows.append(_issue("error", "germline_empty_sequence",
                               f"{allele!r} has an empty germline sequence", segment=seg, allele=allele))
            continue
        if set(seq) - set("ACGT"):
            # Dropped rather than fatal: IUPAC-ambiguous germline cannot be encoded by the native
            # DP, and a handful of IMGT alleles carry one. Reported instead of silently vanishing.
            rows.append(_issue("warn", "germline_ambiguous",
                               f"{allele!r} has non-ACGT bases; it will be dropped from the model",
                               segment=seg, allele=allele))
            continue
        if seg == "V" and translate(seq[:3]) != "C":
            rows.append(_issue("warn", "germline_anchor_frame",
                               f"{allele!r} does not start with a Cys codon ({seq[:3]!r} -> "
                               f"{translate(seq[:3])!r}); the CDR3 anchor looks misplaced",
                               segment=seg, allele=allele))
        if seg == "J" and len(seq) >= 3 and translate(seq[-3:]) not in ("F", "W"):
            rows.append(_issue("warn", "germline_anchor_frame",
                               f"{allele!r} does not end with a Phe/Trp codon ({seq[-3:]!r} -> "
                               f"{translate(seq[-3:])!r}); the CDR3 anchor looks misplaced",
                               segment=seg, allele=allele))
    return pl.DataFrame(rows, schema=_ISSUE_SCHEMA)


def read_fasta(path) -> list[tuple[str, str]]:
    """Read a FASTA into ``(header, sequence)`` pairs, transparently handling ``.gz``.

    Parsing itself is delegated to arda (one FASTA parser in the stack, not two), but arda's reader
    takes a path and plain ``open()``, so a gzipped file reaches it as mojibake and dies on the
    first byte. Gzipped input is therefore decompressed to a temporary file first.
    """
    import gzip
    import tempfile
    from pathlib import Path

    from arda.refbuild.imgt import read_fasta as _read

    path = Path(path)
    if path.suffix != ".gz":
        return _read(path)
    with gzip.open(path, "rt") as src, \
            tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as tmp:
        tmp.write(src.read())
        name = tmp.name
    try:
        return _read(Path(name))
    finally:
        Path(name).unlink(missing_ok=True)


def read_germline_fasta(v, j, d=None, *, anchors=None) -> pl.DataFrame:
    """Build a germline frame from your own FASTA files — the entry point for a custom library.

    The segment comes from **which argument a file was passed as**, so no header convention is
    assumed: a header is either ``>ALLELE`` or ``>ANYTHING|ALLELE`` (arda's D convention), and
    everything after the second ``|`` is ignored.

    Args:
        v: FASTA path for the V alleles.
        j: FASTA path for the J alleles.
        d: Optional FASTA path for the D alleles. Supplying it makes the model ``VDJ``; omitting
            it makes it ``VJ``.
        anchors: Optional CDR3-anchor CSV in OLGA's ``*_gene_CDR3_anchors.csv`` format
            (``gene,anchor_index,function``). When given, the FASTAs are taken to hold
            **full-length** germline and are sliced to the CDR3 region (``full[anchor:]`` for V,
            ``full[:anchor + 3]`` for J). When omitted, the V/J sequences are taken to be
            CDR3-region germline already — :func:`validate_germline` flags it if they are not.

    Returns:
        A germline frame in :func:`load_germline`'s schema, ready for
        :func:`vdjtools.model.io.from_germline`.
    """
    from pathlib import Path

    from .io import _read_anchors

    anchor_map = _read_anchors(Path(anchors)) if anchors else {}
    rows = []
    for path, segment in ((v, "V"), (j, "J"), (d, "D")):
        if path is None:
            continue
        for header, seq in read_fasta(Path(path)):
            parts = header.split("|")
            allele = (parts[1] if len(parts) > 1 else parts[0]).strip()
            seq = seq.strip().upper()
            anchor, functionality = anchor_map.get(allele, (-1, "F"))
            full = ""
            if segment in ("V", "J") and anchor >= 0:
                full, seq = seq, (seq[anchor:] if segment == "V" else seq[:anchor + 3])
            rows.append({
                "allele": allele, "gene": allele.split("*")[0], "segment": segment,
                "sequence": seq, "cdr3_anchor": anchor, "full_germline": full,
                "functionality": functionality,
                "functional": segment == "D" or functionality == "F",
                "status": "ok",
            })
    if not rows:
        raise ValueError("no FASTA records read — check the paths and that the files are not empty")
    return pl.DataFrame(rows)


def reconcile_olga(model) -> pl.DataFrame:
    """Catalog how an OLGA-loaded model's germline relates to arda's (the shared-frame audit).

    OLGA bootstrap models keep OLGA's germline geometry for exact-Pgen fidelity; this reports,
    per V/J allele, whether it resolves in arda and whether the CDR3-region germline matches —
    so divergences (IMGT-version drift, V 3' extent) are flagged, never silent.

    Args:
        model: A :class:`~vdjtools.model.model.Model` loaded via ``from_olga``.

    Returns:
        Per-allele report: ``allele, segment, in_arda, germline_equal, olga_len, arda_len``.
    """
    arda = load_germline(model.locus, model.organism)
    arda_seq = {(r[0], r[1]): r[2] for r in arda.select(["segment", "allele", "sequence"]).iter_rows()}
    rows = []
    for seg in ("v", "j"):
        g = model.genomic[f"genes_{seg}"]
        for allele, olga_seq in zip(g[f"{seg}_allele"], g["cdr3_segment"]):
            key = (seg.upper(), allele)
            a = arda_seq.get(key)
            rows.append(
                {
                    "allele": allele,
                    "segment": seg.upper(),
                    "in_arda": a is not None,
                    "germline_equal": a is not None and a == olga_seq,
                    "olga_len": len(olga_seq),
                    "arda_len": len(a) if a is not None else -1,
                }
            )
    return pl.DataFrame(rows)
