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
        ImportError: If arda (the ``[model]`` extra) is not installed.
        ValueError: If no germline is found for ``locus`` / ``organism``.
    """
    from arda.cdr3fix import load_anchors  # optional dep (the [model] extra)
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
        ImportError: If arda (the ``[model]`` extra) is not installed.
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
