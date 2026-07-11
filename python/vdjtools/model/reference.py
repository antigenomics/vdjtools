"""Germline reference — **arda is the single source of germline truth**.

Every V/D/J germline sequence and CDR3 anchor used anywhere in vdjtools (annotation,
scenario enumeration, contig stitching, generation, and arda-native models) resolves from
arda's germline library by allele name, so the whole pipeline speaks one coordinate frame.

arda's anchor convention is **byte-identical to OLGA's**: ``anchor_nt`` is a 0-based offset
into the full germline marking the conserved Cys104 codon (V) or [FW]118 codon (J); the
CDR3-region germline is ``full[anchor:]`` for V and ``full[:anchor+3]`` for J. So no coordinate
conversion is needed between the two — this module documents and enforces that shared frame
(:func:`reconcile_olga` catalogs any residual sequence differences).

Note (arda coverage today): arda ships the **CDR3-region** germline for V/J (what Pgen needs)
plus **full** D germlines. Full-length V/J germline — needed for contig **stitching** — is not
in arda's shipped data yet and is a prerequisite for Phase 1c (either a new arda helper or the
OLGA model's own full germline for bootstrap models).
"""
from __future__ import annotations

from functools import lru_cache

import polars as pl

_COMP = str.maketrans("ACGT", "TGCA")


def reverse_complement(seq: str) -> str:
    """Reverse complement of a nucleotide string (ACGT)."""
    return seq.translate(_COMP)[::-1]


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


def reconcile_olga(model, *, tol_segment: str = "cdr3_segment") -> pl.DataFrame:
    """Catalog how an OLGA-loaded model's germline relates to arda's (the shared-frame audit).

    OLGA bootstrap models keep OLGA's germline geometry for exact-Pgen fidelity; this reports,
    per V/J allele, whether it resolves in arda and whether the CDR3-region germline matches —
    so divergences (IMGT-version drift, V 3' extent) are flagged, never silent.

    Args:
        model: A :class:`~vdjtools.model.model.Model` loaded via ``from_olga``.
        tol_segment: Which genomic column holds the OLGA CDR3-region germline to compare.

    Returns:
        Per-allele report: ``allele, segment, in_arda, germline_equal, olga_len, arda_len``.
    """
    arda = load_germline(model.locus, model.organism)
    arda_seq = {(r[0], r[1]): r[2] for r in arda.select(["segment", "allele", "sequence"]).iter_rows()}
    rows = []
    for seg in ("v", "j"):
        g = model.genomic[f"genes_{seg}"]
        for allele, olga_seq in zip(g[f"{seg}_allele"], g[tol_segment]):
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
