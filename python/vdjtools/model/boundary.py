"""Where the germline stops, as distinct from which history produced the junction.

:attr:`Scenario.v_end` is a property of the **argmax recombination history** — the most likely
``(V, delV, D, insertions, delJ, J)`` — and it is the right answer to that question. It is the
wrong answer to "where does the V germline stop", because maximising ``P(sequence)`` explains
N-region nucleotides as templated whenever it can: on 8,257 human TRB junctions with external
nucleotide truth it credits the germline two or more residues too far on 136 of them, against 8
for a plain alignment (``antigenomics/vdjtools#182``).

This module answers the boundary question instead, and needs **no DP and no recombination model**
— only the germline and the genetic code. Given the amino-acid junction it reads off how many
residues the germline explains, then hands the boundary codon to :func:`arda.cdr3fix.boundary_nt`,
which decides how far into it the germline still reaches. arda owns that decision because arda
owns the protein alignment and the germline; the whole measurement is in its docstring.

**The germline is the model's own**, never arda's anchor table: an OLGA-bootstrapped model keeps
OLGA germline and an arda-native one keeps arda's, and mixing the two inside one model is how a
Pgen goes quietly wrong (see ``CLAUDE.md``, "arda germline = single source of truth"). So this
reads ``model.genomic["genes_v"]["cdr3_segment"]`` — the germline proper, **not** ``cut_segment``,
whose palindromic P nucleotides are the reverse complement of the germline end and so are exactly
the nucleotides an aligner does *not* credit to the germline.

A V run and a J run may overlap on a short junction. That is not a defect and is not repaired: the
V-end / N / J-start partition is often not identifiable from sequence at all, and each side here
reports only what its own germline supports.
"""

from __future__ import annotations

from .model import Model
from .native import gene_to_allele
from .reference import translate

__all__ = ["germline_boundary"]


def _segments(model: Model, kind: str) -> dict[str, tuple[str, str]]:
    """Per allele, its CDR3-region germline and that germline's translation.

    Resolved once per allele rather than once per row -- a batch holds tens of alleles and
    hundreds of thousands of rows. The J germline is read from its 3' end, so its dangling 5'
    nucleotides (``len % 3``) come off before translating, exactly as arda's ``templated_aa`` does.
    """
    col, frame = (("v_allele", model.genomic["genes_v"]) if kind == "V"
                  else ("j_allele", model.genomic["genes_j"]))
    out: dict[str, tuple[str, str]] = {}
    for allele, seg in zip(frame[col], frame["cdr3_segment"]):
        if seg:
            out[allele] = (seg, translate(seg if kind == "V" else seg[len(seg) % 3:]))
    return out


def germline_boundary(model: Model, cdr3_aas, v=None, j=None):
    """The V/J boundary a germline alignment supports, in CDR3 **nucleotide** space.

    Args:
        model: the model whose germline to align against.
        cdr3_aas: amino-acid CDR3s (junction space: the conserved Cys through the conserved
            Phe/Trp, both included).
        v: per-row V call, or ``None`` to decline the V side for every row. A gene-level name
            resolves to its representative allele (:func:`native.gene_to_allele`) and a
            comma-separated call takes its first allele, which is what the boundary of a
            record with an ambiguous call can honestly mean.
        j: per-row J call, same rules.

    Returns:
        A :class:`polars.DataFrame` with **one row per input row, in input order**, carrying
        ``v_end`` (nucleotides the V germline explains, a half-open end) and ``j_start`` (index of
        the first nucleotide the J germline explains). Either is **null** when there is no germline
        to align against: no call, a call the model does not carry, or a non-functional allele with
        an empty ``cdr3_segment``. Nothing here declines because a junction is unexplainable -- that
        is the DP's business, and the reason to read this rather than :func:`infer_nt`'s output when
        what you want is a boundary.

    Measured against ``isalgo/airr_control``'s ``human.trb.ntvj``, on the **8,132** VDJdb human
    TRB junctions whose boundary every control observation agrees on and whose V and J the
    bundled ``olga:human_T_beta`` model carries, against :func:`infer_nt_batch` on the same rows:

    =========================================  ================  ================
    exact                                      ``Scenario``      this
    =========================================  ================  ================
    ``v.end``, VDJdb residues ``(nt + 1)//3``  7,244 (89.08 %)   7,554 (92.89 %)
    ``j.start``, VDJdb residues ``ceil/3``     7,483 (92.02 %)   7,966 (97.96 %)
    V boundary, nucleotides                    5,293 (65.09 %)   6,536 (80.37 %)
    J boundary, nucleotides                    5,301 (65.19 %)   6,056 (74.47 %)
    =========================================  ================  ================

    The same residue answer from an amino-acid alignment alone -- VDJdb's k-mer scanner, and
    ``arda.cdr3fix`` 2.30.1 -- is 71.79 %, so what the codon decision buys over rounding is 21
    points and what it buys over the argmax history is 4.

    NOTE the one column that gets worse: crediting the V germline **two or more residues** too far
    happens on 61 records here against 39 for the argmax history and 9 for the unextended
    alignment. Every one of them is the protein alignment over-running, not the codon decision --
    a residue past the real boundary that happens to encode the germline's residue extends the
    amino-acid run, and reaching 2 nucleotides into the codon after it then costs a whole residue.
    It is the same 16 % of records on which the residue count is wrong to begin with, and no
    amount of protein evidence can see it: the nucleotides that would give it away are exactly
    the ones not observed.

    Cost is **2.67 us/row** on this batch against :func:`infer_nt_batch`'s 82.81, i.e. 3.2 % --
    the per-allele germline and its translation are resolved once, and what is left per row is two
    string prefix runs and two codon lookups. Do not re-derive this; it is not the bottleneck.
    """
    import polars as pl
    from arda.cdr3fix import boundary_nt

    aas = [(a or "").upper() for a in cdr3_aas]
    n = len(aas)
    if v is not None and len(v) != n:
        raise ValueError("v must have the same length as cdr3_aas")
    if j is not None and len(j) != n:
        raise ValueError("j must have the same length as cdr3_aas")
    alias = gene_to_allele(model)
    germ_v, germ_j = _segments(model, "V"), _segments(model, "J")
    v_end: list[int | None] = [None] * n
    j_start: list[int | None] = [None] * n

    def allele(name):
        if not name:
            return None
        first = name.split(",")[0].strip()
        return alias.get(first, first)

    for i, aa in enumerate(aas):
        if not aa:
            continue
        got = germ_v.get(allele(v[i])) if v is not None else None
        if got is not None:
            seg, prot = got
            m, stop = 0, min(len(prot), len(aa))
            while m < stop and aa[m] == prot[m]:
                m += 1
            v_end[i] = boundary_nt(aa, m, seg, "V")[0]
        got = germ_j.get(allele(j[i])) if j is not None else None
        if got is not None:
            seg, prot = got
            m, stop = 0, min(len(prot), len(aa))
            while m < stop and aa[-1 - m] == prot[-1 - m]:
                m += 1
            j_start[i] = boundary_nt(aa, len(aa) - m, seg, "J")[0]
    return pl.DataFrame({"v_end": pl.Series(v_end, dtype=pl.Int64),
                         "j_start": pl.Series(j_start, dtype=pl.Int64)})
