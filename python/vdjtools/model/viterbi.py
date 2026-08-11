"""Most-likely nucleotide CDR3 for an amino-acid CDR3, and the recombination scenario behind it.

VDJdb-style records carry ``(V, J, CDR3aa)`` and no nucleotides. :func:`infer_nt` returns the
nucleotide CDR3 that maximises the model's generation probability, together with the argmax
scenario — which *is* the V/D/J boundary markup (``v_end``, ``d_start``, ``d_end``, ``j_start``).

**This is the argmax counterpart of a sum that already exists.** :func:`pgen.pgen_aa` marginalises
over every nucleotide sequence encoding a given amino-acid CDR3; the most likely such sequence is
the maximum of the same set, and :func:`pgen.pgen_nt`'s own loops — ``(V, len_v)``, ``(J, len_j)``,
then the D/insertion middle — already enumerate the scenarios. Nothing here re-derives the
recombination model: every probability comes from :func:`pgen.prepare`'s tables.

Two objectives are deliberately separated, because they are not the same number:

* :func:`infer_nt` maximises the **marginal** ``P(nt) = pgen_nt(nt)``, which sums over every
  scenario that could produce that sequence. That is the right target for "the most likely
  nucleotide sequence".
* :func:`best_scenario` then maximises over scenarios **for that fixed sequence**, which is what the
  boundary markup has to come from. The two maxima are over different things and taking the second
  as an answer to the first would be wrong.

⛔ **STATUS — read before using.**

* :func:`best_scenario` is **complete and validated**. It needs no enumeration: it is a max-product
  walk of the ``(V, len_v) x (J, len_j) x middle`` loops :func:`pgen.pgen_nt` already sums over.
  Checked on 500 generated draws per locus: the reported V span is the V germline 500/500, the J
  span the J germline 500/500, the D span the D germline 500/500 (TRB), and
  ``scenario_p <= pgen_nt`` 500/500 — a max cannot exceed the sum it is taken over.
* :func:`infer_nt_bruteforce` is **exact but exponential**, and exists as the TEST ORACLE.
* ⛔ **There is no production ``infer_nt`` yet.** It needs a max-product DP over the
  aa-constrained space (the argmax counterpart of :func:`pgen.pgen_aa`), which is not written.

**Why brute force cannot be the answer, measured on the real VDJdb (79,997 records).** The codon
search space — the product of synonymous-codon counts over the CDR3 — has a median of **5.3e6 for
TRA** and **1.9e7 for TRB**; only **8.9 %** / **1.6 %** of records are at or under 1e5 candidates.
Each candidate costs a full :func:`pgen.pgen_nt`, so enumeration is out by orders of magnitude.

⛔ **And it cannot be rescued by pinning the germline-templated flanks**, which was tried and is
**unsound**: pinning a codon to the V/J germline excludes every sequence in which that germline was
*trimmed*, and the true maximum can be one of them. Caught against the brute-force oracle —
``CAVSDMRF`` under ``TRAV21*01`` gives ``…GTGAGTGAC…`` where the pinned version returns
``…GTGAGCGAT…``. The pin is gone; do not reintroduce it as an optimisation.

⚠ Both functions assume an **in-frame** CDR3, i.e. ``len(nt) == 3 * len(aa)``. That is the VDJdb
population (real, productive receptors) but NOT every draw from :func:`generate` — an out-of-frame
rearrangement has e.g. 8 aa and 25 nt, and no in-frame nt can explain it. Test with
``productive_only=True`` or the comparison is meaningless.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product

from .model import Model
from .pgen import (_CODON_AA, _NUM2NT, _j_candidates, _p_insert, _v_candidates,
                   _Prepared, pgen_nt, prepare)

__all__ = ["Scenario", "best_scenario", "infer_nt_bruteforce", "codon_options"]


def _prep(m) -> _Prepared:
    """Accept a :class:`Model` or an already-:func:`pgen.prepare` -d one.

    `prepare` is not cheap and callers annotating a whole table should hoist it out of the loop, so
    both are accepted rather than forcing a re-prepare per record.
    """
    return m if isinstance(m, _Prepared) else prepare(m)

#: aa -> the codons encoding it as STRINGS, sorted so every result is reproducible.
#:
#: ⚠ `pgen._CODON_AA` is keyed by a TUPLE OF NUCLEOTIDE INDICES (`_NT2NUM`: A0 C1 G2 T3), not by a
#: string — `(3, 2, 1)` is "TGC". Everything here works in strings, so the conversion happens once.
_AA2CODONS: dict[str, list[str]] = {}
for _c, _a in _CODON_AA.items():
    _AA2CODONS.setdefault(_a, []).append("".join(_NUM2NT[i] for i in _c))
for _a in _AA2CODONS:
    _AA2CODONS[_a].sort()


@dataclass(frozen=True)
class Scenario:
    """A nucleotide CDR3 and the recombination that most likely produced it.

    Coordinates are **0-based, half-open, in CDR3 nucleotide space** (the conserved Cys through the
    conserved Phe/Trp, both included) — the same space as ``cdr3_nt`` itself. ``d_start``/``d_end``
    are ``None`` when the model has no D (a VJ chain) or no D was placed.
    """

    cdr3_nt: str
    v_call: str
    j_call: str
    v_end: int                      # nt contributed by the V germline: cdr3_nt[:v_end]
    j_start: int                    # nt contributed by the J germline: cdr3_nt[j_start:]
    d_call: str | None = None
    d_start: int | None = None
    d_end: int | None = None
    pgen: float = 0.0               # marginal P(cdr3_nt) — the quantity infer_nt maximises
    scenario_p: float = 0.0         # P of this single scenario; always <= pgen
    n_candidates: int = 0           # nt sequences actually scored
    truncated: bool = False         # the candidate set hit max_candidates
    runner_up_pgen: float = 0.0     # 2nd-best marginal, for the margin

    @property
    def margin(self) -> float:
        """``pgen / runner_up_pgen``. 1.0 means a tie; large means the argmax is unambiguous.

        ⚠ Report this. With a long non-templated core many sequences are near-equally likely and a
        bare "most likely" claim is then close to meaningless — the margin is what says so.
        """
        return self.pgen / self.runner_up_pgen if self.runner_up_pgen > 0 else float("inf")


def codon_options(aa: str) -> list[list[str]]:
    """Per residue, the codons encoding it. Unknown residues get no options (an empty list)."""
    return [_AA2CODONS.get(a, []) for a in aa.upper()]


def best_scenario(model_or_prep, cdr3_nt: str, v: str | None = None,
                  j: str | None = None) -> Scenario | None:
    """The single most likely recombination for a KNOWN nucleotide CDR3.

    A max-product walk of exactly the loops :func:`pgen.pgen_nt` sums over — ``(V, P(V), len_v)`` ×
    ``(J, P(J), len_j)`` × the middle — so the scenario space, the pruning and every probability are
    the shipped ones. Returns ``None`` when no (V, J) pair can explain the sequence.

    ⚠ The D placement is taken from the middle by the same gapless scan the model uses; ``d_call``
    is ``None`` for a VJ chain.
    """
    prep = _prep(model_or_prep)
    s = cdr3_nt.upper()
    vdj = prep.chain_type == "VDJ"
    v_cands = _v_candidates(prep, s, [v] if v else prep.functional_v)
    j_cands = _j_candidates(prep, s, [j] if j else prep.functional_j)

    best = None
    for V, pv, v_opts in v_cands:
        for J, j_opts in j_cands:
            pj = prep.p_j.get(J if vdj else (V, J), 0.0)
            if pj == 0.0:
                continue
            for len_v, _nv, p_dv in v_opts:
                for len_j, _nj, p_dj in j_opts:
                    if len_v + len_j > len(s):
                        continue
                    middle = s[len_v:len(s) - len_j]
                    # ⛔ MAX throughout, including over the D placement. Ranking (V, J, trims) by the
                    # SUM over middles and only then picking a D would answer a different question:
                    # the most likely (V, J, trims) MARGINALLY, whose best D need not lie on the
                    # single most likely path. This function's contract is one path.
                    if vdj:
                        inner, D, off, dlen = _best_vdj_middle(prep, J, middle)
                    else:
                        inner, D, off, dlen = _vj_middle(prep, middle), None, None, 0
                    if inner == 0.0:
                        continue
                    p = pv * pj * p_dv * p_dj * inner
                    if best is None or p > best[0]:
                        best = (p, V, J, len_v, len_j, D, off, dlen)
    if best is None:
        return None
    p, V, J, len_v, len_j, D, off, dlen = best
    d_call = D
    d_start = d_end = None
    if D is not None and off is not None:
        d_start, d_end = len_v + off, len_v + off + dlen
    return Scenario(cdr3_nt=s, v_call=V, j_call=J, v_end=len_v, j_start=len(s) - len_j,
                    d_call=d_call, d_start=d_start, d_end=d_end, scenario_p=p)


def _vj_middle(prep, middle: str) -> float:
    """P(the whole middle as VJ insertions). A VJ chain has one term, so max and sum coincide."""
    return _p_insert(middle, prep.p_ins["vj"], prep.R["vj"], prep.bias["vj"], from_right=False)


def _best_vdj_middle(prep, j: str, middle: str):
    """``(p, D, offset, length)`` — the single most likely D placement inside ``middle``.

    A max-product mirror of :func:`pgen._d_middle`, which sums
    ``P(D|J)·P(delD|D)·Pins(VD)·Pins(DJ)`` over every D, both trims and every position. Same terms,
    same tables, ``max`` instead of ``+`` — so the winner is the argmax of exactly the distribution
    the model defines, not a heuristic standing next to it.

    ⛔ ``P(D|J) == 0`` prunes the pair, which is where the GENOMIC constraint lives: TRBD2 lies 3' of
    the whole TRBJ1 cluster, so deletional joining can never produce a TRBD2-TRBJ1 pair and the
    model's table already encodes that as a zero. An earlier draft here picked the longest exact D
    substring and ignored ``j`` entirely — it would have happily called an impossible pair.

    ⛔ Not a re-ranking of D by a generative prior over an alignment score, either. That was built
    and measured in a sibling package and changed nothing (gene accuracy 98.9 -> 97.8 % on IGH,
    94.2 -> 94.5 % on human TRB): with 10-18 matched nt the alignment term is 11-20 nats and the
    prior moves it by ~3. Here the model IS the score, so the question does not arise.

    Ties are broken deterministically on ``(D name, offset, length)`` so a run is reproducible.
    """
    m = len(middle)
    maxdl, maxdr = prep.maxpal["d_5"], prep.maxpal["d_3"]
    pins_vd, R_vd, b_vd = prep.p_ins["vd"], prep.R["vd"], prep.bias["vd"]
    pins_dj, R_dj, b_dj = prep.p_ins["dj"], prep.R["dj"], prep.bias["dj"]
    best = (0.0, None, None, 0)
    for d in sorted(prep.functional_d):
        pdj = prep.p_d_given_j.get((j, d), 0.0)
        if pdj == 0.0:
            continue                                   # genomically impossible with this J
        cut = prep.cut["d"][d]
        for idx5 in range(len(cut) + 1):
            for idx3 in range(len(cut) - idx5 + 1):
                d_contrib = cut[idx5:len(cut) - idx3]
                pdel = prep.p_del["d"].get((d, idx5 - maxdl, idx3 - maxdr), 0.0)
                if pdel == 0.0:
                    continue
                ld = len(d_contrib)
                for pos in range(0, m - ld + 1):
                    if d_contrib and middle[pos:pos + ld] != d_contrib:
                        continue
                    w = _p_insert(middle[:pos], pins_vd, R_vd, b_vd, from_right=False)
                    if w == 0.0:
                        continue
                    w *= _p_insert(middle[pos + ld:], pins_dj, R_dj, b_dj, from_right=True)
                    p = pdj * pdel * w
                    if p > best[0]:
                        best = (p, d, pos, ld)
    return best


def infer_nt_bruteforce(model: Model, cdr3_aa: str, v: str | None = None, j: str | None = None, *,
                        max_candidates: int = 200_000) -> Scenario | None:
    """Exact most-likely nucleotide CDR3 by enumerating every codon assignment. **Oracle only.**

    Maximises the **marginal** ``pgen_nt`` over every nucleotide sequence encoding ``cdr3_aa``, then
    reports the argmax scenario for the winner via :func:`best_scenario`. Exact, and exponential in
    the number of residues.

    ⛔ **Not usable on real data** — see the module header: the median VDJdb record has 5.3e6 (TRA)
    to 1.9e7 (TRB) candidates. This exists so the production DP, when written, has an exact
    reference to be tested against on the small cases where both can run. ``max_candidates`` bounds
    the work and sets ``truncated``; a truncated result is the best of a PARTIAL set and is not the
    maximum, which is why it is flagged rather than returned quietly.

    Returns ``None`` when no candidate is explicable under the model — including the out-of-frame
    case, where no in-frame sequence of length ``3 * len(aa)`` exists.
    """
    prep = _prep(model)
    aa = cdr3_aa.upper()
    if not aa:
        return None
    # ⛔ NO germline pinning. It is unsound: pinning a codon to the V/J germline drops every
    # sequence in which that germline was trimmed, and the true maximum can be one of those.
    opts = codon_options(aa)
    if any(not o for o in opts):
        return None
    total = 1
    for o in opts:
        total *= len(o)
        if total > max_candidates:
            break
    truncated = total > max_candidates

    best_p = second_p = 0.0
    best_nt = None
    n = 0
    for combo in product(*opts):
        nt = "".join(combo)
        n += 1
        p = pgen_nt(prep, nt, v, j)
        if p > best_p:
            best_p, second_p, best_nt = p, best_p, nt
        elif p > second_p:
            second_p = p
        if n >= max_candidates:
            truncated = True
            break
    if best_nt is None:
        return None
    sc = best_scenario(prep, best_nt, v, j)
    if sc is None:
        return None
    return Scenario(cdr3_nt=sc.cdr3_nt, v_call=sc.v_call, j_call=sc.j_call, v_end=sc.v_end,
                    j_start=sc.j_start, d_call=sc.d_call, d_start=sc.d_start, d_end=sc.d_end,
                    pgen=best_p, scenario_p=sc.scenario_p, n_candidates=n,
                    truncated=truncated, runner_up_pgen=second_p)


def infer_nt(*_a, **_k):
    """⛔ NOT IMPLEMENTED — the production aa->nt inference.

    Needs a max-product DP with traceback over the aa-constrained space: the argmax counterpart of
    :func:`pgen.pgen_aa`, mirroring ``_aa_dp`` / ``_pgen_aa_vj`` / ``_pgen_aa_vdj`` / ``_d_aa_middle``
    with ``+`` replaced by ``max`` and back-pointers recorded. Use
    :func:`infer_nt_bruteforce` for small cases and as the reference it must reproduce.
    """
    raise NotImplementedError(
        "infer_nt needs the max-product DP over the aa-constrained space; it is not written. "
        "Use infer_nt_bruteforce for small cases (see the module header for why enumeration "
        "cannot scale, and why germline pinning is unsound)."
    )
