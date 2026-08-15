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

**STATUS.**

* :func:`best_scenario` is **complete and validated**. It needs no enumeration: it is a max-product
  walk of the ``(V, len_v) x (J, len_j) x middle`` loops :func:`pgen.pgen_nt` already sums over.
  Checked on 500 generated draws per locus: the reported V span is the V germline 500/500, the J
  span the J germline 500/500, the D span the D germline 500/500 (TRB), and
  ``scenario_p <= pgen_nt`` 500/500 — a max cannot exceed the sum it is taken over.
* :func:`infer_nt` is the **production** path: a codon-constrained max-product DP over the same
  scenario enumeration :func:`pgen.pgen_aa` sums over, with exact branch-and-bound, followed by an
  exact :func:`pgen.pgen_nt` re-score of the top candidates.
* :func:`infer_nt_bruteforce` is **exact but exponential**, and exists as the TEST ORACLE.

**Why brute force cannot be the answer, measured on the real VDJdb (79,997 records).** The codon
search space — the product of synonymous-codon counts over the CDR3 — has a median of **5.3e6 for
TRA** and **1.9e7 for TRB**; only **8.9 %** / **1.6 %** of records are at or under 1e5 candidates.
Each candidate costs a full :func:`pgen.pgen_nt`, so enumeration is out by orders of magnitude.

WARNING: **And it cannot be rescued by pinning the germline-templated flanks**, which was tried and is
**unsound**: pinning a codon to the V/J germline excludes every sequence in which that germline was
*trimmed*, and the true maximum can be one of them. Caught against the brute-force oracle —
``CAVSDMRF`` under ``TRAV21*01`` gives ``…GTGAGTGAC…`` where the pinned version returns
``…GTGAGCGAT…``. The pin is gone; do not reintroduce it as an optimisation.

NOTE: Both functions assume an **in-frame** CDR3, i.e. ``len(nt) == 3 * len(aa)``. That is the VDJdb
population (real, productive receptors) but NOT every draw from :func:`generate` — an out-of-frame
rearrangement has e.g. 8 aa and 25 nt, and no in-frame nt can explain it. Test with
``productive_only=True`` or the comparison is meaningless.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product

from .model import Model
from .pgen import (_CODON_AA, _CODON_TABLE, _NT2NUM, _NUM2NT, _aa_germline_prefix_ok,
                   _d_codons_ok, _j_candidates, _p_insert, _v_candidates,
                   _Prepared, pgen_nt, prepare)

__all__ = ["Scenario", "best_scenario", "infer_nt", "infer_nt_bruteforce", "codon_options"]


#: ``id(model) -> (model, prepared)``. The model reference is stored *and verified*, not just its
#: id: CPython reuses ids, and a stale hit across a TRB->TRD switch in one process is exactly the
#: bug already recorded against ``native.pack``.
_PREP_CACHE: dict = {}


def _prep(m) -> _Prepared:
    """Accept a :class:`Model` or an already-:func:`pgen.prepare` -d one.

    `prepare` is not cheap and callers annotating a whole table should hoist it out of the loop, so
    both are accepted rather than forcing a re-prepare per record.
    """
    if isinstance(m, _Prepared):
        return m
    hit = _PREP_CACHE.get(id(m))
    if hit is not None and hit[0] is m:
        return hit[1]
    p = prepare(m)
    _PREP_CACHE[id(m)] = (m, p)
    return p

#: aa -> the codons encoding it as STRINGS, sorted so every result is reproducible.
#:
#: NOTE: `pgen._CODON_AA` is keyed by a TUPLE OF NUCLEOTIDE INDICES (`_NT2NUM`: A0 C1 G2 T3), not by a
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

        NOTE: Report this. With a long non-templated core many sequences are near-equally likely and a
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

    NOTE: The D placement is taken from the middle by the same gapless scan the model uses; ``d_call``
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
                    # WARNING: MAX throughout, including over the D placement. Ranking (V, J, trims) by the
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

    WARNING: ``P(D|J) == 0`` prunes the pair, which is where the GENOMIC constraint lives: TRBD2 lies 3' of
    the whole TRBJ1 cluster, so deletional joining can never produce a TRBD2-TRBJ1 pair. An earlier
    draft here picked the longest exact D substring and ignored ``j`` entirely — it would have
    happily called an impossible pair.

    NOTE: This inherits the constraint, it does not impose one: it is only as good as the table. An
    earlier revision of this note claimed the model "already encodes that as a zero", and that was
    **not true** — OLGA's own TRB model gives ``P(TRBD2*01 | TRBJ1-6*01) = 0.333``, and EM relearned
    the pair from noisy short-read D calls until the constraint was added to the M-step. Models
    fitted from :mod:`vdjtools.model.infer` now carry the zeros;
    :func:`~vdjtools.model.infer.enforce_dj_order` repairs one that does not, and
    :func:`~vdjtools.model.check.check_model` reports any survivor.

    WARNING: Not a re-ranking of D by a generative prior over an alignment score, either. That was built
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

    WARNING: **Not usable on real data** — see the module header: the median VDJdb record has 5.3e6 (TRA)
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
    # WARNING: NO germline pinning. It is unsound: pinning a codon to the V/J germline drops every
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


def _aa_dp_max(aa: str, template: list[int], specs: list, floor: float = 0.0):
    """Max-product counterpart of :func:`pgen._aa_dp` with traceback: ``(p, nt)`` or ``(0.0, None)``.

    Identical states, identical factors, ``max`` where ``_aa_dp`` sums — so the germline positions
    stay pinned to their segment and every free (N-region) position is scored by exactly the
    insertion model that covers it: ``P(nt_1)`` from ``bias`` then ``P(nt_k | nt_{k-1})`` from the
    VD / DJ / VJ dinucleotide matrix. The ``P(insLen)`` factors are excluded here, as in ``_aa_dp``.

    Args:
        aa: The amino-acid CDR3, upper-case.
        template: CDR3 as nt integers, ``-1`` at free (insertion) positions.
        specs: Per position, ``None`` for germline or ``(kind, R, bias, is_first, is_last)``.
        floor: Prune partial paths at or below this weight. Every factor is a probability ``<= 1``,
            so a prefix can only lose weight — the prune is exact, not a heuristic.

    Returns:
        ``(p, nt)`` with ``p`` the path weight and ``nt`` the nucleotide string, or ``(0.0, None)``
        when no codon assignment survives. Ties break on the nt string, so a run is reproducible.
    """
    dp = {(-1, -1): (1.0, "")}                      # state = (nt[i-1], nt[i-2]) -> (weight, path)
    for i, fixed in enumerate(template):
        spec = specs[i]
        codon_end = i % 3 == 2
        aa_i = aa[i // 3]
        ndp: dict = {}
        for (p1, p2), (w, path) in dp.items():
            for nt in (range(4) if fixed < 0 else (fixed,)):
                ww = w
                if spec is not None:
                    kind, R, bias, is_first, is_last = spec
                    if kind == "L":
                        ww *= bias[nt] if is_first else R[nt, p1]
                    else:                           # "R" — DJ, read from the 3' end
                        if not is_first:
                            ww *= R[p1, nt]
                        if is_last:
                            ww *= bias[nt]
                if ww <= floor:
                    continue
                if codon_end and _CODON_AA[(p2, p1, nt)] != aa_i:
                    continue
                key = (nt, p1)
                cand = (ww, path + _NUM2NT[nt])
                if key not in ndp or cand > ndp[key]:
                    ndp[key] = cand
        if not ndp:
            return 0.0, None
        dp = ndp
    return max(dp.values())


def _mask(spec, default: list[str]) -> list[str]:
    """Resolve a ``v=``/``j=`` argument to the allele list to search.

    Three input modes, because real annotation tables have all three: ``None`` (nothing known),
    one allele name, or **several** — a list, or the comma-separated string an AIRR ``v_call``
    carries when the aligner could not choose (``"TRBV6-2*01,TRBV6-3*01"``). The multi case is not
    a pin: every listed allele is searched and the model picks the most plausible.
    """
    if spec is None:
        return default
    if isinstance(spec, str):
        spec = [s.strip() for s in spec.split(",") if s.strip()]
    return list(spec)


def _shortlist(opts: list, k: int) -> list:
    """Keep only the ``k`` best-scoring alleles (with all of their trims).

    NOTE: **A heuristic, and the only one here.** Searching every V × J pair costs ~23 s per TRB CDR3,
    so an unconstrained call needs a starting guess. Ranking is by the model's own
    ``P(gene)·P(del)`` on codon-compatible trims — not an invented alignment score — but a shortlist
    can still drop the true optimum. It applies *only* where the caller gave nothing; pass the calls
    the record already has and the search is exhaustive over them.
    """
    if k <= 0 or not opts:
        return opts
    best: dict[str, float] = {}
    for o in opts:
        best[o[0]] = max(best.get(o[0], 0.0), o[2])
    keep = {a for a, _ in sorted(best.items(), key=lambda kv: (-kv[1], kv[0]))[:k]}
    return [o for o in opts if o[0] in keep]


def _v_aa_options(prep: _Prepared, aa: str, N: int, v):
    """``[(V, len_v, P(V)·P(delV|V), gv)]`` — 3' V trims whose germline codons translate to ``aa``.

    The codon filter is the aa analogue of :func:`pgen._v_options`' prefix match and is what keeps
    this cheap: a V germline whose *full* codons mistranslate can never appear in any scenario, so
    the whole ``(V, len_v)`` branch dies before a single DP runs.
    """
    out = []
    for V in _mask(v, prep.functional_v):
        pv = prep.p_v.get(V, 0.0)
        if pv == 0.0:
            continue
        cut, maxp = prep.cut["v"][V], prep.maxpal["v_3"]
        for len_v in range(1, min(len(cut), N) + 1):     # V contributes >= 1 nt (matches OLGA)
            pd = prep.p_del["v"].get((V, len(cut) - len_v - maxp), 0.0)
            if pd == 0.0 or not _aa_germline_prefix_ok(cut, aa, len_v):
                continue
            out.append((V, len_v, pv * pd, [_NT2NUM[c] for c in cut[:len_v]]))
    return out


def _j_aa_options(prep: _Prepared, aa: str, N: int, j):
    """``[(J, len_j, P(delJ|J), gj)]`` — 5' J trims whose germline codons translate to ``aa``."""
    out = []
    for J in _mask(j, prep.functional_j):
        cut, maxp = prep.cut["j"][J], prep.maxpal["j_5"]
        for len_j in range(1, min(len(cut), N) + 1):     # J contributes >= 1 nt (matches OLGA)
            idxj = len(cut) - len_j
            pd = prep.p_del["j"].get((J, idxj - maxp), 0.0)
            if pd == 0.0:
                continue
            seq, left = cut[idxj:], N - len_j
            if not all(_CODON_TABLE[seq[3 * c - left:3 * c - left + 3]] == aa[c]
                       for c in range(-(-left // 3), N // 3)):
                continue
            out.append((J, len_j, pd, [_NT2NUM[c] for c in seq]))
    return out


def _vj_scenarios(prep, N, gv, gj, len_v, len_j):
    """The one VJ middle: the whole gap is a VJ insertion. Yields ``(weight, template, specs)``."""
    pins, R, bias = prep.p_ins["vj"], prep.R["vj"], prep.bias["vj"]
    ins = N - len_v - len_j
    if ins < 0 or ins >= len(pins) or pins[ins] == 0.0:
        return
    specs = [None] * N
    for p in range(len_v, len_v + ins):
        specs[p] = ("L", R, bias, p == len_v, p == len_v + ins - 1)
    yield pins[ins], gv + [-1] * ins + gj, specs


def _d_variants(prep, aa: str, N: int):
    """``{D: [(dc, P(delD|D), len, ok_positions)]}`` — every surviving 5'/3' trim of every D.

    Built once per call, because it depends on neither V nor J. ``ok_positions`` is the tuple of
    CDR3 offsets at which this trimmed D's *full* codons translate to ``aa`` — precomputed for the
    same reason: the check is a function of ``(dc, pos)`` alone, but the placement loop would
    otherwise re-run it once per (V trim, J trim) pair. On a 15-residue TRB CDR3 that is ~3e6
    calls against ~2e4, and it was the single largest cost in the search.

    Trims yielding the same ``dc`` are kept separately: :func:`pgen._d_aa_middle` sums them as
    distinct scenarios and this must range over the same set.
    """
    maxdl, maxdr = prep.maxpal["d_5"], prep.maxpal["d_3"]
    out = {}
    for D in prep.functional_d:
        cut = prep.cut["d"][D]
        vs = []
        for idx5 in range(len(cut) + 1):
            for idx3 in range(len(cut) - idx5 + 1):
                pdel = prep.p_del["d"].get((D, idx5 - maxdl, idx3 - maxdr), 0.0)
                if pdel:
                    dc = [_NT2NUM[c] for c in cut[idx5:len(cut) - idx3]]
                    ok = tuple(p for p in range(N - len(dc) + 1) if _d_codons_ok(dc, p, aa))
                    if ok:
                        vs.append((dc, pdel, len(dc), ok))
        out[D] = vs
    return out


def _d_placements(prep, N, J, len_v, len_j, dvars, pins_vd, pins_dj):
    """Yield ``(w_mid, lvd, ldj, pos, dc, ld, D)`` for every single-D middle — the max-product
    mirror of :func:`pgen._d_aa_middle`'s summands, ``P(n_D)·P(D|J)·P(delD|D)·P(insVD)·P(insDJ)``.

    No template is built here: the caller prunes on ``w_mid`` first and most placements never need
    one. ``P(D|J) == 0`` prunes the pair, which is where the genomic D->J locus order lives — see
    :func:`_best_vdj_middle`.
    """
    p1 = prep.p_nd.get(1, 0.0) + prep.p_nd.get(0, 0.0)   # 0-D folds into 1-D, as in pgen
    if p1 == 0.0:
        return
    n_vd, n_dj = len(pins_vd), len(pins_dj)
    right = N - len_j                                    # DJ insertion + J start here
    for D, variants in dvars.items():
        pdg = prep.p_d_given_j.get((J, D), 0.0)
        if pdg == 0.0:
            continue                                     # genomically impossible with this J
        w_d = p1 * pdg
        for dc, pdel, ld, ok in variants:
            w_dd = w_d * pdel
            hi = right - ld
            for pos in ok:
                if pos < len_v:
                    continue
                if pos > hi:
                    break                                # ok_positions is ascending
                lvd, ldj = pos - len_v, hi - pos
                if lvd >= n_vd or ldj >= n_dj:
                    continue
                w = w_dd * pins_vd[lvd] * pins_dj[ldj]
                if w > 0.0:
                    yield w, lvd, ldj, pos, dc, ld, D


def _d_template(prep, N, gv, gj, len_v, len_j, lvd, ldj, pos, dc, ld):
    """``(template, specs)`` for one single-D placement — the layout :func:`_aa_dp_max` consumes."""
    right = N - len_j
    specs = [None] * N
    R_vd, b_vd = prep.R["vd"], prep.bias["vd"]
    R_dj, b_dj = prep.R["dj"], prep.bias["dj"]
    for p in range(len_v, pos):
        specs[p] = ("L", R_vd, b_vd, p == len_v, p == pos - 1)
    for p in range(pos + ld, right):
        specs[p] = ("R", R_dj, b_dj, p == pos + ld, p == right - 1)
    return gv + [-1] * lvd + dc + [-1] * ldj + gj, specs


def _chain_bound(prep, key: str) -> list[float]:
    """``[max Markov weight of an n-nt insertion at ``key``]`` — excluding ``P(n)``, as ``_aa_dp``
    does, so it multiplies straight onto a scenario weight that already carries ``P(insLen)``.

    WARNING: The DP weight cannot be bounded by 1.0. The insertion chain costs ~0.4 per nucleotide, so on
    a 10 nt N-region the true weight is ~1e-6 and a bound of 1.0 leaves the branch-and-bound cutoff
    six orders of magnitude too high — nothing is ever pruned and TRB takes minutes per sequence.
    """
    R, bias = prep.R[key], prep.bias[key]
    rmax, bmax = float(R.max()), float(bias.max())
    return [1.0 if n == 0 else bmax * rmax ** (n - 1) for n in range(len(prep.p_ins[key]))]


def _gap_bound(prep, N: int, vdj: bool, dvars) -> list[float]:
    """``bound[g]`` — an upper bound on any middle filling a ``g`` nt gap, DP weight included.

    Keyed on the gap length because that is what pins down how much insertion chain has to be paid
    for; every other factor is maximised independently, which can only over-estimate (the safe
    direction — a loose bound costs time, never correctness).
    """
    def full(key):                                       # P(n) x the best chain, per length
        pins, ch = prep.p_ins[key], _chain_bound(prep, key)
        return [float(pins[n]) * ch[n] for n in range(len(ch))]

    if not vdj:
        vj = full("vj")
        return [vj[g] if g < len(vj) else 0.0 for g in range(N + 1)]
    vd, dj = full("vd"), full("dj")
    pair = [max((vd[a] * dj[n - a] for a in range(min(n, len(vd) - 1) + 1) if n - a < len(dj)),
                default=0.0) for n in range(N + 1)]
    dmax = [0.0] * (N + 1)
    for variants in dvars.values():
        for _dc, pdel, ld, _ok in variants:
            if ld <= N:
                dmax[ld] = max(dmax[ld], pdel)
    head = ((prep.p_nd.get(1, 0.0) + prep.p_nd.get(0, 0.0))
            * (max(prep.p_d_given_j.values()) if prep.p_d_given_j else 1.0))
    return [head * max((dmax[ld] * pair[g - ld] for ld in range(g + 1)), default=0.0)
            for g in range(N + 1)]


def _scenario_template(prep, N: int, sc: tuple):
    """``(template, specs)`` for one native scenario — the layout :func:`_aa_dp_max` consumes.

    The search returns scenarios, not nucleotide paths, because once the scenario is fixed the free
    positions are just two short insertion blocks and picking their codons is this one cheap DP.
    Carrying back-pointers through the whole sweep would cost far more than the handful of
    reconstructions a caller actually needs.
    """
    _w, V, len_v, J, len_j, D, idx5, idx3, pos = sc
    cutv, cutj = prep.cut["v"][V], prep.cut["j"][J]
    if len_v > len(cutv) or len_j > len(cutj) or len_v + len_j > N:
        return None
    gv = [_NT2NUM[c] for c in cutv[:len_v]]
    gj = [_NT2NUM[c] for c in cutj[len(cutj) - len_j:]]
    if D is None:
        got = list(_vj_scenarios(prep, N, gv, gj, len_v, len_j))
        return (got[0][1], got[0][2]) if got else None
    cutd = prep.cut["d"][D]
    dc = [_NT2NUM[c] for c in cutd[idx5:len(cutd) - idx3]]
    ld = len(dc)
    lvd, ldj = pos - len_v, (N - len_j) - (pos + ld)
    if lvd < 0 or ldj < 0:
        return None
    return _d_template(prep, N, gv, gj, len_v, len_j, lvd, ldj, pos, dc, ld)


def _infer_nt_native(model, prep, aa: str, N: int, v, j, n_best: int):
    """The production path: the native argmax DP for the scenarios, this module for the codons."""
    from . import native

    scen = []
    for vv in _mask(v, [None]):
        for jj in _mask(j, [None]):
            scen += native.best_aa_scenarios(model, aa, vv, jj, max(1, n_best))
    if not scen:
        return None
    best: dict[str, tuple] = {}
    for sc in sorted(scen, key=lambda s: -s[0]):
        built = _scenario_template(prep, N, sc)
        if built is None:
            continue
        p, nt = _aa_dp_max(aa, built[0], built[1])
        if nt is None:
            continue
        # sc[0] is already the full joint weight: the native DP folds in P(V)P(delV)P(J)P(delJ)
        # P(D|J)P(delD)P(insLen) *and* the same insertion Markov product this DP just re-walked.
        if sc[0] > best.get(nt, (0.0,))[0]:
            best[nt] = (sc[0], sc)
    if not best:
        return None

    def pin(spec, chosen):
        return spec if spec is None or (isinstance(spec, str) and "," not in spec) else chosen

    top = sorted(best.items(), key=lambda kv: (-kv[1][0], kv[0]))[:max(1, n_best)]
    rescored = sorted(((native.pgen_nt(model, nt, pin(v, s[1]), pin(j, s[3])), nt, s, w)
                       for nt, (w, s) in top), key=lambda t: (-t[0], t[1]))
    pgen, best_nt, s, w = rescored[0]
    _sw, V, len_v, J, len_j, D, idx5, idx3, pos = s
    ld = len(prep.cut["d"][D]) - idx5 - idx3 if D is not None else 0
    return Scenario(cdr3_nt=best_nt, v_call=V, j_call=J, v_end=len_v, j_start=N - len_j,
                    d_call=D if ld else None, d_start=pos if ld else None,
                    d_end=pos + ld if ld else None, pgen=pgen, scenario_p=w,
                    n_candidates=len(best),
                    runner_up_pgen=rescored[1][0] if len(rescored) > 1 else 0.0)


def infer_nt(model_or_prep, cdr3_aa: str, v=None, j=None, *,
             n_best: int = 8, keep: float = 1e-2, top_vj: int = 8) -> Scenario | None:
    """The most likely nucleotide CDR3 for an amino-acid CDR3, and its recombination markup.

    Two stages, and the second is what makes the first honest:

    1. **Max-product over (scenario, nucleotides) jointly.** Every scenario
       :func:`pgen.pgen_aa` sums over is enumerated — ``(V, delV) x (J, delJ) x (D, delD, position)``
       — and for each, :func:`_aa_dp_max` picks the single best codon assignment: germline positions
       are pinned to their segment, and each free N-region position takes the nucleotide maximising
       ``P(nt_1)·prod P(nt_k | nt_{k-1})`` under the VD / DJ / VJ dinucleotide model covering it.
    2. **The surviving candidates are re-scored with the exact marginal** :func:`pgen.pgen_nt`, and
       the winner is the best of those. Stage 1 maximises the *joint* ``P(nt, scenario)``, which is
       not the *marginal* ``P(nt)`` this function's contract asks for; the re-scoring is what turns
       a Viterbi path into an answer about the marginal, and it is why ``pgen`` on the result is a
       real ``pgen_nt`` rather than a path weight.

    **Both stages earn their place, measured against** :func:`infer_nt_bruteforce` **on generated
    productive draws with V and J pinned (25 TRG, 19 TRA, 10 residues each).**

    .. list-table::
       :header-rows: 1

       * - variant
         - TRG
         - TRA
       * - one scenario, best codon per residue (no search)
         - 9/25
         - 4/19
       * - the joint argmax alone, no re-score
         - 21/25
         - 15/19
       * - argmax + marginal re-score (the defaults)
         - 25/25
         - 19/19

    The cheap shortcut fails because a trim chosen before the codons pins a codon the true optimum
    would have trimmed away — the same unsoundness as pinning the germline flanks.

    **Two implementations.** Pass a :class:`~vdjtools.model.model.Model` and stage 1 runs natively
    (:func:`vdjtools.model.native.best_aa_scenarios`) — **2.5 ms per human TRB CDR3, 0.5 ms per
    TRA**, i.e. all 80k VDJdb records in about 3 minutes. Pass an already-:func:`pgen.prepare` -d
    model and it runs the pure-Python scenario search instead, which is the *reference*
    implementation the native one is validated against and is ~600x slower on TRB. ``keep`` and
    ``top_vj`` apply only to that reference path; the native DP needs neither, because it
    marginalizes over V and J at essentially no cost (0.26 vs 0.23 ms with both free).

    **Three input modes for the calls**, because annotation tables have all three:

    * ``v="TRBV5-1*01", j="TRBJ2-3*01"`` — the normal mode. Exhaustive over the trims of that pair.
    * ``v="TRBV6-2*01,TRBV6-3*01"`` (or a list) — the aligner offered several and could not choose.
      Every listed allele is searched and the model picks the most plausible.
    * ``v=None`` — nothing known; the DP marginalizes over every V (or J).

    Args:
        model_or_prep: A :class:`~vdjtools.model.model.Model` (native, the normal call) or a
            :func:`pgen.prepare` -d one (the pure-Python reference search).
        cdr3_aa: The CDR3 amino-acid sequence (conserved Cys -> conserved Phe/Trp inclusive).
        v: V allele, comma-separated alleles, a list, or ``None`` — see above. Allele-keyed, never
            gene-level: see the ``_gene_idx`` trap.
        j: J allele, comma-separated alleles, a list, or ``None``.
        n_best: How many distinct candidates to re-score with the exact marginal Pgen.
        keep: **Reference path only.** Candidates within this factor of the running best survive
            the prune. ``1.0`` is the pure argmax and the fastest.
        top_vj: **Reference path only.** Shortlist size when a call is ``None``, ranked by the
            model's own ``P(gene)·P(del)``; ``0`` searches everything. It exists because the
            Python search costs ~23 s per unconstrained TRB CDR3. The native path ignores it.

    Returns:
        A :class:`Scenario` with ``pgen`` the marginal of the winner and ``runner_up_pgen`` the
        marginal of the next distinct candidate — report :attr:`Scenario.margin`, because with a
        long non-templated core many sequences are near-equally likely. ``None`` when no nucleotide
        CDR3 encoding ``cdr3_aa`` is explicable under the model.

    Note:
        Tandem-D (``n_D = 2``) scenarios are not enumerated in stage 1. They cannot add candidates —
        a single D trimmed to zero length already reaches every middle — so they can only reorder
        them, and the stage-2 ``pgen_nt`` re-scoring includes the tandem-D mass in full. Raise
        ``n_best`` if a D-D model's ordering matters to you.
    """
    prep = _prep(model_or_prep)
    aa = cdr3_aa.upper()
    if not aa or any(a not in _AA2CODONS for a in aa):
        return None
    N = 3 * len(aa)
    if isinstance(model_or_prep, Model):
        return _infer_nt_native(model_or_prep, prep, aa, N, v, j, n_best)
    v_opts = _v_aa_options(prep, aa, N, v)
    j_opts = _j_aa_options(prep, aa, N, j)
    if v is None:                                        # nothing known -> start from a guess
        v_opts = _shortlist(v_opts, top_vj)
    if j is None:
        j_opts = _shortlist(j_opts, top_vj)
    if not v_opts or not j_opts:
        return None
    vdj = prep.chain_type == "VDJ"

    dvars = _d_variants(prep, aa, N) if vdj else {}
    gap_bound = _gap_bound(prep, N, vdj, dvars)
    # numpy scalar indexing costs ~10x a list index and this runs 10^5 times per sequence
    pins_vd = prep.p_ins["vd"].tolist() if vdj else None
    pins_dj = prep.p_ins["dj"].tolist() if vdj else None
    vd_ch = _chain_bound(prep, "vd") if vdj else None
    dj_ch = _chain_bound(prep, "dj") if vdj else None

    # Priority = the scenario's own weight × the best any middle filling its gap could reach. Both
    # the outer break and the per-placement skip are exact: every factor still to be paid is a
    # probability <= 1, so a partial product can only fall.
    outer = []
    for V, len_v, wv, gv in v_opts:
        for J, len_j, wj, gj in j_opts:
            g = N - len_v - len_j
            if g < 0:
                continue
            pj = prep.p_j.get(J if vdj else (V, J), 0.0)
            if pj and gap_bound[g]:
                outer.append((wv * pj * wj * gap_bound[g], wv * pj * wj,
                              V, J, len_v, len_j, gv, gj))
    outer.sort(key=lambda t: -t[0])

    # nt -> (joint weight, and the path that reached it — no need to re-derive the markup later)
    best: dict[str, tuple] = {}
    best_p, cut = 0.0, 0.0                               # cut = best_p * keep, the live cutoff
    for ub, w_out, V, J, len_v, len_j, gv, gj in outer:
        if ub <= cut:
            break                                        # and so does every remaining scenario
        if vdj:
            for w_mid, lvd, ldj, pos, dc, ld, D in _d_placements(
                    prep, N, J, len_v, len_j, dvars, pins_vd, pins_dj):
                w = w_out * w_mid
                if w * vd_ch[lvd] * dj_ch[ldj] <= cut:
                    continue                             # cannot beat the cutoff, skip the DP
                template, specs = _d_template(prep, N, gv, gj, len_v, len_j, lvd, ldj, pos, dc, ld)
                p, nt = _aa_dp_max(aa, template, specs, floor=cut / w)
                if nt is None:
                    continue
                p *= w
                if p > best.get(nt, (0.0,))[0]:
                    best[nt] = (p, V, J, len_v, len_j, D, pos, ld)
                if p > best_p:
                    best_p, cut = p, p * keep
        else:
            for w_mid, template, specs in _vj_scenarios(prep, N, gv, gj, len_v, len_j):
                w = w_out * w_mid
                if w <= cut:
                    continue
                p, nt = _aa_dp_max(aa, template, specs, floor=cut / w)
                if nt is None:
                    continue
                p *= w
                if p > best.get(nt, (0.0,))[0]:
                    best[nt] = (p, V, J, len_v, len_j, None, None, 0)
                if p > best_p:
                    best_p, cut = p, p * keep
    if not best:
        return None

    # Condition the re-score on the caller's pin when there is one; with several candidate calls
    # there is no single conditioning to marginalise over, so each candidate is scored under the
    # call its own path chose (P(nt, V, J) rather than P(nt) — reported as such).
    def pin(spec, chosen):
        return spec if spec is None or (isinstance(spec, str) and "," not in spec) else chosen

    # The native scorer is bit-identical to the Python one and ~290x faster on TRB, and this is
    # 3/4 of the runtime — so it is worth passing a Model rather than a pre-prepared one.
    if isinstance(model_or_prep, Model):
        from . import native
        score = lambda nt, vv, jj: native.pgen_nt(model_or_prep, nt, vv, jj)   # noqa: E731
    else:
        score = lambda nt, vv, jj: pgen_nt(prep, nt, vv, jj)                   # noqa: E731

    top = sorted(best.items(), key=lambda kv: (-kv[1][0], kv[0]))[:max(1, n_best)]
    rescored = sorted(((score(nt, pin(v, path[1]), pin(j, path[2])), nt, path)
                       for nt, path in top), key=lambda t: (-t[0], t[1]))
    pgen, best_nt, (sp, V, J, len_v, len_j, D, pos, ld) = rescored[0]
    return Scenario(cdr3_nt=best_nt, v_call=V, j_call=J, v_end=len_v, j_start=N - len_j,
                    d_call=D, d_start=pos if ld else None, d_end=pos + ld if ld else None,
                    pgen=pgen, scenario_p=sp, n_candidates=len(best),
                    runner_up_pgen=rescored[1][0] if len(rescored) > 1 else 0.0)
