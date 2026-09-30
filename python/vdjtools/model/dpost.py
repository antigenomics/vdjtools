"""Which D, and where — from an amino-acid junction alone.

Moved here from ``arda.dpost`` in vdjtools 4.8.0 (``antigenomics/vdjtools#183``,
``antigenomics/arda#144``). The division of labour is that arda owns the germline reference and
the markup against it — ``alleles.fasta``, ``cdr3_anchors.tsv``, V/D/J geometry against a *named*
germline — and vdjtools owns the recombination model and everything probabilistic about it. This
module marginalises the generative model's insertion-length and D-trimming distributions and
multiplies in ``P(D | J)``, which are model marginals, so it belongs on this side of that line.
It still reads arda for the germline and the markup, exactly as the rest of :mod:`vdjtools.model`
does.

A VDJdb-style record has no nucleotides, and D segments are short and trimmed at both
ends, so the D is often invisible in the translated junction. Two independent sources of
information remain, and they are complementary:

**Where.** The junction's nucleotide length is known (3x its amino-acid length), so
``insVD + |D surviving| + insDJ`` is *pinned*. Marginalising the generative model's
insertion-length and D-trimming distributions therefore places the D even when the
sequence says nothing at all about it. Measured against OLGA ground truth, the MAP
``d_start`` is a median 1 nt off for mouse TRB, 2 nt for human TRB, and 3 nt for TRD and
IGH.

**Which.** The length constraint is nearly useless for identity — the D length
distributions overlap, so the posterior barely moves off the prior. What identity the
prior does carry is ``P(D | J)``, and for TRB that is mostly genomic order: TRBD2 lies 3'
of the whole TRBJ1 cluster, so a TRBJ1 junction can only have used TRBD1
(:func:`~vdjtools.model.reference.forbidden_dj_pairs`). What otherwise identifies a D is the
amino-acid match, and only where enough D survives: median surviving D is 17 nt for IGH
(~5.7 aa) but 5 nt for human TRB (~1.7 aa).

So neither source alone is enough, and which one dominates flips by locus:

    locus       prior only   aa only   combined     n     (held-out seed, generated)
    human IGH      15 %       81 %       82 %      345
    human TRB      76 %       70 %       82 %      595
    human TRD      86 %       88 %       87 %      699
    mouse TRB      76 %       83 %       85 %      699

"prior only" is ``beta = 0``; "aa only" is argmax of the match score under a uniform prior,
ties broken by marginal usage. Combining wins at IGH and both TRB; TRD is a wash, because
one D gene (TRDD3) accounts for 85 % of rearrangements and the aa match already finds it.

The combination is ``log P(D | M, J) + beta * s_D``: the length-and-J prior, tempered by
the best gapless local alignment score ``s_D`` of the D's three-frame translations against
the non-templated middle of the junction. ``beta`` is fitted per locus and shipped in arda's
``database/vdj/<organism>/d_prior.tsv`` with the distributions themselves (see
:func:`_shipped_prior` for why it is read there rather than vendored), so nothing here needs OLGA
at runtime. It is flat above ~1.25 for TRB, so the shipped values are not delicate.

**Honesty about the numbers.** The table is measured on junctions drawn from the same
generative model that supplies the prior, so the prior's contribution is flattered. The
amino-acid contribution is not — it is germline matching. Rearrangements that genomic order
forbids are excluded from the truth: OLGA's human TRB model emits TRBD2 x TRBJ1 in 8.7 % of
draws, and scoring against those measures agreement with a model artifact.

Out of model, against nucleotide D calls (E <= 0.05) on the real GenBank fixtures: human
TRB 94 %, IGH 85 %, TRD 91 %, mouse TRB 85 %. On TRB, note that both this posterior and the
nucleotide caller enforce the same D2-x-J1 constraint, so their agreement on TRBJ1 records
is guaranteed rather than earned; the TRBJ2 rows, where both D genes stay possible, score
91 % (human) and 81 % (mouse).

Priors exist only for the (organism, locus) pairs with a published model: human IGH, TRB
and TRD, and mouse TRB. Everything else returns ``None`` rather than guessing.

**Call :func:`posterior_d_batch`, never :func:`posterior_d` in a loop.** The per-record form
marks up one junction per call, and markup is most of the cost; the batch form hands the whole
column to arda's batched markup and loads the prior and the D translations once.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

__all__ = ["DPosterior", "DPrior", "posterior_d", "posterior_d_batch", "load_d_prior"]


def _gene(allele: str) -> str:
    return allele.split("*")[0]


@dataclass
class DPrior:
    """The shipped generative-model summaries for one locus."""

    ins_vd: list[float]
    ins_dj: list[float]
    dlen: dict[str, list[float]]                 # allele -> P(surviving nt length)
    d_given_j: dict[str, dict[str, float]]       # j allele -> {d allele: P}
    d_marginal: dict[str, float]
    beta: float

    def __post_init__(self) -> None:
        """Precompute the part of the marginalisation that does not depend on the record.

        The per-record joint is ``P(D) * insVD[a] * sum_L dlen[D][L] * insDJ[n - a - L]``, summed
        over the surviving-D length ``L``, where ``n`` is the pinned nucleotide middle. That inner
        sum is the CONVOLUTION ``dlen[D] * insDJ`` evaluated at ``n - a``, and it is the same for
        every junction -- so it is computed once per prior here, and a record costs one vectorised
        multiply per D allele instead of a Python loop over every (insVD, dlen) pair.
        """
        self._ins_vd = np.asarray(self.ins_vd, dtype=float)
        dj = np.asarray(self.ins_dj, dtype=float)
        self._tail = {a: np.convolve(np.asarray(dl, dtype=float), dj)
                      for a, dl in self.dlen.items() if dl and dj.size}


@dataclass
class DPosterior:
    """Posterior over the D gene, and over where it sits in the junction."""

    locus: str
    d_call: str                                  # MAP D gene
    posterior: float                             # P(MAP gene)
    entropy: float                               # bits, over genes
    by_gene: dict[str, float] = field(default_factory=dict)
    support_aa: int = 0                          # best aa local-align score for the MAP gene
    d_start: int = -1                            # MAP nt offset of D start within the junction
    d_start_ci90: tuple[int, int] = (-1, -1)     # narrowest 90 % credible interval
    n_middle_nt: int = 0                         # insVD + |D| + insDJ, pinned by the length

    @property
    def confident(self) -> bool:
        """A hard call. 0.9 keeps ~the top decile of TRB records and most of IGH."""
        return self.posterior >= 0.9


def _shipped_prior(organism: str) -> Path | None:
    """``database/vdj/<organism>/d_prior.tsv`` from the installed arda, or ``None``.

    ⚠ The table is still read from arda's reference tree, and deliberately not vendored here. It is
    a fitted model artifact and ``antigenomics/arda#144`` is right that it does not belong in a
    germline reference -- but ``arda.hmm`` and ``arda.scenarios`` are parameterised by the same
    file, so a copy here would be a SECOND copy of a fitted artifact, and the two would drift. One
    file, read by whoever needs it, until those two modules follow the posterior across (which is
    `antigenomics/arda#145`); then it moves into the model bundle beside them.
    """
    from arda.paths import vdj_dir

    p = Path(vdj_dir(organism)) / "d_prior.tsv"
    return p if p.exists() else None


@lru_cache(maxsize=8)
def load_d_prior(organism: str, path: Path | None = None) -> dict[str, DPrior]:
    """``{locus: DPrior}``; empty when the organism has no shipped table.

    The table is parsed by :func:`arda.scenarios.load_prior_table`, which sits next to the code that
    WRITES the format, so the layout has exactly one parser. What is added here is the precompute
    the posterior needs (:meth:`DPrior.__post_init__`).

    ``path`` reads a table fitted elsewhere -- ``arda scenarios`` writes this layout -- instead of
    the shipped one. Never: *using* an estimate is not *adopting* it. A caller-supplied ``path``
    that does not exist RAISES; the shipped one is allowed to be absent, since most
    (organism, D-locus) pairs have no table and :func:`posterior_d` returns ``None`` for them.
    """
    from arda.scenarios import load_prior_table

    if path is None and _shipped_prior(organism) is None:
        return {}
    return {locus: DPrior(ins_vd=t.ins_vd, ins_dj=t.ins_dj, dlen=t.dlen,
                          d_given_j=t.d_given_j, d_marginal=t.d_marginal, beta=t.beta)
            for locus, t in load_prior_table(organism, path).items()}


@lru_cache(maxsize=16)
def _d_aa_frames(locus: str, organism: str) -> dict[str, tuple[str, str, str]]:
    """``{D allele: (aa frame 0, 1, 2)}`` from arda's shipped IMGT D germlines."""
    import polars as pl
    from arda.refbuild.translate import translate

    from .reference import load_germline

    d = load_germline(locus, organism).filter(pl.col("segment") == "D")
    return {allele: tuple(translate(seq[f:], 0) for f in (0, 1, 2))
            for allele, seq in zip(d["allele"], d["sequence"])}


def _mask_forbidden(pd: dict[str, float], j_call: str, locus: str) -> dict[str, float]:
    """Zero the D alleles lying 3' of the J, then renormalise.

    Delegates the genomic-order fact itself to
    :func:`~vdjtools.model.reference.forbidden_dj_pairs`, which is where this repository keeps it
    -- TRBD2 sits downstream of the whole TRBJ1 cluster, so no deletional join can reach it. An
    ambiguous J spanning both clusters forbids nothing.
    """
    from .reference import forbidden_dj_pairs

    js = [a.strip() for a in j_call.split(",") if a.strip()]
    if not js or not pd:
        return pd
    banned = forbidden_dj_pairs(list(pd), js, locus)
    # Only a D forbidden against EVERY named J is excluded: an ambiguous call spanning both
    # clusters leaves the join possible.
    dead = {d for d in pd if all((d, j) in banned for j in js)}
    if not dead:
        return pd
    kept = {a: p for a, p in pd.items() if a not in dead}
    total = sum(kept.values())
    return {a: p / total for a, p in kept.items()} if total > 0 else pd


def _logsumexp(xs: list[float]) -> float:
    if not xs:
        return -math.inf
    m = max(xs)
    if m == -math.inf:
        return -math.inf
    return m + math.log(sum(math.exp(x - m) for x in xs))


def _score(prior: DPrior, locus: str, organism: str, j_call: str,
           middle: str, v_end: int) -> DPosterior | None:
    """The posterior itself, given a junction already marked up. No I/O, no markup."""
    from arda import _markup

    n_middle = 3 * len(middle)                   # insVD + |D surviving| + insDJ, pinned
    alleles = sorted(set(prior.dlen) | set(prior.d_marginal))
    if not alleles:
        return None

    # P(D allele | J): the load-bearing prior for TRB, where genomic order forbids TRBD2 x
    # TRBJ1 outright. Back off to the marginal when the J allele is outside the model -- but
    # the marginal pools both J clusters, so re-apply the same mask (the shipped human model
    # has no TRBJ1-6*01 row at all, and would otherwise let TRBD2 back in through the door).
    j_id = j_call.split(",")[0].strip()
    pd_given_j = prior.d_given_j.get(j_id) or _mask_forbidden(prior.d_marginal, j_call, locus)

    # Joint over (D, insVD), marginalising the surviving-D length and insDJ. `prior._tail[allele]`
    # is `dlen[allele]` convolved with `insDJ`, so the whole inner marginalisation is one lookup
    # per insVD bin: `tail[n_middle - a]`, which is the array reversed and sliced.
    ins_vd = prior._ins_vd
    n_a = ins_vd.size
    joint: dict[str, np.ndarray] = {}
    pa = np.zeros(max(1, n_a))
    for allele in alleles:
        base = pd_given_j.get(allele, 0.0)
        tail = prior._tail.get(allele)
        if base <= 0 or tail is None:
            continue
        # tail[n_middle - a] for a in 0..n_a-1, zero where that index leaves the support.
        idx = n_middle - np.arange(n_a)
        vals = np.where((idx >= 0) & (idx < tail.size), tail[np.clip(idx, 0, tail.size - 1)], 0.0)
        row = base * ins_vd * vals
        if row.any():
            joint[allele] = row
            pa[:n_a] += row
    total = float(pa.sum())
    if total <= 0:
        return None

    # Amino-acid evidence: best gapless local alignment of each D's three-frame
    # translations against the non-templated middle.
    frames = _d_aa_frames(locus, organism)
    score: dict[str, int] = {}
    for allele in joint:
        fr = frames.get(allele)
        score[allele] = 0 if not fr or not middle else max(
            (_markup.d_local_align(middle, f)[0] for f in fr if f), default=0)

    # Combine: log-prior + beta * aa score, marginalised over alleles within a gene.
    log_by_allele = {a: math.log(joint[a].sum()) + prior.beta * score[a] for a in joint}
    by_gene_log: dict[str, list[float]] = {}
    for allele, lp in log_by_allele.items():
        by_gene_log.setdefault(_gene(allele), []).append(lp)
    gene_log = {g: _logsumexp(v) for g, v in by_gene_log.items()}
    z = _logsumexp(list(gene_log.values()))
    by_gene = {g: math.exp(lp - z) for g, lp in gene_log.items()}
    best_gene = max(by_gene, key=by_gene.get)
    # `+ 0.0` so a degenerate posterior reports 0.0 rather than -0.0.
    entropy = -sum(p * math.log2(p) for p in by_gene.values() if p > 0) + 0.0

    # Where: d_start = (nt templated by V) + insVD, marginalising D and its trimming.
    pa_norm = pa / total
    order = sorted(range(pa_norm.size), key=lambda i: pa_norm[i], reverse=True)
    cum, chosen = 0.0, []
    for i in order:
        chosen.append(i)
        cum += float(pa_norm[i])
        if cum >= 0.90:
            break
    v_nt = 3 * v_end
    map_a = order[0] if order else -1
    ci = (v_nt + min(chosen), v_nt + max(chosen)) if chosen else (-1, -1)

    best_allele = max((a for a in log_by_allele if _gene(a) == best_gene),
                      key=lambda a: log_by_allele[a])
    return DPosterior(
        locus=locus, d_call=best_gene, posterior=by_gene[best_gene], entropy=entropy,
        by_gene=by_gene, support_aa=score.get(best_allele, 0),
        d_start=(v_nt + map_a) if map_a >= 0 else -1, d_start_ci90=ci,
        n_middle_nt=n_middle)


def posterior_d(junction_aa: str, v_call: str, j_call: str,
                species: str = "human", prior_path: Path | None = None) -> DPosterior | None:
    """Posterior over the D gene (and its position) for one amino-acid junction.

    ``junction_aa`` is junction space (Cys104 .. Phe/Trp118, both included), as in
    :mod:`arda.cdr3fix`. Returns ``None`` when the locus has no D, no shipped model, or
    the junction cannot be marked up.

    ⚠ For more than a handful of records call :func:`posterior_d_batch`: this form marks up one
    junction per call, and the markup is most of the cost.
    """
    from arda.cdr3fix import load_anchors, markup_cdr3, resolve_locus, resolve_species

    organism = resolve_species(species)
    locus = resolve_locus(v_call, j_call)
    prior = load_d_prior(organism, prior_path).get(locus)
    if prior is None:
        return None
    mk = markup_cdr3(junction_aa, v_call, j_call, organism, anchors=load_anchors(organism))
    if mk.v_end < 0 or mk.j_start < 0 or mk.j_start < mk.v_end:
        return None
    return _score(prior, locus, organism, j_call,
                  mk.cdr3_repaired[mk.v_end:mk.j_start], mk.v_end)


def posterior_d_batch(junction_aas: Sequence[str], v_calls: Sequence[str],
                      j_calls: Sequence[str], species: str | Iterable[str] = "human",
                      prior_path: Path | None = None) -> list[DPosterior | None]:
    """:func:`posterior_d` over many records: **one row out per row in, in input order**.

    ``None`` where the record has no D locus, no shipped prior, or a junction that cannot be
    marked up -- never an exception, because a real corpus legitimately contains species and loci
    with no model. ``species`` is one name for every row, or one per row.

    This exists because the per-record form has no plural and a caller with N junctions therefore
    wrote a Python row loop, which on the VDJdb build was the single largest stage of the whole
    build -- 15.96 s over 119,034 keys, 35.8 % of it (``antigenomics/arda#142``). What that loop
    paid for was the markup, once per record; here arda's batched ``markup_records`` gets the whole
    column and the prior and the D translations are loaded once.
    """
    from arda.cdr3fix import markup_records, resolve_species

    import polars as pl

    n = len(junction_aas)
    if not (len(v_calls) == len(j_calls) == n):
        raise ValueError(f"ragged input: {n} junctions, {len(v_calls)} V, {len(j_calls)} J")
    sp = [species] * n if isinstance(species, str) else list(species)
    if len(sp) != n:
        raise ValueError(f"ragged input: {n} junctions, {len(sp)} species")
    if not n:
        return []

    df = pl.DataFrame({"cdr3": list(junction_aas), "v": list(v_calls), "j": list(j_calls),
                       "species": sp})
    out: list[DPosterior | None] = [None] * n
    # One markup call per ORGANISM, because `markup_records` loads one anchor table per call.
    for organism in {resolve_species(s) for s in sp}:
        idx = [i for i, s in enumerate(sp) if resolve_species(s) == organism]
        priors = load_d_prior(organism, prior_path)
        if not priors:
            continue
        recs = markup_records(df[idx], cdr3="cdr3", v="v", j="j", organism=organism)
        for i, mk in zip(idx, recs):
            prior = priors.get(mk.locus)
            if prior is None or mk.v_end < 0 or mk.j_start < 0 or mk.j_start < mk.v_end:
                continue
            out[i] = _score(prior, mk.locus, organism, j_calls[i],
                            mk.cdr3_repaired[mk.v_end:mk.j_start], mk.v_end)
    return out
