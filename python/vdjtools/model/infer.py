"""EM inference of model marginals from nucleotide CDR3 sequences.

Expectation–Maximization over the recombination Bayes net: the **E-step** enumerates every
scenario that could produce each observed nt CDR3 (the same enumeration as ``pgen``), weights
them by the current model, and accumulates *soft counts* per event realization; the **M-step**
re-normalizes those counts in polars to get the next marginals. Trained on out-of-frame reads,
it recovers the raw generation model (no productivity conditioning, so no selection bias).

Closed-loop oracle: generate synthetic sequences from a known model, then ``infer`` must recover
that model's marginals (see the tests). This is the reference driver; the E-step hot loop is a
Phase 1f native-port candidate.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from math import log

import numpy as np
import polars as pl

from .model import Model
from .pgen import _NT2NUM, _common_prefix, _common_suffix, _j_candidates, _v_candidates, prepare
from .schema import normalization_keys, table_columns


@dataclass(slots=True)
class InferenceReport:
    """Per-iteration diagnostics from :func:`infer`."""

    loglik: list[float] = field(default_factory=list)  # mean per-sequence log-Pgen (over scoreable reads)
    n_scoreable: list[int] = field(default_factory=list)
    gene_tv: list[float] = field(default_factory=list)  # V-usage total-variation vs previous iter
    n_iter: int = 0
    converged: bool = False


def _tv(a: pl.DataFrame, b: pl.DataFrame, keycols: list[str]) -> float:
    """Total-variation distance between two tables' ``p``, aligned by ``keycols``."""
    da = {tuple(r[:-1]): r[-1] for r in a.select([*keycols, "p"]).iter_rows()}
    db = {tuple(r[:-1]): r[-1] for r in b.select([*keycols, "p"]).iter_rows()}
    return 0.5 * sum(abs(da.get(k, 0.0) - db.get(k, 0.0)) for k in set(da) | set(db))


def _insert_markov(seq: str, R: np.ndarray, bias: np.ndarray, *, from_right: bool):
    """(markov weight, [(from_nt, to_nt), ...]) for an N-region (excludes the P(len) factor)."""
    n = len(seq)
    if n == 0:
        return 1.0, ()
    nums = [_NT2NUM[c] for c in seq]
    trans = []
    if from_right:
        w = bias[nums[-1]]
        for k in range(n - 2, -1, -1):
            w *= R[nums[k], nums[k + 1]]
            trans.append((nums[k + 1], nums[k]))
    else:
        w = bias[nums[0]]
        for k in range(1, n):
            w *= R[nums[k], nums[k - 1]]
            trans.append((nums[k - 1], nums[k]))
    return w, trans


def _estep_seq(prep, s: str, counts: dict, mask=None, allow_dd: bool = True) -> float:
    """Accumulate one sequence's soft counts into ``counts``; return its Pgen (for log-lik).

    ``mask`` optionally restricts enumeration to a read's aligned genes — ``(v_genes, j_genes,
    d_genes)`` name lists (e.g. from arda). This is what makes VDJ inference tractable: without it
    every V that shares the conserved Cys prefix is a candidate and the D enumeration runs for each.
    ``allow_dd`` gates the tandem (n_D=2) enumeration for this read (arda-anchored D-D learning).
    """
    N = len(s)
    vdj = prep.chain_type == "VDJ"
    local: dict = defaultdict(float)
    total = 0.0

    v_mask = mask[0] if mask else prep.functional_v
    j_mask = mask[1] if mask else prep.functional_j
    d_mask = mask[2] if mask and mask[2] else (prep.functional_d if vdj else None)
    v_cands = _v_candidates(prep, s, v_mask)
    j_cands = _j_candidates(prep, s, j_mask)
    for V, pv, v_opts in v_cands:
        for J, j_opts in j_cands:
            pj = prep.p_j.get(J if vdj else (V, J), 0.0)
            if pj == 0.0:
                continue
            for len_v, nv, p_dv in v_opts:
                for len_j, nj, p_dj in j_opts:
                    if len_v + len_j > N:
                        continue
                    mid = s[len_v:N - len_j]
                    base = pv * pj * p_dv * p_dj
                    if vdj:
                        p1 = prep.p_nd.get(0, 0.0) + prep.p_nd.get(1, 0.0)
                        if p1 > 0.0:
                            total += _accum_vdj(prep, J, V, nv, nj, mid, base * p1, local, d_mask)
                        p2 = prep.p_nd.get(2, 0.0)
                        if p2 > 0.0 and prep.p_d2_given_d1 and allow_dd:
                            total += _accum_dd(prep, J, V, nv, nj, mid, base * p2, local, d_mask)
                    else:
                        total += _accum_vj(prep, V, J, nv, nj, mid, base, local)

    if total > 0.0:
        inv = 1.0 / total
        for (event, key), w in local.items():
            counts[event][key] += w * inv
    return total


def _accum_vj(prep, V, J, nv, nj, mid, base, local) -> float:
    L = len(mid)
    pins = prep.p_ins["vj"]
    if L >= len(pins) or pins[L] == 0.0:
        return 0.0
    mw, trans = _insert_markov(mid, prep.R["vj"], prep.bias["vj"], from_right=False)
    w = base * pins[L] * mw
    if w <= 0.0:
        return 0.0
    local[("v_choice", (V,))] += w
    local[("j_choice", (V, J))] += w
    local[("v_3_del", (V, nv))] += w
    local[("j_5_del", (J, nj))] += w
    local[("vj_ins", (L,))] += w
    for fr, to in trans:
        local[("vj_dinucl", (fr, to))] += w
    return w


def _accum_vdj(prep, J, V, nv, nj, mid, base, local, d_mask=None) -> float:
    pins_vd, pins_dj = prep.p_ins["vd"], prep.p_ins["dj"]
    maxdl, maxdr = prep.maxpal["d_5"], prep.maxpal["d_3"]
    m = len(mid)
    seq_total = 0.0
    for D in (d_mask if d_mask is not None else prep.functional_d):
        pdj = prep.p_d_given_j.get((J, D), 0.0)
        if pdj == 0.0:
            continue
        cut = prep.cut["d"][D]
        for idx5 in range(len(cut) + 1):
            for idx3 in range(len(cut) - idx5 + 1):
                pdel = prep.p_del["d"].get((D, idx5 - maxdl, idx3 - maxdr), 0.0)
                if pdel == 0.0:
                    continue
                dc = cut[idx5:len(cut) - idx3]
                ld = len(dc)
                for pos in range(0, m - ld + 1):
                    if mid[pos:pos + ld] != dc:
                        continue
                    ins_vd, ins_dj = mid[:pos], mid[pos + ld:]
                    lvd, ldj = len(ins_vd), len(ins_dj)
                    if lvd >= len(pins_vd) or pins_vd[lvd] == 0.0:
                        continue
                    if ldj >= len(pins_dj) or pins_dj[ldj] == 0.0:
                        continue
                    mwvd, tvd = _insert_markov(ins_vd, prep.R["vd"], prep.bias["vd"], from_right=False)
                    mwdj, tdj = _insert_markov(ins_dj, prep.R["dj"], prep.bias["dj"], from_right=True)
                    w = base * pdj * pdel * pins_vd[lvd] * mwvd * pins_dj[ldj] * mwdj
                    if w <= 0.0:
                        continue
                    seq_total += w
                    local[("v_choice", (V,))] += w
                    local[("j_choice", (J,))] += w
                    local[("d_gene", (J, D))] += w
                    local[("v_3_del", (V, nv))] += w
                    local[("j_5_del", (J, nj))] += w
                    local[("d_del", (D, idx5 - maxdl, idx3 - maxdr))] += w
                    local[("vd_ins", (lvd,))] += w
                    local[("dj_ins", (ldj,))] += w
                    local[("n_d", (1,))] += w  # this scenario has exactly one D
                    for fr, to in tvd:
                        local[("vd_dinucl", (fr, to))] += w
                    for fr, to in tdj:
                        local[("dj_dinucl", (fr, to))] += w
    return seq_total


def _accum_dd(prep, J, V, nv, nj, mid, base, local, d_mask=None) -> float:
    """Two-D (tandem) soft-count accumulation: ``mid`` = [insVD] D1 [insDD] D2 [insDJ], each D ≥1 nt.

    The reference (naive) enumeration mirroring :func:`~vdjtools.model.pgen._dd_middle`; correctness
    over speed (validated closed-loop on small synthetic models — the native E-step is the fast path).
    """
    maxdl, maxdr = prep.maxpal["d_5"], prep.maxpal["d_3"]
    pins_vd, pins_dd, pins_dj = prep.p_ins["vd"], prep.p_ins["dd"], prep.p_ins["dj"]
    m = len(mid)
    dset = d_mask if d_mask is not None else prep.functional_d
    seq_total = 0.0
    for D1 in dset:
        pd1 = prep.p_d_given_j.get((J, D1), 0.0)
        if pd1 == 0.0:
            continue
        cut1 = prep.cut["d"][D1]
        for i5 in range(len(cut1) + 1):
            for i3 in range(len(cut1) - i5 + 1):
                pdel1 = prep.p_del["d"].get((D1, i5 - maxdl, i3 - maxdr), 0.0)
                if pdel1 == 0.0:
                    continue
                dc1 = cut1[i5:len(cut1) - i3]
                ld1 = len(dc1)
                if ld1 < 1:
                    continue
                for pos1 in range(0, m - ld1):  # leave ≥1 nt for D2
                    if mid[pos1:pos1 + ld1] != dc1:
                        continue
                    lvd = pos1
                    if lvd >= len(pins_vd) or pins_vd[lvd] == 0.0:
                        continue
                    mwvd, tvd = _insert_markov(mid[:pos1], prep.R["vd"], prep.bias["vd"], from_right=False)
                    left = base * pd1 * pdel1 * pins_vd[lvd] * mwvd
                    if left <= 0.0:
                        continue
                    for D2 in dset:
                        pd2g = prep.p_d2_given_d1.get((D1, D2), 0.0)
                        if pd2g == 0.0:
                            continue
                        cut2 = prep.cut["d"][D2]
                        for k5 in range(len(cut2) + 1):
                            for k3 in range(len(cut2) - k5 + 1):
                                pdel2 = prep.p_del_d2.get((D2, k5 - maxdl, k3 - maxdr), 0.0)
                                if pdel2 == 0.0:
                                    continue
                                dc2 = cut2[k5:len(cut2) - k3]
                                ld2 = len(dc2)
                                if ld2 < 1:
                                    continue
                                for pos2 in range(pos1 + ld1, m - ld2 + 1):
                                    if mid[pos2:pos2 + ld2] != dc2:
                                        continue
                                    ins_dd, ins_dj = mid[pos1 + ld1:pos2], mid[pos2 + ld2:]
                                    ldd, ldj = len(ins_dd), len(ins_dj)
                                    if ldd >= len(pins_dd) or pins_dd[ldd] == 0.0:
                                        continue
                                    if ldj >= len(pins_dj) or pins_dj[ldj] == 0.0:
                                        continue
                                    mwdd, tdd = _insert_markov(ins_dd, prep.R["dd"], prep.bias["dd"], from_right=False)
                                    mwdj, tdj = _insert_markov(ins_dj, prep.R["dj"], prep.bias["dj"], from_right=True)
                                    w = left * pd2g * pdel2 * pins_dd[ldd] * mwdd * pins_dj[ldj] * mwdj
                                    if w <= 0.0:
                                        continue
                                    seq_total += w
                                    local[("v_choice", (V,))] += w
                                    local[("j_choice", (J,))] += w
                                    local[("d_gene", (J, D1))] += w
                                    local[("d2_gene", (D1, D2))] += w
                                    local[("v_3_del", (V, nv))] += w
                                    local[("j_5_del", (J, nj))] += w
                                    local[("d_del", (D1, i5 - maxdl, i3 - maxdr))] += w
                                    local[("d2_del", (D2, k5 - maxdl, k3 - maxdr))] += w
                                    local[("vd_ins", (lvd,))] += w
                                    local[("dd_ins", (ldd,))] += w
                                    local[("dj_ins", (ldj,))] += w
                                    local[("n_d", (2,))] += w  # this scenario has two Ds
                                    for fr, to in tvd:
                                        local[("vd_dinucl", (fr, to))] += w
                                    for fr, to in tdd:
                                        local[("dd_dinucl", (fr, to))] += w
                                    for fr, to in tdj:
                                        local[("dj_dinucl", (fr, to))] += w
    return seq_total


# Events re-estimated by EM (the germline stays fixed from the template). n_d is fit too: for a
# single-D model the E-step only ever emits n_D=1 counts, so it renormalizes to δ(1) (a no-op);
# for a D-D model it learns P(n_D=2) from the tandem soft counts.
def _fit_events(manifest) -> list[str]:
    return list(manifest.events)


def _mstep(template: Model, counts: dict, nd_prior: float = 0.0) -> dict[str, pl.DataFrame]:
    if nd_prior and counts.get("n_d"):  # Dirichlet pseudocount on the single-D (n_D=1) bucket
        counts["n_d"][(1,)] += nd_prior
    tables = {}
    for name, event in template.manifest.events.items():
        if not counts.get(name):  # event never accumulated (e.g. n_d on a VJ locus) — keep template
            tables[name] = template.tables[name]
            continue
        cols = list(table_columns(event))  # value cols..., "p"
        value_cols = cols[:-1]
        data: dict = {c: [] for c in cols}
        for key, cnt in counts[name].items():
            for c, val in zip(value_cols, key):
                data[c].append(val)
            data["p"].append(cnt)
        df = pl.DataFrame(data, schema=table_columns(event))
        keys = normalization_keys(event)
        total = pl.col("p").sum().over(keys) if keys else pl.col("p").sum()
        tables[name] = df.with_columns(p=pl.col("p") / total)
    return tables


def _uniform_init(template: Model) -> dict[str, pl.DataFrame]:
    tables = {}
    for name, event in template.manifest.events.items():
        df = template.tables[name].with_columns(p=pl.lit(1.0))
        keys = normalization_keys(event)
        total = pl.col("p").sum().over(keys) if keys else pl.col("p").sum()
        tables[name] = df.with_columns(p=pl.col("p") / total)
    return tables


def _align_init(template: Model, sequences: list[str]) -> dict[str, pl.DataFrame]:
    """Seed gene usage from a best-alignment vote (each read votes its longest-matching V/J).

    Everything else starts uniform. Without this, uniform gene usage makes the E-step enumerate
    every gene (they all share the conserved Cys/Phe anchors) — this concentrates it so pruning
    bites; EM then refines. This is the alignment step a real pipeline (arda) does up front.
    """
    prep = prepare(template)
    vdj = template.chain_type == "VDJ"
    v_votes: dict[str, float] = defaultdict(float)
    j_votes: dict[str, float] = defaultdict(float)
    vj_votes: dict[tuple, float] = defaultdict(float)
    for s in sequences:
        s = s.upper()
        bv = max(prep.functional_v, key=lambda v: _common_prefix(prep.cut["v"][v], s))
        bj = max(prep.functional_j, key=lambda j: _common_suffix(prep.cut["j"][j], s))
        v_votes[bv] += 1.0
        j_votes[bj] += 1.0
        vj_votes[(bv, bj)] += 1.0

    events = template.manifest.events
    tables = _uniform_init(template)
    tables["v_choice"] = _set_p(template.tables["v_choice"], "v_allele", v_votes, normalization_keys(events["v_choice"]))
    if vdj:
        tables["j_choice"] = _set_p(template.tables["j_choice"], "j_allele", j_votes, normalization_keys(events["j_choice"]))
    else:
        tables["j_choice"] = _set_p(template.tables["j_choice"], ("v_allele", "j_allele"), vj_votes, normalization_keys(events["j_choice"]))
    return tables


def _set_p(df: pl.DataFrame, key, votes: dict, norm_keys: list[str]) -> pl.DataFrame:
    """Overwrite ``p`` from a votes dict (keyed by ``key`` col(s)), normalized within ``norm_keys``."""
    if isinstance(key, str):
        p = [votes.get(k, 0.0) for k in df[key]]
    else:
        p = [votes.get(tuple(row), 0.0) for row in df.select(list(key)).iter_rows()]
    out = df.with_columns(p=pl.Series("p", p))
    total = pl.col("p").sum().over(norm_keys) if norm_keys else pl.col("p").sum()
    # groups with no votes stay 0 (undefined conditional — allowed); avoid 0/0.
    return out.with_columns(p=pl.when(total > 0).then(pl.col("p") / total).otherwise(0.0))


# The D-bearing loci where tandem-D (D-D) recombination is biologically documented and modelled by
# default. VJ loci (TRA/TRG/IGK/IGL) have no D and are always single-chain-D-free.
DD_DEFAULT_LOCI = frozenset({"TRB", "TRD", "IGH"})


def _maybe_promote_dd(template: Model, single_d: bool, p_nd2_init: float) -> Model:
    """Promote a single-D VDJ template to a tandem-D model by default (unless ``single_d``).

    D-D is the default for the D-bearing loci; ``single_d=True``, a VJ locus, or an already-tandem
    template all leave it unchanged.
    """
    from .dd import has_tandem, to_dd
    if single_d or template.chain_type != "VDJ" or has_tandem(template):
        return template
    if template.manifest.locus not in DD_DEFAULT_LOCI:
        return template
    return to_dd(template, p_nd2=p_nd2_init)


def infer(
    template: Model,
    sequences: list[str],
    *,
    max_iter: int = 30,
    tol: float = 1e-3,
    init: str = "align",
    masks: list | None = None,
    single_d: bool = False,
    p_nd2_init: float = 0.02,
    dd_allowed: list | None = None,
    nd_prior: float = 0.0,
) -> tuple[Model, InferenceReport]:
    """Re-estimate a model's marginals from nucleotide CDR3s by EM.

    Args:
        template: A model supplying the gene set, germline, and event graph (its marginals are
            replaced). Use one built by ``from_olga`` (or any :class:`Model`).
        sequences: Observed CDR3 nucleotide strings (typically out-of-frame reads).
        max_iter: Maximum EM iterations.
        tol: Stop when the V-usage total-variation between iterations falls below this.
        init: ``"align"`` (seed gene usage from a best-match vote — the default and fastest),
            ``"uniform"`` (each event uniform on its support), or ``"template"`` (warm start).
        masks: Optional per-sequence ``(v_genes, j_genes, d_genes)`` name lists (e.g. from
            :func:`arda_masks`) restricting each read's scenario enumeration to its aligned genes.
            **Strongly recommended for VDJ** — without it the E-step enumerates every Cys-sharing
            V × the full D grid per read (tens of s/seq); with it, VDJ inference is tractable.
        single_d: By default a **tandem-D (D-D)** model is learned for the D-bearing loci
            (IGH/TRD/TRB): a single-D template is promoted with :func:`~vdjtools.model.dd.to_dd`
            (seeding ``P(n_D=2)=p_nd2_init``) and EM learns the true ``P(n_D=2)``. Set ``True`` to
            keep a strict single-D model. No effect on VJ loci or an already-tandem template.
        p_nd2_init: Initial ``P(n_D=2)`` seed when promoting to D-D (ignored if ``single_d``).
        dd_allowed: Optional per-read booleans gating the tandem (``n_D=2``) E-step — a read may be
            tandem only where ``dd_allowed[i]`` is true (e.g. reads arda flags with a ``d2_call``).
            Anchors D-D learning to alignment-detected tandems, countering the tandem-vs-long-insertion
            identifiability that inflates unregularized ``P(n_D=2)`` on real data. ``None`` = all reads.
        nd_prior: Dirichlet/Beta pseudocount added to the single-D (``n_D=1``) soft count each M-step,
            regularizing ``P(n_D=2)`` toward 0. Both anchors combine.

    Returns:
        ``(fitted_model, report)``. For a tandem-D template the E-step enumerates ``n_D=2``
        scenarios and the M-step learns ``P(n_D=2)`` along with the ``d2_gene`` / ``d2_del`` / ``dd``
        events.
    """
    template = _maybe_promote_dd(template, single_d, p_nd2_init)
    upper = [s.upper() for s in sequences]
    if init == "template":
        tables = template.tables
    elif init == "align":
        tables = _align_init(template, upper)
    else:
        tables = _uniform_init(template)
    model = Model(manifest=template.manifest, tables=tables, genomic=template.genomic)
    report = InferenceReport()
    fit = _fit_events(template.manifest)
    seq_masks = masks if masks is not None else [None] * len(upper)
    dd_gate = dd_allowed if dd_allowed is not None else [True] * len(upper)

    for it in range(max_iter):
        prep = prepare(model)
        counts = {name: defaultdict(float) for name in fit}
        ll = 0.0
        n_ok = 0
        for s, mask, allow_dd in zip(upper, seq_masks, dd_gate):
            pg = _estep_seq(prep, s, counts, mask, bool(allow_dd))
            if pg > 0.0:
                ll += log(pg)
                n_ok += 1
        report.loglik.append(ll / n_ok if n_ok else float("-inf"))
        report.n_scoreable.append(n_ok)
        report.n_iter = it + 1

        new_tables = _mstep(template, counts, nd_prior)
        new_model = Model(manifest=template.manifest, tables={**model.tables, **new_tables}, genomic=template.genomic)
        # Converge on marginal stability (V-usage total variation), which is robust to the
        # changing set of scoreable reads that makes raw mean-log-lik non-monotonic.
        tv = _tv(model.tables["v_choice"], new_model.tables["v_choice"], ["v_allele"])
        report.gene_tv.append(tv)
        model = new_model
        if it > 0 and tv < tol:
            report.converged = True
            break

    return model, report


def _gene_to_alleles(model: Model, seg: str) -> dict[str, list[str]]:
    out: dict[str, list[str]] = defaultdict(list)
    for a in model.genomic[f"genes_{seg}"][f"{seg}_allele"]:
        out[a.split("*")[0]].append(a)
    return out


def call_alleles(index: dict[str, list[str]], call: str | None) -> list[str]:
    """All model alleles compatible with one AIRR gene call, ambiguity included.

    Two kinds of ambiguity, and both must widen the mask rather than narrow it:

    * **allele** — a call of ``TRBV20-1*03`` where the truth is ``*01``. Expanding to every model
      allele of the gene keeps the right scenario reachable.
    * **comma-separated genes** — AIRR writes an aligner's tie as ``IGHV3-23*01,IGHV3-23D*01``,
      which means *the aligner could not tell these apart*. Splitting on ``*`` alone keeps only
      ``IGHV3-23`` and silently DROPS ``IGHV3-23D`` — a different gene on a duplicated locus. If
      the truth is the dropped one, its scenario is unreachable and EM misattributes the read.
      Measured on human IGH: 23,176 of 160,324 non-functional clonotypes (14.5%) carry an
      ambiguous V call; TRB 2.0%; TRA/TRD 0%.

    Returns the union over every gene named, deduplicated and order-stable. Unknown genes
    contribute nothing; a call naming no known gene yields ``[]``, which the E-step reads as
    "unrestricted" — the honest degradation, since we know nothing about that read's gene.
    """
    if not call:
        return []
    out: list[str] = []
    seen: set[str] = set()
    for part in call.split(","):
        for a in index.get(part.strip().split("*")[0], []):
            if a not in seen:
                seen.add(a)
                out.append(a)
    return out


def gene_masks(model: Model, v_calls: list[str], j_calls: list[str]) -> list[tuple]:
    """Build per-read ``(v_genes, j_genes, d_genes)`` E-step masks from V/J gene calls.

    Each call is expanded to every model allele of every gene it names — see
    :func:`call_alleles` for why both allele- and comma-ambiguity must widen the mask.
    D is left unrestricted (few D genes, and D calls on the short D germline are unreliable).
    """
    va, ja = _gene_to_alleles(model, "v"), _gene_to_alleles(model, "j")
    return [(call_alleles(va, v), call_alleles(ja, j), None) for v, j in zip(v_calls, j_calls)]


def arda_masks(contigs: list[str], model: Model, *, organism: str = "human") -> tuple[list[str], list[tuple]]:
    """Annotate nt contigs with arda and build ``(junctions, masks)`` for masked :func:`infer`.

    The production path for real reads: ``junctions, masks = arda_masks(contigs, template);
    infer_native(template, junctions, masks=masks)``. arda is a base dependency (ships with vdjtools).
    """
    from .stitch import annotate

    calls = annotate(contigs, organism=organism)
    junctions = calls["junction"].to_list()
    masks = gene_masks(model, calls["v_call"].to_list(), calls["j_call"].to_list())
    return junctions, masks


# --------------------------------------------------------------------------------------------
# Native EM: the E-step runs in C++ (_core.estep_batch); the M-step re-normalizes the returned
# dense soft-count arrays back into polars tables. Same result as pure-Python infer(), much faster.

def _functional_support(template: Model, seg: str) -> set[str]:
    """The allele names the germline says are real and functional for one segment.

    This is the support a ``gene_prior`` is spread over: the germline reference (arda for
    arda-native models, OLGA's own for ``from_olga`` ones) is the authority on which alleles
    exist, and EM should never be able to declare one of them impossible.
    """
    g = template.genomic.get(f"genes_{seg}")
    if g is None:
        return set()
    return set(g.filter(pl.col("functional"))[f"{seg}_allele"].to_list())


def _mstep_native(template: Model, counts, v_alleles, j_alleles, d_alleles, nbins, nd_prior=0.0,
                  gene_prior=0.0) -> dict[str, pl.DataFrame]:
    vdj = template.chain_type == "VDJ"
    mp = template.manifest.palindrome_max
    nV, nJ, nD = len(v_alleles), len(j_alleles), len(d_alleles)

    def norm(df, keys):
        tot = pl.col("p").sum().over(keys) if keys else pl.col("p").sum()
        return df.with_columns(p=pl.when(tot > 0).then(pl.col("p") / tot).otherwise(0.0))

    def gene_choice(df, allele_col, seg, keys):
        """Normalize a gene-choice table, optionally with a Dirichlet prior over the germline.

        ``P(V)=0`` is an ABSORBING STATE of this EM: the E-step weights every scenario by P(V),
        so an allele that reaches zero count can never be re-attributed and is dead for good. One
        unlucky iteration permanently deletes a real gene. Measured on human TRB: the learned
        model kept 30 of the 57 V genes OLGA carries, having seen 54 of them in the training data.

        A pseudocount over the germline's FUNCTIONAL alleles fixes it at the source: the germline
        is the authority on what exists, so nothing real is ever assigned probability zero, while
        the data still decides the actual usage. Non-functional alleles (pseudogenes/ORFs) are
        deliberately NOT given mass -- the model cannot score them anyway.

        The prior only protects alleles the data actually ATTRIBUTED reads to (soft count > 0),
        NOT every functional allele: the E-step commits V-choice to one best-match allele per read,
        so a functional secondary allele (e.g. TRDV2*03 when arda calls TRDV2*01) gets zero soft
        count AND zero deletion/insertion counts. Handing it choice mass anyway makes it
        selectable by the generative sampler while its deletion distribution is all-zero -> the
        sampler draws it and then has nothing to draw a deletion from (IndexError). Guarding on
        ``p > 0`` keeps every mass-bearing allele conditionally complete, which is what
        `absorbing state' protection actually requires: rescue what was seen, do not invent what
        was not (rescale_usage covers the cross-protocol `give me every gene' case separately).

        ``gene_prior=0.0`` (the default) is byte-identical to plain MLE normalization, so the
        exact-Pgen invariant on ``from_olga`` models is untouched.
        """
        if not gene_prior:
            return norm(df, keys)
        ok = _functional_support(template, seg)
        return norm(df.with_columns(
            p=pl.col("p") + pl.when(pl.col(allele_col).is_in(list(ok)) & (pl.col("p") > 0))
                              .then(gene_prior).otherwise(0.0)
        ), keys)

    def deletion(arr, alleles, col, nb, maxpal):
        a = np.repeat(alleles, nb)
        ndel = np.tile(np.arange(nb) - maxpal, len(alleles))
        return norm(pl.DataFrame({col: a, "ndel": ndel.astype(np.int16), "p": list(arr)}), [col])

    def dinucl(arr):
        i = np.arange(16)  # arr index i = to*4 + from
        return norm(pl.DataFrame({"from_nt": (i % 4).astype(np.uint8), "to_nt": (i // 4).astype(np.uint8),
                                  "p": list(arr)}), ["from_nt"])

    t = {}
    t["v_choice"] = gene_choice(pl.DataFrame({"v_allele": v_alleles, "p": list(counts.v_choice)}), "v_allele", "v", [])
    t["v_3_del"] = deletion(counts.v_3_del, v_alleles, "v_allele", nbins["v"], mp["v_3"])
    t["j_5_del"] = deletion(counts.j_5_del, j_alleles, "j_allele", nbins["j"], mp["j_5"])
    if vdj:
        t["j_choice"] = gene_choice(pl.DataFrame({"j_allele": j_alleles, "p": list(counts.j_choice)}), "j_allele", "j", [])
        t["d_gene"] = norm(pl.DataFrame({
            "j_allele": np.repeat(j_alleles, nD), "d_allele": np.tile(d_alleles, nJ),
            "p": list(counts.d_gene)}), ["j_allele"])
        n5, n3 = nbins["d5"], nbins["d3"]
        t["d_del"] = norm(pl.DataFrame({
            "d_allele": np.repeat(d_alleles, n5 * n3),
            "ndel5": np.tile(np.repeat(np.arange(n5) - mp["d_5"], n3), nD).astype(np.int16),
            "ndel3": np.tile(np.arange(n3) - mp["d_3"], n5 * nD).astype(np.int16),
            "p": list(counts.d_del)}), ["d_allele"])
        t["vd_ins"] = norm(pl.DataFrame({"length": np.arange(len(counts.ins_vd), dtype=np.int16), "p": list(counts.ins_vd)}), [])
        t["dj_ins"] = norm(pl.DataFrame({"length": np.arange(len(counts.ins_dj), dtype=np.int16), "p": list(counts.ins_dj)}), [])
        t["vd_dinucl"] = dinucl(counts.dinucl_vd)
        t["dj_dinucl"] = dinucl(counts.dinucl_dj)
        # n_d learned from the soft counts (indexed by the n_D value): single-D emits only n_D=1 mass
        # -> renormalizes to delta(1) (a no-op); D-D emits n_D=1 and n_D=2 -> learns P(n_D=2).
        nd = template.tables["n_d"]
        # Dirichlet/Beta prior: nd_prior pseudocounts added to the single-D (n_D=1) bucket regularize
        # away the tandem over-attribution of unregularized D-D EM on real data.
        nd_counts = [float(counts.n_d[int(k)]) + (nd_prior if int(k) == 1 else 0.0) for k in nd["n_d"]]
        t["n_d"] = norm(nd.with_columns(p=pl.Series("p", nd_counts)), [])
        if "d2_gene" in template.manifest.events:  # tandem-D model — learn the second-D events too
            t["d2_gene"] = norm(pl.DataFrame({
                "d_allele": np.repeat(d_alleles, nD), "d2_allele": np.tile(d_alleles, nD),
                "p": list(counts.d2_gene)}), ["d_allele"])
            t["d2_del"] = norm(pl.DataFrame({
                "d2_allele": np.repeat(d_alleles, n5 * n3),
                "ndel5": np.tile(np.repeat(np.arange(n5) - mp["d_5"], n3), nD).astype(np.int16),
                "ndel3": np.tile(np.arange(n3) - mp["d_3"], n5 * nD).astype(np.int16),
                "p": list(counts.d2_del)}), ["d2_allele"])
            t["dd_ins"] = norm(pl.DataFrame({"length": np.arange(len(counts.ins_dd), dtype=np.int16), "p": list(counts.ins_dd)}), [])
            t["dd_dinucl"] = dinucl(counts.dinucl_dd)
    else:
        t["j_choice"] = norm(pl.DataFrame({
            "v_allele": np.repeat(v_alleles, nJ), "j_allele": np.tile(j_alleles, nV),
            "p": list(counts.j_choice)}), ["v_allele"])
        t["vj_ins"] = norm(pl.DataFrame({"length": np.arange(len(counts.ins_vj), dtype=np.int16), "p": list(counts.ins_vj)}), [])
        t["vj_dinucl"] = dinucl(counts.dinucl_vj)
    return t


def infer_native(
    template: Model,
    sequences: list[str],
    *,
    max_iter: int = 30,
    tol: float = 1e-3,
    init: str = "align",
    masks: list | None = None,
    single_d: bool = False,
    p_nd2_init: float = 0.02,
    dd_allowed: list | None = None,
    nd_prior: float = 0.0,
    gene_prior: float = 0.0,
) -> tuple[Model, InferenceReport]:
    """EM inference with the native C++ E-step — same result as :func:`infer`, much faster.

    Requires the compiled ``_core`` extension. See :func:`infer` for the arguments (including
    ``single_d`` / ``p_nd2_init`` / ``dd_allowed`` / ``nd_prior`` / ``gene_prior``).

    ``gene_prior`` is a Dirichlet pseudocount spread over the germline's **functional** V/J
    alleles in each M-step. ``P(V)=0`` is an absorbing state of this EM — the E-step weights
    scenarios by P(V), so a zeroed allele can never be re-attributed — and on real data that
    silently deletes real genes for good (human TRB: 30 of 57 V genes survived unregularized).
    The germline says which alleles exist; the prior keeps all of them reachable and lets the
    data set the usage. ``0.0`` (default) is byte-identical to plain MLE. Learns tandem-D (``n_D=2``) by
    default on the D-bearing loci: the native E-step accumulates the second-D soft counts via a
    factorized forward/backward pass, read-parallelized across cores.
    """
    from .._core import estep_batch, make_counts
    from .native import _encode, pack

    template = _maybe_promote_dd(template, single_d, p_nd2_init)
    ddflags = [1 if a else 0 for a in dd_allowed] if dd_allowed is not None else []

    upper = [s.upper() for s in sequences]
    if init == "template":
        tables = template.tables
    elif init == "align":
        tables = _align_init(template, upper)
    else:
        tables = _uniform_init(template)
    model = Model(manifest=template.manifest, tables=tables, genomic=template.genomic)

    v_alleles = template.genomic["genes_v"]["v_allele"].to_list()
    j_alleles = template.genomic["genes_j"]["j_allele"].to_list()
    d_alleles = template.genomic["genes_d"]["d_allele"].to_list() if template.chain_type == "VDJ" else []
    vi = {a: i for i, a in enumerate(v_alleles)}
    ji = {a: i for i, a in enumerate(j_alleles)}
    di = {a: i for i, a in enumerate(d_alleles)}
    seqs_enc = [_encode(s) for s in upper]
    if masks is not None:
        vmasks = [[vi[a] for a in mk[0] if a in vi] for mk in masks]
        jmasks = [[ji[a] for a in mk[1] if a in ji] for mk in masks]
        dmasks = [[di[a] for a in (mk[2] or []) if a in di] for mk in masks]
    else:
        vmasks = jmasks = dmasks = []

    report = InferenceReport()
    for it in range(max_iter):
        pm, _, _ = pack(model)
        counts = make_counts(pm)
        ll = estep_batch(pm, seqs_enc, vmasks, jmasks, dmasks, counts, 0, ddflags)
        # Report the per-sequence MEAN, matching infer() (line ~431) -- estep_batch returns the
        # summed log-likelihood, and appending it raw made infer_native's loglik differ from
        # infer's by a factor of n (the two are documented as "same result"). Denominator is the
        # input read count; n_scoreable is populated so the sum can be recovered if needed.
        n = len(seqs_enc)
        report.loglik.append(ll / n if n else float("-inf"))
        report.n_scoreable.append(n)
        nbins = {"v": pm.nbins_v, "j": pm.nbins_j, "d5": pm.nbins_d5, "d3": pm.nbins_d3}
        new_tables = _mstep_native(template, counts, v_alleles, j_alleles, d_alleles, nbins, nd_prior, gene_prior)
        new_model = Model(manifest=template.manifest, tables={**model.tables, **new_tables}, genomic=template.genomic)
        tv = _tv(model.tables["v_choice"], new_model.tables["v_choice"], ["v_allele"])
        report.gene_tv.append(tv)
        report.n_iter = it + 1
        model = new_model
        if it > 0 and tv < tol:
            report.converged = True
            break
    return model, report
