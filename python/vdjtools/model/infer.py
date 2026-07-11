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


def _estep_seq(prep, s: str, counts: dict, mask=None) -> float:
    """Accumulate one sequence's soft counts into ``counts``; return its Pgen (for log-lik).

    ``mask`` optionally restricts enumeration to a read's aligned genes — ``(v_genes, j_genes,
    d_genes)`` name lists (e.g. from arda). This is what makes VDJ inference tractable: without it
    every V that shares the conserved Cys prefix is a candidate and the D enumeration runs for each.
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
                        total += _accum_vdj(prep, J, V, nv, nj, mid, base, local, d_mask)
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
                    for fr, to in tvd:
                        local[("vd_dinucl", (fr, to))] += w
                    for fr, to in tdj:
                        local[("dj_dinucl", (fr, to))] += w
    return seq_total


# Events re-estimated by EM (n_d and the germline stay fixed from the template).
def _fit_events(manifest) -> list[str]:
    return [n for n, ev in manifest.events.items() if ev.kind.value != "n_d"]


def _mstep(template: Model, counts: dict) -> dict[str, pl.DataFrame]:
    tables = {}
    for name, event in template.manifest.events.items():
        if event.kind.value == "n_d":
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


def infer(
    template: Model,
    sequences: list[str],
    *,
    max_iter: int = 30,
    tol: float = 1e-3,
    init: str = "align",
    masks: list | None = None,
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

    Returns:
        ``(fitted_model, report)``.
    """
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

    for it in range(max_iter):
        prep = prepare(model)
        counts = {name: defaultdict(float) for name in fit}
        ll = 0.0
        n_ok = 0
        for s, mask in zip(upper, seq_masks):
            pg = _estep_seq(prep, s, counts, mask)
            if pg > 0.0:
                ll += log(pg)
                n_ok += 1
        report.loglik.append(ll / n_ok if n_ok else float("-inf"))
        report.n_scoreable.append(n_ok)
        report.n_iter = it + 1

        new_tables = _mstep(template, counts)
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


def gene_masks(model: Model, v_calls: list[str], j_calls: list[str]) -> list[tuple]:
    """Build per-read ``(v_genes, j_genes, d_genes)`` E-step masks from V/J gene calls.

    Each call is expanded from an allele to *all* model alleles of its gene, so allele-level
    ambiguity (a call of ``TRBV20-1*03`` vs the true ``*01``) never excludes the right scenario.
    D is left unrestricted (few D genes, and D calls on the short D germline are unreliable).
    """
    va, ja = _gene_to_alleles(model, "v"), _gene_to_alleles(model, "j")

    def alleles(m, call):
        return m.get(call.split("*")[0], []) if call else []

    return [(alleles(va, v), alleles(ja, j), None) for v, j in zip(v_calls, j_calls)]


def arda_masks(contigs: list[str], model: Model, *, organism: str = "human") -> tuple[list[str], list[tuple]]:
    """Annotate nt contigs with arda and build ``(junctions, masks)`` for masked :func:`infer`.

    The production path for real reads: ``junctions, masks = arda_masks(contigs, template);
    infer(template, junctions, masks=masks)``. arda is the ``[model]`` extra.
    """
    from .stitch import annotate

    calls = annotate(contigs, organism=organism)
    junctions = calls["junction"].to_list()
    masks = gene_masks(model, calls["v_call"].to_list(), calls["j_call"].to_list())
    return junctions, masks
