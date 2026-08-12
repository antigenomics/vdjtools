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

import warnings
from collections import defaultdict
from dataclasses import dataclass, field, fields
from math import log

import numpy as np
import polars as pl

from .model import Model
from .pgen import _NT2NUM, _common_prefix, _common_suffix, _j_candidates, _v_candidates, prepare
from .schema import normalization_keys, table_columns


@dataclass(slots=True)
class InferenceReport:
    """Per-iteration diagnostics from :func:`infer` — the training log.

    Every field after ``converged`` is metadata recorded so a saved model can say what it was
    fitted on and how; all are defaulted, so constructing a bare report still works.
    """

    loglik: list[float] = field(default_factory=list)  # mean per-sequence log-Pgen (over scoreable reads)
    n_scoreable: list[int] = field(default_factory=list)
    gene_tv: list[float] = field(default_factory=list)  # relative log-likelihood change vs previous iter (the convergence signal)
    n_iter: int = 0
    converged: bool = False
    n_sequences: int = 0
    max_iter: int = 0
    tol: float = 0.0
    init: str = ""
    native: bool = False
    elapsed_s: float = 0.0
    finished_at: str = ""
    template_source: str = ""

    def to_dict(self) -> dict:
        """The report as a JSON-serializable dict (one entry of a model's ``training["runs"]``)."""
        from dataclasses import asdict

        return asdict(self)

    @classmethod
    def from_dict(cls, obj: dict) -> "InferenceReport":
        """Rebuild a report from :meth:`to_dict`, ignoring keys this version does not know."""
        known = {f.name for f in fields(cls)}
        return cls(**{k: v for k, v in obj.items() if k in known})

    def to_frame(self) -> pl.DataFrame:
        """Per-iteration training log: ``iter, loglik, n_scoreable, rel_change``."""
        n = len(self.loglik)
        return pl.DataFrame({
            "iter": list(range(1, n + 1)),
            "loglik": self.loglik,
            "n_scoreable": (self.n_scoreable + [None] * n)[:n],
            "rel_change": (self.gene_tv + [None] * n)[:n],
        })


def _record(report: InferenceReport, model: Model, template: Model, *, started: float,
            max_iter: int, tol: float, init: str, native: bool, n_sequences: int) -> Model:
    """Stamp a finished run onto the report and append it to the model's training log.

    Appends rather than overwrites, so a warm-start refit of an already-fitted model keeps both
    runs and the history reads as the sequence of things that were actually done to the model.
    """
    import time as _time

    report.n_sequences = n_sequences
    report.max_iter = max_iter
    report.tol = tol
    report.init = init
    report.native = native
    report.elapsed_s = round(_time.monotonic() - started, 3)
    report.finished_at = _time.strftime("%Y-%m-%dT%H:%M:%S", _time.gmtime())
    report.template_source = template.manifest.source
    runs = list((template.training or {}).get("runs", []))
    runs.append(report.to_dict())
    return Model(manifest=model.manifest, tables=model.tables, genomic=model.genomic,
                 training={"runs": runs})


def _write_checkpoint(path, model: Model, template: Model, report: InferenceReport,
                      meta: dict) -> None:
    """Save the model as it stands, with the training log up to this iteration.

    Written atomically-ish (to a sibling directory, then swapped) so a job killed mid-write leaves
    the previous checkpoint intact rather than a half-written model.
    """
    import shutil
    from pathlib import Path

    from .io import save_model

    path = Path(path)
    snapshot = _record(report, model, template, **meta)
    tmp = path.with_name(path.name + ".partial")
    if tmp.exists():
        shutil.rmtree(tmp)
    save_model(snapshot, tmp)
    if path.exists():
        shutil.rmtree(path)
    tmp.rename(path)


def resume(path, sequences, **kw):
    """Continue EM from a checkpoint written by ``infer_native(checkpoint=...)``.

    A warm start from the saved marginals, so the fit picks up where it stopped rather than
    realigning from scratch. The checkpoint carries its own training log and this run **appends**
    to it, so the full history survives across as many interruptions as it takes — which is what
    makes a multi-hour fit on a time-limited queue practical.

    Args:
        path: The checkpoint directory (or an already-loaded :class:`Model`).
        sequences: The same sequences the interrupted run was fitting.
        **kw: Passed to :func:`infer_native`. ``init`` is forced to ``"template"``, and
            ``checkpoint`` defaults to ``path`` so the continued run keeps saving to the same
            place — pass ``checkpoint=None`` to stop checkpointing.

    Returns:
        ``(model, report)`` — ``report`` covers this run only; ``training_frame(model)`` shows
        every run.

    Example:
        >>> model, rep = resume("ckpt/IGH", seqs, max_iter=10)
    """
    from .io import load_model

    base = path if isinstance(path, Model) else load_model(path)
    if not isinstance(path, Model):
        kw.setdefault("checkpoint", path)
    kw["init"] = "template"
    return infer_native(base, sequences, **kw)


def print_progress(stream=None, prefix: str = ""):
    """A ready-made ``progress=`` callback that reports each EM iteration as it happens.

    EM on a large D-bearing locus runs for tens of minutes with nothing to show for it, and the
    training log only becomes readable once the fit *returns* — so a long run is indistinguishable
    from a hung one. This prints the log-likelihood and its relative change per iteration, which is
    exactly the quantity the convergence test uses, so you can watch it approach ``tol``.

    Args:
        stream: Where to write; defaults to ``sys.stderr`` so it never pollutes piped output.
        prefix: Prepended to each line, e.g. the locus being built.

    Returns:
        A callable suitable for ``infer(progress=...)`` / ``infer_native(progress=...)``.

    Example:
        >>> infer_native(template, seqs, progress=print_progress(prefix="[TRB] "))
        [TRB] iter  1  loglik -37.0067  rel      inf  n=122703
        [TRB] iter  2  loglik -33.9738  rel 8.20e-02  n=122703
    """
    import sys

    out = stream if stream is not None else sys.stderr

    def _report(iteration: int, loglik: float, rel: float, n_scoreable: int) -> None:
        rel_s = "     inf" if rel == float("inf") else f"{rel:.2e}"
        print(f"{prefix}iter {iteration:2d}  loglik {loglik:9.4f}  rel {rel_s}  n={n_scoreable}",
              file=out, flush=True)

    return _report


def training_frame(obj) -> pl.DataFrame:
    """The training log of a model (or a single report) as one tidy frame across all runs.

    Args:
        obj: A :class:`~vdjtools.model.model.Model` carrying a ``training`` log, or an
            :class:`InferenceReport`.

    Returns:
        ``run, iter, loglik, n_scoreable, rel_change`` — empty when the model was never fitted here
        (every bundled model, and anything imported straight from OLGA).

    Example:
        >>> m, rep = infer_native(template, seqs, max_iter=10)
        >>> training_frame(m)     # loglik per iteration, ready to plot
    """
    if isinstance(obj, InferenceReport):
        return obj.to_frame().with_columns(run=pl.lit(0, dtype=pl.Int64)).select(
            ["run", "iter", "loglik", "n_scoreable", "rel_change"])
    runs = (getattr(obj, "training", None) or {}).get("runs", [])
    if not runs:
        return pl.DataFrame(schema={"run": pl.Int64, "iter": pl.Int64, "loglik": pl.Float64,
                                    "n_scoreable": pl.Int64, "rel_change": pl.Float64})
    parts = [InferenceReport.from_dict(r).to_frame().with_columns(run=pl.lit(i, dtype=pl.Int64))
             for i, r in enumerate(runs)]
    return pl.concat(parts).select(["run", "iter", "loglik", "n_scoreable", "rel_change"])


def _loglik_rel(loglik: list[float]) -> float:
    """Relative change in mean log-likelihood between the last two iterations (``inf`` before iter 1).

    The natural EM convergence signal, and it reflects the WHOLE model (trims/insertions/D-D), not just
    V usage — which, arda-masked, settles in ~2 iterations while those are still moving. It is usable
    because the fixes keep the scoreable-read set stable across iterations, so mean log-lik is now
    monotone (the non-monotonicity that once forced a V-usage criterion is gone). Scale-free, so one
    ``tol`` works across loci whose absolute log-lik differs (TRA ≈ −19, IGK ≈ −13, TRD ≈ −35).
    """
    if len(loglik) < 2:
        return float("inf")
    prev, cur = loglik[-2], loglik[-1]
    return abs(cur - prev) / (abs(prev) + 1e-12)


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
        # Split each read's vote across ALL genes tied for the longest germline match, not the first.
        # Germline-identical paralogs (TRBV6-2/6-5/6-6, IGKV2-28/2D-28) tie EXACTLY, so a single-winner
        # vote (Python ``max`` returns the first) hands the whole family to one representative and seeds
        # the rest at P(V)=0 — which the E-step's ``if pv==0: continue`` then makes an absorbing state
        # no amount of data escapes. Sharing the tie keeps every indistinguishable sibling reachable.
        vsc = [(_common_prefix(prep.cut["v"][v], s), v) for v in prep.functional_v]
        jsc = [(_common_suffix(prep.cut["j"][j], s), j) for j in prep.functional_j]
        bv_best = max(sc for sc, _ in vsc)
        bj_best = max(sc for sc, _ in jsc)
        bvs = [v for sc, v in vsc if sc == bv_best]
        bjs = [j for sc, j in jsc if sc == bj_best]
        wv, wj = 1.0 / len(bvs), 1.0 / len(bjs)
        for v in bvs:
            v_votes[v] += wv
        for j in bjs:
            j_votes[j] += wj
        for v in bvs:
            for j in bjs:
                vj_votes[(v, j)] += wv * wj

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
    progress=None,
    checkpoint=None,
    checkpoint_every: int = 1,
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
        progress: Optional ``callable(iteration, loglik, rel_change, n_scoreable)`` invoked after
            every iteration — use :func:`print_progress` to watch a long fit converge live.
        checkpoint: Directory to save the model into after each iteration, so a long fit survives
            being interrupted. Continue it with :func:`resume`. The checkpoint carries the training
            log so far, and the resumed run appends to it.
        checkpoint_every: Write a checkpoint every N iterations (default every one).

    Returns:
        ``(fitted_model, report)``. For a tandem-D template the E-step enumerates ``n_D=2``
        scenarios and the M-step learns ``P(n_D=2)`` along with the ``d2_gene`` / ``d2_del`` / ``dd``
        events.
    """
    import time as _time

    started = _time.monotonic()
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
    meta = dict(started=started, max_iter=max_iter, tol=tol, init=init, native=False,
                n_sequences=len(upper))
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
        model = new_model
        # Converge on the relative log-likelihood improvement (whole-model, monotone post-fix), not V
        # usage alone: V is arda-masked so it settles in ~2 iters while trims/insertions/D-D still move.
        rel = _loglik_rel(report.loglik)
        report.gene_tv.append(rel)
        if progress is not None:
            progress(report.n_iter, report.loglik[-1], rel, report.n_scoreable[-1])
        if checkpoint is not None and report.n_iter % checkpoint_every == 0:
            _write_checkpoint(checkpoint, model, template, report, meta)
        if it > 0 and rel < tol:
            report.converged = True
            break

    model = _record(report, model, template, **meta)
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


#: Columns `infer_frame` will look for as the nucleotide junction, in preference order.
_JUNCTION_COLS = ("junction", "junction_nt", "cdr3_nt", "cdr3nt", "sequence")


def sanitize_junctions(df: pl.DataFrame, col: str, *, ambiguous: str | None = "A",
                       where: str = "infer_frame") -> pl.DataFrame:
    """Make a junction column safe for the native encoder, which knows only A/C/G/T.

    Real annotated reads carry the occasional ambiguous base (``N``, or an IUPAC code), and an
    unhandled one surfaces as a ``KeyError`` from deep inside the E-step.

    Args:
        df: Clonotype frame.
        col: Junction column.
        ambiguous: A single base to substitute for every non-ACGT character (default ``"A"``), or
            ``None`` to drop those clonotypes instead. Substituting is the default because it keeps
            the clonotype: an ambiguous base is one uncertain position in a junction that is
            otherwise perfectly good evidence, and on these reads it affects ~0.01% of rows, so
            dropping costs sample size for no gain in correctness. It is a *substitution*, not a
            marginalization — the base is treated as read, so a run with many ambiguous positions
            will bias the insertion model toward the substituted base and should use ``None``.
        where: Caller name, used in the warning.

    Returns:
        A frame with the column uppercased and cleaned; rows are dropped only when
        ``ambiguous is None``.
    """
    df = df.with_columns(pl.col(col).str.to_uppercase())
    bad = df.filter(~pl.col(col).str.contains(r"^[ACGT]+$")).height
    if not bad:
        return df
    if ambiguous is None:
        warnings.warn(
            f"{where}: dropped {bad} of {df.height} clonotypes whose {col!r} contains non-ACGT "
            f"bases (the recombination model is defined over A/C/G/T only)", stacklevel=3)
        return df.filter(pl.col(col).str.contains(r"^[ACGT]+$"))
    if ambiguous not in ("A", "C", "G", "T"):
        raise ValueError(f"ambiguous must be one of A, C, G, T or None, got {ambiguous!r}")
    warnings.warn(
        f"{where}: {bad} of {df.height} clonotypes have non-ACGT bases in {col!r}; substituting "
        f"{ambiguous!r} (pass ambiguous=None to drop them instead)", stacklevel=3)
    return df.with_columns(pl.col(col).str.replace_all(r"[^ACGT]", ambiguous))


def infer_frame(template, clones: pl.DataFrame, *, seq_col: str | None = None,
                v_col: str = "v_call", j_col: str = "j_call", use_calls: bool = True,
                native: bool = True, ambiguous: str | None = "A", **kw):
    """Fit a model from a **clonotype frame** — the ergonomic entry point to EM.

    Wraps :func:`infer_native` with the two steps every caller otherwise repeats: find the
    nucleotide junction column, and turn the frame's V/J calls into per-read E-step masks with
    :func:`gene_masks`. The masks matter enormously on a D-bearing locus — without them the E-step
    enumerates every Cys-sharing V against the full D grid for every read.

    Args:
        template: A :class:`~vdjtools.model.model.Model` supplying the gene set, germline and event
            graph, **or** a locus string (e.g. ``"TRB"``) to build one with
            :func:`~vdjtools.model.io.from_arda`.
        clones: Clonotype frame. Needs a nucleotide junction column and, for ``use_calls``,
            ``v_call``/``j_call``.
        seq_col: Explicit junction column; auto-detected from :data:`_JUNCTION_COLS` otherwise.
        v_col: V-call column.
        j_col: J-call column.
        use_calls: Build per-read masks from the V/J calls. Turn off only if the frame's calls are
            untrustworthy — inference then enumerates every gene and gets much slower.
        native: Use :func:`infer_native` (default). ``False`` runs the pure-Python :func:`infer`.
        ambiguous: What to do with a junction holding a non-ACGT base — substitute this base
            (default ``"A"``), or ``None`` to drop the clonotype. See :func:`sanitize_junctions`.
        **kw: Passed through (``max_iter``, ``tol``, ``init``, ``single_d``, ``nd_prior``,
            ``gene_prior``, ...).

    Returns:
        ``(fitted_model, report)`` — the model carries the run in its ``training`` log.

    Raises:
        ValueError: If no junction column is found, or it holds no usable sequences.

    Example:
        >>> m, rep = infer_frame("TRB", clones, max_iter=10)
        >>> training_frame(m)
    """
    if isinstance(template, str):
        from .io import from_arda

        template = from_arda(template)
    col = seq_col or next((c for c in _JUNCTION_COLS if c in clones.columns), None)
    if col is None:
        raise ValueError(
            f"no nucleotide junction column in {clones.columns}; pass seq_col=")
    df = clones.filter(pl.col(col).is_not_null() & (pl.col(col).str.len_bytes() > 0))
    df = sanitize_junctions(df, col, ambiguous=ambiguous, where="infer_frame")
    if not df.height:
        raise ValueError(f"column {col!r} has no usable sequences")
    seqs = df[col].to_list()
    masks = None
    if use_calls and v_col in df.columns and j_col in df.columns:
        masks = gene_masks(template, df[v_col].to_list(), df[j_col].to_list())
    fn = infer_native if native else infer
    return fn(template, seqs, masks=masks, **kw)


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


def augment_from_oracle(learned: Model, oracle: Model) -> Model:
    """Fill functional genes the learned model left at P=0 with the ORACLE's own usage and conditionals.

    A learned model carries only the genes arda saw *producibly* in the training repertoire; a user's
    library (different protocol/tissue) can be full of genes 5'RACE or this cohort never amplified, or
    that arda called but the germline can't emit (a source mismatch, or a hard-call tie it lost). The
    OLGA oracle models every functional gene it knows with a real per-gene usage and deletion profile,
    and differs from the learned model ONLY in D/D-D handling — orthogonal to V/J gene identity — so its
    V/J genes transplant cleanly. For each functional (scoreable) gene present in the oracle but absent
    from the learned model, this copies the oracle's choice mass AND every child table the gene parents
    (its deletion profile, P(J|V) on a VJ locus, P(D|J) for a J), then renormalizes — so no functional
    gene is silently missing. Usage is the oracle's here; :func:`~vdjtools.model.rescale.rescale_usage`
    adapts it to the user's actual library (its cross-protocol job). Idempotent.
    """
    from difflib import SequenceMatcher

    tables = dict(learned.tables)
    for seg in ("v", "j"):
        acol, choice_ev = f"{seg}_allele", f"{seg}_choice"
        g = learned.genomic[f"genes_{seg}"]
        cut = {r[acol]: r["cut_segment"] for r in g.iter_rows(named=True)}
        by_gene: dict[str, list[str]] = defaultdict(list)
        for r in g.filter(pl.col("functional")).iter_rows(named=True):
            by_gene[r["gene"]].append(r[acol])
        lmass = {r[acol]: r["p"] for r in tables[choice_ev].group_by(acol).agg(pl.col("p").sum()).iter_rows(named=True)}
        omass = {r[acol]: r["p"] for r in oracle.tables[choice_ev].group_by(acol).agg(pl.col("p").sum()).iter_rows(named=True)}
        oracle_donors = [a for a in cut if omass.get(a, 0) > 0 and cut.get(a)]
        if not oracle_donors:
            continue
        floor = min(p for p in omass.values() if p > 0) * 0.5      # usage for genes no reference models
        direct: set[str] = set()          # gene the oracle has -> copy its own alleles verbatim
        proxy: dict[str, str] = {}        # gene neither has -> rep allele <- germline-nearest oracle allele
        for alleles in by_gene.values():
            if any(lmass.get(a, 0) > 0 for a in alleles):          # already present in learned
                continue
            own = [a for a in alleles if omass.get(a, 0) > 0 and cut.get(a)]
            if own:
                direct |= set(own)
                continue
            rep = next((a for a in alleles if a.endswith("*01") and cut.get(a)),
                       next((a for a in alleles if cut.get(a)), None))
            if rep is not None:                                    # else: no germline (unscoreable ORF)
                proxy[rep] = max(oracle_donors, key=lambda b: SequenceMatcher(None, cut[rep], cut[b]).ratio())
        if not direct and not proxy:
            continue
        moved = list(direct) + list(proxy)
        for ev in list(tables):        # choice table + every child (deletion, P(J|V), P(D|J))
            if acol not in tables[ev].columns or ev not in oracle.tables or acol not in oracle.tables[ev].columns:
                continue
            parts = [tables[ev].filter(~pl.col(acol).is_in(moved))]
            if direct:                 # oracle carries the gene: copy its rows verbatim (own usage + profile)
                parts.append(oracle.tables[ev].filter(pl.col(acol).is_in(list(direct))))
            for rep, donor in proxy.items():   # neither has it: nearest oracle gene's profile, relabelled
                parts.append(oracle.tables[ev].filter(pl.col(acol) == donor).with_columns(pl.lit(rep).alias(acol)))
            tables[ev] = pl.concat(parts)
        nk = normalization_keys(learned.manifest.events[choice_ev])
        if proxy:                      # proxy genes inherited the donor's usage — override with the floor
            tables[choice_ev] = tables[choice_ev].with_columns(
                p=pl.when(pl.col(acol).is_in(list(proxy))).then(floor).otherwise(pl.col("p")))
        tot = pl.col("p").sum().over(nk) if nk else pl.col("p").sum()
        tables[choice_ev] = tables[choice_ev].with_columns(p=pl.when(tot > 0).then(pl.col("p") / tot).otherwise(0.0))
    return Model(manifest=learned.manifest, tables=tables, genomic=learned.genomic,
                 training=learned.training)


def _nearest_donor(seq: str, candidates: dict[str, str]) -> str | None:
    """The candidate allele whose germline is most similar to ``seq`` (``{allele: germline}``)."""
    from difflib import SequenceMatcher

    if not candidates or not seq:
        return None
    return max(candidates, key=lambda b: SequenceMatcher(None, seq, candidates[b]).ratio())


def _renormalize(table: pl.DataFrame, keys: list[str]) -> pl.DataFrame:
    tot = pl.col("p").sum().over(keys) if keys else pl.col("p").sum()
    return table.with_columns(p=pl.when(tot > 0).then(pl.col("p") / tot).otherwise(0.0))


def extend_alleles(model: Model, germline: pl.DataFrame, *, weight: float = 1.0) -> Model:
    """Add alleles from a larger germline library to an existing model, seeded from what it knows.

    The use case is a model fitted against one reference meeting a richer one — a newer IMGT
    release, a population-specific library, your own genotyped alleles. Every new allele needs a
    germline row, a choice probability and a full set of child conditionals (its deletion profile,
    and ``P(J|V)`` / ``P(D|J)`` where it is a parent), none of which the library supplies.

    Seeding uses the strongest evidence available for each case:

    - **A new allele of a gene the model already has.** Its choice mass is ``weight ×`` the mean
      mass of that gene's existing alleles, and its child tables are copied from a gene-mate. A new
      IMGT allele of a known gene is a polymorphism whose carriers use it about as often as the
      ``*01``, so the gene's own level is the right prior.
    - **A brand-new gene.** Child tables come from the germline-nearest existing allele, and the
      choice mass is a **floor** of half the smallest non-zero mass in the model. Plausible shape,
      deliberately tiny mass — there is no evidence at all for how often it is used.

    Deletion rows copied from a donor are clipped to the new allele's own germline length, so an
    extension can never introduce the unreachable mass
    :func:`~vdjtools.model.check.check_model` flags.

    Existing alleles are **never modified**, including their germline: silently swapping the
    sequence under an allele the model was fitted on would invalidate every conditional that
    references it. A library that disagrees about an existing allele is reported by
    ``check_model``'s ``germline_source`` check, not fixed here.

    Args:
        model: The model to extend.
        germline: A germline frame (see :func:`~vdjtools.model.io.from_germline` for the schema),
            typically a superset of the model's own.
        weight: Scales the seeded mass for new alleles of known genes. ``1.0`` gives a new allele
            the gene's average; ``0.5`` is a more conservative half of it.

    Returns:
        A new, validated and renormalized :class:`Model`. Idempotent — extending twice with the
        same library changes nothing the second time.

    Note:
        This *seeds*, it does not estimate. Follow it with
        ``infer_native(extended, seqs, init="template")`` to let data set the new probabilities.

    Example:
        >>> bigger = extend_alleles(m, load_germline("TRB", "human"))
        >>> bigger, rep = infer_native(bigger, seqs, init="template", max_iter=5)
    """
    from . import reference as ref
    from .io import _d_genomic, _vj_genomic

    gl = ref.normalize_germline(germline)
    pal = model.manifest.palindrome_max
    tables = dict(model.tables)
    genomic = dict(model.genomic)

    for seg in ("v", "j", "d"):
        frame_name = f"genes_{seg}"
        if frame_name not in genomic:
            continue
        acol, choice_ev = f"{seg}_allele", {"v": "v_choice", "j": "j_choice", "d": "d_gene"}[seg]
        if choice_ev not in tables:
            continue
        existing = genomic[frame_name]
        known = set(existing[acol].to_list())
        add = gl.filter((pl.col("segment") == seg.upper()) & ~pl.col("allele").is_in(list(known)))
        if not add.height:
            continue

        new_rows = (_vj_genomic(add, seg, pal[f"{seg}_3" if seg == "v" else "j_5"])
                    if seg in ("v", "j") else _d_genomic(add, pal["d_5"], pal["d_3"]))
        if not new_rows.height:
            continue
        genomic[frame_name] = pl.concat([existing, new_rows.select(existing.columns)])

        cut = {r[acol]: r["cut_segment"] for r in existing.iter_rows(named=True) if r["cut_segment"]}
        mass = {r[acol]: r["p"] for r in
                tables[choice_ev].group_by(acol).agg(pl.col("p").sum()).iter_rows(named=True)}
        by_gene: dict[str, list[str]] = defaultdict(list)
        for allele in cut:
            by_gene[allele.split("*")[0]].append(allele)
        # Per-ROW floor, not per-allele: on a conditioned choice table (P(J|V), P(D|J)) every
        # parent group sums to 1, so a single row's probability is the comparable scale. Summing
        # over parents would be n_parents too large.
        row_p = tables[choice_ev].filter(pl.col("p") > 0)["p"]
        floor = float(row_p.min()) * 0.5 if row_p.len() else 1e-6
        # Snapshot the per-gene totals BEFORE any donor rows are copied in, so step (2) below
        # restores what the model actually had rather than what the copy just inflated it to.
        gene_expr = pl.col(acol).str.split("*").list.first()
        nk = normalization_keys(model.manifest.events[choice_ev])
        old_totals = (tables[choice_ev].with_columns(_g=gene_expr)
                      .group_by([*nk, "_g"]).agg(pl.col("p").sum().alias("_old")))

        donors: dict[str, str] = {}      # new allele -> existing allele to copy conditionals from
        new_genes: set[str] = set()      # new alleles whose GENE is also new to the model
        for r in new_rows.iter_rows(named=True):
            allele, gene, seq = r[acol], r["gene"], r["cut_segment"]
            mates = [a for a in by_gene.get(gene, []) if mass.get(a, 0) > 0]
            if mates:
                donors[allele] = mates[0]
            else:
                nearest = _nearest_donor(seq, {a: s for a, s in cut.items() if mass.get(a, 0) > 0})
                if nearest is None:
                    continue
                donors[allele] = nearest
                new_genes.add(allele)
        if not donors:
            continue

        new_cut_len = {r[acol]: len(r["cut_segment"]) for r in new_rows.iter_rows(named=True)}
        for ev_name in list(tables):
            if acol not in tables[ev_name].columns:
                continue
            parts = [tables[ev_name]]
            for allele, donor in donors.items():
                block = (tables[ev_name].filter(pl.col(acol) == donor)
                         .with_columns(pl.lit(allele).alias(acol)))
                if not block.height:
                    continue
                block = _clip_deletions(block, ev_name, model.manifest.events.get(ev_name),
                                        new_cut_len[allele], pal)
                parts.append(block)
            tables[ev_name] = pl.concat(parts)

        # Each new allele's choice rows were copied verbatim from its donor above, which already
        # puts them on the right per-parent scale. Two corrections follow.
        seeded = tables[choice_ev].with_columns(_g=gene_expr)
        # (1) A brand-new GENE has no evidence for its usage at all, so it gets a floor rather
        #     than its donor's real mass -- a plausible shape, a deliberately tiny weight.
        #     `weight` scales a new allele of a KNOWN gene relative to its gene-mate.
        seeded = seeded.with_columns(
            p=pl.when(pl.col(acol).is_in(list(new_genes))).then(pl.lit(floor))
            .when(pl.col(acol).is_in(list(donors))).then(pl.col("p") * weight)
            .otherwise(pl.col("p")))
        # (2) Preserve each pre-existing GENE's total usage. Alleles of one gene are alternative
        #     versions of the same gene, not extra genes -- a diploid carries at most two -- so a
        #     richer library must SPLIT a gene's mass more finely, never multiply it. Without this,
        #     extending human TRB from 1 to ~3 alleles per gene inflated gene-level V usage by up
        #     to 6 points, silently reweighting every Pgen through those genes.
        seeded = (seeded.join(old_totals, on=[*nk, "_g"], how="left")
                  .with_columns(_new=pl.col("p").sum().over([*nk, "_g"])))
        seeded = seeded.with_columns(
            p=pl.when((pl.col("_new") > 0) & (pl.col("_old") > 0))
            .then(pl.col("p") * pl.col("_old") / pl.col("_new"))
            .otherwise(pl.col("p"))
        ).drop(["_g", "_old", "_new"])
        tables[choice_ev] = _renormalize(seeded, nk)
        # A donor's child tables were copied wholesale, so each new allele's own conditionals must
        # be re-normalized within their own group (clipping deletions removed some of their mass).
        for ev_name, ev_obj in model.manifest.events.items():
            if ev_name == choice_ev or acol not in tables[ev_name].columns:
                continue
            nk = normalization_keys(ev_obj)
            if acol in nk:
                tables[ev_name] = _renormalize(tables[ev_name], nk)

    return Model(manifest=model.manifest, tables=tables, genomic=genomic,
                 training=model.training).validate()


#: Deletion event -> the ``palindrome_max`` keys bounding it (mirrors ``check._DELETION_ENDS``).
_DEL_ENDS = {"v_3_del": ("v_3",), "j_5_del": ("j_5",),
             "d_del": ("d_5", "d_3"), "d2_del": ("d_5", "d_3")}


def _clip_deletions(block: pl.DataFrame, ev_name: str, event, cut_len: int,
                    pal: dict) -> pl.DataFrame:
    """Drop deletion rows a new allele's own germline cannot reach (see ``check.check_model``).

    Reachability matches the Pgen DP: ``ndel <= len(cut) - Σ palindrome_max - 1`` for V/J (which
    must each leave one nt) and ``ndel5 + ndel3 <= len(cut) - Σ palindrome_max`` for D.
    """
    ends = _DEL_ENDS.get(ev_name)
    if event is None or ends is None:
        return block
    pal_total = sum(pal.get(e, 0) for e in ends)
    if event.kind.value == "deletion" and "ndel" in block.columns:
        return block.filter(pl.col("ndel") <= cut_len - pal_total - 1)
    if event.kind.value == "deletion_2d" and "ndel5" in block.columns:
        return block.filter(pl.col("ndel5") + pl.col("ndel3") <= cut_len - pal_total)
    return block


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
    progress=None,
    checkpoint=None,
    checkpoint_every: int = 1,
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
    import time as _time

    from .._core import estep_batch, make_counts
    from .native import _encode, pack

    started = _time.monotonic()
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
    meta = dict(started=started, max_iter=max_iter, tol=tol, init=init, native=True,
                n_sequences=len(upper))
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
        report.n_iter = it + 1
        model = new_model
        # Converge on the relative log-likelihood improvement (whole-model, monotone post-fix), not V
        # usage alone: V is arda-masked so it settles in ~2 iters while trims/insertions/D-D still move.
        rel = _loglik_rel(report.loglik)
        report.gene_tv.append(rel)
        if progress is not None:
            progress(report.n_iter, report.loglik[-1], rel, report.n_scoreable[-1])
        if checkpoint is not None and report.n_iter % checkpoint_every == 0:
            _write_checkpoint(checkpoint, model, template, report, meta)
        if it > 0 and rel < tol:
            report.converged = True
            break
    # Completion pass: EM only keeps genes whose CDR3s were producibly observed, so functional genes
    # arda never saw (or saw but couldn't emit) end at P=0 — yet a user's library may be full of them.
    # ``template`` is the OLGA oracle (from_olga seed), differing only in D/D-D, so transplant its own
    # usage + conditionals for every functional gene the learned model is missing. Only in the prior
    # mode (``gene_prior > 0``), i.e. "keep every real gene reachable".
    if gene_prior > 0:
        model = augment_from_oracle(model, template)
    model = _record(report, model, template, **meta)
    return model, report
