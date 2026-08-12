"""Score sequences under a model: likelihood, BIC, Pgen distributions, entropy and diversity.

Everything here is **data-conditioned** — it needs sequences, or sequences drawn from the model —
which is what separates it from :mod:`vdjtools.model.analyze` (model-only information theory).

Three questions this answers:

- **How well does this model explain these sequences?** :func:`model_fit` — log-likelihood, AIC, BIC.
- **Do two models score the same repertoire the same way?** :func:`compare_pgen` +
  :func:`pgen_summary`.
- **How much diversity does this model actually generate?** :func:`diversity` — the entropy of the
  generated sequence distribution and the ``~10^x`` effective-diversity number that follows from it.

Two conventions run through the module and are worth stating once.

**Likelihoods use nucleotide Pgen.** ``Σ Pgen_nt`` over all nt CDR3s is 1, so ``log Pgen_nt`` is a
proper log-likelihood and BIC is meaningful. ``Pgen_aa`` sums only the in-frame, stop-free
nucleotide fiber of a translation, so ``Σ Pgen_aa < 1``: an amino-acid log-likelihood is
unnormalized, and the missing constant **differs between models**. Amino-acid scoring is supported
(real clonotype tables are often aa-only) but is a relative score on one fixed sequence set, never
an absolute one.

**Pgen = 0 never becomes ``-inf``.** A sequence the model cannot generate gets ``pgen = 0.0`` and a
null log-probability, and aggregates are taken over the scoreable subset — the same convention
:class:`~vdjtools.model.infer.InferenceReport` already uses. Every aggregate reports ``n`` beside
``n_scoreable`` so a flattering log-likelihood earned on 10% of the data is visible rather than
hidden.
"""
from __future__ import annotations

import math

import numpy as np
import polars as pl

from . import native
from .model import Model
from .schema import _allele_col, normalization_keys

_LOG10 = math.log(10.0)
_LOG2 = math.log(2.0)


def _is_nt(seq: str) -> bool:
    return bool(seq) and set(seq.upper()) <= set("ACGT")


def _resolve_kind(seqs: list[str], kind: str) -> str:
    """Decide nt vs aa for the whole set; ``"auto"`` requires the set to be homogeneous."""
    if kind in ("nt", "aa"):
        return kind
    if kind != "auto":
        raise ValueError(f"kind must be 'nt', 'aa' or 'auto', got {kind!r}")
    if not seqs:
        return "aa"
    nt = [_is_nt(s) for s in seqs]
    if all(nt):
        return "nt"
    if not any(nt):
        return "aa"
    raise ValueError(
        "kind='auto' cannot type a mixed set of nucleotide and amino-acid sequences; "
        "pass kind='nt' or kind='aa' explicitly"
    )


def _allele_index(model: Model, seg: str) -> tuple[set[str], dict[str, list[str]]]:
    """(all alleles, gene -> its alleles) for one segment, from the model's germline frame."""
    g = model.genomic.get(f"genes_{seg}")
    if g is None:
        return set(), {}
    alleles = g[f"{seg}_allele"].to_list()
    by_gene: dict[str, list[str]] = {}
    for a in alleles:
        by_gene.setdefault(a.split("*")[0], []).append(a)
    return set(alleles), by_gene


def _resolve_call(call, alleles: set[str], by_gene: dict[str, list[str]], seg: str,
                  on_unknown: str) -> str | None:
    """Map an AIRR V/J call onto a model allele, or ``None`` to marginalize.

    Real repertoires carry **gene-level** calls (``TRBV9``) while the model is keyed by allele, and
    ``native._gene_idx`` raises rather than silently returning the V/J-agnostic Pgen (which was
    2.38x too high — see the trap note in CLAUDE.md). So resolve here, explicitly: a gene with
    exactly one allele in the model (the ``collapse=True`` default of ``load_bundled``) resolves;
    anything ambiguous or unknown obeys ``on_unknown``.
    """
    if call is None or call == "" or not alleles:
        return None
    candidates: set[str] = set()
    for part in str(call).split(","):
        part = part.strip()
        if not part:
            continue
        if part in alleles:
            candidates.add(part)
        else:
            candidates.update(by_gene.get(part.split("*")[0], []))
    if len(candidates) == 1:
        return candidates.pop()
    if on_unknown == "marginalize":
        return None
    reason = "is ambiguous in" if candidates else "is not in"
    raise KeyError(
        f"{seg} call {call!r} {reason} this model's allele set "
        f"({len(candidates)} match(es)); pass on_unknown='marginalize' to score it "
        f"{seg}-agnostically, or use a model whose alleles match your calls"
    )


def _unpack(sequences, v, j, seq_col, v_col, j_col) -> tuple[list[str], list, list]:
    """Accept either parallel lists or a clonotype frame; return (seqs, v_calls, j_calls)."""
    if isinstance(sequences, (pl.DataFrame, pl.LazyFrame)):
        df = sequences.collect() if isinstance(sequences, pl.LazyFrame) else sequences
        col = seq_col or next(
            (c for c in ("junction", "junction_nt", "junction_aa", "cdr3_nt", "cdr3_aa", "sequence")
             if c in df.columns), None)
        if col is None:
            raise ValueError(
                f"no sequence column found in {df.columns}; pass seq_col=")
        seqs = df[col].to_list()
        vs = df[v_col].to_list() if v_col and v_col in df.columns else (
            df["v_call"].to_list() if "v_call" in df.columns else [None] * len(seqs))
        js = df[j_col].to_list() if j_col and j_col in df.columns else (
            df["j_call"].to_list() if "j_call" in df.columns else [None] * len(seqs))
        return seqs, vs, js
    seqs = list(sequences)
    vs = list(v) if v is not None else [None] * len(seqs)
    js = list(j) if j is not None else [None] * len(seqs)
    if len(vs) != len(seqs) or len(js) != len(seqs):
        raise ValueError("v and j must have the same length as sequences")
    return seqs, vs, js


#: Below this many sequences, thread startup costs more than the nt Pgen it saves.
_NT_THREAD_MIN = 32


def _pgen_nt_many(model: Model, seqs: list[str], vres: list, jres: list,
                  threads: int) -> list[float]:
    """Nucleotide Pgen over many sequences, threaded.

    There is no native nt batch entry point, but the ``pgen_nt`` binding releases the GIL, so a
    plain thread pool gets the same parallelism. Worth it: a V/J-marginalized nt Pgen sums over
    every V and J and costs tens of milliseconds, which is minutes for a Monte-Carlo diversity
    estimate on one core.
    """
    if threads == 1 or len(seqs) < _NT_THREAD_MIN:
        return [native.pgen_nt(model, s, a, b) for s, a, b in zip(seqs, vres, jres)]
    import os
    from concurrent.futures import ThreadPoolExecutor

    native.pack(model)  # populate the pack cache once, before the workers race for it
    n_workers = threads if threads > 0 else max(1, (os.cpu_count() or 2) - 2)
    with ThreadPoolExecutor(max_workers=n_workers) as pool:
        return list(pool.map(lambda t: native.pgen_nt(model, t[0], t[1], t[2]),
                             zip(seqs, vres, jres)))


def pgen_frame(model: Model, sequences, *, v=None, j=None, kind: str = "auto",
               use_calls: bool = True, on_unknown: str = "error", threads: int = 0,
               seq_col: str | None = None, v_col: str | None = None,
               j_col: str | None = None) -> pl.DataFrame:
    """Generation probability of each sequence under ``model``, as a tidy frame.

    Args:
        model: The recombination model to score with.
        sequences: A list of junction/CDR3 strings, **or** a clonotype ``pl.DataFrame`` (the
            sequence column and ``v_call``/``j_call`` are auto-detected).
        v: Optional per-sequence V calls, when ``sequences`` is a list.
        j: Optional per-sequence J calls, when ``sequences`` is a list.
        kind: ``"nt"``, ``"aa"``, or ``"auto"`` (default) to detect from the sequences. ``"auto"``
            requires a homogeneous set.
        use_calls: Condition each sequence's Pgen on its own V/J call. ``False`` marginalizes over
            all V/J — note this changes the quantity from ``P(junction, V, J)`` to ``P(junction)``.
        on_unknown: What to do with a call that does not resolve to exactly one model allele —
            ``"error"`` (default) or ``"marginalize"``.
        threads: Worker threads for amino-acid scoring (``0`` = auto). Nucleotide scoring is serial.
        seq_col: Explicit sequence column, when ``sequences`` is a frame.
        v_col: Explicit V-call column, when ``sequences`` is a frame.
        j_col: Explicit J-call column, when ``sequences`` is a frame.

    Returns:
        ``sequence, v_call, j_call, kind, pgen, log_pgen, log10_pgen, scoreable`` — one row per
        input sequence, in input order. ``log_pgen`` is null where ``pgen`` is 0.

    Raises:
        KeyError: If a V/J call does not resolve and ``on_unknown="error"``.
        ValueError: On a mixed nt/aa set under ``kind="auto"``, or a missing sequence column.
    """
    seqs, vs, js = _unpack(sequences, v, j, seq_col, v_col, j_col)
    kind = _resolve_kind(seqs, kind)

    if use_calls:
        v_alleles, v_by_gene = _allele_index(model, "v")
        j_alleles, j_by_gene = _allele_index(model, "j")
        vres = [_resolve_call(c, v_alleles, v_by_gene, "V", on_unknown) for c in vs]
        jres = [_resolve_call(c, j_alleles, j_by_gene, "J", on_unknown) for c in js]
    else:
        vres = jres = [None] * len(seqs)

    if kind == "aa":
        pgens = native.pgen_aa_batch(model, seqs, v=vres, j=jres, threads=threads)
    else:
        pgens = _pgen_nt_many(model, seqs, vres, jres, threads)

    p = np.asarray(pgens, dtype=float)
    ok = p > 0
    logp = np.full(len(p), np.nan)
    logp[ok] = np.log(p[ok])
    return pl.DataFrame({
        "sequence": seqs,
        "v_call": vres,
        "j_call": jres,
        "kind": [kind] * len(seqs),
        "pgen": p,
        "log_pgen": logp,
        "log10_pgen": logp / _LOG10,
        "scoreable": ok,
    }).with_columns(
        pl.when(pl.col("scoreable")).then(pl.col("log_pgen")).otherwise(None).alias("log_pgen"),
        pl.when(pl.col("scoreable")).then(pl.col("log10_pgen")).otherwise(None).alias("log10_pgen"),
    )


def free_params(model: Model, *, by_event: bool = False, eps: float = 0.0, tol: float = 1e-9,
                reachable_only: bool = True) -> int | pl.DataFrame:
    """Number of free parameters ``k`` in a model — the penalty term of AIC / BIC.

    ``k = Σ_events Σ_groups max(support − 1, 0)``, where a *group* is one normalization group
    (:func:`~vdjtools.model.schema.normalization_keys`, which already splits a dinucleotide table
    by ``from_nt``, giving the right ``4 × (4 − 1) = 12``), the ``−1`` is the simplex constraint,
    and **support is the number of cells carrying probability**, not the number of rows.

    Counting rows instead of support would be badly wrong: deletion bins past a germline's length,
    alleles pinned to zero and insertion tail bins are **structural** zeros that parameterize
    nothing. On human TRB that is the difference between ~3,600 and ~700 parameters for ``v_3_del``
    alone.

    Two group kinds contribute nothing and are dropped:

    - **undefined conditionals** — the all-zero groups the schema explicitly permits, kept only for
      gene-index alignment;
    - **unreachable conditionals** — a group whose parent allele has zero marginal probability
      (``from_olga`` fills ``P(D|J)`` uniformly even where ``P(J) = 0``). These sum to 1 but are
      not estimable, and counting them inflates ``k`` by ``n_D − 1`` for every dead J.

    Caveat, and it is a real one: a support-based count cannot tell a structural zero from a
    parameter EM happened to drive to exactly zero, so ``k`` is a **lower bound**. BIC is therefore
    only comparable between models counted the same way — which is the case for any two models
    compared through this function.

    Args:
        model: The model to count.
        by_event: Return a per-event breakdown instead of the total.
        eps: A cell counts toward the support when ``p > eps``.
        tol: A group whose probabilities sum to at most this is treated as an undefined conditional.
        reachable_only: Drop groups whose parent allele has zero marginal probability.

    Returns:
        ``k`` as an ``int``, or (with ``by_event``) a frame ``event, n_groups, k``.
    """
    from .analyze import gene_marginal

    rows = []
    for name, event in model.manifest.events.items():
        df = model.tables[name]
        keys = normalization_keys(event)
        if keys:
            agg = df.group_by(keys).agg(total=pl.col("p").sum(),
                                        support=(pl.col("p") > eps).sum())
        else:
            agg = df.select(total=pl.col("p").sum(), support=(pl.col("p") > eps).sum())
        agg = agg.filter(pl.col("total") > tol)
        if reachable_only and len(event.given) == 1:
            parent_col = _allele_col(model.manifest.events[event.given[0]])
            if parent_col in agg.columns:
                marg = gene_marginal(model, parent_col.split("_")[0])
                agg = agg.filter(
                    pl.col(parent_col).replace_strict(marg, default=0.0,
                                                      return_dtype=pl.Float64) > 0)
        k = agg.select((pl.col("support") - 1).clip(lower_bound=0).sum()).item() if agg.height else 0
        rows.append({"event": name, "n_groups": agg.height, "k": int(k)})

    if by_event:
        return pl.DataFrame(rows, schema={"event": pl.Utf8, "n_groups": pl.Int64, "k": pl.Int64})
    return int(sum(r["k"] for r in rows))


def model_fit(model: Model, sequences, *, weights=None, k: int | None = None,
              **kw) -> pl.DataFrame:
    """Log-likelihood, AIC and BIC of a sequence set under a model.

    Args:
        model: The model to score with.
        sequences: Sequences or a clonotype frame — see :func:`pgen_frame`.
        weights: Optional per-sequence weights (e.g. ``duplicate_count``), or the name of a column
            when ``sequences`` is a frame. The log-likelihood becomes ``Σ w·log p`` and ``n``
            becomes ``Σ w``, i.e. exactly as if each clonotype were repeated ``w`` times.
        k: Override the free-parameter count; defaults to :func:`free_params`.
        **kw: Passed to :func:`pgen_frame` (``kind``, ``use_calls``, ``on_unknown``, ``threads``,
            ``seq_col``, ``v_col``, ``j_col``, ``v``, ``j``).

    Returns:
        A one-row frame: ``kind, conditioned, n, n_scoreable, frac_scoreable, loglik_sum,
        loglik_mean, k, aic, bic``. ``conditioned`` records whether V/J were conditioned on, because
        that changes the sample space (``P(junction, V, J)`` vs ``P(junction)``) and two fits are
        only comparable when it matches.

    Example:
        >>> model_fit(load_bundled("TRB", "learned"), held_out_junctions)
    """
    if isinstance(weights, str):
        if not isinstance(sequences, pl.DataFrame):
            raise ValueError("a string `weights` names a column, so `sequences` must be a frame")
        w = sequences[weights].to_numpy().astype(float)
    elif weights is not None:
        w = np.asarray(weights, dtype=float)
    else:
        w = None

    scored = pgen_frame(model, sequences, **kw)
    ok = scored["scoreable"].to_numpy()
    logp = scored["log_pgen"].fill_null(0.0).to_numpy()
    if w is None:
        w = np.ones(len(logp))
    if len(w) != len(logp):
        raise ValueError("weights must have the same length as sequences")

    n_total = float(w.sum())
    n_ok = float(w[ok].sum())
    ll = float((w[ok] * logp[ok]).sum())
    kk = free_params(model) if k is None else int(k)
    # n_ok, not n_total: the likelihood is only over what the model could score, so the BIC sample
    # size must match. frac_scoreable sits next to it so the trade-off is never invisible.
    bic = kk * math.log(n_ok) - 2 * ll if n_ok > 0 else float("nan")
    return pl.DataFrame([{
        "kind": scored["kind"][0] if scored.height else None,
        "conditioned": bool(kw.get("use_calls", True)),
        "n": n_total,
        "n_scoreable": n_ok,
        "frac_scoreable": n_ok / n_total if n_total else float("nan"),
        "loglik_sum": ll,
        "loglik_mean": ll / n_ok if n_ok else float("nan"),
        "k": kk,
        "aic": 2 * kk - 2 * ll,
        "bic": bic,
    }])


def compare_pgen(a: Model, b: Model, sequences, *, labels: tuple[str, str] = ("a", "b"),
                 **kw) -> pl.DataFrame:
    """Score one sequence set under two models and pair the results up.

    Args:
        a: First model.
        b: Second model.
        sequences: Sequences or a clonotype frame — see :func:`pgen_frame`.
        labels: Names for the two models, used in the output column suffixes.
        **kw: Passed to :func:`pgen_frame` for both models, so the two are scored identically.

    Returns:
        ``sequence, v_call_<a>, v_call_<b>, j_call_<a>, j_call_<b>, kind, pgen_<a>, pgen_<b>,
        log10_<a>, log10_<b>, delta_log10`` — one row per sequence. The two V/J columns differ when
        the models resolve a gene-level call to different alleles. ``delta_log10`` is null unless
        both models scored the sequence.

    Example:
        >>> summary = pgen_summary(compare_pgen(olga, learned, seqs, labels=("olga", "learned")))
    """
    la, lb = labels
    if la == lb:
        raise ValueError(f"labels must differ, got {labels!r}")
    fa = pgen_frame(a, sequences, **kw)
    fb = pgen_frame(b, sequences, **kw)
    return pl.DataFrame({
        "sequence": fa["sequence"],
        f"v_call_{la}": fa["v_call"], f"v_call_{lb}": fb["v_call"],
        f"j_call_{la}": fa["j_call"], f"j_call_{lb}": fb["j_call"],
        "kind": fa["kind"],
        f"pgen_{la}": fa["pgen"], f"pgen_{lb}": fb["pgen"],
        f"log10_{la}": fa["log10_pgen"], f"log10_{lb}": fb["log10_pgen"],
    }).with_columns(
        delta_log10=pl.col(f"log10_{la}") - pl.col(f"log10_{lb}")
    )


def pgen_summary(cmp: pl.DataFrame, *, labels: tuple[str, str] = ("a", "b")) -> pl.DataFrame:
    """Summarise a :func:`compare_pgen` frame: agreement, offset, and coverage.

    Args:
        cmp: The output of :func:`compare_pgen`.
        labels: The same labels that call used.

    Returns:
        A one-row frame: ``n, n_scoreable_a, n_scoreable_b, n_scoreable_both, only_a_scoreable,
        only_b_scoreable, mean_log10_a, mean_log10_b, median_log10_a, median_log10_b, mean_delta,
        median_delta, sd_delta, q05_delta, q95_delta, pearson_log10, spearman_log10, ks_stat,
        ks_p``.

    Note:
        ``only_a_scoreable`` / ``only_b_scoreable`` are the headline numbers, not the correlations:
        one model assigning Pgen 0 to thousands of sequences the other scores fine is the finding,
        and a mean-delta-only report would hide it entirely. The KS statistic compares the two
        **marginal** log10 distributions; the samples are paired, so its p-value is anticonservative
        — lean on the statistic.
    """
    from scipy.stats import ks_2samp, spearmanr

    la, lb = labels
    a = cmp[f"log10_{la}"].to_numpy()
    b = cmp[f"log10_{lb}"].to_numpy()
    ok_a, ok_b = ~np.isnan(a), ~np.isnan(b)
    both = ok_a & ok_b
    d = cmp["delta_log10"].drop_nulls().to_numpy()

    ks_stat = ks_p = float("nan")
    if ok_a.sum() and ok_b.sum():
        ks = ks_2samp(a[ok_a], b[ok_b])
        ks_stat, ks_p = float(ks.statistic), float(ks.pvalue)
    pearson = spearman = float("nan")
    if both.sum() > 1:
        av, bv = a[both], b[both]
        if av.std() > 0 and bv.std() > 0:
            pearson = float(np.corrcoef(av, bv)[0, 1])
            spearman = float(spearmanr(av, bv).statistic)

    return pl.DataFrame([{
        "n": cmp.height,
        "n_scoreable_a": int(ok_a.sum()), "n_scoreable_b": int(ok_b.sum()),
        "n_scoreable_both": int(both.sum()),
        "only_a_scoreable": int((ok_a & ~ok_b).sum()),
        "only_b_scoreable": int((ok_b & ~ok_a).sum()),
        "mean_log10_a": float(np.mean(a[ok_a])) if ok_a.any() else float("nan"),
        "mean_log10_b": float(np.mean(b[ok_b])) if ok_b.any() else float("nan"),
        "median_log10_a": float(np.median(a[ok_a])) if ok_a.any() else float("nan"),
        "median_log10_b": float(np.median(b[ok_b])) if ok_b.any() else float("nan"),
        "mean_delta": float(np.mean(d)) if d.size else float("nan"),
        "median_delta": float(np.median(d)) if d.size else float("nan"),
        "sd_delta": float(np.std(d, ddof=1)) if d.size > 1 else float("nan"),
        "q05_delta": float(np.quantile(d, 0.05)) if d.size else float("nan"),
        "q95_delta": float(np.quantile(d, 0.95)) if d.size else float("nan"),
        "pearson_log10": pearson, "spearman_log10": spearman,
        "ks_stat": ks_stat, "ks_p": ks_p,
    }])


# --- information content and diversity ---------------------------------------------------------

def pgen_spectrum(model: Model, *, n: int = 10_000, seed: int = 0, bins: int = 40,
                  productive_only: bool = False, sequences=None, **kw) -> pl.DataFrame:
    """The model's Pgen distribution, as a histogram table ready to plot.

    Args:
        model: The model whose Pgen spectrum to take.
        n: Sequences to generate when ``sequences`` is not given.
        seed: Generation seed.
        bins: Number of equal-width log10 bins.
        productive_only: Restrict generation to productive rearrangements.
        sequences: Score these instead of generating (e.g. a real repertoire), so an observed
            spectrum can be overlaid on the model's own.
        **kw: Passed to :func:`pgen_frame`.

    Returns:
        ``bin_left, bin_right, bin_mid, count, frac`` over ``log10(Pgen)``.
    """
    scored = _model_pgens(model, n=n, seed=seed, productive_only=productive_only,
                          sequences=sequences, **kw)
    x = scored["log10_pgen"].drop_nulls().to_numpy()
    if not x.size:
        raise ValueError("no scoreable sequences — cannot build a Pgen spectrum")
    counts, edges = np.histogram(x, bins=bins)
    return pl.DataFrame({
        "bin_left": edges[:-1], "bin_right": edges[1:],
        "bin_mid": (edges[:-1] + edges[1:]) / 2.0,
        "count": counts.astype(np.int64),
        "frac": counts / counts.sum(),
    })


def _model_pgens(model: Model, *, n: int, seed: int, productive_only: bool, sequences=None,
                 **kw) -> pl.DataFrame:
    """Pgen of ``n`` sequences drawn from the model (or of a supplied set)."""
    if sequences is not None:
        return pgen_frame(model, sequences, **kw)
    from .generate import generate

    gen = generate(model, n, seed=seed, productive_only=productive_only)
    # Score nt, unconditioned: the Monte-Carlo entropy estimator below needs P(junction), the
    # same quantity the sequences were drawn from.
    kw.setdefault("use_calls", False)
    kw.setdefault("kind", "nt")
    return pgen_frame(model, gen["junction_nt"].to_list(), **kw)


def diversity(model: Model, *, n: int = 5_000, seed: int = 0,
              productive_only: bool = False) -> pl.DataFrame:
    """How much diversity this model generates — entropy in bits, and effective diversity.

    Two independent readings of "total diversity", reported side by side because they answer
    different questions and differ by orders of magnitude:

    - **Scenario entropy** ``H_scenario`` — the information in one *recombination event*, summed
      over the Bayes net (:func:`~vdjtools.model.analyze.total_entropy`). This is an upper bound on
      the sequence entropy, because different scenarios can produce the same junction.
    - **Sequence entropy** ``H_sequence`` — the entropy of the junction distribution itself,
      estimated by Monte Carlo. Sequences drawn from the model *are* distributed as ``Pgen``, so
      ``E[−log₂ Pgen]`` over generated sequences is an unbiased estimator of ``H`` and the sample
      standard error comes free.

    From those, two Hill numbers:

    - ``diversity_shannon = 2^H_sequence`` (Hill ``q = 1``) — the classic "effective number of
      distinct sequences" figure, the one usually quoted as ``~10^x`` for a locus.
    - ``diversity_simpson = 1 / E[Pgen]`` (Hill ``q = 2``) — the inverse coincidence probability,
      i.e. how many sequences you would need for two independent draws to collide. Exact, because
      ``Σ Pgen² = E_{s∼Pgen}[Pgen]``, so the generated sample estimates it directly. It is always
      the smaller number: it weights the common sequences more heavily.

    Args:
        model: The model to characterise.
        n: Sequences to generate. The Shannon estimate converges quickly; ``diversity_simpson`` is
            driven by the few highest-Pgen draws and needs more samples for a tight interval —
            check ``pgen_mean_se`` against ``pgen_mean``.
        seed: Generation seed, so the estimate is reproducible.
        productive_only: Estimate over productive rearrangements only. Off by default: the
            unrestricted distribution is the one ``Pgen`` normalizes over, so the estimator is
            unbiased there and merely descriptive here.

    Returns:
        A one-row frame: ``n, scenario_entropy_bits, scenario_diversity, sequence_entropy_bits,
        sequence_entropy_se_bits, diversity_shannon, diversity_simpson, pgen_mean, pgen_mean_se,
        median_log10_pgen, q05_log10_pgen, q95_log10_pgen``.

    Example:
        >>> diversity(load_bundled("TRB", "olga"), n=20_000)
    """
    from .analyze import total_entropy

    h_scenario = float(total_entropy(model)["contribution_bits"].sum())
    scored = _model_pgens(model, n=n, seed=seed, productive_only=productive_only)
    p = scored["pgen"].to_numpy()
    p = p[p > 0]
    if p.size < 2:
        raise ValueError("fewer than two scoreable generated sequences — cannot estimate diversity")

    log2p = np.log(p) / _LOG2
    h_seq = float(-log2p.mean())
    h_se = float(log2p.std(ddof=1) / math.sqrt(p.size))
    pmean = float(p.mean())
    return pl.DataFrame([{
        "n": int(p.size),
        "scenario_entropy_bits": h_scenario,
        "scenario_diversity": 2.0 ** h_scenario,
        "sequence_entropy_bits": h_seq,
        "sequence_entropy_se_bits": h_se,
        "diversity_shannon": 2.0 ** h_seq,
        "diversity_simpson": 1.0 / pmean,
        "pgen_mean": pmean,
        "pgen_mean_se": float(p.std(ddof=1) / math.sqrt(p.size)),
        "median_log10_pgen": float(np.median(np.log(p) / _LOG10)),
        "q05_log10_pgen": float(np.quantile(np.log(p) / _LOG10, 0.05)),
        "q95_log10_pgen": float(np.quantile(np.log(p) / _LOG10, 0.95)),
    }])
