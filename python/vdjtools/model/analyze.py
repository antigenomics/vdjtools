"""Information-theoretic diagnostics for a recombination :class:`~vdjtools.model.model.Model`.

Turns a model's declared Bayes net (:mod:`~vdjtools.model.events`) and its marginal tables into
three views used for validation and for the appendix figures:

- :func:`entropy_table` — per-event Shannon entropy of that part of the rearrangement (bits):
  the marginal entropy ``H(X)`` of the event's realization and, where the event is conditioned,
  the expected conditional entropy ``H(X | parents)``.
- :func:`mutual_information` — the information each declared edge carries, ``I(child; parent)`` =
  ``H(child) − H(child | parent)`` (bits), plus ``I(V; J)`` and the within-D ``I(delD5; delD3)``.
- :func:`bayes_net_dot` / :func:`render_bayes_net` — a graphviz DAG (bnlearn style) with nodes
  annotated by ``H(X)`` and edges by ``I``; rendered to PDF/PNG via the ``dot`` CLI.

Everything is read straight from the polars tables, so it works identically on a legacy OLGA
bootstrap model and an EM-inferred native model — the two are directly comparable.
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import numpy as np
import polars as pl

from .model import Model
from .schema import _allele_col, normalization_keys

_LOG2 = np.log(2.0)


def _H(p: np.ndarray) -> float:
    """Shannon entropy (bits) of a probability vector; unnormalized-safe, 0·log0 = 0."""
    p = np.asarray(p, dtype=float)
    s = p.sum()
    if s <= 0:
        return 0.0
    p = p[p > 0] / s
    return float(-(p * np.log(p)).sum() / _LOG2) + 0.0  # +0.0 folds -0.0 → 0.0


def gene_marginal(model: Model, seg: str) -> dict[str, float]:
    """Marginal ``P(seg)`` as {allele: prob}, forward-propagated over the Bayes net.

    V is always a root. VDJ J is a root (``P(J)``); VJ J is ``P(J|V)`` marginalized over V.
    D is ``P(D|J)`` marginalized over J.
    """
    t = model.tables
    if seg == "v":
        return dict(zip(t["v_choice"]["v_allele"], t["v_choice"]["p"]))
    if seg == "j":
        jt = t["j_choice"]
        if "v_allele" not in jt.columns:  # VDJ: root marginal
            return dict(zip(jt["j_allele"], jt["p"]))
        pv = gene_marginal(model, "v")  # VJ: Σ_v P(v) P(j|v)
        out: dict[str, float] = {}
        for v, j, p in jt.select(["v_allele", "j_allele", "p"]).iter_rows():
            out[j] = out.get(j, 0.0) + pv.get(v, 0.0) * p
        return out
    if seg == "d":
        pj = gene_marginal(model, "j")
        out = {}
        for j, d, p in t["d_gene"].select(["j_allele", "d_allele", "p"]).iter_rows():
            out[d] = out.get(d, 0.0) + pj.get(j, 0.0) * p
        return out
    if seg == "d2":  # tandem second D: Σ_d1 P(D1) P(D2|D1)
        pd1 = gene_marginal(model, "d")
        out = {}
        for d1, d2, p in t["d2_gene"].select(["d_allele", "d2_allele", "p"]).iter_rows():
            out[d2] = out.get(d2, 0.0) + pd1.get(d1, 0.0) * p
        return out
    raise ValueError(seg)


def _stationary(R: np.ndarray) -> np.ndarray:
    """Stationary distribution of a column-stochastic 4×4 dinucleotide matrix ``R[next, prev]``."""
    w, v = np.linalg.eig(R)
    k = int(np.argmin(np.abs(w - 1.0)))
    pi = np.real(v[:, k])
    pi = np.abs(pi)
    return pi / pi.sum() if pi.sum() else np.full(4, 0.25)


def _dinucl_R(df: pl.DataFrame) -> np.ndarray:
    R = np.zeros((4, 4))
    for frm, to, p in df.select(["from_nt", "to_nt", "p"]).iter_rows():
        R[to, frm] = p
    return R


def _cond_entropy(df: pl.DataFrame, parent_col: str, pmarg: dict[str, float]) -> tuple[float, float]:
    """Return ``(H_marginal, H_conditional)`` in bits for a table ``P(X | parent)``.

    ``X`` is every non-parent, non-``p`` column jointly. ``pmarg`` is the parent's marginal.
    """
    xcols = [c for c in df.columns if c not in (parent_col, "p")]
    hcond = 0.0
    for (pv,), grp in df.group_by([parent_col], maintain_order=True):
        w = pmarg.get(pv, 0.0)
        if w > 0:
            hcond += w * _H(grp["p"].to_numpy())
    marg = (
        df.with_columns(pl.col("p") * pl.col(parent_col).replace_strict(pmarg, default=0.0))
        .group_by(xcols)
        .agg(pl.col("p").sum())
    )
    return _H(marg["p"].to_numpy()), hcond


def entropy_table(model: Model) -> pl.DataFrame:
    """Per-event entropy (bits): marginal ``H(X)`` and conditional ``H(X | parents)``.

    Returns a tidy frame ``(event, kind, given, n_states, H_bits, H_cond_bits)`` — one row per
    event of the model's declared graph, in graph order.
    """
    t = model.tables
    rows: list[dict] = []

    def add(event: str, h: float, hc: float, n: int) -> None:
        ev = model.manifest.events[event]
        rows.append({
            "event": event, "kind": ev.kind.value, "given": ",".join(ev.given) or "-",
            "n_states": n, "H_bits": round(h, 4), "H_cond_bits": round(hc, 4),
        })

    for name, ev in model.manifest.events.items():
        df = t[name]
        given = ev.given
        if len(given) > 1:
            # single-parent factorizations only (V/D/J-conditioned). A multi-parent event
            # (e.g. P(D|V,J)) needs a joint parent marginal — not yet handled here.
            raise NotImplementedError(f"analyze: event {name!r} has >1 parent; multi-parent MI unsupported")
        if ev.kind.value == "gene_choice":
            seg = _allele_col(ev).split("_")[0]
            marg = gene_marginal(model, seg)
            h = _H(np.fromiter(marg.values(), float))
            if given:
                _, hc = _cond_entropy(df, _allele_col(model.manifest.events[given[0]]), gene_marginal(model, given[0].split("_")[0]))
            else:
                hc = h
            add(name, h, hc, len(marg))
        elif ev.kind.value in ("deletion", "deletion_2d"):
            pcol = _allele_col(model.manifest.events[given[0]])
            h, hc = _cond_entropy(df, pcol, gene_marginal(model, given[0].split("_")[0]))
            realiz = ["ndel"] if ev.kind.value == "deletion" else ["ndel5", "ndel3"]
            add(name, h, hc, df.select(realiz).n_unique())  # distinct deletion states, not rows
        elif ev.kind.value == "ins_length":
            h = _H(df["p"].to_numpy())
            add(name, h, h, df.height)
        elif ev.kind.value == "dinucleotide":
            # H(X) here is the entropy of the Markov *stationary* base composition, and H(X|prev)
            # the per-step conditional entropy — a composition summary of the N-region, independent
            # of the insertion-length distribution (which the ins_length event carries separately).
            R = _dinucl_R(df)
            pi = _stationary(R)
            hc = float(sum(pi[frm] * _H(R[:, frm]) for frm in range(4)))
            add(name, _H(pi), hc, 4)
        elif ev.kind.value == "n_d":
            add(name, _H(df["p"].to_numpy()), _H(df["p"].to_numpy()), df.height)
    return pl.DataFrame(rows)


def mutual_information(model: Model) -> pl.DataFrame:
    """Mutual information (bits) carried by informative pairs of the model.

    One row per declared parent→child edge (``I(child; parent) = H(child) − H(child|parent)``),
    plus ``I(V; J)`` (0 by construction for a VDJ model — V, J are independent roots) and the
    within-D deletion coupling ``I(delD5; delD3 | D)`` — for the second D too on a tandem model.
    """
    t = model.tables
    rows: list[dict] = []
    for name, ev in model.manifest.events.items():
        if len(ev.given) > 1:
            raise NotImplementedError(f"analyze: event {name!r} has >1 parent; multi-parent MI unsupported")
        if not ev.given:
            continue
        parent = ev.given[0]
        pcol = _allele_col(model.manifest.events[parent])
        hm, hc = _cond_entropy(t[name], pcol, gene_marginal(model, parent.split("_")[0]))
        rows.append({"a": name, "b": parent, "mi_bits": round(hm - hc, 4)})

    # I(V; J): explicit — 0 for VDJ (independent roots), >0 for VJ (encoded in P(J|V)).
    if model.chain_type == "VDJ":
        rows.append({"a": "v_choice", "b": "j_choice", "mi_bits": 0.0})

    # Within-D 5'/3' deletion coupling: the genuine conditional MI, E_D[ I(delD5; delD3 | D) ]
    # (averaging over D, not marginalizing — a D-marginal joint would inflate MI via Simpson mixing).
    for tbl, acol, seg, lab in (("d_del", "d_allele", "d", "delD5"), ("d2_del", "d2_allele", "d2", "delD2_5")):
        if tbl in t:
            mi = _within_d_deletion_mi(t[tbl], acol, gene_marginal(model, seg))
            rows.append({"a": lab, "b": lab.replace("5", "3"), "mi_bits": round(mi, 4)})
    return pl.DataFrame(rows)


def _within_d_deletion_mi(df: pl.DataFrame, allele_col: str, dmarg: dict[str, float]) -> float:
    """E_D[ I(delD5; delD3 | D) ] (bits) — 5'/3' trim coupling within a D, averaged over D usage."""
    total = 0.0
    for (d,), grp in df.group_by([allele_col], maintain_order=True):
        w = dmarg.get(d, 0.0)
        if w <= 0:
            continue
        m5: dict[int, float] = {}
        m3: dict[int, float] = {}
        cells: dict[tuple[int, int], float] = {}
        for n5, n3, p in grp.select(["ndel5", "ndel3", "p"]).iter_rows():
            cells[(n5, n3)] = cells.get((n5, n3), 0.0) + p
            m5[n5] = m5.get(n5, 0.0) + p
            m3[n3] = m3.get(n3, 0.0) + p
        mi_d = _H(np.fromiter(m5.values(), float)) + _H(np.fromiter(m3.values(), float)) - _H(np.fromiter(cells.values(), float))
        total += w * mi_d
    return total


def total_entropy(model: Model) -> pl.DataFrame:
    """Per-event contribution to the **scenario entropy** ``H(recombination event)``, in bits.

    For a Bayes net the joint entropy is the sum of each node's conditional entropy given its
    parents, so this frame's ``contribution_bits`` sums to the entropy of one whole recombination
    scenario — the information content of the process, and the basis of the ``2^H`` diversity
    figure in :func:`vdjtools.model.score.diversity`.

    The insertion regions need care, and this is the only place the accounting is non-obvious.
    :func:`entropy_table` reports a dinucleotide event's **per-step** conditional entropy (a
    composition summary, independent of how long the N-region is). Its contribution to the scenario
    is that per-step entropy times the expected number of steps, so the dinucleotide row here is
    ``E[length] · H_step`` using its paired ``*_ins`` event's mean length. The first inserted base
    is treated as drawn from the chain's stationary distribution rather than the model's separate
    first-base bias — the standard decomposition, and worth a fraction of a bit at most.

    Args:
        model: The model to measure.

    Returns:
        ``event, kind, contribution_bits`` — one row per event, in graph order.

    Example:
        >>> total_entropy(m)["contribution_bits"].sum()   # bits per rearrangement
    """
    ent = entropy_table(model)
    cond = {r["event"]: r["H_cond_bits"] for r in ent.to_dicts()}
    rows = []
    for name, ev in model.manifest.events.items():
        kind = ev.kind.value
        bits = cond.get(name, 0.0)
        if kind == "dinucleotide":
            ins = name.replace("_dinucl", "_ins")
            t = model.tables.get(ins)
            mean_len = 0.0
            if t is not None:
                mean_len = float((t["length"].cast(pl.Float64) * t["p"]).sum())
            bits = bits * mean_len
        rows.append({"event": name, "kind": kind, "contribution_bits": round(bits, 6)})
    return pl.DataFrame(rows, schema={"event": pl.Utf8, "kind": pl.Utf8,
                                      "contribution_bits": pl.Float64})


# --- graphviz (bnlearn-style) --------------------------------------------------------------

_KIND_COLOR = {
    "gene_choice": "#cfe8ff", "n_d": "#ffd6a5", "deletion": "#d7f0d7",
    "deletion_2d": "#d7f0d7", "ins_length": "#f0e0ff", "dinucleotide": "#f0e0ff",
}


def bayes_net_dot(model: Model, *, title: str | None = None) -> str:
    """Graphviz DOT for the model's Bayes net: nodes labelled with ``H(X)``, edges with ``I``."""
    ent = {r["event"]: r["H_bits"] for r in entropy_table(model).to_dicts()}
    mi = {(r["a"], r["b"]): r["mi_bits"] for r in mutual_information(model).to_dicts()}
    lab = title or f"{model.organism} {model.locus} ({model.chain_type})  ·  {model.manifest.source}"
    out = ["digraph bn {", '  rankdir=LR;', '  node [style=filled, fontname="Helvetica", shape=ellipse];',
           '  edge [fontname="Helvetica", fontsize=9];', f'  labelloc="t"; label="{lab}";']
    for name, ev in model.manifest.events.items():
        h = ent.get(name, 0.0)
        color = _KIND_COLOR.get(ev.kind.value, "#eeeeee")
        out.append(f'  "{name}" [fillcolor="{color}", label="{name}\\nH={h:.2f} bits"];')
    for name, ev in model.manifest.events.items():
        for parent in ev.given:
            w = mi.get((name, parent))
            edge_lab = f' [label="I={w:.2f}"]' if w is not None else ""
            out.append(f'  "{parent}" -> "{name}"{edge_lab};')
    out.append("}")
    return "\n".join(out)


def render_dot(dot: str, path: str | Path, *, fmt: str = "pdf") -> Path:
    """Render any DOT source to ``path`` via the graphviz ``dot`` CLI; returns the output path.

    Args:
        dot: DOT source, e.g. from :func:`bayes_net_dot` or :func:`compare_net_dot`.
        path: Output path; the suffix is replaced with ``fmt``.
        fmt: Any format ``dot`` supports (``pdf``, ``png``, ``svg``, ...).

    Raises:
        RuntimeError: If the ``dot`` CLI is not on ``PATH``.
    """
    if not shutil.which("dot"):
        raise RuntimeError("graphviz 'dot' CLI not found on PATH (brew install graphviz)")
    out = Path(path).with_suffix(f".{fmt}")
    subprocess.run(["dot", f"-T{fmt}", "-o", str(out)], input=dot, text=True, check=True)
    return out


def render_bayes_net(model: Model, path: str | Path, *, fmt: str = "pdf") -> Path:
    """Render :func:`bayes_net_dot` to ``path`` via the ``dot`` CLI; returns the output path."""
    return render_dot(bayes_net_dot(model), path, fmt=fmt)


def compare_entropy(models: dict[str, Model]) -> pl.DataFrame:
    """Stack :func:`entropy_table` across models into a wide ``event × model`` H(X) matrix."""
    wide: dict[str, dict[str, float]] = {}
    for label, m in models.items():
        for r in entropy_table(m).to_dicts():
            wide.setdefault(r["event"], {})[label] = r["H_bits"]
    rows = [{"event": ev, **cols} for ev, cols in wide.items()]
    return pl.DataFrame(rows)


# --- model vs model -------------------------------------------------------------------------

def _to_gene(df: pl.DataFrame) -> pl.DataFrame:
    """Collapse every ``*_allele`` column to gene level and re-sum ``p`` within the new key."""
    allele_cols = [c for c in df.columns if c.endswith("_allele")]
    if not allele_cols:
        return df
    out = df.with_columns([pl.col(c).str.split("*").list.first().alias(c) for c in allele_cols])
    keys = [c for c in out.columns if c != "p"]
    return out.group_by(keys).agg(pl.col("p").sum())


def _tv_jsd(pa: np.ndarray, pb: np.ndarray) -> tuple[float, float]:
    """(total variation, Jensen-Shannon divergence in bits) between two aligned distributions."""
    sa, sb = pa.sum(), pb.sum()
    if sa <= 0 or sb <= 0:
        return float("nan"), float("nan")
    pa, pb = pa / sa, pb / sb
    tv = float(0.5 * np.abs(pa - pb).sum())
    m = 0.5 * (pa + pb)
    ok = m > 0
    with np.errstate(divide="ignore", invalid="ignore"):
        ka = np.where(pa[ok] > 0, pa[ok] * np.log(pa[ok] / m[ok]), 0.0).sum()
        kb = np.where(pb[ok] > 0, pb[ok] * np.log(pb[ok] / m[ok]), 0.0).sum()
    return tv, float(0.5 * (ka + kb) / _LOG2)


def compare_models(a: Model, b: Model, *, labels: tuple[str, str] = ("a", "b"),
                   by: str = "allele") -> pl.DataFrame:
    """Per-event distance between two models — the parameter-level ``compare_networks``.

    The two models' tables are aligned on the **union** of their realization keys with zero fill,
    so a gene one model knows and the other does not contributes to the distance instead of being
    dropped. For a conditioned event the distance is computed per parent group and averaged
    **weighted by the parent's marginal** (the same weighting the conditional-entropy code uses),
    so a rarely-used V's deletion profile cannot dominate the number; ``tv_max`` reports the worst
    single group, which is what finds the one broken gene an average hides.

    Args:
        a: First model.
        b: Second model.
        labels: Names for the two models (used in error messages and the DOT title).
        by: ``"allele"`` (default) or ``"gene"``. Use ``"gene"`` to compare models built on
            different germline vintages or sources — an OLGA-namespace model against an
            arda-namespace one only lines up at gene level.

    Returns:
        One row per event of the union of both graphs: ``event, kind, given, status, n_groups,
        support_a, support_b, support_shared, support_only_a, support_only_b, tv, tv_max,
        jsd_bits``. ``status`` is ``shared``, ``only_a``, ``only_b``, or ``schema_differs`` (the
        event exists in both but is factorized differently, e.g. ``j_choice`` is ``P(J|V)`` on a VJ
        locus and a root ``P(J)`` on a VDJ one). Distances are null unless the status is ``shared``.

    Note:
        Jensen-Shannon is the primary metric: it is symmetric, bounded by 1 bit, and **finite when
        the supports differ**, which is exactly the case here. KL is deliberately not reported — it
        is infinite whenever one model assigns zero to something the other does not.

    Example:
        >>> compare_models(load_bundled("TRB", "olga"), load_bundled("TRB", "learned"), by="gene")
    """
    if by not in ("allele", "gene"):
        raise ValueError(f"by must be 'allele' or 'gene', got {by!r}")
    rows = []
    for name in dict.fromkeys([*a.manifest.events, *b.manifest.events]):
        ev = a.manifest.events.get(name) or b.manifest.events[name]
        status = "shared" if name in a.manifest.events and name in b.manifest.events else (
            "only_a" if name in a.manifest.events else "only_b")
        base = {"event": name, "kind": ev.kind.value, "given": ",".join(ev.given) or "-",
                "status": status}
        if status != "shared":
            rows.append({**base, "n_groups": 0, "support_a": 0, "support_b": 0,
                         "support_shared": 0, "support_only_a": 0, "support_only_b": 0,
                         "tv": None, "tv_max": None, "jsd_bits": None})
            continue

        ta, tb = a.tables[name], b.tables[name]
        if by == "gene":
            ta, tb = _to_gene(ta), _to_gene(tb)
        if set(ta.columns) != set(tb.columns):
            # Same event name, different factorization -- e.g. j_choice is P(J|V) on a VJ locus
            # but a root P(J) on a VDJ one. The two parameterize different things, so there is no
            # meaningful distance; say so instead of joining on a column one side lacks.
            rows.append({**base, "status": "schema_differs", "n_groups": 0,
                         "support_a": int((ta["p"] > 0).sum()),
                         "support_b": int((tb["p"] > 0).sum()),
                         "support_shared": 0, "support_only_a": 0, "support_only_b": 0,
                         "tv": None, "tv_max": None, "jsd_bits": None})
            continue
        keys = [c for c in ta.columns if c != "p"]
        joined = ta.join(tb, on=keys, how="full", coalesce=True, suffix="_b").with_columns(
            pl.col("p").fill_null(0.0), pl.col("p_b").fill_null(0.0))

        group_keys = [c for c in normalization_keys(ev) if c in joined.columns]
        weights = _parent_weights(a, b, ev, group_keys, by)
        tvs, jsds, ws = [], [], []
        groups = ([(tuple(k), g) for k, g in joined.group_by(group_keys, maintain_order=True)]
                  if group_keys else [((), joined)])
        for key, grp in groups:
            tv, jsd = _tv_jsd(grp["p"].to_numpy(), grp["p_b"].to_numpy())
            if np.isnan(tv):
                continue
            w = weights.get(key[0] if len(key) == 1 else key, 1.0) if weights else 1.0
            tvs.append(tv); jsds.append(jsd); ws.append(w)

        wsum = sum(ws)
        tv_w = sum(t * w for t, w in zip(tvs, ws)) / wsum if wsum > 0 else None
        jsd_w = sum(t * w for t, w in zip(jsds, ws)) / wsum if wsum > 0 else None
        sa = int((joined["p"] > 0).sum())
        sb = int((joined["p_b"] > 0).sum())
        shared = int(((joined["p"] > 0) & (joined["p_b"] > 0)).sum())
        rows.append({**base, "n_groups": len(tvs), "support_a": sa, "support_b": sb,
                     "support_shared": shared, "support_only_a": sa - shared,
                     "support_only_b": sb - shared,
                     "tv": None if tv_w is None else round(tv_w, 6),
                     "tv_max": round(max(tvs), 6) if tvs else None,
                     "jsd_bits": None if jsd_w is None else round(jsd_w, 6)})
    return pl.DataFrame(rows)


def _parent_weights(a: Model, b: Model, ev, group_keys: list[str], by: str) -> dict:
    """Marginal weight per conditioning group, averaged over the two models."""
    if len(group_keys) != 1 or not ev.given:
        return {}
    seg = group_keys[0].split("_")[0]
    if seg not in ("v", "j", "d", "d2"):
        return {}
    try:
        ma, mb = gene_marginal(a, seg), gene_marginal(b, seg)
    except (KeyError, ValueError):
        return {}
    out: dict[str, float] = {}
    for key in set(ma) | set(mb):
        out[key.split("*")[0] if by == "gene" else key] = (
            out.get(key.split("*")[0] if by == "gene" else key, 0.0)
            + 0.5 * (ma.get(key, 0.0) + mb.get(key, 0.0)))
    return out


def compare_usage(a: Model, b: Model, seg: str = "v", *, by: str = "gene") -> pl.DataFrame:
    """Side-by-side gene (or allele) usage of two models — the protocol-bias view.

    Args:
        a: First model.
        b: Second model.
        seg: ``"v"``, ``"j"``, ``"d"`` or ``"d2"``.
        by: ``"gene"`` (default) or ``"allele"``. Gene level is the meaningful comparison: allele
            calls on short reads are mismapping-prone, so allele-resolution usage is noise.

    Returns:
        ``name, p_a, p_b, log2_ratio`` sorted by descending ``p_a``, over the union of both models'
        genes. ``log2_ratio`` is null where either side is zero.
    """
    ma, mb = gene_marginal(a, seg), gene_marginal(b, seg)
    if by == "gene":
        ma, mb = _collapse_marginal(ma), _collapse_marginal(mb)
    rows = [{"name": k, "p_a": ma.get(k, 0.0), "p_b": mb.get(k, 0.0)}
            for k in sorted(set(ma) | set(mb))]
    return (pl.DataFrame(rows, schema={"name": pl.Utf8, "p_a": pl.Float64, "p_b": pl.Float64})
            .with_columns(log2_ratio=pl.when((pl.col("p_a") > 0) & (pl.col("p_b") > 0))
                          .then((pl.col("p_b") / pl.col("p_a")).log(2))
                          .otherwise(None))
            .sort("p_a", descending=True))


def _collapse_marginal(marg: dict[str, float]) -> dict[str, float]:
    out: dict[str, float] = {}
    for allele, p in marg.items():
        gene = allele.split("*")[0]
        out[gene] = out.get(gene, 0.0) + p
    return out


def compare_net_dot(a: Model, b: Model, *, labels: tuple[str, str] = ("a", "b"),
                    title: str | None = None) -> str:
    """Graphviz DOT contrasting two models' Bayes nets — the ``compare_networks`` picture.

    One DAG over the union of both graphs. Edges present in both are solid black; an edge only in
    ``a`` is blue and dashed, only in ``b`` red and dotted. Node fill intensity scales with that
    event's Jensen-Shannon divergence between the two models, and node labels carry
    ``ΔH = H_a − H_b`` so the structural and the quantitative differences are visible at once.

    Args:
        a: First model.
        b: Second model.
        labels: Names for the two models, shown in the title and legend.
        title: Graph title; defaults to a summary of the two models.

    Returns:
        DOT source — render it with :func:`render_dot`.
    """
    la, lb = labels
    diff = {r["event"]: r for r in compare_models(a, b).to_dicts()}
    ha = {r["event"]: r["H_bits"] for r in entropy_table(a).to_dicts()}
    hb = {r["event"]: r["H_bits"] for r in entropy_table(b).to_dicts()}
    lab = title or (f"{la} vs {lb}  ·  {a.organism} {a.locus} ({a.chain_type}) vs "
                    f"{b.organism} {b.locus} ({b.chain_type})")

    out = ["digraph compare {", '  rankdir=LR;',
           '  node [style=filled, fontname="Helvetica", shape=ellipse];',
           '  edge [fontname="Helvetica", fontsize=9];',
           f'  labelloc="t"; label="{lab}";']
    for name, r in diff.items():
        jsd = r["jsd_bits"]
        if r["status"] != "shared":
            fill, extra = ("#cfe0ff" if r["status"] == "only_a" else "#ffd6d6"), r["status"]
            node_lab = f"{name}\\n({extra.replace('only_a', la + ' only').replace('only_b', lb + ' only')})"
        else:
            # White at JSD 0 -> saturated orange at 1 bit, the JSD maximum, so the scale is absolute.
            t = min(max(jsd or 0.0, 0.0), 1.0)
            fill = "#{:02x}{:02x}{:02x}".format(255, int(255 - 90 * t), int(255 - 200 * t))
            node_lab = (f"{name}\\nJSD={jsd:.3f} bits\\n"
                        f"dH={ha.get(name, 0.0) - hb.get(name, 0.0):+.2f}")
        out.append(f'  "{name}" [fillcolor="{fill}", label="{node_lab}"];')

    edges_a = {(p, c) for c, ev in a.manifest.events.items() for p in ev.given}
    edges_b = {(p, c) for c, ev in b.manifest.events.items() for p in ev.given}
    for parent, child in sorted(edges_a | edges_b):
        if (parent, child) in edges_a and (parent, child) in edges_b:
            style = 'color="black"'
        elif (parent, child) in edges_a:
            style = f'color="blue", style=dashed, label="{la} only"'
        else:
            style = f'color="red", style=dotted, label="{lb} only"'
        out.append(f'  "{parent}" -> "{child}" [{style}];')
    out.append("}")
    return "\n".join(out)
