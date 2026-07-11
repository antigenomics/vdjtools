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
from .schema import _allele_col

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


def render_bayes_net(model: Model, path: str | Path, *, fmt: str = "pdf") -> Path:
    """Render :func:`bayes_net_dot` to ``path`` via the ``dot`` CLI; returns the output path."""
    if not shutil.which("dot"):
        raise RuntimeError("graphviz 'dot' CLI not found on PATH (brew install graphviz)")
    out = Path(path).with_suffix(f".{fmt}")
    subprocess.run(["dot", f"-T{fmt}", "-o", str(out)], input=bayes_net_dot(model),
                   text=True, check=True)
    return out


def compare_entropy(models: dict[str, Model]) -> pl.DataFrame:
    """Stack :func:`entropy_table` across models into a wide ``event × model`` H(X) matrix."""
    wide: dict[str, dict[str, float]] = {}
    for label, m in models.items():
        for r in entropy_table(m).to_dicts():
            wide.setdefault(r["event"], {})[label] = r["H_bits"]
    rows = [{"event": ev, **cols} for ev, cols in wide.items()]
    return pl.DataFrame(rows)
