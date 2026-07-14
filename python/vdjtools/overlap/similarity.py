"""Sequence-similarity-weighted repertoire overlap (TINA / Leinster-Cobbold).

Where :mod:`vdjtools.overlap.metrics` scores overlap on an *exact* clonotype match and
:mod:`vdjtools.overlap.fuzzy` on a *within-edit-scope* match, this module scores it
through a **continuous CDR3 similarity kernel** ``Z`` (Leinster-Cobbold / Nei; Schmidt
et al. 2016 *TINA*). Two repertoires are relative-abundance vectors ``p, q`` over their
clonotypes and the overlap is built on the bilinear form ``pᵀZq = Σᵢⱼ pᵢ Zᵢⱼ qⱼ``:

- **cosine** (TINA_w): ``S = pᵀZq / sqrt((pᵀZp)(qᵀZq))``
- **Morisita-Horn**: ``S = 2·pᵀZq / (pᵀZp + qᵀZq)``
- ``distance = 1 − S``.

``Z`` is a **kernel of a seqtree gap-block alignment penalty** ``P`` (``≥0``, ``0``
identical), symmetric, with ``Zᵢᵢ = 1``:

- ``kernel="exp"`` — ``Zᵢⱼ = exp(−Pᵢⱼ/τ)`` (Leinster-Cobbold / Nei); ``τ`` defaults to
  ``SubstitutionMatrix.blosum62().scale()`` (= 14).
- ``kernel="step"`` — ``Zᵢⱼ = 1[Pᵢⱼ ≤ max_penalty]``. On unit cost with indels prohibited
  this is exactly vdjmatch's fuzzy edit-distance overlap.
- ``kernel="identity"`` — ``Z = I`` (match on the clonotype key). This recovers the exact
  frequency overlap: cosine collapses to the classical cosine of the shared-frequency
  vectors and Morisita to the classical Morisita-Horn index.

The three kernels are the same code on the same spine — identity and step are the exact
special cases (exact / fuzzy overlap) of the continuous ``exp`` form. The CDR3 kernel is
**block-diagonalised by the non-cdr3 key fields** (V/J/locus): exp/step never connect
clonotypes that differ outside the CDR3, so identity is the exact ``τ→0`` limit of exp on
the same key.

The penalty is built with **seqtree** (a base dependency): the **dense** path scores
every clonotype pair via :func:`seqtree.score_matrix` (``O(N²)``, for small ``N`` and the
tests); the **sparse** path uses :func:`seqtree.pairwise_batch` only to *find* near
candidates, then **re-scores each with the same gap-block model**
(:func:`seqtree.gapblock_score`, same matrix / gap_open / gap_prior as dense) and keeps
those with penalty ``≤`` the threshold — so dense and sparse agree on every retained pair.
It assembles a ``scipy.sparse`` ``Z`` (the practical path at scale). Both within-sample
blocks (``Z_AA``, ``Z_BB``) are needed in full — a diagonal-only approximation of ``pᵀZp``
is wrong.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import polars as pl

from ..io.schema import JUNCTION_AA, COUNT

_SEQTREE_HINT = (
    "seqtree is required for vdjtools.overlap.similarity; install the extra with "
    "seqtree is a base dependency of vdjtools -- reinstall with `pip install --force-reinstall vdjtools`."
)

#: Standard 20 amino acids; CDR3s with any other symbol (``*``, ``X``, digits) are dropped
#: before indexing (seqtree's ``"aa"`` alphabet rejects them).
_STANDARD_AA = r"^[ACDEFGHIKLMNPQRSTVWY]+$"

#: Above this clonotype count the dense ``O(N²)`` path auto-switches to the sparse path.
_DENSE_MAX_N = 1500


def _require_seqtree():
    """Import and return the ``seqtree`` module; raise a helpful error if missing."""
    try:
        import seqtree  # noqa: F401
    except ImportError as exc:  # pragma: no cover - exercised only without seqtree
        raise ImportError(_SEQTREE_HINT) from exc
    return seqtree


@dataclass
class SimilarityMatrices:
    """The three similarity blocks plus the aligned weight vectors for a sample pair.

    Attributes:
        z_ab: Cross-sample kernel, shape ``(n_a, n_b)`` (dense ``np.ndarray`` or a
            ``scipy.sparse`` matrix).
        z_aa: Within-``a`` kernel, shape ``(n_a, n_a)``; diagonal 1, symmetric.
        z_bb: Within-``b`` kernel, shape ``(n_b, n_b)``.
        keys_a: Clonotype key tuples for ``a`` in row order.
        keys_b: Clonotype key tuples for ``b`` in row order.
        freq_a: Within-sample relative abundance of each ``a`` clonotype (sums to 1).
        freq_b: Within-sample relative abundance of each ``b`` clonotype.
        kernel: The kernel used (``"exp"``, ``"step"``, ``"identity"``).
        tau: The kernel bandwidth ``τ`` (``None`` for non-``exp`` kernels).
        sparse: Whether the blocks are ``scipy.sparse`` matrices.
    """

    z_ab: Any
    z_aa: Any
    z_bb: Any
    keys_a: list
    keys_b: list
    freq_a: np.ndarray
    freq_b: np.ndarray
    kernel: str
    tau: float | None
    sparse: bool


def _aggregate(df: pl.DataFrame, key: list[str]) -> pl.DataFrame:
    """Collapse to unique clonotype keys (summing counts), drop non-standard-AA CDR3s,
    recompute within-sample frequency, and attach a 0-based row index."""
    if JUNCTION_AA not in key:
        raise ValueError(f"key must contain {JUNCTION_AA!r} (the similarity match unit); "
                         f"got {tuple(key)!r}")
    agg = df.group_by(key, maintain_order=True).agg(pl.col(COUNT).sum().alias("_count"))
    agg = agg.filter(pl.col(JUNCTION_AA).str.contains(_STANDARD_AA))
    total = agg["_count"].sum() or 1
    return agg.with_columns(
        (pl.col("_count") / total).alias("_freq"),
        pl.int_range(pl.len(), dtype=pl.Int64).alias("_idx"),
    )


def _keys(agg: pl.DataFrame, key: list[str]) -> list[tuple]:
    """Row-ordered clonotype key tuples (for identity matching)."""
    return list(agg.select(key).iter_rows())


def _block_labels(keys_a: list, keys_b: list, key: list[str]):
    """Integer block ids from the NON-cdr3 key fields (V/J/locus/…), over a shared space.

    Two clonotypes share a block iff their non-cdr3 key columns are all equal. The exp/step
    CDR3 kernels are zeroed across blocks (block-diagonalisation) so identity — which matches
    the full key — is the exact ``τ→0`` limit of exp on the same spine. Returns
    ``(labels_a, labels_b)`` as ``int`` arrays, or ``(None, None)`` when the key is cdr3-only
    (no gating).
    """
    if len(key) <= 1:
        return None, None
    cut = key.index(JUNCTION_AA)
    idx: dict = {}

    def label(k):
        return idx.setdefault(k[:cut] + k[cut + 1:], len(idx))

    la = np.array([label(k) for k in keys_a], dtype=np.int64)
    lb = np.array([label(k) for k in keys_b], dtype=np.int64)
    return la, lb


def _resolve_matrix(seqtree, matrix, kernel: str):
    """Resolve the substitution matrix: caller value, else blosum62 (exp) / unit (step)."""
    if matrix is not None:
        return matrix
    if kernel == "exp":
        return seqtree.SubstitutionMatrix.blosum62()
    return None  # step: unit cost (edit-count penalty), matches vdjmatch


def _gap_prior(seqtree, gap_prior, matrix):
    """Map a ``gap_prior`` spec to a seqtree ``GapPrior`` (or ``None``).

    Placement of the single indel block only matters for unequal-length CDR3 pairs.
    ``"frame"`` pins the block to a common frame column (the only transitive rule);
    ``"central"`` biases it to the loop centre; ``"none"`` lets the score decide.
    """
    if gap_prior is None or callable(gap_prior):
        return gap_prior
    scale = matrix.scale() if matrix is not None else 1
    lam = round(1.5 * scale)
    if gap_prior == "frame":
        return seqtree.frame_prior(lam, 0)  # left-anchored common frame
    if gap_prior == "central":
        return seqtree.central_prior(lam)
    if gap_prior == "none":
        return None
    raise ValueError(f"gap_prior must be 'frame', 'central', 'none', a GapPrior, or None; "
                     f"got {gap_prior!r}")


def _dense_penalty(seqtree, qs, rs, matrix, gap_open, gap_prior, threads) -> np.ndarray:
    """Full dense penalty matrix of ``qs`` (rows) against ``rs`` (cols)."""
    sm = seqtree.score_matrix(qs, rs, matrix=matrix, gap_open=gap_open,
                              gap_prior=gap_prior, threads=threads)
    return np.asarray(sm).astype(float)


def _kernel_dense(P: np.ndarray, kernel: str, tau: float, max_penalty) -> np.ndarray:
    """Apply the kernel to a dense penalty matrix."""
    if kernel == "exp":
        return np.exp(-P / tau)
    if kernel == "step":
        return (P <= max_penalty).astype(float)
    raise ValueError(f"kernel must be 'exp', 'step' or 'identity'; got {kernel!r}")


def _identity_block(keys_x: list, keys_y: list, sparse: bool):
    """Indicator ``Z[i,j] = 1[keys_x[i] == keys_y[j]]`` (dense or sparse)."""
    index = {}
    for j, k in enumerate(keys_y):
        index.setdefault(k, []).append(j)
    rows, cols = [], []
    for i, k in enumerate(keys_x):
        for j in index.get(k, ()):
            rows.append(i)
            cols.append(j)
    if sparse:
        from scipy.sparse import coo_matrix
        data = np.ones(len(rows))
        return coo_matrix((data, (rows, cols)), shape=(len(keys_x), len(keys_y))).tocsr()
    Z = np.zeros((len(keys_x), len(keys_y)))
    Z[rows, cols] = 1.0
    return Z


def _sparse_candidate_pairs(seqtree, qs, rs, matrix, gap_open, max_penalty, threads):
    """Near-neighbour ``(i, j)`` candidate pairs of ``qs`` vs ``rs`` (for later re-scoring).

    Uses :func:`seqtree.pairwise_batch` purely to *find* candidates. Its (prior-free)
    penalty is a lower bound on the gap-block penalty **with** a (non-negative) prior, so
    this candidate set is a superset of the pairs the re-scoring keeps — nothing scoring
    ``≤ max_penalty`` under the dense model is missed.
    """
    params = seqtree.SearchParams(
        max_subs=max_penalty, max_ins=max_penalty, max_dels=max_penalty,
        max_total_edits=2 * max_penalty, max_penalty=max_penalty,
        matrix=matrix if matrix is not None else "",
        gap_open=gap_open, gap_extend=1, engine="auto", mode="all")
    res = seqtree.pairwise_batch(qs, rs, params, "aa", threads)
    for i, hits in enumerate(res):
        for h in hits:
            yield i, h.ref_id


def _sparse_kernel_block(seqtree, qs, rs, matrix, gap_open, gap_prior, kernel, tau,
                         max_penalty, threads, block_q, block_r, self_block: bool):
    """Sparse kernel block: a genuine sparsification of the dense gap-block model.

    :func:`pairwise_batch` finds candidate near pairs; each is **re-scored with the same
    gap-block model** (:func:`seqtree.gapblock_score`, same ``matrix`` / ``gap_open`` /
    ``gap_prior`` as dense) and kept iff its penalty ``≤ max_penalty``. V/J gating zeroes
    any pair whose non-cdr3 blocks differ. The diagonal is forced to 1 on a self block.
    """
    from scipy.sparse import coo_matrix
    cells: dict = {}
    for i, j in _sparse_candidate_pairs(seqtree, qs, rs, matrix, gap_open,
                                        max_penalty, threads):
        if block_q is not None and block_q[i] != block_r[j]:
            continue  # V/J/locus gate: differ outside the CDR3 → no similarity
        pen, _ = seqtree.gapblock_score(qs[i], rs[j], matrix=matrix,
                                        gap_open=gap_open, gap_prior=gap_prior)
        if pen > max_penalty:
            continue
        cells[(i, j)] = float(np.exp(-pen / tau)) if kernel == "exp" else 1.0
    if cells:
        rows, cols = zip(*cells.keys())
        data = list(cells.values())
    else:
        rows, cols, data = (), (), ()
    Z = coo_matrix((data, (rows, cols)), shape=(len(qs), len(rs))).tocsr()
    if self_block:
        Z.setdiag(1.0)  # exact self-similarity even if the search missed the self hit
    return Z


def similarity_matrix(a: pl.DataFrame, b: pl.DataFrame, *,
                      key: "tuple[str, ...]" = (JUNCTION_AA,), kernel: str = "exp",
                      tau: float | None = None, matrix=None, max_penalty: int | None = None,
                      gap_prior="central",
                      gap_open: int | None = None, dense: bool | None = None,
                      threads: int = 0) -> SimilarityMatrices:
    """Build the three CDR3-similarity blocks ``Z_AB, Z_AA, Z_BB`` for a sample pair.

    Args:
        a: First clonotype frame (canonical schema).
        b: Second clonotype frame.
        key: Columns forming the clonotype identity; **must include** ``junction_aa`` (the
            similarity match unit). Default ``("junction_aa",)``.
        kernel: ``"exp"`` (``exp(−P/τ)``), ``"step"`` (``1[P ≤ max_penalty]``), or
            ``"identity"`` (``Z = I`` on the key — exact overlap, no seqtree).
        tau: Kernel bandwidth for ``"exp"``; defaults to the matrix scale
            (``blosum62().scale() == 14``).
        matrix: A ``seqtree.SubstitutionMatrix``; defaults to BLOSUM62 for ``"exp"`` and
            to unit cost (edit-count penalty) for ``"step"``.
        max_penalty: Penalty ceiling. Required (and the step threshold) for ``"step"``
            (default ``1`` = one edit); for the sparse ``"exp"`` path it is the
            neighbourhood cutoff (default from ``τ`` so ``exp(−P/τ) ≳ 1e-3``).
        gap_prior: Single-indel block placement rule: ``"central"`` (default — biases the
            indel to the loop centre, where CDR3 length variation sits; seqtree's own
            gapblock recommendation for pairwise scoring), ``"frame"`` (left-anchored
            common frame), ``"none"``, or a ``seqtree.GapPrior``.
        gap_open: Block-opening cost. Defaults to the seqtree score default for ``"exp"``
            and to a gap-prohibiting value for ``"step"`` (substitution-only, matching
            vdjmatch's default scope).
        dense: Force the dense (``True``) or sparse (``False``) path; ``None`` auto-selects
            (dense while ``max(n_a, n_b) ≤ 1500``).
        threads: Worker threads for the native search (``0`` = all cores).

    Returns:
        A :class:`SimilarityMatrices`.

    Raises:
        ImportError: If seqtree (or, for the sparse path, scipy) is missing.
        ValueError: On an unknown ``kernel``/``gap_prior`` or a ``key`` without ``junction_aa``.
    """
    key = list(key)
    a_agg = _aggregate(a, key)
    b_agg = _aggregate(b, key)
    keys_a, keys_b = _keys(a_agg, key), _keys(b_agg, key)
    block_a, block_b = _block_labels(keys_a, keys_b, key)
    freq_a = a_agg["_freq"].to_numpy().astype(float)
    freq_b = b_agg["_freq"].to_numpy().astype(float)
    n_a, n_b = len(keys_a), len(keys_b)

    if kernel == "identity":
        # Exact-overlap spine: no seqtree, Z = I on the clonotype key.
        use_sparse = bool(dense is False)
        z_ab = _identity_block(keys_a, keys_b, use_sparse)
        z_aa = _identity_block(keys_a, keys_a, use_sparse)
        z_bb = _identity_block(keys_b, keys_b, use_sparse)
        return SimilarityMatrices(z_ab, z_aa, z_bb, keys_a, keys_b,
                                  freq_a, freq_b, kernel, None, use_sparse)

    seqtree = _require_seqtree()
    sub_matrix = _resolve_matrix(seqtree, matrix, kernel)
    if tau is None:
        tau = float(sub_matrix.scale()) if sub_matrix is not None else 14.0
    if kernel == "step":
        if max_penalty is None:
            max_penalty = 1
        if gap_open is None:
            gap_open = 10 ** 6  # prohibit indels: substitution-only, = vdjmatch scope
    else:  # exp
        if max_penalty is None:
            # Sparse exp neighbourhood cutoff: drop kernel weights below ~1e-3.
            max_penalty = int(np.ceil(tau * np.log(1e3)))
        if gap_open is None:
            # seqtree's own default (2·scale); make it explicit so SearchParams can use it.
            gap_open = 2 * (sub_matrix.scale() if sub_matrix is not None else 1)
    prior = _gap_prior(seqtree, gap_prior, sub_matrix)

    use_dense = dense if dense is not None else max(n_a, n_b) <= _DENSE_MAX_N
    cdr3_a = a_agg[JUNCTION_AA].to_list()
    cdr3_b = b_agg[JUNCTION_AA].to_list()

    if use_dense:
        p_ab = _dense_penalty(seqtree, cdr3_a, cdr3_b, sub_matrix, gap_open, prior, threads)
        p_aa = _dense_penalty(seqtree, cdr3_a, cdr3_a, sub_matrix, gap_open, prior, threads)
        p_bb = _dense_penalty(seqtree, cdr3_b, cdr3_b, sub_matrix, gap_open, prior, threads)
        p_aa = np.minimum(p_aa, p_aa.T)  # defensive symmetrisation of the square blocks
        p_bb = np.minimum(p_bb, p_bb.T)
        z_ab = _kernel_dense(p_ab, kernel, tau, max_penalty)
        z_aa = _kernel_dense(p_aa, kernel, tau, max_penalty)
        z_bb = _kernel_dense(p_bb, kernel, tau, max_penalty)
        if block_a is not None:
            # Block-diagonalise by the non-cdr3 key fields (V/J/locus): zero every entry
            # whose clonotypes differ outside the CDR3.
            z_ab = z_ab * (block_a[:, None] == block_b[None, :])
            z_aa = z_aa * (block_a[:, None] == block_a[None, :])
            z_bb = z_bb * (block_b[:, None] == block_b[None, :])
        return SimilarityMatrices(z_ab, z_aa, z_bb, keys_a, keys_b,
                                  freq_a, freq_b, kernel, tau, False)

    z_ab = _sparse_kernel_block(seqtree, cdr3_a, cdr3_b, sub_matrix, gap_open, prior,
                                kernel, tau, max_penalty, threads, block_a, block_b,
                                self_block=False)
    z_aa = _sparse_kernel_block(seqtree, cdr3_a, cdr3_a, sub_matrix, gap_open, prior,
                                kernel, tau, max_penalty, threads, block_a, block_a,
                                self_block=True)
    z_bb = _sparse_kernel_block(seqtree, cdr3_b, cdr3_b, sub_matrix, gap_open, prior,
                                kernel, tau, max_penalty, threads, block_b, block_b,
                                self_block=True)
    return SimilarityMatrices(z_ab, z_aa, z_bb, keys_a, keys_b,
                              freq_a, freq_b, kernel, tau, True)


def _weights(freq: np.ndarray, weight: str) -> np.ndarray:
    """Per-clonotype weight vector for the requested mode."""
    n = len(freq)
    if weight in ("freq", "frequency"):
        return freq
    if weight == "presence":
        return np.full(n, 1.0 / n) if n else np.zeros(0)
    raise ValueError(f"weight must be 'freq' or 'presence'; got {weight!r}")


def _quad(x: np.ndarray, Z, y: np.ndarray) -> float:
    """Bilinear form ``xᵀ Z y`` for a dense or sparse ``Z``."""
    if x.size == 0 or y.size == 0:
        return 0.0
    return float(x @ (Z @ y))


def similarity_overlap(a: pl.DataFrame, b: pl.DataFrame, *,
                       key: "tuple[str, ...]" = (JUNCTION_AA,), metric: str = "cosine",
                       weight: str = "freq", kernel: str = "exp",
                       tau: float | None = None, matrix=None,
                       max_penalty: int | None = None, dense: bool | None = None,
                       threads: int = 0) -> dict:
    """Similarity-weighted overlap between two repertoires.

    Builds ``Z`` via :func:`similarity_matrix`, forms the weight vectors ``p, q``, and
    returns the cosine (TINA_w) or Morisita-Horn similarity on ``pᵀZq``.

    Args:
        a: First clonotype frame.
        b: Second clonotype frame.
        key: Clonotype identity key (must include ``junction_aa``).
        metric: ``"cosine"`` (``pᵀZq / sqrt(pᵀZp·qᵀZq)``) or ``"morisita"``
            (``2·pᵀZq / (pᵀZp + qᵀZq)``).
        weight: ``"freq"`` (relative abundance) or ``"presence"`` (uniform per clonotype,
            TINA-unweighted).
        kernel: ``"exp"``, ``"step"``, or ``"identity"`` (see :func:`similarity_matrix`).
            ``"identity"`` recovers the classical (exact) cosine / Morisita-Horn.
        tau: Kernel bandwidth for ``"exp"`` (default ``14``).
        matrix: Substitution matrix (see :func:`similarity_matrix`).
        max_penalty: Penalty ceiling / step threshold (see :func:`similarity_matrix`).
        dense: Force the dense (``True``) or sparse (``False``) kernel path; ``None``
            auto-selects (see :func:`similarity_matrix`). Both agree on retained pairs.
        threads: Worker threads for the native search.

    Returns:
        Dict with keys ``similarity, distance, pTZq, pTZp, qTZq, metric, kernel``.

    Raises:
        ValueError: On an unknown ``metric``/``weight``/``kernel``.
    """
    sm = similarity_matrix(a, b, key=key, kernel=kernel, tau=tau, matrix=matrix,
                           max_penalty=max_penalty, dense=dense, threads=threads)
    p = _weights(sm.freq_a, weight)
    q = _weights(sm.freq_b, weight)
    pTZq = _quad(p, sm.z_ab, q)
    pTZp = _quad(p, sm.z_aa, p)
    qTZq = _quad(q, sm.z_bb, q)

    if metric == "cosine":
        denom = np.sqrt(pTZp * qTZq)
        sim = pTZq / denom if denom > 0 else 0.0
    elif metric == "morisita":
        denom = pTZp + qTZq
        sim = 2.0 * pTZq / denom if denom > 0 else 0.0
    else:
        raise ValueError(f"metric must be 'cosine' or 'morisita'; got {metric!r}")

    return {"similarity": float(sim), "distance": float(1.0 - sim),
            "pTZq": pTZq, "pTZp": pTZp, "qTZq": qTZq,
            "metric": metric, "kernel": kernel}
