"""Clustering-evaluation metrics (purity / homogeneity / parsimony / q-measure).

Given a ground-truth labelling (e.g. the antigen/epitope a clonotype binds) and a
predicted clustering (cluster ids), these functions score how well the clustering
recovers the truth. They are the ``clustereval`` family of information-theoretic and
set-overlap metrics, reimplemented here on a contingency matrix so they can grade any
clonotype clustering — TCRnet components, metaclonotype groups, GLIPH motifs, …

The contingency matrix ``n[i, j]`` counts items with true class ``i`` and predicted
cluster ``j`` (rows = true classes, columns = clusters), built with
:func:`sklearn.metrics.cluster.contingency_matrix`. All entropies use natural logs.

**Singleton convention.** An unclustered item must get its **own** unique cluster id
rather than being dropped or lumped into a shared "noise" cluster — otherwise purity
and homogeneity are silently inflated. :func:`assign_singleton_ids` maps a sentinel
(``None`` / ``-1``) onto distinct negative ids for exactly this.
"""
from __future__ import annotations

import numpy as np

_SKLEARN_HINT = (
    "scikit-learn is required for vdjtools.sc.cluster_eval; install the extra with "
    "`pip install 'vdjtools[sc]'` (or `pip install scikit-learn`)."
)


def _contingency(labels_true, labels_pred) -> np.ndarray:
    """Return the dense ``n[i, j]`` contingency matrix (rows=true, cols=pred)."""
    try:
        from sklearn.metrics.cluster import contingency_matrix
    except ImportError as exc:  # pragma: no cover - exercised only without sklearn
        raise ImportError(_SKLEARN_HINT) from exc
    if len(labels_true) != len(labels_pred):
        raise ValueError(
            f"labels_true and labels_pred must be equal length; got "
            f"{len(labels_true)} and {len(labels_pred)}"
        )
    if len(labels_true) == 0:
        raise ValueError("labels must be non-empty")
    return np.asarray(contingency_matrix(labels_true, labels_pred), dtype=np.float64)


def _entropy(counts: np.ndarray) -> float:
    """Shannon entropy (natural log) of a count vector; 0 for an empty/degenerate one."""
    counts = counts[counts > 0]
    total = counts.sum()
    if total <= 0:
        return 0.0
    p = counts / total
    return float(-(p * np.log(p)).sum())


def assign_singleton_ids(pred, *, sentinel=None):
    """Give every unclustered item (``sentinel``) its own unique negative cluster id.

    Args:
        pred: Iterable of predicted cluster ids; ``sentinel`` marks unclustered items.
        sentinel: The value flagging "no cluster" (default ``None``; ``-1`` is common).

    Returns:
        A ``list`` of cluster ids with each sentinel replaced by a distinct negative
        integer (``-1, -2, …``), so no two unclustered items share a cluster.

    Example:
        >>> assign_singleton_ids([5, None, 5, None])
        [5, -1, 5, -2]
    """
    out = []
    nxt = -1
    for value in pred:
        if value == sentinel or (sentinel is None and value is None):
            out.append(nxt)
            nxt -= 1
        else:
            out.append(value)
    return out


# --------------------------------------------------------------------------- #
# metric cores — operate on a precomputed contingency matrix ``n``
# --------------------------------------------------------------------------- #
def _purity(n: np.ndarray) -> float:
    return float(n.max(axis=0).sum() / n.sum())


def _normalized_purity(n: np.ndarray) -> float:
    total = n.sum()
    pur = n.max(axis=0).sum() / total
    pmin = n.sum(axis=1).max() / total
    if pmin >= 1.0:
        return 1.0
    return float((pur - pmin) / (1.0 - pmin))


def _inverse_purity(n: np.ndarray) -> float:
    return float(n.max(axis=1).sum() / n.sum())


def _normalized_inverse_purity(n: np.ndarray) -> float:
    total = n.sum()
    inv = n.max(axis=1).sum() / total
    imin = n.shape[0] / total
    if imin >= 1.0:
        return 1.0
    return float((inv - imin) / (1.0 - imin))


def _homogeneity(n: np.ndarray) -> float:
    h_c = _entropy(n.sum(axis=1))
    if h_c <= 0.0:
        return 1.0
    h_k = _entropy(n.sum(axis=0))
    h_ck = _entropy(n.ravel())
    h_c_given_k = max(0.0, h_ck - h_k)
    return float(1.0 - h_c_given_k / h_c)


def _parsimony(n: np.ndarray) -> float:
    h_c = _entropy(n.sum(axis=1))
    denom = np.log(n.sum()) - h_c
    if denom <= 0.0:
        return 1.0
    h_ck = _entropy(n.ravel())
    h_k_given_c = max(0.0, h_ck - h_c)
    return float(1.0 - h_k_given_c / denom)


def _q_measure(n: np.ndarray, beta: float = 1.0) -> float:
    h = _homogeneity(n)
    p = _parsimony(n)
    if h <= 0.0 or p <= 0.0:
        return 0.0
    return float((1.0 + beta) * h * p / (beta * h + p))


# --------------------------------------------------------------------------- #
# public metrics — each builds its own contingency from the label arrays
# --------------------------------------------------------------------------- #
def purity(labels_true, labels_pred) -> float:
    """Cluster purity: mean over clusters of the dominant true class fraction.

    ``purity = (1/N) · Σ_j max_i n[i, j]`` — 1.0 when every cluster is class-pure.
    """
    return _purity(_contingency(labels_true, labels_pred))


def normalized_purity(labels_true, labels_pred) -> float:
    """Purity rescaled against its one-cluster floor ``pmin = max_i n_i. / N``.

    ``(purity - pmin) / (1 - pmin)``; returns 1.0 when ``pmin == 1`` (a single true
    class). This maps the trivial "everything in one cluster" purity to 0 and the
    perfect clustering to 1.
    """
    return _normalized_purity(_contingency(labels_true, labels_pred))


def inverse_purity(labels_true, labels_pred) -> float:
    """Inverse purity: mean over true classes of the dominant cluster fraction.

    ``(1/N) · Σ_i max_j n[i, j]`` — the completeness-flavoured dual of purity.
    """
    return _inverse_purity(_contingency(labels_true, labels_pred))


def normalized_inverse_purity(labels_true, labels_pred) -> float:
    """Inverse purity rescaled against its all-singletons floor.

    ``(inv - imin) / (1 - imin)`` with ``imin = |unique(true)| / N`` (the inverse
    purity you get when every item is its own cluster); returns 1.0 when ``imin == 1``.
    """
    return _normalized_inverse_purity(_contingency(labels_true, labels_pred))


def homogeneity(labels_true, labels_pred) -> float:
    """Homogeneity: ``1 - H(C|K) / H(C)`` — do clusters contain a single true class?

    Returns 1.0 when ``H(C) == 0`` (a single true class, trivially homogeneous).
    ``H(C|K) = max(0, H(C,K) - H(K))`` with all entropies in natural logs.
    """
    return _homogeneity(_contingency(labels_true, labels_pred))


def parsimony(labels_true, labels_pred) -> float:
    """Parsimony: ``1 - H(K|C) / (ln N - H(C))`` — penalises fragmenting a class.

    Returns 1.0 when the denominator ``ln N - H(C)`` is 0. ``H(K|C) = max(0, H(C,K)
    - H(C))``. Falls to 0 when every item is its own cluster (maximal fragmentation).
    """
    return _parsimony(_contingency(labels_true, labels_pred))


def q_measure(labels_true, labels_pred, beta: float = 1.0) -> float:
    """Weighted harmonic mean of homogeneity ``h`` and parsimony ``p``.

    ``(1 + beta) · h · p / (beta · h + p)``; 0 when either ``h`` or ``p`` is ``≤ 0``
    (there is nothing to balance). ``beta`` weights homogeneity relative to parsimony.
    """
    return _q_measure(_contingency(labels_true, labels_pred), beta=beta)


def cluster_eval(labels_true, labels_pred, *, beta: float = 1.0) -> dict:
    """Compute the full clustering-evaluation metric suite in one pass.

    Args:
        labels_true: Ground-truth class label per item (e.g. bound epitope).
        labels_pred: Predicted cluster id per item (same length as ``labels_true``).
            Unclustered items should already carry unique ids — see
            :func:`assign_singleton_ids`.
        beta: Homogeneity weight for :func:`q_measure`.

    Returns:
        Dict with keys ``purity, normalized_purity, inverse_purity,
        normalized_inverse_purity, homogeneity, parsimony, q_measure``.

    Raises:
        ImportError: If scikit-learn is not installed (see the ``sc`` extra).
        ValueError: If the label arrays are empty or unequal length.
    """
    n = _contingency(labels_true, labels_pred)  # built once, shared by every metric
    return {
        "purity": _purity(n),
        "normalized_purity": _normalized_purity(n),
        "inverse_purity": _inverse_purity(n),
        "normalized_inverse_purity": _normalized_inverse_purity(n),
        "homogeneity": _homogeneity(n),
        "parsimony": _parsimony(n),
        "q_measure": _q_measure(n, beta=beta),
    }
