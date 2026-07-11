"""Pairwise sample distances and low-dimensional clustering.

Formalises the legacy ``CalcPairwiseDistances`` + ``ClusterSamples`` workflow (and
the ad-hoc MDS the aging example notebook did inline): compute an all-pairs distance
matrix from a repertoire-overlap metric, then embed it in 2-D (MDS) or build a
hierarchy (hclust).

An overlap *similarity* is turned into a *distance* with the legacy per-metric
normalisation (``OverlapMetricNormalization``):

- ``F``, ``F2``, ``D`` (frequency/diversity overlaps, ``(0, 1]``) → ``-log10(x + 1e-9)``;
- ``R`` (correlation, ``[-1, 1]``) → ``(1 - x) / 2``;
- ``jaccard`` (similarity index, ``[0, 1]``) → ``1 - x``.

The diagonal is forced to ``0`` and the matrix is symmetric. When a ``scope`` is
passed the (fuzzy) :func:`vdjtools.overlap.fuzzy.fuzzy_overlap_metrics` ``fuzzy_F`` is
used instead of the exact metric.
"""
from __future__ import annotations

import math

import numpy as np
import polars as pl

from ..io.schema import CDR3_AA, J_CALL, V_CALL
from .fuzzy import fuzzy_overlap_metrics
from .metrics import overlap_pair

#: Default clonotype match key (CDR3 aa + V + J), matching the exact-overlap default.
DEFAULT_KEY = (CDR3_AA, V_CALL, J_CALL)

_SKLEARN_HINT = (
    "scikit-learn is required for cluster_samples(method='mds'); install the extra "
    "with `pip install 'vdjtools[overlap]'` (or `pip install scikit-learn`)."
)

#: similarity-metric name -> distance transform (legacy OverlapMetricNormalization).
_TRANSFORM = {
    "F": lambda x: -math.log10(x + 1e-9),
    "F2": lambda x: -math.log10(x + 1e-9),
    "D": lambda x: -math.log10(x + 1e-9),
    "R": lambda x: (1.0 - x) / 2.0,
    "jaccard": lambda x: 1.0 - x,
}


def _named(samples) -> "list[tuple[str, pl.DataFrame]]":
    """Normalise ``list | dict`` of samples to an ordered ``(name, frame)`` list."""
    if isinstance(samples, dict):
        return list(samples.items())
    return [(str(i), df) for i, df in enumerate(samples)]


def _similarity(a: pl.DataFrame, b: pl.DataFrame, metric: str,
                key, scope) -> float:
    """Raw overlap *similarity* between two samples for the requested metric."""
    if scope is not None:
        # Fuzzy path: only the frequency-weighted fuzzy-F is defined here.
        if metric != "F":
            raise ValueError("fuzzy distances (scope=) support metric='F' only")
        return fuzzy_overlap_metrics(a, b, scope=scope)["fuzzy_F"]
    _, m = overlap_pair(a, b, key=key)
    if metric == "jaccard":
        denom = m["d1"] + m["d2"] - m["d12"]
        return m["d12"] / denom if denom else 0.0
    if metric == "R":
        return m["R"] if m["R"] is not None else 0.0  # legacy coerces undefined R -> 0
    if metric in ("F", "F2", "D"):
        return m[metric]
    raise ValueError(f"unknown metric {metric!r}; expected one of {sorted(_TRANSFORM)}")


def pairwise_distances(samples, metric: str = "F", key=DEFAULT_KEY,
                       scope: str | None = None, form: str = "matrix") -> pl.DataFrame:
    """All-pairs distance matrix over a collection of samples.

    Args:
        samples: A ``list`` of clonotype frames (named ``"0".."N-1"``) or a ``dict``
            mapping sample name to frame.
        metric: Overlap similarity to base the distance on: ``"F"``, ``"F2"``,
            ``"D"`` (→ ``-log10``), ``"R"`` (→ ``(1-x)/2``), or ``"jaccard"``
            (→ ``1-x``). See the module docstring.
        key: Exact-match clonotype key (default CDR3 aa + V + J); ignored when
            ``scope`` is given.
        scope: If set, use fuzzy overlap within this vdjmatch edit scope
            (``"subs,ins,dels,total"``) and the ``fuzzy_F`` similarity instead of the
            exact ``metric`` (only ``metric="F"`` is valid then).
        form: ``"matrix"`` for a wide frame (a ``sample`` column plus one column per
            sample) or ``"long"`` for a ``sample_a, sample_b, distance`` frame.

    Returns:
        A symmetric distance matrix with a zero diagonal, in the requested ``form``.
    """
    named = _named(samples)
    names = [n for n, _ in named]
    n = len(named)
    dist = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            sim = _similarity(named[i][1], named[j][1], metric, key, scope)
            d = _TRANSFORM[metric](sim)
            dist[i, j] = dist[j, i] = d

    if form == "long":
        rows = [(names[i], names[j], float(dist[i, j]))
                for i in range(n) for j in range(n)]
        return pl.DataFrame(rows, orient="row",
                            schema=["sample_a", "sample_b", "distance"])
    if form == "matrix":
        data = {"sample": names}
        for j, name in enumerate(names):
            data[name] = dist[:, j].tolist()
        return pl.DataFrame(data)
    raise ValueError(f"form must be 'matrix' or 'long'; got {form!r}")


def _matrix(dist: pl.DataFrame) -> "tuple[list[str], np.ndarray]":
    """Extract ``(sample_names, N×N array)`` from a matrix-form distance frame."""
    names = dist["sample"].to_list()
    mat = dist.select(names).to_numpy()
    return names, mat


def cluster_samples(dist: pl.DataFrame, method: str = "mds", n_components: int = 2,
                    metadata: pl.DataFrame | None = None) -> pl.DataFrame:
    """Embed / cluster samples from a precomputed distance matrix.

    Args:
        dist: A matrix-form distance frame from :func:`pairwise_distances` (a
            ``sample`` column plus one column per sample).
        method: ``"mds"`` — metric MDS (``sklearn.manifold.MDS`` with
            ``dissimilarity="precomputed"``) → ``n_components`` coordinate columns
            ``mds1..mdsK``; or ``"hclust"`` — average-linkage hierarchy
            (``scipy.cluster.hierarchy``) → a dendrogram ``leaf_order`` and a flat
            ``cluster`` label (``fcluster`` into ``n_components`` clusters).
        n_components: MDS output dimensionality (``method="mds"``) or the number of
            flat clusters (``method="hclust"``).
        metadata: Optional frame carrying a ``sample`` column plus per-sample columns
            (e.g. ``age``, ``group``) to left-join onto the result for colouring.

    Returns:
        A ``pl.DataFrame`` with one row per sample and the embedding / cluster
        columns, plus any joined ``metadata``.

    Raises:
        ImportError: If ``method="mds"`` and scikit-learn is not installed.
        ValueError: If ``method`` is not ``"mds"`` or ``"hclust"``.
    """
    names, mat = _matrix(dist)

    if method == "mds":
        try:
            from sklearn.manifold import MDS
        except ImportError as exc:  # pragma: no cover - exercised only without sklearn
            raise ImportError(_SKLEARN_HINT) from exc
        coords = MDS(n_components=n_components, dissimilarity="precomputed",
                     random_state=0, normalized_stress="auto").fit_transform(mat)
        out = pl.DataFrame({"sample": names} | {
            f"mds{k + 1}": coords[:, k].tolist() for k in range(n_components)
        })
    elif method == "hclust":
        from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
        from scipy.spatial.distance import squareform
        z = linkage(squareform(mat, checks=False), method="average")
        labels = fcluster(z, t=n_components, criterion="maxclust")
        order = dendrogram(z, no_plot=True)["leaves"]
        leaf_order = [order.index(i) for i in range(len(names))]
        out = pl.DataFrame({"sample": names, "cluster": labels.tolist(),
                            "leaf_order": leaf_order})
    else:
        raise ValueError(f"method must be 'mds' or 'hclust'; got {method!r}")

    if metadata is not None:
        out = out.join(metadata, on="sample", how="left")
    return out
