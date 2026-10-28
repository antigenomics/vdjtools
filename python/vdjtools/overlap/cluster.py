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

from ..io.schema import JUNCTION_AA, J_CALL, V_CALL
from .fuzzy import fuzzy_overlap_metrics
from .metrics import _aggregate, _overlap_from_agg
from .similarity import similarity_overlap

#: Default clonotype match key (CDR3 aa + V + J), matching the exact-overlap default.
DEFAULT_KEY = (JUNCTION_AA, V_CALL, J_CALL)

_SKLEARN_HINT = (
    "scikit-learn is required for cluster_samples(method='mds'); install the extra "
    "with `pip install 'vdjtools[overlap]'` (or `pip install scikit-learn`)."
)

#: similarity-metric name -> distance transform (legacy OverlapMetricNormalization).
#: The similarity-weighted overlaps (TINA cosine / Morisita) are already in ``[0, 1]``,
#: so their distance is simply ``1 - similarity``.
_TRANSFORM = {
    "F": lambda x: -math.log10(x + 1e-9),
    "F2": lambda x: -math.log10(x + 1e-9),
    "D": lambda x: -math.log10(x + 1e-9),
    "R": lambda x: (1.0 - x) / 2.0,
    "jaccard": lambda x: 1.0 - x,
    "similarity_cosine": lambda x: 1.0 - x,
    "similarity_morisita": lambda x: 1.0 - x,
}


def _named(samples) -> "list[tuple[str, pl.DataFrame]]":
    """Normalise ``list | dict`` of samples to an ordered ``(name, frame)`` list."""
    if isinstance(samples, dict):
        return list(samples.items())
    return [(str(i), df) for i, df in enumerate(samples)]


#: Metrics computed exactly from a shared-clonotype join (as opposed to the
#: engine-delegated fuzzy / sequence-similarity paths).
_EXACT = frozenset({"F", "F2", "D", "R", "jaccard"})


def _exact_similarity(m: dict, metric: str) -> float:
    """Overlap *similarity* for an exact metric, from an ``overlap_pair`` metrics dict."""
    if metric == "jaccard":
        denom = m["d1"] + m["d2"] - m["d12"]
        return m["d12"] / denom if denom else 0.0
    if metric == "R":
        return m["R"] if m["R"] is not None else 0.0  # legacy coerces undefined R -> 0
    return m[metric]  # F, F2, D


def _delegated_similarity(a: pl.DataFrame, b: pl.DataFrame, metric: str,
                          key, scope) -> float:
    """Similarity for the engine-delegated paths: fuzzy (``scope=``) and TINA cosine/morisita."""
    if scope is not None:
        # Fuzzy path: only the frequency-weighted fuzzy-F is defined here.
        if metric != "F":
            raise ValueError("fuzzy distances (scope=) support metric='F' only")
        return fuzzy_overlap_metrics(a, b, scope=scope)["fuzzy_F"]
    sub = "cosine" if metric == "similarity_cosine" else "morisita"
    return similarity_overlap(a, b, key=key, metric=sub)["similarity"]


def _exact_matrix(aggs, key, metric):
    """Batch exact overlaps through a shared sparse clonotype vocabulary."""
    from scipy.sparse import csr_matrix
    from ..cores import require_memory

    n = len(aggs)
    require_memory(48 * n * n + 64 * sum(a.height for a in aggs),
                   "exact-overlap numerical workspace (reduce samples or depth)")
    blocks = pl.concat([a.select(*key, "_freq", "_count").with_columns(pl.lit(i).alias("_sample"))
                        for i, a in enumerate(aggs)], how="vertical_relaxed").drop_nulls(key)
    values = blocks["_freq"].to_numpy()
    if not np.all(np.isfinite(values) & (values >= 0) & (values <= 1)):
        return None  # Preserve the pairwise arithmetic for signed/non-finite input frequencies.
    totals = None
    if metric in {"F", "F2"}:
        if any(a.schema["_count"] != pl.Int64 for a in aggs):
            return None
        totals = np.array([a["_count"].cast(pl.Float64).sum() or 1 for a in aggs])
        if np.any((totals <= 0) | (totals >= 2**53)):
            return None  # Keep legacy arithmetic beyond exact float64 integer accumulation.
    vocab = blocks.select(key).unique(maintain_order=True).with_row_index("_clone")
    encoded = blocks.join(vocab, on=key)
    rows, cols = encoded["_sample"].to_numpy(), encoded["_clone"].to_numpy()
    shape = (n, vocab.height)
    present = csr_matrix((np.ones(len(rows), dtype=np.int64), (rows, cols)), shape=shape)
    sizes = np.array([a.height for a in aggs], dtype=float)
    if metric in {"D", "jaccard"}:
        shared = (present @ present.T).toarray()
        denominator = (sizes[:, None] * sizes[None, :] if metric == "D" else
                       sizes[:, None] + sizes[None, :] - shared)
        similarity = np.divide(shared, denominator, out=np.zeros(shape=(n, n)),
                               where=denominator > 0)
    else:
        counts = csr_matrix((encoded["_count"].cast(pl.Float64).to_numpy(), (rows, cols)), shape=shape)
        if metric == "F":
            mass = (counts @ present.T).toarray() / totals[:, None]
            similarity = np.sqrt(mass * mass.T)
        else:
            counts.data = np.sqrt(counts.data)
            root_totals = np.sqrt(totals)
            similarity = (counts @ counts.T).toarray() / (root_totals[:, None] * root_totals[None, :])
    distance = 1 - similarity if metric == "jaccard" else -np.log10(similarity + 1e-9)
    np.fill_diagonal(distance, 0)
    return distance


def pairwise_distances(samples, metric: str = "F", key=DEFAULT_KEY,
                       scope: str | None = None, form: str = "matrix") -> pl.DataFrame:
    """All-pairs distance matrix over a collection of samples.

    Args:
        samples: A ``list`` of clonotype frames (named ``"0".."N-1"``) or a ``dict``
            mapping sample name to frame.
        metric: Overlap similarity to base the distance on: ``"F"``, ``"F2"``,
            ``"D"`` (→ ``-log10``), ``"R"`` (→ ``(1-x)/2``), ``"jaccard"`` (→ ``1-x``),
            or the sequence-similarity-weighted ``"similarity_cosine"`` /
            ``"similarity_morisita"`` (TINA / Leinster-Cobbold, → ``1-x``). See the
            module docstring.
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
    if metric not in _TRANSFORM:
        raise ValueError(f"unknown metric {metric!r}; expected one of {sorted(_TRANSFORM)}")
    if form not in {"matrix", "long"}:
        raise ValueError(f"form must be 'matrix' or 'long'; got {form!r}")
    named = _named(samples)
    names = [n for n, _ in named]
    n = len(named)
    from ..cores import require_memory
    workspace = (64 if form == "long" else 16) * n * n
    if form == "long":
        workspace += 2 * n * sum(len(str(name).encode()) for name in names)
    require_memory(workspace, "pairwise-distance output (reduce sample count)")
    # Aggregate once and batch the nonnegative frequency metrics. Pearson keeps its
    # centered pairwise calculation; raw-moment matrix products can lose its variance.
    if scope is None and metric in _EXACT:
        klist = list(key)
        aggs = [_aggregate(df, klist) for _, df in named]
        # With at most three samples there are at most three joins; avoid a global
        # vocabulary that regresses a few deep samples.
        batched = _exact_matrix(aggs, klist, metric) if n > 3 and metric != "R" else None
        if batched is not None:
            dist = batched
        else:
            dist = np.zeros((n, n), dtype=float)
            for i in range(n):
                for j in range(i + 1, n):
                    _, m = _overlap_from_agg(aggs[i], aggs[j], klist)
                    d = _TRANSFORM[metric](_exact_similarity(m, metric))
                    dist[i, j] = dist[j, i] = d
    else:
        dist = np.zeros((n, n), dtype=float)
        for i in range(n):
            for j in range(i + 1, n):
                sim = _delegated_similarity(named[i][1], named[j][1], metric, key, scope)
                d = _TRANSFORM[metric](sim)
                dist[i, j] = dist[j, i] = d

    if form == "long":
        if not n:
            return pl.DataFrame(schema=["sample_a", "sample_b", "distance"])
        return (pl.DataFrame({"sample_a": names})
                .join(pl.DataFrame({"sample_b": names}), how="cross", maintain_order="left_right")
                .with_columns(pl.Series("distance", dist.ravel())))
    return pl.DataFrame({"sample": names, **{name: dist[:, j] for j, name in enumerate(names)}})


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
