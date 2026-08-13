"""A frozen (V gene x k-mer) feature space: reduced alphabets, TF-IDF, and a sparse basis.

:mod:`vdjtools.features.kmer` produces tidy per-sample profiles. This module produces the thing a
*model* wants: a fixed-width vector per sample, on a vocabulary and a rotation that were decided
once and then never move, so two collaborators' matrices are comparable.

Three problems stand between a k-mer count and that vector, and each has a specific answer here.

**Sparsity.** A V-pinned 4-mer over 20 residues is a code space of ~8e6, of which one repertoire
touches ~4e5 and any two repertoires share far fewer. Two independent reductions apply: the
alphabet (:func:`reduced_alphabet` collapses the 20 residues into BLOSUM62 groups, so chemically
equivalent substitutions stop being different features) and the vocabulary (a document-frequency
window, below).

**Scale.** Raw k-mer weight tracks sequencing depth and CDR3 length before it tracks anything
immunological. TF-IDF plus an L2 row norm removes both: IDF down-weights the germline-adjacent
k-mers that appear in every repertoire, and the row norm makes a deep and a shallow sample
comparable.

**Noise.** The surviving matrix is still tens of thousands of columns wide and mostly zero, so it
is projected onto a small basis by truncated SVD. ``scipy.sparse.linalg.svds`` operates on the
sparse matrix directly and never densifies it -- and, deliberately, **scikit-learn is not used**:
scipy is a base dependency of vdjtools while sklearn is an extra, and a portable feature block
that imports an extra fails late, after every other block of a sample has already been computed.

The document-frequency window is what makes the vocabulary robust rather than merely small. A
k-mer present in nearly every sample carries no contrast (it is germline or near it), and one
present in a handful cannot be estimated; both ends are cut, which is the same reasoning that
bounds the public-clonotype panel's incidence window.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import polars as pl

from ..io.schema import JUNCTION_AA, V_CALL, strip_allele, weight_expr

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def _native():
    from .. import _core

    return _core


def reduced_alphabet(n_groups: int = 8) -> dict[str, int]:
    """Partition the 20 amino acids into ``n_groups`` by BLOSUM62 similarity.

    ``seqtree``'s ``SubstitutionMatrix.penalty(a, b) = s(a,a) + s(b,b) - 2 s(a,b)`` is the Gram
    transform of the log-odds matrix -- that is, a **squared Euclidean distance**. So the spectral
    embedding needs no Laplacian and no kernel choice: classical MDS (double-centre, eigendecompose)
    recovers coordinates whose pairwise distances are exactly that penalty, and the clustering
    happens there. Building an affinity matrix and a normalised Laplacian would be the same
    computation with an arbitrary bandwidth bolted on.

    Ward linkage on those coordinates, via ``scipy.cluster.hierarchy`` -- a base dependency, and
    for 20 points the choice of clusterer is not where the answer comes from.

    Args:
        n_groups: Number of groups. ``20`` returns the identity partition (each residue its own
            group), which is the un-reduced alphabet.

    Returns:
        ``{residue: group index}`` over the 20 standard amino acids, group indices contiguous
        from 0. Deterministic: no random initialisation anywhere in the path.
    """
    if not 2 <= n_groups <= 20:
        raise ValueError(f"n_groups must be in [2, 20]; got {n_groups}")
    if n_groups == 20:
        return {a: i for i, a in enumerate(AMINO_ACIDS)}

    import seqtree
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform

    m = seqtree.SubstitutionMatrix.blosum62()
    n = len(AMINO_ACIDS)
    d2 = np.array([[float(m.penalty(a, b)) for b in AMINO_ACIDS] for a in AMINO_ACIDS])

    # Classical MDS on the squared-distance matrix: B = -0.5 * J d2 J is the Gram matrix, and its
    # positive eigenpairs are the coordinates. Negative eigenvalues are clipped -- BLOSUM62's Gram
    # transform is only approximately Euclidean, and a negative eigenvalue is an imaginary axis,
    # not a real direction to cluster along.
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ d2 @ J
    w, V = np.linalg.eigh(B)
    keep = w > 1e-9
    coords = V[:, keep] * np.sqrt(w[keep])

    Z = linkage(squareform(_pairwise(coords), checks=False), method="ward")
    labels = fcluster(Z, t=n_groups, criterion="maxclust") - 1
    # Relabel in order of first appearance so the mapping does not depend on scipy's cluster
    # numbering, which is stable but undocumented.
    order, seen = {}, 0
    out = {}
    for a, lab in zip(AMINO_ACIDS, labels):
        if lab not in order:
            order[lab] = seen
            seen += 1
        out[a] = order[lab]
    return out


def _pairwise(X: np.ndarray) -> np.ndarray:
    """Euclidean distance matrix of the rows of ``X``."""
    d2 = np.maximum((X ** 2).sum(1)[:, None] + (X ** 2).sum(1)[None, :] - 2 * X @ X.T, 0.0)
    return np.sqrt(d2)


def alphabet_table(groups: dict[str, int]) -> tuple[np.ndarray, int]:
    """A 256-entry char -> group table for the native kernel, and the group count.

    Everything not in ``groups`` maps to ``-1``, which voids any window containing it. That
    includes ``X``, ``*`` and lowercase: a k-mer spanning an unknown residue is not a k-mer, and
    mapping it to a wildcard group would invent a count that was never observed.
    """
    table = np.full(256, -1, dtype=np.int8)
    for residue, g in groups.items():
        table[ord(residue)] = g
    return table, (max(groups.values()) + 1 if groups else 0)


def pattern_of(spec: str) -> list[int]:
    """Parse a k-mer shape: ``"xxx"`` is an ungapped 3-mer, ``"xx.x"`` a gapped 3-mer over 4.

    ``x`` keeps a position, ``.`` (or ``_``) skips it. A pattern must start and end with ``x`` --
    a leading or trailing gap is the same feature as the shorter pattern, shifted, and admitting
    both would put two names on one column.
    """
    keep = [1 if c == "x" else 0 for c in spec]
    if not keep or any(c not in "x._" for c in spec):
        raise ValueError(f"pattern spec must be 'x'/'.'/'_' characters; got {spec!r}")
    if not keep[0] or not keep[-1]:
        raise ValueError(f"pattern must start and end with 'x'; got {spec!r}")
    if sum(keep) < 1:
        raise ValueError(f"pattern must keep at least one position; got {spec!r}")
    return keep


@dataclass
class KmerSpace:
    """A frozen (V, k-mer) vocabulary with its IDF and optional SVD basis.

    Everything needed to turn a repertoire into the same vector next year, in another lab. Fitted
    once by :func:`fit_kmer_space`; applied by :meth:`transform`.
    """

    pattern: list[int]
    groups: dict[str, int]
    v_genes: list[str]
    flank: int
    codes: np.ndarray                 # selected code ids, ascending
    idf: np.ndarray                   # per selected column
    components: np.ndarray | None = None   # (n_components, n_columns) SVD basis, rows orthonormal
    meta: dict = field(default_factory=dict)
    _lookup: np.ndarray | None = field(default=None, repr=False, compare=False)

    @property
    def n_columns(self) -> int:
        return len(self.codes)

    @property
    def n_components(self) -> int:
        return 0 if self.components is None else self.components.shape[0]

    def lookup(self) -> np.ndarray:
        """code -> column index (``-1`` outside the vocabulary), built once and reused.

        Held as one int32 buffer and passed to the kernel by pointer. Materialised lazily because
        it is the size of the whole code space -- ~8e6 int32 = 33 MB for a V-pinned 4-mer over 20
        residues -- which is cheap once per corpus and ruinous once per sample.
        """
        if self._lookup is None:
            table, n_alpha = alphabet_table(self.groups)
            space = _native().kmer_code_space(self.pattern, n_alpha, len(self.v_genes))
            lk = np.full(space, -1, dtype=np.int32)
            lk[self.codes] = np.arange(len(self.codes), dtype=np.int32)
            self._lookup = lk
        return self._lookup

    def _encode(self, df: pl.DataFrame, weight: str):
        vidx = {v: i for i, v in enumerate(self.v_genes)}
        d = df.with_columns(strip_allele(pl.col(V_CALL).cast(pl.Utf8)).alias("_v"),
                            weight_expr(weight).alias("_w"))
        junctions = d[JUNCTION_AA].to_list()
        # An unmapped V gene becomes -1, which the kernel folds into its own bucket rather than
        # dropping. A collaborator on a different IMGT release or on Adaptive nomenclature has
        # genes we never saw; losing their k-mers silently is the failure mode this avoids.
        v_codes = np.fromiter((vidx.get(v, -1) for v in d["_v"]), dtype=np.int32, count=d.height)
        return junctions, v_codes, d["_w"].to_numpy().astype(float)

    def counts(self, df: pl.DataFrame, *, weight: str = "freq") -> np.ndarray:
        """Raw clone-weighted counts on this vocabulary, one entry per column."""
        junctions, v_codes, w = self._encode(df, weight)
        table, n_alpha = alphabet_table(self.groups)
        return np.asarray(_native().kmer_gather(
            junctions, v_codes, w, self.pattern, table.tolist(), self.lookup(),
            n_alpha, len(self.v_genes), self.flank, self.n_columns))

    def transform(self, df: pl.DataFrame, *, weight: str = "freq",
                  sublinear: bool = True) -> np.ndarray:
        """TF-IDF vector for one repertoire, L2-normalised, projected if a basis was fitted.

        Args:
            df: One repertoire's clonotype frame, already filtered. Pass the *work frame* if the
                caller uses one -- ``weight_expr`` reads whatever ``frequency`` holds.
            weight: Clone weight ladder, as elsewhere in vdjtools.
            sublinear: ``log1p`` the term frequency before IDF. On by default: a clone expanded
                1000-fold is not 1000 times more informative about which k-mers a repertoire
                contains, and without it one clonal expansion dominates the whole vector.

        Returns:
            ``(n_components,)`` if a basis was fitted, else ``(n_columns,)``.
        """
        tf = self.counts(df, weight=weight)
        if sublinear:
            tf = np.log1p(tf)
        x = tf * self.idf
        norm = np.linalg.norm(x)
        if norm > 0:
            x = x / norm
        return x if self.components is None else self.components @ x


def fit_kmer_space(frames, *, pattern: str = "xxxx", n_groups: int = 8, flank: int = 4,
                   v_genes=None, weight: str = "freq", min_df: float = 0.02,
                   max_df: float = 0.80, max_columns: int = 50_000,
                   n_components: int = 32, sublinear: bool = True,
                   threads: int = 0) -> KmerSpace:
    """Fit the vocabulary, IDF and SVD basis on a reference corpus.

    This is the only corpus-fitted part of the k-mer block, and it is fitted **once**, on studies
    disjoint from anything it will later be scored on. Nothing here is re-estimated per cohort;
    that is the whole point.

    Args:
        frames: Iterable of per-sample clonotype frames (one locus).
        pattern: Shape spec -- ``"xxx"``/``"xxxx"`` ungapped, ``"xx.x"``/``"x.x.x"`` gapped.
        n_groups: BLOSUM62 alphabet groups; ``20`` disables the reduction.
        flank: Residues trimmed from each end. ``4`` is the measured germline core (an N-terminal
            4-mer is shared by 31% of clonotypes, a central one by 0.08%).
        v_genes: V vocabulary; inferred from the corpus when ``None``.
        min_df: Drop k-mers seen in fewer than this fraction of samples -- too rare to estimate.
        max_df: Drop k-mers seen in more than this fraction -- present everywhere is no contrast,
            which for a CDR3 core means germline or near it.
        max_columns: Keep at most this many columns, highest document frequency first among the
            survivors. A cap, and reported when it binds rather than silently truncating.
        n_components: SVD rank; ``0`` ships the TF-IDF columns unprojected.
        sublinear: ``log1p`` the term frequency (see :meth:`KmerSpace.transform`).
        threads: Native worker threads across samples; ``0`` = auto.

    Returns:
        A fitted :class:`KmerSpace`.
    """
    frames = list(frames)
    if not frames:
        raise ValueError("need at least one frame to fit a k-mer space")
    pat = pattern_of(pattern)
    groups = reduced_alphabet(n_groups)
    table, n_alpha = alphabet_table(groups)

    if v_genes is None:
        seen = set()
        for df in frames:
            seen |= set(df.select(strip_allele(pl.col(V_CALL).cast(pl.Utf8)).alias("v"))["v"]
                        .drop_nulls().to_list())
        v_genes = sorted(seen)
    v_genes = list(v_genes)
    vidx = {v: i for i, v in enumerate(v_genes)}

    core = _native()
    space = core.kmer_code_space(pat, n_alpha, len(v_genes))

    junctions, v_codes, weights = [], [], []
    for df in frames:
        d = df.with_columns(strip_allele(pl.col(V_CALL).cast(pl.Utf8)).alias("_v"),
                            weight_expr(weight).alias("_w"))
        junctions.append(d[JUNCTION_AA].to_list())
        v_codes.append([vidx.get(v, -1) for v in d["_v"]])
        weights.append(d["_w"].to_numpy().astype(float).tolist())

    rows = core.kmer_rows(junctions, v_codes, weights, pat, table.tolist(),
                          n_alpha, len(v_genes), flank, threads)
    df_counts = np.asarray(core.kmer_document_frequency(rows, space), dtype=np.int64)

    n = len(frames)
    lo, hi = min_df * n, max_df * n
    keep = np.flatnonzero((df_counts >= lo) & (df_counts <= hi))
    if keep.size == 0:
        raise ValueError(f"no k-mer survives the document-frequency window "
                         f"[{min_df}, {max_df}] over {n} samples")
    dropped_cap = 0
    if keep.size > max_columns:
        dropped_cap = keep.size - max_columns
        keep = keep[np.argsort(df_counts[keep])[::-1][:max_columns]]
        keep.sort()

    # Smoothed IDF: +1 inside the log and on both counts, so a column at the window edge cannot
    # produce a zero or infinite weight.
    idf = np.log((1.0 + n) / (1.0 + df_counts[keep])) + 1.0
    sp = KmerSpace(pattern=pat, groups=groups, v_genes=v_genes, flank=flank,
                   codes=keep.astype(np.int64), idf=idf)

    if n_components > 0:
        X = _corpus_matrix(rows, sp, sublinear)
        k = min(n_components, min(X.shape) - 1)
        if k < 1:
            raise ValueError(f"cannot fit {n_components} components from a {X.shape} matrix")
        from scipy.sparse.linalg import svds

        # svds returns ascending singular values; flip so component 0 is the leading direction.
        # v0 is pinned: ARPACK starts from a random vector otherwise and the basis would differ
        # run to run, which for a *frozen* artifact is the one thing it must not do.
        m = min(X.shape)
        _, s, Vt = svds(X, k=k, v0=np.ones(m) / np.sqrt(m))
        order = np.argsort(s)[::-1]
        comp = Vt[order]
        # Sign is arbitrary in any SVD; pin it by the sign of each component's largest-magnitude
        # loading so a refit does not silently flip a column's meaning.
        lead = np.abs(comp).argmax(axis=1)
        signs = np.sign(comp[np.arange(comp.shape[0]), lead])
        signs[signs == 0] = 1.0
        sp.components = comp * signs[:, None]

    sp.meta = {"n_samples": n, "code_space": int(space), "surviving": int(keep.size),
               "dropped_by_cap": int(dropped_cap), "n_alphabet": int(n_alpha)}
    return sp


def _corpus_matrix(rows, space: KmerSpace, sublinear: bool):
    """The fitted corpus as a sparse TF-IDF matrix, on the selected vocabulary only.

    Built from the aggregated rows already in hand rather than re-reading the frames, and never
    densified: a 23k x 50k dense matrix is 9 GB, while the sparse form is the ~5% that is real.
    """
    from scipy.sparse import csr_matrix

    col_of = space.lookup()
    indptr, indices, data = [0], [], []
    for row in rows:
        codes = np.asarray(row.codes, dtype=np.int64)
        vals = np.asarray(row.weights, dtype=float)
        cols = col_of[codes]
        sel = cols >= 0
        indices.append(cols[sel])
        data.append(vals[sel])
        indptr.append(indptr[-1] + int(sel.sum()))
    idx = np.concatenate(indices) if indices else np.empty(0, dtype=np.int32)
    val = np.concatenate(data) if data else np.empty(0, dtype=float)
    if sublinear:
        val = np.log1p(val)
    val = val * space.idf[idx]
    X = csr_matrix((val, idx, np.asarray(indptr)), shape=(len(rows), space.n_columns))
    # L2-normalise rows: the same depth correction transform() applies, so the basis is fitted in
    # the geometry it will be used in.
    norms = np.sqrt(X.multiply(X).sum(axis=1)).A.ravel()
    norms[norms == 0] = 1.0
    from scipy.sparse import diags

    return diags(1.0 / norms) @ X


def save_kmer_spaces(spaces: dict[str, KmerSpace], path) -> None:
    """Freeze a per-locus set of spaces to one ``.npz``.

    A space that cannot be written down is not frozen, whatever the docstring says: the whole
    claim is that a collaborator recomputes the same columns, which requires the vocabulary, the
    IDF and the basis to travel. Stored per locus with the locus in the key.

    The alphabet is stored as the group id of each of the 20 residues in :data:`AMINO_ACIDS`
    order, not as the ``n_groups`` that produced it -- a future change to the clustering must not
    silently re-partition an already-frozen space.
    """
    import numpy as np

    out: dict[str, np.ndarray] = {"loci": np.array(sorted(spaces), dtype=object)}
    for locus, sp in sorted(spaces.items()):
        out[f"{locus}/pattern"] = np.asarray(sp.pattern, dtype=np.int8)
        out[f"{locus}/groups"] = np.asarray([sp.groups[a] for a in AMINO_ACIDS], dtype=np.int8)
        out[f"{locus}/v_genes"] = np.array(sp.v_genes, dtype=object)
        out[f"{locus}/flank"] = np.asarray(sp.flank, dtype=np.int32)
        out[f"{locus}/codes"] = sp.codes.astype(np.int64)
        out[f"{locus}/idf"] = sp.idf.astype(np.float64)
        if sp.components is not None:
            out[f"{locus}/components"] = sp.components.astype(np.float64)
    np.savez_compressed(path, **out)


def load_kmer_spaces(path) -> dict[str, KmerSpace]:
    """Read back what :func:`save_kmer_spaces` wrote."""
    import numpy as np

    z = np.load(path, allow_pickle=True)
    spaces: dict[str, KmerSpace] = {}
    for locus in [str(x) for x in z["loci"]]:
        groups = {a: int(g) for a, g in zip(AMINO_ACIDS, z[f"{locus}/groups"])}
        comp = z[f"{locus}/components"] if f"{locus}/components" in z.files else None
        spaces[locus] = KmerSpace(
            pattern=[int(x) for x in z[f"{locus}/pattern"]],
            groups=groups,
            v_genes=[str(v) for v in z[f"{locus}/v_genes"]],
            flank=int(z[f"{locus}/flank"]),
            codes=z[f"{locus}/codes"],
            idf=z[f"{locus}/idf"],
            components=comp,
        )
    return spaces
