"""Sparse feature rows and bounded dense column workspaces for fitting only.

Zeros are implicit; unobserved row/block pairs have a separate mask. No feature
matrix or standardized sample-by-feature array is written to disk. Fitted parameters
use the existing Corpus/TablePCA formats.
"""
from dataclasses import dataclass
from ctypes import CDLL
import sys

import numpy as np
from scipy import sparse, linalg
from scipy.sparse.linalg import LinearOperator

# Byte bound for one dense column chunk; up to eight such temporaries coexist.
CHUNK_BYTES = 8 * 1024**2


def _width(n):
    return max(1, CHUNK_BYTES // (8 * max(n, 1)))


@dataclass
class Block:
    name: str
    values: object
    present: np.ndarray


class BlockMatrix:
    """Float64 sparse blocks with implicit zeros and explicit missing observations."""
    def __init__(self, blocks):
        self.blocks = blocks
        self.shape = (blocks[0].values.shape[0], sum(b.values.shape[1] for b in blocks))

    @property
    def nbytes(self):
        return sum(b.values.data.nbytes + b.values.indices.nbytes + b.values.indptr.nbytes
                   + b.present.nbytes for b in self.blocks)

    @classmethod
    def from_dense(cls, x, columns):
        from .layout import parse
        blocks, start = [], 0
        while start < len(columns):
            group = parse(columns[start])[1]
            end = start + 1
            while end < len(columns) and parse(columns[end])[1] == group:
                end += 1
            for j in range(start, end, _width(x.shape[0])):
                part = np.array(x[:, j:min(end, j + _width(x.shape[0]))], dtype=float, copy=True)
                present = np.isfinite(part).any(1)
                part[~present] = 0
                blocks.append(Block(group, sparse.csc_matrix(part), present))
            start = end
        return cls(blocks)

    def observed(self):
        return np.logical_or.reduce([b.present for b in self.blocks])

    def take_rows(self, keep):
        if keep.all():
            return self
        from ..cores import require_memory
        if (self.observed() & ~keep).any():
            require_memory(self.nbytes, 'sparse observed-row selection')
            return BlockMatrix([Block(b.name, b.values[keep].tocsc(), b.present[keep]) for b in self.blocks])
        # Absent rows contain no entries: only row indices change, not values/column pointers.
        require_memory(sum(b.values.indices.nbytes for b in self.blocks) + 16*len(keep),
                       'sparse observed-row indices')
        mapping = np.cumsum(keep, dtype=np.int64) - 1
        mapping[~keep] = -1
        blocks = []
        for b in self.blocks:
            indices = np.empty_like(b.values.indices)
            local = mapping.astype(indices.dtype, copy=False)
            # np.take casts source indices to intp even when an output is supplied.
            for start in range(0, len(indices), max(1, CHUNK_BYTES // 16)):
                stop = start + max(1, CHUNK_BYTES // 16)
                np.take(local, b.values.indices[start:stop], out=indices[start:stop], mode='clip')
            if indices.min(initial=0) < 0:
                raise ValueError('unobserved sparse rows contain stored entries')
            values = sparse.csc_matrix((b.values.data, indices, b.values.indptr),
                                       shape=(int(keep.sum()), b.values.shape[1]), copy=False)
            blocks.append(Block(b.name, values, b.present[keep]))
        return BlockMatrix(blocks)

    def chunks(self):
        offset = 0
        for b in self.blocks:
            for j in range(0, b.values.shape[1], _width(self.shape[0])):
                end = min(b.values.shape[1], j + _width(self.shape[0]))
                part = b.values[:, j:end].toarray()
                part[~b.present] = np.nan
                yield offset + j, offset + end, part
            offset += b.values.shape[1]


class _Row:
    """The existing observed-bin writer's tiny array sink, one group of one sample."""
    def __init__(self, width):
        self.width, self.present = width, False
        self.indices = np.empty(0, dtype=np.int32)
        self.values = np.empty(0, dtype=float)

    def fill(self, value):
        self.present = bool(np.isfinite(value))
        self.indices = np.empty(0, dtype=np.int32)
        self.values = np.empty(0, dtype=float)

    def __setitem__(self, key, values):
        values = np.asarray(values, dtype=float)
        if isinstance(key, slice):
            if key != slice(None):
                raise ValueError('sparse row accepts only a complete slice or explicit indices')
            assert values.size == self.width
            self.present = bool(np.isfinite(values).any())
            indices = np.flatnonzero(values != 0)
            values = values[indices]
        else:
            indices = np.asarray(key)
        self.indices = indices.astype(np.int32, copy=False)
        self.values = values


class MatrixBuilder:
    """Independent row slots; workers never mutate another sample's sparse storage."""
    def __init__(self, n, groups):
        self.n, self.groups = n, groups
        self.rows = {name: [_Row(width) for _ in range(n)] for name, width in groups}

    def row(self, i):
        return {name: rows[i] for name, rows in self.rows.items()}

    def finish(self):
        from ..cores import require_memory
        blocks = []
        for name, width in self.groups:
            rows = self.rows.pop(name)
            present = np.array([r.present for r in rows])
            lengths = np.array([r.values.size if r.present else 0 for r in rows], dtype=np.int64)
            count = int(lengths.sum())
            require_memory(24 * count + 16 * (self.n + width + 1), 'sparse feature assembly')
            data = np.empty(count)
            indices = np.empty(count, dtype=np.int32)
            indptr = np.r_[0, lengths.cumsum()]
            for i, row in enumerate(rows):
                if present[i]:
                    data[indptr[i]:indptr[i+1]] = row.values
                    indices[indptr[i]:indptr[i+1]] = row.indices
            del rows
            values = sparse.csr_matrix((data, indices, indptr), shape=(self.n, width)).tocsc()
            blocks.append(Block(name, values, present))
        return BlockMatrix(blocks)


class CenteredOperator(LinearOperator):
    """Sparse deviations + one observed-row outer product per block + column offset."""
    def __init__(self, blocks, offset):
        self.blocks, self.offset = blocks, offset
        super().__init__(dtype=np.dtype(float), shape=(blocks[0][0].shape[0], len(offset)))

    def _matmat(self, q):
        out = np.zeros((self.shape[0], q.shape[1]))
        start = 0
        for data, present, baseline in self.blocks:
            end = start + data.shape[1]
            out += data @ q[start:end]
            out += present[:, None] * (baseline @ q[start:end])[None, :]
            start = end
        out -= (self.offset @ q)[None, :]
        return out

    def _rmatmat(self, q):
        out = np.empty((self.shape[1], q.shape[1]))
        start, sums = 0, q.sum(0)
        for data, present, baseline in self.blocks:
            end = start + data.shape[1]
            out[start:end] = data.T @ q
            out[start:end] += baseline[:, None] * (present @ q)[None, :]
            out[start:end] -= self.offset[start:end, None] * sums[None, :]
            start = end
        return out

    def _matvec(self, q):
        return self._matmat(np.asarray(q).reshape(-1, 1))[:, 0]

    def _rmatvec(self, q):
        return self._rmatmat(np.asarray(q).reshape(-1, 1))[:, 0]


def randomized_svd(matrix, k, *, seed, n_iter=4, oversamples=10):
    """Seeded LU-normalized range algorithm on an operator.

    sklearn's public solver requires an array/sparse matrix, so its operator form
    uses the same NumPy RNG and public SciPy LU/QR/SVD kernels here. Oversampling,
    transpose choice and sign policy match its established defaults.
    """
    from sklearn.utils.extmath import svd_flip
    transpose = matrix.shape[1] > matrix.shape[0]
    a = matrix.T if transpose else matrix
    q = np.random.RandomState(seed).normal(size=(a.shape[1], k + oversamples))
    for _ in range(n_iter):
        q = linalg.lu(a @ q, permute_l=True)[0]
        q = linalg.lu(a.T @ q, permute_l=True)[0]
    q = linalg.qr(a @ q, mode='economic')[0]
    u, s, vt = linalg.svd((a.T @ q).T, full_matrices=False)
    u = q @ u
    u, vt = svd_flip(u, vt, u_based_decision=not transpose)
    return (vt[:k].T, s[:k], u[:, :k].T) if transpose else (u[:, :k], s[:k], vt[:k])


def _spans(values):
    """Bound sparse index/value temporaries by bytes, including empty columns."""
    start, width = 0, values.shape[1]
    while start < width:
        limit = values.indptr[start] + max(1, CHUNK_BYTES // 64)
        end = np.searchsorted(values.indptr, limit, side='right') - 1
        end = max(start + 1, min(width, end, start + max(1, CHUNK_BYTES // 256)))
        yield start, end
        start = end


def _indices(values, start, end):
    return np.repeat(np.arange(end-start), np.diff(values.indptr[start:end+1]))


def sparse_quantiles(block, probabilities):
    """Exact linear quantiles with implicit zeros and observed-only denominators."""
    values = block.values
    width, present = values.shape[1], int(block.present.sum())
    counts = np.zeros(width, dtype=np.int64)
    result = np.empty((len(probabilities), width))
    for start, end in _spans(values):
        ids = _indices(values, start, end)
        data = values.data[values.indptr[start]:values.indptr[end]]
        finite = np.isfinite(data)
        count = present - np.bincount(ids[~finite], minlength=end-start)
        counts[start:end] = count
        nonzero = finite & (data != 0)
        ids, data = ids[nonzero], data[nonzero]
        number = np.bincount(ids, minlength=end-start)
        negative = np.bincount(ids[data < 0], minlength=end-start)
        zeros = count - number
        offsets = np.r_[0, number.cumsum()[:-1]]
        ordered = np.r_[data[np.lexsort((data, ids))], 0.0]
        def kth(rank):
            in_zero = (rank >= negative) & (rank < negative + zeros)
            index = np.where(rank < negative, rank, rank - zeros)
            index = offsets + np.maximum(index, 0)
            index = np.minimum(index, ordered.size-1)
            return np.where(in_zero, 0.0, ordered[index])
        for i, probability in enumerate(probabilities):
            rank = probability * np.maximum(count-1, 0)
            fraction = rank - np.floor(rank)
            a, b = kth(np.floor(rank).astype(np.int64)), kth(np.ceil(rank).astype(np.int64))
            # NumPy's linear interpolation chooses the nearer endpoint for stability.
            interpolated = np.where(fraction < .5, a+(b-a)*fraction, b-(b-a)*(1-fraction))
            result[i, start:end] = np.where(count > 0, interpolated, np.nan)
    return counts, result


def _moments(values, present, baseline, *, fill=None, deviations=False):
    """Observed means/variances, or all-row moments after column-wise imputation."""
    width, n = values.shape[1], values.shape[0]
    counts, means, variances = np.zeros(width, dtype=np.int64), np.zeros(width), np.zeros(width)
    observed = int(present.sum())
    for start, end in _spans(values):
        ids = _indices(values, start, end)
        data = values.data[values.indptr[start]:values.indptr[end]]
        if deviations:
            data = data + baseline[start:end][ids]
        valid = np.isfinite(data)
        zero = observed - np.diff(values.indptr[start:end+1])
        count = zero + np.bincount(ids[valid], minlength=end-start)
        base = baseline[start:end]
        sums = zero*base + np.bincount(ids[valid], weights=data[valid], minlength=end-start)
        missing = n-count
        denominator = count if fill is None else np.full(end-start, n)
        if fill is not None:
            sums += missing*fill[start:end]
        mean = np.divide(sums, denominator, out=np.zeros(end-start), where=denominator > 0)
        squares = zero*(base-mean)**2
        squares += np.bincount(ids[valid], weights=(data[valid]-mean[ids[valid]])**2, minlength=end-start)
        if fill is not None:
            squares += missing*(fill[start:end]-mean)**2
        counts[start:end], means[start:end] = count, mean
        variances[start:end] = np.divide(squares, denominator, out=np.zeros(end-start), where=denominator > 0)
    return counts, means, variances


def _transformed(block, lo, hi, composition):
    values = block.values.copy()
    baseline = np.clip(np.zeros(values.shape[1]), lo, hi)
    if composition:
        baseline = np.sqrt(np.clip(baseline, 0, 1))
    for start, end in _spans(values):
        ids = _indices(values, start, end)
        begin, stop = values.indptr[start], values.indptr[end]
        data = np.clip(values.data[begin:stop], lo[start:end][ids], hi[start:end][ids])
        finite = data[np.isfinite(data)]
        if composition and ((finite < -1e-12) | (finite > 1+1e-12)).any():
            raise ValueError('composition frequencies must be in [0,1]')
        values.data[begin:stop] = np.sqrt(np.clip(data, 0, 1)) if composition else data
    return values, baseline


def fit_locus(x, columns, *, feature_definition, n_components, seed, max_components, winsor_p):
    from . import corpus as C, layout as L
    from ..cores import require_memory
    if not isinstance(x, BlockMatrix):
        x = BlockMatrix.from_dense(x, columns)
    n, p = x.shape
    cap = min(n_components if isinstance(n_components, int) else max_components, n-1, p)
    if cap < 1:
        raise ValueError('component cap must be positive')
    require_memory(3*x.nbytes + 64*(n+p)*(cap+10) + 8*CHUNK_BYTES + 160*p,
                   'bounded corpus PCA workspace')
    declarations = {g.name: g for g in L._feature_groups('vsig', feature_definition)}
    group_names = np.array([L.parse(c)[1] for c in columns])
    supports = [declarations[group].dynamic[1] if group in declarations and declarations[group].dynamic
                else L.support_of(c, feature_definition=feature_definition)
                for c, group in zip(columns, group_names)]
    bounds = {q: (np.full(p, -np.inf), np.full(p, np.inf)) for q in C.WINSOR_PS}
    loc, scale, variances, n_obs = np.zeros(p), np.zeros(p), np.zeros(p), np.zeros(p, dtype=np.int64)
    probabilities = [prob for q in C.WINSOR_PS for prob in (q, 1-q)]
    start = 0
    for block in x.blocks:
        end = start + block.values.shape[1]
        count, quantiles = sparse_quantiles(block, probabilities)
        n_obs[start:end] = count
        sides = np.array([L.SUPPORTS[s] for s in supports[start:end]])
        for i, q in enumerate(C.WINSOR_PS):
            lo, hi = quantiles[2*i:2*i+2]
            bounds[q][0][start:end] = np.where(sides[:, 0] & np.isfinite(lo), lo, -np.inf)
            bounds[q][1][start:end] = np.where(sides[:, 1] & np.isfinite(hi), hi, np.inf)
        work, baseline = _transformed(block, bounds[winsor_p][0][start:end],
                                      bounds[winsor_p][1][start:end], block.name in L.COMPOSITION_BLOCKS)
        if block.name in L.COMPOSITION_BLOCKS:
            _, loc[start:end], variances[start:end] = _moments(work, block.present, baseline)
        else:
            # The small diversity block keeps the original robust estimator exactly.
            for j in range(0, end-start, _width(n)):
                stop = min(end-start, j+_width(n))
                chunk = block.values[:, j:stop].toarray()
                np.clip(chunk, bounds[winsor_p][0][start+j:start+stop],
                        bounds[winsor_p][1][start+j:start+stop], out=chunk)
                chunk[~block.present] = np.nan
                loc[start+j:start+stop], scale[start+j:start+stop] = C.robust_loc_scale(chunk)
        start = end
        del work
    for group in L.COMPOSITION_BLOCKS:
        mask = group_names == group
        scale[mask] = np.sqrt(np.nansum(variances[mask]))
    if not (scale > C.MIN_SPREAD).any():
        return None
    blocks, offset, variance, start = [], np.zeros(p), np.zeros(p), 0
    for block in x.blocks:
        width = block.values.shape[1]
        end = start+width
        values, baseline = _transformed(block, bounds[winsor_p][0][start:end],
                                        bounds[winsor_p][1][start:end], block.name in L.COMPOSITION_BLOCKS)
        safe = np.where(scale[start:end] > C.MIN_SPREAD, scale[start:end], 1)
        base_z = C._standardize(baseline, loc[start:end], scale[start:end])
        for a, b in _spans(values):
            ids = _indices(values, a, b) + a
            begin, stop = values.indptr[a], values.indptr[b]
            z = (values.data[begin:stop] - loc[start:end][ids]) / safe[ids]
            z = np.where(scale[start:end][ids] > C.MIN_SPREAD, z, 0)
            values.data[begin:stop] = np.where(np.isfinite(z), z, 0) - base_z[ids]
        values.eliminate_zeros()
        _, offset[start:end], variance[start:end] = _moments(
            values, block.present, base_z, fill=np.zeros(width), deviations=True)
        blocks.append((values.tocsr(), block.present, base_z))
        start = end
    operator = CenteredOperator(blocks, offset)
    u, singular, vt = randomized_svd(operator, cap, seed=seed)
    del u
    return C._finish_locus(columns, supports, loc, scale, n_obs, bounds, operator,
                           singular**2/max(n-1, 1), vt.T, n, p,
                           feature_definition=feature_definition, n_components=n_components,
                           total_variance=float(variance.sum()*n/(n-1)), offset=offset)


def feature_matrices(items, vocab, definition, feature_kw, n_jobs):
    """Use the existing sparse motif writer; keep only independent numeric row slots."""
    from . import corpus as C, layout as L
    from .cohort import parallel_rows
    builders, columns, maps = {}, {}, {}
    for locus, vocabulary in vocab.items():
        groups = []
        maps[locus] = {}
        for group in L._feature_groups('vsig', definition):
            if locus not in group.emitted_loci:
                continue
            names = vocabulary[group.name] if group.dynamic else tuple(group.features)
            groups.append((group.name, len(names)))
            maps[locus][group.name] = {name: j for j, name in enumerate(names)}
        builders[locus] = MatrixBuilder(len(items), groups)
        columns[locus] = L.raw_columns('vsig', locus, vocabulary, feature_definition=definition)

    def fill(indexed):
        i, item = indexed
        output = {loc: {name: (row, maps[loc][name]) for name, row in builder.row(i).items()}
                  for loc, builder in builders.items()}
        C._row_for(item, 'vsig', vocab, {**feature_kw, '_raw_array': output})

    parallel_rows(enumerate(items), fill, n_jobs)
    matrices = {loc: (builder.finish(), columns[loc]) for loc, builder in builders.items()}
    # glibc retains freed per-sample buffers in worker arenas after assembly.
    # Return those pages before fitting; other allocators still use the same RAM guard.
    if sys.platform == 'linux':
        trim = getattr(CDLL(None), 'malloc_trim', None)
        if trim is not None:
            trim(0)
    return matrices


def trimmed_fraction(x, lo, hi):
    moved = seen = start = 0
    for block in x.blocks:
        values = block.values
        observed = int(block.present.sum())
        for a, b in _spans(values):
            ids = _indices(values, a, b)
            data = values.data[values.indptr[a]:values.indptr[b]]
            lower, upper = lo[start+a:start+b], hi[start+a:start+b]
            finite = np.isfinite(data)
            zero = observed - np.diff(values.indptr[a:b+1])
            seen += int(observed*(b-a) - (~finite).sum())
            moved += int(zero[np.clip(np.zeros(b-a), lower, upper) != 0].sum())
            moved += int((finite & (np.clip(data, lower[ids], upper[ids]) != data)).sum())
        start += values.shape[1]
    return float(moved/seen) if seen else 0.0


def table_pca(x, k, *, seed):
    """Median-impute/StandardScaler semantics, with a bounded centered operator."""
    from ..cores import require_memory
    from sklearn import __version__
    from sklearn.utils.extmath import svd_flip
    n, p = x.shape
    k = min(k, n, p)
    # Match PCA(auto): full small/near-full fits, or covariance fits since sklearn 1.5.
    modern = tuple(map(int, __version__.split('.')[:2])) >= (1, 5)
    covariance = modern and p <= 1000 and n >= 10*p
    full = covariance or max(n, p) <= 500 or k >= .8*min(n, p)
    sketch_width = k if covariance else (min(n, p) if full else k+10)
    iterations = 0 if full else (7 if k < .1*min(n, p) else 4)
    require_memory(3*x.nbytes + 64*(n+p)*sketch_width + 8*CHUNK_BYTES + 96*p + (16*p*p if covariance else 0),
                   'bounded sample feature PCA workspace')
    fill, mean, scale, center, variance = (np.zeros(p) for _ in range(5))
    blocks, start = [], 0
    eps = np.finfo(float).eps
    for block in x.blocks:
        width = block.values.shape[1]
        end = start+width
        _, quantiles = sparse_quantiles(block, [.5])
        median = np.nan_to_num(quantiles[0])
        fill[start:end] = median
        _, means, vs = _moments(block.values, block.present, np.zeros(width), fill=median)
        # StandardScaler's roundoff bound identifies constant columns.
        scales = np.sqrt(vs)
        scales[vs <= n*eps*vs + (n*means*eps)**2] = 1
        mean[start:end], scale[start:end] = means, scales
        values = block.values.copy()
        for a, b in _spans(values):
            ids = _indices(values, a, b) + a
            begin, stop = values.indptr[a], values.indptr[b]
            data = values.data[begin:stop]
            values.data[begin:stop] = np.where(np.isfinite(data), data, median[ids]) / scales[ids]
        values.eliminate_zeros()
        # All rows begin at the imputed value; presence removes that baseline.
        base = -median/scales
        missing = (median-means)/scales
        _, local_mean, local_variance = _moments(values, block.present, base,
                                                fill=np.zeros(width), deviations=True)
        center[start:end] = local_mean + missing
        variance[start:end] = local_variance*n/(n-1)
        blocks.append((values.tocsr(), block.present, base))
        start = end
    operator = CenteredOperator(blocks, center-(fill-mean)/scale)
    if covariance:
        # PCA's tall, <=1000-column covariance regime needs only a small Gram matrix.
        gram = np.empty((p, p))
        for start in range(0, p, _width(n)):
            stop = min(p, start+_width(n))
            basis = np.zeros((p, stop-start))
            basis[np.arange(start, stop), np.arange(stop-start)] = 1
            gram[:, start:stop] = operator.T @ (operator @ basis)
        eigenvalues, vectors = linalg.eigh(gram)
        singular = np.sqrt(np.maximum(eigenvalues[-k:][::-1], 0))
        from .corpus import _fix_signs
        vt = _fix_signs(vectors[:, -k:][:, ::-1]).T
    else:
        u, singular, vt = randomized_svd(operator, k, seed=seed, n_iter=iterations, oversamples=sketch_width-k)
        if modern or full:
            _, vt = svd_flip(u, vt, u_based_decision=not modern)
        del u
    total = variance.sum()
    ratio = singular**2/(n-1)/total if total > 0 else np.zeros(k)
    return fill, mean, scale, center, vt, ratio
