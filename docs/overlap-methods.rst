Exact cohort distance calculations
==================================

For four or more samples, ``pairwise_distances`` batches F/F2/D/Jaccard through
one shared clonotype vocabulary and sparse products. Up to three samples, it keeps
pairwise joins: there are at most three pairs, and a global vocabulary adds
unnecessary work for a few deep samples. Both paths aggregate each sample once.
The batch uses existing SciPy kernels, with no new native extension or fitted parameters. Pearson ``R`` retains the centered pairwise
calculation to avoid cancellation in a raw-moment covariance formula. Signed or
non-finite frequency tables also retain their existing pairwise arithmetic.
F/F2 batching uses canonical Int64 counts with effective read totals below
:math:`2^{53}`; other count types and larger totals retain the pairwise path.

For :math:`N` samples and :math:`K` fully specified clonotype keys, let
:math:`B_{ik}` be one when sample :math:`i` contains key :math:`k`, and zero otherwise.
Let :math:`W_{ik}` be that key's within-sample read frequency, including an explicit
zero for a zero-count clonotype. Frequencies retain their original sample totals;
null-key rows cannot match another sample but still contribute to those totals and
the sample richness :math:`n_i`. Thus missing calls are not silently matched.

The implementation in ``overlap.cluster._exact_matrix`` computes

.. math::

   C = BB^T,\qquad A = WB^T,\qquad H = \sqrt{W}\sqrt{W}^{\,T},

where the square root is elementwise. Entry :math:`C_{ij}` counts shared keys,
:math:`A_{ij}` sums sample :math:`i`'s frequencies on those keys, and
:math:`H_{ij}` sums their geometric mean frequencies. Therefore the similarities are

.. math::

   D_{ij}=\frac{C_{ij}}{n_i n_j},\qquad
   F_{ij}=\sqrt{A_{ij} A_{ji}},\qquad
   F2_{ij}=H_{ij},\qquad
   J_{ij}=\frac{C_{ij}}{n_i+n_j-C_{ij}}.

A zero denominator gives similarity zero. The existing distance transforms remain
:math:`-\log_{10}(s+10^{-9})` for D/F/F2 and :math:`1-s` for Jaccard, where
:math:`s` denotes the corresponding similarity. The diagonal is set to zero and
both output forms preserve sample order.

The code accumulates read counts before normalization. Let :math:`Q_{ik}` be the
nonnegative integer read count of fully specified key :math:`k`, and let
:math:`t_i` be the sample's total read count, including null-key rows. A zero total
uses :math:`t_i=1`, matching the single-pair normalization. Then

.. math::

   A_{ij}=\frac{(QB^T)_{ij}}{t_i},\qquad
   H_{ij}=\frac{(\sqrt{Q}\sqrt{Q}^{\,T})_{ij}}{\sqrt{t_i}\sqrt{t_j}}.

For nonnegative integer totals below :math:`2^{53}`, float64 shared-count
accumulation is exact. This avoids accumulating millions of small normalized
frequencies for F. F2 still rounds its square roots and sums. Vocabulary indices
follow first occurrence, keeping accumulation order repeatable for the same inputs.

These are algebraic identities, not new estimators. Floating-point summation order
can differ. For a sum of :math:`m` nonnegative terms in the normal floating-point
range, the usual accumulation bound is :math:`\gamma_m=mu/(1-mu)` when
:math:`mu<1`, with float64 unit roundoff :math:`u=2^{-53}`; products and square roots
add their own rounding. A similarity perturbation from :math:`x` to :math:`y`
changes its log distance by at most
:math:`|x-y|/[\ln(10)(\min(x,y)+10^{-9})]`. Bitwise equality with pairwise joins is
not promised. Boundary tests cover empty samples, zero counts, duplicate keys,
null calls, custom keys and signed inputs.

With :math:`P` aggregated sample/key rows, :math:`p_i` fully specified keys in
sample :math:`i`, and :math:`d_k` samples containing key :math:`k`, expected
fixed-length-key construction, sparse row sorting and product work is
:math:`O(P+\sum_i p_i\log(p_i+1)+\sum_k d_k^2+N^2)`. The returned dense distances inherently require
:math:`O(N^2)` storage; numerical workspace is checked against the available memory
allocation. Matrix and long-form outputs use native columns rather than constructing
one Python float or tuple per result cell. Large all-pairs outputs still require
adequate memory or a smaller sample collection.

