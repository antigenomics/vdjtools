/*
 * Copyright (c) 2014-2024, OOO «MiLaboratory»
 *
 * IN NO EVENT SHALL THE INVENTORS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
 * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
 * ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF THE INVENTORS HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE LICENSOR HAS NO
 * OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
 * MODIFICATIONS. THE LICENSOR MAKES NO REPRESENTATIONS AND EXTENDS NO
 * WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR
 * PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY PATENT,
 * TRADEMARK OR OTHER RIGHTS.
 */

package com.antigenomics.vdjtools.diversity

import com.antigenomics.vdjtools.misc.AtomicDouble
import com.antigenomics.vdjtools.misc.ExecUtil
import com.antigenomics.vdjtools.misc.MathUtil
import groovyx.gpars.GParsPool

import static com.antigenomics.vdjtools.diversity.RichnessEstimateType.*

/**
 * Species richness estimates based on works of Anne Chao et al. Implemented according to the following paper:
 * {@url http://viceroy.eeb.uconn.edu/estimates/EstimateSPages/EstSUsersGuide/References/ColwellEtAl2012.pdf}
 */
class ChaoEstimator {
    private final FrequencyTable frequencyTable
    private final long n
    private final double F0, Sobs, F1, F2

    /**
     * Creates an instance of class that computes Chao richness estiamtes.
     * @param frequencyTable a {@link com.antigenomics.vdjtools.diversity.FrequencyTable} summary for sample of interest.
     */
    ChaoEstimator(FrequencyTable frequencyTable) {
        this.frequencyTable = frequencyTable
        this.n = frequencyTable.count
        this.Sobs = frequencyTable.diversity
        this.F1 = frequencyTable.singletons
        this.F2 = frequencyTable.doubletons
        this.F0 = F1 * (F1 - 1) / 2 / (F2 + 1)
    }

    /**
     * Computes the estimate for lower bound of total repertoire richness,
     * i.e. the total number of clonotypes in individual that was sampled.
     * @return Chao lower bound estimate for total richness.
     */
    SpeciesRichness chao1() {
        new SpeciesRichness(
                Sobs + F0,
                Math.sqrt(
                        F0 + F1 * (2 * F1 - 1) * (2 * F1 - 1) / 4 / (F2 + 1) / (F2 + 1) +
                                F1 * F1 * F2 * (F1 - 1) * (F1 - 1) / 4 / (F2 + 1) / (F2 + 1) / (F2 + 1) / (F2 + 1)
                ),
                frequencyTable.count,
                TotalDiversityLowerBoundEstimate)
    }

    /**
     * Extrapolates observed richness according to multinomial model.
     * @param extrapolateTo number of reads, should be greater than the total number of reads in a sample
     * @return extrapolated richness estimate.
     */
    SpeciesRichness chaoE(long extrapolateTo) {
        double mStar = extrapolateTo - n

        if (mStar < 0)
            throw new IllegalArgumentException("Should extrapolate farther than the total count of the sample")

        // here is the derivative to compute variance
        double dF0dF1 = (2 * F1 - 1) / 2 / (F2 + 1),
               dF0dF2 = -F1 * (F1 - 1) / 2 / (F2 + 1) / (F2 + 1)
        double brackets = 1.0 - Math.pow(1.0 - F1 / n / F0, mStar),
               dBrackets = -mStar * Math.pow(1.0 - F1 / n / F0, mStar - 1),
               dBracketsdF1 = dBrackets * (-1 / n / F0 + F1 / n / F0 / F0 * dF0dF1),
               dBracketsdF2 = dBrackets * (F1 / n / F0 / F0 * dF0dF2)
        double dSdF1 = dF0dF1 * brackets + F0 * dBracketsdF1,
               dSdF2 = dF0dF2 * brackets + F0 * dBracketsdF2,
               cov11 = F1 * (1 - F1 / (Sobs + F0)), cov22 = F2 * (1 - F2 / (Sobs + F0)),
               cov12 = -F1 * F2 / (Sobs + F0)

        new SpeciesRichness(
                Sobs + F0 * brackets,
                Math.sqrt(Sobs * (1.0 - Sobs / (Sobs + F0)) +
                        dSdF1 * dSdF1 * cov11 + dSdF2 * dSdF2 * cov22 + 2 * dSdF1 * dSdF2 * cov12),
                extrapolateTo,
                Extrapolated)
    }

    /**
     * Interpolates observed richness based on multinomial model.
     * @param interpolateTo number of reads, should be less than the total number of reads in a sample
     * @return interpolated richness estimate.
     */
    SpeciesRichness chaoI(long interpolateTo) {
        if (n > Integer.MAX_VALUE)
            throw new UnsupportedOperationException()

        if (interpolateTo < 0 || interpolateTo > n)
            throw new IllegalArgumentException("Should interpolate within the size of sample")

        //double denom = CombinatoricsUtils.binomialCoefficientLog((int) n, (int) interpolateTo)
        final int n = (int) n, m = (int) interpolateTo
        final double denom = MathUtil.logFactorialRatio(n, m)
        final AtomicDouble sum1 = new AtomicDouble(), sum2 = new AtomicDouble()

        GParsPool.withPool ExecUtil.THREADS, { // this is quite time consuming
            frequencyTable.bins.eachParallel { FrequencyTable.FrequencyTableBin bin ->
                int k = bin.count,
                    f = bin.diversity

                if (k <= n - m) {
                    double alpha = Math.exp(MathUtil.logFactorialRatio(n - k, m) - denom)
                    //def alpha = Math.exp(CombinatoricsUtils.binomialCoefficientLog((int) n - k,
                    //        (int) interpolateTo) - denom)
                    sum1.addAndGet(f * alpha)
                    sum2.addAndGet(f * (1 - alpha) * (1 - alpha))
                } else {
                    sum2.addAndGet(f)
                }
            }
        }

        double Sind = Sobs - sum1.get()

        new SpeciesRichness(
                Sind,
                Math.sqrt(sum2.get() - Sind * Sind / (Sobs + F0)),
                m,
                m == n ? Observed : Interpolated)
    }
}
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
