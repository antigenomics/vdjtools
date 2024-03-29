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

import com.antigenomics.vdjtools.overlap.OverlapType
import com.antigenomics.vdjtools.pool.PooledSample
import com.antigenomics.vdjtools.pool.SampleAggregator
import com.antigenomics.vdjtools.sample.Sample
import groovy.transform.PackageScope
import org.apache.commons.math3.util.CombinatoricsUtils

/**
 * Class that computes richness estimates and diversity indices. 
 * All computations are performed via the {@link com.antigenomics.vdjtools.diversity.FrequencyTable} object.
 * Note that the clonotype is always computed based on certain {@link OverlapType},
 * which tells how to collapse clonotypes, e.g. consider identical CDR3 nucleotide or amino acid sequences, etc.
 * Therefore in some cases the result will be different from one obtained using raw clonotype frequencies.
 * This will not happen in case {@link OverlapType#Strict} is used,
 * which is recommended for most purposes.
 *
 * @see ResamplingEstimator
 * @see com.antigenomics.vdjtools.diversity.DiversityEstimate
 */
class ExactEstimator extends DiversityEstimator {
    private final long extrapolateTo
    private final ChaoEstimator chaoEstimator

    /**
     * INTERNAL, constructor used for {@link ResamplingEstimator.}
     */
    @PackageScope
    static ExactEstimator basicDiversityEstimates(FrequencyTable frequencyTable) {
        new ExactEstimator(frequencyTable, -1)
    }

    /**
     * Creates an instance of individual-based diversity estimates class.
     * @param frequencyTable a {@link com.antigenomics.vdjtools.diversity.FrequencyTable} summary for sample of interest.
     * @param extrapolateTo desired extrapolation extent. Used to compute the
     * {@link com.antigenomics.vdjtools.diversity.ChaoEstimator#chaoE} estimate.
     *                      For most cases, it should be set to the size of the largest sample when several samples are to be compared.
     */
    ExactEstimator(FrequencyTable frequencyTable, long extrapolateTo) {
        super(frequencyTable, EstimationMethod.Exact)
        this.extrapolateTo = extrapolateTo
        this.chaoEstimator = new ChaoEstimator(frequencyTable)
    }

    /**
     * Creates an instance of individual-based diversity estimates class.
     * Will compute {@link com.antigenomics.vdjtools.diversity.FrequencyTable} for a given sample.
     * @param sample a sample to analyze.
     * @param intersectionType {@link OverlapType} that will be used
     *                         to collapse sample during {@link com.antigenomics.vdjtools.diversity.FrequencyTable} computation
     * @param extrapolateTo desired extrapolation extent. Used to compute the
     * {@link com.antigenomics.vdjtools.diversity.ChaoEstimator#chaoE} estimate.
     *                      For most cases, it should be set to the size of the largest sample when several samples are to be compared.
     */
    ExactEstimator(Sample sample, OverlapType intersectionType, long extrapolateTo) {
        this(new FrequencyTable(sample, intersectionType), extrapolateTo)
    }

    /**
     * Creates an instance of sample-based diversity estimates class.
     * Will summarize clonotype occurrences to build {@link com.antigenomics.vdjtools.diversity.FrequencyTable}
     * @param pool a pool of several samples to analyze.
     * @param extrapolateTo desired extrapolation extent. Used to compute the
     * {@link com.antigenomics.vdjtools.diversity.ChaoEstimator#chaoE} estimate.
     *                      For most cases, it should be set to the number of samples in the largest sample pool if several are to be compared.
     */
    ExactEstimator(PooledSample pool, long extrapolateTo) {
        this(new FrequencyTable(pool), extrapolateTo)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    DiversityIndex getShannonWienerIndex() {
        // todo: std computaiton ?
        def mean = 0
        frequencyTable.bins.each { FrequencyTable.FrequencyTableBin bin ->
            def h = bin.diversity * bin.freq * Math.log(bin.freq)
            mean -= h
        }
        new DiversityIndex(Math.exp(mean), 0, frequencyTable.count)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    DiversityIndex getNormalizedShannonWienerIndex() {
        def mean = 0
        def totalDiversity = 0

        frequencyTable.bins.each { FrequencyTable.FrequencyTableBin bin ->
            def h = bin.diversity * bin.freq * Math.log(bin.freq)
            mean -= h
            totalDiversity += bin.diversity
        }

        new DiversityIndex(mean / Math.log(totalDiversity), 0, frequencyTable.count)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    DiversityIndex getInverseSimpsonIndex() {
        // todo: std computaiton ?
        def mean = 0
        frequencyTable.bins.each { FrequencyTable.FrequencyTableBin bin ->
            mean += bin.diversity * bin.freq * bin.freq
        }
        //std -= mean * mean
        new DiversityIndex(1.0 / mean, 0, frequencyTable.count)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    DiversityIndex getDxxIndex(double fraction) {
        if (fraction < 0 || fraction > 1)
            throw new IllegalArgumentException("Fraction value should be within [0,1] bounds.")

        def div = 0, freqSum = 0

        frequencyTable.bins.find {
            freqSum += it.freq
            div++

            freqSum >= fraction
        }

        new DiversityIndex(1.0 - div / (double) frequencyTable.diversity, 0, frequencyTable.count)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    SpeciesRichness getEfronThisted(int maxDepth, double cvThreshold) {
        double S = -1, D = -1, CV

        for (int depth = 1; depth <= maxDepth; depth++) {
            final double[] h = new double[depth], nx = new double[depth]
            for (int y = 1; y <= depth; y++) {
                nx[y - 1] = frequencyTable[y]

                // Calculate Euler coefficients
                for (int x = 1; x <= y; x++) {
                    def coef = CombinatoricsUtils.binomialCoefficientDouble((int) (y - 1), (int) (x - 1))
                    if (x % 2 == 1)
                        h[x - 1] += coef
                    else
                        h[x - 1] -= coef
                }
            }

            // Extrapolate count
            S = frequencyTable.diversity + (double) (0..<depth).sum { int i -> h[i] * nx[i] }
            D = Math.sqrt((double) (0..<depth).sum { int i -> h[i] * h[i] * nx[i] })
            CV = D / S

            // Go to maximum count depth, but balance that STD doesn't get too high
            if (CV >= cvThreshold)
                break
        }

        new SpeciesRichness((long) S, (long) D, frequencyTable.count, RichnessEstimateType.TotalDiversityLowerBoundEstimate)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    SpeciesRichness getChao1() {
        chaoEstimator.chao1()
    }

    /**
     * {@inheritDoc}
     */
    @Override
    DiversityEstimate getChaoE() {
        extrapolateTo < frequencyTable.count ? DiversityEstimate.DUMMY : chaoEstimator.chaoE(extrapolateTo)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    SpeciesRichness getObservedDiversity() {
        new SpeciesRichness(frequencyTable.diversity, 0, frequencyTable.count, RichnessEstimateType.Observed)
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
