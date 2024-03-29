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
/**
 * Base class for implementations that compute various species richness estimates and diversity indices.
 * @see com.antigenomics.vdjtools.diversity.DiversityIndex
 * @see com.antigenomics.vdjtools.diversity.SpeciesRichness
 */
abstract class DiversityEstimator {
    protected final FrequencyTable frequencyTable
    protected final EstimationMethod estimationMethod

    /**
     * Protected constructor.
     * @param frequencyTable a {@link com.antigenomics.vdjtools.diversity.FrequencyTable} summary for sample of interest.
     * @param estimationMethod specifies the method that is used in this estimator implementation. 
     */
    protected DiversityEstimator(FrequencyTable frequencyTable, EstimationMethod estimationMethod) {
        this.frequencyTable = frequencyTable
        this.estimationMethod = estimationMethod
    }

    /*
     * Diversity indices
     */

    /**
     * Computes the Shannon-Wiener diversity index.
     * @return {@link com.antigenomics.vdjtools.diversity.DiversityEstimate} object handling the result.
     */
    abstract DiversityIndex getShannonWienerIndex()

    /**
     * Computes the normalized Shannon-Wiener diversity index.
     * @return {@link com.antigenomics.vdjtools.diversity.DiversityEstimate} object handling the result.
     */
    abstract DiversityIndex getNormalizedShannonWienerIndex()

    /**
     * Computes the Inverse Simpson diversity index.
     * @return {@link com.antigenomics.vdjtools.diversity.DiversityEstimate} object handling the result.
     */
    abstract DiversityIndex getInverseSimpsonIndex()

    /**
     * Computes the diversity index DXX, where XX is a specified fraction.
     * The estimate equals to {@code 1 - n / N}, where {@code n} is the minimum number of clonotypes accounting 
     * for at least XX% of the total reads and {@code N} is the total number of clonotypes in the sample.
     * @param fraction a fraction of reads at which the diversity estimate will be computed.
     * @return {@link com.antigenomics.vdjtools.diversity.DiversityEstimate} object handling the result.
     */
    abstract DiversityIndex getDxxIndex(double fraction)

    /**
     * Computes diversity index that equals to one minus the minimum fraction of clonotypes accounting 
     * for at least 50% of the total reads.
     * @return {@link com.antigenomics.vdjtools.diversity.DiversityEstimate} object handling the result.
     * @see com.antigenomics.vdjtools.diversity.DiversityEstimator#getDxxIndex
     */
    DiversityEstimate getD50Index() {
        getDxxIndex(0.5)
    }

    /*
     * Richness estimates
     */

    /**
     * Computes the Efron-Thisted lower bound estimate of total diversity.
     * @param maxDepth max number of terms in Euler series to compute.
     * @param cvThreshold coefficient of variance threshold, after reaching it higher order terms will be ignored.
     * @return {@link com.antigenomics.vdjtools.diversity.DiversityEstimate} object handling the result.
     */
    abstract SpeciesRichness getEfronThisted(int maxDepth, double cvThreshold)

    /**
     * Computes the Efron-Thisted lower bound estimate of total diversity.
     * @return {@link com.antigenomics.vdjtools.diversity.DiversityEstimate} object handling the result.
     */
    SpeciesRichness getEfronThisted() {
        getEfronThisted(20, 0.05)
    }

    /**
     * Computes the Chao1 lower bound estimate of total diversity.
     * @return {@link com.antigenomics.vdjtools.diversity.DiversityEstimate} object handling the result.
     */
    abstract SpeciesRichness getChao1()

    /**
     * Computes the extrapolated Chao diversity estimate.
     * @return {@link com.antigenomics.vdjtools.diversity.DiversityEstimate} object handling the result.
     */
    abstract DiversityEstimate getChaoE()

    /**
     * Computes the observed diversity, i.e. the number of clonotypes in the sample.
     * @return {@link com.antigenomics.vdjtools.diversity.DiversityEstimate} object handling the result.
     */
    abstract SpeciesRichness getObservedDiversity()

    /*
     * end estimates
     */

    /**
     * Gets the underlying frequency table.
     * @return {@link com.antigenomics.vdjtools.diversity.DiversityEstimate} object handling the result.
     */
    FrequencyTable getFrequencyTable() {
        frequencyTable
    }

    /**
     * Gets the computation method that is used in this estimator implementation.
     * @return estimation method (exact/resampling).
     */
    EstimationMethod getEstimationMethod() {
        return estimationMethod
    }

    /**
     * List of fields that will be included in tabular output .
     */
    static final String[] ESTIMATE_NAMES = ["observedDiversity",
                                            "chaoE",
                                            "efronThisted", "chao1",
                                            "d50Index",
                                            "shannonWienerIndex",
                                            "normalizedShannonWienerIndex",
                                            "inverseSimpsonIndex"]

    /**
     * Computes all the diversity estimates.
     * @return a map {@code [name -> value]} of computed diversity estimates.
     */
    Map<String, DiversityEstimate> computeAll() {
        ESTIMATE_NAMES.collectEntries { [(it): this."$it"] }
    }

    /**
     * Header string, used for tabular output.
     */
    static final String HEADER = ESTIMATE_NAMES.collect { [it + "_mean", it + "_std"] }.flatten().join("\t")

    /**
     * Plain text row for tabular output.
     */
    @Override
    String toString() {
        ESTIMATE_NAMES.collect { this."$it" }.join("\t")
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
