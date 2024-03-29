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

import com.antigenomics.vdjtools.ClonotypeWrapperContainer
import com.antigenomics.vdjtools.misc.AtomicDouble
import com.antigenomics.vdjtools.misc.AtomicDoubleArray

/**
 * A class that computes summary statistics of repertoire clonality divided into several levels:
 *
 * 1. singletons (encountered once), doubletons (encountered twice) and high-order clonotypes - 
 *    those are base quantities to estimate the lower bound on total repertoire diversity.
 *
 * 2. cumulative frequency for several quantiles (e.g. top 25%, next 25%, ...) of high-order clonotypes
 *
 * 3. details for top N clonotypes.
 *
 */
class QuantileStats {
    private final int numberOfQuantiles
    private final AtomicDoubleArray quantileFreqs
    private final AtomicDouble highOrderFreq = new AtomicDouble(),
                               doubletonFreq = new AtomicDouble(),
                               singletonFreq = new AtomicDouble()

    /**
     * Summarizes quantile statistics for a given sample.
     * @param clonotypeContainer a set of clonotypes.
     * @param numberOfQuantiles number of quantiles for 2nd level of detalizaiton.
     */
    QuantileStats(ClonotypeWrapperContainer clonotypeContainer, int numberOfQuantiles) {
        this.numberOfQuantiles = numberOfQuantiles
        this.quantileFreqs = new AtomicDoubleArray(numberOfQuantiles)

        if (!clonotypeContainer.isSorted())
            throw new RuntimeException("Clonotype container should be sorted to be used as input for this statistic")

        update(clonotypeContainer)
    }

    /**
     * Summarizes quantile statisitcs for a given sample.
     * @param clonotypeContainer a set of clonotypes.
     */
    QuantileStats(ClonotypeWrapperContainer clonotypeContainer) {
        this(clonotypeContainer, 5)
    }

    /**
     * Internal - adds more clonotyps to stats.
     */
    private void update(ClonotypeWrapperContainer clonotypeContainer) {
        int n = clonotypeContainer.diversity, m = -1

        for (int i = n - 1; i >= 0; i--) {
            def clonotype = clonotypeContainer[i]
            def count = clonotype.getCount(),
                freq = clonotype.getFreq()
            boolean highOrderFlag = false
            switch (count) {
                case 1:
                    singletonFreq.addAndGet(freq)
                    break
                case 2:
                    doubletonFreq.addAndGet(freq)
                    break
                default:
                    highOrderFlag = true
                    break
            }
            if (highOrderFlag) {
                m = i
                break
            }
        }

        for (int i = 0; i <= m; i++) {
            int q = (i * numberOfQuantiles) / (m + 1)
            def clonotype = clonotypeContainer[i]
            def freq = clonotype.getFreq()
            highOrderFreq.addAndGet(freq)
            quantileFreqs.addAndGet(q, freq)
        }
    }

    /**
     * Gets the number of 2nd level summary quantiles.
     * @return number of 2nd level quantiles.
     */
    int getNumberOfQuantiles() {
        return numberOfQuantiles
    }

    /**
     * Gets frequency for a given quantile.
     * @param quantile quantile index, should be less than {@link #numberOfQuantiles} and greater or equal than {@code 0}.
     * @return selected quantile frequency.
     * @throws IndexOutOfBoundsException wrong quantile index.
     */
    double getQuantileFrequency(int quantile) {
        if (quantile < 0 || quantile >= numberOfQuantiles)
            throw new IndexOutOfBoundsException()
        quantileFreqs.get(quantile)
    }

    /**
     * Gets the frequency of singletons, i.e. clonotypes represented by a single read.
     * @return singleton frequency.
     */
    double getSingletonFreq() {
        singletonFreq.get()
    }

    /**
     * Gets the frequency of doubletons, i.e. clonotypes represented by two reads.
     * @return doubleton frequency.
     */
    double getDoubletonFreq() {
        doubletonFreq.get()
    }

    /**
     * Gets the frequency of high-order clonotypes, i.e. clonotypes represented by more than two reads.
     * @return high-order clonotype frequency.
     */
    double getHighOrderFreq() {
        highOrderFreq.get()
    }

    /**
     * Header string, used for tabular output.
     */
    static final String HEADER = "type\tname\tvalue"

    /**
     * Plain text row for tabular output.
     */
    @Override
    String toString() {
        ["set\t3+\t$highOrderFreq",
         "set\t2\t$doubletonFreq",
         "set\t1\t$singletonFreq",
         (0..<numberOfQuantiles).collect { "quantile\tQ${it + 1}\t${getQuantileFrequency(it)}" }].flatten().join("\n")
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
