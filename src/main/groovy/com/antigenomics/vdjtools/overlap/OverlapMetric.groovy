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

package com.antigenomics.vdjtools.overlap

import static com.antigenomics.vdjtools.overlap.OverlapMetricNormalization.*

/**
 * An enum that defines overlap metric, a function that characterizes the extent of overlap between a pair of samples.
 */
public enum OverlapMetric {
    /**
     * Correlation between sample frequencies of _overlapping_ clonotypes
     */
    /*   */ Correlation("R", R),
    /**
     * Spearman Correlation between sample frequencies of _overlapping_ clonotypes
     */
    /*   */ CorrelationS("Rs", R),
    /**
     * Ratio of observed to expected numbers of unique overlapping clonotypes, {@code div12 / div1 / div2}
     */
            Diversity("D", NegLog),
    /**
     * Geometric mean of sums of frequencies of overlapping clonotypes, {@code sqrt ( freq12 * freq21 )}
     */
            Frequency("F", NegLog),
    /**
     * Sum of geometric means of frequencies of overlapping clonotypes, {@code sum ( sqrt ( freq12(i) * freq21(i) ) ), i = 1..div12}
     */
            Frequency2("F2", NegLog),
    /**
     * Jensen-Shannon divergence between Variable segment usage vectors
     */
            vJSD("vJSD", None),
    /**
     * Jensen-Shannon divergence between concatenated Variable and Joining segment usage vectors
     */
            vjJSD("vjJSD", None),
    /**
     * Jensen-Shannon divergence between flattened Variable-Joining segment pairing matrices
     */
            vj2JSD("vj2JSD", None),
    /**
     * Jensen-Shannon divergence between spectratypes
     */
            sJSD("sJSD", None),
    /**
     * Jaccard index
     */
            Jaccard("Jaccard", Index),
    /*
    ChaoJaccard("ChaoJaccard", Index),
    ChaoSorensen("ChaoSorensen", Index),
    */
            /**
             * MorisitaHorn index
             */
            MorisitaHorn("MorisitaHorn", Index)

    public final String shortName
    public final OverlapMetricNormalization normalization

    /**
     * Defines a new overlap metric
     * @param shortName short name
     * @param normalization normalization type
     */
    public OverlapMetric(String shortName, OverlapMetricNormalization normalization) {
        this.shortName = shortName
        this.normalization = normalization
    }

    /**
     * Gets {@code IntersectMetric} by short name
     * @param shortName short name
     * @return
     */
    public static OverlapMetric getByShortName(String name) {
        name = name.toUpperCase()
        values().find { it.shortName.toUpperCase() == name }
    }

    /**
     * A list of existing {@code IntersectMetric} short names
     */
    public static String allowedNames = values().collect { it.shortName }.join(",")
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
