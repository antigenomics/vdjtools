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
 * A species richness estimate that reflects the total number of clonotypes in an immune repertoire.
 */
class SpeciesRichness extends DiversityEstimate {
    private final RichnessEstimateType type

    /**
     * Creates a structure holding diversity estimate summary.
     * @param mean expected value of a diversity estimate.
     * @param std standard deviation of a diversity estimate.
     * @param numberOfReads number of reads in the sample that was analyzed.
     * @param type richness estimate type. 
     */
    SpeciesRichness(long mean, long std, long numberOfReads, RichnessEstimateType type) {
        super(mean, std, numberOfReads)
        this.type = type
    }

    /**
     * Creates a structure holding diversity estimate summary.
     * @param mean expected value of a diversity estimate.
     * @param std standard deviation of a diversity estimate.
     * @param numberOfReads number of reads in the sample that was analyzed.
     * @param type richness estimate type.
     */
    SpeciesRichness(double mean, double std, long numberOfReads, RichnessEstimateType type) {
        this((long) mean, (long) std, numberOfReads, type)
    }

    /**
     * Gets the richness estimate type.
     * @return estimate type (interpolated/observed/extrapolated/lower bound estimate on total diversity).
     */
    RichnessEstimateType getType() {
        type
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
