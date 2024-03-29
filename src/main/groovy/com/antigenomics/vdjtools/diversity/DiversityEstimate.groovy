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
 * An estimate of repertoire diversity computed based on the abundances of sampled clonotypes.
 */
class DiversityEstimate {
    protected final mean, std
    protected final long size

    public static DiversityEstimate DUMMY = new DiversityEstimate("NA", "NA", -1)

    /**
     * Creates a structure holding diversity estimate summary.
     * @param mean expected value of a diversity estimate.
     * @param std standard deviation of a diversity estimate.
     * @param size size of the sample that was analyzed.
     */
    DiversityEstimate(mean, std, long size) {
        this.mean = mean
        this.std = std
        this.size = size
    }

    /**
     * Gets the mean value of diversity estimate. 
     * @return mean value.
     */
    def getMean() {
        mean
    }

    /**
     * Gets the standard deviation of diversity estimate. 
     * @return standard deviation.
     */
    def getStd() {
        std
    }

    /**
     * Gets the size of the sample that was analyzed.
     * @return sample size.
     */
    long getSize() {
        size
    }

    /**
     * Plain text row for tabular output
     */
    @Override
    public String toString() {
        [mean, std].join("\t")
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
