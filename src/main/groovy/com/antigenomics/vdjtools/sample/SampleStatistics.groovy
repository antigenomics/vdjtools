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

package com.antigenomics.vdjtools.sample

/**
 * An object holding some basic statistics for a sample collection 
 */
public class SampleStatistics {
    private final long minCount, maxCount
    private final double minFreq, maxFreq
    private final int minDiversity, maxDiversity

    /**
     * Initializes a new sample statistics object with pre-computed values 
     * @param minCount minimal number of reads in a sample from given sample collection
     * @param maxCount maximal number of reads in a sample from given sample collection
     * @param minFreq minimal total frequency of clonotypes in a sample from given sample collection
     * @param maxFreq maximal total frequency of clonotypes in a sample from given sample collection
     * @param minDiversity minimal number of clonotypes in a sample from given sample collection
     * @param maxDiversity maximal number of clonotypes in a sample from given sample collection
     */
    public SampleStatistics(long minCount, long maxCount,
                            double minFreq, double maxFreq,
                            int minDiversity, int maxDiversity) {
        this.minCount = minCount
        this.maxCount = maxCount
        this.minFreq = minFreq
        this.maxFreq = maxFreq
        this.minDiversity = minDiversity
        this.maxDiversity = maxDiversity
    }

    /**
     * Gets the minimal number of reads in a sample from given sample collection
     * @return size of the smallest sample in sample collection
     */
    public long getMinCount() {
        minCount
    }

    /**
     * Gets the maximal number of reads in a sample from given sample collection
     * @return size of the largest sample in sample collection
     */
    public long getMaxCount() {
        maxCount
    }

    public double getMinFreq() {
        minFreq
    }

    public double getMaxFreq() {
        maxFreq
    }

    public int getMinDiversity() {
        minDiversity
    }

    public int getMaxDiversity() {
        maxDiversity
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
