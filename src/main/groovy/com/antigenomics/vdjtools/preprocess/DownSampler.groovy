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

package com.antigenomics.vdjtools.preprocess

import com.antigenomics.vdjtools.sample.Clonotype
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.misc.MathUtil

/**
 * A class that implements down-sampling procedure, i.e.
 * selecting {@code n < N} reads from a given sample with {@code N} reads
 *
 * TODO: rewrite method in plain Java. Known issue: this will not work for samples with getCount() > Integer.MAX_VALUE
 */
public class DownSampler implements Sampler{
    private final Clonotype[] flattenedClonotypes
    private final Sample sample
    private final boolean unweighted

    /**
     * Create a down-sampler for the specified sample 
     * @param sample sample that would be down-sampled
     */
    public DownSampler(Sample sample) {
        this(sample, false)
    }

    /**
     * Create a down-sampler for the specified sample 
     * @param sample sample that would be down-sampled
     * @param unweighted don't weight clonotypes by frequency during sampling 
     */
    public DownSampler(Sample sample, boolean unweighted) {
        if (!unweighted && sample.count > Integer.MAX_VALUE) {
            throw new RuntimeException("Couldn't downsample samples with > ${Integer.MAX_VALUE} cells")
        }

        this.sample = sample
        this.flattenedClonotypes = new Clonotype[unweighted ? sample.diversity : sample.count]
        this.unweighted = unweighted

        int counter = 0
        sample.each {
            if (unweighted) {
                flattenedClonotypes[counter++] = it
            } else {
                for (int i = 0; i < it.count; i++) {
                    flattenedClonotypes[counter++] = it
                }
            }
        }
    }

    /**
     * Gets a specified number of reads from a given sample
     * @param count number of reads (weighted) or clonotypes (unweighted) to take
     * @return a newly create down-sampled sample, or the underlying sample if the number of reads is greater or equal to the sample size
     */
    public Sample reSample(int count) {
        if (unweighted ? count >= sample.diversity : count >= sample.count) {
            return new Sample(sample)
        } else {
            MathUtil.shuffle(flattenedClonotypes)

            def countMap = new HashMap<Clonotype, Integer>() // same as with strict overlap

            if (unweighted) {
                for (int i = 0; i < count; i++) {
                    def clonotype = flattenedClonotypes[i]
                    countMap.put(clonotype, clonotype.count)
                }
            } else {
                for (int i = 0; i < count; i++) {
                    def clonotype = flattenedClonotypes[i]
                    countMap.put(clonotype, (countMap[clonotype] ?: 0) + 1)
                }
            }

            return new Sample(sample, countMap)
        }
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
