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

package com.antigenomics.vdjtools.pool;

import com.antigenomics.vdjtools.misc.MathUtil;
import com.antigenomics.vdjtools.sample.Clonotype;

/**
 * A base class for all clonotype aggregators used for summarizing clonotype counts and occurrences.
 * Representative clonotype for the clonotype aggregatior is implementation-dependent.
 */
public abstract class ClonotypeAggregator {
    private int sampleId;
    private int incidenceCount, count;
    private double freq, freqRem = 0;

    /**
     * Creates a new conotype aggregator.
     *
     * @param clonotype base clonotype instance.
     * @param sampleId  index of clonotype's parent sample.
     */
    public ClonotypeAggregator(Clonotype clonotype, int sampleId) {
        this.incidenceCount = 1;
        this.freqRem = clonotype.getFreq();
        this.count = (int) clonotype.getCount();
        this.sampleId = sampleId;
    }

    final void combine(Clonotype other, int sampleId) {
        if (this.sampleId != sampleId) {
            this.sampleId = sampleId;
            incidenceCount++;
            freq += freqRem > 0 ? Math.log10(freqRem) : 0;  // add remainder
            freqRem = other.getFreq();                      // start accumulating freq
        } else {
            freqRem += other.getFreq();                     // continue accumulating freq
        }

        count += other.getCount();

        tryReplace(other, sampleId);
    }

    protected abstract boolean tryReplace(Clonotype other, int sampleId);

    /**
     * Gets the number of samples this clonotype was detected in.
     *
     * @return number of samples.
     */
    public int getIncidenceCount() {
        return incidenceCount;
    }

    /**
     * Gets the total count of this clonotype in all samples so far.
     *
     * @return clonotype count.
     */
    public long getCount() {
        return count;
    }

    /**
     * INTERNAL used for computing pooled clonotype frequencies.
     */
    public double getFreqLogSum() {
        return freq + (freqRem > 0 ? Math.log10(freqRem) : 0);
    }

    /**
     * INTERNAL used for computing pooled clonotype frequencies.
     */
    public double getFreqGeomMean() {
        return Math.pow(10.0, getFreqLogSum() / incidenceCount);
    }

    /**
     * INTERNAL used for computing pooled clonotype frequencies.
     */
    public double getFreqGeomMean(int numberOfSamples) {
        return Math.pow(10.0, (getFreqLogSum() + (numberOfSamples - incidenceCount) * MathUtil.JITTER_LOG10) / numberOfSamples);
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
