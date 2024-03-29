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

package com.antigenomics.vdjtools.basic

import com.antigenomics.vdjtools.sample.Clonotype

abstract class Histogram {
    protected final boolean unweighted
    protected final int min, max, span

    protected final double[] innerHist
    protected final int[] lengths

    protected int count = 0
    protected double freq = 0

    public final String HEADER

    protected Histogram(boolean unweighted, int min, int max) {
        this.unweighted = unweighted

        this.min = min
        this.max = max

        this.span = max - min + 1
        this.innerHist = new double[span]
        this.lengths = (min..max) as double[]
        this.HEADER = lengths.collect().join("\t")
    }

    protected abstract int getValue(Clonotype clonotype)

    public int bin(Clonotype clonotype) {
        Math.min(max, Math.max(min, getValue(clonotype))) - min
    }

    /**
     * Update histogram with a single clonotype 
     * @param clonotype clonotype to add
     */
    public void add(Clonotype clonotype) {
        if (unweighted) {
            this.innerHist[bin(clonotype)]++
            freq++
        } else {
            this.innerHist[bin(clonotype)] += clonotype.freq
            freq += clonotype.freq
        }
        count++
    }

    /**
     * Update histogram with a set of clonotypes 
     * @param sample sample to add
     */
    public void addAll(Iterable<Clonotype> sample) {
        sample.each { Clonotype clonotype ->
            add(clonotype)
        }
    }

    /**
     * Combines two histograms
     * @param other histogram to update with
     */
    public void addAll(Histogram other) {
        if (other.unweighted != this.unweighted ||
                other.min != this.min || other.max != this.max)
            throw new Exception("Can't combine histograms of different type")

        this.count += other.count
        this.freq += other.freq
        other.innerHist.eachWithIndex { it, ind ->
            innerHist[ind] += it
        }
    }

    /**
     * Clears all spectratype counters
     */
    public void clear() {
        for (int i = 0; i < span; i++)
            this.innerHist[i] = 0
        count = 0
        freq = 0
    }

    /**
     * Normalize and get the histogram. Values correspond to specified clonotype feature, ranging from {@code min} to {@code max}
     * @return
     */
    public double[] getHistogram() {
        getHistogram(true)
    }

    /**
     * Normalize and get the histogram array. Values correspond to specified clonotype feature, ranging from {@code min} to {@code max}
     * @param normalized if true will normalize clonotype to total sum of 1.
     * @return
     */
    public double[] getHistogram(boolean normalized) {
        def histogram = new double[span]

        double _freq = (normalized && freq > 0) ? freq : 1.0

        for (int i = 0; i < span; i++)
            histogram[i] = this.innerHist[i] / _freq

        histogram
    }

    /**
     * Gets min CDR3 length of clonotype (lower boundary on spectratype length)
     * @return
     */
    int getMin() {
        min
    }

    /**
     * Gets max CDR3 length of clonotype (upper boundary on spectratype length)
     * @return
     */
    int getMax() {
        max
    }

    /**
     * Gets the number of bins in spectratype histogram, i.e. {@code max - min + 1}
     * @return
     */
    public int getSpan() {
        span
    }

    /**
     * Gets total frequency of summarized clonotypes
     * @return
     */
    public double getFreq() {
        freq
    }

    /**
     * Gets total read count of summarized clonotypes
     * @return
     */
    public int getCount() {
        count
    }

    /**
     * An array of CDR3 lengths that correspond to spectratype bins
     */
    public int[] getLengths() {
        return lengths
    }

    boolean isUnweighted() {
        unweighted
    }

    /**
     * Plain text row for tabular output
     */
    @Override
    public String toString() {
        histogram.collect().join("\t")
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
