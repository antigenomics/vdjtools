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
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.misc.ExecUtil

/**
 * Class that represents Variable and Joining segment frequency (usage) vector and V-J pairing matrix
 */
public class SegmentUsage {
    public static boolean VERBOSE = true

    private final Map<String, double[]> vSegmentUsage = new HashMap<>(),
                                        jSegmentUsage = new HashMap<>(), vjSegmentUsage = new HashMap<>()
    private final Map<String, Integer> sampleIndex = new HashMap<>()
    private Map<String, Double> sortedVSegmTotal = new HashMap<>(), sortedJSegmTotal = new HashMap<>()
    private final int n
    private final boolean unweighted

    /**
     * Creates a SegmentUsage for a sample collection. 
     * Initializes with a sample collection to obtain all segment names that could be encountered 
     * in order to provide same table layout for different samples.
     * All computations are done within constructor.
     * @param sampleCollection sample collection to analyze
     * @param unweighted will count each unique clonotype once if set to true. Will weight each clonotype by its frequency otherwise
     */
    public SegmentUsage(SampleCollection sampleCollection, boolean unweighted) {
        this.n = sampleCollection.size()
        this.unweighted = unweighted
        // ok to use in batch intersect as all samples are pre-loaded
        sampleCollection.eachWithIndex { it, ind -> process(it, ind) }
        summarize()
    }

    /**
     * Creates a SegmentUsage for a sample collection.
     * Initializes with a set of samples to obtain all segment names that could be encountered
     * in order to provide same table layout for different samples.
     * All computations are done within constructor.
     * Internal constructor
     * @param samples a set of sample to analyze
     * @param unweighted will count each unique clonotype once if set to true. Will weight each clonotype by its frequency otherwise
     */
    public SegmentUsage(Sample[] samples, boolean unweighted) {
        this.n = samples.length
        this.unweighted = unweighted
        // ok to use in batch intersect as all samples are pre-loaded
        samples.eachWithIndex { it, ind -> process(it, ind) }
        summarize()
    }

    /**
     * Process a single sample
     * @param sample sample
     * @param index sample index
     */
    private void process(Sample sample, int index) {
        //todo: everywhere
        ExecUtil.report(this, "Processing sample ${sample.sampleMetadata.sampleId}", VERBOSE)
        sample.each { Clonotype clonotype ->
            def vArray = vSegmentUsage[clonotype.v],
                jArray = jSegmentUsage[clonotype.j],
                vjArray = vjSegmentUsage[clonotype.v + "\t" + clonotype.j]

            if (!vArray)
                vSegmentUsage.put(clonotype.v, vArray = new double[n])

            if (!jArray)
                jSegmentUsage.put(clonotype.j, jArray = new double[n])

            if (!vjArray)
                vjSegmentUsage.put(clonotype.v + "\t" + clonotype.j, vjArray = new double[n])

            def increment = unweighted ? 1 : clonotype.freq

            vArray[index] += increment
            jArray[index] += increment
            vjArray[index] += increment
        }
        sampleIndex.put(sample.sampleMetadata.sampleId, index)
    }

    /**
     * Summarize (calculate sums) 
     */
    private void summarize() {
        vSegmentUsage.each {
            sortedVSegmTotal.put(it.key, (double) it.value.collect().sum() ?: 0)
        }
        sortedVSegmTotal = sortedVSegmTotal.sort()
        jSegmentUsage.each {
            sortedJSegmTotal.put(it.key, (double) it.value.collect().sum() ?: 0)
        }
        sortedJSegmTotal = sortedJSegmTotal.sort()
    }

    /**
     * Gets J usage vector, J segment names are ordered as {@code jUsageHeader()}
     * @param sampleId sample id (as in sample collection / sample name)
     * @return
     */
    public double[] jUsageVector(String sampleId) {
        jUsageVector(getSampleIndex(sampleId))
    }

    /**
     * Gets J usage vector, J segment names are ordered as {@code jUsageHeader()}
     * @param sampleIndex sample index, as in initial sample collection/set
     * @return
     */
    public double[] jUsageVector(int sampleIndex) {
        usageVector(jSegmentUsage, sortedJSegmTotal, sampleIndex)
    }

    /**
     * Gets V usage vector, V segment names are ordered as {@code vUsageHeader()}
     * @param sampleId sample id (as in sample collection / sample name)
     * @return
     */
    public double[] vUsageVector(String sampleId) {
        vUsageVector(getSampleIndex(sampleId))
    }

    /**
     * Gets V usage vector, V segment names are ordered as {@code vUsageHeader()}
     * @param sampleIndex sample index, as in initial sample collection/set
     * @return
     */
    public double[] vUsageVector(int sampleIndex) {
        usageVector(vSegmentUsage, sortedVSegmTotal, sampleIndex)
    }

    /**
     * Gets sample index by sample id 
     * @param sampleId sample id (as in sample collection / sample name)
     * @return
     */
    private int getSampleIndex(String sampleId) {
        if (!this.sampleIndex.containsKey(sampleId))
            throw new IllegalArgumentException("$sampleId is not in the sample collection used to build segment usage matrix")
        this.sampleIndex[sampleId]
    }

    /**
     * INTERNAL
     */
    private static double[] usageVector(Map<String, double[]> usageMap, Map<String, Double> totalMap, int sampleIndex) {
        def sampleTotal = usageMap.values().collect { it[sampleIndex] }.sum() ?: 0.0
        totalMap.collect {
            usageMap[it.key][sampleIndex] / (double) (sampleTotal + 1e-7)
        } as double[]
    }

    /**
     * Gets V-J pairing matrix, J segment (row) names are ordered as {@code jUsageHeader()} and
     * V segment (column) names are ordered as {@code vUsageHeader()}
     * @param sampleId sample id (as in sample collection / sample name)
     * @return
     */
    public double[][] vjUsageMatrix(String sampleId) {
        vjUsageMatrix(getSampleIndex(sampleId))
    }

    /**
     * Gets V-J pairing matrix, J segment (row) names are ordered as {@code jUsageHeader()} and
     * V segment (column) names are ordered as {@code vUsageHeader()}
     * @param sampleIndex sample index, as in initial sample collection/set
     * @return
     */
    public double[][] vjUsageMatrix(int sampleIndex) {
        double sampleTotal = (double) (jSegmentUsage.values().collect { it[sampleIndex] }.sum() ?: 0.0)

        def matrix = new double[jSegmentUsage.size()][vSegmentUsage.size()]

        sortedJSegmTotal.eachWithIndex { jEntry, ii ->
            sortedVSegmTotal.eachWithIndex { vEntry, jj ->
                def v = vEntry.key, j = jEntry.key, vj = v + "\t" + j
                matrix[ii][jj] = vjSegmentUsage.containsKey(vj) ? (vjSegmentUsage[vj][sampleIndex] / sampleTotal) : 0d
            }
        }

        matrix
    }

    /**
     * Gets an ordered list of J segment names. V/J/V-J vectors and matrices are ordered correspondingly
     * @return
     */
    public String[] jUsageHeader() {
        sortedJSegmTotal.collect { it.key }
    }

    /**
     * Gets an ordered list of V segment names. V/J/V-J vectors and matrices are ordered correspondingly
     * @return
     */
    public String[] vUsageHeader() {
        sortedVSegmTotal.collect { it.key }
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
