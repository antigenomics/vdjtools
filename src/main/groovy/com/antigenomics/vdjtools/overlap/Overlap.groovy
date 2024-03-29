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

import com.antigenomics.vdjtools.join.JointSample
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SamplePair
import com.antigenomics.vdjtools.sample.metadata.MetadataTable
import com.antigenomics.vdjtools.misc.ExecUtil

/**
 * A class that performs an overlap between a pair of samples,
 * holds an exhaustive information on the extent of an overlap and
 * computes a set of overlap metrics
 */
public class Overlap {
    public static boolean VERBOSE = true

    private final SamplePair samplePair
    private final JointSample jointSample
    private final OverlapEvaluator intersectionEvaluator
    private final Map<OverlapMetric, Double> intersectMetricCache
    private final Collection<OverlapMetric> intersectMetrics
    private final int div1, div2, div12, div21, count1, count2, count12, count21
    private final double freq1, freq2, freq12, freq21
    private final boolean store
    private final String header1, header2,
                         id1, id2, meta1, meta2

    /**
     * INTERNAL, just sets up all fields 
     * @param samplePair
     * @param jointSample
     * @param intersectionEvaluator
     * @param intersectMetrics
     * @param intersectMetricCache
     * @param div1
     * @param div2
     * @param div12
     * @param div21
     * @param count1
     * @param count2
     * @param count12
     * @param count21
     * @param freq1
     * @param freq2
     * @param freq12
     * @param freq21
     * @param header1
     * @param header2
     * @param id1
     * @param id2
     * @param meta1
     * @param meta2
     * @param store
     */
    private Overlap(SamplePair samplePair, JointSample jointSample, OverlapEvaluator intersectionEvaluator,
                               Collection<OverlapMetric> intersectMetrics, Map<OverlapMetric, Double> intersectMetricCache,
                               int div1, int div2, int div12, int div21,
                               int count1, int count2, int count12, int count21,
                               double freq1, double freq2, double freq12, double freq21,
                               String header1, String header2, String id1, String id2, String meta1, String meta2,
                               boolean store) {
        this.samplePair = samplePair
        this.jointSample = jointSample
        this.intersectionEvaluator = intersectionEvaluator
        this.intersectMetrics = intersectMetrics
        this.intersectMetricCache = intersectMetricCache
        this.div1 = div1
        this.div2 = div2
        this.div12 = div12
        this.div21 = div21
        this.count1 = count1
        this.count2 = count2
        this.count12 = count12
        this.count21 = count21
        this.freq1 = freq1
        this.freq2 = freq2
        this.freq12 = freq12
        this.freq21 = freq21
        this.header1 = header1
        this.header2 = header2
        this.id1 = id1
        this.id2 = id2
        this.meta1 = meta1
        this.meta2 = meta2
        this.store = store
    }

    /**
     * Intersects a pair of samples and stores all results.
     * Will load both samples into memory for the initialization step.
     * Pre-computes all overlap metrics.
     * @param sample1 first sample to be intersected.
     * @param sample2 second sample to be intersected.
     * @param intersectionType clonotype matching rule
     */
    public Overlap(Sample sample1, Sample sample2, OverlapType intersectionType) {
        this(new SamplePair(sample1, sample2), intersectionType, false)
    }

    /**
     * Intersects a pair of samples and stores all results.
     * Will load both samples into memory for the initialization step.
     * Pre-computes all overlap metrics.
     * @param samplePair an object holding samples to be intersected
     * @param intersectionType clonotype matching rule
     */
    public Overlap(SamplePair samplePair, OverlapType intersectionType) {
        this(samplePair, intersectionType, false)
    }

    /**
     * Intersects a pair of samples and stores all results. 
     * Pre-computes all overlap metrics.
     * @param samplePair an object holding samples to be intersected
     * @param intersectionType clonotype matching rule
     * @param store holds all samples in memory if set to {@code true}
     */
    public Overlap(SamplePair samplePair, OverlapType intersectionType, boolean store) {
        this(samplePair, intersectionType, store, OverlapMetric.values())
    }

    /**
     * Intersects a pair of samples and stores all results.
     * @param samplePair an object holding samples to be intersected
     * @param intersectionType clonotype matching rule
     * @param store holds all samples in memory if set to {@code true}
     * @param intersectMetrics a list of overlap metrics that should be pre-computed
     */
    public Overlap(SamplePair samplePair,
                              OverlapType intersectionType,
                              boolean store,
                              Collection<OverlapMetric> intersectMetrics) {
        this.store = store
        this.samplePair = store ? samplePair : null
        ExecUtil.report(this, "Intersecting samples #${samplePair.i} and ${samplePair.j}", VERBOSE)
        def jointSample = new JointSample(intersectionType, [samplePair[0], samplePair[1]] as Sample[])
        this.jointSample = store ? jointSample : null
        def intersectionEvaluator = new OverlapEvaluator(jointSample)
        this.intersectionEvaluator = store ? intersectionEvaluator : null
        this.intersectMetrics = intersectMetrics
        this.intersectMetricCache = new HashMap<>()

        intersectMetrics.each {
            intersectMetricCache.put(it, intersectionEvaluator.computeIntersectionMetric(it))
        }

        this.div1 = samplePair[0].diversity
        this.div2 = samplePair[1].diversity
        this.div12 = jointSample.diversity
        this.div21 = div12
        this.count1 = samplePair[0].count
        this.count2 = samplePair[1].count
        this.count12 = jointSample.getIntersectionCount(0, 1)
        this.count21 = jointSample.getIntersectionCount(1, 0)
        this.freq1 = samplePair[0].freqAsInInput
        this.freq2 = samplePair[1].freqAsInInput
        this.freq12 = jointSample.getIntersectionFreq(0, 1)
        this.freq21 = jointSample.getIntersectionFreq(1, 0)

        this.header1 = samplePair[0].sampleMetadata.parent.columnHeader1
        this.header2 = samplePair[0].sampleMetadata.parent.columnHeader2
        this.id1 = samplePair[0].sampleMetadata.sampleId
        this.id2 = samplePair[1].sampleMetadata.sampleId
        this.meta1 = samplePair[0].sampleMetadata.toString()
        this.meta2 = samplePair[1].sampleMetadata.toString()
    }

    /**
     * Gets the value of a specified overlap metric. Uses cache.
     * @param intersectMetric overlap metric type
     * @return value of overlap metric which can lie in {@code ( - inf , + inf )}
     * @throws Exception if metric is not pre-computed and {@code store=false}
     */
    public double getMetricValue(OverlapMetric intersectMetric) {
        def value = intersectMetricCache[intersectMetric]
        if (value == null) {
            if (!store)
                throw new Exception("Cannot provided value for ${intersectMetric.shortName} as " +
                        "\$store=false and the value is not precomputed")
            else
                intersectMetricCache.put(intersectMetric,
                        value = intersectionEvaluator.computeIntersectionMetric(intersectMetric))
        }
        value
    }

    /**
     * Gets the first sample in overlap
     * @return sample object
     * @throws Exception if {@code store=false}
     */
    public Sample getSample1() {
        if (!store)
            throw new Exception("Cannot access this property as \$store=false")
        samplePair[0]
    }

    /**
     * Gets the second sample in overlap
     * @return sample object
     * @throws Exception if {@code store=false}
     */
    public Sample getSample2() {
        if (!store)
            throw new Exception("Cannot access this property as \$store=false")
        samplePair[1]
    }

    /**
     * Gets the joint sample that contains all shared clonotypes
     * @return joint sample object
     */
    public JointSample getJointSample() {
        if (!store)
            throw new Exception("Cannot access this property as \$store=false")
        jointSample
    }

    /**
     * Gets the number of unique (up to matching rule) clonotypes in first sample
     * @return the diversity of first sample
     */
    public int getDiv1() {
        div1
    }

    /**
     * Gets the number of unique (up to matching rule) clonotypes in second sample
     * @return the diversity of second sample
     */
    public int getDiv2() {
        div2
    }

    /**
     * Gets the number of unique (up to matching rule) clonotypes that overlap between samples
     * Same as {@code getDiv21 ( )}
     * @return sample overlap diversity
     */
    public int getDiv12() {
        div12
    }

    /**
     * Gets the number of unique (up to matching rule) clonotypes that overlap between samples. 
     * Same as {@code getDiv12 ( )}
     * @return sample overlap diversity
     */
    public int getDiv21() {
        div21
    }

    /**
     * Gets the number of reads in the first sample 
     * @return read count
     */
    public int getCount1() {
        count1
    }

    /**
     * Gets the number of reads in the second sample 
     * @return read count
     */
    public int getCount2() {
        count2
    }

    /**
     * Gets the number of reads in the first sample that belong to clonotypes overlapping between samples
     * @return read count of overlapping clonotypes according to first sample
     */
    public int getCount12() {
        count12
    }

    /**
     * Gets the number of reads in the second sample that belong to clonotypes overlapping between samples
     * @return read count of overlapping clonotypes according to second sample
     */
    public int getCount21() {
        count21
    }

    /**
     * Gets the total frequency of the first sample
     * @return frequency , should return {@code 1.0} in most cases
     */
    public double getFreq1() {
        freq1
    }

    /**
     * Gets the total frequency of the second sample 
     * @return frequency , should return {@code 1.0} in most cases
     */
    public double getFreq2() {
        freq2
    }

    /**
     * Gets the frequency of reads in the first sample that belong to clonotypes overlapping between samples
     * @return frequency of overlapping clonotypes according to first sample
     */
    public double getFreq12() {
        freq12
    }

    /**
     * Gets the frequency of reads in the second sample that belong to clonotypes overlapping between samples
     * @return frequency of overlapping clonotypes according to second sample
     */
    public double getFreq21() {
        freq21
    }

    /**
     * List of fields that will be included in tabular output 
     */
    public static final String[] OUTPUT_FIELDS = ["div1", "div2", "div12", "div21",
                                                  "count1", "count2", "count12", "count21",
                                                  "freq1", "freq2", "freq12", "freq21"]

    /**
     * Swaps samples and all fields: {@code div1} is swapped with {@code div2}, etc.
     * Method is mostly used for output as by default only the lower triangular of overlap matrix is stored
     * @return a paired overlap instance for swapped pair of samples
     */
    public Overlap getReverse() {
        new Overlap(store ? samplePair.reverse : null, store ? jointSample.reverse : null,
                intersectionEvaluator, intersectMetrics, intersectMetricCache,
                div2, div1, div21, div12,
                count2, count1, count21, count12,
                freq2, freq1, freq21, freq12,
                header1, header2, id2, id1, meta2, meta1,
                store)
    }

    /**
     * Header string, used for tabular output
     */
    public String getHeader() {
        ["1_$MetadataTable.SAMPLE_ID_COLUMN", "2_$MetadataTable.SAMPLE_ID_COLUMN",
         OUTPUT_FIELDS.collect(), intersectMetrics.collect { it.shortName },
         header1, header2].flatten().join("\t")
    }

    /**
     * Generic header string
     */
    public static String HEADER =
            ["1_$MetadataTable.GENERIC_METADATA_TABLE.SAMPLE_ID_COLUMN",
             "2_$MetadataTable.GENERIC_METADATA_TABLE.SAMPLE_ID_COLUMN",
             OUTPUT_FIELDS.collect(), OverlapMetric.collect { it.shortName },
             MetadataTable.GENERIC_METADATA_TABLE.columnHeader1,
             MetadataTable.GENERIC_METADATA_TABLE.columnHeader2].flatten().join("\t")

    /**
     * Plain text row for tabular output
     */
    @Override
    public String toString() {
        [id1, id2,
         OUTPUT_FIELDS.collect { this."$it" }, intersectMetrics.collect { intersectMetricCache[it] },
         meta1, meta2].flatten().join("\t")
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
