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

import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.sample.SamplePair
import com.antigenomics.vdjtools.sample.metadata.MetadataTable
import com.antigenomics.vdjtools.misc.ExecUtil
import groovyx.gpars.GParsPool

import java.util.concurrent.atomic.AtomicInteger

/**
 * A class that implements all-vs-all paired overlap for a sample collection and
 * holds information on those intersections in a manner it could be easily accessed
 */
public class PairwiseOverlap {
    private final SampleCollection sampleCollection
    private final OverlapType intersectionType
    private final Collection<OverlapMetric> intersectMetrics
    private final Overlap[][] pairedIntersectionCache
    private final int numberOfSamples

    /**
     * Intersects clonotype lists for all unique pairs of samples in a given sample collection.
     * Will load both samples into memory for the initialization step. 
     * Pre-computes all overlap metrics.
     * @param sampleCollection a list of samples
     * @param intersectionType clonotype matching rule
     */
    public PairwiseOverlap(SampleCollection sampleCollection,
                                   OverlapType intersectionType) {
        this(sampleCollection, intersectionType, false, false)
    }

    /**
     * Intersects clonotype lists for all unique pairs of samples in a given sample collection.
     * Pre-computes all overlap metrics.
     * @param sampleCollection a list of samples
     * @param intersectionType clonotype matching rule
     * @param store holds all samples in memory if set to {@code true}
     * @param lowMem if set to {@code true}, will not load all samples in memory, but rather load a sample pair at each step
     */
    public PairwiseOverlap(SampleCollection sampleCollection,
                                   OverlapType intersectionType,
                                   boolean store, boolean lowMem) {
        this(sampleCollection, intersectionType, store, lowMem, OverlapMetric.values())
    }

    /**
     * Intersects clonotype lists for all unique pairs of samples in a given sample collection
     * @param sampleCollection a list of samples
     * @param intersectionType clonotype matching rule
     * @param store holds all samples in memory if set to {@code true}
     * @param lowMem if set to {@code true}, will not load all samples in memory, but rather load a sample pair at each step
     * @param intersectMetrics a list of overlap metrics that should be pre-computed
     */
    public PairwiseOverlap(SampleCollection sampleCollection,
                                   OverlapType intersectionType,
                                   boolean store, boolean lowMem,
                                   Collection<OverlapMetric> intersectMetrics) {
        if (store && lowMem)
            throw new Exception("Isn't it illogical to use 'store' and 'lowMem' options simultaneously?")

        this.sampleCollection = sampleCollection
        this.intersectionType = intersectionType
        this.intersectMetrics = intersectMetrics
        this.numberOfSamples = sampleCollection.size()
        this.pairedIntersectionCache = new Overlap[numberOfSamples][numberOfSamples]

        int totalPairs = numberOfSamples * (numberOfSamples - 1) / 2
        def progressCounter = new AtomicInteger()

        def intersect = { SamplePair pair ->
            pairedIntersectionCache[pair.i][pair.j] =
                    new Overlap(pair, intersectionType, store, intersectMetrics)
            int progr
            if ((progr = progressCounter.incrementAndGet()) % 10 == 0) {
                ExecUtil.report(this, "Processed $progr of $totalPairs pairs. " + ExecUtil.memoryFootprint())
            }
        }

        ExecUtil.report(this, "Started batch overlap for $numberOfSamples samples ($totalPairs pairs)")

        if (lowMem) {
            for (int i = 0; i < numberOfSamples - 1; i++) {
                def pairs = sampleCollection.listPairs(i)
                pairs.each(intersect)
            }
        } else {
            def pairs = sampleCollection.listPairs()

            GParsPool.withPool ExecUtil.THREADS, {
                pairs.eachParallel(intersect)
            }
        }
    }

    /**
     * Gets a paired overlap for a given pair of samples
     * @param i first sample index
     * @param j second sample index
     * @return {@code PairedIntersection} for samples ordered as {@code [i , j]}. Will return {@code null} if {@code i == j}
     * @throws {@code IndexOutOfBoundsException} if {@code i} or {@code j} are not in {@code [0 , numberOfSamples)}
     */
    public Overlap getAt(int i, int j) {
        if (i == j)
            return null // todo: re-implement with dummy batch overlap

        boolean reverse = i > j
        if (reverse)
            (i, j) = [j, i]

        if (i >= numberOfSamples || j < 0)
            throw new IndexOutOfBoundsException()

        reverse ? pairedIntersectionCache[i][j] : pairedIntersectionCache[i][j].reverse
    }

    /**
     * Header string, used for tabular output
     */
    public String getHeader() {
        ["1_$MetadataTable.SAMPLE_ID_COLUMN", "2_$MetadataTable.SAMPLE_ID_COLUMN",
         Overlap.OUTPUT_FIELDS.collect(), intersectMetrics.collect { it.shortName },
         sampleCollection.metadataTable.columnHeader1,
         sampleCollection.metadataTable.columnHeader2].flatten().join("\t")
    }

    /**
     * Plain text row for tabular output
     */
    @Override
    public String toString() {
        (0..<(numberOfSamples - 1)).collect { int i ->
            ((i + 1)..<numberOfSamples).collect { int j ->
                pairedIntersectionCache[i][j].toString()
            }
        }.flatten().join("\n")
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
