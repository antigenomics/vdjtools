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

import com.antigenomics.vdjtools.Countable
import com.antigenomics.vdjtools.join.ClonotypeKeyGen
import com.antigenomics.vdjtools.join.key.ClonotypeKey
import com.antigenomics.vdjtools.overlap.OverlapType
import com.antigenomics.vdjtools.pool.ClonotypeAggregator
import com.antigenomics.vdjtools.pool.PooledSample
import com.antigenomics.vdjtools.pool.SampleAggregator
import com.antigenomics.vdjtools.pool.StoringClonotypeAggregator
import com.antigenomics.vdjtools.sample.Sample

/**
 * A base class for providing info on the frequencies of rare (singletons, doubletons)
 * and tabulating frequencies of abundant clonotypes. The class is a wrapper for mapping between 
 * clonotype frequency and the number of clonotypes with a given frequency table.
 */
class FrequencyTable {
    private final long count, diversity
    private final Map<Long, FrequencyTableBin> frequencyMap = new HashMap<>()

    /**
     * Creates a frequency table from cache.
     * @param frequencyTableCache a {@code [clonotype count -> number of clonotypes]} mapping.
     */
    FrequencyTable(Map<Long, Long> frequencyTableCache) {
        long count = 0, diversity = 0
        frequencyTableCache.each {
            frequencyMap.put((long) it.key, // absolutely necessary conversion here, time to switch to Scala :)
                    new FrequencyTableBin(it.key, it.value))
            count += it.key * it.value
            diversity += it.value
        }
        this.count = count
        this.diversity = diversity
    }

    /**
     * Gets the frequency table cache.
     * @return {@code [clonotype count -> number of clonotypes]} mapping.
     */
    Map<Long, Long> getCache() {
        frequencyMap.values().collectEntries {
            [(it.count): it.diversity]
        }
    }

    /**
     * Creates frequency table that bins clonotypes according to their frequency.
     * @param sample sample to tabulate.
     */
    FrequencyTable(Sample sample) {
        this(sample, OverlapType.Strict)
    }

    /**
     * Creates frequency table that bins clonotypes according to their frequency.
     * @param sample sample to tabulate.
     * @param intersectionType overlap type used to collapse clonotypes.
     */
    FrequencyTable(Sample sample, OverlapType intersectionType) {

        // collapse clonotypes by a specific key
        def clonotypeKeyGen = new ClonotypeKeyGen(intersectionType)

        // todo: parallelization possible
        def hashedCounts = new HashMap<ClonotypeKey, Counter>()

        sample.each {
            def key = clonotypeKeyGen.generateKey(it)
            def counter = hashedCounts[key]
            if (!counter)
                hashedCounts.put(key, counter = new Counter())
            counter.add(it)
        }

        this.diversity = hashedCounts.size()
        def counters = hashedCounts.values()

        // compute frequency table
        counters.each {
            long x = it.count
            def bin = frequencyMap[x]
            if (!bin)
                frequencyMap.put(x, bin = new FrequencyTableBin(x))
            bin.increment()
        }

        this.count = sample.count
    }

    /**
     * Creates frequency table that bins clonotypes according to their occurrence in a set of samples.
     * @param pool pooled samples.
     */
    FrequencyTable(PooledSample pool) {
        long diversity = 0, count = 0

        // compute frequency table
        pool.each { StoringClonotypeAggregator it ->
            long x = it.incidenceCount
            def bin = frequencyMap[x]
            if (!bin)
                frequencyMap.put(x, bin = new FrequencyTableBin(x))
            bin.increment()
            diversity++
            count += x
        }

        this.count = count
        this.diversity = diversity
    }

    /**
     * Gets the total number of clonotypes in this table, 
     * in accordance with {@link OverlapType} used to collapse the sample.
     * @return total number of clonotypes in the table after they were collapsed.
     */
    long getDiversity() {
        diversity
    }

    /**
     * Gets the total read count in this table.
     * @return total number of reads.
     */
    long getCount() {
        this.count
    }

    /**
     * Gets the number of singleton clonotypes, i.e. those met only once.
     * @return singleton frequency.
     */
    int getSingletons() {
        this[1]
    }

    /**
     * Gets the number of doubleton clonotypes, i.e. those met only twice
     * @return doubleton frequency.
     */
    int getDoubletons() {
        this[2]
    }

    /**
     * Gets the number of clonotypes with the given count.
     * @param count number of reads.
     * @return number of clonotypes with a given read count.
     */
    int getAt(long count) {
        def bin = frequencyMap[count]
        bin ? bin.diversity : 0
    }

    /**
     * Gets the bins that have at least one clonotype in them.
     * @return a collection of non-empty bins.
     */
    Collection<FrequencyTableBin> getBins() {
        Collections.unmodifiableCollection(frequencyMap.values().sort { -it.count })
    }

    /**
     * Header for plain-text output. 
     */
    static final String HEADER = "count_bin\tclonotypes_in_bin"

    @Override
    String toString() {
        bins.collect() {
            it.count + "\t" + it.diversity
        }.join("\n")
    }

    /**
     * Frequency table bin.
     */
    class FrequencyTableBin {
        private final long count
        private int diversity = 0

        /**
         * Creates a new frequency table bin.
         * @param count number of reads, used as key.
         */
        FrequencyTableBin(long count) {
            this.count = count
        }

        /**
         * Creates a new frequency table bin.
         * @param count number of reads, used as key.
         * @param diversity number of clonotypes in this bin.
         */
        FrequencyTableBin(long count, long diversity) {
            this.count = count
            this.diversity = diversity
        }

        /**
         * Increment the number of clonotypes in this bin.
         */
        void increment() {
            this.diversity++
        }

        /**
         * Decrement the number of clonotypes in this bin.
         */
        void decrement() {
            this.diversity--
        }

        /**
         * Get the number of reads that specify this bin.
         * @return read count for this bin.

         public long getCount() {
         this.count
         }/**
         * Gets the number of clonotypes in this bin.
         * @return clonotype diversity in this bin.
         * @see com.antigenomics.vdjtools.diversity.FrequencyTable#getDiversity()
         */
        int getDiversity() {
            this.diversity
        }

        /**
         * Gets the ratio of clonotypes in this bin and the total diversity of the sample.
         * @return ratio of clonotypes in this bin.
         * @see com.antigenomics.vdjtools.diversity.FrequencyTable#getDiversity()
         */
        double getRelativeDiversity() {
            this.diversity / (double) FrequencyTable.this.diversity
        }

        /**
         * Get the clonotype frequency that specify this bin.
         * @return frequency in this bin.
         */
        double getFreq() {
            this.count / (double) FrequencyTable.this.count
        }

        @Override
        boolean equals(o) {
            if (this.is(o)) return true
            if (getClass() != o.class) return false

            FrequencyTableBin that = (FrequencyTableBin) o

            this.count == that.count
        }

        @Override
        int hashCode() {
            (int) (this.count ^ (this.count >>> 32))
        }
    }

    private static class Counter {
        int count = 0

        void add(Countable other) {
            count += other.count
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
