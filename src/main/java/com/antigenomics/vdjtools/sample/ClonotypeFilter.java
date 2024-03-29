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

package com.antigenomics.vdjtools.sample;

import com.antigenomics.vdjtools.misc.AtomicDouble;

import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;

/**
 * A base class for clonotype filtering rule. This class also provides summary statistics.
 */
public abstract class ClonotypeFilter {
    private final AtomicInteger passedClonotypes = new AtomicInteger(),
            totalClonotypes = new AtomicInteger();
    private final AtomicLong passedCount = new AtomicLong(),
            totalCount = new AtomicLong();
    private final AtomicDouble passedFreq = new AtomicDouble(),
            totalFreq = new AtomicDouble();
    private final boolean negative;

    protected ClonotypeFilter(boolean negative) {
        this.negative = negative;
    }

    protected ClonotypeFilter() {
        this(false);
    }

    /**
     * Checks whether a given clonotype passes the filter.
     *
     * @param clonotype a clonotype.
     * @return clonotype passes the filter and should be retained.
     */
    public boolean pass(Clonotype clonotype) {
        boolean pass = negative ^ checkPass(clonotype);

        if (pass) {
            passedClonotypes.incrementAndGet();
            passedCount.addAndGet(clonotype.getCount());
            passedFreq.addAndGet(clonotype.getFreq());
        }

        totalClonotypes.incrementAndGet();
        totalCount.addAndGet(clonotype.getCount());
        totalFreq.addAndGet(clonotype.getFreq());

        return pass;
    }

    protected abstract boolean checkPass(Clonotype clonotype);

    /**
     * Gets filtering statistics.
     *
     * @return clonotype filtering statistics.
     */
    public ClonotypeFilterStats getStats() {
        return new ClonotypeFilterStats(passedClonotypes.get(), totalClonotypes.get(),
                passedCount.get(), totalCount.get(),
                passedFreq.get(), totalFreq.get());
    }

    /**
     * Gets filtering statistics and resets counters.
     *
     * @return clonotype filtering statistics.
     */
    public ClonotypeFilterStats getStatsAndFlush() {
        ClonotypeFilterStats stats = getStats();
        passedClonotypes.set(0);
        totalClonotypes.set(0);
        passedCount.set(0L);
        totalCount.set(0L);
        passedFreq.set(0.0);
        totalFreq.set(0.0);
        return stats;
    }

    /**
     * Clonotype filtering statistics.
     */
    public static class ClonotypeFilterStats {
        private final int passedClonotypes, totalClonotypes;
        private final long passedCount, totalCount;
        private final double passedFreq, totalFreq;

        public ClonotypeFilterStats(int passedClonotypes, int totalClonotypes,
                                    long passedCount, long totalCount,
                                    double passedFreq, double totalFreq) {
            this.passedClonotypes = passedClonotypes;
            this.totalClonotypes = totalClonotypes;
            this.passedCount = passedCount;
            this.totalCount = totalCount;
            this.passedFreq = passedFreq;
            this.totalFreq = totalFreq;
        }

        public int getPassedClonotypes() {
            return passedClonotypes;
        }

        public int getTotalClonotypes() {
            return totalClonotypes;
        }

        public long getPassedCount() {
            return passedCount;
        }

        public long getTotalCount() {
            return totalCount;
        }

        public double getPassedFreq() {
            return passedFreq;
        }

        public double getTotalFreq() {
            return totalFreq;
        }

        public static final String HEADER = "passed_clones\ttotal_clones\t" +
                "passed_count\ttotal_count\t" +
                "passed_freq\ttotal_freq";

        @Override
        public String toString() {
            return passedClonotypes + "\t" + totalClonotypes + "\t" +
                    passedCount + "\t" + totalCount + "\t" +
                    passedFreq + "\t" + totalFreq;
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
