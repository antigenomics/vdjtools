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

package com.antigenomics.vdjtools.preprocess;

import com.antigenomics.vdjtools.misc.Segment;
import com.antigenomics.vdjtools.sample.Clonotype;
import com.antigenomics.vdjtools.sample.Sample;
import com.milaboratory.core.sequence.NucleotideSequence;
import com.milaboratory.core.tree.NeighborhoodIterator;
import com.milaboratory.core.tree.SequenceTreeMap;
import com.milaboratory.core.tree.TreeSearchParameters;

import java.util.HashMap;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * A class implementing simple frequency-based error correction.
 */
public class Corrector {
    private final TreeSearchParameters treeSearchParameters;
    private final double logRatioThreshold;
    private final boolean requireSegmentMatch;

    /**
     * Creates a sample corrector. Clonotypes that are no more than 2 mismatches difference from each other,
     * have a more than 20 times ratio threshold will be subject to correction.
     */
    public Corrector() {
        this(2, 0.05f, false);
    }

    /**
     * Creates a sample corrector. Clonotypes that are in a certain edit distance from each other,
     * have a certain ratio threshold and matching segments (if specified) will be subject to correction.
     *
     * @param maxMismatches       maximum number of mismatches allowed.
     * @param ratioThreshold      clonotype abundance ratio.
     * @param requireSegmentMatch requires matching of V/J segment.
     */
    public Corrector(int maxMismatches, float ratioThreshold, boolean requireSegmentMatch) {
        this.treeSearchParameters = new TreeSearchParameters(maxMismatches, 0, 0);
        this.logRatioThreshold = -Math.log10(ratioThreshold);
        this.requireSegmentMatch = requireSegmentMatch;
    }

    /**
     * Corrects specified sample for clonotypes that result from PCR/sequencing errors.
     *
     * @param sample a sample.
     * @return corrected sample, a deep copy.
     */
    public Sample correct(Sample sample) {
        final SequenceTreeMapProvider sequenceTreeMapProvider =
                requireSegmentMatch ? new SegmentSequenceTreeMapProvider(sample) :
                        new SimpleSequenceTreeMapProvider(sample);

        Map<Clonotype, Integer> samplerMap =
                sample.getClonotypes().parallelStream().collect(
                        Collectors.toMap(
                                it -> it,
                                it -> computeCorrectedCount(it, sequenceTreeMapProvider.get(it)),
                                (count1, count2) -> count1 + count2
                        )
                );

        return new Sample(sample, samplerMap);
    }

    private int computeCorrectedCount(Clonotype it, SequenceTreeMap<NucleotideSequence, Clonotype> sequenceTreeMap) {
        NeighborhoodIterator<NucleotideSequence, Clonotype> neighborhoodIterator =
                sequenceTreeMap.getNeighborhoodIterator(it.getCdr3ntBinary(), treeSearchParameters);

        int totalCount = (int) it.getCount();

        for (Clonotype other : neighborhoodIterator.it()) {
            if (!other.equals(it)) {
                int mismatches = neighborhoodIterator.getMismatches();
                double logRatio = Math.log10(it.getCount() / (double) other.getCount());

                if (logRatio > mismatches * logRatioThreshold) {
                    totalCount += other.getCount();
                } else if (logRatio < -mismatches * logRatioThreshold) {

                    return 0;
                }
            }
        }

        return totalCount;
    }

    public TreeSearchParameters getTreeSearchParameters() {
        return treeSearchParameters;
    }

    public double getRatioThreshold() {
        return Math.pow(10, -logRatioThreshold);
    }

    private static interface SequenceTreeMapProvider {
        public SequenceTreeMap<NucleotideSequence, Clonotype> get(Clonotype clonotype);
    }

    private static class SegmentSequenceTreeMapProvider implements SequenceTreeMapProvider {
        private final Map<VJBin, SequenceTreeMap<NucleotideSequence, Clonotype>> sequenceTreeMapBySegment =
                new HashMap<>();

        public SegmentSequenceTreeMapProvider(Sample sample) {
            for (int i = sample.getDiversity() - 1; i >= 0; i--) {
                Clonotype clonotype = sample.getAt(i);
                VJBin vjBin = new VJBin(clonotype);
                SequenceTreeMap<NucleotideSequence, Clonotype> sequenceTreeMap = sequenceTreeMapBySegment.get(vjBin);

                if (sequenceTreeMap == null) {
                    sequenceTreeMapBySegment.put(vjBin,
                            sequenceTreeMap = new SequenceTreeMap<>(NucleotideSequence.ALPHABET));
                }

                sequenceTreeMap.put(clonotype.getCdr3ntBinary(), clonotype);
            }
        }

        public SequenceTreeMap<NucleotideSequence, Clonotype> get(Clonotype clonotype) {
            return sequenceTreeMapBySegment.get(new VJBin(clonotype));
        }

        private static class VJBin {
            private final Segment v, j;

            public VJBin(Clonotype clonotype) {
                this(clonotype.getVBinary(), clonotype.getJBinary());
            }

            public VJBin(Segment v, Segment j) {
                this.v = v;
                this.j = j;
            }

            @Override
            public boolean equals(Object o) {
                if (this == o) return true;
                if (o == null || getClass() != o.getClass()) return false;

                VJBin vjBin = (VJBin) o;

                return j.equals(vjBin.j) && v.equals(vjBin.v);
            }

            @Override
            public int hashCode() {
                return 31 * v.hashCode() + j.hashCode();
            }
        }
    }

    private static class SimpleSequenceTreeMapProvider implements SequenceTreeMapProvider {
        private final SequenceTreeMap<NucleotideSequence, Clonotype> sequenceTreeMap =
                new SequenceTreeMap<>(NucleotideSequence.ALPHABET);

        public SimpleSequenceTreeMapProvider(Sample sample) {
            for (int i = sample.getDiversity() - 1; i >= 0; i--) {
                Clonotype clonotype = sample.getAt(i);

                sequenceTreeMap.put(clonotype.getCdr3ntBinary(), clonotype);
            }
        }

        public SequenceTreeMap<NucleotideSequence, Clonotype> get(Clonotype clonotype) {
            return sequenceTreeMap;
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
