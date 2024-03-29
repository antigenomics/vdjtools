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

package com.antigenomics.vdjtools.join;

import com.antigenomics.vdjtools.ClonotypeWrapperContainer;
import com.antigenomics.vdjtools.join.key.ClonotypeKey;
import com.antigenomics.vdjtools.overlap.OverlapType;
import com.antigenomics.vdjtools.sample.Clonotype;
import com.antigenomics.vdjtools.sample.Sample;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;

import java.util.*;

/**
 * A join between several clonotype tables, performed using a specific clonotype matching rule {@link #getOverlapType}.
 * Note that clonotypes are grouped to joint clonotypes with a representative clonotype and a list of
 * convergent variants for each sample under the same matching rule. Joint clonotypes can be pre-filtered using
 * a specified {@see JoinFilter}.
 */
public class JointSample implements ClonotypeWrapperContainer<JointClonotype> {
    private final Sample[] samples;
    private final double[] transformedCountSum;
    private final double[] intersectionFreq;
    private final double[][] intersectionFreqMatrix;
    private final long[] intersectionCount;
    private final long[][] intersectionCountMatrix;
    private final int[] intersectionDiv;
    private final int[] totalDiv;
    private final int[][] intersectionDivMatrix;
    private final List<JointClonotype> jointClonotypes;
    private final double totalMeanFreq, minMeanFreq;
    private final int numberOfSamples;
    private final long count;
    private final OverlapType overlapType;
    private final boolean reverse;

    private JointSample(Sample[] samples, double[] transformedCountSum,
                        double[] intersectionFreq, double[][] intersectionFreqMatrix,
                        long[] intersectionCount, long[][] intersectionCountMatrix,
                        int[] totalDiv, int[] intersectionDiv, int[][] intersectionDivMatrix,
                        List<JointClonotype> jointClonotypes,
                        double totalMeanFreq, double minMeanFreq,
                        int numberOfSamples, long count,
                        OverlapType overlapType, boolean reverse) {
        this.samples = samples;
        this.transformedCountSum = transformedCountSum;
        this.intersectionFreq = intersectionFreq;
        this.intersectionFreqMatrix = intersectionFreqMatrix;
        this.intersectionCount = intersectionCount;
        this.intersectionCountMatrix = intersectionCountMatrix;
        this.totalDiv = totalDiv;
        this.intersectionDiv = intersectionDiv;
        this.intersectionDivMatrix = intersectionDivMatrix;
        this.jointClonotypes = jointClonotypes;
        this.totalMeanFreq = totalMeanFreq;
        this.minMeanFreq = minMeanFreq;
        this.numberOfSamples = numberOfSamples;
        this.count = count;
        this.overlapType = overlapType;
        this.reverse = reverse;
    }

    /**
     * Joins clonotype tables from several samples together.
     * Only joint clonotypes detected in two or more samples are retained
     *
     * @param overlapType clonotype matching rule.
     * @param samples     a list of samples to join.
     */
    public JointSample(OverlapType overlapType, Sample[] samples) {
        this(overlapType, samples, new OccurrenceJoinFilter());
    }

    /**
     * Joins clonotype tables from several samples together.
     *
     * @param overlapType clonotype matching rule.
     * @param samples     a list of samples to join.
     * @param joinFilter  specified a filter for joint clonotypes, e.g. by number of occurrences.
     */
    public JointSample(OverlapType overlapType, Sample[] samples,
                       JoinFilter joinFilter) {
        this.numberOfSamples = samples.length;
        this.samples = samples;
        this.transformedCountSum = new double[numberOfSamples];
        this.totalDiv = new int[numberOfSamples];
        this.intersectionDiv = new int[numberOfSamples];
        this.intersectionFreq = new double[numberOfSamples];
        this.intersectionFreqMatrix = new double[numberOfSamples][numberOfSamples];
        this.intersectionCount = new long[numberOfSamples];
        this.intersectionCountMatrix = new long[numberOfSamples][numberOfSamples];
        this.intersectionDivMatrix = new int[numberOfSamples][numberOfSamples];
        this.overlapType = overlapType;
        this.reverse = false;

        ClonotypeKeyGen clonotypeKeyGen = new ClonotypeKeyGen(overlapType);

        Map<ClonotypeKey, JointClonotype> clonotypeMap = new HashMap<>();
        int sampleIndex = 0;
        for (Sample sample : samples) {
            for (Clonotype clonotype : sample) {
                ClonotypeKey key = clonotypeKeyGen.generateKey(clonotype);

                JointClonotype jointClonotype = clonotypeMap.get(key);

                if (jointClonotype == null) {
                    clonotypeMap.put(key, jointClonotype = new JointClonotype(this));
                }

                jointClonotype.addVariant(clonotype, sampleIndex);
            }
            sampleIndex++;
        }

        this.jointClonotypes = new ArrayList<>(clonotypeMap.size() / 2);

        double totalMeanFreq = 0, minMeanFreq = 1;
        int count = 0;
        for (JointClonotype jointClonotype : clonotypeMap.values()) {
            for (int i = 0; i < numberOfSamples; i++) {
                if (jointClonotype.present(i)) {
                    totalDiv[i]++;
                }
            }

            if (joinFilter.pass(jointClonotype)) {
                jointClonotypes.add(jointClonotype);

                double meanFreq = jointClonotype.getBaseFreq();
                totalMeanFreq += meanFreq;
                minMeanFreq = Math.min(minMeanFreq, meanFreq);

                for (int i = 0; i < numberOfSamples; i++) {
                    if (jointClonotype.present(i)) {
                        count += jointClonotype.getCount(i);

                        double freq1 = jointClonotype.getFreq(i);
                        int count1 = jointClonotype.getCount(i);

                        intersectionCount[i] += count1;
                        intersectionFreq[i] += freq1;
                        intersectionDiv[i]++;

                        for (int j = i + 1; j < numberOfSamples; j++) {
                            if (jointClonotype.present(j)) {
                                double freq2 = jointClonotype.getFreq(j);
                                int count2 = jointClonotype.getCount(j);

                                intersectionFreqMatrix[i][j] += freq1;
                                intersectionFreqMatrix[j][i] += freq2;
                                intersectionCountMatrix[i][j] += count1;
                                intersectionCountMatrix[j][i] += count2;
                                intersectionDivMatrix[i][j]++;
                            }
                        }
                    }
                }
            }
            for (int i = 0; i < numberOfSamples; i++) {
                transformedCountSum[i] += transformCount(jointClonotype.getCount(i));
            }
        }
        this.totalMeanFreq = totalMeanFreq;
        this.minMeanFreq = minMeanFreq;
        this.count = count;

        Collections.sort(jointClonotypes);
    }

    private static double transformCount(int count) {
        // www.jstor.org/stable/2332343
        double transformFactor = 3 / (double) 8;
        return Math.sqrt(count + transformFactor);
    }

    private double goodnessOfFit(final JointClonotype jointClonotype) {

        // compute expected counts and chi-square
        double G = 0, expected = 0, expectedD = 0;

        for (int i = 0; i < numberOfSamples; i++) {
            expected += transformCount(jointClonotype.getCount(i));
            expectedD += getTransformedCount(i);
        }
        expected /= expectedD;
        for (int i = 0; i < numberOfSamples; i++) {
            double observed = transformCount(jointClonotype.getCount(i));

            G += observed > 0 ? observed * Math.log(observed / expected / getTransformedCount(i)) : 0;
        }
        return 2 * G;
    }

    /**
     * EXPERIMENTAL Computes the probability that the variance of joint clonotype's abundance can be
     * explained by sampling stochastics for each clonotype.
     */
    public void computeAndCorrectSamplingPValues() {
        List<JointClonotype> jointClonotypes = new ArrayList<>(this.jointClonotypes);

        final ChiSquaredDistribution distribution = new ChiSquaredDistribution(numberOfSamples - 1);

        for (JointClonotype jointClonotype : jointClonotypes) {
            double G = goodnessOfFit(jointClonotype);
            jointClonotype.samplingPValue = 1.0d - distribution.cumulativeProbability(G);
        }

        Collections.sort(jointClonotypes, new Comparator<JointClonotype>() {
            @Override
            public int compare(JointClonotype o1, JointClonotype o2) {
                return -Double.compare(o1.getSamplingPValue(), o2.getSamplingPValue());
            }
        });

        // todo: http://www.biomedcentral.com/1471-2105/9/303 ?
        for (int i = 0; i < jointClonotypes.size(); i++) {
            JointClonotype jointClonotype = jointClonotypes.get(i);
            jointClonotype.samplingPValue *= jointClonotypes.size() / (i + 1);
            jointClonotype.samplingPValue = Math.min(1.0, jointClonotype.samplingPValue);
        }
    }

    /**
     * Gets the number of samples that were joined.
     *
     * @return number of samples.
     */
    public int getNumberOfSamples() {
        return numberOfSamples;
    }

    /**
     * Gets the original sample that was joined by index.
     *
     * @param sampleIndex sample index.
     * @return original sample.
     */
    public Sample getSample(int sampleIndex) {
        return samples[getIndex(sampleIndex)];
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getFreq() {
        return 1.0;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getFreqAsInInput() {
        return 1.0;
    }

    /**
     * Total number of joint clonotypes in this joint sample, not counting convergent variants.
     *
     * @return total number of joint conotypes.
     */
    @Override
    public int getDiversity() {
        return jointClonotypes.size();
    }

    /**
     * Gets the total read count in this joint sample.
     *
     * @return read count sum for all joint clonotypes.
     */
    @Override
    public long getCount() {
        return count;
    }

    private double getTransformedCount(int sampleIndex) {
        return transformedCountSum[getIndex(sampleIndex)];
    }

    protected long getCount(int sampleIndex) {
        return getSample(sampleIndex).getCount();
    }

    /**
     * Gets joint clonotype by index.
     *
     * @param index joint clonotype index.
     * @return joint clonotype.
     */
    @Override
    public JointClonotype getAt(int index) {
        JointClonotype jointClonotype = jointClonotypes.get(index);
        return reverse ? jointClonotype.changeParent(this) : jointClonotype;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean isSorted() {
        return true;
    }

    /**
     * Gets the original number of clonotypes in a given sample, not accounting for sub-variants.
     *
     * @param sampleIndex sample index.
     * @return diversity of a given sample, not accounting for sub-variants.
     */
    public int getTotalDiv(int sampleIndex) {
        return totalDiv[getIndex(sampleIndex)];
    }

    /**
     * Gets the number of overlapping clonotypes in a given sample, not accounting for sub-variants.
     *
     * @param sampleIndex sample index.
     * @return number of clonotypes in intersection coming from a given sample, not accounting for sub-variants.
     */
    public int getIntersectionDiv(int sampleIndex) {
        return intersectionDiv[getIndex(sampleIndex)];
    }

    /**
     * Gets the number of intersecting clonotypes for a given pair of samples, not accounting for sub-variants.
     *
     * @param sampleIndex1 index of first sample.
     * @param sampleIndex2 index of second sample.
     * @return number of clonotypes in intersection, not accounting for sub-variants.
     */
    public int getIntersectionDiv(int sampleIndex1, int sampleIndex2) {
        return sampleIndex1 < sampleIndex2 ? intersectionDivMatrix[sampleIndex1][sampleIndex2] :
                intersectionDivMatrix[sampleIndex2][sampleIndex1];
    }

    /**
     * Gets the frequency sum of overlapping clonotypes in a given sample.
     *
     * @param sampleIndex sample index.
     * @return sum of overlapping clonotype frequencies in a given sample.
     */
    public double getIntersectionFreq(int sampleIndex) {
        return intersectionFreq[getIndex(sampleIndex)];
    }

    /**
     * Gets the frequency sum of clonotypes overlapping specifies between samples.
     * The sum is computed based on clonotype frequencies in the first sample.
     *
     * @param sampleIndex1 index of first sample.
     * @param sampleIndex2 index of second sample.
     * @return sum of overlapping clonotype frequncies, according to first sample.
     */
    public double getIntersectionFreq(int sampleIndex1, int sampleIndex2) {
        return intersectionFreqMatrix[getIndex(sampleIndex1)][getIndex(sampleIndex2)];
    }

    /**
     * Gets the read count sum of overlapping clonotypes in a given sample.
     *
     * @param sampleIndex sample index.
     * @return sum of overlapping clonotype read counts in a given sample.
     */
    public long getIntersectionCount(int sampleIndex) {
        return intersectionCount[getIndex(sampleIndex)];
    }

    /**
     * Gets the read count sum of clonotypes overlapping specifies between samples.
     * The sum is computed based on clonotype frequencies in the first sample.
     *
     * @param sampleIndex1 index of first sample.
     * @param sampleIndex2 index of second sample.
     * @return sum of overlapping clonotype read counts, according to first sample.
     */
    public long getIntersectionCount(int sampleIndex1, int sampleIndex2) {
        return intersectionCountMatrix[getIndex(sampleIndex1)][getIndex(sampleIndex2)];
    }

    /**
     * Gets the clonotype matching rule used to construct a given joint sample.
     *
     * @return clonotype matching rule.
     */
    public OverlapType getOverlapType() {
        return overlapType;
    }

    /**
     * INTERNAL used for plain-text table output of joint sample output.
     */
    public double getTotalMeanFreq() {
        return totalMeanFreq;
    }

    /**
     * INTERNAL used for plain-text table output of joint sample output.
     */
    public double calcFreq(double meanFreq) {
        return meanFreq / totalMeanFreq;
    }

    /**
     * INTERNAL used for plain-text table output of joint sample output.
     */
    public int calcCount(double freq) {
        return (int) (freq / minMeanFreq);
    }

    int getIndex(int sampleIndex) {
        return reverse ? (numberOfSamples - sampleIndex - 1) : sampleIndex;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        JointSample that = (JointSample) o;

        return Arrays.equals(samples, that.samples) && reverse == that.reverse;
    }

    @Override
    public int hashCode() {
        return 31 * Arrays.hashCode(samples) + (reverse ? 1231 : 1237);
    }

    @Override
    public Iterator<JointClonotype> iterator() {
        final Iterator<JointClonotype> iterator = jointClonotypes.iterator();

        if (reverse) {
            return new Iterator<JointClonotype>() {
                @Override
                public boolean hasNext() {
                    return iterator.hasNext();
                }

                @Override
                public JointClonotype next() {
                    return iterator.next().changeParent(JointSample.this);
                }

                @Override
                public void remove() {
                    throw new UnsupportedOperationException("remove");
                }
            };
        }

        return iterator;
    }

    /**
     * INTERNAL Reverses the sample order given joint sample.
     *
     * @return
     */
    public JointSample getReverse() throws Exception {
        if (reverse)
            throw new RuntimeException("Already reversed, for the sake of performance multiple reverse is disabled");

        return new JointSample(samples, transformedCountSum,
                intersectionFreq, intersectionFreqMatrix,
                intersectionCount, intersectionCountMatrix,
                totalDiv, intersectionDiv, intersectionDivMatrix,
                jointClonotypes,
                totalMeanFreq, minMeanFreq,
                numberOfSamples, count, overlapType, true);
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
