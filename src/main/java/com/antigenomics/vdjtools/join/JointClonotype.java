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

import com.antigenomics.vdjtools.ClonotypeWrapper;
import com.antigenomics.vdjtools.join.key.ClonotypeKey;
import com.antigenomics.vdjtools.misc.MathUtil;
import com.antigenomics.vdjtools.sample.Clonotype;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * A joint clonotype, including all convergent variants (up to clonotype matching rule, {@link com.antigenomics.vdjtools.overlap.OverlapType})
 * and all occurrences of the specified clonotype in joined samples. Convergent variants are hereafter termed "variants" for simplicity.
 */
public class JointClonotype implements Comparable<JointClonotype>, ClonotypeWrapper {
    private final JointSample parent;
    private final List[] variantsBySample;
    private final int[] counts;
    private int peak = -1, occurrences = -1;
    private Clonotype representative = null;
    private double meanFreq = -1;
    protected double samplingPValue = 1;

    private JointClonotype(JointSample parent,
                           List[] variantsBySample, int[] counts,
                           int peak, Clonotype representative, double meanFreq) {
        this.parent = parent;
        this.variantsBySample = variantsBySample;
        this.counts = counts;
        this.peak = peak;
        this.representative = representative;
        this.meanFreq = meanFreq;
    }

    /**
     * Creates an empty joint clonotype.
     *
     * @param parent parent joint sample.
     */
    public JointClonotype(JointSample parent) {
        this.parent = parent;
        this.variantsBySample = new List[parent.getNumberOfSamples()];
        this.counts = new int[parent.getNumberOfSamples()];
    }

    /**
     * Adds a given clonotype variant to the joint sample.
     *
     * @param variant     clonotype variant.
     * @param sampleIndex index of the sample where the specified clonotype variant was detected.
     */
    @SuppressWarnings("unchecked")
    void addVariant(Clonotype variant, int sampleIndex) {
        List<Clonotype> variants = variantsBySample[sampleIndex];
        if (variants == null)
            variantsBySample[sampleIndex] = (variants = new ArrayList<>());
        variants.add(variant);
        counts[sampleIndex] += variant.getCount();
    }

    /**
     * Gets parent joint sample.
     *
     * @return joint sample this joint clonotype belongs to.
     */
    @Override
    public JointSample getParent() {
        return parent;
    }

    public int getPeak() {
        if (peak < 0) {
            int peakCount = -1;
            for (int i = 0; i < counts.length; i++) {
                if (counts[i] > peakCount) {
                    peak = i;
                    peakCount = counts[i];
                }
            }
        }
        return parent.getIndex(peak);
    }

    /**
     * Gets the number of occurrences of this clonotype. The count is not weighted by the number convergent variants.
     *
     * @return number of samples this clonotype was detected in.
     */
    public int getOccurrences() {
        if (occurrences < 0) {
            occurrences = 0;
            for (int i = 0; i < counts.length; i++) {
                if (present(i))
                    occurrences++;
            }
        }
        return occurrences;
    }

    /**
     * EXPERIMENTAL Gets the probability that the variance of joint clonotype abundance can
     * be explained by sampling stochastics.
     *
     * @return p-value for clonotype abundance variance.
     */
    public double getSamplingPValue() {
        return samplingPValue;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public Clonotype getClonotype() {
        if (representative == null) {
            int max = 0;
            getPeak();
            for (Object o : variantsBySample[peak]) {
                Clonotype clonotype = (Clonotype) o;
                if (clonotype.getCount() > max) {
                    max = (int) clonotype.getCount();
                    representative = clonotype;
                }
            }
        }
        return representative;
    }

    /**
     * Gets the number of variants associated with a given clonotype in specified sample.
     *
     * @param sampleIndex sample index.
     * @return number of variants associated with the clonotype.
     */
    public int getNumberOfVariants(int sampleIndex) {
        return variantsBySample[parent.getIndex(sampleIndex)].size();
    }

    /**
     * Gets the number reads associated with the clonotype and all its variants in the specified sample.
     *
     * @param sampleIndex sample index.
     * @return number of reads associated with the clonotypes and all its variants in a given sample.
     */
    public int getCount(int sampleIndex) {
        return counts[parent.getIndex(sampleIndex)];
    }

    /**
     * Gets the clonotype frequency in a given sample, all variants are accounted.
     *
     * @param sampleIndex sample index.
     * @return frequency of clonotype and all its variants in a given sample.
     */
    public double getFreq(int sampleIndex) {
        return getCount(sampleIndex) / (double) parent.getCount(sampleIndex);
    }

    /**
     * Gets clonotype frequency at a given sample relative to the total frequency of overlapping
     * (i.e. all that passed {@link com.antigenomics.vdjtools.join.JoinFilter}) clonotypes in the specified sample.
     *
     * @param sampleIndex sample index.
     * @return frequency of clonotype and all its variants among all overlapping clonotypes in a given sample.
     */
    public double getFreqWithinIntersection(int sampleIndex) {
        return getFreq(sampleIndex) / parent.getIntersectionFreq(sampleIndex);
    }

    /**
     * Gets a representative (geometric mean) frequency of a given joint clonotype.
     *
     * @return non-normalized joint clonotype frequency.
     */
    public double getBaseFreq() {
        if (meanFreq < 0) {
            meanFreq = 1;
            for (int i = 0; i < parent.getNumberOfSamples(); i++) {
                meanFreq *= (getFreq(i) + MathUtil.JITTER);
            }
            meanFreq = Math.pow(meanFreq, 1.0 / (double) parent.getNumberOfSamples());
        }
        return meanFreq;
    }

    /**
     * Gets a representative (geometric mean) frequency of a given joint clonotype.
     * The frequency is normalized so that the total mean frequency will sum to 1.0 in JointSample,
     * thus allowing to use it during sample output.
     *
     * @return normalized joint clonotype frequency.
     */
    @Override
    public double getFreq() {
        return parent.calcFreq(getBaseFreq());
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getFreqAsInInput() {
        return getFreq();
    }

    /**
     * Gets the total number of clonotype variants associated with a given clonotype in all joined samples,
     * accounting for specified clonotype matching rule.
     *
     * @return total number of clonotype variants.
     */
    @Override
    @SuppressWarnings("unchecked")
    public int getDiversity() {
        Set<ClonotypeKey> variantKeys = new HashSet<>();
        ClonotypeKeyGen clonotypeKeyGen = new ClonotypeKeyGen(parent.getOverlapType());
        for (List variantList : variantsBySample) {
            for (Object variant : variantList) {
                variantKeys.add(clonotypeKeyGen.generateKey((Clonotype) variant));
            }
        }
        return variantKeys.size();
    }

    /**
     * Gets a representative count of a given joint clonotype, which is proportional to representative
     * (geometric mean) frequency. The count is scaled so it is equal to 1 for the smallest JointClonotype,
     * thus allowing to use it during sample output.
     *
     * @return normalized joint clonotype count.
     */
    @Override
    public long getCount() {
        return parent.calcCount(getBaseFreq());
    }

    /**
     * Checks whether the clonotype was detected in a given sample (even with at least one read).
     *
     * @param sampleIndex sample index.
     * @return true if the clonotype was detected in a given sample, false otherwise.
     */
    public boolean present(int sampleIndex) {
        return counts[parent.getIndex(sampleIndex)] > 0;
    }

    JointClonotype changeParent(JointSample newParent) {
        return new JointClonotype(newParent, variantsBySample, counts, peak, representative, meanFreq);
    }

    @Override
    public int compareTo(JointClonotype o) {
        return Double.compare(o.getBaseFreq(), this.getBaseFreq());
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
