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

import com.antigenomics.vdjtools.overlap.OverlapType
import com.antigenomics.vdjtools.sample.Sample

/**
 * A class that is used to compute a rarefaction curve based on {@link com.antigenomics.vdjtools.diversity.ChaoEstimator}.
 * Rarefaction curve shows the relation between sample richness and sample size.
 * Implementation is based on Colwell et al 2012 paper.
 *
 * @url http://viceroy.eeb.uconn.edu/estimates/EstimateSPages/EstSUsersGuide/References/ColwellEtAl2012.pdf
 */
class Rarefaction {
    private final FrequencyTable frequencyTable
    private final ChaoEstimator chaoEstimator
    private final long n

    /**
     * Computes the interpolated/exact/extrapolated richness estimate for a given sample size.
     * @param coord sample size.
     * @return species richness estimate.
     */
    SpeciesRichness getAt(long coord) {
        coord > n ? chaoEstimator.chaoE(coord) : chaoEstimator.chaoI(coord)
    }

    /**
     * Creates an instance that will perform rarefaction for a given sample.
     * @param sample sample to analyze.
     * @param intersectionType specifies a method that will be used to collapse clonotypes when computing a {@link com.antigenomics.vdjtools.diversity.FrequencyTable}.
     */
    Rarefaction(Sample sample, OverlapType intersectionType) {
        this.frequencyTable = new FrequencyTable(sample, intersectionType)
        this.n = frequencyTable.count
        this.chaoEstimator = new ChaoEstimator(frequencyTable)
    }

    /**
     * Creates an instance that will perform rarefaction for a given frequency table.
     * @param frequencyTable a {@link com.antigenomics.vdjtools.diversity.FrequencyTable} summarizing the number of clonotypes met once, twice, ... in a given clonotype set.
     */
    Rarefaction(FrequencyTable frequencyTable) {
        this.frequencyTable = frequencyTable
        this.n = frequencyTable.count
        this.chaoEstimator = new ChaoEstimator(frequencyTable)
    }

    /**
     * Creates an instance that will perform rarefaction for a given sample. 
     * Will use {@code com.antigenomics.vdjtools.overlap.IntersectionType # Strict} to build a {@link com.antigenomics.vdjtools.diversity.FrequencyTable}.
     * @param sample sample to analyze.
     */
    Rarefaction(Sample sample) {
        this(sample, OverlapType.Strict)
    }

    /**
     * Builds a rarefaction curve interpolation up to given sample size.
     * @return create the interpolated part of rarefaction curve.
     */
    ArrayList<RarefactionPoint> interpolate() {
        build(n)
    }

    /**
     * Extrapolates rarefaction curve from the size of a given sample up to a specified size.
     * @param to where to extrapolate.
     * @return the extrapolated part of rarefaction curve.
     * @throws IllegalArgumentException if {@code to} is less than {@code sample size + 1}.
     */
    ArrayList<RarefactionPoint> extrapolate(long to) {
        build(n + 1, to)
    }

    /**
     * Builds a rarefaction curve interpolation up to given sample size.
     * @param numberOfPoints number of size steps.
     * @return create the interpolated part of rarefaction curve.
     */
    ArrayList<RarefactionPoint> interpolate(int numberOfPoints) {
        build(n, numberOfPoints)
    }

    /**
     * Extrapolates rarefaction curve starting from the size of a given sample up to a specified size.
     * @param to where to extrapolate.
     * @param numberOfPoints number of points in rarefaction curve.
     * @return
     * @throws IllegalArgumentException if {@code to} is less than {@code sample size + 1}.
     */
    ArrayList<RarefactionPoint> extrapolate(long to, int numberOfPoints) {
        build(n + 1, to, numberOfPoints)
    }

    /**
     * Builds a rarefaction curve from sample size {@code 0} up to sample size {@code to}.
     * @param to sample size for last rarefaction point.
     * @return rarefaction curve.
     */
    ArrayList<RarefactionPoint> build(long to) {
        build(0L, to)
    }

    /**
     * Builds a rarefaction curve starting from sample size {@code from} up to sample size {@code to}.
     * @param from sample size for first rarefaction point.
     * @param to sample size for last rarefaction point.
     * @return rarefaction curve.
     */
    ArrayList<RarefactionPoint> build(long from, long to) {
        build(from, to, Math.min(101, (int) n))
    }

    /**
     * Builds a rarefaction curve starting from sample size {@code from} up to sample size {@code to}.
     * @param from sample size for first rarefaction point.
     * @param to sample size for last rarefaction point.
     * @param numberOfPoints number of points in rarefaction curve.
     * @return rarefaction curve.
     */
    ArrayList<RarefactionPoint> build(long from, long to, int numberOfPoints) {
        if (from > to)
            throw new IllegalArgumentException("from should be less than to")
        def rarefactionCurve = new ArrayList<RarefactionPoint>()

        double step = (to - from) / (double) (numberOfPoints - 1)
        boolean hasExact = false

        for (int i = 0; i < numberOfPoints - 1; i++) {
            long m = from + i * step

            if (m == n)
                hasExact = true

            rarefactionCurve.add(new RarefactionPoint(this[m]))
        }

        if (!hasExact)
            rarefactionCurve.add(new RarefactionPoint(chaoEstimator.chaoI(n)))

        // add the last point (more robust & label placing in plotting)
        rarefactionCurve.add(new RarefactionPoint(chaoEstimator.chaoE(to)))

        rarefactionCurve.sort { it.x }
    }

    /**
     * Holds summary for rarefied richness estimate.
     */
    static class RarefactionPoint {
        private final double x, mean, ciU, ciL
        private final RichnessEstimateType richnessType

        public RarefactionPoint(SpeciesRichness richness) {
            this.richnessType = richness.type
            this.x = richnessType == RichnessEstimateType.TotalDiversityLowerBoundEstimate ? 1e20 : richness.size
            this.mean = richness.mean
            this.ciL = mean - 1.96 * richness.std
            this.ciU = mean + 1.96 * richness.std
        }

        /**
         * Gets the coordinate, i.e. sample size.
         * @return {@code x} coordinate of rarefaction curve.
         */
        double getX() {
            x
        }

        /**
         * Gets the mean of rarefied richness estimate.
         * @return {@code y} coordinate of rarefaction curve.
         */
        double getMean() {
            mean
        }

        /**
         * Gets the upper bound of 95% confidence interval of rarefied richness estimate.
         * @return upper bound coordinate of rarefaction curve.
         */
        double getCiU() {
            ciU
        }

        /**
         * Gets the lower bound of 95% confidence interval of rarefied richness estimate.
         * @return lower bound coordinate of rarefaction curve.
         */
        double getCiL() {
            ciL
        }

        /**
         * Gets the richness estimate type: interpolated, exact or extrapolated.
         * @return richness type.
         */
        RichnessEstimateType getRichnessType() {
            richnessType
        }

        /**
         * Header string, used for tabular output.
         */
        static final String HEADER = "x\tmean\tciL\tciU\ttype"

        /**
         * Plain text row for tabular output.
         */
        @Override
        public String toString() {
            [x, mean, ciL, ciU, richnessType.id].join("\t")
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
