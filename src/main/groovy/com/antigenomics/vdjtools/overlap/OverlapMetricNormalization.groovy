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

/**
 * Normalization type that is recommended for a given IntersectMetric,
 * should transform overlap metric value to {@code [0 , +inf)} scale.
 */
public enum OverlapMetricNormalization {
    /**
     * Negative logarithm normalization {@code -log10(x + 1e-9)}.
     */
    /*   */ NegLog(0),
    /**
     * Correlation metrics, normalized as {@code ( 1 - x ) /2}.
     */
            R(1),
    /**
     * Simialrity index, normalized as {@code ( 1 - x )}.
     */
            Index(2),
    /**
     * Metrics for which normalization is not required.
     */
            None(-1)

    /**
     * Normalization type ID, used for passing to R scripts as argument.
     */
    public final int id

    private OverlapMetricNormalization(int id) {
        this.id = id
    }

    /**
     * Normalize raw overlap metric accordingly.
     * @param x overlap metric value.
     * @return normalized value.
     */
    public double normalize(double x) {
        switch (this) {
            case NegLog:
                return -Math.log10(x + 1e-9)
            case R:
                return (1.0 - x) / 2
            case Index:
                return 1.0 - x
        }
        x
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
