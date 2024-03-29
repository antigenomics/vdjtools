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
public final class StatUtil {
    private StatUtil() {

    }

    /// POISSON

    // P(X==k|lambda)
    public static double poissonPdf(int k, double lambda) {
        double ret;
        if (k < 0 || k == Integer.MAX_VALUE) {
            return 0;
        } else if (k == 0) {
            ret = -lambda;
        } else {
            ret = -getStirlingError(k) -
                    getDeviancePart(k, lambda) -
                    0.5 * FastMath.log(MathUtils.TWO_PI) - 0.5 * FastMath.log(k);
        }
        return FastMath.exp(ret);
    }

    // P(X<=k|lambda)
    public static double poissonCdf(int k, double lambda) {
        if (k < 0) {
            return 0;
        }
        if (k == Integer.MAX_VALUE) {
            return 1;
        }
        return Gamma.regularizedGammaQ((double) k + 1.0,
                lambda, 1e-12,
                10000000);
    }

    // P(X>=k|lambda)
    public static double poissonPValue(int k, double lambda) {
        return 1.0 - poissonCdf(k - 1, lambda);
    }

    // P'(X>=k|lambda), P' means special treatment of the X=k case (0.5 factor)
    public static double poissonPValueSymm(int k, double lambda) {
        return 1.0 - poissonCdf(k, lambda) + 0.5 * poissonPdf(k, lambda);
    }

    // two-sided P'(X?k|lambda)
    public static double poissonPValue2(int k, double lambda) {
        double p = poissonPValue(k, lambda);
        return Math.min(p, 1.0 - p);
    }

    /// BINOMIAL

    // P(X==k|n,p)
    public static double binomialPdf(int k, int n, double p) {
        return Math.exp(binomialCoefficientLog(n, k) +
                (n - k) * Math.log(1.0 - p) + k * Math.log(p));
    }

    // P(X<=k|n,p)
    public static double binomialCdf(int k, int n, double p) {
        return regularizedBeta(1.0 - p, n - k, k + 1);
    }

    public static double binomialCdf(int k, double n, double p) {
        return regularizedBeta(1.0 - p, n - k, k + 1);
    }

    // P(X>=k|n,p)
    public static double binomialPValue(int k, int n, double p) {
        return 1.0 - binomialCdf(k - 1, n, p);
    }

    public static double binomialPValue(int k, double n, double p) {
        return 1.0 - binomialCdf(k - 1, n, p);
    }

    // P'(X>=k|n,p), P' means special treatment of the X=k case (0.5 factor)
    public static double binomialPValueSymm(int k, int n, double p) {
        return 1.0 - binomialCdf(k, n, p) + 0.5 * binomialPdf(k, n, p);
    }

    // two-sided P'(X?k|n,p)
    public static double binomialPValue2(int k, int n, double p) {
        double pp = binomialPValue(k, n, p);
        return Math.min(pp, 1.0 - pp);
    }

    /// NEGATIVE BINOMIAL

    /**
     * Computes the PDF of the distribution of number of successes in a sequence of
     * iid Bernoulli trials before a specified number of failures occurs.
     *
     * @param k number of successes
     * @param r number of failures until the experiment is stopped
     * @param p success probability in each experiment
     * @return probability
     */
    public static double negativeBinomialPdf(int k, int r, double p) {
        return Math.exp(binomialCoefficientLog(k + r - 1, k) +
                r * Math.log(1.0 - p) + k * Math.log(p));
    }

    /**
     * Computes the CDF of the distribution of number of successes in a sequence of
     * iid Bernoulli trials before a specified number of failures occurs.
     *
     * @param k number of successes
     * @param r number of failures until the experiment is stopped
     * @param p success probability in each experiment
     * @return cumulative
     */
    public static double negativeBinomialCdf(int k, double r, double p) {
        return 1.0 - regularizedBeta(p, k + 1, r);
    }

    public static double betaBinomialPdf(int k, int n, double alpha, double beta) {
        return Math.exp(
                (Beta.logBeta(k + alpha, n - k + beta) - Beta.logBeta(alpha, beta)) +
                        (Gamma.logGamma(n + 1) - Gamma.logGamma(k + 1) - Gamma.logGamma(n - k + 1))
        );
    }

    /// BETA BINOMIAL

    public static double betaBinomialCdf(int k, int n, double alpha, double beta) {
        double sum = 0;

        for (int i = 0; i < k; i++) {
            sum += betaBinomialPdf(i, n, alpha, beta);
        }

        return Math.min(1.0, sum);
    }

    public static double betaBinomialPvalueFast(int k, int n, double alpha, double beta) {
        return betaBinomialPvalueFast(k, n, alpha, beta, 1e-10);
    }

    public static double betaBinomialPvalueFast(int k, int n, double alpha, double beta, double pThreshold) {
        double sum = 1.0 + 0.5 * betaBinomialPdf(k, n, alpha, beta);

        for (int i = 0; i < k; i++) {
            sum -= betaBinomialPdf(i, n, alpha, beta);

            if (sum <= pThreshold) {
                return pThreshold;
            }
        }

        return sum;
    }

    /// NORMAL

    public static double normalCdf(double x, double mean, double sd) {
        final double dev = x - mean;
        if (FastMath.abs(dev) > 40 * sd) {
            return dev < 0 ? 0.0d : 1.0d;
        }
        return 0.5 * Erf.erfc(-dev / (sd * SQRT2));
    }

    /// HYPERGEOMETRIC

    public static double hypgeomCdf(int n12, int n1, int n2, int nSamples) {
        // http://journals.sagepub.com.sci-hub.cc/doi/pdf/10.2466/pms.1998.87.1.51

        int v = Math.max(0, n1 + n2 - nSamples),
                w = Math.min(n1, n2);

        assert n12 >= v && n12 <= w;

        double pPrev = 1, T = pPrev, S = pPrev;

        for (int i = v + 1; i <= w; i++) {
            pPrev *= ((double) (n1 - i + 1) * (n2 - i + 1)) / i / (nSamples - n1 - n2 + i);
            T += pPrev;

            if (i == n12) {
                S = T - 0.5 * pPrev;
            }
        }

        return S / T;
    }

    public static double hypgeomPValue(int n12, int n1, int n2, int nSamples) {
        double p = hypgeomCdf(n12, n1, n2, nSamples);

        return Math.min(p, 1 - p);
    }

    /// misc from Apache


    private static final double SQRT2 = Math.sqrt(2.0);
    /**
     * 1/2 * log(2 &#960;).
     */
    private static final double HALF_LOG_2_PI = 0.5 * FastMath.log(MathUtils.TWO_PI);

    /**
     * exact Stirling expansion error for certain values.
     */
    private static final double[] EXACT_STIRLING_ERRORS = {0.0, /* 0.0 */
            0.1534264097200273452913848, /* 0.5 */
            0.0810614667953272582196702, /* 1.0 */
            0.0548141210519176538961390, /* 1.5 */
            0.0413406959554092940938221, /* 2.0 */
            0.03316287351993628748511048, /* 2.5 */
            0.02767792568499833914878929, /* 3.0 */
            0.02374616365629749597132920, /* 3.5 */
            0.02079067210376509311152277, /* 4.0 */
            0.01848845053267318523077934, /* 4.5 */
            0.01664469118982119216319487, /* 5.0 */
            0.01513497322191737887351255, /* 5.5 */
            0.01387612882307074799874573, /* 6.0 */
            0.01281046524292022692424986, /* 6.5 */
            0.01189670994589177009505572, /* 7.0 */
            0.01110455975820691732662991, /* 7.5 */
            0.010411265261972096497478567, /* 8.0 */
            0.009799416126158803298389475, /* 8.5 */
            0.009255462182712732917728637, /* 9.0 */
            0.008768700134139385462952823, /* 9.5 */
            0.008330563433362871256469318, /* 10.0 */
            0.007934114564314020547248100, /* 10.5 */
            0.007573675487951840794972024, /* 11.0 */
            0.007244554301320383179543912, /* 11.5 */
            0.006942840107209529865664152, /* 12.0 */
            0.006665247032707682442354394, /* 12.5 */
            0.006408994188004207068439631, /* 13.0 */
            0.006171712263039457647532867, /* 13.5 */
            0.005951370112758847735624416, /* 14.0 */
            0.005746216513010115682023589, /* 14.5 */
            0.005554733551962801371038690 /* 15.0 */
    };


    /**
     * Compute the error of Stirling's series at the given value.
     * <p>
     * References:
     * <ol>
     * <li>Eric W. Weisstein. "Stirling's Series." From MathWorld--A Wolfram Web
     * Resource. <a target="_blank"
     * href="http://mathworld.wolfram.com/StirlingsSeries.html">
     * http://mathworld.wolfram.com/StirlingsSeries.html</a></li>
     * </ol>
     * </p>
     *
     * @param z the value.
     * @return the Striling's series error.
     */
    static double getStirlingError(double z) {
        double ret;
        if (z < 15.0) {
            double z2 = 2.0 * z;
            if (FastMath.floor(z2) == z2) {
                ret = EXACT_STIRLING_ERRORS[(int) z2];
            } else {
                ret = Gamma.logGamma(z + 1.0) - (z + 0.5) * FastMath.log(z) +
                        z - HALF_LOG_2_PI;
            }
        } else {
            double z2 = z * z;
            ret = (0.083333333333333333333 -
                    (0.00277777777777777777778 -
                            (0.00079365079365079365079365 -
                                    (0.000595238095238095238095238 -
                                            0.0008417508417508417508417508 /
                                                    z2) / z2) / z2) / z2) / z;
        }
        return ret;
    }

    /**
     * A part of the deviance portion of the saddle point approximation.
     * <p>
     * References:
     * <ol>
     * <li>Catherine Loader (2000). "Fast and Accurate Computation of Binomial
     * Probabilities.". <a target="_blank"
     * href="http://www.herine.net/stat/papers/dbinom.pdf">
     * http://www.herine.net/stat/papers/dbinom.pdf</a></li>
     * </ol>
     * </p>
     *
     * @param x  the x value.
     * @param mu the average.
     * @return a part of the deviance.
     */
    static double getDeviancePart(double x, double mu) {
        double ret;
        if (FastMath.abs(x - mu) < 0.1 * (x + mu)) {
            double d = x - mu;
            double v = d / (x + mu);
            double s1 = v * d;
            double s = Double.NaN;
            double ej = 2.0 * x * v;
            v *= v;
            int j = 1;
            while (s1 != s) {
                s = s1;
                ej *= v;
                s1 = s + ej / ((j * 2) + 1);
                ++j;
            }
            ret = s1;
        } else {
            ret = x * FastMath.log(x / mu) + mu - x;
        }
        return ret;
    }

    /**
     * Compute the logarithm of the PMF for a binomial distribution
     * using the saddle point expansion.
     *
     * @param x the value at which the probability is evaluated.
     * @param n the number of trials.
     * @param p the probability of success.
     * @param q the probability of failure (1 - p).
     * @return log(p(x)).
     */
    static double logBinomialProbability(int x, int n, double p, double q) {
        double ret;
        if (x == 0) {
            if (p < 0.1) {
                ret = -getDeviancePart(n, n * q) - n * p;
            } else {
                ret = n * FastMath.log(q);
            }
        } else if (x == n) {
            if (q < 0.1) {
                ret = -getDeviancePart(n, n * p) - n * q;
            } else {
                ret = n * FastMath.log(p);
            }
        } else {
            ret = getStirlingError(n) - getStirlingError(x) -
                    getStirlingError(n - x) - getDeviancePart(x, n * p) -
                    getDeviancePart(n - x, n * q);
            double f = (MathUtils.TWO_PI * x * (n - x)) / n;
            ret = -0.5 * FastMath.log(f) + ret;
        }
        return ret;
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
