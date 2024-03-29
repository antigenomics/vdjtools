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

import com.antigenomics.vdjtools.misc.Software
import com.antigenomics.vdjtools.preprocess.DownSampler
import com.antigenomics.vdjtools.misc.ExecUtil
import org.junit.Test

import static com.antigenomics.vdjtools.TestUtil.getResource
import static com.antigenomics.vdjtools.io.SampleStreamConnection.load
import static com.antigenomics.vdjtools.overlap.OverlapMetric.*

class DistanceMeasureTest {
    static void checkMetricValue(OverlapMetric intersectMetric, double val) {
        switch (intersectMetric) {
            case CorrelationS:
            case Correlation:
                assert val >= -1 && val <= 1
                break

            case Frequency:
            case Frequency2:
            case Diversity:
            case MorisitaHorn:
            case Jaccard:
                //case OverlapMetric.ChaoJaccard:
                //case OverlapMetric.ChaoSorensen:
                assert val >= 0 && val <= 1
                break

            case sJSD:
            case vJSD:
            case vj2JSD:
            case vjJSD:
                assert val >= 0
                break
        }

    }

    @Test
    public void test() {
        def software = Software.VDJtools
        def resStream = getResource("samples/${software.toString().toLowerCase()}.txt.gz")
        def sample = load(resStream, software)

        // downsample
        def downSampler = new DownSampler(sample)
        def nResamples = 200
        int smallCount = sample.count * 0.1, largeCount = sample.count * 0.5

        ExecUtil.quiet()

        def largeSampleOverlapWorse = new HashMap(), selfOverlapWorse = new HashMap()

        values().each {
            largeSampleOverlapWorse.put(it, 0)
            selfOverlapWorse.put(it, 0)
        }

        nResamples.times {
            def s1 = downSampler.reSample(smallCount),
                s2 = downSampler.reSample(smallCount)

            def smallIntersection = new Overlap(s1, s2, OverlapType.AminoAcid)

            s1 = downSampler.reSample(largeCount)
            s2 = downSampler.reSample(largeCount)

            def largeIntersection = new Overlap(s1, s2, OverlapType.AminoAcid),
                selfIntersection = new Overlap(s1, s1, OverlapType.AminoAcid)

            values().each {
                def val1 = smallIntersection.getMetricValue(it),
                    val2 = largeIntersection.getMetricValue(it),
                    val3 = selfIntersection.getMetricValue(it)

                // check values
                checkMetricValue(it, val1)
                checkMetricValue(it, val2)
                checkMetricValue(it, val3)

                // Assure that distance has decreased
                if (it.normalization.normalize(val1) < it.normalization.normalize(val2)) {
                    largeSampleOverlapWorse.put(it, largeSampleOverlapWorse[it] + 1.0)
                }

                // self overlap is always the closest one
                if (it.normalization.normalize(val2) < it.normalization.normalize(val3)) {
                    selfOverlapWorse.put(it, selfOverlapWorse[it] + 1.0)
                }
            }
        }

        double failureFreq

        println "Larger sample has worse overlap"
        largeSampleOverlapWorse.each {
            failureFreq = (it.value / nResamples)
            println(it.key.toString() + "\t" + failureFreq)
            // Note that due to the scaling used in diversity measure (overlap size / (sample1 size * sample2 size))
            // when two equally-well intersecting sub-sample pairs are drawn,
            // the pair with the smallest samples is considered the closest one
            // For RepSeq data, the other possible normalization (overlap size / sqrt(sample1 size * sample2 size))
            // will actually bias towards larger samples, due to higher probability of grabbing similar variants
            // It was empirically estimated that pow(sample1 size * sample2 size, 0.8) is the optimal choice
            if (it.key != Diversity & it.key != CorrelationS)
                assert failureFreq <= 0.01
        }

        println "Self overlap is worse than resampled"
        selfOverlapWorse.each {
            failureFreq = (it.value / nResamples)
            println(it.key.toString() + "\t" + failureFreq)
            assert failureFreq <= 0.01
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
