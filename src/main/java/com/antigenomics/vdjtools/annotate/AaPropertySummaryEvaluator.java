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

package com.antigenomics.vdjtools.annotate;

import com.antigenomics.vdjtools.annotate.partitioning.FullCdr3;
import com.antigenomics.vdjtools.annotate.partitioning.SequenceRegion;
import com.antigenomics.vdjtools.preprocess.DownSampler;
import com.antigenomics.vdjtools.sample.Clonotype;
import com.antigenomics.vdjtools.sample.Sample;
import com.milaboratory.core.sequence.AminoAcidSequence;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;

import java.util.*;

public class AaPropertySummaryEvaluator {
    private final AaProperty aaProperty;
    private final SequenceRegion sequenceRegion;
    private final boolean averageByLength;
    private final boolean weightByFrequency;

    public AaPropertySummaryEvaluator(AaProperty aaProperty, SequenceRegion sequenceRegion,
                                      boolean averageByLength, boolean weightByFrequency) {
        this.aaProperty = aaProperty;
        this.sequenceRegion = sequenceRegion;
        this.averageByLength = averageByLength;
        this.weightByFrequency = weightByFrequency;

        if (aaProperty instanceof Cdr3ContactEstimate && !(sequenceRegion instanceof FullCdr3)) {
            throw new IllegalArgumentException("Cdr3ContactEstimate property can only be used " +
                    "with FullCdr3 sequence region.");
        }
    }

    public float compute(Clonotype clonotype) {
        if (!clonotype.isCoding()) {
            throw new IllegalArgumentException("Cannot compute amino acid properties for non-coding clonotypes.");
        }

        AminoAcidSequence aaSeq = sequenceRegion.extractAminoAcid(clonotype);

        float value = 0;

        for (int i = 0; i < aaSeq.size(); i++) {
            value += aaProperty.compute(aaSeq, i);
        }

        if (averageByLength) {
            value /= aaSeq.size();
        }

        return value;
    }

    public AaPropertySummary compute(Sample sample) {
        double[] values;
        float sum = 0;
        int k = 0;

        if (weightByFrequency) {
            if (sample.getCount() > 1000000L) {
                DownSampler ds = new DownSampler(sample);
                sample = ds.reSample(1000000);
            }

            values = new double[(int) sample.getCount()];

            for (Clonotype clonotype : sample) {
                if (clonotype.isCoding()) {
                    float value = compute(clonotype);
                    for (int i = 0; i < clonotype.getCount(); i++) {
                        values[k++] = value;
                        sum += value;
                    }
                }
            }
        } else {
            values = new double[(int) sample.getCount()];

            for (Clonotype clonotype : sample) {
                if (clonotype.isCoding()) {
                    float value = compute(clonotype);
                    values[k++] = value;
                    sum += value;
                }
            }
        }

        Percentile percentile = new Percentile();

        percentile.setData(values, 0, k);

        return new AaPropertySummary(sum / k,
                (float) percentile.evaluate(50),
                (float) percentile.evaluate(25),
                (float) percentile.evaluate(75));
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
