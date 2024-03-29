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

package com.antigenomics.vdjtools.pool;

import com.antigenomics.vdjtools.sample.Clonotype;
import com.antigenomics.vdjtools.sample.ClonotypeFilter;
import com.antigenomics.vdjtools.sample.Sample;

/**
 * A clonotype filter that filters out all clonotypes that do not pass a certain ratio threshold
 * when compared to the matching most abundant clonotype in another sample.
 * Used in {@link com.antigenomics.vdjtools.preprocess.Decontaminate}.
 */
public class RatioFilter extends ClonotypeFilter {
    private final SampleAggregator<MaxClonotypeAggregator> sampleAggregator;
    private final double thresholdRatio;

    public RatioFilter(Iterable<Sample> samples, double thresholdRatio, boolean negative) {
        super(negative);
        this.sampleAggregator = new SampleAggregator<>(samples, new MaxClonotypeAggregatorFactory());
        this.thresholdRatio = thresholdRatio;
    }

    public RatioFilter(Iterable<Sample> samples, double thresholdRatio) {
        this(samples, thresholdRatio, false);
    }

    public RatioFilter(Iterable<Sample> samples) {
        this(samples, 20.0);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected boolean checkPass(Clonotype clonotype) {
        MaxClonotypeAggregator aggregator = sampleAggregator.getAt(clonotype);
        return aggregator == null ||
                aggregator.getMaxFreq() < clonotype.getFreq() * thresholdRatio;
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
