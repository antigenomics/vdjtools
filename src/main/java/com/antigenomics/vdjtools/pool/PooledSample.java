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

import com.antigenomics.vdjtools.ClonotypeWrapperContainer;
import com.antigenomics.vdjtools.overlap.OverlapType;
import com.antigenomics.vdjtools.sample.SampleCollection;

import java.util.*;

public class PooledSample implements ClonotypeWrapperContainer<StoringClonotypeAggregator> {
    private final List<StoringClonotypeAggregator> clonotypes;
    private final long count;

    @SuppressWarnings("unchecked")
    public PooledSample(SampleCollection samples) {
        this(new SampleAggregator(samples,
                new StoringClonotypeAggregatorFactory(), OverlapType.Strict));
    }

    public PooledSample(SampleAggregator<StoringClonotypeAggregator> sampleAggregator) {
        this.clonotypes = new ArrayList<>(sampleAggregator.getDiversity());

        long count = 0;

        for (StoringClonotypeAggregator clonotypeAggregator : sampleAggregator) {
            clonotypeAggregator.setParent(this);
            int x = (int) clonotypeAggregator.getCount();
            count += x;
            clonotypes.add(clonotypeAggregator);
        }

        this.count = count;

        Collections.sort(clonotypes,
                (o1, o2) -> {
                    return Long.compare(o2.getCount(), o1.getCount()); // inverse - sort descending
                });
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
     * {@inheritDoc}
     */
    @Override
    public long getCount() {
        return count;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getDiversity() {
        return clonotypes.size();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public StoringClonotypeAggregator getAt(int index) {
        return clonotypes.get(index);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean isSorted() {
        return true;
    }

    @Override
    public Iterator<StoringClonotypeAggregator> iterator() {
        return clonotypes.iterator();
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
