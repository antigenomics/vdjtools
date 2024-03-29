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

import com.antigenomics.vdjtools.join.ClonotypeKeyGen;
import com.antigenomics.vdjtools.join.key.ClonotypeKey;
import com.antigenomics.vdjtools.overlap.OverlapType;
import com.antigenomics.vdjtools.sample.Clonotype;
import com.antigenomics.vdjtools.sample.Sample;

import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 * A sample aggregator used in constructing {@link com.antigenomics.vdjtools.pool.PooledSample}.
 *
 * @param <T> clonotype aggregator type.
 */
public class SampleAggregator<T extends ClonotypeAggregator> implements Iterable<T> {
    private final Map<ClonotypeKey, T> innerMap = new HashMap<>();
    private final ClonotypeKeyGen clonotypeKeyGen;
    private final long count;

    /**
     * Aggregates clonotypes from the specified set of samples. In case clonotypes are matching
     * exactly (see {@link com.antigenomics.vdjtools.overlap.OverlapType#Strict}) in several samples
     * they are pooled in a single "pooled clonotype".
     *
     * @param samples                    a set of samples.
     * @param clonotypeAggregatorFactory clonotype aggregation rule.
     */
    public SampleAggregator(Iterable<Sample> samples,
                            ClonotypeAggregatorFactory<T> clonotypeAggregatorFactory) {
        this(samples, clonotypeAggregatorFactory, OverlapType.Strict);
    }

    /**
     * Aggregates clonotypes from the specified set of samples. Both clonotypes matching between samples
     * and convergent variants of each clonotype are pooled according to specified clonotype matching rule.
     *
     * @param samples                    a set of samples.
     * @param clonotypeAggregatorFactory clonotype aggregation rule.
     * @param overlapType                clonotype matching rule.
     */
    public SampleAggregator(Iterable<Sample> samples,
                            ClonotypeAggregatorFactory<T> clonotypeAggregatorFactory,
                            OverlapType overlapType) {
        this.clonotypeKeyGen = new ClonotypeKeyGen(overlapType);
        int sampleId = 0;
        long count = 0;

        for (Sample sample : samples) {
            System.out.println("[" + (new Date().toString()) + " " + "SamplePool] " +
                    "Pooling sample " + sample.getSampleMetadata().getSampleId());

            for (Clonotype clonotype : sample) {
                ClonotypeKey clonotypeKey = clonotypeKeyGen.generateKey(clonotype);

                ClonotypeAggregator clonotypeAggregator = getAt(clonotypeKey);
                if (clonotypeAggregator == null) {
                    innerMap.put(clonotypeKey, clonotypeAggregatorFactory.create(clonotype, sampleId));
                } else {
                    clonotypeAggregator.combine(clonotype, sampleId);
                }
            }

            count += sample.getCount();
            sampleId++;
        }

        // todo: sort

        this.count = count;
    }

    private T getAt(ClonotypeKey clonotypeKey) {
        return innerMap.get(clonotypeKey);
    }

    public T getAt(Clonotype clonotype) {
        return getAt(clonotypeKeyGen.generateKey(clonotype));
    }

    public int getDiversity() {
        return innerMap.size();
    }

    public long getCount() {
        return count;
    }

    @Override
    public Iterator<T> iterator() {
        return innerMap.values().iterator();
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
