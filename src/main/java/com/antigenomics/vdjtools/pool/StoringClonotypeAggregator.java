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

import com.antigenomics.vdjtools.ClonotypeWrapper;
import com.antigenomics.vdjtools.join.ClonotypeKeyGen;
import com.antigenomics.vdjtools.join.key.ClonotypeKey;
import com.antigenomics.vdjtools.sample.Clonotype;

import java.util.HashSet;
import java.util.Set;

/**
 * A clonotype aggregator that stores the representative clonotype and acts as a clonotype wrapper,
 * a.k.a. pooled clonotype.
 */
public class StoringClonotypeAggregator extends MaxClonotypeAggregator implements ClonotypeWrapper {
    private Clonotype clonotype;
    private int convergence, occurrences;
    private PooledSample parent;
    private final Set<ClonotypeKey> variants = new HashSet<>();
    private final ClonotypeKeyGen clonotypeKeyGen = new ClonotypeKeyGen();

    protected StoringClonotypeAggregator(Clonotype clonotype, int sampleId) {
        super(clonotype, sampleId);
        this.clonotype = clonotype;
        this.convergence = 1;
        this.occurrences = 1;
        this.variants.add(clonotypeKeyGen.generateKey(clonotype));
    }

    @Override
    protected boolean tryReplace(Clonotype clonotype, int sampleId) {
        ClonotypeKey clonotypeKey = clonotypeKeyGen.generateKey(clonotype);

        if (!variants.contains(clonotypeKey)) {
            convergence++;
            variants.add(clonotypeKey);
        }
        occurrences++;

        if (super.tryReplace(clonotype, sampleId)) {
            this.clonotype = clonotype;
            return true;
        }
        return false;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public Clonotype getClonotype() {
        return clonotype;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public PooledSample getParent() {
        return parent;
    }

    public void setParent(PooledSample parent) {
        this.parent = parent;
    }

    /**
     * Gets the total number of convergent variants for this pooled clonotype, i.e.
     * the total number of unique nucleotide variants in all samples.
     *
     * @return number of convergent variants.
     */
    @Override
    public int getDiversity() {
        return convergence;
    }

    /**
     * Gets the total number of clonotype occurrences - counting the total number of
     * nucleotide variants in all samples. Differs from diversity/convergence as the
     * same nucleotide variant found in two samples will be counted twice.
     *
     * @return number of occurrences
     */
    public int getOccurrences() {
        return occurrences;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getFreq() {
        return parent != null ? (getCount() / (double) parent.getCount()) : Double.NaN;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getFreqAsInInput() {
        return 1.0;
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
