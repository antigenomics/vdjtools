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
import com.antigenomics.vdjtools.ClonotypeWrapperContainer;
import com.antigenomics.vdjtools.join.key.*;
import com.antigenomics.vdjtools.overlap.OverlapType;
import com.antigenomics.vdjtools.sample.Clonotype;

import java.util.HashSet;
import java.util.Set;

/**
 * Clonotype key generator implementing a certain clonotype matching rule.
 */
public class ClonotypeKeyGen {
    private final OverlapType overlapType;

    /**
     * Creates a clonotype key generator with {@link OverlapType#Strict} clonotype matching rule.
     */
    public ClonotypeKeyGen() {
        this(OverlapType.Strict);
    }

    /**
     * Creates a clonotype key generator for a specified clonotype matching rule.
     *
     * @param overlapType clonotype matching rule.
     */
    public ClonotypeKeyGen(OverlapType overlapType) {
        this.overlapType = overlapType;
    }

    /**
     * Generates keys for all clonotypes in a given sample.
     *
     * @param clonotypeWrapperContainer a sample.
     * @return a set of clonotype keys.
     */
    public Set<ClonotypeKey> generateKeySet(ClonotypeWrapperContainer<? extends ClonotypeWrapper> clonotypeWrapperContainer) {
        Set<ClonotypeKey> keySet = new HashSet<>();
        for (ClonotypeWrapper clonotypeWrapper : clonotypeWrapperContainer) {
            keySet.add(generateKey(clonotypeWrapper));
        }
        return keySet;
    }

    /**
     * Generates a key for a given clonotype wrapper under specified matching rule.
     *
     * @param clonotypeWrapper a clonotype wrapper.
     * @return clonotype key.
     */
    public ClonotypeKey generateKey(ClonotypeWrapper clonotypeWrapper) {
        return generateKey(clonotypeWrapper.getClonotype());
    }

    /**
     * Generates a key for a given clonotype under specified matching rule.
     *
     * @param clonotype a clonotype.
     * @return clonotype key.
     */
    public ClonotypeKey generateKey(Clonotype clonotype) {
        switch (overlapType) {
            case Nucleotide:
                return new NtKey(clonotype);

            case NucleotideV:
                return new NtVKey(clonotype);

            case NucleotideVJ:
                return new NtVJKey(clonotype);

            case AminoAcid:
                return new AaKey(clonotype);

            case AminoAcidV:
                return new AaVKey(clonotype);

            case AminoAcidVJ:
                return new AaVJKey(clonotype);

            case AminoAcidNonNucleotide:
                return new AaNotNtKey(clonotype);

            case Strict:
                return new StrictKey(clonotype);

            default:
                throw new UnsupportedOperationException();
        }
    }

    /**
     * Gets the clonotype matching rule for this key generator.
     *
     * @return clonotype matching rule.
     */
    public OverlapType getOverlapType() {
        return overlapType;
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
