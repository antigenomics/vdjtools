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

import com.antigenomics.vdjtools.sample.Clonotype;
import com.milaboratory.core.sequence.AminoAcidSequence;

public class Cdr3ContactEstimate extends AaProperty {
    private final float[][] contactProbsTRA, contactProbsTRB;
    private final int l;

    public Cdr3ContactEstimate(float[][] contactProbsTRA, float[][] contactProbsTRB) {
        this.contactProbsTRA = contactProbsTRA;
        this.contactProbsTRB = contactProbsTRB;

        if (contactProbsTRA.length != AminoAcidSequence.ALPHABET.size()) {
            throw new IllegalArgumentException("Number of rows of contactProbs TRA " +
                    "be equal to AminoAcidSequence.ALPHABET size.");
        }

        if (contactProbsTRB.length != AminoAcidSequence.ALPHABET.size()) {
            throw new IllegalArgumentException("Number of rows of contactProbs TRB " +
                    "be equal to AminoAcidSequence.ALPHABET size.");
        }

        l = contactProbsTRA[0].length;

        for (int i = 1; i < contactProbsTRA.length; i++) {
            if (contactProbsTRA[i].length != l) {
                throw new IllegalArgumentException("Number of columns of " +
                        "contactProbs TRA should be the same for all amino acids");
            }
        }

        for (int i = 0; i < contactProbsTRB.length; i++) {
            if (contactProbsTRB[i].length != l) {
                throw new IllegalArgumentException("Number of columns of " +
                        "contactProbs TRB should be the same for all amino acids");
            }
        }
    }

    @Override
    public String getName() {
        return "cdr3contact";
    }

    @Override
    public float compute(AminoAcidSequence sequence, int pos) {
        return Float.NaN;
    }

    @Override
    public float compute(Clonotype clonotype, int pos) {
        AminoAcidSequence sequence = clonotype.getCdr3aaBinary();
        byte aa = sequence.codeAt(pos);

        int posBin = (int) ((l - 1) * (pos / (float) (sequence.size() - 1)));

        switch (clonotype.getVBinary().getChain()) {
            case TRA:
                return contactProbsTRA[aa][posBin];
            case TRB:
                return contactProbsTRB[aa][posBin];
        }

        return compute(sequence, pos);
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
