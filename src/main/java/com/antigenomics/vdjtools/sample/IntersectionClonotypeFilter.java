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

package com.antigenomics.vdjtools.sample;

import com.antigenomics.vdjtools.join.ClonotypeKeyGen;
import com.antigenomics.vdjtools.join.key.ClonotypeKey;
import com.antigenomics.vdjtools.overlap.OverlapType;

import java.util.Set;

/**
 * Filter based on clonotype intersection with the list of clonotypes from the specified sample according to
 * the specified clonotype matching rule.
 */
public class IntersectionClonotypeFilter extends ClonotypeFilter {
    private final ClonotypeKeyGen clonotypeKeyGen;
    private final Set<ClonotypeKey> keySet;

    public IntersectionClonotypeFilter(OverlapType overlapType, Sample sample, boolean negative) {
        super(negative);
        this.clonotypeKeyGen = new ClonotypeKeyGen(overlapType);
        this.keySet = new ClonotypeKeyGen(overlapType).generateKeySet(sample);
    }

    public IntersectionClonotypeFilter(OverlapType overlapType, Sample sample) {
        this(overlapType, sample, false);
    }

    @Override
    protected boolean checkPass(Clonotype clonotype) {
        return keySet.contains(clonotypeKeyGen.generateKey(clonotype));
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
