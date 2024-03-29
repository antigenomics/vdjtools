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

package com.antigenomics.vdjtools.annotate

import com.antigenomics.vdjtools.annotate.partitioning.FullCdr3
import com.antigenomics.vdjtools.sample.Clonotype

class AAPropertyAnnotator implements ClonotypeAnnotator {
    static final String NORMALIZED_SUFFIX = ".avg"

    private final name
    private final AaPropertySummaryEvaluator propertyCalculator

    AAPropertyAnnotator(String name, boolean normalized) {
        this.propertyCalculator = new AaPropertySummaryEvaluator(KnownAminoAcidProperties.INSTANCE.getByName(name),
                new FullCdr3(), normalized, false)
        this.name = normalized ? (name + NORMALIZED_SUFFIX) : name
    }

    @Override
    String getName() {
        name
    }

    @Override
    String getCategory() {
        "aaprop"
    }

    @Override
    String annotate(Clonotype clonotype) {
        clonotype.coding ? (propertyCalculator.compute(clonotype)) : ""
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
