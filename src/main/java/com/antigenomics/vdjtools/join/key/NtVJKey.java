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

package com.antigenomics.vdjtools.join.key;

import com.antigenomics.vdjtools.sample.Clonotype;

/**
 * {@inheritDoc}
 */
public class NtVJKey extends ClonotypeKey {
    public NtVJKey(Clonotype clonotype) {
        super(clonotype);
    }

    @Override
    public boolean equals(Clonotype other) {
        return clonotype.getCdr3ntBinary().equals(other.getCdr3ntBinary()) &&
                clonotype.getVBinary().equals(other.getVBinary()) &&
                clonotype.getJBinary().equals(other.getJBinary());
    }

    @Override
    public int hashCode() {
        return 31 * (clonotype.getCdr3ntBinary().hashCode() * 31 + clonotype.getVBinary().hashCode()) +
                clonotype.getJBinary().hashCode();
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
