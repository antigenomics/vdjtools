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
 * A clonotype key, which implements {@link #equals} and {@link #hashCode}
 * according to specified {@link com.antigenomics.vdjtools.overlap.OverlapType}.
 * {@see com.antigenomics.vdjtools.join.key.ClonotypeKey}
 */
public abstract class ClonotypeKey {
    protected final Clonotype clonotype;

    public ClonotypeKey(Clonotype clonotype) {
        this.clonotype = clonotype;
    }

    abstract boolean equals(Clonotype other);

    @Override
    public abstract int hashCode();

    @Override
    public boolean equals(Object o) {
        // no it shouldn't
        return o != null && this.equals(((ClonotypeKey) o).clonotype);
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
