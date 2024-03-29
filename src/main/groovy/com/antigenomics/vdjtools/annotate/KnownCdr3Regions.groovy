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

import com.antigenomics.vdjtools.annotate.partitioning.Cdr3Center
import com.antigenomics.vdjtools.annotate.partitioning.Cdr3Region
import com.antigenomics.vdjtools.annotate.partitioning.DGermline
import com.antigenomics.vdjtools.annotate.partitioning.DJJunction
import com.antigenomics.vdjtools.annotate.partitioning.FullCdr3
import com.antigenomics.vdjtools.annotate.partitioning.JGermline
import com.antigenomics.vdjtools.annotate.partitioning.VDJunction
import com.antigenomics.vdjtools.annotate.partitioning.VGermline
import com.antigenomics.vdjtools.annotate.partitioning.VJJunction

class KnownCdr3Regions {
    private final Map<String, Cdr3Region> regionsByName

    static final KnownCdr3Regions INSTANCE = new KnownCdr3Regions()

    private KnownCdr3Regions() {
        this.regionsByName = [new VGermline(), new DGermline(), new JGermline(),
                              new VDJunction(), new DJJunction(),
                              new VJJunction(), new FullCdr3(),
                              Cdr3Center.CDR3_CENTER_5,
                              Cdr3Center.CDR3_CENTER_3].collectEntries {
            [(it.name.toLowerCase()): it]
        }
    }

    Cdr3Region getByName(String name) {
        name = name.toLowerCase()
        if (!regionsByName.containsKey(name))
            throw new IllegalArgumentException("Bad CDR3 region name '$name', allowed values: ${allowedNames}")
        regionsByName[name]
    }

    List<String> getAllowedNames() {
        regionsByName.keySet().collect()
    }

    List<Cdr3Region> getAll() {
        regionsByName.values().collect()
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
