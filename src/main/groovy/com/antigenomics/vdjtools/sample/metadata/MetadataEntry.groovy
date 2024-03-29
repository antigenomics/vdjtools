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

package com.antigenomics.vdjtools.sample.metadata

import groovy.transform.PackageScope

/**
 * Metadata entry, an object holding metadata table value that is associated with a given sample and column 
 */
public class MetadataEntry {
    public final MetadataTable grandParent
    public final SampleMetadata parent
    public final String columnId
    public String value

    /**
     * Creates new metadata entry
     * @param grandParent metadata table
     * @param parent sample metadata, a set of column values for a given sample (e.g. metadata table row)
     * @param columnId id of metadata table column
     * @param value metadata table value that is associated with a given sample and column
     */
    public MetadataEntry(MetadataTable grandParent, SampleMetadata parent, String columnId, String value) {
        this.grandParent = grandParent
        this.parent = parent
        this.columnId = columnId
        this.value = value
    }

    /**
     * INTERNAL re-assignes a given entry to new metadata table and sample 
     * @param grandParent
     * @param parent
     * @return
     */
    @PackageScope
    MetadataEntry changeParent(MetadataTable grandParent, SampleMetadata parent) {
        new MetadataEntry(grandParent, parent, columnId, value)
    }

    /**
     * Checks if a given metadata entry contains a numeric value
     * @return {@code true} if value could be converted to number, {@code false} otherwise
     */
    public boolean isNumeric() {
        value.isDouble()
    }

    /**
     * Converts the metadata entry value to double-precision number
     * @return a double precision value if sample is numeric, {@code Double.NaN} otherwise
     */
    public double asNumeric() {
        numeric ? value.toDouble() : Double.NaN
    }

    @Override
    public boolean equals(o) {
        if (this.is(o)) return true
        if (getClass() != o.class) return false

        MetadataEntry that = (MetadataEntry) o

        if (columnId != that.columnId) return false
        if (value != that.value) return false

        return true
    }

    @Override
    public int hashCode() {
        int result
        result = columnId.hashCode()
        result = 31 * result + value.hashCode()
        return result
    }

    @Override
    public String toString() {
        value
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
