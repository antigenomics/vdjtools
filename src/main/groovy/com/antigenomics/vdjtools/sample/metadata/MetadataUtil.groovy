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

import org.apache.commons.io.FilenameUtils

/**
 * Some useful utils for metadata manipulation 
 */
public class MetadataUtil {
    private static final Map<String, Integer> sampleHash = new HashMap<>()

    /**
     * Converts a file name to sample id 
     * @param fileName file name to convert
     * @return sample id, a shortcut for file name without any path and extension
     */
    public static String fileName2id(String fileName) {
        FilenameUtils.getBaseName(
                fileName.endsWith(".gz") ?
                        FilenameUtils.getBaseName(fileName) :
                        fileName)
    }

    /**
     * Creates sample metadata object and assigns it to a generic metadata table 
     * @param sampleId short unique identifier of a sample
     * @return sample metadata object assigned to a generic metadata table
     */
    public static SampleMetadata createSampleMetadata(String sampleId) {
        def idCount = (sampleHash[sampleId] ?: 0) + 1
        sampleHash.put(sampleId, idCount)
        defaultMetadataTable.createRow((idCount > 0 ? "$idCount." : "") + sampleId, new ArrayList<String>())
    }

    /**
     * Gets a generic metadata table
     * @return metadata table which contains all statically-created metadata entries
     */
    public static MetadataTable getDefaultMetadataTable() {
        MetadataTable.GENERIC_METADATA_TABLE
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
