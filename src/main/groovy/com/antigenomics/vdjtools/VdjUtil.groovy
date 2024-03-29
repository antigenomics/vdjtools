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

package com.antigenomics.vdjtools

import com.antigenomics.vdjtools.io.FileInputStreamFactory
import com.antigenomics.vdjtools.io.SampleStreamConnection
import com.antigenomics.vdjtools.misc.Software
import com.antigenomics.vdjtools.sample.ClonotypeFilter
import com.antigenomics.vdjtools.sample.CompositeClonotypeFilter
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SampleCollection

class VdjUtil {
    /**
     * Loads a sample collection using custom metadata file.
     * File should contain two columns: first with file path and second with sample id.
     * Additional columns will be stored as metadata entries.
     * First line of file should contain header that includes metadata field names.
     * Samples will be ordered as they appear in file.
     * @param sampleMetadataFileName metadata file path.
     * @param software software used to get processed samples.
     * @param store if set to true, all loaded samples will be stored in memory (only has effect if lazy is set to true).
     * @param lazy if set to true, all samples will be immediately loaded, otherwise samples will be loaded by request.
     * @param strict if set to false, will ignore samples with missing files, otherwise will throw an exception in such case.
     * @param sort whether to sort sample metadata by sample id.
     */
    static SampleCollection loadAll(String sampleMetadataFileName, Software software = Software.VDJtools,
                                    boolean store = true, boolean lazy = true, boolean strict = true, boolean sort = false
    ) {
        new SampleCollection(sampleMetadataFileName, software, store, lazy, strict, sort)
    }

    /**
     * Builds a sample collection from a pre-defined list of sample file names.
     * Samples will be assigned to generic metadata table, sample order will be preserved.
     * @param sampleFileNames list of sample file names
     * @param store if set to true, all loaded samples will be stored in memory (only has effect if lazy is set to true)
     * @param lazy if set to true, all samples will be immediately loaded, otherwise samples will be loaded by request
     * @param strict if set to false, will ignore samples with missing files, otherwise will throw an exception in such case
     * @param sort not sort sample metadata by sample id 
     */
    static SampleCollection loadAll(List<String> sampleFileNames, Software software = Software.VDJtools,
                                    boolean store = true, boolean lazy = true, boolean strict = true, boolean sort = false
    ) {
        new SampleCollection(sampleFileNames, software, store, lazy, strict, sort)
    }

    static SampleCollection toCollection(List<Sample> samples) {
        SampleCollection.fromSampleList(samples)
    }

    static Sample load(String fileName, Software software) {
        SampleStreamConnection.load(new FileInputStreamFactory(fileName), software)
    }

    static Sample filter(Sample sample, ClonotypeFilter... clonotypeFilter) {
        new Sample(sample, new CompositeClonotypeFilter(clonotypeFilter))
    }

    static SampleCollection filter(SampleCollection sampleCollection, ClonotypeFilter... clonotypeFilter) {
        SampleCollection.fromSampleList(
                sampleCollection.collect { new Sample(it, new CompositeClonotypeFilter(clonotypeFilter)) }
        )
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
