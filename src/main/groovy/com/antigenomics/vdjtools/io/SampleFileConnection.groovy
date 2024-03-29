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

package com.antigenomics.vdjtools.io

import com.antigenomics.vdjtools.misc.Software
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.metadata.SampleMetadata

/**
 * A wrapper for plain-text clonotype table stored in a file.
 * This is a semi-internal class to provide lazy-loading support for SampleCollection
 */
public class SampleFileConnection extends SampleStreamConnection {
    /**
     * Loads a sample from the specified file.
     * @param fileName path to file containing the sample.
     * @return sample object filled with clonotypes.
     */
    public static Sample load(String fileName) {
        load(fileName, Software.VDJtools)
    }

    /**
     * Loads a sample from the specified file.
     * @param fileName path to file containing the sample.
     * @param software type of software used to create clonotype table. Specifies how the plain-text input will be parsed.
     * @return sample object filled with clonotypes.
     */
    public static Sample load(String fileName, Software software) {
        load(new FileInputStreamFactory(fileName), software)
    }

    /**
     * Loads a sample from the specified file.
     * @param fileName path to file containing the sample.
     * @param software type of software used to create clonotype table. Specifies how the plain-text input will be parsed.
     * @param sampleMetadata a metadata object that will be associated with a given sample.
     * @return sample object filled with clonotypes.
     */
    public static Sample load(String fileName, Software software, SampleMetadata sampleMetadata) {
        load(new FileInputStreamFactory(fileName), software, sampleMetadata)
    }

    /**
     * Creates a sample connection, an object that could be used to access (load to memory, store, etc) a sample stored as plain text.
     * Will load the sample upon initialization and store it into memory. Generic sample metadata will be associated with the underlying sample.
     * @param fileName path to file containing the sample.
     * @param software type of software used to create clonotype table. Specifies how the plain-text input will be parsed.
     */
    public SampleFileConnection(String fileName, Software software) {
        super(new FileInputStreamFactory(fileName), software)
    }

    /**
     * Creates a sample connection, an object that could be used to access (load to memory, store, etc) a sample stored as plain text.
     * Will load the sample upon initialization and store it into memory.
     * @param fileName path to file containing the sample.
     * @param software type of software used to create clonotype table. Specifies how the plain-text input will be parsed.
     * @param sampleMetadata a metadata object that will be associated with a given sample.
     */
    public SampleFileConnection(String fileName, Software software, SampleMetadata sampleMetadata) {
        this(fileName, software, sampleMetadata, false, true)
    }

    /**
     * Creates a sample connection, an object that could be used to access (load to memory, store, etc) a sample stored as plain text.
     * @param fileName path to file containing the sample.
     * @param software type of software used to create clonotype table. Specifies how the plain-text input will be parsed.
     * @param sampleMetadata a metadata object that will be associated with a given sample.
     * @param lazy will load the sample only when {@code getSample ( )} is called.
     * @param store sample will be stored into memory after loading.
     */
    public SampleFileConnection(String fileName, Software software, SampleMetadata sampleMetadata,
                                boolean lazy, boolean store) {
        super(new FileInputStreamFactory(fileName), software, sampleMetadata, lazy, store)
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
