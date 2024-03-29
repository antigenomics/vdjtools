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
import com.antigenomics.vdjtools.misc.ExecUtil

import static com.antigenomics.vdjtools.sample.metadata.MetadataUtil.createSampleMetadata

/**
 * A wrapper for plain-text clonotype table that could be accessed using {@code InputStream}.
 * This is a semi-internal class to provide lazy-loading support for SampleCollection.
 */
public class SampleStreamConnection implements SampleConnection {
    private final InputStreamFactory inputStreamFactory
    private final SampleMetadata sampleMetadata
    private final Software software
    private Sample _sample = null
    private final boolean lazy, store

    /**
     * Loads a sample from a given input stream provider.
     * @param inputStreamFactory an input stream provider.
     * @return sample object filled with clonotypes.
     */
    public static Sample load(InputStreamFactory inputStreamFactory) {
        load(inputStreamFactory, Software.VDJtools)
    }

    /**
     * Loads a sample from a given input stream provider.
     * @param inputStreamFactory an input stream provider.
     * @param software type of software used to create clonotype table. Specifies how the plain-text input will be parsed.
     * @return sample object filled with clonotypes.
     */
    public static Sample load(InputStreamFactory inputStreamFactory, Software software) {
        new SampleStreamConnection(inputStreamFactory, software).sample
    }

    /**
     * Loads a sample from a given input stream provider.
     * @param inputStreamFactory an input stream provider.
     * @param software type of software used to create clonotype table. Specifies how the plain-text input will be parsed.
     * @param sampleMetadata a metadata object that will be associated with a given sample.
     * @return sample object filled with clonotypes.
     */
    public static Sample load(InputStreamFactory inputStreamFactory, Software software, SampleMetadata sampleMetadata) {
        new SampleStreamConnection(inputStreamFactory, software, sampleMetadata).sample
    }

    /**
     * Creates a sample connection, an object that could be used to access (load to memory, store, etc) a sample stored as plain text.
     * Will load the sample upon initialization and store it into memory. Generic sample metadata will be associated with the underlying sample.
     * @param inputStreamFactory an input stream provider.
     * @param software type of software used to create clonotype table. Specifies how the plain-text input will be parsed.
     */
    public SampleStreamConnection(InputStreamFactory inputStreamFactory, Software software) {
        this(inputStreamFactory, software, createSampleMetadata(inputStreamFactory.getId()))
    }

    /**
     * Creates a sample connection, an object that could be used to access (load to memory, store, etc) a sample stored as plain text.
     * Will load the sample upon initialization and store it into memory.
     * @param inputStreamFactory an input stream provider.
     * @param software type of software used to create clonotype table. Specifies how the plain-text input will be parsed.
     * @param sampleMetadata a metadata object that will be associated with a given sample.
     */
    public SampleStreamConnection(InputStreamFactory inputStreamFactory, Software software, SampleMetadata sampleMetadata) {
        this(inputStreamFactory, software, sampleMetadata, false, true)
    }

    /**
     * Creates a sample connection, an object that could be used to access (load to memory, store, etc) a sample stored as plain text.
     * @param inputStreamFactory an input stream provider.
     * @param software type of software used to create clonotype table. Specifies how the plain-text input will be parsed.
     * @param sampleMetadata a metadata object that will be associated with a given sample.
     * @param lazy will load the sample only when {@code getSample ( )} is called.
     * @param store sample will be stored into memory after loading.
     */
    public SampleStreamConnection(InputStreamFactory inputStreamFactory, Software software, SampleMetadata sampleMetadata,
                                  boolean lazy, boolean store) {
        this.inputStreamFactory = inputStreamFactory
        this.sampleMetadata = sampleMetadata
        this.software = software
        this.lazy = lazy
        this.store = store

        if (!lazy)
            _sample = _load()
    }

    /**
     * INTERNAL loads a given sample into memory.
     * @return loaded sample filled with clonotypes.
     */
    private Sample _load() {
        println "[${new Date()} SampleStreamConnection] Loading sample $sampleMetadata.sampleId"
        def inputStream = inputStreamFactory.create()
        def sample = Sample.fromInputStream(inputStream, sampleMetadata, software, -1, true, software.collapseRequired)
        println "[${new Date()} SampleStreamConnection] Loaded sample $sampleMetadata.sampleId with " +
                "$sample.diversity clonotypes and $sample.count cells. " + ExecUtil.memoryFootprint()
        sample
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public Sample getSample() {
        _sample ?: (store ? (_sample = _load()) : _load())
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public Sample haveAGlance() {
        _sample ?: Sample.fromInputStream(inputStreamFactory.create(), sampleMetadata, software, -1, false, false)
    }

    @Override
    public String toString() {
        "SampleFileConnection{$inputStreamFactory.id>${sampleMetadata.sampleId},storing=${store},loaded=${_sample != null}"
    }

    InputStreamFactory getInputStreamFactory() {
        return inputStreamFactory
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
