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

import com.antigenomics.vdjtools.io.FileInputStreamFactory
import com.antigenomics.vdjtools.io.SampleFileConnection
import com.antigenomics.vdjtools.sample.SampleCollection
import groovy.transform.PackageScope

import static com.antigenomics.vdjtools.misc.ExecUtil.*

/**
 * Base class that handles metadata for a sample collection
 */
class MetadataTable implements Iterable<SampleMetadata> {
    /**
     * Generic metadata table. All statically created SampleMetadata objects are internally dependent on it
     */
    @PackageScope
    static final MetadataTable GENERIC_METADATA_TABLE = new MetadataTable()

    private final HashMap<String, Integer> columnId2Index = new HashMap<>()
    private final List<String> columnIds
    private List<String> sampleOrder = new ArrayList<>()
    private final Map<String, SampleMetadata> metadataBySample = new HashMap<>()

    /**
     * Creates an empty metadata table
     */
    MetadataTable() {
        this.columnIds = new ArrayList<>()
    }

    /**
     * Creates a metadata table with pre-defined list of columns
     * @param columnIds column names
     */
    MetadataTable(List<String> columnIds) {
        this.columnIds = columnIds

        columnIds.eachWithIndex { it, ind ->
            if (columnId2Index.containsKey(it))
                throw new IllegalArgumentException("Duplicate metadata columns not allowed")

            columnId2Index.put(it, ind)
        }
    }

    /**
     * Creates a deep copy of this metadata
     * @return
     */
    public MetadataTable copy() {
        def metadataTable = new MetadataTable(columnIds.collect())
        this.each {
            metadataTable.addSample(it)
        }
        metadataTable
    }

    /**
     * Selects a subset of samples from this metadata table based on a specified rule.
     * This method creates a deep copy of metadata table.
     * @param entryFilter metadata entry filtering rule.
     * @param sampleIds sample (row) IDs to keep.
     * @return a copy of metadata table containing selected samples (rows) only.
     */
    public MetadataTable select(MetadataEntryFilter entryFilter = BlankMetadataEntryFilter.INSTANCE,
                                Set<String> sampleIds = sampleIds) {
        def metadataTable = new MetadataTable(columnIds.collect())
        def sampleFilter = new SampleMetadataFilter(entryFilter)
        this.findAll {
            if (!this.sampleIds.contains(it.sampleId)) {
                throw new RuntimeException("Sample $it.sampleId is missing in sample collection." +
                        "Current sample list: ${this.sampleIds.join(",")}.")
            }
            sampleIds.contains(it.sampleId) && sampleFilter.passes(it)
        }.each {
            metadataTable.addSample(it)
        }
        metadataTable
    }

    /**
     * Reorders sample list according to provided order
     * @param sampleOrder list of sample ids in new order
     */
    void sort(List<String> sampleOrder) {
        if (sampleOrder.size() != sampleCount)
            throw new IllegalArgumentException("Bad sample order, " +
                    "number of sample ids and sample collection size don't match")

        sampleOrder.each {
            if (!columnId2Index.containsKey(it))
                throw new IllegalArgumentException("Bad sample order, unknown sample id $it")
        }

        this.sampleOrder = sampleOrder
    }

    /**
     * Sorts samples in metadata according to sample id.
     */
    void sort() {
        sampleOrder.sort()
    }

    /**
     * Sorts samples in metadata according to a given column in ascending order.
     * Sorting is performed according to column type (Factor/Numeric/Semi-Numeric).
     * NaNs for numeric columns are always put last
     * @param columnId column identifier
     */
    void sort(String columnId) {
        sort(columnId, true)
    }

    /**
     * Sorts samples in metadata according to a given column.
     * Sorting is performed according to column type (Factor/Numeric/Semi-Numeric).
     * NaNs for numeric columns are always put last
     * @param columnId column identifier
     * @param ascending true for ascending, false for descending sort
     */
    void sort(String columnId, boolean ascending) {
        checkColumnId(columnId)

        def columnInfo = getInfo(columnId)

        switch (columnInfo.metadataColumnType) {
            case MetadataColumnType.Numeric:
            case MetadataColumnType.SemiNumeric:
                sampleOrder.sort {
                    def value = this[it, columnId].asNumeric()
                    Double.isNaN(value) ?
                            (ascending ? Double.MIN_VALUE : Double.MAX_VALUE) :
                            (ascending ? value : -value)
                }
                break
            case MetadataColumnType.Factor:
                sampleOrder.sort(true, new Comparator<String>() {
                    @Override
                    int compare(String o1, String o2) {
                        ascending ? o1.compareTo(o2) : -o1.compareTo(o2)
                    }
                })
                break
        }
    }

    /**
     * Add metadata column
     * @param columnId name of column
     * @param values column values for all samples in collection, in corresponding order
     */
    void addColumn(String columnId, List<String> values) {
        if (values.size() != sampleCount)
            throw new Exception("Number of values should be equal to number of samples")

        addColumnId(columnId)

        sampleOrder.eachWithIndex { it, ind ->
            metadataBySample[it].addEntry(columnId, values[ind])
        }
    }

    /**
     * Add metadata column
     * @param columnId name of column
     * @param value value used to fill the column
     */
    void addColumn(String columnId, String value) {
        addColumnId(columnId)

        sampleOrder.eachWithIndex { it, ind ->
            metadataBySample[it].addEntry(columnId, value)
        }
    }

    /**
     * Add metadata column
     * @param columnId name of column
     * @param sample2value map containing column values for
     * @param strict if true all samples should be matched. Otherwise non-existing samples will be ignored,
     *               while missing samples will take "NA" value
     */
    void addColumn(String columnId, Map<String, String> sample2value, boolean strict) {
        if (strict && sampleCount != sample2value.size())
            throw new Exception("Number of values should be equal to number of samples")

        addColumnId(columnId)

        sample2value.each {
            if (metadataBySample.containsKey(it.key)) {
                metadataBySample[it.key].addEntry(columnId, it.value)
            } else if (strict) {
                throw new Exception("Unknown sample $it.key")
            }
        }

        // fill missing values
        if (!strict) {
            metadataBySample.keySet().each {
                if (!sample2value.containsKey(it)) {
                    metadataBySample[it].addEntry(columnId, "NA")
                }
            }
        }
    }

    private void addColumnId(String columnId) {
        columnId2Index.put(columnId, columnIds.size())
        columnIds.add(columnId)
    }

    private void checkColumnId(String columnId) {
        if (!containsColumn(columnId))
            throw new IllegalArgumentException("Column $columnId doesn't exist")
    }

    public int getColumnIndex(String columnId) {
        columnId = columnId.toUpperCase()
        columnIds.findIndexOf { it.toUpperCase() == columnId }
    }

    /**
     * Check if metadata table contains given column
     * @param columnId specific column
     * @return true if metadata table contains given column, otherwise false
     */
    boolean containsColumn(String columnId) {
        columnIds.contains(columnId)
    }

    /**
     * Creates a sample row, appends the metadata to table and binds sample to a given table
     * @param sampleId sample name
     * @param rowValues values of corresponding metadata fields
     * @return sample metadata instance
     */
    SampleMetadata createRow(String sampleId, List<String> rowValues) {
        // Create metadata entries, assign current metadata as parent
        def entries = new ArrayList<MetadataEntry>()

        def sampleMetadata = SampleMetadata.predefined(sampleId, entries, this)

        rowValues.eachWithIndex { String value, int i ->
            entries.add(new MetadataEntry(this, sampleMetadata, columnIds[i], value))
        }

        // Record sample to metadata if all is fine

        addSample(sampleMetadata)

        sampleMetadata
    }

    /**
     * Creates an empty sample row, appends the metadata to table and binds sample to a given table
     * @param sampleId sample name
     * @return sample metadata instance
     */
    SampleMetadata createRow(String sampleId) {
        createRow(sampleId, new ArrayList<String>())
    }

    /**
     * Adds a single sample metadata to given metadata table 
     * @param row sample metadata
     */
    void addSample(SampleMetadata row) {
        // Extensive checks
        checkSample(row)

        // Copy & change ownership
        row.changeParent(this)

        // Record this sample
        sampleOrder.add(row.sampleId)

        metadataBySample.put(row.sampleId, row)
    }

    private void checkSample(SampleMetadata row) {
        if (metadataBySample.containsKey(row.sampleId))
            throw new IllegalArgumentException("Duplicate samples not allowed")

        if (row.entries.size() != columnIds.size())
            throw new Exception("Bad metadata row\n$row\n- number of entries in row " +
                    "and in parent metadata should agree")

        for (int i = 0; i < columnIds.size(); i++) {
            def rowValue = row.entries[i].columnId, metadataValue = columnIds[i]
            if (rowValue != metadataValue)
                throw new Exception("Bad metadata row\n$row\n- column#${i} don't match between row " +
                        "($rowValue) and parent metadata ($metadataValue)")
        }
    }

    /**
     * Write metadata table copy to file assuming that one-to-one sample output will also be
     * placed to the same directory. Used by routines like Decontaminate, ApplySampleAsFilter or Downsample. 
     * @param outputPrefix output prefix
     * @param compress indicates whether samples will be stored as compressed
     * @param filters list of filter names applied to data this time
     */
    public void storeWithOutput(String outputPrefix, boolean compress, String... filters) {
        def metadataPath = formMetadataPath(outputPrefix)

        def metadataTableCopy = this.copy()

        filters.each { filter ->
            if (metadataTableCopy.containsColumn("..filter..")) {
                metadataTableCopy.getColumn("..filter..").each {
                    it.value += ",$filter"
                }
            } else {
                metadataTableCopy.addColumn("..filter..", filter)
            }
        }

        new File(metadataPath).withPrintWriter { pw ->
            pw.println("$FILE_NAME_COLUMN\t$SAMPLE_ID_COLUMN\t" + metadataTableCopy.columnHeader)
            metadataTableCopy.each {
                def sampleOutputPath = formOutputPath(outputPrefix, it.sampleId)

                pw.println([relativeSamplePath(metadataPath, sampleOutputPath) + (compress ? ".gz" : ""),
                            it.sampleId,
                            it.toString()].join("\t"))
            }
        }
    }

    /**
     * Write the metadata table to file, inheriting sample paths from provided sample collection.
     * Sample connections should be provided by {@link FileInputStreamFactory} for this method to 
     * work without throwing exception
     * @param outputPrefix output prefix
     * @param sampleCollection parent sample collection
     */
    public void storeWithOutput(String outputPrefix, SampleCollection sampleCollection, String splitterValue = null) {
        def metadataPath = formMetadataPath(outputPrefix, splitterValue)

        new File(metadataPath).createNewFile()

        new File(metadataPath).withPrintWriter { pw ->
            pw.println("$FILE_NAME_COLUMN\t$SAMPLE_ID_COLUMN\t" + this.columnHeader)
            this.eachWithIndex { it, ind ->
                def connection = sampleCollection.sampleMap[it.sampleId] as SampleFileConnection
                def fileInputStream = connection.inputStreamFactory as FileInputStreamFactory

                pw.println([relativeSamplePath(metadataPath, fileInputStream.fileName),
                            it.sampleId,
                            it.toString()].join("\t"))
            }
        }
    }

    /**
     * Gets info on values contained in a given metadata column
     * @param columnId column name
     * @return column info
     */
    MetadataColumnInfo getInfo(String columnId) {
        new MetadataColumnInfo(this, columnId)
    }

    /**
     * Gets all metadata entries in a given column
     * @param columnId column name
     * @return metadata entries
     */
    List<MetadataEntry> getColumn(String columnId) {
        if (columnId == SAMPLE_ID_COLUMN) {
            return sampleOrder.collect { new MetadataEntry(this, getRow(it), SAMPLE_ID_COLUMN, it) }
        }

        checkColumnId(columnId)

        int index = columnId2Index[columnId]

        getColumn(index)
    }

    /**
     * Gets all metadata entries in a given column
     * @param columnIndex column index
     * @return metadata entries
     */
    List<MetadataEntry> getColumn(int columnIndex) {
        if (columnIndex < 0 || columnIndex >= columnIds.size())
            throw new IndexOutOfBoundsException()

        sampleOrder.collect {
            metadataBySample[it].entries[columnIndex]
        }
    }

    /**
     * Gets row from metadata table
     * @param sampleId sample name
     * @return metadata entries for a given sample
     */
    SampleMetadata getRow(String sampleId) {
        metadataBySample[sampleId]
    }

    /**
     * Gets row from metadata table
     * @param sampleIndex sample index
     * @return metadata entries for a given sample
     */
    SampleMetadata getRow(int sampleIndex) {
        getRow(sampleOrder[sampleIndex])
    }

    /**
     * Iterate through sample names (row names)
     * @return
     */
    Iterator<String> getSampleIterator() {
        sampleOrder.iterator()
    }

    /**
     * Iterate through metadata entry ids (column names)
     * @return
     */
    Iterator<String> getColumnIterator() {
        columnIds.iterator()
    }

    /**
     * Gets a specific metadata entry
     * @param sampleId sample name
     * @param columnId metadata id
     * @return metadata entry
     */
    MetadataEntry getAt(String sampleId, String columnId) {
        if (!columnId2Index.containsKey(columnId))
            return null

        metadataBySample[sampleId].entries[columnId2Index[columnId]]
    }

    /**
     * Gets the number of different metadata entry types (columns) in this table
     * @return number of columns in the table
     */
    int getNumberOfColumns() {
        columnIds.size()
    }

    /**
     * Gets the number of different samples (rows) in this table
     * @return number of rows in the table
     */
    int getSampleCount() {
        metadataBySample.size()
    }

    Set<String> getSampleIds() {
        Collections.unmodifiableSet(metadataBySample.keySet())
    }

    public static final String SAMPLE_ID_COLUMN = "sample_id", FILE_NAME_COLUMN = "file_name"

    String getColumnHeader() {
        columnIds.size() > 0 ? columnIds.collect().join("\t") : "metadata_blank"
    }

    @PackageScope
    String getColumnHeader1() {
        columnIds.size() > 0 ? columnIds.collect { "1_$it" }.join("\t") : "1_metadata_blank"
    }

    @PackageScope
    String getColumnHeader2() {
        columnIds.size() > 0 ? columnIds.collect { "2_$it" }.join("\t") : "2_metadata_blank"
    }

    @Override
    String toString() {
        "\t" + columnIds.collect().join("\t") + "\n" +
                metadataBySample.collect { it.key + "\t" + it.value.entries.collect() }.join("\n")
    }

    @Override
    Iterator<SampleMetadata> iterator() {
        def orderIter = sampleOrder.iterator()
        [hasNext: { orderIter.hasNext() }, next: { metadataBySample[orderIter.next()] }] as Iterator
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
