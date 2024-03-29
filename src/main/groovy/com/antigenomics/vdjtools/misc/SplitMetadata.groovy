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

package com.antigenomics.vdjtools.misc

import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.sample.metadata.BlankMetadataEntryFilter

def cli = new CliBuilder(usage: "SplitMetadata [options] metadata.txt output_dir")
cli.h("display help message")
cli.c(longOpt: "columns", argName: "string1,string2,...", args: 1, required: true,
        "Column name(s) to split metadata by.")

def opt = cli.parse(args)

if (opt == null) {
    System.exit(2)
}

if (opt.h || opt.arguments().size() != 2) {
    cli.usage()
    System.exit(2)
}

def metadataFileName = opt.arguments()[0], columnIds = ((String) opt.c).split(","),
    outputDir = ExecUtil.toDirPath(opt.arguments()[1])

def scriptName = getClass().canonicalName.split("\\.")[-1]

// Lazy load sample list, need to get absolute paths
println "[${new Date()} $scriptName] Checking sample(s)"
def sampleCollection = new SampleCollection((String) metadataFileName, Software.VDJtools, false)

println "[${new Date()} $scriptName] Splitting metadata by $columnIds"

def sampleIdByMetadataValue = new HashMap<String, List<String>>()

sampleCollection.metadataTable.each { sampleMetadata ->
    def key = columnIds.collect { sampleMetadata[it].value }.join(".")
    def sampleList = sampleIdByMetadataValue[key]
    if (sampleList == null) {
        sampleIdByMetadataValue.put(key, sampleList = new ArrayList<String>())
    }
    sampleList.add(sampleMetadata.sampleId)
}

sampleIdByMetadataValue.each {
    def filteredMetadataTable = sampleCollection.metadataTable.select(BlankMetadataEntryFilter.INSTANCE,
            new HashSet<String>(it.value))
    filteredMetadataTable.storeWithOutput(outputDir, sampleCollection, it.key)
}

println "[${new Date()} $scriptName] Finished"
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
