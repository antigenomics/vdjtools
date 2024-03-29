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
import com.antigenomics.vdjtools.sample.metadata.MetadataEntryExpressionFilter

def cli = new CliBuilder(usage: "FilterMetadata [options] metadata.txt output_dir output_suffix")
cli.h("display help message")
cli.f(longOpt: "filter", argName: "string", args: 1, required: true,
        "Filter expression, metadata column names should be marked with ${MetadataEntryExpressionFilter.FILTER_MARK}, " +
                "e.g. \"__chain__=~/TR[AB]/\" or \"__chain__=='TRA'||__chain__=='TRB'\"")

def opt = cli.parse(args)

if (opt == null) {
    System.exit(2)
}

if (opt.h || opt.arguments().size() != 3) {
    cli.usage()
    System.exit(2)
}

def metadataFileName = opt.arguments()[0], filter = (String) opt.f,
    outputDir = ExecUtil.toDirPath(opt.arguments()[1]), outputSuffix = opt.arguments()[2]

def scriptName = getClass().canonicalName.split("\\.")[-1]

// Lazy load sample list, need to get absolute paths
println "[${new Date()} $scriptName] Checking sample(s)"
def sampleCollection = new SampleCollection((String) metadataFileName, Software.VDJtools, false)

println "[${new Date()} $scriptName] Filtering metadata by $filter"
def filteredMetadataTable = sampleCollection.metadataTable.select(new MetadataEntryExpressionFilter(filter))
filteredMetadataTable.storeWithOutput(outputDir, sampleCollection, outputSuffix)
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
