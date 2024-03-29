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

package com.antigenomics.vdjtools.preprocess

import com.antigenomics.vdjtools.io.SampleWriter
import com.antigenomics.vdjtools.sample.*
import com.antigenomics.vdjtools.sample.metadata.MetadataTable

import static com.antigenomics.vdjtools.misc.ExecUtil.formOutputPath

def cli = new CliBuilder(usage: "FilterBySegment [options] " +
        "[sample1 sample2 sample3 ... if -m is not specified] output_prefix")
cli.h("Display help message.")
cli.m(longOpt: "metadata", argName: "filename", args: 1,
        "Metadata file. First and second columns should contain file name and sample id. " +
                "Header is mandatory and will be used to assign column names for metadata.")
cli.e(longOpt: "negative", "Will retain only clonotypes that lack specified V/D/J segments.")
cli.v(longOpt: "v-segments", argName: "v1,v2,...", args: 1, "A comma-separated list of Variable segments. " +
        "Incomplete names will act as wildcards.")
cli.d(longOpt: "d-segments", argName: "d1,d2,...", args: 1, "A comma-separated list of Diversity segments. " +
        "Incomplete names will act as wildcards.")
cli.j(longOpt: "j-segments", argName: "j1,j2,...", args: 1, "A comma-separated list of Joining segments. " +
        "Incomplete names will act as wildcards.")
cli._(longOpt: "save-freqs", "Preserve clonotype frequencies as in original sample. " +
        "By default, output sample(s) will be re-normalized to have a sum of clonotype frequencies equal to 1.")
cli.c(longOpt: "compress", "Compress output sample files.")

def opt = cli.parse(args)

if (opt == null)
    System.exit(2)

if (opt.h || opt.arguments().size() == 0) {
    cli.usage()
    System.exit(2)
}

// Check if metadata is provided

def metadataFileName = opt.m

if (metadataFileName ? opt.arguments().size() != 1 : opt.arguments().size() < 2) {
    if (metadataFileName)
        println "Only output prefix should be provided in case of -m"
    else
        println "At least 1 sample files should be provided if not using -m"
    cli.usage()
    System.exit(2)
}

// Remaining arguments

def outputFilePrefix = opt.arguments()[-1],
    saveFreqs = (boolean) opt.'save-freqs',
    compress = (boolean) opt.c,
    negative = (boolean) opt.e

def vFilter = opt.v ? new VFilter(((String) opt.v).split(",")) : BlankClonotypeFilter.INSTANCE,
    dFilter = opt.d ? new DFilter(((String) opt.d).split(",")) : BlankClonotypeFilter.INSTANCE,
    jFilter = opt.j ? new JFilter(((String) opt.j).split(",")) : BlankClonotypeFilter.INSTANCE

def scriptName = getClass().canonicalName.split("\\.")[-1]

//
// Batch load all samples (lazy)
//

println "[${new Date()} $scriptName] Reading sample(s)"

def sampleCollection = metadataFileName ?
        new SampleCollection((String) metadataFileName) :
        new SampleCollection(opt.arguments()[0..-2])

println "[${new Date()} $scriptName] ${sampleCollection.size()} sample(s) loaded"

//
// Iterate over samples & filter
//

def writer = new SampleWriter(compress, !saveFreqs)

def filter = new CompositeClonotypeFilter(negative, vFilter, dFilter, jFilter)
new File(formOutputPath(outputFilePrefix, "segfilter", "summary")).withPrintWriter { pw ->
    def header = "$MetadataTable.SAMPLE_ID_COLUMN\t" +
            sampleCollection.metadataTable.columnHeader + "\t" +
            ClonotypeFilter.ClonotypeFilterStats.HEADER

    pw.println(header)

    sampleCollection.eachWithIndex { sample, ind ->
        def sampleId = sample.sampleMetadata.sampleId
        println "[${new Date()} $scriptName] Filtering $sampleId.."

        def filteredSample = new Sample(sample, filter)

        // print output
        writer.writeConventional(filteredSample, outputFilePrefix)

        def stats = filter.getStatsAndFlush()

        pw.println([sampleId, sample.sampleMetadata, stats].join("\t"))
    }
}

sampleCollection.metadataTable.storeWithOutput(outputFilePrefix, compress,
        "segfilter:${negative ? "remove" : "keep"}:${opt.v ?: "."}:${opt.d ?: "."}:${opt.j ?: "."}")

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
