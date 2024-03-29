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
import com.antigenomics.vdjtools.sample.CountFilter
import com.antigenomics.vdjtools.sample.*
import com.antigenomics.vdjtools.sample.metadata.MetadataTable

import static com.antigenomics.vdjtools.misc.ExecUtil.formOutputPath

def DEFAULT_FREQ_THRESHOLD = "0",
    DEFAULT_QUANTILE_THRESHOLD = "1",
    DEFAULT_COUNT_THRESHOLD = "0"

def cli = new CliBuilder(usage: "FilterByFrequency [options] " +
        "[sample1 sample2 sample3 ... if -m is not specified] output_prefix")
cli.h("display help message")
cli.m(longOpt: "metadata", argName: "filename", args: 1,
        "Metadata file. First and second columns should contain file name and sample id. " +
                "Header is mandatory and will be used to assign column names for metadata.")
cli.f(longOpt: "freq-threshold", argName: "double", args: 1,
        "Clonotype frequency threshold. Set it to 0 to disable. " +
                "[default = $DEFAULT_FREQ_THRESHOLD]")
cli.x(longOpt: "count-threshold", argName: "int", args: 1,
        "Clonotype count (number of reads) threshold. Set it to 0 to disable. [default = $DEFAULT_COUNT_THRESHOLD]")
cli.q(longOpt: "quantile-threshold", argName: "double", args: 1,
        "Quantile threshold. Will retain a set of top N clonotypes " +
                "so that their total frequency is equal or less to the specified threshold, " +
                "e.g. '0.25' will select top 25% clonotypes by frequency. " +
                "Set it to 1 to disable." +
                "[default = $DEFAULT_QUANTILE_THRESHOLD]")
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
    freqThreshold = (opt.f ?: DEFAULT_FREQ_THRESHOLD).toDouble(),
    countThreshold = (opt.x ?: DEFAULT_COUNT_THRESHOLD).toInteger(),
    quantileThreshold = (opt.q ?: DEFAULT_QUANTILE_THRESHOLD).toDouble(),
    saveFreqs = (boolean) opt.'save-freqs',
    compress = (boolean) opt.c

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

def filter = new CompositeClonotypeFilter(
        new CountFilter(countThreshold),
        new FrequencyFilter(freqThreshold),
        new QuantileFilter(quantileThreshold)
)

new File(formOutputPath(outputFilePrefix, "freqfilter", "summary")).withPrintWriter { pw ->
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
        "freqfilter:$freqThreshold:$quantileThreshold:$countThreshold")

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
