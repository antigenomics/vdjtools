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

import com.antigenomics.vdjtools.io.SampleFileConnection
import com.antigenomics.vdjtools.io.SampleWriter
import com.antigenomics.vdjtools.overlap.OverlapType
import com.antigenomics.vdjtools.sample.ClonotypeFilter
import com.antigenomics.vdjtools.sample.IntersectionClonotypeFilter
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.sample.metadata.MetadataTable

import static com.antigenomics.vdjtools.misc.ExecUtil.formOutputPath

def scriptName = getClass().canonicalName.split("\\.")[-1]

def I_TYPE_DEFAULT = "strict"
def cli = new CliBuilder(usage: "ApplySampleAsFilter [options] " +
        "[sample1 sample2 ... if not -m] filter_sample output_prefix")
cli.h("display help message")
cli.m(longOpt: "metadata", argName: "filename", args: 1,
        "Metadata file. First and second columns should contain file name and sample id. " +
                "Header is mandatory and will be used to assign column names for metadata.")
cli.i(longOpt: "intersect-type", argName: "string", args: 1,
        "Intersection type to apply. " +
                "Allowed values: $OverlapType.allowedNames. " +
                "Will use '$I_TYPE_DEFAULT' by default.")
cli.e(longOpt: "negative", "Will report clonotypes that are not present in filter_sample. " +
        "The default action is to retain only them.")
cli._(longOpt: "save-freqs", "Preserve clonotype frequencies as in original sample. " +
        "By default, output sample(s) will be re-normalized to have a sum of clonotype frequencies equal to 1.")
cli.c(longOpt: "compress", "Compress output sample files.")

def opt = cli.parse(args)

if (opt == null) {
    //cli.usage()
    System.exit(2)
}

if (opt.h) {
    cli.usage()
    System.exit(2)
}

// Check if enough arguments are provided

def metadataFileName = opt.m

if (metadataFileName ? opt.arguments().size() != 2 : opt.arguments().size() < 3) {
    if (metadataFileName)
        println "Output prefix and filter sample should be provided in case of -m"
    else
        println "At least 1 sample, filter sample and output path should be provided if not using -m"
    cli.usage()
    System.exit(2)
}

// IO stuff

def filterFileName = opt.arguments()[-2],
    saveFreqs = (boolean) opt.'save-freqs',
    compress = (boolean) opt.c,
    outputFilePrefix = opt.arguments()[-1]

// Parameters

def intersectionType = OverlapType.getByShortName((opt.i ?: I_TYPE_DEFAULT)),
    negative = (boolean) opt.e

if (!intersectionType) {
    println "[ERROR] Bad overlap type specified ($opt.i). " +
            "Allowed values are: $OverlapType.allowedNames"
    System.exit(2)
}

//
// Load samples
//

println "[${new Date()} $scriptName] Reading input samples & filter sample"

def sampleCollection = metadataFileName ?
        new SampleCollection((String) metadataFileName) :
        new SampleCollection(opt.arguments()[0..-3])

println "[${new Date()} $scriptName] Loading filter sample"

def filterSample = SampleFileConnection.load(filterFileName)

def clonotypeFilter = new IntersectionClonotypeFilter(intersectionType, filterSample, negative)

//
// Filter samples
//

println "[${new Date()} $scriptName] Filtering (${negative ? "negative" : "positive"}) and writing output"

def sw = new SampleWriter(compress, !saveFreqs)

new File(formOutputPath(outputFilePrefix, "asaf", "summary")).withPrintWriter { pw ->
    def header = "$MetadataTable.SAMPLE_ID_COLUMN\t" +
            sampleCollection.metadataTable.columnHeader + "\t" +
            ClonotypeFilter.ClonotypeFilterStats.HEADER

    pw.println(header)

    sampleCollection.each { Sample sample ->
        // Filter
        def sampleId = sample.sampleMetadata.sampleId

        println "[${new Date()} $scriptName] Filtering $sampleId sample."
        def filteredSample = new Sample(sample, clonotypeFilter)

        // print filter stats
        def stats = clonotypeFilter.getStatsAndFlush()
        pw.println([sampleId, sample.sampleMetadata, stats].join("\t"))

        // print output
        sw.writeConventional(filteredSample, outputFilePrefix)
    }
}

sampleCollection.metadataTable.storeWithOutput(outputFilePrefix, compress,
        "asaf:$filterSample.sampleMetadata.sampleId:${negative ? "-" : "+"}:$intersectionType.shortName")

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
