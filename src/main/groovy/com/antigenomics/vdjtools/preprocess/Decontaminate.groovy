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

import com.antigenomics.vdjtools.misc.Software
import com.antigenomics.vdjtools.io.SampleWriter
import com.antigenomics.vdjtools.pool.RatioFilter
import com.antigenomics.vdjtools.sample.ClonotypeFilter
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.sample.metadata.MetadataTable

import static com.antigenomics.vdjtools.misc.ExecUtil.formOutputPath

def DEFAULT_CONT_RATIO = "20"

def cli = new CliBuilder(usage: "Decontaminate [options] " +
        "[sample1 sample2 ... if not -m] output_prefix")
cli.h("display help message")
cli.m(longOpt: "metadata", argName: "filename", args: 1,
        "Metadata file. First and second columns should contain file name and sample id. " +
                "Header is mandatory and will be used to assign column names for metadata.")
cli._(longOpt: "read-based", "Remove contamination based on read ratios, not fractions. " +
        "Useful when dealing with contaminations on the same sequencing lane, not that " +
        "good for FACS-related contaminations.")
cli.r(longOpt: "ratio", argName: "double", args: 1,
        "Parent-to-child clonotype frequency ratio for contamination filtering [default = $DEFAULT_CONT_RATIO]")
cli._(longOpt: "low-mem", "Will process all sample pairs sequentially, avoiding" +
        " loading all of them into memory. Slower but memory-efficient mode. [UNUSED]")
cli.c(longOpt: "compress", "Compress output sample files.")

def opt = cli.parse(args)

if (opt == null) {
    System.exit(2)
}

if (opt.h) {
    cli.usage()
    System.exit(2)
}

// Check if metadata is provided

def metadataFileName = opt.m

if (metadataFileName ? opt.arguments().size() != 1 : opt.arguments().size() < 3) {
    if (metadataFileName)
        println "Output prefix should be provided in case of -m"
    else
        println "At least 2 sample files should be provided if not using -m"
    cli.usage()
    System.exit(2)
}

def outputFilePrefix = opt.arguments()[-1],
// todo: implement low-mem version
// this will lead eventually to the problem of removing Clonotype ref from dynamic clonotype
    lowMem = (boolean) opt.'low-mem',
    ratio = (opt.r ?: DEFAULT_CONT_RATIO).toDouble(),
    readBased = (boolean) opt.'read-based',
    compress = (boolean) opt.c

def scriptName = getClass().canonicalName.split("\\.")[-1]

//
// Batch load samples
//

println "[${new Date()} $scriptName] Reading samples"

def sampleCollection = metadataFileName ?
        new SampleCollection((String) metadataFileName, Software.VDJtools, false, true) :
        new SampleCollection(opt.arguments()[0..-2], Software.VDJtools, false, true)

println "[${new Date()} $scriptName] ${sampleCollection.size()} samples prepared"

//
// Vamos a generar un sample aunado
//

println "[${new Date()} $scriptName] Creating sample pool for filtering"

def ratioFilter = readBased ? new RatioFilter(sampleCollection, ratio) :
        new RatioFilter(sampleCollection, ratio)

//
// Go through all sample once more and perform freq-based filtering
//
def sw = new SampleWriter(compress, true)

new File(formOutputPath(outputFilePrefix, "dec", "summary")).withPrintWriter { pw ->
    def header = "$MetadataTable.SAMPLE_ID_COLUMN\t" +
            sampleCollection.metadataTable.columnHeader + "\t" +
            ClonotypeFilter.ClonotypeFilterStats.HEADER

    pw.println(header)

    sampleCollection.each { sample ->
        def sampleId = sample.sampleMetadata.sampleId
        println "[${new Date()} $scriptName] Filtering $sampleId"
        def newSample = new Sample(sample, ratioFilter)

        sw.writeConventional(newSample, outputFilePrefix)

        def stats = ratioFilter.getStatsAndFlush()

        pw.println([sampleId, sample.sampleMetadata, stats].join("\t"))
    }
}

sampleCollection.metadataTable.storeWithOutput(outputFilePrefix, compress, "dec:$ratio")

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
