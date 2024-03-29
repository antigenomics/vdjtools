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

package com.antigenomics.vdjtools.annotate

import com.antigenomics.vdjtools.io.SampleWriter
import com.antigenomics.vdjtools.sample.SampleCollection

def DEFAULT_ANNOTATORS = ["cdr3length", "ndnsize", "insertsize",
                          "hydropathy", "charge", "polarity", "strength", "cdr3contact"].join(","),
    ALLOWED_ANNOTATORS = KnownAnnotators.INSTANCE.allowedNames.join(",")

def cli = new CliBuilder(usage: "Annotate [options] " +
        "[sample1 sample2 ... if not -m] output_prefix")
cli.h("display help message")
cli.m(longOpt: "metadata", argName: "filename", args: 1,
        "Metadata file. First and second columns should contain file name and sample id. " +
                "Header is mandatory and will be used to assign column names for metadata.")
cli.a(longOpt: "annotators", argName: "param1,param2,...", args: 1,
        "Comma-separated list annotator names. The can contain basic annotators (e.g. CDR3 insert size) and " +
                "annotators computing sum of CDR3 amino acid property values. " +
                "Amino acid annotator will report average property values for each CDR3 (divide them by " +
                "CDR3 length) if suffix '${AAPropertyAnnotator.NORMALIZED_SUFFIX}' is appended to its name. " +
                "Allowed values: $ALLOWED_ANNOTATORS. " +
                "[default=$DEFAULT_ANNOTATORS]")
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

if (metadataFileName ? opt.arguments().size() != 1 : opt.arguments().size() < 2) {
    if (metadataFileName)
        println "Only output prefix should be provided in case of -m"
    else
        println "At least 1 sample files should be provided if not using -m"
    cli.usage()
    System.exit(2)
}

def outputFilePrefix = opt.arguments()[-1],
    annotatorNames = (opt.a ?: DEFAULT_ANNOTATORS).split(",").collect { it.toLowerCase() },
    compress = (boolean) opt.c

def scriptName = getClass().canonicalName.split("\\.")[-1]

def badAnnotatorNames = annotatorNames.findAll { !KnownAnnotators.INSTANCE.allowedNames.contains(it) }

if (!badAnnotatorNames.empty) {
    println "The following annotator names are invalid: $badAnnotatorNames. Allowed values: $ALLOWED_ANNOTATORS."
    System.exit(2)
}

//
// Batch load samples
//

println "[${new Date()} $scriptName] Reading samples"

def sampleCollection = metadataFileName ?
        new SampleCollection((String) metadataFileName) :
        new SampleCollection(opt.arguments()[0..-2])

println "[${new Date()} $scriptName] ${sampleCollection.size()} samples prepared"

//
// Create annotators
//

def sampleAnnotator = new SampleAnnotator(annotatorNames.collect { KnownAnnotators.INSTANCE.getByName(it) })

//
// Iterate over samples and annotate
//
def sw = new SampleWriter(compress)

sampleCollection.eachWithIndex { sample, ind ->
    def sampleId = sample.sampleMetadata.sampleId
    println "[${new Date()} $scriptName] Annotating $sampleId.."

    sampleAnnotator.annotate(sample)

    // print output
    sw.writeConventional(sample, outputFilePrefix)
}

sampleCollection.metadataTable.storeWithOutput(outputFilePrefix, compress,
        "annot:${annotatorNames.join(",")}")

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
