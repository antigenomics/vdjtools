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

import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.sample.metadata.MetadataTable

import static com.antigenomics.vdjtools.misc.ExecUtil.formOutputPath

def DEFAULT_AA_PROPERTIES = ["hydropathy", "charge", "polarity", "strength", "mjenergy", "count"].join(","),
    ALLOWED_AA_PROPERTIES = KnownAminoAcidProperties.INSTANCE.allowedNames.join(","),
    DEFAULT_REGIONS = ["CDR3-full", "VJ-junc", "V-germ", "J-germ", "CDR3-center-5"].join(","),
    ALLOWED_REGIONS = KnownCdr3Regions.INSTANCE.allowedNames.join(",")

def cli = new CliBuilder(usage: "CalcCdrAaStats [options] " +
        "[sample1 sample2 sample3 ... if -m is not specified] output_prefix")
cli.h("display help message")
cli.m(longOpt: "metadata", argName: "filename", args: 1,
        "Metadata file. First and second columns should contain file name and sample id. " +
                "Header is mandatory and will be used to assign column names for metadata.")
cli.w(longOpt: "weighted", "Weight all statistics by clonotype frequency.")
cli.n(longOpt: "normalize", "Normalize property value by dividing them by corresponding CDR3 region length.")
cli.a(longOpt: "aa-properties", argName: "property1,...", args: 1,
        "Comma-separated list of amino-acid properties to summarize. " +
                "Allowed values: $ALLOWED_AA_PROPERTIES. " +
                "[default = $DEFAULT_AA_PROPERTIES]")
cli.r(longOpt: "region-list", argName: "region1,...", args: 1,
        "List of CDR3 regions to analyze. " +
                "Allowed segments: $ALLOWED_REGIONS. " +
                "[default = $DEFAULT_REGIONS]")
/*cli.p(longOpt: "plot", "Plot amino acid property distributions for a specified list of segments.")
cli.f(longOpt: "factor", argName: "string", args: 1, "Metadata entry used to group samples in plot.")
cli._(longOpt: "plot-normalized", "Will normalize regions by the total number of AAs in them.")
cli._(longOpt: "plot-type", argName: "pdf|png", args: 1, "Plot output format [default=pdf]")
cli._(longOpt: "include-cfw", "Consider first and last AAs of CDR3, which are normally conserved C and F/W")*/


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
    normalize = (boolean) opt.n, weighted = (boolean) opt.w,
    regionNames = (opt.r ?: DEFAULT_REGIONS).split(",").collect { it.toLowerCase() },
    propertyNames = (opt.a ?: DEFAULT_AA_PROPERTIES).split(",").collect { it.toLowerCase() }

/*plot = (boolean) opt.p,
    plotType = (opt.'plot-type' ?: "pdf").toString(),
    includeCFW = (boolean) opt.'include-cfw'*/

def badRegionNames = regionNames.findAll { !KnownCdr3Regions.INSTANCE.allowedNames.contains(it) }

if (badRegionNames.size() > 0) {
    println "[ERROR] Unknown amino acid properties: ${badRegionNames.join(",")}. " +
            "Allowed values are $ALLOWED_REGIONS"
    System.exit(2)
}

def badProperties = propertyNames.findAll { !KnownAminoAcidProperties.INSTANCE.allowedNames.contains(it) }

if (badProperties.size() > 0) {
    println "[ERROR] Unknown amino acid properties: ${badProperties.join(",")}. " +
            "Allowed values are $DEFAULT_AA_PROPERTIES"
    System.exit(2)
}

def scriptName = getClass().canonicalName.split("\\.")[-1]

//
// Batch load all samples (lazy)
//

println "[${new Date()} $scriptName] Reading sample(s)"

def sampleCollection = metadataFileName ?
        new SampleCollection((String) metadataFileName) :
        new SampleCollection(opt.arguments()[0..-2])

println "[${new Date()} $scriptName] ${sampleCollection.size()} sample(s) prepared"

//
// Prepare summary calculators
//

def summarizerMap = new HashMap<String, AaPropertySummaryEvaluator>()

regionNames.each { regionName ->
    propertyNames.each { propertyName ->
        try {
            summarizerMap.put(regionName + "\t" + propertyName,
                    new AaPropertySummaryEvaluator(KnownAminoAcidProperties.INSTANCE.getByName(propertyName),
                            KnownCdr3Regions.INSTANCE.getByName(regionName),
                            normalize, weighted))
        } catch (IllegalArgumentException e) {
            // do nothing - incompatible region <> property pair
        }
    }
}

//
// Compute and write summary for each region<>property pair and each sample
//

def outputFileName = formOutputPath(outputFilePrefix, "cdr3aa", "stat",
        (weighted ? "wt" : "unwt"),
        (normalize ? "norm" : "unnorm"))

new File(outputFileName).withPrintWriter { pw ->
    def header = "$MetadataTable.SAMPLE_ID_COLUMN\t" +
            sampleCollection.metadataTable.columnHeader + "\t" +
            "region\tproperty\tmean\tq25\tmedian\tq75"

    pw.println(header)

    def sampleCounter = 0

    sampleCollection.each { Sample sample ->
        println "[${new Date()} $scriptName] Running $sample.sampleMetadata.sampleId"

        summarizerMap.each {
            println "[${new Date()} $scriptName] Summarizing for $it.key"

            def summary = it.value.compute(sample)

            pw.println([sample.sampleMetadata.sampleId, sample.sampleMetadata,
                        it.key,
                        summary.mean, summary.q25, summary.median, summary.q75].join("\t"))
        }

        println "[${new Date()} $scriptName] ${++sampleCounter} sample(s) processed"
    }
}

/*
if (plot) {
    RUtil.execute("cdr3aa_profile.r",
            outputFileName,
            toPlotPath(outputFileName, plotType),
            opt.f ? (sampleCollection.metadataTable.getColumnIndex(opt.f) + 2).toString() : "0",
            RUtil.logical(opt.'plot-normalized')
    )
}
*/

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
