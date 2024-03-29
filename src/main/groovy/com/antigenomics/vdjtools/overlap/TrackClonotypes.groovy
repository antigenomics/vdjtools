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

package com.antigenomics.vdjtools.overlap

import com.antigenomics.vdjtools.misc.Software
import com.antigenomics.vdjtools.io.SampleWriter
import com.antigenomics.vdjtools.join.JointSample
import com.antigenomics.vdjtools.join.OccurrenceJoinFilter
import com.antigenomics.vdjtools.join.SampleSpecificJoinFilter
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.sample.metadata.MetadataTable
import com.antigenomics.vdjtools.misc.RUtil

import static com.antigenomics.vdjtools.misc.ExecUtil.formOutputPath
import static com.antigenomics.vdjtools.misc.ExecUtil.toPlotPath

def I_TYPE_DEFAULT = "strict", TOP_DEFAULT = "100", TOP_MAX = 200
def cli = new CliBuilder(usage: "TrackClonotypes [options] " +
        "[sample1 sample2 sample3 ... if -m is not specified] output_prefix")
cli.h("display help message")

// General

cli.m(longOpt: "metadata", argName: "filename", args: 1,
        "Metadata file. First and second columns should contain file name and sample id. " +
                "Header is mandatory and will be used to assign column names for metadata." +
                "If column named 'time' is present, it will be used to specify time point sequence.")
cli.i(longOpt: "intersect-type", argName: "string", args: 1,
        "Intersection rule to apply. Allowed values: $OverlapType.allowedNames. " +
                "Will use '$I_TYPE_DEFAULT' by default.")
// Tracking and output

cli.x(longOpt: "track-sample", argName: "int", args: 1,
        "A zero-based index of time point to track. " +
                "Will otherwise consider all clonotypes that were detected in 2+ samples")
cli.t(longOpt: "top", args: 1, "Number of top clonotypes which will be provided in the collapsed joint table " +
        "and shown on the summary stacked area plot. " +
        "Values > $TOP_MAX are not allowed, as they would make the plot unreadable. [default = $TOP_DEFAULT]")
cli.c(longOpt: "compress", "Compress output sample files.")

// Plotting

cli.s(longOpt: "sequence", argName: "[t1,t2,t3,...]", args: 1,
        "[plotting] Time point sequence. Unused if -m is specified. " +
                "If not specified, either time values from metadata, " +
                "or sample indexes (as in command line) are used.")
cli.f(longOpt: "factor", argName: "string", args: 1,
        "[plotting] Column name, as in metadata. Factor to be treated as time variable. [default = \"time\"]")
cli.p(longOpt: "plot", "[plotting] Turns on plotting.")
cli._(longOpt: "plot-type", argName: "pdf|png", args: 1, "Plot output format [default=pdf]")

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

if (metadataFileName ? opt.arguments().size() != 1 : opt.arguments().size() < 4) {
    if (metadataFileName)
        println "Only output prefix should be provided in case of -m"
    else
        println "At least 3 sample files should be provided if not using -m"
    cli.usage()
    System.exit(2)
}

// Other parameters

def trackSample = (opt.x ?: "-1").toInteger(),
    plot = (boolean) opt.p, compress = (boolean) opt.c,
    timeFactor = opt.f ?: "time",
    outputPrefix = opt.arguments()[-1],
    plotType = (opt.'plot-type' ?: "pdf").toString()

def scriptName = getClass().canonicalName.split("\\.")[-1]

// Select overlap type

def iName = opt.i ?: I_TYPE_DEFAULT
def intersectionType = OverlapType.getByShortName(iName)

if (!intersectionType) {
    println "[ERROR] Bad overlap type specified ($iName). " +
            "Allowed values are: $OverlapType.allowedNames"
    System.exit(2)
}

// Define number of clonotypes to show explicitly

def top = (opt.t ?: TOP_DEFAULT).toInteger()

if (top > TOP_MAX) {
    println "[ERROR] Specified number of top clonotypes should not exceed $TOP_MAX"
    System.exit(2)
}

//
// Batch load samples
//

println "[${new Date()} $scriptName] Reading in all samples"

def sampleCollection = metadataFileName ?
        new SampleCollection((String) metadataFileName, Software.VDJtools, true, false) :
        new SampleCollection(opt.arguments()[0..-2], Software.VDJtools, true, false)

def metadataTable = sampleCollection.metadataTable

if (sampleCollection.size() < 3) {
    println "[ERROR] Metadata file should contain at least 3 samples"
    System.exit(2)
}

println "[${new Date()} $scriptName] ${sampleCollection.size()} samples loaded"

//
// Get time points & sort samples
//

// get time points (this will be unused if metadata already contains $timeFactor column)
def timePoints = opt.s ? (opt.s.toString())[1..-2].split(",").collect() : // [-10,0.5,1,20]
        (0..<sampleCollection.size()).collect { it.toString() } // value provided with argument or uniform

int n = timePoints.size()

if (n != sampleCollection.size()) {
    println "[ERROR] $n time points provided, " +
            "while sample collection contains ${sampleCollection.size()} samples"
    System.exit(2)
}

// sort
if (!metadataTable.containsColumn(timeFactor)) {
    // add corresponding column to metadata
    metadataTable.addColumn(timeFactor, timePoints as List<String>)
}

metadataTable.sort(timeFactor)

// Take sorted values from metadata
timePoints = metadataTable.getColumn(timeFactor).collect { it.asNumeric().toString() }

println "[${new Date()} $scriptName] Time points " +
        "${timePoints.join(",")} and corresponding samples ${metadataTable.sampleIterator.collect().join(",")} " +
        "will be processed"

//
// Join samples
//

println "[${new Date()} $scriptName] Joining samples by ${intersectionType}"

def jointSample = new JointSample(intersectionType, sampleCollection.collect() as Sample[],
        trackSample > -1 ? new SampleSpecificJoinFilter(trackSample) : new OccurrenceJoinFilter(2))

jointSample.computeAndCorrectSamplingPValues()

//
// Write output tables
//

println "[${new Date()} $scriptName] Writing tabular output"

def sampleWriter = new SampleWriter(compress)
sampleWriter.write(jointSample, formOutputPath(outputPrefix, "tracking", intersectionType.shortName, "table"))

def tableCollapsedOutputPath = formOutputPath(outputPrefix, "tracking", intersectionType.shortName, "table", "collapsed")
if (top >= 0) {
    sampleWriter = new SampleWriter(false)
    sampleWriter.write(jointSample, tableCollapsedOutputPath, top, true)
}

//
// Write summary output
//

println "[${new Date()} $scriptName] Writing output"

def summaryOutputPath = formOutputPath(outputPrefix, "tracking", intersectionType.shortName, "summary")
new File(summaryOutputPath).withPrintWriter { pw ->
    pw.println("1_$MetadataTable.SAMPLE_ID_COLUMN\t2_$MetadataTable.SAMPLE_ID_COLUMN\t" +
            "value\tmetric\t1_time\t2_time\t" +
            sampleCollection.metadataTable.columnHeader1 + "\t" +
            sampleCollection.metadataTable.columnHeader2)

    double[][] freqTable = new double[n][n],
               divTable = new double[n][n],
               countTable = new double[n][n]

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            freqTable[i][j] = jointSample.getIntersectionFreq(i, j)
            countTable[i][j] = jointSample.getIntersectionCount(i, j)
            divTable[i][j] = jointSample.getIntersectionDiv(i, j)
        }
    }

    for (int i = 0; i < jointSample.numberOfSamples; i++) {
        for (int j = 0; j < jointSample.numberOfSamples; j++) {
            def ids = [i, j].collect { jointSample.getSample(it).sampleMetadata.sampleId }.join("\t"),
                metadata = [i, j].collect { jointSample.getSample(it).sampleMetadata.toString() }.join("\t")
            def pointPair = [timePoints[i], timePoints[j]]

            if (i != j) {
                pw.println([ids, RUtil.asNumeric(freqTable[i][j]), "frequency", pointPair, metadata].flatten().join("\t"))
                pw.println([ids, RUtil.asNumeric(divTable[i][j]), "diversity", pointPair, metadata].flatten().join("\t"))
                pw.println([ids, RUtil.asNumeric(countTable[i][j]), "count", pointPair, metadata].flatten().join("\t"))
            } else {
                pw.println([ids, "NA", "frequency", pointPair, metadata].flatten().join("\t"))
                pw.println([ids, "NA", "diversity", pointPair, metadata].flatten().join("\t"))
                pw.println([ids, "NA", "count", pointPair, metadata].flatten().join("\t"))
            }
        }
    }
}

//
// Plotting via R
//

if (plot) {
    println "[${new Date()} $scriptName] Writing plots"

    // Plot all the heatmaps
    RUtil.execute("tracking_similarity_map.r",
            summaryOutputPath,
            toPlotPath(summaryOutputPath, plotType))

    // Plot a stack plot of top X clonotype abundances
    RUtil.execute("tracking_stackplot.r",
            timeFactor,
            timePoints.join(","),
            tableCollapsedOutputPath,
            formOutputPath(outputPrefix, "tracking", intersectionType.shortName, "stackplot", "pdf"))

    // Plot a "heatcourse" plot of top X clonotype abundances
    RUtil.execute("tracking_heatcourse.r",
            timeFactor,
            timePoints.join(","),
            tableCollapsedOutputPath,
            formOutputPath(outputPrefix, "tracking", intersectionType.shortName, "heatplot", "pdf"))
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
