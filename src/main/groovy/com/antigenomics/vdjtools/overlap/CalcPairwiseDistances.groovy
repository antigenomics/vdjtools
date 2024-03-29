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
import com.antigenomics.vdjtools.basic.SegmentUsage
import com.antigenomics.vdjtools.sample.SampleCollection

import static com.antigenomics.vdjtools.misc.ExecUtil.formOutputPath
import static com.antigenomics.vdjtools.misc.ExecUtil.toPlotPath
import static com.antigenomics.vdjtools.misc.RUtil.execute

def I_TYPE_DEFAULT = "aa"
def cli = new CliBuilder(usage: "CalcPairwiseDistances [options] " +
        "[sample1 sample2 sample3 ... if not -m] output_prefix")
cli.h("display help message")
cli.m(longOpt: "metadata", argName: "filename", args: 1,
        "Metadata file. First and second columns should contain file name and sample id. " +
                "Header is mandatory and will be used to assign column names for metadata.")
cli._(longOpt: "low-mem", "Will process all sample pairs sequentially, avoiding" +
        " loading all of them into memory. Slower but memory-efficient mode.")
cli.i(longOpt: "intersect-type", argName: "string", args: 1,
        "Intersection rule to apply. Allowed values: $OverlapType.allowedNames. " +
                "Will use '$I_TYPE_DEFAULT' by default.")
cli._(longOpt: "plot-type", argName: "pdf|png", args: 1, "Plot output format [default=pdf]")
cli.p(longOpt: "plot", "[plotting] Turns plotting on.")

def opt = cli.parse(args)

if (opt == null) {
    //cli.usage()
    System.exit(2)
}

if (opt.h) {
    cli.usage()
    System.exit(0)
}

// Check if metadata is provided

def metadataFileName = opt.m

if (metadataFileName ? opt.arguments().size() != 1 : opt.arguments().size() < 4) {
    if (metadataFileName)
        println "[ERROR] Only output prefix should be provided in case of -m"
    else
        println "[ERROR] At least 3 sample file names should be provided if not using -m"
    cli.usage()
    System.exit(2)
}

def outputPrefix = opt.arguments()[-1],
    lowMem = (boolean) opt.'low-mem',
    plot = (boolean) opt.p,
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

//
// Batch load all samples
//

println "[${new Date()} $scriptName] Reading samples"

boolean store, lazy
(store, lazy) = lowMem ? [false, true] : [true, false]

def sampleCollection = metadataFileName ?
        new SampleCollection((String) metadataFileName, Software.VDJtools, store, lazy) :
        new SampleCollection(opt.arguments()[0..-2], Software.VDJtools, store, lazy)

if (sampleCollection.size() < 3) {
    println "[ERROR] Metadata file should contain at least 3 samples"
    System.exit(2)
}

println "[${new Date()} $scriptName] ${sampleCollection.size()} samples loaded"

//
// Perform overlap for all specified overlap types
//

println "[${new Date()} $scriptName] Intersecting by $intersectionType"

SegmentUsage.VERBOSE = false
Overlap.VERBOSE = false
OverlapEvaluator.VERBOSE = false

def pairedIntersectionBatch = new PairwiseOverlap(sampleCollection, intersectionType)

println "[${new Date()} $scriptName] Writing results"

def outputFileName = formOutputPath(outputPrefix, "intersect", "batch", intersectionType.shortName)

new File(outputFileName).withPrintWriter { pw ->
    pw.println(pairedIntersectionBatch.header)
    pw.println(pairedIntersectionBatch.toString())
}

if (plot) {
    println "[${new Date()} $scriptName] Plotting"
    execute("pairwise_distance_plot.r", outputFileName, toPlotPath(outputFileName, plotType))
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
