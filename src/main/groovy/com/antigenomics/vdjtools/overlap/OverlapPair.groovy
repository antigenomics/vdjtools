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
import com.antigenomics.vdjtools.sample.SampleCollection

import static com.antigenomics.vdjtools.misc.ExecUtil.formOutputPath
import static com.antigenomics.vdjtools.misc.ExecUtil.toPlotPath
import static com.antigenomics.vdjtools.misc.RUtil.execute

def I_TYPE_DEFAULT = "strict", TOP_DEFAULT = "20", TOP_MAX = 100
def cli = new CliBuilder(usage: "OverlapPair [options] sample1 sample2 output_prefix")
cli.h("display help message")
cli.i(longOpt: "intersect-type", argName: "string", args: 1,
        "Intersection rule to apply. Allowed values: $OverlapType.allowedNames. " +
                "Will use '$I_TYPE_DEFAULT' by default.")
cli.t(longOpt: "top", args: 1, "Number of top clonotypes which will be provided in the collapsed joint table " +
        "and shown on the summary stacked area plot. " +
        "Values > $TOP_MAX are not allowed, as they would make the plot unreadable. [default = $TOP_DEFAULT]")
cli.p(longOpt: "plot", "Generate a scatterplot to characterize overlapping clonotypes. " +
        "Also generate abundance difference plot if -c option is specified. " +
        "(R installation with ggplot2, grid and gridExtra packages required).")
cli._(longOpt: "plot-type", argName: "pdf|png", args: 1, "Plot output format [default=pdf]")
cli.c(longOpt: "compress", "Compress output sample files.")
cli._(longOpt: "plot-area-v2", "Use alternative stacked area plot.")

def opt = cli.parse(args)

if (opt == null)
    System.exit(2)

if (opt.h || opt.arguments().size() < 3) {
    cli.usage()
    System.exit(2)
}

def sample1FileName = opt.arguments()[0], sample2FileName = opt.arguments()[1],
    outputPrefix = opt.arguments()[2],
    compress = (boolean) opt.c,
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
// Load samples
//

println "[${new Date()} $scriptName] Reading samples $sample1FileName and $sample2FileName"

def sampleCollection = new SampleCollection([sample1FileName, sample2FileName], Software.VDJtools, true, false)

//
// Perform an overlap by CDR3NT & V segment
//

println "[${new Date()} $scriptName] Intersecting"

def pairedIntersection = new Overlap(sampleCollection.listPairs()[0], intersectionType, true)
def jointSample = pairedIntersection.jointSample
jointSample.computeAndCorrectSamplingPValues()

//
// Generate and write output
//

println "[${new Date()} $scriptName] Writing output"

new File(formOutputPath(outputPrefix, "paired", intersectionType.shortName, "summary")).withPrintWriter { pw ->
    pw.println(pairedIntersection.header)
    pw.println(pairedIntersection.toString())
}


def sampleWriter = new SampleWriter(compress)
sampleWriter.write(jointSample, formOutputPath(outputPrefix, "paired", intersectionType.shortName, "table"))

def tableCollapsedOutputPath = formOutputPath(outputPrefix, "paired", intersectionType.shortName, "table", "collapsed")
if (top >= 0) {
    sampleWriter = new SampleWriter(false)
    sampleWriter.write(jointSample, tableCollapsedOutputPath, top, true)
}

if (opt.p) {
    println "[${new Date()} $scriptName] Plotting"

    def sample1 = sampleCollection[0],
        sample2 = sampleCollection[1]

    def xyFile = new File(outputPrefix + ".xy.txt")
    xyFile.withPrintWriter { pw ->
        pw.println("x\ty")
        jointSample.each { jointClone ->
            pw.println((0..1).collect { sampleIndex ->
                jointClone.getFreq(sampleIndex)
            }.join("\t"))
        }
    }
    xyFile.deleteOnExit()

    def xxFile = new File(outputPrefix + ".xx.txt")
    xxFile.withPrintWriter { pw ->
        pw.println("xx")
        sample1.each { pw.println(it.freq) }
    }
    xxFile.deleteOnExit()

    def yyFile = new File(outputPrefix + ".yy.txt")
    yyFile.withPrintWriter { pw ->
        pw.println("yy")
        sample2.each { pw.println(it.freq) }
    }
    yyFile.deleteOnExit()

    execute("intersect_pair_scatter.r", sample1.sampleMetadata.sampleId, sample2.sampleMetadata.sampleId,
            outputPrefix + ".xy.txt", outputPrefix + ".xx.txt", outputPrefix + ".yy.txt",
            formOutputPath(outputPrefix, intersectionType.shortName, "paired", "scatter", plotType))

    def areaPlot = opt.'plot-area-v2' ? "intersect_pair_area_v2.r" : "intersect_pair_area.r"
    execute(areaPlot, sample1.sampleMetadata.sampleId, sample2.sampleMetadata.sampleId,
            tableCollapsedOutputPath, toPlotPath(tableCollapsedOutputPath, plotType))
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
