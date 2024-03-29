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

package com.antigenomics.vdjtools.basic

import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.misc.RUtil

import static com.antigenomics.vdjtools.misc.ExecUtil.formOutputPath
import static com.antigenomics.vdjtools.misc.ExecUtil.toPlotPath

def cli = new CliBuilder(usage: "PlotFancyVJUsage [options] input_name output_prefix")
cli.h("display help message")
cli.u(longOpt: "unweighted", "Will count each clonotype only once, apart from conventional frequency-weighted histogram.")
cli._(longOpt: "plot-type", argName: "pdf|png", args: 1, "Plot output format [default=pdf]")

def opt = cli.parse(args)

if (opt == null)
    System.exit(2)

if (opt.h || opt.arguments().size() != 2) {
    cli.usage()
    System.exit(2)
}

def unweighted = (boolean) opt.u,
    outputFilePrefix = opt.arguments()[1],
    plotType = (opt.'plot-type' ?: "pdf").toString()

def scriptName = getClass().canonicalName.split("\\.")[-1]

//
// Read the sample
//

println "[${new Date()} $scriptName] Reading sample"

def sampleCollection = new SampleCollection([opt.arguments()[0]])

def sampleId = sampleCollection.metadataTable.getRow(0).sampleId

// Calculate segment usage
def segmentUsage = new SegmentUsage(sampleCollection, unweighted)

// Output and plotting
println "[${new Date()} $scriptName] Writing output"

def outputFileName = formOutputPath(outputFilePrefix, "fancyvj", (unweighted ? "unwt" : "wt"))

new File(outputFileName).withPrintWriter { pw ->
    pw.println(".\t" + segmentUsage.vUsageHeader().collect().join("\t"))
    def vjMatrix = segmentUsage.vjUsageMatrix(sampleId)
    vjMatrix.eachWithIndex { double[] vVectorByJ, int i ->
        pw.println(segmentUsage.jUsageHeader()[i] + "\t" + vVectorByJ.collect().join("\t"))
    }
}

println "[${new Date()} $scriptName] Plotting data (be patient, complex graphics)"

RUtil.execute("vj_pairing_plot.r",
        outputFileName, toPlotPath(outputFileName, plotType)
)


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
