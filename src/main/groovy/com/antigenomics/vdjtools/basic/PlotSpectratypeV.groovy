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

def TOP_DEFAULT = "12", TOP_MAX = 12
def cli = new CliBuilder(usage: "PlotSpectratypeV [options] input_name output_prefix")
cli.h("display help message")
cli.t(longOpt: "top", args: 1, "Number of top V segments to present on the histogram. " +
        "Values > $TOP_MAX are not allowed, as they would make the plot unreadable. [default = $TOP_DEFAULT]")
cli.u(longOpt: "unweighted", "Will count each clonotype only once, apart from conventional frequency-weighted histogram.")
cli._(longOpt: "plot-type", argName: "pdf|png", args: 1, "Plot output format [default=pdf]")

def opt = cli.parse(args)

if (opt == null)
    System.exit(2)

if (opt.h || opt.arguments().size() != 2) {
    cli.usage()
    System.exit(2)
}

def outputFilePrefix = opt.arguments()[1],
    top = (opt.t ?: TOP_DEFAULT).toInteger(),
    unweighted = (boolean) opt.u,
    plotType = (opt.'plot-type' ?: "pdf").toString()

if (top > TOP_MAX) {
    println "[ERROR] Specified number of top V segments should not exceed 20"
    System.exit(2)
}

def scriptName = getClass().canonicalName.split("\\.")[-1]

//
// Read the sample
//

println "[${new Date()} $scriptName] Reading sample"

def sampleCollection = new SampleCollection([opt.arguments()[0]])

def sample = sampleCollection[0]

// Calculate spectratype

def spectratypeV = new SpectratypeV(false, unweighted)

spectratypeV.addAll(sample)

def collapsedSpectratypes = spectratypeV.collapse(top)

// Prepare output table

def spectraMatrix = new double[spectratypeV.span][top + 1]

def otherHistogram = collapsedSpectratypes["other"].getHistogram(false)
for (int i = 0; i < spectratypeV.span; i++) {
    spectraMatrix[i][0] = otherHistogram[i]
}

collapsedSpectratypes.findAll { it.key != "other" }.eachWithIndex { it, ind ->
    def histogram = it.value.getHistogram(false)
    for (int i = 0; i < spectratypeV.span; i++) {
        spectraMatrix[i][top - ind] = histogram[i]
    }
}

def table = "Len\tOther\t" + collapsedSpectratypes.findAll { it.key != "other" }.collect { it.key }.reverse().join("\t")
for (int i = 0; i < spectratypeV.span; i++) {
    table += "\n" + spectratypeV.lengths[i] + "\t" + spectraMatrix[i].collect().join("\t")
}

// Output

println "[${new Date()} $scriptName] Writing output and plotting data"

def outputFileName = formOutputPath(outputFilePrefix, "spectraV", (unweighted ? "unwt" : "wt"))

new File(outputFileName).withPrintWriter { pw ->
    pw.println(table)
}

RUtil.execute("fancy_spectratype.r",
        outputFileName, toPlotPath(outputFileName, plotType), "Variable segment", RUtil.logical(false)
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
