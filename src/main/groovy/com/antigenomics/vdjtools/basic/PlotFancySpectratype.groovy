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

def TOP_DEFAULT = "10", TOP_MAX = 20
def cli = new CliBuilder(usage: "PlotFancySpectratype [options] input_name output_prefix")
cli.h("display help message")
cli.t(longOpt: "top", args: 1, "Number of top clonotypes to present on the histogram. " +
        "Values > $TOP_MAX are not allowed, as they would make the plot legend unreadable. [default = $TOP_DEFAULT]")
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
    plotType = (opt.'plot-type' ?: "pdf").toString()

if (top > TOP_MAX) {
    println "[ERROR] Specified number of top clonotypes should not exceed $TOP_MAX"
    System.exit(2)
}

def scriptName = getClass().canonicalName.split("\\.")[-1]

//
// Read the sample
//

println "[${new Date()} $scriptName] Reading sample"

def sampleCollection = new SampleCollection([opt.arguments()[0]])

def sample = sampleCollection[0]

top = Math.min(top, sample.diversity)

// Calculate spectratype

def fancySpectratype = new FancySpectratype(sample, top)

def table = "Len\tOther\t" +
        fancySpectratype.topClonotypes.reverse().collect { it.cdr3aa }.join("\t")
for (int i = 0; i < fancySpectratype.spectratype.span; i++) {
    table += "\n" + fancySpectratype.spectratype.lengths[i] + "\t" +
            fancySpectratype.spectraMatrix[i].collect().join("\t")
}

// Output

println "[${new Date()} $scriptName] Writing output and plotting data"

def outputFileName = formOutputPath(outputFilePrefix, "fancyspectra")

new File(outputFileName).withPrintWriter { pw ->
    pw.println(table)
}

RUtil.execute("fancy_spectratype.r",
        outputFileName, toPlotPath(outputFileName, plotType), "Clonotype", "TRUE"
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
