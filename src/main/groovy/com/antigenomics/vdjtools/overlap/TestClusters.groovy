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

import com.antigenomics.vdjtools.overlap.permutations.DiscreteFactorClusterStats

import static com.antigenomics.vdjtools.misc.ExecUtil.formOutputPath
import static com.antigenomics.vdjtools.misc.ExecUtil.toPlotPath
import static com.antigenomics.vdjtools.misc.RUtil.execute

def MEASURE_DEFAULT = "F", I_TYPE_DEFAULT = "aa"

def cli = new CliBuilder(usage: "TestClusters [options] input_prefix [output_prefix]\n" +
        "NOTE: input_prefix should be equal to output_prefix specified for" +
        "ClusterSamples execution, -i, -e parameters should also match.")

cli.h("display help message")
cli.e(longOpt: "measure", argName: "string", args: 1,
        "Distance measure to use, allowed values are ${OverlapMetric.allowedNames}. " +
                "[default = $MEASURE_DEFAULT]")
cli.i(longOpt: "intersect-type", argName: "string", args: 1,
        "Intersection rule, as used in CalcPairwiseDistances." +
                "Allowed values: $OverlapType.allowedNames. " +
                "Will use '$I_TYPE_DEFAULT' by default.")
cli.n(longOpt: "num-factor", "Treat factor as numeric")
cli._(longOpt: "plot-type", argName: "pdf|png", args: 1, "Plot output format [default=pdf]")

def opt = cli.parse(args)

if (opt == null)
    System.exit(2)

if (opt.h || opt.arguments().size() < 1) {
    cli.usage()
    System.exit(2)
}

def scriptName = getClass().canonicalName.split("\\.")[-1]

def iName = opt.i ?: I_TYPE_DEFAULT
def intersectionType = OverlapType.getByShortName(iName)

if (!intersectionType) {
    println "[ERROR] Bad overlap type specified ($iName). " +
            "Allowed values are: $OverlapType.allowedNames"
    System.exit(2)
}

intersectionType = intersectionType.shortName

def measureName = (opt.e ?: MEASURE_DEFAULT).toUpperCase(),
    numFactor = (boolean) opt.n,
    inputPrefix = opt.arguments()[0],
    mdsFileName = formOutputPath(inputPrefix, "mds", intersectionType, measureName),
    outputPrefix = opt.arguments().size() > 1 ? opt.arguments()[1] : inputPrefix,
    plotType = (opt.'plot-type' ?: "pdf").toString()

//
// Permutation testing
//

if (numFactor) {
    // todo: finish with distance correlation
} else {
    println "[${new Date()} $scriptName] Running permutation testing for factor ~ cluster dependence"
    def permsOutputPath = formOutputPath(outputPrefix, "perms", intersectionType, measureName)
    def summary = new DiscreteFactorClusterStats(mdsFileName).performPermutations(10000)
    if (summary) {
        DiscreteFactorClusterStats.writeSummary(summary, permsOutputPath)
        execute("cluster_permutations_plot.r",
                permsOutputPath,
                toPlotPath(permsOutputPath, plotType)
        )
        new File(permsOutputPath).delete()
    } else {
        println "[${new Date()} $scriptName] No way - less than 2 factor levels are present."
    }
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
