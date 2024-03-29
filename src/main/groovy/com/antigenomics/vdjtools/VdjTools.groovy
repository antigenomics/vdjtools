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


package com.antigenomics.vdjtools

import com.antigenomics.vdjtools.annotate.Annotate
import com.antigenomics.vdjtools.annotate.CalcCdrAaStats
import com.antigenomics.vdjtools.annotate.CalcDegreeStats
import com.antigenomics.vdjtools.annotate.SegmentsToFamilies
import com.antigenomics.vdjtools.basic.*
import com.antigenomics.vdjtools.diversity.CalcDiversityStats
import com.antigenomics.vdjtools.diversity.PlotQuantileStats
import com.antigenomics.vdjtools.diversity.RarefactionPlot
import com.antigenomics.vdjtools.operate.JoinSamples
import com.antigenomics.vdjtools.operate.PoolSamples
import com.antigenomics.vdjtools.overlap.*
import com.antigenomics.vdjtools.preprocess.*
import com.antigenomics.vdjtools.misc.*

import java.util.jar.JarFile

def getVersion = {
    (getClass().classLoader.findResource(JarFile.MANIFEST_NAME).text =~
            /Implementation-Version: (.+)/)[0][1]
}

def printHelp = {
    def version = getVersion()

    println "VDJtools V$version"
    println ""
    println "Run as \$java -jar vdjtools-${version}.jar ROUTINE_NAME arguments"
    println ""
    println "[Basic]"
    println "CalcBasicStats"
    println "CalcSpectratype"
    println "CalcSegmentUsage"
    println "PlotFancySpectratype"
    println "PlotSpectratypeV"
    println "PlotFancyVJUsage"
    println ""
    println "[Diversity]"
    println "CalcDiversityStats"
    println "RarefactionPlot"
    println "PlotQuantileStats"
    println ""
    println "[Overlap]"
    println "OverlapPair"
    println "CalcPairwiseDistances"
    println "ClusterSamples"
    println "TestClusters"
    println "TrackClonotypes"
    println ""
    println "[Preprocessing]"
    println "ApplySampleAsFilter"
    println "FilterNonFunctional"
    println "FilterByFrequency"
    println "DownSample"
    println "Decontaminate"
    println "FilterBySegment"
    println "SelectTop"
    println ""
    println "[Operation]"
    println "PoolSamples"
    println "JoinSamples"
    println "(Enrichment) -> Deprecated"
    println ""
    println "[Annotation]"
    println "(ScanDatabase) -> moved to VDJdb since 1.0.5, please visit vdjdb.cdr3.net"
    println "CalcCdrAaStats"
    println "CalcDegreeStats"
    println "Annotate"
    println "SegmentsToFamilies"
    println ""
    println "[Util]"
    println "FilterMetadata"
    println "SplitMetadata"
    println "Convert"
    println "RInstall"
    println ""
    println "* Run with 'discard_scripts' option prior to ROUTINE_NAME to clean up R scripts upon execution"
}

def getScript = { String scriptName ->
    switch (scriptName.toUpperCase()) {
        case "CALCBASICSTATS":
            return new CalcBasicStats()
        case "CALCSPECTRATYPE":
            return new CalcSpectratype()
        case "CALCSEGMENTUSAGE":
            return new CalcSegmentUsage()
        case "PLOTFANCYSPECTRATYPE":
            return new PlotFancySpectratype()
        case "PLOTSPECTRATYPEV":
            return new PlotSpectratypeV()
        case "PLOTFANCYVJUSAGE":
            return new PlotFancyVJUsage()

        case "CALCDIVERSITYSTATS":
            return new CalcDiversityStats()
        case "RAREFACTIONPLOT":
            return new RarefactionPlot()
        case "PLOTQUANTILESTATS":
            return new PlotQuantileStats()

        case "OVERLAPPAIR":
            return new OverlapPair()
        case "CALCPAIRWISEDISTANCES":
            return new CalcPairwiseDistances()
        case "CLUSTERSAMPLES":
            return new ClusterSamples()
        case "TESTCLUSTERS":
            return new TestClusters()
        case "TRACKCLONOTYPES":
            return new TrackClonotypes()

        case "APPLYSAMPLEASFILTER":
            return new ApplySampleAsFilter()
        case "DECONTAMINATE":
            return new Decontaminate()
        case "FILTERBYSEGMENT":
            return new FilterBySegment()
        case "FILTERBYFREQUENCY":
            return new FilterByFrequency()
        case "FILTERNONFUNCTIONAL":
            return new FilterNonFunctional()
        case "DOWNSAMPLE":
            return new DownSample()
        case "SELECTTOP":
            return new SelectTop()
        case "CORRECT":
            return new Correct()

        case "POOLSAMPLES":
            return new PoolSamples()
        case "JOINSAMPLES":
            return new JoinSamples()
        case "ENRICHMENT":
            println "Deprecated"
            System.exit(0)
            break

        case "CALCCDRAASTATS":
            return new CalcCdrAaStats()
        case "CALCDEGREESTATS":
            return new CalcDegreeStats()
        case "ANNOTATE":
            return new Annotate()
        case "SEGMENTSTOFAMILIES":
            return new SegmentsToFamilies()
        case "SCANDATABASE":
            println "Moved to VDJdb since 1.0.5, please visit vdjdb.cdr3.net"
            System.exit(0)
            break

        case "FILTERMETADATA":
            return new FilterMetadata()
        case "SPLITMETADATA":
            return new SplitMetadata()
        case "CONVERT":
            return new Convert()
        case "RINSTALL":
            return new RInstall()
        case "-H":
        case "H":
        case "-HELP":
        case "HELP":
        case "":
            printHelp()
            println ""
            System.exit(0)
            break

        case "DUMMY":
            System.exit(0)
            break

        default:
            printHelp()
            println ""
            println "Unknown routine name $scriptName"
            System.exit(0)
    }
}

if (args.length == 0) {
    printHelp()
} else {
    if (args[0].toLowerCase() == "discard_scripts") {
        if (args.length == 1) {
            printHelp()
            System.exit(0)
        }
        RUtil.REMOVE_R_SCRIPTS = true
        args = args[1..-1]
    }
    def script = getScript(args[0])
    try {
        ExecUtil.run(script, args.size() > 1 ? args[1..-1] : [""])
    } catch (Exception e) {
        def version = getVersion()
        println "[ERROR] ${e.toString()}, see _vdjtools_error.log for details"
        new File("_vdjtools_error.log").withWriterAppend { writer ->
            writer.println("[${new Date()} BEGIN]")
            writer.println("[Script]")
            writer.println(args[0])
            writer.println("[CommandLine]")
            writer.println("executing vdjtools-${version}.jar ${args.join(" ")}")
            writer.println("[Message]")
            writer.println(e.toString())
            writer.println("[StackTrace-Short]")
            writer.println(e.stackTrace.findAll { it.toString().contains("com.antigenomics.vdjtools") }.join("\n"))
            writer.println("[StackTrace-Full]")
            e.printStackTrace(new PrintWriter(writer))
            writer.println("[END]")
        }
        System.exit(1)
    }
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
