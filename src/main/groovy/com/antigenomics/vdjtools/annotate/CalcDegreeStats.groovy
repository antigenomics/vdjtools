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

import com.antigenomics.vdjtools.graph.DegreeStatisticsCalculator
import com.antigenomics.vdjtools.graph.DummyClonotypeGroupingFactory
import com.antigenomics.vdjtools.graph.VClonotypeGroupingFactory
import com.antigenomics.vdjtools.graph.VJClonotypeGroupingFactory
import com.antigenomics.vdjtools.graph.VJLClonotypeGroupingFactory
import com.antigenomics.vdjtools.io.SampleFileConnection
import com.antigenomics.vdjtools.io.SampleWriter
import com.antigenomics.vdjtools.overlap.OverlapType
import com.antigenomics.vdjtools.pool.PooledSample
import com.antigenomics.vdjtools.pool.SampleAggregator
import com.antigenomics.vdjtools.pool.StoringClonotypeAggregatorFactory
import com.antigenomics.vdjtools.sample.SampleCollection


def DEFAULT_SEARCH_SCOPE = "1,0,1"

def cli = new CliBuilder(usage: "CalcDegreeStats [options] " +
        "[sample1 sample2 ... if not -m] output_prefix")
cli.h("display help message")
cli.m(longOpt: "metadata", argName: "filename", args: 1,
        "Metadata file. First and second columns should contain file name and sample id. " +
                "Header is mandatory and will be used to assign column names for metadata.")
cli.b(longOpt: "background", argName: "filename", args: 1,
        "Background (control) sample file to compute random graph statistics. In case not provided," +
                "will pool input samples and use them as control.")
cli.o(longOpt: "search-scope", argName: "s,id,t", args: 1,
        "Search scope: number of substitutions (s), indels (id) and total number of mismatches (t) allowed. " +
                "Default is $DEFAULT_SEARCH_SCOPE")
cli.g(longOpt: "grouping", argName: "type", args: 1,
        "Primary grouping type, limits the scope of possible clonotype comparisons when computing clonotype degree. " +
                "Group counts are used in primary statistical test (p.value in output)." +
                "Allowed values: 'dummy' (no grouping, default), " +
                "'v' (same V), 'vj' (same V and J) and 'vjl' (same V, J and CDR3 length).")
cli.g2(longOpt: "grouping2", argName: "type", args: 1,
        "Secondary grouping type, doesn't limit possible clonotype comparisons. " +
                "Group counts are used in secondary statistical test (p.value.2 in output). " +
                "Allowed values: 'dummy', 'v', 'vj' and 'vjl'. " +
                "By default will select 'vjl' if no indels are allowed and " +
                "'vj' otherwise.")
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

// General options

def outputFilePrefix = opt.arguments()[-1],
    backgroundSample = opt.b,
    compress = (boolean) opt.c,
    optSearchScope = (opt.o ?: DEFAULT_SEARCH_SCOPE).split(",")

// Search scope

if (optSearchScope.length != 3 || optSearchScope.any { !it.isInteger() || it.toInteger() < 0 }) {
    println "[ERROR] Bad search scope $optSearchScope"
    System.exit(2)
}
def searchScope = optSearchScope.collect { it.toInteger() } as int[]

// Grouping

def grouping = (opt.g ?: "dummy").toLowerCase(),
    grouping2 = (opt.g2 ?: (searchScope[1] == 0 ? "vjl" : "vj")).toLowerCase()

def getGroupingFactory = { name ->
    switch (name) {
        case "vjl":
            return new VJLClonotypeGroupingFactory()
        case "vj":
            return new VJClonotypeGroupingFactory()
        case "v":
            return new VClonotypeGroupingFactory()
        default:
            return DummyClonotypeGroupingFactory.INSTANCE
    }
}

def groupingFactory = getGroupingFactory(grouping),
    groupingFactory2 = getGroupingFactory(grouping2)

def scriptName = getClass().canonicalName.split("\\.")[-1]

// Prepare samples

println "[${new Date()} $scriptName] Reading samples"

def sampleCollection = metadataFileName ?
        new SampleCollection((String) metadataFileName) : // todo: if no control, perhaps store samples?
        new SampleCollection(opt.arguments()[0..-2])

println "[${new Date()} $scriptName] ${sampleCollection.size()} samples prepared"

// Compute control degree statistics

def bgDegreeStatCalc = new DegreeStatisticsCalculator(searchScope[0],
        searchScope[1], searchScope[2], groupingFactory, groupingFactory2)

if (backgroundSample) {
    // Load control sample
    println "[${new Date()} $scriptName] Loading control sample"
    bgDegreeStatCalc.inititalize(SampleFileConnection.load((String) backgroundSample))
} else {
    // Create control by aggregating all samples on the fly
    println "[${new Date()} $scriptName] No control sample provided. Creating pooled sample to be used as control"
    def sampleAggr = new SampleAggregator(sampleCollection,
            new StoringClonotypeAggregatorFactory(),
            OverlapType.Strict)
    bgDegreeStatCalc.inititalize(new PooledSample(sampleAggr))
}

// Compute degree statistics and write output

def sw = new SampleWriter(compress)

sampleCollection.eachWithIndex { sample, ind ->
    def sampleId = sample.sampleMetadata.sampleId
    println "[${new Date()} $scriptName] Computing degree statistics for $sampleId.."

    // Initialize sample degree statistics
    def sampleDegreeStatCalc = new DegreeStatisticsCalculator(searchScope[0],
            searchScope[1], searchScope[2], groupingFactory, groupingFactory2)
    sampleDegreeStatCalc.inititalize(sample)

    // Create annotator and annotate sample with both background and observed degree statistics
    def annotator = new DegreeStatisticsAnnotator(sampleDegreeStatCalc, bgDegreeStatCalc)
    annotator.annotate(sample)

    // print output
    sw.writeConventional(sample, outputFilePrefix)
}

sampleCollection.metadataTable.storeWithOutput(outputFilePrefix, compress,
        "degstat")

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
