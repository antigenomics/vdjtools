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

package com.antigenomics.vdjtools.operate

import com.antigenomics.vdjtools.diversity.FrequencyTable
import com.antigenomics.vdjtools.misc.Software
import com.antigenomics.vdjtools.io.SampleWriter
import com.antigenomics.vdjtools.overlap.OverlapType
import com.antigenomics.vdjtools.pool.PooledSample
import com.antigenomics.vdjtools.pool.SampleAggregator
import com.antigenomics.vdjtools.pool.StoringClonotypeAggregatorFactory
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.misc.ExecUtil

import static com.antigenomics.vdjtools.misc.ExecUtil.formOutputPath

def I_TYPE_DEFAULT = "aa"
def cli = new CliBuilder(usage: "PoolSamples [options] " +
        "[sample1 sample2 sample3 ... if not -m] output_prefix")
cli.h("display help message")
cli.m(longOpt: "metadata", argName: "filename", args: 1,
        "Metadata file. First and second columns should contain file name and sample id. " +
                "Header is mandatory and will be used to assign column names for metadata.")
cli.i(longOpt: "intersect-type", argName: "string", args: 1,
        "Comma-separated list of overlap types to apply. " +
                "Allowed values: $OverlapType.allowedNames. " +
                "Will use '$I_TYPE_DEFAULT' by default.")
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

def outputPrefix = opt.arguments()[-1],
    compress = (boolean) opt.c

ExecUtil.ensureDir(outputPrefix)

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

def sampleCollection = metadataFileName ?
        new SampleCollection((String) metadataFileName, Software.VDJtools, false, true) :
        new SampleCollection(opt.arguments()[0..-2], Software.VDJtools, false, true)

println "[${new Date()} $scriptName] ${sampleCollection.size()} samples loaded"

//
// Pool samples
//

println "[${new Date()} $scriptName] Pooling with $intersectionType, this may take a while"

def cloneAggrFact = new StoringClonotypeAggregatorFactory()

def sampleAggr = new SampleAggregator(sampleCollection, cloneAggrFact, intersectionType)

//
// Sort pooled sample and write output
//

println "[${new Date()} $scriptName] Normalizing and sorting pooled clonotypes"
def pooledSample = new PooledSample(sampleAggr)

println "[${new Date()} $scriptName] Writing output"
def writer = new SampleWriter(compress, true)
writer.write(pooledSample, formOutputPath(outputPrefix, "pool", intersectionType.shortName, "table"))

def incidenceTable = new FrequencyTable(pooledSample)

new File(formOutputPath(outputPrefix, "pool", intersectionType.shortName, "summary")).withPrintWriter { pw ->
    pw.println("clonotype.incidence\tnumber.of.clonotypes")
    pw.println(incidenceTable.toString())
}

println "[${new Date()} $scriptName] Finished."


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
