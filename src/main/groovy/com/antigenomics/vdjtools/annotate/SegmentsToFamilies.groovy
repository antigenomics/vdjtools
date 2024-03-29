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

import com.antigenomics.vdjtools.io.SampleWriter
import com.antigenomics.vdjtools.misc.CommonUtil
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.sample.SegmentConverter

def cli = new CliBuilder(usage: "SegmentsToFamilies [options] " +
        "[sample1 sample2 ... if not -m] output_prefix")
cli.h("display help message")
cli.m(longOpt: "metadata", argName: "filename", args: 1,
        "Metadata file. First and second columns should contain file name and sample id. " +
                "Header is mandatory and will be used to assign column names for metadata.")
cli.s(longOpt: "species", argName: "name", args: 1,
        "Dataset species, 'human' or 'mouse'.", required: true)
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

def outputFilePrefix = opt.arguments()[-1],
    species = (String) opt.s,
    compress = (boolean) opt.c

if (!["human", "mouse"].any { it.equalsIgnoreCase(species) }) {
    println "Should specify either human or mouse as species."
    System.exit(2)
}

def scriptName = getClass().canonicalName.split("\\.")[-1]

//
// Batch load samples
//

println "[${new Date()} $scriptName] Reading samples"

def sampleCollection = metadataFileName ?
        new SampleCollection((String) metadataFileName) :
        new SampleCollection(opt.arguments()[0..-2])

println "[${new Date()} $scriptName] ${sampleCollection.size()} samples prepared"

// Load segment conversions
def vSegmentMap = new HashMap<String, String>(),
    jSegmentMap = new HashMap<String, String>()
CommonUtil.resourceStreamReader("vj_families.txt").splitEachLine("\t") {
    if (it[0].equalsIgnoreCase(species)) {
        if (it[2].equalsIgnoreCase("v")) {
            vSegmentMap.put(it[3], it[4])
        } else {
            jSegmentMap.put(it[3], it[4])
        }
    }
}

def converter = new SegmentConverter(vSegmentMap, jSegmentMap)

//
// Iterate over samples and change V segments
//
def sw = new SampleWriter(compress)

sampleCollection.eachWithIndex { sample, ind ->
    def sampleId = sample.sampleMetadata.sampleId
    println "[${new Date()} $scriptName] Changing segments for $sampleId.."

    // print output
    sw.writeConventional(new Sample(sample, converter), outputFilePrefix)
}

sampleCollection.metadataTable.storeWithOutput(outputFilePrefix, compress,
        "segm2fam")

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
