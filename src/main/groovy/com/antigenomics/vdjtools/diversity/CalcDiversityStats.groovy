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

package com.antigenomics.vdjtools.diversity

import com.antigenomics.vdjtools.overlap.OverlapType
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.sample.metadata.MetadataTable

import static com.antigenomics.vdjtools.misc.ExecUtil.formOutputPath

def I_TYPE_DEFAULT = OverlapType.Strict, RESAMPLES_DEFAULT = "3"
def cli = new CliBuilder(usage: "CalcDiversityStats [options] " +
        "[sample1 sample2 sample3 ... if -m is not specified] output_prefix")
cli.h("display help message")
cli.m(longOpt: "metadata", argName: "filename", args: 1,
        "Metadata file. First and second columns should contain file name and sample id. " +
                "Header is mandatory and will be used to assign column names for metadata.")
cli.x(longOpt: "downsample-to", argName: "integer", args: 1,
        "Number of reads to take for down-sampled sample diversity estimate." +
                "[default=number of reads in the smallest sample]")
cli.X(longOpt: "extrapolate-to", argName: "integer", args: 1,
        "Number of reads to take for extrapolated sample diversity estimate." +
                "[default=number of reads in the largest sample]")
cli.i(longOpt: "intersect-type", argName: "string", args: 1,
        "Intersection rule to apply. Allowed values: $OverlapType.allowedNames. " +
                "Will use '$I_TYPE_DEFAULT' by default.")
cli._(longOpt: "resample-trials", argName: "integer", args: 1,
        "Number of resamples for corresponding estimator. [default = $RESAMPLES_DEFAULT]")

def opt = cli.parse(args)

if (opt == null)
    System.exit(2)

if (opt.h || opt.arguments().size() == 0) {
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

// Other arguments

def intersectionType = opt.i ? OverlapType.getByShortName((String) opt.i) : I_TYPE_DEFAULT,
    resampleCount = (opt."resample-trials" ?: RESAMPLES_DEFAULT).toInteger(),
    outputPrefix = opt.arguments()[-1]

def scriptName = getClass().canonicalName.split("\\.")[-1]

//
// Batch load all samples (lazy)
//

println "[${new Date()} $scriptName] Reading samples"

def sampleCollection = metadataFileName ?
        new SampleCollection((String) metadataFileName) :
        new SampleCollection(opt.arguments()[0..-2])

println "[${new Date()} $scriptName] ${sampleCollection.size()} samples to analyze"

//
// Set up downsample/extrapolate read counts
//

def sampleStats = (!opt.x || !opt.X) ? sampleCollection.sampleStatistics : null
def minReads = (opt.x ?: "$sampleStats.minCount").toInteger(),
    maxReads = (opt.X ?: "$sampleStats.maxCount").toInteger()

//
// Compute and output diversity measures
//

def headerBase = "$MetadataTable.SAMPLE_ID_COLUMN\t" +
        sampleCollection.metadataTable.columnHeader + "\treads\tdiversity"

def exactOutputPath = formOutputPath(outputPrefix, "diversity",
        intersectionType.shortName, EstimationMethod.Exact.name),
    resampledOutputPath = formOutputPath(outputPrefix, "diversity",
            intersectionType.shortName, EstimationMethod.Resampled.name)

new File(exactOutputPath).withPrintWriter { pwExact ->
    new File(resampledOutputPath).withPrintWriter { pwResampling ->
        pwExact.println(headerBase + "\textrapolate_reads\t" + ExactEstimator.HEADER)
        pwResampling.println(headerBase + "\tresample_reads\t" + ResamplingEstimator.HEADER)

        sampleCollection.each { Sample sample ->
            println "[${new Date()} $scriptName] Analyzing $sample.sampleMetadata.sampleId"

            def rowBase = [sample.sampleMetadata.sampleId, sample.sampleMetadata,
                           sample.count, sample.diversity].join("\t")

            def exactEstimator = new ExactEstimator(sample, intersectionType, maxReads),
                resamplingEstimator = new ResamplingEstimator(sample, intersectionType, minReads, resampleCount)

            pwExact.println(rowBase + "\t" + maxReads + "\t" + exactEstimator)
            pwResampling.println(rowBase + "\t" + minReads + "\t" + resamplingEstimator)
        }
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
