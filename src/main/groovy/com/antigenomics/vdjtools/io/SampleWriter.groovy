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

package com.antigenomics.vdjtools.io

import com.antigenomics.vdjtools.misc.Software
import com.antigenomics.vdjtools.join.JointSample
import com.antigenomics.vdjtools.pool.PooledSample
import com.antigenomics.vdjtools.sample.Clonotype
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.misc.ExecUtil

import java.util.zip.GZIPOutputStream

/**
 * A class implementing output of Sample and JointSample to plain-text file
 */
public class SampleWriter {
    private final Software software
    private final String header
    private final boolean compress
    private final List<String> printFields

    public String getHeader() {
        header
    }

    public String getClonotypeString(Clonotype clonotype) {
        printFields.collect {
            clonotype."$it"
        }.join("\t")
    }

    /**
     * Creates a sample writer capable to output samples in plain-text files according to specified software format.
     * @param compress specifies whether to compress resulting output file
     * @param renormalize tells whether to perform re-normalization (compute frequency by dividing read count by total
     *                    number of reads in sample) or preserve original frequencies as in input
     */
    public SampleWriter(boolean compress = false, boolean renormalize = false) {
        this(Software.VDJtools, compress, renormalize)
    }

    /**
     * Creates a sample writer capable to output samples in plain-text files according to specified software format.
     * Will compress the resulting output file if {@code compress = true} and append ".gz" to the file name.
     * @param software table layout that will be used during output
     * @param compress specifies whether to compress resulting output file
     * @param renormalize tells whether to perform re-normalization (compute frequency by dividing read count by total
     *                    number of reads in sample) or preserve original frequencies as in input
     *
     * @deprecated writing back in various output formats is not supported by CLI
     */
    @Deprecated
    public SampleWriter(Software software, boolean compress, boolean renormalize) {
        this.software = software
        this.header = (software.headerLineCount > 1 ?
                "${software.name()}-header-blank\n" * (software.headerLineCount - 1) : "") +
                (software.comment ?: "") +
                software.printFields.join("\t")
        this.compress = compress
        this.printFields = renormalize ? software.printFields :
                software.printFields.collect { it.replace("freq", "freqAsInInput") }
    }

    /**
     * Gets buffered writer 
     * @param fileName
     * @return
     */
    public BufferedWriter getWriter(String fileName) {
        if (compress)
            fileName += ".gz"

        def fos = new FileOutputStream(fileName)

        new BufferedWriter(new OutputStreamWriter(compress ?
                new GZIPOutputStream(fos) : fos))
    }

    /**
     * Writes a sample to a given directory/path using conventional sample naming suffix. 
     * Will assume that output to a directory is required if path ends with "/".
     * @param sample sample to write
     * @param outputPrefix output path prefix or output directory (if end with "/")
     */
    public void writeConventional(Sample sample, String outputPrefix) {
        write(sample, ExecUtil.formOutputPath(outputPrefix, sample))
    }

    /**
     * Writes a sample to a given directory/path using conventional sample naming suffix. 
     * Will assume that output to a directory is required if path ends with "/".
     * @param sample sample to write
     * @param outputPrefix output path prefix or output directory (if end with "/")
     * @param top number of top clonotype to write, -1 to write all clonotypes
     * @param collapse specifies whether to store the information on clonotypes that not got it to {@code top}
     *        ones as a separate single entry put at the end of the file
     */
    public void writeConventional(Sample sample, String outputPrefix, int top, boolean collapse) {
        write(sample, ExecUtil.formOutputPath(outputPrefix, sample), top, collapse)
    }

    /**
     * Writes a sample as a plain-text table to the specified path.
     * @param sample sample to write
     * @param fileName output path
     */
    public void write(Sample sample, String fileName) {
        write(sample, fileName, -1, false)
    }

    private static String appendAnnotation(String annotation) {
        annotation ? ("\t" + annotation) : ""
    }

    private static String appendAnnotationDummy(String annotation) {
        annotation ? ("\t" + annotation.split("\t").collect { "" }.join("\t")) : ""
    }

    String getFullHeader(Sample sample) {
        header + appendAnnotation(sample.annotationHeader)
    }

    String getFullClonotypeString(Clonotype clonotype) {
        getClonotypeString(clonotype) + appendAnnotation(clonotype.annotation)
    }

    /**
     * Writes a sample as a plain-text table to the specified path.
     * @param sample sample to write
     * @param fileName output path
     * @param top number of top clonotype to write, -1 to write all clonotypes
     * @param collapse specifies whether to store the information on clonotypes that not got it to {@code top}
     *        ones as a separate single entry put at the end of the file
     */
    public void write(Sample sample, String fileName, int top, boolean collapse) {
        def printWriter = getWriter(fileName)

        top = top > sample.diversity || top < 0 ? sample.diversity : top
        printWriter.println(getFullHeader(sample))

        long count = 0
        double freq = 0.0

        for (int i = 0; i < top; i++) {
            def clonotype = sample[i]

            if (collapse) {
                count += clonotype.count
                freq += clonotype.freq
            }

            printWriter.println(getFullClonotypeString(clonotype))
        }

        if (collapse && top < sample.diversity) {
            // Collapsed
            printWriter.println(software.printFields.collect {
                if (it == "count")
                    sample.count - count
                else if (it == "freq")
                    sample.freqAsInInput - freq
                else
                    "NotShown"
            }.join("\t") + appendAnnotationDummy(sample.annotationHeader))
        }

        printWriter.close()
    }

    /**
     * Writes a joint sample as a plain-text table to the specified path.
     * @param sample sample to write
     * @param fileName output path
     * @param top number of top clonotype to write, -1 to write all clonotypes
     * @param collapse specifies whether to store the information on clonotypes that not got it to {@code top}
     *        ones as a separate single entry put at the end of the file
     */
    public void write(JointSample jointSample, String fileName, int top, boolean collapse) {
        def printWriter = getWriter(fileName)
        top = top > jointSample.diversity || top < 0 ? jointSample.diversity : top

        def sampleIndices = (0..<jointSample.numberOfSamples)

        printWriter.println(header + "\tpeak\toccurences\tsampling.p\t" +
                sampleIndices.collect { jointSample.getSample(it).sampleMetadata.sampleId }.join("\t"))

        double collapsedMeanFreq = 0.0
        double[] freqArr = new double[jointSample.numberOfSamples]

        for (int i = 0; i < top; i++) {
            def jointClonotype = jointSample[i],
                clonotype = jointClonotype.clonotype

            if (collapse) {
                collapsedMeanFreq += jointClonotype.baseFreq
                sampleIndices.each { int j ->
                    freqArr[j] += jointClonotype.getFreq(j)
                }
            }

            printWriter.println(
                    [software.printFields.collect {
                        if (it == "count")
                            jointClonotype.count
                        else if (it == "freq")
                            jointClonotype.freq
                        else
                            clonotype."$it"
                    },
                     jointClonotype.peak,
                     jointClonotype.occurrences,
                     jointClonotype.samplingPValue,
                     sampleIndices.collect { int j ->
                         jointClonotype.getFreq(j)
                     }].flatten().join("\t"))
        }

        if (collapse) {
            if (top < jointSample.diversity) {
                // Collapsed
                def collapsedFreq = jointSample.totalMeanFreq - collapsedMeanFreq,
                    collapsedFreqArr = sampleIndices.collect { int j ->
                        jointSample.getIntersectionFreq(j) - freqArr[j]
                    }
                printWriter.println(
                        [software.printFields.collect {
                            if (it == "count")
                                jointSample.calcCount(collapsedFreq)
                            else if (it == "freq")
                                jointSample.calcFreq(collapsedFreq)
                            else
                                "NotShown"
                        },
                         collapsedFreqArr.findIndexOf { it == collapsedFreqArr.max() },
                         collapsedFreqArr.findAll { it > 0 }.sum(),
                         0,
                         collapsedFreqArr
                        ].flatten().join("\t"))
            }

            // Not in the overlap
            def nonOverlappingFreqArr = sampleIndices.collect { int j ->
                jointSample.getSample(j).freqAsInInput - jointSample.getIntersectionFreq(j)
            }
            printWriter.println(
                    [software.printFields.collect {
                        if (it == "count")
                            0
                        else if (it == "freq")
                            0.0
                        else
                            "NonOverlapping"
                    },
                     nonOverlappingFreqArr.findIndexOf { it == nonOverlappingFreqArr.max() },
                     nonOverlappingFreqArr.findAll { it > 0 }.sum(),
                     0,
                     nonOverlappingFreqArr
                    ].flatten().join("\t"))
        }

        printWriter.close()
    }

    /**
     * Writes a pooled sample as a plain-text table to the specified path.
     * @param sample sample to write.
     * @param fileName output path.
     */
    public void write(PooledSample pooledSample, String fileName) {
        def printWriter = getWriter(fileName)

        printWriter.println(header + "\tincidence\tconvergence\toccurrences")

        pooledSample.each { pooledClonotype ->
            printWriter.println(
                    [software.printFields.collect {
                        if (it == "count")
                            pooledClonotype.count
                        else if (it == "freq")
                            pooledClonotype.count / (double) pooledSample.count
                        else
                            pooledClonotype.clonotype."$it"
                    },
                     pooledClonotype.incidenceCount,
                     pooledClonotype.diversity,
                     pooledClonotype.occurrences
                    ].flatten().join("\t"))
        }

        printWriter.close()
    }

    /**
     * Writes a joint sample as a plain-text table to the specified path.
     * @param sample sample to write
     * @param fileName output path
     */
    public void write(JointSample jointSample, String fileName) {
        write(jointSample, fileName, -1, false)
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
