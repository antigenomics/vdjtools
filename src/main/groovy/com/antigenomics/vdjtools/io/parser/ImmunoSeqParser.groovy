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

package com.antigenomics.vdjtools.io.parser

import com.antigenomics.vdjtools.misc.CommonUtil
import com.antigenomics.vdjtools.misc.Software
import com.antigenomics.vdjtools.sample.Clonotype
import com.antigenomics.vdjtools.sample.Sample

import static com.antigenomics.vdjtools.misc.CommonUtil.*

/**
 * A clonotype parser implementation that handles Adaptive Biotechnologies (tm) immunoSEQ (tm) assay
 * output format, see
 * {@url http://www.adaptivebiotech.com/content/immunoseq-0}
 *
 * This parser is intended for samples obtained using "Export sample" option from ImmunoSEQ analyzer,
 * not "Export sample V2"
 */
class ImmunoSeqParser extends ClonotypeStreamParser {
    protected boolean initialized = false
    protected int countColumn, countColumn2, freqColumn, cdr3StartColumn,
                  cdr3ntColumn, cdr3aaColumn, cdr3LenColumn,
                  jStartColumn, inFrameColumn,
            vColumn0, dColumn0, jColumn0,
            vColumn1, dColumn1, jColumn1,
                  vColumn2, dColumn2, jColumn2,
                  vEndColumn, dStartColumn, dEndColumn

    /**
     * {@inheritDoc}
     */
    protected ImmunoSeqParser(Iterator<String> innerIter, Sample sample) {
        super(innerIter, Software.ImmunoSeq, sample)
    }

    /**
     * {@inheritDoc}
     */
    protected ImmunoSeqParser(Iterator<String> innerIter, Software software, Sample sample) {
        super(innerIter, software, sample)
    }

    /**
     * Performs parser initialization based on the header string.
     *
     * @throws RuntimeException if header line doesn't contain required columns
     */
    protected synchronized void ensureInitialized() {
        if (initialized)
            return

        // Parsing header line to determine positions of certain columns with clones properties

        String headerLine = this.header[0];
        String[] splitHeaderLine = headerLine.split(software.delimiter)

        countColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("reads") }
        countColumn2 = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("templates") }
        freqColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("frequency") }
        cdr3StartColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("v_index") }
        cdr3LenColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("cdr3_length") }
        cdr3ntColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("rearrangement") }
        inFrameColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("frame_type") }
        cdr3aaColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("amino_acid") }
        vColumn0 = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("v_family_ties") }
        dColumn0 = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("d_family_ties") }
        jColumn0 = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("j_family_ties") }
        vColumn1 = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("v_family") }
        dColumn1 = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("d_family") }
        jColumn1 = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("j_family") }
        vColumn2 = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("v_gene") }
        dColumn2 = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("d_gene") }
        jColumn2 = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("j_gene") }
        vEndColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("n1_index") }
        dStartColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("d_index") }
        dEndColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("n2_index") }
        jStartColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("j_index") }

        if ([countColumn, countColumn2, freqColumn,
             cdr3StartColumn, cdr3LenColumn,
             cdr3ntColumn, cdr3aaColumn,
             vColumn0, dColumn0, jColumn0,
             vColumn1, dColumn1, jColumn1,
             vColumn2, dColumn2, jColumn2,
             vEndColumn, dStartColumn, dEndColumn, jStartColumn,
             inFrameColumn].any { it < 0 })
            throw new RuntimeException("Some mandatory columns are absent in the input file.")

        // Initialized
        initialized = true
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected Clonotype innerParse(String clonotypeString) {
        ensureInitialized()

        def splitString = clonotypeString.split("\t")

        // As-is data
        def count = splitString[countColumn].isInteger() ? splitString[countColumn].toInteger() : null

        if (count == null) {
            // try to rescue count from "templates" column
            count = splitString[countColumn2].toInteger()
        }

        def freq = splitString[freqColumn].toDouble()

        int cdr3start = splitString[cdr3StartColumn].toInteger(),
            cdr3Len = splitString[cdr3LenColumn].toInteger()

        def inFrame = splitString[inFrameColumn].equalsIgnoreCase("in")

        def cdr3End = cdr3start + cdr3Len,
            readLen = splitString[cdr3ntColumn].length(),
            padding = 0

        if (cdr3End > readLen) {
            // handling short reads
            padding = cdr3End - readLen
            cdr3End = readLen
        }

        def cdr3nt = splitString[cdr3ntColumn][cdr3start..<cdr3End] + ("N" * padding)
        def cdr3aa = toUnifiedCdr3Aa(inFrame ? splitString[cdr3aaColumn] : translate(cdr3nt))

        String v, d, j
        (v, d, j) = extractVDJImmunoSeq(
                splitString[[vColumn0, dColumn0, jColumn0]],
                splitString[[vColumn1, dColumn1, jColumn1]],
                splitString[[vColumn2, dColumn2, jColumn2]])

        // Fixing mess with CDR3s that are failed to be extracted

        def jStart = splitString[jStartColumn].toInteger()

        boolean isComplete = true
        if (cdr3start >= 0 &&
                (cdr3aa.length() == 0 || cdr3nt.length() != 3 * cdr3aa.length())) {
            cdr3nt = splitString[cdr3ntColumn]
            if (cdr3aa.length() > 0) {
                // see https://github.com/mikessh/vdjtools/issues/30 for the reason for workaround
                int to = cdr3start + cdr3aa.length() * 3
                isComplete = to <= cdr3nt.length()
                cdr3nt = isComplete ? cdr3nt.substring(cdr3start, to) : cdr3nt.substring(cdr3start) // in-frame
                cdr3aa = toUnifiedCdr3Aa(inFrame ? splitString[cdr3aaColumn] : translate(cdr3nt))
            } else {
                // it seems to be hard to get conventional out-of-frame translation here
                // but we'll try to reconstruct it
                if (jStart >= 0) {
                    def jRef = getJReferencePoint(cdr3nt.substring(jStart))
                    if (jRef >= 0) {
                        cdr3nt = cdr3nt.substring(cdr3start, jStart + jRef + 4)
                        cdr3aa = translate(cdr3nt).replaceAll(/([atgc#\?])+/, "~") // to unified look
                    }
                }
                isComplete = cdr3aa.length() > 0
            }
        }

        inFrame = inFrame && cdr3aa.length() > 0 && CommonUtil.inFrame(cdr3aa)
        boolean noStop = noStop(cdr3aa)

        // Correctly record segment points
        cdr3start = cdr3start < 0 ? 0 : cdr3start

        def segmPoints = [
                splitString[vEndColumn].toInteger() - 1 - cdr3start,
                splitString[dStartColumn].toInteger() - cdr3start,
                splitString[dEndColumn].toInteger() - 1 - cdr3start,
                splitString[jStartColumn].toInteger() - cdr3start].collect { it < 0 ? -1 : it } as int[]

        new Clonotype(sample, count, freq,
                segmPoints, v, d, j,
                cdr3nt, cdr3aa,
                inFrame, noStop, isComplete)
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
