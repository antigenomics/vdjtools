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

import com.antigenomics.vdjtools.misc.Software
import com.antigenomics.vdjtools.sample.Clonotype
import com.antigenomics.vdjtools.sample.Sample

import static com.antigenomics.vdjtools.misc.CommonUtil.*

/**
 * A clonotype parser implementation that handles Adaptive Biotechnologies (tm) immunoSEQ (tm) assay
 * output format, see
 * {@url http://www.adaptivebiotech.com/content/immunoseq-0}
 *
 * This parser is intended for samples obtained using "Export sample V2" option from ImmunoSEQ analyzer,
 * not "Export sample"
 */
class ImmunoSeqV2Parser extends ImmunoSeqParser {
    /**
     * {@inheritDoc}
     */
    protected ImmunoSeqV2Parser(Iterator<String> innerIter, Sample sample) {
        super(innerIter, Software.ImmunoSeqV2, sample)
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

        countColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("count (templates/reads)") }
        countColumn2 = -1
        freqColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("frequencyCount (%)") }
        cdr3StartColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("vIndex") }
        cdr3LenColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("cdr3Length") }
        cdr3ntColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("nucleotide") }
        inFrameColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("sequenceStatus") }
        cdr3aaColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("aminoAcid") }
        vColumn0 = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("vFamilyTies") }
        dColumn0 = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("dFamilyTies") }
        jColumn0 = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("jFamilyTies") }
        vColumn1 = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("vFamilyName") }
        dColumn1 = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("dFamilyName") }
        jColumn1 = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("jFamilyName") }
        vColumn2 = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("vGeneName") }
        dColumn2 = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("dGeneName") }
        jColumn2 = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("jGeneName") }
        vEndColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("n1Index") }
        dStartColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("dIndex") }
        dEndColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("n2Index") }
        jStartColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("jIndex") }

        if ([countColumn, freqColumn,
             cdr3StartColumn,  cdr3LenColumn,
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
