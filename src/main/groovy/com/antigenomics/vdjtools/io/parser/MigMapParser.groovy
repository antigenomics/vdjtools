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

import static com.antigenomics.vdjtools.misc.CommonUtil.toUnifiedCdr3Aa

/**
 * A clonotype parser implementation that handles output from IgBlastWrapper-derived software.
 * {@url https://github.com/mikessh/igblastwrp}
 */
public class MigMapParser extends ClonotypeStreamParser {
    private boolean initialized = false
    private int countColumn, freqColumn, cdr3ntColumn, cdr3aaColumn,
                vColumn, dColumn, jColumn,
                vEndColumn, dStartColumn, dEndColumn, jStartColumn,
                inFrameColumn, noStopColumn, isCompleteColumn

    /**
     * {@inheritDoc}
     */
    protected MigMapParser(Iterator<String> innerIter, Sample sample) {
        super(innerIter, Software.MigMap, sample)
    }

    /**
     * Performs parser initialization based on the header string.
     *
     * @throws RuntimeException if header line doesn't contain required columns
     */
    private synchronized void ensureInitialized() {
        if (initialized)
            return

        // Parsing header line to determine positions of certain columns with clones properties

        String headerLine = this.header[0];
        String[] splitHeaderLine = headerLine.split(software.delimiter)

        countColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("count") }
        freqColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("freq") }
        cdr3ntColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("cdr3nt") }
        cdr3aaColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("cdr3aa") }
        vColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("v") }
        dColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("d") }
        jColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("j") }
        vEndColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("v.end.in.cdr3") }
        dStartColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("d.start.in.cdr3") }
        dEndColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("d.end.in.cdr3") }
        jStartColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("j.start.in.cdr3") }
        inFrameColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("in.frame") }
        noStopColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("no.stop") }
        isCompleteColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("complete") }

        if ([countColumn, freqColumn, cdr3ntColumn, cdr3aaColumn,
             vColumn, dColumn, jColumn,
             vEndColumn, dStartColumn, dEndColumn, jStartColumn,
             inFrameColumn, noStopColumn, isCompleteColumn].any { it < 0 })
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

        /*
            Old layout:
             0     1      2  3  4  5       6           7    8     9    10    11   12    13
            "freq\tcount\tv\td\tj\tcdr3nt\tcdr3aa\t" + FR1, CDR1, FR2, CDR2, FR3, CDR3, FR4 +
               14               15                  16    17      18    19 
            "\tcdr.insert.qual\tmutations.qual\t" + vEnd, dStart, dEnd, jStart + "\t" +
            20    21     22     23
            vDel, dDel5, dDel3, jDel +
            24      25       26       27
            "pol.v\tpol.d.5\tpol.d.3\tpol.j" +
               28        29        30       31        32
            "\thas.cdr3\tin.frame\tno.stop\tcomplete\tcanonical"
         */

        def splitString = clonotypeString.split("\t")

        def freq = splitString[freqColumn].toDouble()
        def count = splitString[countColumn].toInteger()
        def cdr3nt = splitString[cdr3ntColumn] == "." || splitString[cdr3ntColumn].contains("N") ?
                "" : splitString[cdr3ntColumn],
            cdr3aa = splitString[cdr3aaColumn] == "." || splitString[cdr3aaColumn].contains("X") ?
                    "" : toUnifiedCdr3Aa(splitString[cdr3aaColumn])

        String v, d, j
        (v, d, j) = CommonUtil.extractVDJ(splitString[[vColumn, dColumn, jColumn]])

        boolean inFrame, noStop, isComplete
        (inFrame, noStop, isComplete) = splitString[[inFrameColumn, noStopColumn, isCompleteColumn]].collect {
            it.toBoolean()
        }

        if (cdr3nt == "") {
            isComplete = false // handling Ns in CDR3
        }

        def segmPoints = splitString[[vEndColumn, dStartColumn, dEndColumn, jStartColumn]].collect {
            it.toInteger()
        } as int[]

        def clonotype = new Clonotype(sample,
                count, freq,
                segmPoints, v, d, j,
                cdr3nt, cdr3aa,
                inFrame, noStop, isComplete)

        return clonotype
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
