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
 * A clonotype parser implementation that handles simple tab-delimited input, see
 * {@url https://github.com/mikessh/vdjtools/wiki/Input#simple}
 */
public class BaseParser extends ClonotypeStreamParser {
    /**
     * {@inheritDoc}
     */
    protected BaseParser(Iterator<String> innerIter, Sample sample) {
        super(innerIter, Software.VDJtools, sample)
    }

    /**
     * {@inheritDoc}
     */
    protected BaseParser(Iterator<String> innerIter, Software software, Sample sample) {
        super(innerIter, software, sample)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected Clonotype innerParse(String clonotypeString) {
        def splitString = clonotypeString.split(software.delimiter)

        def count = splitString[0].toInteger()
        def freq = splitString[1].toDouble()

        def cdr3nt = splitString[2]
        def cdr3aa = splitString[3].length() == 0 || splitString[3] == PLACEHOLDER ?
                translate(splitString[3]) : splitString[3]
        cdr3aa = toUnifiedCdr3Aa(cdr3aa)

        String v, d, j
        (v, d, j) = extractVDJ(splitString[4..6])

        boolean inFrame = inFrame(cdr3aa),
                noStop = noStop(cdr3aa),
                isComplete = true


        def segmPoints = (7..10).collect {
            (splitString.size() <= it || !splitString[it].isInteger()) ? -1 : splitString[it].toInteger()
        } as int[]

        new Clonotype(sample, count, freq,
                segmPoints, v, d, j,
                cdr3nt, cdr3aa,
                inFrame, noStop, isComplete,
                extractAnnotation(splitString))
    }

    private static String extractAnnotation(String[] splitLine) {
        splitLine.size() > 11 ? splitLine[11..-1].join("\t") : null
    }

    @Override
    String getAnnotationHeader() {
        extractAnnotation(header[0].split("\t"))
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
