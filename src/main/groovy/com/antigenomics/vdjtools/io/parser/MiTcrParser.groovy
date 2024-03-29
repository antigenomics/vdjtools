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
 * A clonotype parser implementation that handles output from MiTCR software, see
 * {@url http://mitcr.milaboratory.com/}
 */
public class MiTcrParser extends ClonotypeStreamParser {
    /**
     * {@inheritDoc}
     */
    protected MiTcrParser(Iterator<String> innerIter, Sample sample) {
        super(innerIter, Software.MiTcr, sample)
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
        def cdr3aa = toUnifiedCdr3Aa(splitString[5]) // replace ~

        String v, d, j
        (v, d, j) = extractVDJ(splitString[[7, 11, 9]])

        def segmPoints = [splitString[12].toInteger(),
                          splitString[13].isInteger() ? splitString[13].toInteger() : -1,
                          splitString[14].isInteger() ? splitString[14].toInteger() : -1,
                          splitString[15].toInteger()] as int[]

        boolean inFrame = inFrame(cdr3aa),
                noStop = noStop(cdr3aa),
                isComplete = true

        new Clonotype(sample, count, freq,
                segmPoints, v, d, j, cdr3nt, cdr3aa,
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
