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

class ImgtHighVQuestParser extends ClonotypeStreamParser {
    /**
     * {@inheritDoc}
     */
    protected ImgtHighVQuestParser(Iterator<String> innerIter, Sample sample) {
        super(innerIter, Software.ImgtHighVQuest, sample)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected Clonotype innerParse(String clonotypeString) {
        def splitString = clonotypeString.split("\t")

        if (splitString.length < 107)
            return null

        def count = 1
        def freq = 0

        def cdr3start = splitString[60].isInteger() ?
                splitString[60].toInteger() :
                -1 // this is called "junction start" here. Junction = CDR3 + conserved C, F/W

        def cdr3nt = splitString[15].toUpperCase()

        if (!(cdr3nt =~ /^[ATGCatgc]+$/))
            return null // no N's allowed

        def cdr3aa = toUnifiedCdr3Aa(translate(cdr3nt))

        String v, d, j
        (v, j, d) = extractVDJ(splitString[3..5]).collect {
            def splitRecord = it.split(" ")
            splitRecord.length > 1 ? splitRecord[1] : splitRecord[0]
        }

        def segmPoints = [
                splitString[63],
                splitString[76],
                splitString[77],
                splitString[106]
        ].collect {
            it.isInteger() ? (it.toInteger() - cdr3start) : -1 // subtract cdr3start
        }.collect {
            (it >= 0 && it < cdr3nt.length()) ? it : -1 // sometimes segment bounds appear out of junction region
        } as int[]

        boolean inFrame = inFrame(cdr3aa),
                noStop = noStop(cdr3aa), isComplete = cdr3aa.length() > 0


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
