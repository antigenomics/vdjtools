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
 * A clonotype parser implementation that handles output from MiGEC software, see
 * {@url https://github.com/mikessh/migec}
 */
public class MiGecParser extends ClonotypeStreamParser {
    /**
     * {@inheritDoc}
     */
    protected MiGecParser(Iterator<String> innerIter, Sample sample) {
        super(innerIter, Software.MiGec, sample)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected Clonotype innerParse(String clonotypeString) {
        /*
             0  Count
             1  Fraction
             2  CDR3 nucleotide sequence
             3  CDR3 amino acid sequence
             4  V segments
             5  J segments
             6  D segments
             7  Last V nucleotide position
             8  First D nucleotide position
             9  Last D nucleotide position
             10 First J nucleotide position
             11 Good events
             12 Total events
             13 Good reads
             14 Total reads
          */

        def splitString = clonotypeString.split("\t")

        def count = splitString[0].toInteger()
        def freq = splitString[1].toDouble()

        def cdr3nt = splitString[2]
        def cdr3aa = toUnifiedCdr3Aa(splitString[3])


        String v, j, d
        (v, j, d) = extractVDJ(splitString[4..6])

        boolean inFrame = inFrame(cdr3aa), noStop = noStop(cdr3aa), isComplete = true

        def segmPoints = [
                splitString[7].toInteger(),
                splitString[8].isInteger() ? splitString[8].toInteger() : -1,
                splitString[9].isInteger() ? splitString[9].toInteger() : -1,
                splitString[10].toInteger()] as int[]

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
