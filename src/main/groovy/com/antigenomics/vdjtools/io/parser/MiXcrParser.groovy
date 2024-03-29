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
 * A clonotype parser implementation that handles output from MiXCR software, see
 * {@url http://mixcr.milaboratory.com/}
 */
public class MiXcrParser extends ClonotypeStreamParser {
    private boolean initialized = false
    private int countColumn, freqColumn, cdr3ntColumn, cdr3aaColumn,
                vHitsColumn, dHitsColumn, jHitsColumn,
                vAlignmentsColumn, dAlignmentsColumn, jAlignmentsColumn,
                numberOfColumns

    /**
     * {@inheritDoc}
     */
    protected MiXcrParser(Iterator<String> innerIter, Sample sample) {
        super(innerIter, Software.MiXcr, sample)
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

        countColumn = splitHeaderLine.findIndexOf {
            it.equalsIgnoreCase("Clone count") || it.equalsIgnoreCase("cloneCount")
        }
        freqColumn = splitHeaderLine.findIndexOf {
            it.equalsIgnoreCase("Clone fraction") || it.equalsIgnoreCase("cloneFraction")
        }
        cdr3ntColumn = splitHeaderLine.findIndexOf {
            it.equalsIgnoreCase("N. Seq. CDR3") || it.equalsIgnoreCase("nSeqCDR3") || it.equalsIgnoreCase("nSeqImputedCDR3")
        }
        cdr3aaColumn = splitHeaderLine.findIndexOf {
            it.equalsIgnoreCase("AA. Seq. CDR3") || it.equalsIgnoreCase("aaSeqCDR3") || it.equalsIgnoreCase("aaSeqImputedCDR3")
        }
        vAlignmentsColumn = splitHeaderLine.findIndexOf {
            it =~ /(?i)V alignment/ || it =~ /(?i)VAlignment/ || it =~ /(?i)allVAlignments/
        }
        dAlignmentsColumn = splitHeaderLine.findIndexOf {
            it =~ /(?i)D alignment/ || it =~ /(?i)DAlignment/ || it =~ /(?i)allDAlignments/
        }
        jAlignmentsColumn = splitHeaderLine.findIndexOf {
            it =~ /(?i)J alignment/ || it =~ /(?i)JAlignment/ || it =~ /(?i)allJAlignments/
        }
        vHitsColumn = splitHeaderLine.findIndexOf {
            it =~ /(?i)V hits/ || it =~ /(?i)VHits/ || it =~ /(?i)allVHitsWithScore/
        }
        dHitsColumn = splitHeaderLine.findIndexOf {
            it =~ /(?i)D hits/ || it =~ /(?i)DHits/ || it =~ /(?i)allDHitsWithScore/
        }
        jHitsColumn = splitHeaderLine.findIndexOf {
            it =~ /(?i)J hits/ || it =~ /(?i)JHits/ || it =~ /(?i)allJHitsWithScore/
        }
        if (countColumn == -1 || freqColumn == -1 || cdr3ntColumn == -1 || cdr3aaColumn == -1 ||
                vAlignmentsColumn == -1 || dAlignmentsColumn == -1 || jAlignmentsColumn == -1)
            throw new RuntimeException("Some mandatory columns are absent in the input file.")

        //println([countColumn, freqColumn,cdr3ntColumn,cdr3aaColumn,vAlignmentsColumn,dAlignmentsColumn,jAlignmentsColumn])

        numberOfColumns = splitHeaderLine.size()

        // Initialized
        initialized = true
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected Clonotype innerParse(String clonotypeString) {
        ensureInitialized()

        def splitString = clonotypeString.split(software.delimiter, numberOfColumns)

        def count = (int) (splitString[countColumn].toDouble())
        def freq = splitString[freqColumn].toDouble()

        def cdr3nt = splitString[cdr3ntColumn]

        def cdr3aa = splitString[cdr3aaColumn] // no need to unify, MiXCR is based on milib

        String v, d, j
        (v, d, j) = extractVDJ(splitString[[vHitsColumn, dHitsColumn, jHitsColumn]])

        List<Alignment> vAlignemtns = parseAlignments(splitString[vAlignmentsColumn])
        List<Alignment> dAlignemtns = parseAlignments(splitString[dAlignmentsColumn])
        List<Alignment> jAlignemtns = parseAlignments(splitString[jAlignmentsColumn])

        def segmPoints = [vAlignemtns.size() > 0 && vAlignemtns[0] != null ?
                                  vAlignemtns[0].seq2End - 1 : 0,
                          dAlignemtns.size() > 0 && dAlignemtns[0] != null ?
                                  dAlignemtns[0].seq2Begin : -1,
                          dAlignemtns.size() > 0 && dAlignemtns[0] != null ?
                                  dAlignemtns[0].seq2End - 1 : -1,
                          jAlignemtns.size() > 0 && jAlignemtns[0] != null ?
                                  jAlignemtns[0].seq2Begin : cdr3nt.size() - 1] as int[]

        boolean inFrame = inFrame(cdr3aa),
                noStop = noStop(cdr3aa),
                isComplete = true

        new Clonotype(sample, count, freq,
                segmPoints, v, d, j, cdr3nt, cdr3aa,
                inFrame, noStop, isComplete)
    }

    private static List<Alignment> parseAlignments(String alignmentsLine) {
        if (alignmentsLine.isEmpty())
            return Collections.EMPTY_LIST

        String[] splitByAlignments = alignmentsLine.split(";", -1)
        List<Alignment> ret = new ArrayList<>()
        for (String alignmentString : splitByAlignments) {
            if (alignmentString.trim().empty) {
                ret.add(null)
            } else {
                String[] splitByFields = alignmentString.split("\\|")
                ret.add(new Alignment(splitByFields[0].toInteger(), splitByFields[1].toInteger(),
                        splitByFields[3].toInteger(), splitByFields[4].toInteger()));
            }
        }
        ret
    }

    /**
     * Now all positions are in coordinates of clonal sequence (not CDR3 coordinates). In case of non-CDR3 clonal
     * sequence there is currently no way to convert positions to CDR3 positions, so non-CDR3 positions will be
     * forbidden before appropriate information will be made available in MiXCR output.
     *
     * Seq1 - reference sequence
     *
     * Seq2 - clonal sequence
     */
    private static class Alignment {
        int seq1Begin, seq1End, seq2Begin, seq2End

        Alignment(int seq1Begin, int seq1End, int seq2Begin, int seq2End) {
            this.seq1Begin = seq1Begin
            this.seq1End = seq1End
            this.seq2Begin = seq2Begin
            this.seq2End = seq2End
        }
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
