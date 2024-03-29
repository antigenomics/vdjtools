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

import com.antigenomics.vdjtools.misc.Segment
import com.antigenomics.vdjtools.misc.Software
import com.antigenomics.vdjtools.sample.Clonotype
import com.antigenomics.vdjtools.sample.Sample
import com.milaboratory.core.sequence.Sequence

/**
 * Base class for providing parsing of various RepSeq software output.
 * The stream parser is not thread-safe. 
 */
public abstract class ClonotypeStreamParser implements Iterable<Clonotype> {
    private static final int WARNINGS_TO_DISPLAY = 5
    protected final List<String> header = new ArrayList<>()
    protected final Software software
    protected final Iterator<String> innerIter
    protected final Sample sample
    private int skippedLineCount = 0, commentLineCount
    private boolean printedWarning = false
    private final boolean hasComment
    protected int totalLines = 0
    protected final Map<String, Integer> badClonotypeStats = ["NO_CDR3NT" : 0,
                                                              "NO_CDR3AA" : 0,
                                                              "NO_V"      : 0,
                                                              "NO_J"      : 0,
                                                              "ZERO_COUNT": 0,
                                                              "ZERO_FREQ" : 0,
                                                              "BAD_LINE"  : 0];

    /**
     * Creates a new instance of clonotype parser. 
     * It is a clonotype factory wrapped around the specified input stream
     * @param innerIter object that iterates over file lines, i.e. rows in the clonotype table
     * @param software software used to create the clonotype table. Specifies the parser
     * @param sample a blank sample to fill up with clonotypes
     */
    protected ClonotypeStreamParser(Iterator<String> innerIter, Software software, Sample sample) {
        this.software = software
        this.innerIter = innerIter
        this.sample = sample
        this.hasComment = software.comment && software.comment.length() > 0
        this.commentLineCount = software.headerLineCount
    }

    /**
     * Gets an instance of clonotype parser.
     * It is a clonotype factory wrapped around the specified input stream
     * @param inputStream input stream that will be read
     * @param software software used to create the clonotype table. Specifies the parser
     * @param sample a blank sample to fill up with clonotypes
     * @return clonotype parser object
     */
    public static ClonotypeStreamParser create(InputStream inputStream, Software software, Sample sample) {
        ClonotypeStreamParser parser
        def reader = new BufferedReader(new InputStreamReader(inputStream))
        def innerIter = reader.iterator()

        switch (software) {
            case Software.MiTcr:
                parser = new MiTcrParser(innerIter, sample)
                break
            case Software.MigMap:
                parser = new MigMapParser(innerIter, sample)
                break
            case Software.VDJtools:
            case Software.VDJtoolsRenorm:
                parser = new BaseParser(innerIter, sample)
                break
            case Software.MiGec:
                parser = new MiGecParser(innerIter, sample)
                break
            case Software.ImmunoSeq:
                parser = new ImmunoSeqParser(innerIter, sample)
                break
            case Software.ImmunoSeqV2:
                parser = new ImmunoSeqV2Parser(innerIter, sample)
                break
            case Software.ImgtHighVQuest:
                parser = new ImgtHighVQuestParser(innerIter, sample)
                break
            case Software.MiXcr:
                parser = new MiXcrParser(innerIter, sample)
                break
            case Software.ImSeq:
                parser = new ImSeqParser(innerIter, sample)
                break
            case Software.Vidjil:
                parser = new VidjilParser(innerIter, sample)
                break
            case Software.RTCR:
                parser = new RtcrParser(innerIter, sample)
                break
            default:
                throw new UnsupportedOperationException("Don't know how to parse $software data")
        }

        (0..<software.headerLineCount).each { parser.header.add(innerIter.next()) }

        return parser
    }

    /**
     * Parses a string into clonotype in a {@code Software}-dependent manner 
     * @param clonotypeString string to parse
     * @return string to parse
     */
    protected abstract Clonotype innerParse(String clonotypeString)

    /**
     * Parses a string into clonotype in a {@code Software}-dependent manner and performs some consistency checks.
     * Skips input strings that result in incomplete / bad clonotypes.
     * @param clonotypeString string to parse
     * @return a clonotype instance or {@code null} if input string was skipped
     */
    public Clonotype parse(String clonotypeString) {
        try {
            if (hasComment && clonotypeString.startsWith(software.comment)) {
                commentLineCount++
                return null
            }

            if (clonotypeString.trim().empty) {
                return null // ignore blank lines
            }

            totalLines++

            def clonotype = innerParse(clonotypeString)

            def badFieldMap = clonotype ? ["NO_CDR3NT" : missingEntry(clonotype.cdr3ntBinary),
                                           "NO_CDR3AA" : missingEntry(clonotype.cdr3aaBinary),
                                           "NO_V"      : missingEntry(clonotype.VBinary),
                                           "NO_J"      : missingEntry(clonotype.JBinary),
                                           "ZERO_COUNT": clonotype.count == 0,
                                           "ZERO_FREQ" : !software.perReadOutput && clonotype.freqAsInInput == 0] :
                    ["BAD_LINE": true]

            if (badFieldMap.any { it.value }) {
                badFieldMap.each {
                    if (it.value) {
                        badClonotypeStats[it.key]++
                    }
                }
                if (!printedWarning) {
                    printedWarning = true
                    println "[WARNING] Some of the essential fields are bad/missing " +
                            "for the following clonotype string (displaying first $WARNINGS_TO_DISPLAY warnings)"
                }
                if (skippedLineCount++ < WARNINGS_TO_DISPLAY) {
                    println badFieldMap.findAll { it.value }.collect { it.key }.join(",") + ":"
                    println "$clonotypeString"
                }
                return null
            }

            return clonotype
        } catch (Exception e) {
            throw new RuntimeException("Unable to parse clonotype string $clonotypeString " +
                    "for $software input type: ${e.message}", e)
        }
    }

    /**
     * Gets the header string for clonotype annotation lines (if available)
     * @return annotation line header or null if not available
     */
    public String getAnnotationHeader() {
        null
    }

    /**
     * Gets the number of skipped input lines
     * @return number of input lines that resulted in bad/incomplete clonotypes
     */
    public int getSkippedLineCount() {
        skippedLineCount
    }

    /**
     * Gets the number of comment input lines 
     * @return number of input lines that have a comment character at their beginning
     */
    public int getCommentLineCount() {
        commentLineCount
    }

    /**
     * As for now, just reports summary statistics to {@code stdout}
     */
    public void finish() {
        println "[${new Date()} ClonotypeStreamParser] Finished parsing $totalLines lines. " +
                "$commentLineCount header and $skippedLineCount bad line(s) were skipped. " +
                "Bad line statistics: ${badClonotypeStats.collect { it.key + "=" + it.value }.join(", ")}"
    }

    /**
     * INTERNAL checks if a given entry string is missing
     * @param entry
     * @return
     */
    private static boolean missingEntry(Sequence entry) {
        !entry || entry.size() == 0
    }

    private static boolean missingEntry(Segment segment) {
        segment == Segment.MISSING
    }

    /**
     * Implementation that simplifies usage syntax {@code parser.each{ Clonotype c -> ...}}
     * @return a clonotype iterator wrapped around the inner plain text table iterator
     */
    @Override
    public Iterator<Clonotype> iterator() {
        [hasNext: {
            innerIter.hasNext()
        }, next : {
            parse(innerIter.next())
        }] as Iterator
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
