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

package com.antigenomics.vdjtools.misc

import java.util.regex.Pattern
import java.util.zip.GZIPInputStream

/**
 * Class containing commonly used static functions for sequence manipulation and I/O
 */
public class CommonUtil {
    public static final String PLACEHOLDER = "."

    public static final String AA_LIST = /[FLSYCWPHQRIMTNKVADEGX\\*\\?]/ // Amino-acid sequence regex pattern

    public static final char[] NTS = ['A', 'T', 'G', 'C'], // list of allowed nucleotide characters (excluding N)
                               AAS = ['F', 'L', 'S', 'Y',  // list of allowed amino acid characters (excluding X and stop codon)
                                      'C', 'W', 'P', 'H',
                                      'Q', 'R', 'I', 'M',
                                      'T', 'N', 'K', 'V',
                                      'A', 'D', 'E', 'G']

    /**
     * Converts amino acid character to corresponding byte code
     * @param aa amino acid character (or stop codon), should be upper case
     * @return byte code in 0..20, or -1 if unknown character was provided
     */
    public static final byte aa2code(char aa) {
        switch (aa) {
            case 'F':
                return 0
            case 'L':
                return 1
            case 'S':
                return 2
            case 'Y':
                return 3
            case 'C':
                return 4
            case 'W':
                return 5
            case 'P':
                return 6
            case 'H':
                return 7
            case 'Q':
                return 8
            case 'R':
                return 9
            case 'I':
                return 10
            case 'M':
                return 11
            case 'T':
                return 12
            case 'N':
                return 13
            case 'K':
                return 14
            case 'V':
                return 15
            case 'A':
                return 16
            case 'D':
                return 17
            case 'E':
                return 18
            case 'G':
                return 19
            case '*':
                return 20
            default:
                return -1
        }
    }

    /**
     * Converts byte code to corresponding amino acid character
     * @param code byte code, should be in 0..20
     * @return amino acid character (or stop codon), or X for unknown byte code
     */
    public static final char code2aa(byte code) {
        switch (code) {
            case 0:
                return 'F'
            case 1:
                return 'L'
            case 2:
                return 'S'
            case 3:
                return 'Y'
            case 4:
                return 'C'
            case 5:
                return 'W'
            case 6:
                return 'P'
            case 7:
                return 'H'
            case 8:
                return 'Q'
            case 9:
                return 'R'
            case 10:
                return 'I'
            case 11:
                return 'M'
            case 12:
                return 'T'
            case 13:
                return 'N'
            case 14:
                return 'K'
            case 15:
                return 'V'
            case 16:
                return 'A'
            case 17:
                return 'D'
            case 18:
                return 'E'
            case 19:
                return 'G'
            case 20:
                return '*'
            default:
                return 'X'
        }
    }

    /**
     * Gets a byte code that corresponds to a given nucleotide
     * @param nt nucleotide, should be upper case
     * @return 0..3 for "A", "T", "G" or "C", -1 otherwise
     */
    public static final byte nt2code(char nt) {
        switch (nt) {
            case 'A':
                return 0
            case 'T':
                return 1
            case 'G':
                return 2
            case 'C':
                return 3
            default:
                return -1
        }
    }

    /**
     * Gets a nucleotide that corresponds to a given byte code
     * @param code nucleotide byte code
     * @return nucleotide (upper-case) for 0..3, "N" otherwise
     */
    public static final char code2nt(byte code) {
        switch (code) {
            case 0:
                return 'A'
            case 1:
                return 'T'
            case 2:
                return 'G'
            case 3:
                return 'C'
            default:
                return 'N'
        }
    }

    /**
     * Gets the complementary nucleotide (case-insensitive)
     * @param nt nucleotide
     * @return complementary nucleotide or N if a character other than A, T, G or C is passed
     */
    public static final char rcNt(char nt) {
        switch (nt) {
            case 'A':
                return 'T'
            case 'G':
                return 'C'
            case 'T':
                return 'A'
            case 'C':
                return 'G'
            case 'a':
                return 't'
            case 'g':
                return 'c'
            case 't':
                return 'a'
            case 'c':
                return 'g'
            default:
                return 'N'
        }
    }

    /**
     * Gets the reverse complement of a given nucleotide sequence
     * @param seq a nucleotide sequence
     * @return reverse complement, sequence on a complementary strand
     */
    public static String rc(String seq) {
        def chars = seq.reverse().toCharArray()
        for (int i = 0; i < chars.length; i++) {
            chars[i] = rcNt(chars[i])
        }
        return new String(chars)
    }

    /**
     * Translates a given codon 
     * @param codon a string of three nucleotides
     * @return corresponding amino acid, "*" for stop codon, "X" if codon contains undefined nucleotide "N", "?" for incomplete codon
     */
    public static char codon2aa(String codon) {
        String codonUpper = codon.toUpperCase()
        switch (codonUpper) {
            case 'TTT': return 'F'
            case 'TTC': return 'F'
            case 'TTA': return 'L'
            case 'TTG': return 'L'
            case 'TCT': return 'S'
            case 'TCC': return 'S'
            case 'TCA': return 'S'
            case 'TCG': return 'S'
            case 'TAT': return 'Y'
            case 'TAC': return 'Y'
            case 'TAA': return '*'
            case 'TAG': return '*'
            case 'TGT': return 'C'
            case 'TGC': return 'C'
            case 'TGA': return '*'
            case 'TGG': return 'W'
            case 'CTT': return 'L'
            case 'CTC': return 'L'
            case 'CTA': return 'L'
            case 'CTG': return 'L'
            case 'CCT': return 'P'
            case 'CCC': return 'P'
            case 'CCA': return 'P'
            case 'CCG': return 'P'
            case 'CAT': return 'H'
            case 'CAC': return 'H'
            case 'CAA': return 'Q'
            case 'CAG': return 'Q'
            case 'CGT': return 'R'
            case 'CGC': return 'R'
            case 'CGA': return 'R'
            case 'CGG': return 'R'
            case 'ATT': return 'I'
            case 'ATC': return 'I'
            case 'ATA': return 'I'
            case 'ATG': return 'M'
            case 'ACT': return 'T'
            case 'ACC': return 'T'
            case 'ACA': return 'T'
            case 'ACG': return 'T'
            case 'AAT': return 'N'
            case 'AAC': return 'N'
            case 'AAA': return 'K'
            case 'AAG': return 'K'
            case 'AGT': return 'S'
            case 'AGC': return 'S'
            case 'AGA': return 'R'
            case 'AGG': return 'R'
            case 'GTT': return 'V'
            case 'GTC': return 'V'
            case 'GTA': return 'V'
            case 'GTG': return 'V'
            case 'GCT': return 'A'
            case 'GCC': return 'A'
            case 'GCA': return 'A'
            case 'GCG': return 'A'
            case 'GAT': return 'D'
            case 'GAC': return 'D'
            case 'GAA': return 'E'
            case 'GAG': return 'E'
            case 'GGT': return 'G'
            case 'GGC': return 'G'
            case 'GGA': return 'G'
            case 'GGG': return 'G'
            default:
                if (codonUpper.contains("N") && codonUpper.length() == 3)
                    return "X" // undefined
                else
                    return '?' // incomplete
        }
    }

    /**
     * INTERNAL gets resource stream reader for a given resource
     * @param resourceName relative path to a given resource, e.g. "/imgt/segments.all.txt"
     * @return resource stream reader
     */
    public static InputStreamReader resourceStreamReader(String resourceName) {
        new InputStreamReader(resourceStream(resourceName))
    }

    public static InputStream resourceStream(String resourceName) {
        CommonUtil.class.classLoader.getResourceAsStream(resourceName)
    }

    /**
     * Gets file input stream for the specified file
     * @param fileName path to file, will assume gzipped file if ends with ".gz"
     * @return file input stream, or corresponding wrapper for gzipped files
     */
    public static InputStream getFileStream(String fileName) {
        fileName.endsWith(".gz") ? new GZIPInputStream(new FileInputStream(fileName)) : new FileInputStream(fileName)
    }

    public static String getSubSequence(String sequence, int from, int to) {
        String left = "", right = ""
        if (from < 0) {
            left = "N" * (-from)
            from = 0
        }
        if (to > sequence.length()) {
            right = "N" * (to - sequence.length())
            to = sequence.length()
        }

        left + sequence.substring(from, to) + right
    }

    /**
     * INTERNAL Bring V/D/J alleles to unified look
     * todo: consider allowing minor alleles
     */
    public static List<String> extractVDJ(List<String> vdj) {
        vdj.collect {
            def major = it.split(",")[0] // Taking only first variant in case of ties
            major = major.split("\\*")[0] // trim allele if present
            major = major.replaceAll("\"", "").trim() // zap characters introduced by file editing in external software (Excel,etc)
            major.length() > 0 ? major : PLACEHOLDER
        }
    }

    public static List<String> extractVDJImmunoSeq(List<String> vdj) {
        extractVDJ(vdj).collect {
            if (it.equalsIgnoreCase("unresolved")) {
                return PLACEHOLDER
            } else {
                // Immunoseq uses TCRBV instead of TRBV
                def res = it.replaceAll("TCR", "TR")
                // TCRBV12-03 -> TRBV12-3
                (1..9).each { int i ->
                    res = res.replaceAll("0$i".toString(), "$i".toString())
                }
                return res
            }
        }
    }

    public static List<String> extractVDJImmunoSeq(List<String> vdjFamilyTies,
                                                   List<String> vdjFamilies,
                                                   List<String> vdjGenes) {
        vdjFamilyTies = extractVDJImmunoSeq(vdjFamilyTies)
        vdjFamilies = extractVDJImmunoSeq(vdjFamilies)
        vdjGenes = extractVDJImmunoSeq(vdjGenes)
        def res = []
        vdjGenes.eachWithIndex { String gene, int i ->
            res << (gene == PLACEHOLDER ?
                    (vdjFamilies[i] == PLACEHOLDER ? vdjFamilyTies[i] : vdjFamilies[i]) : gene)
        }
        res
    }

    /*
     * CDR3 mapping & translation
     */

    final static String CYS_REGEX = /TG[TC]/,
                        PHE_TRP_REGEX = /(?:TGG|TT[TC])(?:GG[ATGC]|GC[ATGC])...GG[ATGC]/,
                        PHE_TRP_SHORT_REGEX = /(?:TGG|TT[TC])(?:GG[ATGC]|GC[ATGC])/

    private final static J_REF_PATTERN = Pattern.compile(PHE_TRP_REGEX),
                         J_REF_SHORT_PATTERN = Pattern.compile(PHE_TRP_SHORT_REGEX)

    static int getJReferencePoint(String sequence) {
        def matchers = [J_REF_PATTERN, J_REF_SHORT_PATTERN].collect { it.matcher(sequence) }
        def matcher = matchers.find { it.find() }
        matcher ? (matcher.start() - 1) : -1
    }

    final static String OOF_SYMBOLS_POSSIBLE = /([atgc#~_\?])+/, OOF_CHAR = "_", STOP_CHAR = "*"

    /**
     * Translates CDR3 sequence in V><J directions simultaneously, supports non-coding sequences
     * @param seq sequence to translate
     * @return amino acid sequence, as String object
     */
    static String translate(String seq) {
        if (seq.length() == 0)
            return ""

        def aaSeq = ""
        def oof = seq.size() % 3
        if (oof > 0) {
            def mid = (int) (seq.size() / 2)
            seq = seq.substring(0, mid) + ("?" * (3 - oof)) + seq.substring(mid, seq.length())
        }

        def leftEnd = -1, rightEnd = -1
        for (int i = 0; i <= seq.size() - 3; i += 3) {
            def codon = seq.substring(i, i + 3)
            if (codon.contains("?")) {
                leftEnd = i
                break
            }
            aaSeq += codon2aa(codon)
        }

        if (oof == 0)
            return aaSeq

        def aaRight = ""
        for (int i = seq.size(); i >= 3; i -= 3) {
            def codon = seq.substring(i - 3, i)
            if (codon.contains("?")) {
                rightEnd = i
                break
            }
            aaRight += codon2aa(codon)
        }

        return aaSeq + seq.substring(leftEnd, rightEnd).toLowerCase() + aaRight.reverse()
    }

    static String toUnifiedCdr3Aa(String seq) {
        seq.replaceAll(OOF_SYMBOLS_POSSIBLE, OOF_CHAR)
    }

    static boolean noStop(String seq) {
        !(seq.contains(STOP_CHAR))
    }

    static boolean inFrame(String seq) {
        !((boolean) (seq =~ OOF_SYMBOLS_POSSIBLE))
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
