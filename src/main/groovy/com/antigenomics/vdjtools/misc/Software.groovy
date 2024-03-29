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

/**
 * An enum that holds the list of supported Rep-seq processing software and parsing rules for their output.
 */
enum Software {
    MiTcr("\t", null, 2, false, false, ["count", "freq",
                                        "cdr3nt", _, _, "cdr3aa",
                                        _, "v", _, "j", _, "d",
                                        "VEnd", "DStart", "DEnd", "JStart"]),
    MiXcr("\t", null, 1, true, false, ["count", "freq",
                                       "cdr3nt", "cdr3aa",
                                       "v", "d", "j",
                                       "VEnd", "DStart", "DEnd", "JStart"]),
    MiGec("\t", null, 1, false, false, ["count", "freq",
                                        "cdr3nt", "cdr3aa",
                                        "v", "j", "d",
                                        "VEnd", "DStart", "DEnd", "JStart"]),
    MigMap("\t", null, 1, true, false, [_, _, "count", "freq",
                                        "cdr1nt", "cdr2nt", "cdr3nt", "cdr1aa", "cdr2aa", "cdr3aa",
                                        "inFrame", "noStop", "complete",
                                        4.times { _ }].flatten()),
    ImmunoSeq("\t", null, 1, true, true, ["cdr3nt", "cdr3aa", "count", "freq", "cdr3Length",
                                           _, "v", 5.times { _ },
                                           _, "d", 5.times { _ },
                                           _, "j", 5.times { _ },
                                           6.times { _ },
                                           _, "VEnd", "DStart", "DEnd", "JStart"].flatten()),
    ImmunoSeqV2("\t", null, 1, true, true, ["cdr3nt", "cdr3aa", "count", "freq", "cdr3Length",
                                           _, "v", 5.times { _ },
                                           _, "d", 5.times { _ },
                                           _, "j", 5.times { _ },
                                           6.times { _ },
                                           _, "VEnd", "DStart", "DEnd", "JStart"].flatten()),
    ImgtHighVQuest("\t", null, 1, true, true, [3.times { _ }, "v", "j", "d",
                                               9.times { _ }, "cdr3nt", 47.times { _ }, // comprehensive output
                                               "VEnd", 12.times { _ }, "DStart", "DEnd", 28.times {
        _
    }, "JStart"].flatten()),
    // Tricky for ImSEQ:
    // Output is not tab-delimited table, no freq, incomplete V/J names, no D, etc 'per read output'
    // is specified to re-calculate clonotype frequencies
    ImSeq("[\t:]", null, 0, true, true, ["v", "cdr3nt", "j", "count"]),
    Vidjil("\t", null, 0, true, true, []),
    // Same for RTCR - need to recompute frequencies
    RTCR("\t", null, 1, true, true, ["count", "cdr3aa", "v", "j",
                                       "cdr3nt", "VEnd", "JStart"]),
    VDJtools("\t", null, 1, false, false, ["count", "freq",
                                           "cdr3nt", "cdr3aa",
                                           "v", "d", "j",
                                           "VEnd", "DStart", "DEnd", "JStart"]),
    VDJtoolsRenorm("\t", null, 1, true, true, ["count", "freq",
                                           "cdr3nt", "cdr3aa",
                                           "v", "d", "j",
                                           "VEnd", "DStart", "DEnd", "JStart"])

    static final String _ = "blank"

    final String delimiter, comment
    final boolean collapseRequired, perReadOutput
    final int headerLineCount
    @Deprecated
    final List<String> printFields

    /**
     * Parsing rule for clonotype table produced by specified V-(D)-J mapping software.
     * @param delimiter delimited used in clonotype table.
     * @param comment comment character.
     * @param headerLineCount number of header lines.
     * @param collapseRequired if true, will additionally collapse clonotypes with identical V-CDR3nt-J (should be set to true if per read output is true).
     * @param perReadOutput if set to true, assumes per-read output and accumulates clonotype counts.
     * @param printFields fields that should be printed when writing a sample in specific format. (deprecated, default output is VDJtools_.
     */
    Software(String delimiter, String comment, int headerLineCount, boolean collapseRequired, boolean perReadOutput,
             List<String> printFields) {
        this.delimiter = delimiter
        this.comment = comment
        this.headerLineCount = headerLineCount
        this.collapseRequired = collapseRequired
        this.perReadOutput = perReadOutput
        this.printFields = printFields

        if (perReadOutput && !collapseRequired)
            throw new RuntimeException("Collapsing should always be performed for a software with per-read output")
    }

    /**
     * Gets a software by name, case-insensitive. 
     * @param name software name.
     * @return software type.
     */
    static Software byName(String name) {
        def software = values().find {
            it.toString().toLowerCase() == name.toLowerCase()
        }
        if (!software)
            throw new IllegalArgumentException("Unknown software $name")
        software
    }

    /**
     * List of allowed software names.
     */
    static String allowedNames = values().collect { it.toString().toLowerCase() }.join(",")
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
