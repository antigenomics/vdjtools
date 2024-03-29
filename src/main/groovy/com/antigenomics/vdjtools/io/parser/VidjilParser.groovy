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
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.metadata.MetadataEntry
import groovy.json.JsonSlurper

class VidjilParser extends BaseParser {
    public static final String VIDJIL_SAMPLE_ID_COL = "vidjil_id"

    /**
     * {@inheritDoc}
     */
    protected VidjilParser(Iterator<String> innerIter, Sample sample) {
        super(jsonToTabular(innerIter, sample), Software.Vidjil, sample)
    }

    /* Only used fields shown here
    {
    "clones": [
    {
      "germline": "IGH",
      "id": "CACGGCCTTGTATTACTGTGCACCCGGAGGTATGGACGTCTGGGGCCAAG",
      "name": "IGHV3-9*01 7/CCCGGA/17 IGHJ6*02",
      "reads": [ 189991 ],
      "seg": {
        "3": "IGHJ6*02",
        "3del": 17,
        "3start": 282,
        "4": "IGHD6-13*01",
        "4end": 231,
        "4start": 220,
        "5": "IGHV3-9*01",
        "5del": 7,
        "5end": 275,
        "N": 6,
        "cdr3": {
          "aa": "APGGMDV",
          "start": 274,
          "stop": 294
        },
        "junction": {
          "aa": "CAPGGMDVW",
          "productive": true,
          "start": 271,
          "stop": 297
        }
      },
      "seg_stat": {
        "2": 55411,
        "3": 134580
      },
      "sequence": "GGAGTCGGGGGAGGCTTGGTACAGCCTGGCAGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTTGATGATTATGCCATGCACTGGGTCCGGCAAGCTCCAGGGAAGGGCCTGGAGTGGGTCTCAGGTATTAGTTGGAATAGTGGTAGCATAGGCTATGCGGACTCTGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCCCTGTATCTGCAAATGAACAGTCTGAGAGCTGAGGACACGGCCTTGTATTACTGTGCACCCGGAGGTATGGACGTCTGGGGCCAAGGGACCCTGGTCACC",
    },
    ...
     */

    private static Iterator<String> jsonToTabular(Iterator<String> innerIter, Sample sample) {
        def jsonSb = new StringBuilder()

        while (innerIter.hasNext()) {
            jsonSb.append(innerIter.next()).append("\n")
        }

        def jsonClones = new JsonSlurper().parseText(jsonSb.toString()).clones

        def tabulatedOutput = new ArrayList<String>(jsonClones.size())

        def sampleId = (sample.sampleMetadata[VIDJIL_SAMPLE_ID_COL] ?:
                new MetadataEntry(null, null, null, "0")).value.toInteger()

        jsonClones.each { cloneObj ->
            def segmentationInfo = cloneObj.seg,
                junctionInfo = segmentationInfo.junction

            if (segmentationInfo && junctionInfo) {
                tabulatedOutput << [
                        cloneObj.reads[sampleId],                                                     // count
                        0,                                                                            // freq
                        cloneObj.sequence[(junctionInfo.start - 1)..<junctionInfo.stop],              // CDR3nt
                        junctionInfo.aa,                                                              // CDR3aa
                        segmentationInfo."5",                                                         // V
                        segmentationInfo."4" ?: ".",                                                  // D
                        segmentationInfo."3",                                                         // J
                        segmentationInfo."5end" - junctionInfo.start,                                 // Vend
                        (segmentationInfo."4start" ?: (junctionInfo.start - 1)) - junctionInfo.start, // Dstart
                        (segmentationInfo."4end" ?: (junctionInfo.start - 1)) - junctionInfo.start,   // Dend
                        segmentationInfo."3start" - junctionInfo.start                                // Jstart
                ].join("\t")
            }
        }

        tabulatedOutput.iterator()
    }

    @Override
    String getAnnotationHeader() {
        null
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
