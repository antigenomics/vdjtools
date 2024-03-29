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

package com.antigenomics.vdjtools.annotate

import com.antigenomics.vdjtools.misc.CommonUtil
import com.milaboratory.core.sequence.AminoAcidAlphabet
import com.milaboratory.core.sequence.AminoAcidSequence

class KnownAminoAcidProperties {
    private final Map<String, AaProperty> aminoAcidPropertiesByName = new HashMap<>()

    public static final KnownAminoAcidProperties INSTANCE = new KnownAminoAcidProperties()

    private KnownAminoAcidProperties() {
        // first load simple ones
        def lines = CommonUtil.resourceStream("profile/aa_property_table.txt").readLines().findAll {
            !it.startsWith("#")
        }

        def header = lines[0].split("[\t ]+")

        if (header[0] != "amino_acid") {
            throw new RuntimeException("Simple amino acid property file header should start with 'amino_acid'.")
        }

        def propertyNames = header[1..-1]

        int n = AminoAcidSequence.ALPHABET.size()

        def propertyValues = new float[propertyNames.size()][n]

        lines[1..-1].each { line ->
            def splitLine = line.split("[\t ]+")
            byte aa = AminoAcidSequence.ALPHABET.symbolToCode(splitLine[0].charAt(0))

            splitLine[1..-1].eachWithIndex { String it, int ind ->
                propertyValues[ind][aa] = it.toFloat()
            }
        }

        propertyNames.eachWithIndex { String name, int ind ->
            name = name.toLowerCase()
            aminoAcidPropertiesByName.put(name,
                    new SimpleAaProperty(name, propertyValues[ind]))
        }

        // load contact model
        loadContactProbs()
    }

    private void loadContactProbs() {
        def lines = CommonUtil.resourceStream("profile/cdr3contact.txt").readLines().collect { it.split("[ \t]") }
        def header = lines[0]

        def tcrChainColumn = header.findIndexOf { it == "tcr_chain" },
            posRelTcrColumn = header.findIndexOf { it == "pos_rel_tcr" },
            aaTcrColumn = header.findIndexOf { it == "aa_tcr" },
            pColumn = header.findIndexOf { it == "P" },
            posMaxColumn = header.findIndexOf { it == "pos_max" }

        int r = AminoAcidSequence.ALPHABET.size(), c = lines[1][posMaxColumn].toInteger()

        float[][] tcrAContactProbs = new float[r][c],
                  tcrBContactProbs = new float[r][c]

        lines[1..<lines.size()].each {
            int pos = it[posRelTcrColumn].toInteger(),
                aa = AminoAcidSequence.ALPHABET.symbolToCode(it[aaTcrColumn].charAt(0))
            float prob = (float) it[pColumn].toDouble()

            if (it[tcrChainColumn].equalsIgnoreCase("TRA")) {
                tcrAContactProbs[aa][pos] = prob
            } else if (it[tcrChainColumn].equalsIgnoreCase("TRB")) {
                tcrBContactProbs[aa][pos] = prob
            }
        }

        def contactEstimateProperty = new Cdr3ContactEstimate(tcrAContactProbs, tcrBContactProbs)

        aminoAcidPropertiesByName.put(contactEstimateProperty.name.toLowerCase(), contactEstimateProperty)
    }

    AaProperty getByName(String name) {
        name = name.toLowerCase()
        if (!aminoAcidPropertiesByName.containsKey(name))
            throw new IllegalArgumentException("Bad AaProperty name '$name', allowed values: ${allowedNames}")
        aminoAcidPropertiesByName[name]
    }

    List<String> getAllowedNames() {
        aminoAcidPropertiesByName.keySet().collect()
    }

    List<AaProperty> getAll() {
        aminoAcidPropertiesByName.values().collect()
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
