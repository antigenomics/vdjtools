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

class KnownAnnotators {
    private final Map<String, ClonotypeAnnotator> annotatorsByName = new HashMap<>()

    public static final KnownAnnotators INSTANCE = new KnownAnnotators()

    private KnownAnnotators() {
        ["cdr3Length", "NDNSize", "insertSize", "VDIns", "DJIns"].each {
            annotatorsByName.put(it.toLowerCase(), new BaseAnnotator(it))
        }

        KnownAminoAcidProperties.INSTANCE.allowedNames.each {
            def annotator = new AAPropertyAnnotator(it, true)
            annotatorsByName.put(annotator.name.toLowerCase(), annotator)
            annotator = new AAPropertyAnnotator(it, false)
            annotatorsByName.put(annotator.name.toLowerCase(), annotator)
        }
    }


    ClonotypeAnnotator getByName(String name) {
        name = name.toLowerCase()
        if (!annotatorsByName.containsKey(name))
            throw new IllegalArgumentException("Bad annotator name '$name', allowed values: ${allowedNames}")
        annotatorsByName[name]
    }

    List<String> getAllowedNames() {
        annotatorsByName.keySet().collect()
    }

    List<ClonotypeAnnotator> getAll() {
        annotatorsByName.values().collect()
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
