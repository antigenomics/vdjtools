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

import com.antigenomics.vdjtools.TestUtil
import com.antigenomics.vdjtools.misc.CommonUtil
import com.antigenomics.vdjtools.sample.Sample
import com.milaboratory.core.sequence.AminoAcidSequence
import org.junit.Test

class AnnotateTest {
    @Test
    void aaPropertyTest() {
        // Simple consistency test sum of values <> result of annotation of AA sequence that contains them all

        def seq1 = new AminoAcidSequence(CommonUtil.AAS.collect().join(""))
        def n = CommonUtil.AAS.size()

        KnownAminoAcidProperties.INSTANCE.getAll().
                findAll { it instanceof SimpleAaProperty }.each { SimpleAaProperty prop ->

            assert prop.values.collect().sum() == (0..<n).sum { prop.compute(seq1, it) }
        }
    }

    @Test
    void annotatorTest() {
        // Try all possible annotators

        def annotator = new SampleAnnotator(KnownAnnotators.INSTANCE.getAll())

        [TestUtil.DEFAULT_SAMPLE_COLLECTION, TestUtil.SINGLE_EMPTY_SAMPLE].each { samples ->
            samples.each { sample1 ->
                def sample = new Sample(sample1) // clone
                annotator.annotate(sample)

                assert sample.annotationHeader
                sample.each {
                    assert it.annotation
                }
            }
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
