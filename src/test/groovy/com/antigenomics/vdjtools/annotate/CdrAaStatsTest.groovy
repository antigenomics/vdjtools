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
import com.antigenomics.vdjtools.annotate.partitioning.Cdr3Center
import com.antigenomics.vdjtools.annotate.partitioning.FullCdr3
import com.antigenomics.vdjtools.sample.Clonotype
import com.antigenomics.vdjtools.misc.CommonUtil
import org.junit.Test

class CdrAaStatsTest {
    @Test
    void summaryTest() {
        // Test statistics

        def propertySummary = new AaPropertySummaryEvaluator(
                KnownAminoAcidProperties.INSTANCE.getByName("count"),
                new FullCdr3(),
                false,
                true
        )

        [TestUtil.DEFAULT_SAMPLE_COLLECTION, TestUtil.SINGLE_EMPTY_SAMPLE].each { samples ->
            samples.each { sample ->
                def summary = propertySummary.compute(sample)

                ["mean", "q25", "median", "q75"].each {
                    println "$it=${summary."$it"}"
                }
            }
        }
    }

    def clonotypes = [
            new Clonotype(null, 1, 1.0d,
                    [2, -1, -1, -10] as int[], "TRAV5", CommonUtil.PLACEHOLDER, "TRAJ48",
                    "TGTCATGAGAAATTAACCTTT",
                    "CHEKLTF", true, true, true),
            new Clonotype(null, 1, 1.0d,
                    [14, -1, -1, 19] as int[], "TRAV5", CommonUtil.PLACEHOLDER, "TRAJ48",
                    "TGCCTCGTGGGTGACTCGTACACGGGCAGGAGAGCACTTACTTTT",
                    "CLVGDSYTGRRALTF", true, true, true),
            new Clonotype(null, 1, 1.0d,
                    [-1, -1, -1, -1] as int[], "TRAV5", CommonUtil.PLACEHOLDER, "TRAJ48",
                    "A" * 45,
                    "A" * 6 + "F" * 3 + "A" * 6, true, true, true),
            new Clonotype(null, 1, 1.0d,
                    [-1, -1, -1, -1] as int[], "TRAV5", CommonUtil.PLACEHOLDER, "TRAJ48",
                    "A" * 45,
                    "F" * 6 + "A" * 3 + "F" * 6, true, true, true)]

    @Test
    void regionTest() {
        clonotypes.each { clonotype ->
            KnownCdr3Regions.INSTANCE.getAll().each { region ->
                println region.name + "\t" + region.extractAminoAcid(clonotype)
                println region.name + "\t" + region.extractNucleotide(clonotype)
            }
        }
    }

    @Test
    void  propertyTest() {
        def propertySummary = new AaPropertySummaryEvaluator(
                KnownAminoAcidProperties.INSTANCE.getByName("strength"),
                Cdr3Center.CDR3_CENTER_5,
                false,
                false
        )

        assert propertySummary.compute(clonotypes[2]) == 3
        assert propertySummary.compute(clonotypes[3]) == 2
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
