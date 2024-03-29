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
package com.antigenomics.vdjtools.basic

import com.antigenomics.vdjtools.TestUtil
import org.junit.Test

class FancySpectratypeTest {
    @Test
    void test0() {
        [TestUtil.DEFAULT_SAMPLE_COLLECTION, TestUtil.SINGLE_EMPTY_SAMPLE].each { samples ->
            samples.each {
                if (it.diversity > 0) {
                    def fancySpectratype = new FancySpectratype(it, 10)

                    def values = fancySpectratype.spectraMatrix.collect { it.toList() }.flatten()

                    assert values.every { it >= 0 }
                    assert Math.abs(values.sum() - 1.0) < 1e-8
                } else {
                    // Test no error for empty
                    assert new FancySpectratype(TestUtil.SINGLE_EMPTY_SAMPLE.first(), 10).
                            spectraMatrix.collect { it.toList() }.flatten().every { it == 0 }
                }
            }
        }
    }
}
