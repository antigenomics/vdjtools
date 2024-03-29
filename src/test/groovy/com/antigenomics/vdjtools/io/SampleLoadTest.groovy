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

package com.antigenomics.vdjtools.io

import com.antigenomics.vdjtools.misc.Software
import org.junit.Test

import static Software.*
import static com.antigenomics.vdjtools.TestUtil.getResource
import static com.antigenomics.vdjtools.io.SampleStreamConnection.load

class SampleLoadTest {
    private static void loadTest(Software software, int count, int diversity, String suffix = '') {
        def resStream = getResource("samples/${software.toString().toLowerCase()}${suffix}.txt.gz")
        def sample = load(resStream, software)

        assert sample.count == count
        assert sample.diversity == diversity
        if (software.perReadOutput) {
            // check if frequency is recalculated
            assert Math.abs(sample.freqAsInInput - 1.0f) < 1e-5f
            assert Math.abs(sample.freq - 1.0f) < 1e-5f
        }

        assert sample.any { it.freq <= 1.0f }
    }

    @Test
    public void mitcrTest() {
        loadTest(MiTcr, 10000, 6493)
    }

    @Test
    public void migecTest() {
        loadTest(MiGec, 6147, 2420)
    }

    @Test
    public void migmapTest() {
        loadTest(MigMap, 720, 703)
    }

    @Test
    public void simpleTest() {
        loadTest(VDJtools, 11878, 2257)
    }

    @Test
    public void immunoseqTest() {
        loadTest(ImmunoSeq, 1108335, 42876)
    }

    @Test
    public void immunoseqTemplateTest() {
        loadTest(ImmunoSeq, 654, 1, "_1")
    }

    @Test
    public void immunoseqV2Test() {
        loadTest(ImmunoSeqV2, 1108335, 42876)
    }

    @Test
    public void imgtTest() {
        loadTest(ImgtHighVQuest, 9681, 7199)
    }

    @Test
    public void mixcrTest() {
        loadTest(MiXcr, 96132, 262)
    }

    @Test
    public void mixcrFullLengthTest() {
        loadTest(MiXcr, 14156, 20, ".fl")
    }

    @Test
    public void mixcr17LengthTest() {
        loadTest(MiXcr, 900, 873, ".171")
    }

    @Test
    public void mixcr21Test() {
        loadTest(MiXcr, 900, 873, ".21")
    }

    @Test
    public void mixcr3Test() {
        loadTest(MiXcr, 4509, 96, ".3")
    }

    @Test
    public void imseqTest() {
        loadTest(ImSeq, 647, 86)
    }


    @Test
    public void vidjilTest() {
        loadTest(Vidjil, 1667310, 165)
    }

    @Test
    public void rtcrTest() {
        loadTest(RTCR, 875, 693)
    }

    @Test
    public void emptyTest() {
        def resStream = getResource("samples/empty.txt")
        def sample = load(resStream)

        assert sample.count == 0
        assert sample.diversity == 0
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
