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


package com.antigenomics.vdjtools

import com.antigenomics.vdjtools.io.InputStreamFactory
import com.antigenomics.vdjtools.misc.Software
import com.antigenomics.vdjtools.sample.SampleCollection

import java.util.zip.GZIPInputStream

import static com.antigenomics.vdjtools.io.SampleStreamConnection.load

class TestUtil {
    static final SampleCollection DEFAULT_SAMPLE_COLLECTION = loadSamples(),
            SINGLE_EMPTY_SAMPLE = SampleCollection.fromSampleList([load(getResource("samples/empty.txt"), Software.VDJtools)])

    private static SampleCollection loadSamples() {
        def samples = Software.values().collect {
            load(getResource("samples/${it.toString().toLowerCase()}.txt.gz"), it)
        }
        samples.add(load(getResource("samples/empty.txt"), Software.VDJtools))

        SampleCollection.fromSampleList(samples)
    }


    public static InputStreamFactory getResource(String resourceName) {
        [
                create: {
                    def is = TestUtil.class.classLoader.getResourceAsStream(resourceName)
                    //print(resourceName)
                    resourceName.endsWith(".gz") ? new GZIPInputStream(is) : is
                },
                getId : { resourceName.split("/")[-1] }
        ] as InputStreamFactory
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
