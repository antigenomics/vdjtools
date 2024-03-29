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

package com.antigenomics.vdjtools.diversity

import org.apache.commons.math3.distribution.ZipfDistribution
import org.apache.commons.math3.random.Well19937c

class FrequencyTableGenerator {
    int numberOfSpecies = 500, observedSpecies = 100
    double exponent = 3.0

    FrequencyTableGenerator() {


    }

    FrequencyTableGenerator(int numberOfSpecies, double exponent) {
        this.numberOfSpecies = numberOfSpecies
        this.exponent = exponent
    }

    public FrequencyTable create() {
        def distr = new ZipfDistribution(new Well19937c(21051102L),
                numberOfSpecies, exponent)

        def cache = new HashMap<Long, Long>()

        observedSpecies.times {
            def count = distr.sample()
            cache.put(count, (cache[count] ?: 0) + 1)
        }

        new FrequencyTable(cache)
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
