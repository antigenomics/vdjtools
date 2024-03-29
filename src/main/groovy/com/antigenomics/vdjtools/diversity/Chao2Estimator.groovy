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

import com.antigenomics.vdjtools.join.JointSample

class Chao2Estimator {
    final JointSample jointSample

    Chao2Estimator(JointSample jointSample) {
        this.jointSample = jointSample
    }

    DiversityEstimate compute() {
        def sObs = jointSample.diversity,
            q1 = 0, q2 = 0, m = jointSample.numberOfSamples

        jointSample.each {
            def occurences = it.occurrences
            if (occurences == 1)
                q1++
            else if (occurences == 2)
                q2++
        }

        def q1q2 = q1 / q2

        new DiversityEstimate(sObs + (m - 1) / m * q1 * (q1 - 1) / 2 / (q2 + 1),
                Math.sqrt(q2 * (0.5 * Math.pow(q1q2, 2) + Math.pow(q1q2, 3) + 0.25 * Math.pow(q1q2, 4))),
                m)
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
