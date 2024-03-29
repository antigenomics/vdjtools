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

import com.antigenomics.vdjtools.misc.ExecUtil
import com.antigenomics.vdjtools.sample.Clonotype
import com.antigenomics.vdjtools.sample.Sample
import groovyx.gpars.GParsPool

class SampleAnnotator {
    final List<ClonotypeAnnotator> annotators

    SampleAnnotator(List<ClonotypeAnnotator> annotators) {
        this.annotators = annotators
    }

    void annotate(Sample sample) {
        if (sample.annotationHeader) {
            sample.annotationHeader += "\t" + annotators.collect { it.category + "." + it.name }.join("\t")
            GParsPool.withPool ExecUtil.THREADS, {
                sample.eachParallel { Clonotype clonotype ->
                    clonotype.annotation += "\t" + annotators.collect { it.annotate(clonotype) }.join("\t")
                }
            }
        } else {
            sample.annotationHeader = annotators.collect { it.category + "." + it.name }.join("\t")
            GParsPool.withPool ExecUtil.THREADS, {
                sample.eachParallel { Clonotype clonotype ->
                    clonotype.annotation = annotators.collect { it.annotate(clonotype) }.join("\t")
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
