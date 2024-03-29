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
package com.antigenomics.vdjtools.sample;

import com.antigenomics.vdjtools.misc.Segment;
import com.antigenomics.vdjtools.misc.SegmentFactory;

import java.util.HashMap;
import java.util.Map;

public class SegmentConverter implements ClonotypeConverter {
    private final Map<String, Segment> vSegmentMap = new HashMap<>(),
            jSegmentMap = new HashMap<>();

    public SegmentConverter(Map<String, String> vSegmentMap, Map<String, String> jSegmentMap) {
        for (Map.Entry<String, String> conv : vSegmentMap.entrySet()) {
            this.vSegmentMap.put(conv.getKey(),
                    SegmentFactory.INSTANCE.create(conv.getValue()));
        }

        for (Map.Entry<String, String> conv : jSegmentMap.entrySet()) {
            this.jSegmentMap.put(conv.getKey(),
                    SegmentFactory.INSTANCE.create(conv.getValue()));
        }
    }

    @Override
    public Clonotype convert(Clonotype clonotype) {
        return clonotype.withSegments(
                vSegmentMap.getOrDefault(clonotype.getV(), clonotype.getVBinary()),
                jSegmentMap.getOrDefault(clonotype.getJ(), clonotype.getJBinary()));
    }
}
