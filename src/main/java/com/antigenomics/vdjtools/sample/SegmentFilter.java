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

import java.util.Collections;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;

public abstract class SegmentFilter extends ClonotypeFilter {
    private int mySegmentSetSize = 0;
    private final String[] segmentNames;
    private final Set<String> segmentSet = new HashSet<>();

    public SegmentFilter(boolean negative, String... segmentNames) {
        super(negative);
        this.segmentNames = segmentNames;
    }

    public SegmentFilter(String... segmentNames) {
        this(false, segmentNames);
    }

    protected abstract String getSegmentName(Clonotype clonotype);

    private void refreshLazy() {
        if (mySegmentSetSize != SegmentFactory.INSTANCE.size()) {
            for (String name : segmentNames) {
                segmentSet.addAll(SegmentFactory.INSTANCE.getAtFuzzy(name)
                        .stream()
                        .map(Segment::getName)
                        .collect(Collectors.toList()));
            }
            mySegmentSetSize = SegmentFactory.INSTANCE.size();
        }
    }

    @Override
    protected boolean checkPass(Clonotype clonotype) {
        refreshLazy();
        return segmentSet.contains(getSegmentName(clonotype));
    }

    public Set<String> getSegmentSet() {
        return Collections.unmodifiableSet(segmentSet);
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
