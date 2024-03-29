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

package com.antigenomics.vdjtools.misc;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * An immune receptor segment factory, acting as a cache and ensuring
 * compatibility of segment annotations. This class is a singleton.
 */
public class SegmentFactory {
    public static final SegmentFactory INSTANCE = new SegmentFactory();

    protected final Map<String, Segment> segmentCache = new HashMap<>();

    private SegmentFactory() {
        segmentCache.put(Segment.MISSING.name, Segment.MISSING);
    }

    /**
     * Creates a new segment by identifier or returns corresponding segment object if it already exists.
     *
     * @param name segment identifier.
     * @return segment.
     */
    public Segment create(String name) {
        Segment segment = segmentCache.get(name);

        if (segment == null) {
            segmentCache.put(name, segment = new Segment(name));
        }

        return segment;
    }

    /**
     * Gets the size of segment cache.
     *
     * @return number of segments.
     */
    public int size() {
        return segmentCache.size();
    }

    /**
     * Gets segment by identifier.
     *
     * @param name segment identifier.
     * @return segment.
     */
    public Segment getAt(String name) {
        return segmentCache.get(name);
    }

    public List<Segment> getAtFuzzy(String namePart) {
        return segmentCache.entrySet()
                .stream()
                .filter(segmentEntry -> segmentEntry.getKey().startsWith(namePart))
                .map(Map.Entry::getValue)
                .collect(Collectors.toList());
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
