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

import com.antigenomics.vdjtools.sample.Chain;

/**
 * An immune receptor segment.
 */
public class Segment {
    public static final Segment MISSING = new Segment(CommonUtil.PLACEHOLDER);

    protected final String name;

    Segment(String name) {
        this.name = name;
    }

    /**
     * Gets segment identifier.
     *
     * @return segment identifier.
     */
    public String getName() {
        return name;
    }

    @Override
    public String toString() {
        return name;
    }

    public Chain getChain() {
        if (name.length() < 3) {
            return Chain.NA;
        }

        switch (name.substring(0, 3).toUpperCase()) {
            case "TRA":
                return Chain.TRA;
            case "TRB":
                return Chain.TRB;
            case "IGH":
                return Chain.IGH;
            case "IGL":
                return Chain.IGL;
            case "IGK":
                return Chain.IGK;
            case "TRG":
                return Chain.TRG;
            case "TRD":
                return Chain.TRD;
        }

        return Chain.NA;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Segment segment = (Segment) o;

        if (!name.equals(segment.name)) return false;

        return true;
    }

    @Override
    public int hashCode() {
        return name.hashCode();
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
