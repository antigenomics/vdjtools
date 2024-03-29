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
package com.antigenomics.vdjtools.graph;

public class DegreeStatistics {
    private final int degree;
    private final long primaryGroupCount, secondaryGroupCount;

    public static final DegreeStatistics UNDEF = new DegreeStatistics(-1, -1, -1);

    public DegreeStatistics(int degree, long primaryGroupCount, long secondaryGroupCount) {
        this.degree = degree;
        this.primaryGroupCount = primaryGroupCount;
        this.secondaryGroupCount = secondaryGroupCount;
    }

    public int getDegree() {
        return degree;
    }

    public long getPrimaryGroupCount() {
        return primaryGroupCount;
    }

    public long getSecondaryGroupCount() {
        return secondaryGroupCount;
    }

    @Override
    public String toString() {
        return degree + "\t" + primaryGroupCount + "\t" + secondaryGroupCount;
    }
}
