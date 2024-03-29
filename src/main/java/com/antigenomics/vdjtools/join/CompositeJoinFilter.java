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

package com.antigenomics.vdjtools.join;

/**
 * A set of filters all of which should be passed by a given joint clonotype.
 */
public class CompositeJoinFilter implements JoinFilter {
    private final JoinFilter[] filters;

    /**
     * Creates a set of joint that all should be passed.
     *
     * @param filters joint clonotype filters.
     */
    public CompositeJoinFilter(JoinFilter... filters) {
        this.filters = filters;
    }

    /**
     * Checks whether a given joint clonotype passes all specified filters.
     *
     * @param jointClonotype a joint clonotype to check.
     * @return true if all filters are passed and joint clonotype should be retained, false otherwise.
     */
    @Override
    public boolean pass(JointClonotype jointClonotype) {
        for (JoinFilter filter : filters)
            if (!filter.pass(jointClonotype))
                return false;

        return true;
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
