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
 * Filters joint clonotypes according to the number of samples they were detected.
 */
public class OccurrenceJoinFilter implements JoinFilter {
    private final int occurrenceThreshold;

    /**
     * Creates a filter that retains all joint clonotypes detected two or more times.
     */
    public OccurrenceJoinFilter() {
        this(2);
    }

    /**
     * Creates a filter that retains all joint clonotypes detected the specified number of times or more.
     *
     * @param occurrenceThreshold threshold for the number of samples this clonotype was detected (inclusive).
     */
    public OccurrenceJoinFilter(int occurrenceThreshold) {
        this.occurrenceThreshold = occurrenceThreshold;
    }

    public int getOccurrenceThreshold() {
        return occurrenceThreshold;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean pass(JointClonotype jointClonotype) {
        int detectionCounter = 0;
        for (int i = 0; i < jointClonotype.getParent().getNumberOfSamples(); i++) {
            if (jointClonotype.present(i) && ++detectionCounter == occurrenceThreshold) return true;
        }
        return false;
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
