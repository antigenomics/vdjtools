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

package com.antigenomics.vdjtools;

/**
 * Something that has abundance data associated with it.
 */
public interface Countable {
    /**
     * Gets the number of variants associated with a given object.
     * It can be either number of convergent sub-variants associated
     * with a given clonotype or number of clonotypes in sample, or
     * the number of composite clonotypes that group several
     * clonotype sub-variants for joint/pooled samples.
     *
     * @return number of variants.
     */
    public int getDiversity();

    /**
     * Gets the number of reads associated with a given object.
     *
     * @return number of reads.
     */
    public long getCount();

    /**
     * Gets the share of reads associated with a given object.
     * Should return 1 for clonotype containers.
     *
     * @return share of reads.
     */
    public double getFreq();

    /**
     * Gets the non-normalized frequency (share of reads) for a given object.
     * For the {@link com.antigenomics.vdjtools.sample.Sample} the sum of
     * frequencies as specified in plain-text input file before normalization is returned.
     *
     * @return non-normalized share of reads.
     */
    public double getFreqAsInInput();
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
