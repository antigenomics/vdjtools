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

package com.antigenomics.vdjtools.basic

import com.antigenomics.vdjtools.overlap.OverlapType
import com.antigenomics.vdjtools.sample.Clonotype

/**
 * Class that represents spectratype collection grouped by Variable segments,
 * i.e. distribution of clonotypes / clonotype frequency with a given Variable segment by
 * CDR3 amino acid / nucleotide sequence length
 */
public class SpectratypeV {
    private final Map<String, Spectratype> spectratypes = new HashMap<>()
    private final Spectratype dummy

    /**
     * Creates a blank Variable segment spectratype instance.
     */
    public SpectratypeV() {
        this(false, false)
    }

    /**
     * Creates a blank Variable segment spectratype instance.
     * @param intersectionType overlap type to deduce whether amino acid or nucleotide sequence CDR3 should be used to determine the histogram bin
     */
    public SpectratypeV(OverlapType intersectionType) {
        this(intersectionType.aminoAcid, false)
    }

    /**
     * Creates a blank Variable segment spectratype instance.
     * @param intersectionType overlap type to deduce whether amino acid or nucleotide sequence CDR3 should be used to determine the histogram bin
     * @param unweighted will count each unique clonotype once if set to true. Will weight each clonotype by its frequency otherwise
     */
    public SpectratypeV(OverlapType intersectionType, boolean unweighted) {
        this(intersectionType.aminoAcid, unweighted)
    }

    /**
     * Creates a blank Variable segment spectratype instance.
     * @param aminoAcid will use amino acid CDR3 sequence to determine the histogram bin if true. Will use nucleotide sequence otherwise
     * @param unweighted will count each unique clonotype once if set to true. Will weight each clonotype by its frequency otherwise
     */
    public SpectratypeV(boolean aminoAcid, boolean unweighted) {
        this.dummy = new Spectratype(aminoAcid, unweighted)
    }

    /**
     * Clears all Variable segment spectratype counters 
     */
    public void clear() {
        spectratypes.values().each { it.clear() }
    }

    /**
     * Update Variable segment spectratype with a set of clonotypes 
     * @param sample sample to add
     */
    public void addAll(Iterable<Clonotype> sample) {
        sample.each { clonotype ->
            def vSpectra = spectratypes[clonotype.v]
            if (!vSpectra)
                spectratypes.put(clonotype.v, vSpectra = new Spectratype(dummy.aminoAcid, dummy.unweighted))
            vSpectra.add(clonotype)
        }
    }

    /**
     * Gets the current Variable segment list 
     * @return
     */
    public List<String> vSegmentList() {
        spectratypes.sort { -it.value.freq }.collect { it.key }
    }

    /**
     * Gets a spectratype for a given Variable segment
     * @param vSegmentName Variable segment for sub-setting
     * @return
     */
    public Spectratype getAt(String vSegmentName) {
        spectratypes[vSegmentName]
    }

    /**
     * Selects spectratypes for {@code 0..top-1} Variable segments. 
     * Note that you should not normalize the resulting histograms if you want to see the correct picture,
     * so call {@code spectratype.getHistogram ( false )}
     * @param top number of most frequent Variable segments to take
     * @return map with Variable segment name {@code String} as key and {@code Spectratype} as value
     */
    public Map<String, Spectratype> collapse(int top) {
        def collapsedSpectratypes = new HashMap<String, Spectratype>()
        def otherSpectratype
        collapsedSpectratypes.put("other", otherSpectratype = new Spectratype(dummy.aminoAcid, dummy.unweighted))

        spectratypes.sort { -it.value.freq }.eachWithIndex { it, ind ->
            if (ind < top)
                collapsedSpectratypes.put(it.key, it.value)
            else
                otherSpectratype.addAll(it.value)
        }

        collapsedSpectratypes.sort { -it.value.freq }
    }

    /**
     * Gets the number of bins in a spectratype histogram. 
     * It is the same for all Variable segments
     * @return
     */
    public int getSpan() {
        dummy.span
    }

    /**
     * Gets an array of CDR3 lengths that correspond to spectratype bins.
     * It is the same for all Variable segments
     */
    public int[] getLengths() {
        dummy.lengths
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
