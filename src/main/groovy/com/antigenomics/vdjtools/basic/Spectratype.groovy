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
import com.antigenomics.vdjtools.sample.Sample

/**
 * Class that represents spectratype, i.e. distribution of clonotypes / clonotype frequency by 
 * CDR3 amino acid / nucleotide sequence length
 */
public class Spectratype extends Histogram {
    protected final boolean aminoAcid

    /**
     * Creates a spectratype instance with a given sample. All calculations are performed within constructor.
     * @param sample sample to initialize with
     * @param aminoAcid will use amino acid CDR3 sequence if true. Will use nucleotide sequence otherwise
     * @param unweighted will count each unique clonotype once if set to true. Will weight each clonotype by its frequency otherwise
     */
    public Spectratype(Sample sample, boolean aminoAcid, boolean unweighted) {
        this(aminoAcid, unweighted)
        addAll(sample)
    }

    /**
     * Creates a spectratype instance with a given sample. All calculations are performed within constructor.
     * Internal constructor 
     * @param sample sample sample to initialize with
     * @param intersectionType overlap type to deduce whether amino acid or nucleotide sequence CDR3 should be used
     * @param unweighted will count each unique clonotype once if set to true. Will weight each clonotype by its frequency otherwise
     */
    public Spectratype(Sample sample, OverlapType intersectionType, boolean unweighted) {
        this(intersectionType.aminoAcid, unweighted)
        addAll(sample)
    }

    /**
     * Creates a blank spectratype instance.
     * @param intersectionType overlap type to deduce whether amino acid or nucleotide sequence CDR3 should be used
     * @param unweighted will count each unique clonotype once if set to true. Will weight each clonotype by its frequency otherwise
     */
    public Spectratype(OverlapType intersectionType, boolean unweighted) {
        this(intersectionType.aminoAcid, unweighted)
    }

    /**
     * Creates a blank spectratype instance. 
     * @param aminoAcid will use amino acid CDR3 sequence if true. Will use nucleotide sequence otherwise
     * @param unweighted will count each unique clonotype once if set to true. Will weight each clonotype by its frequency otherwise
     */
    public Spectratype(boolean aminoAcid, boolean unweighted) {
        super(unweighted, aminoAcid ? 7 : 21, aminoAcid ? 23 : 69)
        this.aminoAcid = aminoAcid
    }

    /**
     * Update spectratype with a set of clonotypes, but omits several top clonotypes.
     * Internal, used for fancy (detalized) spectratype.
     * @param sample sample to add
     * @param top top clonotypes which will be excluded from histogram, and returned as list
     * @return list of top clonotypes
     */
    public List<Clonotype> addAllFancy(Iterable<Clonotype> sample, int top) {
        def topClonotypes = new LinkedList<Clonotype>();
        int counter = 0
        sample.each { Clonotype clonotype ->
            if (counter++ < top)
                topClonotypes.add(clonotype)
            else
                add(clonotype)
        }
        topClonotypes.sort { -it.freq }
    }

    @Override
    protected int getValue(Clonotype clonotype) {
        aminoAcid ? clonotype.cdr3nt.length() / 3 : clonotype.cdr3nt.length()
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
