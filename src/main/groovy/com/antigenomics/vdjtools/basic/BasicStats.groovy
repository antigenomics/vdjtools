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
import com.antigenomics.vdjtools.join.ClonotypeKeyGen
import com.antigenomics.vdjtools.join.key.ClonotypeKey
import com.antigenomics.vdjtools.sample.Sample
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics

/**
 * This class provides some basic stats that could be computed for a RepSeq sample:
 * read count, observed diversity, characteristic clone size, cdr3/N-N insert/N-D-N sequence size,
 * fraction of non-coding clonotypes and repertoire convergence 
 */
public class BasicStats {
    private final DescriptiveStatistics cloneSize, cloneSizeGeom
    private final double cdr3ntLength, insertSize, ndnSize
    private int ncDiversity
    private double ncFrequency
    private final long count
    private final int diversity
    private final double convergence
    private final boolean weighted

    /**
     * Create an instance of BasicStats class. All computations are done within constructor.
     * @param sample sample to be analyzed
     */
    public BasicStats(Sample sample) {
        this(sample, true)
    }

    /**
     * Create an instance of BasicStats class. All computations are done within constructor.
     * @param sample sample to be analyzed
     * @param weighted if set to {@code true}, will use clonotype frequency to weight mean_cdr3nt_length, mean_insert_size and mean_ndn_size
     */
    public BasicStats(Sample sample, boolean weighted) {
        this.count = sample.count
        this.diversity = sample.diversity

        this.cloneSize = new DescriptiveStatistics()
        this.cloneSizeGeom = new DescriptiveStatistics()

        double cdr3ntLength = 0, insertSize = 0, ndnSize = 0

        this.weighted = weighted

        def keyGenAA = new ClonotypeKeyGen(OverlapType.AminoAcid),
            keyGenNT = new ClonotypeKeyGen(OverlapType.Nucleotide)
        def aaSet = new HashSet<ClonotypeKey>(),
            ntSet = new HashSet<ClonotypeKey>()

        int weight = 1
        double denom = 0

        sample.each {
            cloneSize.addValue(it.freq)
            cloneSizeGeom.addValue(Math.log10(it.freq))

            if (weighted)
                weight = it.count

            cdr3ntLength += weight * it.cdr3nt.length()

            def x = it.insertSize
            if (x > -1)
                insertSize += weight * x

            x = it.NDNSize
            if (x > -1)
                ndnSize += weight * x

            if (!it.coding) {
                ncDiversity++
                ncFrequency += it.freq
            }

            aaSet.add(keyGenAA.generateKey(it))
            ntSet.add(keyGenNT.generateKey(it))

            denom += weight
        }

        this.convergence = ntSet.size() / (double) aaSet.size()
        this.cdr3ntLength = cdr3ntLength / denom
        this.ndnSize = ndnSize / denom
        this.insertSize = insertSize / denom
    }

    /**
     * Gets mean clonotype frequency 
     * @return
     */
    public double getMeanFrequency() {
        cloneSize.mean
    }

    /**
     * Gets geometric mean of clonotype frequency 
     * @return
     */
    public double getGeomMeanFrequency() {
        Math.pow(10, cloneSizeGeom.mean)
    }

    /**
     * Gets mean CDR3 nucleotide sequence length. 
     * If {@code weighted = true} the mean will be weighted by clonotype count 
     * @return
     */
    public double getMeanCdr3ntLength() {
        cdr3ntLength
    }

    /**
     * Gets mean NDN region size, i.e. mean number of nucleotides between V segment end and J segment start 
     * If {@code weighted = true} the mean will be weighted by clonotype count
     * @return
     */
    public double getMeanNDNSize() {
        ndnSize
    }

    /**
     * Gets mean insert size, i.e. mean number of nucleotides in V-D and D-J regions if D segment is determined or
     * V-J regions if not. If {@code weighted = true} the mean will be weighted by clonotype count
     * @return
     */
    public double getMeanInsertSize() {
        insertSize
    }

    /**
     * Gets the number of reads in sample 
     * @return
     */
    public long getCount() {
        count
    }

    /**
     * Gets the total number of clonotypes in sample
     * @return
     */
    public int getDiversity() {
        diversity
    }

    /**
     * Gets the total number of non-coding (containing a stop or frameshift) clonotypes in sample 
     * @return
     */
    public int getNcDiversity() {
        ncDiversity
    }

    /**
     * Gets the cumulative frequency of non-coding (containing a stop or frameshift) clonotypes in sample 
     * @return
     */
    public double getNcFrequency() {
        ncFrequency
    }

    /**
     * Gets the convergence level of sample, i.e. number of CDR3 nucleotide sequences per CDR3 amino acid sequence 
     * @return
     */
    public double getConvergence() {
        convergence
    }

    /**
     * Header used for tabular output 
     */
    public static final String HEADER = "count\tdiversity\t" +
            "mean_frequency\tgeomean_frequency\t" +
            "nc_diversity\tnc_frequency\t" +
            "mean_cdr3nt_length\tmean_insert_size\tmean_ndn_size\t" +
            "convergence"

    /**
     * Plain text row for tabular output 
     */
    @Override
    public String toString() {
        [count, diversity,
         meanFrequency, geomMeanFrequency,
         ncDiversity, ncFrequency,
         meanCdr3ntLength, meanInsertSize, meanNDNSize,
         convergence
        ].flatten().join("\t")
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
