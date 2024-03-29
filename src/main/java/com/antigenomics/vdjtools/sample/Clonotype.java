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

import com.antigenomics.vdjtools.ClonotypeWrapper;
import com.antigenomics.vdjtools.Countable;
import com.antigenomics.vdjtools.misc.CommonUtil;
import com.antigenomics.vdjtools.misc.Segment;
import com.antigenomics.vdjtools.misc.SegmentFactory;
import com.milaboratory.core.sequence.AminoAcidSequence;
import com.milaboratory.core.sequence.NucleotideSequence;

import java.util.Arrays;
import java.util.List;

/**
 * A class holding comprehensive info on a T- or B-cell clonotype.
 * CDR stands for Complementarity Determining Region
 */
public class Clonotype implements Comparable<Clonotype>, Countable, ClonotypeWrapper {
    private Sample parent;
    private int count;
    private double freq;

    private String annotation;

    private final int[] segmPoints;
    private final Segment v, d, j;
    private final NucleotideSequence cdr3nt;
    private final AminoAcidSequence cdr3aa;

    private final boolean inFrame, isComplete, noStop;

    /**
     * Creates a new clonotype explicitly setting all fields
     *
     * @param parent     parent sample
     * @param count      clonotype count, number of reads associated with this clonotype
     * @param freq       clonotype frequency, fraction of reads associated with this clonotype in the parent sample
     * @param segmPoints an array containing Variable segment end, Diversity segment start, Diversity segment end and Joining segment start in CDR3 coordinates
     * @param v          Variable segment identifier
     * @param d          Diversity segment identifier
     * @param j          Joining segment identifier
     * @param cdr3nt     nucleotide sequence of CDR3
     * @param cdr3aa     amino acid sequence of CDR3
     * @param inFrame    tells if this clonotype is in frame, i.e. it's sequence doesn't contain frameshifts
     * @param noStop     tells if the sequence of this clonotype doesn't contain any stop codon
     * @param isComplete tells if this clonotype is complete, i.e. CDR3 sequence for this clonotype is completely determined
     * @param annotation clonotype annotation string (can be null)
     */
    public Clonotype(Sample parent, int count, double freq,
                     int[] segmPoints, Segment v, Segment d, Segment j,
                     NucleotideSequence cdr3nt, AminoAcidSequence cdr3aa,
                     boolean inFrame, boolean noStop, boolean isComplete,
                     String annotation) {
        this.parent = parent;
        this.count = count;
        this.freq = freq;
        this.segmPoints = segmPoints;
        this.v = v;
        this.d = d;
        this.j = j;
        this.cdr3nt = cdr3nt;
        this.cdr3aa = cdr3aa;
        this.inFrame = inFrame;
        this.isComplete = isComplete;
        this.noStop = noStop;
        this.annotation = annotation;
    }

    /**
     * Creates a new clonotype explicitly setting all fields
     *
     * @param parent     parent sample
     * @param count      clonotype count, number of reads associated with this clonotype
     * @param freq       clonotype frequency, fraction of reads associated with this clonotype in the parent sample
     * @param segmPoints an array containing Variable segment end, Diversity segment start, Diversity segment end and Joining segment start in CDR3 coordinates
     * @param v          Variable segment identifier
     * @param d          Diversity segment identifier
     * @param j          Joining segment identifier
     * @param cdr3nt     nucleotide sequence of CDR3
     * @param cdr3aa     amino acid sequence of CDR3
     * @param inFrame    tells if this clonotype is in frame, i.e. it's sequence doesn't contain frameshifts
     * @param noStop     tells if the sequence of this clonotype doesn't contain any stop codon
     * @param isComplete tells if this clonotype is complete, i.e. CDR3 sequence for this clonotype is completely determined
     */
    public Clonotype(Sample parent, int count, double freq,
                     int[] segmPoints, Segment v, Segment d, Segment j,
                     NucleotideSequence cdr3nt, AminoAcidSequence cdr3aa,
                     boolean inFrame, boolean noStop, boolean isComplete) {
        this(parent, count, freq, segmPoints, v, d, j, cdr3nt, cdr3aa, inFrame, noStop, isComplete, null);
    }

    /**
     * Creates a new clonotype explicitly setting all fields
     *
     * @param parent     parent sample
     * @param count      clonotype count, number of reads associated with this clonotype
     * @param freq       clonotype frequency, fraction of reads associated with this clonotype in the parent sample
     * @param segmPoints an array containing Variable segment end, Diversity segment start, Diversity segment end and Joining segment start in CDR3 coordinates
     * @param v          Variable segment identifier
     * @param d          Diversity segment identifier
     * @param j          Joining segment identifier
     * @param cdr3nt     nucleotide sequence of CDR3
     * @param cdr3aa     amino acid sequence of CDR3
     * @param inFrame    tells if this clonotype is in frame, i.e. it's sequence doesn't contain frameshifts
     * @param noStop     tells if the sequence of this clonotype doesn't contain any stop codon
     * @param isComplete tells if this clonotype is complete, i.e. CDR3 sequence for this clonotype is completely determined
     * @param annotation clonotype annotation string (can be null)
     */
    public Clonotype(Sample parent, int count, double freq,
                     int[] segmPoints, String v, String d, String j,
                     String cdr3nt, String cdr3aa,
                     boolean inFrame, boolean noStop, boolean isComplete,
                     String annotation) {
        this(parent, count, freq, segmPoints,
                SegmentFactory.INSTANCE.create(v), SegmentFactory.INSTANCE.create(d), SegmentFactory.INSTANCE.create(j),
                new NucleotideSequence(cdr3nt), new AminoAcidSequence(cdr3aa),
                inFrame, noStop, isComplete, annotation);
    }

    /**
     * Creates a new clonotype explicitly setting all fields
     *
     * @param parent     parent sample
     * @param count      clonotype count, number of reads associated with this clonotype
     * @param freq       clonotype frequency, fraction of reads associated with this clonotype in the parent sample
     * @param segmPoints an array containing Variable segment end, Diversity segment start, Diversity segment end and Joining segment start in CDR3 coordinates
     * @param v          Variable segment identifier
     * @param d          Diversity segment identifier
     * @param j          Joining segment identifier
     * @param cdr3nt     nucleotide sequence of CDR3
     * @param cdr3aa     amino acid sequence of CDR3
     * @param inFrame    tells if this clonotype is in frame, i.e. it's sequence doesn't contain frameshifts
     * @param noStop     tells if the sequence of this clonotype doesn't contain any stop codon
     * @param isComplete tells if this clonotype is complete, i.e. CDR3 sequence for this clonotype is completely determined
     */
    public Clonotype(Sample parent, int count, double freq,
                     int[] segmPoints, String v, String d, String j,
                     String cdr3nt, String cdr3aa,
                     boolean inFrame, boolean noStop, boolean isComplete) {
        this(parent, count, freq, segmPoints, v, d, j, cdr3nt, cdr3aa, inFrame, noStop, isComplete, null);
    }

    /**
     * A copying constructor
     *
     * @param toCopy clonotype to copy all fields from
     */
    public Clonotype(Clonotype toCopy) {
        this(toCopy, toCopy.parent, toCopy.count);
    }

    /**
     * A copying constructor
     *
     * @param toCopy    clonotype to copy all fields from
     * @param newParent sample this clonotype will be assigned to
     */
    public Clonotype(Clonotype toCopy, Sample newParent) {
        this(toCopy, newParent, toCopy.count);
    }

    /**
     * A copying constructor
     *
     * @param toCopy    clonotype to copy all fields from
     * @param newParent sample this clonotype will be assigned to
     * @param newCount  number of reads that will be associated with this clonotype in new sample
     */
    public Clonotype(Clonotype toCopy, Sample newParent, int newCount) {
        this(newParent, newCount, toCopy.freq,
                toCopy.segmPoints, toCopy.v, toCopy.d, toCopy.j,
                toCopy.cdr3nt, toCopy.cdr3aa,
                toCopy.inFrame, toCopy.noStop, toCopy.isComplete,
                toCopy.annotation);
    }

    /**
     * Gets the annotation string associated with this clonotype.
     *
     * @return clonotype annotation or null if not available
     */
    public String getAnnotation() {
        return annotation;
    }

    /**
     * Sets the annotation string associated with this clonotype.
     */
    public void setAnnotation(String annotation) {
        this.annotation = annotation;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getDiversity() {
        return 1;
    }

    /**
     * Gets the clonotype count
     *
     * @return number of reads that are associated with this clonotype
     */
    public long getCount() {
        return count;
    }

    /**
     * Gets the frequency of this clonotype as specified in plain-text input file
     *
     * @return fraction of reads that are associated with this clonotype according to plain-text input file
     */
    public double getFreqAsInInput() {
        return freq;
    }

    /**
     * Gets the frequency of this clonotype in its parent sample
     *
     * @return the fraction of reads that are associated with this clonotype in its parent sample
     */
    public double getFreq() {
        return count / (double) parent.getCount();
    }

    /**
     * Gets the Variable segment of this clonotype
     *
     * @return string identifier of Variable segment
     */
    public String getV() {
        return v.toString();
    }

    public Segment getVBinary() {
        return v;
    }

    /**
     * Gets the Diversity segment of this clonotype
     *
     * @return string identifier of Diversity segment
     */
    public String getD() {
        return d.toString();
    }

    public Segment getDBinary() {
        return d;
    }

    /**
     * Gets the Joining segment of this clonotype
     *
     * @return string identifier of Joining segment
     */
    public String getJ() {
        return j.toString();
    }

    public Segment getJBinary() {
        return j;
    }

    /**
     * Gets the nucleotide sequence of CDR3 region of this clonotype
     *
     * @return
     */
    public String getCdr3nt() {
        return cdr3nt.toString();
    }

    public NucleotideSequence getCdr3ntBinary() {
        return cdr3nt;
    }

    /**
     * Gets the amino acid sequence of CDR3 region of this clonotype
     *
     * @return
     */
    public String getCdr3aa() {
        return cdr3aa.toString();
    }

    public AminoAcidSequence getCdr3aaBinary() {
        return cdr3aa;
    }

    /**
     * Returns {@code true} if this clonotype is in frame, {@code false} otherwise.
     * In case Variable segment sequence is not fully known, this usually involves checking if
     * Variable and Joining segment parts of CDR3 sequence are in-frame
     *
     * @return {@code true} if the clonotype sequence doesn't contain a frameshift, {@code false} otherwise
     */
    public boolean isInFrame() {
        return inFrame;
    }

    /**
     * Returns {@code true} if this clonotype is complete, {@code false} otherwise
     *
     * @return {@code true} if whole CDR3 sequence is known, {@code false} otherwise
     */
    public boolean isComplete() {
        return isComplete;
    }

    /**
     * Returns {@code true} if this clonotype sequence doesn't contain a stop codon, {@code false} otherwise
     *
     * @return {@code true} if the clonotype sequence doesn't contain a stop codon, {@code false} otherwise
     */
    public boolean isNoStop() {
        return noStop;
    }

    /**
     * Returns {@code true} if this clonotype sequence is coding, {@code false} otherwise
     *
     * @return {@code true} if the clonotype sequence doesn't contain a stop codon or a frameshift, {@code false} otherwise
     */
    public boolean isCoding() {
        return noStop && inFrame;
    }

    /**
     * Gets the position of the last nucleotide of Variable segment in CDR3 coordinates
     *
     * @return zero-based position of last Variable segment nucleotide, relative to first CDR3 nucleotide (first bp of conserved Cys codon)
     */
    public int getVEnd() {
        return segmPoints[0];
    }

    /**
     * Gets the position of the first nucleotide of Diversity segment in CDR3 coordinates
     *
     * @return zero-based position of first Diversity segment nucleotide, relative to first CDR3 nucleotide (first bp of conserved Cys codon)
     */
    public int getDStart() {
        return segmPoints[1];
    }

    /**
     * Gets the position of the last nucleotide of Diversity segment in CDR3 coordinates
     *
     * @return zero-based position of last Diversity segment nucleotide, relative to first CDR3 nucleotide (first bp of conserved Cys codon)
     */
    public int getDEnd() {
        return segmPoints[2];
    }

    /**
     * Gets the position of the first nucleotide of Joining segment in CDR3 coordinates
     *
     * @return zero-based position of first Joining segment nucleotide, relative to first CDR3 nucleotide (first bp of conserved Cys codon)
     */
    public int getJStart() {
        return segmPoints[3];
    }

    /**
     * Gets the length of CDR3 region.
     *
     * @return length of CDR3 nucleotide sequence.
     */
    public int getCdr3Length() {
        return cdr3nt.size();
    }

    /**
     * Gets the number of nucleotides added between the Variable segment end and Diversity segment start
     *
     * @return number of nucleotides inserted into Variable-Diversity segment junction, or {@code -1} if Diversity segment is not defined
     */
    public int getVDIns() {
        return (segmPoints[0] >= 0 && segmPoints[1] >= 0) ? segmPoints[1] - segmPoints[0] - 1 : -1;
    }

    /**
     * Gets the number of nucleotides added between the Diversity segment end and Joining segment start
     *
     * @return number of nucleotides inserted into Diversity-Joining segment junction, or {@code -1} if Diversity segment is not defined
     */
    public int getDJIns() {
        return (segmPoints[3] >= 0 && segmPoints[2] >= 0) ? segmPoints[3] - segmPoints[2] - 1 : -1;
    }

    /**
     * Gets the total number of added nucleotides
     *
     * @return number of nucleotides inserted into Variable-Diversity and Diversity-Joining segment junctions, or {@code -1} if Diversity segment is not defined
     */
    public int getInsertSize() {
        return (segmPoints[0] >= 0 && segmPoints[1] >= 0 && segmPoints[2] >= 0 && segmPoints[3] >= 0) ? getVDIns() + getDJIns() : -1;
    }

    /**
     * Gets the number of nucleotides that are between Variable and Joining segments
     *
     * @return number of nucleotides between Variable and Joining segments, including those of Diversity segment (if present)
     */
    public int getNDNSize() {
        return (segmPoints[0] >= 0 && segmPoints[3] >= 0) ? segmPoints[3] - segmPoints[0] - 1 : 0;
    }

    /**
     * INTERNAL used for output
     *
     * @return
     */
    public String getBlank() {
        return CommonUtil.PLACEHOLDER;
    }

    /**
     * Gets the parent sample
     *
     * @return parent clonotype container
     */
    @Override
    public Sample getParent() {
        return parent;
    }

    void append(Clonotype other) {
        this.count += other.count;
        this.freq += other.freq;
    }

    double recalculateFrequency() {
        this.freq = getFreq();
        return freq;
    }

    @Override
    public int compareTo(Clonotype o) {
        return -Integer.compare(this.count, o.count);
    }


    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Clonotype clonotype = (Clonotype) o;

        if (!cdr3nt.equals(clonotype.cdr3nt)) return false;
        if (!j.equals(clonotype.j)) return false;
        if (!parent.equals(clonotype.parent)) return false;
        if (!v.equals(clonotype.v)) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = parent.hashCode();
        result = 31 * result + v.hashCode();
        result = 31 * result + j.hashCode();
        result = 31 * result + cdr3nt.hashCode();
        return result;
    }

    Clonotype withSegments(Segment v, Segment j) {
        return new Clonotype(parent, count, freq,
                segmPoints, v, d, j,
                cdr3nt, cdr3aa,
                inFrame, noStop, isComplete,
                annotation
        );
    }

    @Override
    public Clonotype getClonotype() {
        return this;
    }

    @Override
    public String toString() {
        return "Clonotype{" +
                "parent=" + parent +
                ", count=" + count +
                ", freq=" + freq +
                ", annotation='" + annotation + '\'' +
                ", segmPoints=" + Arrays.toString(segmPoints) +
                ", v=" + v +
                ", d=" + d +
                ", j=" + j +
                ", cdr3nt=" + cdr3nt +
                ", cdr3aa=" + cdr3aa +
                ", inFrame=" + inFrame +
                ", isComplete=" + isComplete +
                ", noStop=" + noStop +
                '}';
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
