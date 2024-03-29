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

import com.antigenomics.vdjtools.ClonotypeWrapperContainer;
import com.antigenomics.vdjtools.io.parser.ClonotypeStreamParser;
import com.antigenomics.vdjtools.misc.Software;
import com.antigenomics.vdjtools.sample.metadata.SampleMetadata;

import java.io.InputStream;
import java.util.*;

/**
 * An implementation of Rep-Seq sample.
 */
public class Sample implements ClonotypeWrapperContainer<Clonotype> {
    private final List<Clonotype> clonotypes = new ArrayList<>();
    private final SampleMetadata sampleMetadata;
    private double frequency = 0;
    private long count = 0;
    private int diversity = 0;
    private String annotationHeader = null;

    public Sample(SampleMetadata sampleMetadata) {
        this.sampleMetadata = sampleMetadata;
    }

    /**
     * Creates a deep copy of a given sample, re-assigning sample metadata.
     *
     * @param other          sample to copy.
     * @param sampleMetadata new sample metadata.
     */
    public Sample(Sample other, SampleMetadata sampleMetadata) {
        this.sampleMetadata = sampleMetadata;
        this.annotationHeader = other.annotationHeader;

        for (Clonotype clonotype : other.clonotypes) {
            this.addClonotype(new Clonotype(clonotype, this));
        }
    }

    /**
     * Creates a new sample by selecting a subset of reads from a given sample. Number of reads
     * to be taken for each clonotypes is specified explicitly.
     *
     * @param other      sample to subset from.
     * @param samplerMap number of reads to take for each clonotype.
     */
    public Sample(Sample other, Map<Clonotype, Integer> samplerMap) {
        this.sampleMetadata = other.sampleMetadata;
        this.annotationHeader = other.annotationHeader;

        for (Clonotype clonotype : other.clonotypes) {
            Integer newCount = samplerMap.get(clonotype);

            if (newCount != null && newCount > 0)
                this.addClonotype(new Clonotype(clonotype, this, newCount));
        }

        Collections.sort(clonotypes);
    }

    public Sample(Sample other, ClonotypeConverter clonotypeConverter) {
        this.sampleMetadata = other.sampleMetadata;
        this.annotationHeader = other.annotationHeader;

        for (Clonotype clonotype : other.clonotypes) {
            this.addClonotype(clonotypeConverter.convert(clonotype));
        }
    }

    /**
     * Creates a new sample by filtering and selecting top N clonotypes from the specified sample.
     *
     * @param other  sample to filter and select from.
     * @param filter a clonotype filter.
     * @param top    if set to value other than -1 will select only top N most abundant matching clonotypes.
     */
    public Sample(Sample other, ClonotypeFilter filter, int top) {
        this.sampleMetadata = other.sampleMetadata;
        this.annotationHeader = other.annotationHeader;

        for (Clonotype clonotype : other.clonotypes) {
            if (top > -1 && this.getDiversity() == top)
                break;

            if (filter.pass(clonotype))
                this.addClonotype(new Clonotype(clonotype, this));
        }
    }

    /**
     * Creates a new sample by selecting top N clonotypes from the specified sample.
     *
     * @param other sample to select from.
     * @param top   if set to value other than -1 will select only top N most abundant matching clonotypes.
     */
    public Sample(Sample other, int top) {
        this(other, BlankClonotypeFilter.INSTANCE, top);
    }

    /**
     * Creates a new sample by filtering clonotypes from the specified sample.
     *
     * @param other  sample to filter.
     * @param filter a clonotype filter.
     */
    public Sample(Sample other, ClonotypeFilter filter) {
        this(other, filter, -1);
    }

    /**
     * Clones a given sample.
     *
     * @param other sample to clone.
     */
    public Sample(Sample other) {
        this(other, BlankClonotypeFilter.INSTANCE, -1);
    }

    /**
     * Reads sample from input stream.
     *
     * @param inputStream    input stream containing plain-text clonotype table.
     * @param sampleMetadata sample metadata.
     * @param software       software, used for parsing.
     * @param top            select top N clonotypes only. Set to -1 to select all clonotypes.
     * @param store          if set to true, will store sample to memory. Otherwise will create an instance of the sample that will be read on demand.
     * @param collapse       if set to true, will collapse the sample combining duplicate clonotypes.
     * @return sample instance.
     */
    public static Sample fromInputStream(InputStream inputStream,
                                         SampleMetadata sampleMetadata,
                                         Software software,
                                         int top, boolean store, boolean collapse) {
        Sample sample = new Sample(sampleMetadata);
        ClonotypeStreamParser clonotypeStreamParser = ClonotypeStreamParser.create(inputStream, software, sample);

        sample.annotationHeader = clonotypeStreamParser.getAnnotationHeader();

        boolean sorted = !collapse;
        int prevCount = Integer.MAX_VALUE;

        Map<Clonotype, Clonotype> existingClonotypes = new HashMap<>();

        for (Clonotype clonotype : clonotypeStreamParser) {
            if (top > -1 && sample.getDiversity() == top)
                break;

            if (clonotype != null) {
                int count = (int) clonotype.getCount();

                if (sorted && count > prevCount) {
                    sorted = false;
                    prevCount = count;
                }

                Clonotype existing = null;

                if (collapse) {
                    existing = existingClonotypes.get(clonotype);

                    if (existing != null) {
                        existing.append(clonotype);
                    } else {
                        existingClonotypes.put(clonotype, clonotype);
                    }
                }

                sample.addClonotype(clonotype, store, existing);
            }
        }

        clonotypeStreamParser.finish(); // report progress

        // on-demand sorting
        if (!sorted)
            Collections.sort(sample.clonotypes);

        // Re-calculate frequencies for per read storing software
        if (software.isPerReadOutput()) {
            sample.frequency = 0;
            for (Clonotype clonotype : sample) {
                sample.frequency += clonotype.recalculateFrequency();
            }
        }

        return sample;
    }

    /**
     * Reads sample from input stream. Stores the sample into memory.
     *
     * @param inputStream    input stream containing plain-text clonotype table.
     * @param sampleMetadata sample metadata.
     * @param software       software, used for parsing.
     * @return sample instance.
     */
    public static Sample fromInputStream(InputStream inputStream,
                                         SampleMetadata sampleMetadata,
                                         Software software) {
        return fromInputStream(inputStream, sampleMetadata, software, -1, true, software.isCollapseRequired());
    }

    private void addClonotype(Clonotype clonotype) {
        addClonotype(clonotype, true, null);
    }

    private void addClonotype(Clonotype clonotype, boolean store, Clonotype existingClonotype) {
        count += clonotype.getCount();
        frequency += clonotype.getFreqAsInInput();

        if (existingClonotype == null) {
            diversity++;
            if (store)
                clonotypes.add(clonotype);
        }
    }

    /**
     * Gets the clonotype annotation header string
     *
     * @return clonotype annotation header or null if not available
     */
    public String getAnnotationHeader() {
        return annotationHeader;
    }


    /**
     * Gets the clonotype annotation header string
     */
    public void setAnnotationHeader(String annotationHeader) {
        this.annotationHeader = annotationHeader;
    }

    /**
     * Gets the metadata associted with the sample.
     *
     * @return sample metadata.
     */
    public SampleMetadata getSampleMetadata() {
        return sampleMetadata;
    }

    /**
     * Gets the list of clonotypes in a given sample.
     * Added for compatibility with 1.8 stream operations.
     *
     * @return clonotype list.
     */
    public List<Clonotype> getClonotypes() {
        return Collections.unmodifiableList(clonotypes);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getFreq() {
        return 1.0;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getFreqAsInInput() {
        return frequency;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public long getCount() {
        return count;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getDiversity() {
        return diversity;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public Clonotype getAt(int index) {
        if (index < 0 || index >= clonotypes.size())
            throw new IndexOutOfBoundsException();
        return clonotypes.get(index);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean isSorted() {
        return true;
    }

    @Override
    public Iterator<Clonotype> iterator() {
        return clonotypes.iterator();
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
