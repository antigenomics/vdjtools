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

package com.antigenomics.vdjtools.overlap

/**
 * An enum that defines the clonotype matching rule. 
 * Used both when collapsing "identical" clonotypes in a sample and when matching clonotypes between samples 
 */
public enum OverlapType {

    /**
     * Intersection rule. Clonotypes match if their CDR3 nucleotide sequences match 
     */
    /*   */ Nucleotide("nt", false),
    /**
     * Intersection rule. Clonotypes match if both V segments and CDR3 nucleotide sequences match 
     */
            NucleotideV("ntV", false),
    /**
     * Intersection rule. Clonotypes match if both V, J segments and CDR3 nucleotide sequences match 
     */
            NucleotideVJ("ntVJ", false),
    /**
     * Intersection rule. Clonotypes match if their CDR3 amino acid sequences match  
     */
            AminoAcid("aa", true),
    /**
     * Intersection rule. Clonotypes match if both V segments and CDR3 amino acid sequences match 
     */
            AminoAcidV("aaV", true),
    /**
     * Intersection rule. Clonotypes match if both V, J segments and CDR3 amino acid sequences match 
     */
            AminoAcidVJ("aaVJ", true),
    /**
     * Intersection rule. Clonotypes match if their CDR3 amino acid sequences match, 
     * but their CDR3 nucleotide sequences do not match. Contamination-proof matching
     */
            AminoAcidNonNucleotide("aa!nt", true),
    /**
     * Intersection rule. Clonotypes match if both V, J segments, somatic hypermutations and CDR3 nucleotide sequences match
     */
            Strict("strict", false)

    final String shortName
    final boolean aminoAcid

    /**
     * Defines a new clonotype matching rule
     * @param shortName short name
     * @param aminoAcid if {@code true} matching relies on amino acid sequences; {@code false} for nucleotide sequences
     */
    public OverlapType(String shortName, boolean aminoAcid) {
        this.shortName = shortName
        this.aminoAcid = aminoAcid
    }

    /**
     * Gets {@code IntersectionType} by short name
     * @param shortName short name
     * @return
     */
    public static OverlapType getByShortName(String shortName) {
        values().find { it.shortName.toUpperCase() == shortName.toUpperCase() }
    }

    /**
     * A list of existing {@code IntersectionType} short names
     */
    public static String allowedNames = values().collect { it.shortName }.join(",")
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
