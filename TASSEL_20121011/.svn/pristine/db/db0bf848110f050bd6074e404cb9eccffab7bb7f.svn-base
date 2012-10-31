/*
 * Alignment
 */
package net.maizegenetics.pal.alignment;

import java.io.PrintWriter;
import java.io.Serializable;

import net.maizegenetics.pal.datatype.DataType;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.report.Report;

/**
 *
 */
public interface Alignment extends Serializable, Report {

    public static enum SITE_SCORE_TYPE {

        None, QualityScore, ImputedProbablity, Dosage
    };

    /**
     * Returns value for given taxon and site.
     *
     * @param taxon taxon
     * @param site site
     *
     * @return value
     */
    public char getBaseChar(int taxon, int site);

    /**
     * Returns value for given taxon and site.
     *
     * @param taxon taxon
     * @param site site
     *
     * @return value
     */
    public String getBaseString(int taxon, int site);

    /**
     * Return sequence for given taxon in specified range
     * (end site included).
     *
     * @param taxon taxon
     * @param startSite start site
     * @param endSite end site
     *
     * @return sequence
     */
    public char[] getBaseCharRange(int taxon, int startSite, int endSite);

    /**
     * Returns value for given taxon and site.
     *
     * @param taxon taxon
     * @param site site
     *
     * @return value
     */
    public byte getBase(int taxon, int site);

    /**
     * Returns value for given taxon , site, and allele.
     *
     * @param taxon taxon
     * @param site site
     * @param allele allele
     *
     * @return value
     */
    public byte getBase(int taxon, int site, int allele);

    /**
     * Returns value for given taxon, locus, physical
     * position, and allele.
     *
     * @param taxon taxon
     * @param locus locus
     * @param physicalPosition physical position
     * @param allele allele
     * 
     * @return value
     */
    public byte getBase(int taxon, Locus locus, int physicalPosition, int allele);

    /**
     * Return sequence for given taxon in specified range
     * (end site included).
     *
     * @param taxon taxon
     * @param startSite start site
     * @param endSite end site
     *
     * @return sequence
     */
    public byte[] getBaseRange(int taxon, int startSite, int endSite);

    /**
     * Return reference sequence for given reference in specified range
     * (end site not included).
     *
     * @param startSite start site
     * @param endSite end site
     *
     * @return reference sequence
     */
    public byte[] getReference(int startSite, int endSite);

    /**
     * Return reference sequence.
     *
     * @return reference sequence
     */
    public byte[] getReference();

    /**
     * Return whether this alignment has defined reference sequence.
     *
     * @return true if this alignment has reference sequence.
     */
    public boolean hasReference();

    /**
     * Get site summary for specified site.
     *
     * @param site site
     * @return site summary
     */
    public SiteSummary getSiteSummary(int site);

    /**
     *  Get SNP IDs.
     *
     * @return site names.
     */
    public String[] getSNPIDs();

    /**
     * Get SNP ID for specified site.
     *
     * @param site site
     * @return site name
     */
    public String getSNPID(int site);

    /**
     * Returns total number of sites for all alignments.
     *
     * @return number of sites
     */
    public int getSiteCount();

    /**
     * Return number of sites for given locus.
     *
     * @param locus locus
     *
     * @return number of sites
     */
    public int getLocusSiteCount(Locus locus);

    /**
     * Returns number of sequences (taxa).
     *
     * @return number of sequences
     */
    public int getSequenceCount();

    /**
     * Return DataType of this alignment.
     *
     * @return data type
     */
    public DataType getDataType();

    /**
     * Returns string representation of single sequence in
     * alignment with gap characters included.
     *
     * @return aligned sequence string
     */
    public String getAlignedSequenceString(int sequence);

    /**
     * Returns string representation of single sequence in
     * alignment with gap characters included.
     *
     * @param sequence sequence
     * @param fromSite start site
     * @param toSite end site (exclusive)
     *
     * @return formatted string
     */
    public String getAlignedSequenceString(int sequence, int fromSite, int toSite);

    /**
     * Returns a string representing a single sequence (including gaps)
     * from this alignment for easy visual reading.  This is most important
     * for TextDataType formats.
     *
     * @param sequence sequence
     *
     * @return formatted string
     */
    public String getFormattedSequenceString(int sequence);

    /**
     * Returns a string representing a single sequence (including gaps)
     * from this alignment for easy visual reading.  This is most important
     * for TextDataType formats.
     *
     * @param sequence sequence
     * @param fromSite start site
     * @param toSite end site (exclusive)
     *
     * @return formatted string
     */
    public String getFormattedSequenceString(int sequence, int fromSite, int toSite);

    /**
     * Fills a [numsequences][length] matrix with indices.
     * Each index represents the sequence state (values defined by the
     * dataType), -1 means a gap.
     *
     * @return states
     */
    public int[][] getStates();

    /**
     * Returns the physical position at given site.
     *
     * @param site site
     *
     * @return physical position
     */
    public int getPositionInLocus(int site);

    /**
     * Return site of given physical position in locus.
     *
     * @param physicalPosition physical position
     * @param locus locus
     *
     * @return index
     */
    public int getSiteOfPhysicalPosition(int physicalPosition, Locus locus);

    /**
     * Returns position type for given site.
     * (eg.  I=intron, E=exon, P=promoter, 1=first, 2=second, 3=third, etc.)
     *
     * @param site site
     *
     * @return position type
     */
    public byte getPositionType(int site);

    /**
     * Returns position types.
     * (eg.  I=intron, E=exon, P=promoter, 1=first, 2=second, 3=third, etc.)
     *
     * @return position types
     */
    public byte[] getPositionTypes();

    /**
     * Return Locus Name for given site.
     *
     * @param site site
     *
     * @return Locus Name
     */
    public String getLocusName(int site);

    /**
     * Return Locus for given site.
     *
     * @param site site
     *
     * @return Locus
     */
    public Locus getLocus(int site);

    /**
     * Return all loci.
     *
     * @return loci
     */
    public Locus[] getLoci();

    /**
     * Return number of loci.
     *
     * @return number of loci
     */
    public int getNumLoci();

    /**
     * Writes a report of this alignment to give writer.
     *
     * @param out writer
     */
    public void report(PrintWriter out);

    /**
     * Returns the site score of the given sequence and site.
     *
     * @param seq seqence index
     * @param site site
     *
     * @return site score.
     */
    public float getSiteScore(int seq, int site);

    /**
     * Returns the site scores.
     *
     * @return site scores.
     */
    public float[][] getSiteScores();

    /**
     * Returns true if this alignment has site scores.
     *
     * @return true if this alignment has site scores.
     */
    public boolean hasSiteScores();

    /**
     * Return what type of site scores this alignment has.
     *
     * @return site score type.
     */
    public SITE_SCORE_TYPE getSiteScoreType();

    /**
     * Return size of indel at given site.
     *
     * @param site site
     *
     * @return indel size
     */
    public int getIndelSize(int site);

    /**
     * Returns whether give site is an indel.
     *
     * @param site site
     *
     * @return true if indel
     */
    public boolean isIndel(int site);

    /**
     * Returns whether all sites are polymorphic.
     *
     * @return true if all sites are polymorphic.
     */
    public boolean isAllPolymorphic();

    /**
     * Return whether given site is polymorphic.
     *
     * @param site site
     *
     * @return true if given site is polymorphic.
     */
    public boolean isPolymorphic(int site);

    /**
     * Return reference allele at given site.
     *
     * @param site site
     *
     * @return reference allele
     */
    public byte getReferenceAllele(int site);

    /**
     * Return most common allele at given site.
     * Gap is included as state.
     *
     * @param site site
     *
     * @return most common allele
     */
    public byte getMajorAllele(int site);

    /**
     * Return most common minor allele at given site.
     * Gap is included as state.
     *
     * @param site site
     *
     * @return most common minor allele
     */
    public byte getMinorAllele(int site);

    /**
     * Return all minor alleles at given site.
     * Gap is included as state.
     *
     * @param site site
     *
     * @return all minor alleles
     */
    public byte[] getMinorAlleles(int site);

    /**
     * Return all alleles at given site in order of frequency.
     * Gap is included as state.
     *
     * @param site site
     *
     * @return all alleles
     */
    public byte[] getAlleles(int site);

    /**
     * Return frequency for most common minor allele at given site.
     * Gap is included as state.
     *
     * @param site site
     *
     * @return frequency
     */
    public double getMinorAlleleFrequency(int site);

    /**
     * Return this alignment's id group.
     *
     * @return id group.
     */
    public IdGroup getIdGroup();

    /**
     * Return taxa name at given index.
     *
     * @param index
     * @return taxa name
     */
    public String getTaxaName(int index);
    
    /**
     * Return full taxa name at given index.
     * 
     * @param index
     * @return full taxa name
     */
    public String getFullTaxaName(int index);

    /**
     * Gets the Genome Assembly.
     *
     * @return the genome assembly.
     */
    public String getGenomeAssembly();

    /**
     * Return whether is positive strand at given site.
     *
     * @param site site
     *
     * @return whether is positive strand.
     */
    public boolean isPositiveStrand(int site);

    /**
     * Returns individual alignments within this alignment.
     *
     * @return list of alignments.
     */
    public Alignment[] getAlignments();

    /**
     * Return sorted list of alleles from highest frequency to lowest
     * at given site in alignment. Resulting double dimension array
     * holds alleles (bytes or chars in some cases) in result[0].  And the counts
     * are in result[1]. Counts haploid values twice and diploid
     * values once. Higher ploids are considered unknown and not counted.
     *
     * @param a alignment
     * @param site site
     * @param includeGAP whether to include GAP
     * @return sorted list of alleles and counts
     */
    public int[][] getAllelesSortedByFrequency(int site, boolean includeGAP);
}
