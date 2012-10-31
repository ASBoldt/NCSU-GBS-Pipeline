/*
 */
package net.maizegenetics.pal.alignment;

import java.io.PrintWriter;
import java.io.StringWriter;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import net.maizegenetics.pal.datatype.DataType;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.SimpleIdGroup;

/**
 *
 */
abstract public class AbstractAlignment implements Alignment {

    private static final long serialVersionUID = -5197800047652332969L;
    private final DataType myDataType;
    private final IdGroup myIdGroup;
    /** The Reference */
    private byte[] myReference = null;

    public AbstractAlignment(Alignment a) {
        this(a.getIdGroup(), a.getDataType());
    }

    public AbstractAlignment(IdGroup idGroup, DataType dataType) {
        myIdGroup = SimpleIdGroup.getInstance(idGroup);
        myDataType = dataType;
    }

    public char[] getBaseCharRange(int taxon, int startSite, int endSite) {

        char[] result = new char[endSite - startSite + 1];
        for (int i = startSite; i < endSite + 1; i++) {
            result[i] = getBaseChar(taxon, i);
        }
        return result;

    }

    public byte[] getBaseRange(int taxon, int startSite, int endSite) {

        byte[] result = new byte[endSite - startSite + 1];
        for (int i = startSite; i < endSite + 1; i++) {
            result[i] = getBase(taxon, i);
        }
        return result;

    }

    public String getBaseString(int taxon, int site) {
        return myDataType.getFormattedString(getBaseChar(taxon, site));
    }

    public byte[] getReference(int startSite, int endSite) {

        if ((myReference == null) || (myReference.length == 0)) {
            return null;
        }

        byte[] result = new byte[endSite - startSite];
        for (int i = startSite; i < endSite; i++) {
            result[i] = getReferenceAllele(i);
        }
        return result;

    }

    public int getSequenceCount() {
        return myIdGroup.getIdCount();
    }

    public DataType getDataType() {
        return myDataType;
    }

    public String getAlignedSequenceString(int sequence) {

        char[] data = new char[getSiteCount()];
        for (int i = 0; i < getSiteCount(); i++) {
            data[i] = getBaseChar(sequence, i);
        }
        return new String(data);

    }

    public String getAlignedSequenceString(int sequence, int fromSite, int toSite) {

        char[] data = new char[toSite - fromSite];
        for (int i = fromSite; i < toSite; i++) {
            data[i] = getBaseChar(sequence, i);
        }
        return new String(data);

    }

    public String getFormattedSequenceString(int sequence) {
        return getAlignedSequenceString(sequence);
    }

    public String getFormattedSequenceString(int sequence, int fromSite, int toSite) {
        return getAlignedSequenceString(sequence, fromSite, toSite);
    }

    public int[][] getStates() {

        int[][] indices = new int[getSequenceCount()][getSiteCount()];

        for (int i = 0; i < getSequenceCount(); i++) {

            for (int j = 0; j < getSiteCount(); j++) {

                indices[i][j] = myDataType.getState(getBaseChar(i, j));

                if (indices[i][j] >= myDataType.getNumStates()) {
                    indices[i][j] = -1;
                }
            }
        }

        return indices;

    }

    public String getLocusName(int site) {
        return getLocus(site).getName();
    }

    public void report(PrintWriter out) {
    }

    public float getSiteScore(int seq, int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public float[][] getSiteScores() {
        if (hasSiteScores() == false) {
            return null;
        }
        float[][] f = new float[getSequenceCount()][getSiteCount()];
        for (int i = 0; i < getSequenceCount(); i++) {
            for (int j = 0; j < getSiteCount(); j++) {
                f[i][j] = getSiteScore(i, j);
            }
        }
        return f;
    }

    public boolean hasSiteScores() {
        return false;
    }

    public SITE_SCORE_TYPE getSiteScoreType() {
        return Alignment.SITE_SCORE_TYPE.None;
    }

    public int getIndelSize(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public boolean isIndel(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public boolean isAllPolymorphic() {

        for (int i = 0, n = getSiteCount(); i < n; i++) {
            if (!isPolymorphic(i)) {
                return false;
            }
        }

        return true;

    }

    public boolean isPolymorphic(int site) {

        char first = DataType.UNKNOWN_CHARACTER;
        for (int i = 0, n = getSequenceCount(); i < n; i++) {
            char current = getBaseChar(i, site);
            if ((current != DataType.UNKNOWN_CHARACTER) && (current != DataType.GAP_CHARACTER)) {
                if (first == DataType.UNKNOWN_CHARACTER) {
                    first = current;
                } else if (first != current) {
                    return true;
                }
            }
        }

        return false;

    }

    public byte[] getMinorAlleles(int site) {
        int[][] alleles = getAllelesSortedByFrequency(site, true);
        int resultSize = alleles[0].length - 1;
        byte[] result = new byte[resultSize];
        for (int i = 0; i < resultSize; i++) {
            result[i] = (byte) alleles[0][i + 1];
        }
        return result;
    }

    public byte[] getAlleles(int site) {
        int[][] alleles = getAllelesSortedByFrequency(site, true);
        int resultSize = alleles[0].length;
        byte[] result = new byte[resultSize];
        for (int i = 0; i < resultSize; i++) {
            result[i] = (byte) alleles[0][i];
        }
        return result;
    }

    public double getMinorAlleleFrequency(int site) {

        int[][] alleles = getAllelesSortedByFrequency(site, true);

        int numAlleles = alleles[0].length;
        if (numAlleles >= 2) {
            int totalNonMissing = 0;
            for (int i = 0; i < numAlleles; i++) {
                totalNonMissing = totalNonMissing + alleles[1][i];
            }
            return (double) alleles[1][1] / totalNonMissing;
        } else {
            return 0.0;
        }

    }

    public byte getMajorAllele(int site) {

        int[][] alleles = getAllelesSortedByFrequency(site, true);

        if (alleles[0].length >= 1) {
            return (byte) alleles[0][0];
        } else {
            return DataType.UNKNOWN_BYTE;
        }

    }

    public byte getMinorAllele(int site) {

        int[][] alleles = getAllelesSortedByFrequency(site, true);

        if (alleles[0].length >= 2) {
            return (byte) alleles[0][1];
        } else {
            return DataType.UNKNOWN_BYTE;
        }

    }

    public SiteSummary getSiteSummary(int site) {
        return SiteSummary.getInstance(this, site);
    }

    public String toString() {

        StringWriter sw = new StringWriter();
        AlignmentUtils.print(this, new PrintWriter(sw));

        return sw.toString();
    }

    public IdGroup getIdGroup() {
        return myIdGroup;
    }

    public String getTaxaName(int index) {
        return myIdGroup.getIdentifier(index).getName();
    }

    public String getFullTaxaName(int index) {
        return myIdGroup.getIdentifier(index).getFullName();
    }

    /**
     * Sets reference sequence.  Number of sites should be
     * set before calling this.
     * 
     * @param reference
     */
    protected void setReference(byte[] reference) {
        if ((reference == null) || (reference.length == 0)) {
            myReference = null;
        } else {
            myReference = reference;
            AlignmentUtils.cleanSequenceString(myReference);
            if (getSiteCount() != myReference.length) {
                throw new IllegalArgumentException("AbstractAlignment: init: reference differs in length with sequences.");
            }
        }
    }

    public boolean hasReference() {
        if ((myReference != null) && (myReference.length != 0)) {
            return true;
        }
        return false;
    }

    public byte getReferenceAllele(int site) {
        return myReference[site];
    }

    /**
     * Return reference sequence.
     *
     * @return reference sequence
     */
    public byte[] getReference() {
        return myReference;
    }

    /**
     * Gets the Genome Assembly.
     *
     * @return the genome assembly.
     */
    public String getGenomeAssembly() {
        throw new UnsupportedOperationException("Not supported.");
    }

    /**
     * Return whether is positive strand at given site.
     *
     * @param site site
     *
     * @return whether is positive strand.
     */
    public boolean isPositiveStrand(int site) {
        throw new UnsupportedOperationException("Not supported.");
    }

    public Alignment[] getAlignments() {
        return new Alignment[]{this};
    }

    /**
     * Return sorted list of alleles from highest frequency to lowest
     * at given site in alignment. Resulting double dimension array
     * holds alleles (bytes or chars in some cases) in result[0].  And the counts
     * are in result[1]. This counts every value twice
     * as this default implementation does support diploids.
     *
     * @param a alignment
     * @param site site
     * @param includeGAP whether to include GAP
     * @return sorted list of alleles and counts
     */
    public int[][] getAllelesSortedByFrequency(int site, boolean includeGAP) {

        Map alleleCounts = new HashMap();
        for (int j = 0, n = getSequenceCount(); j < n; j++) {
            char current = getBaseChar(j, site);
            Integer count = (Integer) alleleCounts.get(current);
            if (count != null) {
                alleleCounts.put(current, new Integer(count.intValue() + 1));
            } else {
                alleleCounts.put(current, 1);
            }
        }

        int numAlleles = 0;
        Iterator itr = alleleCounts.keySet().iterator();
        while (itr.hasNext()) {
            char current = ((Character) itr.next()).charValue();
            if ((GdpdmBLOBUtils.UNKNOWN_CHARACTER != current) && (GdpdmBLOBUtils.GAP_CHARACTER != current)) {
                numAlleles++;
            } else if (includeGAP && (GdpdmBLOBUtils.GAP_CHARACTER == current)) {
                numAlleles++;
            }
        }

        int[][] result = new int[2][numAlleles];
        int count = 0;
        itr = alleleCounts.keySet().iterator();
        while (itr.hasNext()) {
            char current = ((Character) itr.next()).charValue();
            if ((GdpdmBLOBUtils.UNKNOWN_CHARACTER != current) && (GdpdmBLOBUtils.GAP_CHARACTER != current)) {
                result[1][count] = ((Integer) alleleCounts.get(current)).intValue() * 2;
                result[0][count] = current;
                count++;
            } else if (includeGAP && (GdpdmBLOBUtils.GAP_CHARACTER == current)) {
                result[1][count] = ((Integer) alleleCounts.get(current)).intValue() * 2;
                result[0][count] = current;
                count++;
            }
        }

        boolean change = true;
        while (change) {

            change = false;

            for (int k = 0; k < numAlleles - 1; k++) {

                if (result[1][k] < result[1][k + 1]) {

                    int temp = result[0][k];
                    result[0][k] = result[0][k + 1];
                    result[0][k + 1] = temp;

                    int tempCount = result[1][k];
                    result[1][k] = result[1][k + 1];
                    result[1][k + 1] = tempCount;

                    change = true;
                }
            }

        }

        return result;

    }
}