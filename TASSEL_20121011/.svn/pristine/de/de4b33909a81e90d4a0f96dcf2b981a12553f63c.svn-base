/*
 * SiteSummary
 */
package net.maizegenetics.pal.alignment;

/**
 *
 * @author terry
 */
public class SiteSummary {

    private final String mySiteName;
    private final int myPosition;
    private final Locus myLocus;
    private final char[] myAlleles;
    private final int[] myAlleleCounts;
    private final int myTotalSeqCount;

    public SiteSummary(String name, int position, Locus locus, int numSeqences, char[] alleles, int[] alleleCounts) {

        mySiteName = name;
        myPosition = position;
        myLocus = locus;
        myTotalSeqCount = numSeqences;

        if (alleles.length != alleleCounts.length) {
            throw new IllegalArgumentException("SiteSummary: init: number of alleles should equal number allele counts.");
        }
        myAlleles = alleles;
        myAlleleCounts = alleleCounts;

    }

    public static SiteSummary getInstance(Alignment a, int site) {
        return getInstance(a, site, true);
    }

    public static SiteSummary getInstance(Alignment a, int site, boolean includeGap) {

        String name = a.getSNPID(site);
        int position = a.getPositionInLocus(site);
        Locus locus = a.getLocus(site);
        int numSequences = a.getSequenceCount();
        int[][] alleles = a.getAllelesSortedByFrequency(site, includeGap);
        int numAlleles = alleles[0].length;
        int[] alleleCounts = new int[numAlleles];
        char[] alleleValues = new char[numAlleles];
        for (int i = 0, n = numAlleles; i < n; i++) {
            alleleCounts[i] = alleles[1][i];
            alleleValues[i] = (char) alleles[0][i];
        }

        return new SiteSummary(name, position, locus, numSequences, alleleValues, alleleCounts);

    }

    public double getAlleleFrequency(int allele) {
        return (double) myAlleleCounts[allele] / (double) myTotalSeqCount;
    }

    /**
     * This returns the number of sites with missing data.
     * GAP is considered a state by default. (See getInstance methods)
     *
     * @return number of missing alleles.
     */
    public int getNumberMissing() {
        int result = myTotalSeqCount * 2;
        for (int i = 0; i < myAlleles.length; i++) {
            result = result - myAlleleCounts[i];
        }
        return result / 2;
    }

    public String getName() {
        return mySiteName;
    }

    public int getPosition() {
        return myPosition;
    }

    public Locus getLocus() {
        return myLocus;
    }

    public char[] getAlleles() {
        return myAlleles;
    }

    public int[] getAlleleCounts() {
        return myAlleleCounts;
    }
}
