    /*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.pal.alignment;

import net.maizegenetics.pal.datatype.IUPACNucleotides;

import java.nio.ByteBuffer;
import java.util.Arrays;
import net.maizegenetics.pal.datatype.DataType;

/**
 *
 * @author ed
 */
public class Pack1Alignment extends AbstractAlignment {

    protected byte[][] alleleBLOB;
    protected int[] variableSites;  //positions of the variable sites
    protected byte[] SNPid; //SNP IDs
    protected byte[] variableSitesBLOBHeader;  //header is need to clone this object or output
    protected byte[] SNPidBLOBHeader; //header for SNP id BLOB
    protected byte[] majorAlleles,  minorAlleles;
    protected Locus locus;  //this will be changed to locus object and eventually an list of loci
    private boolean isOnlyPolymorphic = true;
    protected int numSites;
    protected int SNPidLength; //length of each SNPid in the BLOB
    protected byte[] qualitySiteScore = null;

    //One thought is to have a inflated cache of 1M elements.  We can inflate 300M elements in 2s, so this would not
    //be very expensive and provide high performance for analyses in smaller ranges.
    public Pack1Alignment(byte[][] alleleBLOB, byte[] variableSitesBLOB, byte[] SNPidBLOB, byte[] qualityBLOB) {
        super(AllelePositionBLOBUtils.getTaxa(alleleBLOB), new IUPACNucleotides());
        this.alleleBLOB = alleleBLOB;
        String locusName=(new String(variableSitesBLOB, GdpdmBLOBUtils.locusField[0], GdpdmBLOBUtils.locusField[1])).trim();
        String chromoName=(locusName.length()<3)?locusName:"?";
        locus = new Locus(locusName, chromoName, -1, -1, null, null);
        initVariableSites(variableSitesBLOB);
        initSNPids(SNPidBLOB);
        if (qualityBLOB != null) {
            qualitySiteScore = qualityBLOB;
        }
        initMajorMinorAlleles();
        validateInstance();
    }

    public Pack1Alignment(byte[][] alleleBLOB, byte[] variableSitesBLOB, byte[] SNPidBLOB) {
        this(alleleBLOB, variableSitesBLOB, SNPidBLOB, null);
    }

    public Pack1Alignment(byte[][] alleleBLOB) {
        super(AllelePositionBLOBUtils.getTaxa(alleleBLOB), new IUPACNucleotides());
        this.alleleBLOB = alleleBLOB;
        locus = new Locus((new String(alleleBLOB[0], GdpdmBLOBUtils.locusField[0],
                GdpdmBLOBUtils.locusField[1])).trim(), null, -1, -1, null, null);

        ByteBuffer bb = ByteBuffer.wrap(alleleBLOB[0]);
        bb.position(GdpdmBLOBUtils.numSitesField[0]);
        numSites = bb.getInt();
        System.out.println("num: " + numSites);
        initMajorMinorAlleles();
        validateInstance();
    }

    private void validateInstance() {
        if (getPositionInLocus(0) > getPositionInLocus(getSiteCount() - 1)) {
            throw new IllegalStateException("Pack1Alignment: validateInstance: Start position is larger than end position.");
        }
    }

    private void initVariableSites(byte[] variableSitesBLOB) {
        ByteBuffer bb = ByteBuffer.wrap(variableSitesBLOB);
        bb.position(GdpdmBLOBUtils.numSitesField[0]);
        this.numSites = bb.getInt();
        variableSites = new int[getSiteCount()];
        bb.position(GdpdmBLOBUtils.totalHeaderPadding);
        for (int i = 0; i < getSiteCount(); i++) {
            variableSites[i] = bb.getInt();
        }
        variableSitesBLOBHeader = new byte[GdpdmBLOBUtils.totalHeaderPadding];
        bb.position(0);
        bb.get(variableSitesBLOBHeader, 0, GdpdmBLOBUtils.totalHeaderPadding);
    }

    private void initSNPids(byte[] SNPidBLOB) {
        ByteBuffer bb = ByteBuffer.wrap(SNPidBLOB);
        SNPidBLOBHeader = new byte[GdpdmBLOBUtils.totalHeaderPadding];
        bb.get(SNPidBLOBHeader, 0, GdpdmBLOBUtils.totalHeaderPadding);
        bb.position(GdpdmBLOBUtils.idLengthField);
        this.SNPidLength = bb.getInt()/8;
        SNPid = new byte[SNPidBLOB.length - GdpdmBLOBUtils.totalHeaderPadding];
        bb.position(GdpdmBLOBUtils.totalHeaderPadding);
        bb.get(SNPid);
    }

    private void initMajorMinorAlleles() {
        majorAlleles = new byte[getSiteCount()];
        minorAlleles = new byte[getSiteCount()];
        for (int i = 0; i < getSiteCount(); i++) {
            byte[] alleles = getAlleles(i);
            majorAlleles[i] = (alleles.length > 0) ? alleles[0] : DataType.UNKNOWN_BYTE;
            minorAlleles[i] = (alleles.length > 1) ? alleles[1] : DataType.UNKNOWN_BYTE;
        }
    }

    public byte[] getVariableSitesBLOB() {
        byte[] vsBLOB = new byte[GdpdmBLOBUtils.totalHeaderPadding + (4 * variableSites.length)];
        ByteBuffer bb = ByteBuffer.wrap(vsBLOB);
        bb.put(variableSitesBLOBHeader);
        for (int v : variableSites) {
            bb.putInt(v);
        }
        return vsBLOB;
    }

    public byte[] getSNPidBLOB() {
        byte[] SNPidBLOB = new byte[SNPidBLOBHeader.length + SNPid.length];
        ByteBuffer bb = ByteBuffer.wrap(SNPidBLOB);
        bb.put(SNPidBLOBHeader);
        bb.put(SNPid);
        return SNPidBLOB;
    }

    public byte[][] getAlleleBLOBs() {
        return alleleBLOB;
    }

    public byte[] getAlleleBLOBs(int seq) {
        return alleleBLOB[seq];
    }

    public char getBaseChar(int seq, int site) {
        return (char) getBase(seq, site);
    }

    public byte getBase(int seq, int site) {
        if (isOnlyPolymorphic) {
            return AllelePositionBLOBUtils.getBaseFromAlleleBLOB(alleleBLOB[seq], site);
        }
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public byte getBaseCode(int seq, int site) {
        return AlignmentUtils.getAlleleValue(alleleBLOB[seq], site, Pack2Alignment.NUM_BITS_PER_VALUE.four, GdpdmBLOBUtils.totalHeaderPadding);
    }

    @Override
    public byte getMajorAllele(int site) {
        return majorAlleles[site];
    }

    @Override
    public byte getMinorAllele(int site) {
        return minorAlleles[site];
    }

    /**
     * Returns value for given taxon , site, and allele.
     *
     * @param taxon taxon
     * @param site site
     * @param allele allele
     *
     * @return value
     */
    public byte getBase(int taxon, int site, int allele) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

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
    public byte getBase(int taxon, Locus locus, int physicalPosition, int allele) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public byte[] getBaseRange(int seq, int startSite, int endSite) {
        if (isOnlyPolymorphic) {
            return AllelePositionBLOBUtils.getBaseRangeFromAlleleBLOB(alleleBLOB[seq], startSite, endSite);
        }
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public int getSiteCount() {
        return numSites;
    }

    public int getSNPidLength() {
        return SNPidLength;
    }

    public byte[] getSNPid() {
        return SNPid;
    }

    public String[] getSNPIDs() {
        int sites = getSiteCount();
        String[] SNPids = new String[sites];
        for (int i = 0; i < sites; i++) {
            SNPids[i] = getSNPID(i);
        }
        return SNPids;
    }

    public String getSNPID(int site) {
        ByteBuffer bb = ByteBuffer.wrap(getSNPid());
        bb.position(site * getSNPidLength());
        StringBuilder sb = new StringBuilder(getSNPidLength());
        for (int i = 0; i < getSNPidLength(); i++) {
            sb.append((char) bb.get());
        }
        String snpID = sb.toString().trim();
        if (snpID.length() == 0) {
            return getLocusName(site) + "_" + getPositionInLocus(site);
        }
        return snpID;
    }

    public int getLocusSiteCount(Locus locus) {
        return getSiteCount();
    }

    public int getPositionInLocus(int site) {
        try {
            return variableSites[site];
        } catch (Exception e) {
            return site;
        }
    }

    public int getSiteOfPhysicalPosition(int physicalPosition, Locus locus) {
        try {
            return Arrays.binarySearch(variableSites, physicalPosition);
        } catch (Exception e) {
            return physicalPosition;
        }
    }

    public byte getPositionType(int site) {
        return PositionType.ALL_GROUP;
    //todo need to check for indels and return that if necessary.
    }

    public byte[] getPositionTypes() {
        //throw new UnsupportedOperationException("Not supported yet.");
        byte[] positionTypes = new byte[this.getSiteCount()];
        for (int i = 0; i < positionTypes.length; i++) {
            positionTypes[i] = this.getPositionType(i);
        }
        return positionTypes;
    }

    public Locus getLocus(int site) {
        return locus;
    }

    public Locus[] getLoci() {
        return new Locus[]{locus};
    }

    public int getNumLoci() {
        return 1;
    }

    public float getSiteScore(int seq, int site) {
        if (qualitySiteScore == null) {
            return Float.NaN;
        }
        return (float) qualitySiteScore[site + GdpdmBLOBUtils.totalHeaderPadding];
    }

    public boolean hasSiteScores() {
        return (qualitySiteScore != null);
    }

    public SITE_SCORE_TYPE getSiteScoreType() {
        if (qualitySiteScore == null) {
            return Alignment.SITE_SCORE_TYPE.None;
        }
        return Alignment.SITE_SCORE_TYPE.QualityScore;
    }

    /**
     * Return sorted list of alleles from highest frequency to lowest
     * at given site in alignment. Resulting double dimension array
     * holds alleles (actually bytes) in result[0].  And the counts
     * are in result[1]. Counts haploid values twice and diploid
     * values once. Higher ploids are considered unknown and not counted.
     *
     * @param a alignment
     * @param site site
     * @param includeGAP whether to include GAP
     * @return sorted list of alleles and counts
     */
    public int[][] getAllelesSortedByFrequency(int site, boolean includeGAP) {

        int[] alleleCounts = new int[256];
        for (int j = 0, n = getSequenceCount(); j < n; j++) {
            byte[] currentByteArray = AllelePositionBLOBUtils.getSNPValueFromHalfByte(getBaseCode(j, site));

            // Counts haploid values twice and diploid values once.
            // Higher ploids are considered unknown and not counted.
            int incrementAmount = 0;
            if (currentByteArray.length == 1) {
                incrementAmount = 2;
            } else if (currentByteArray.length == 2) {
                incrementAmount = 1;
            }

            // cycles through for each possible allele.
            for (int i = 0; i < currentByteArray.length; i++) {
                byte current = currentByteArray[i];
                alleleCounts[current] = alleleCounts[current] + incrementAmount;
            }
        }

        int numAlleles = 0;
        for (int j = 0, n = alleleCounts.length; j < n; j++) {
            if ((alleleCounts[j] != 0) && (GdpdmBLOBUtils.UNKNOWN_CHARACTER != j) && (GdpdmBLOBUtils.GAP_CHARACTER != j)) {
                numAlleles++;
            } else if (includeGAP && (alleleCounts[j] != 0) && (GdpdmBLOBUtils.GAP_CHARACTER == j)) {
                numAlleles++;
            }
        }

        int[][] result = new int[2][numAlleles];
        int count = 0;
        for (int k = 0, n = alleleCounts.length; k < n; k++) {
            if ((alleleCounts[k] != 0) && (GdpdmBLOBUtils.UNKNOWN_CHARACTER != k) && (GdpdmBLOBUtils.GAP_CHARACTER != k)) {
                result[1][count] = alleleCounts[k];
                result[0][count] = k;
                count++;
            } else if ((GdpdmBLOBUtils.GAP_CHARACTER == k) && (includeGAP) && (alleleCounts[k] != 0)) {
                result[1][count] = alleleCounts[k];
                result[0][count] = k;
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
