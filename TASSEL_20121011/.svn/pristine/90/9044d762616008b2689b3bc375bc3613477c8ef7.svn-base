/*
 * Pack2Alignment
 */
package net.maizegenetics.pal.alignment;

import net.maizegenetics.pal.datatype.TextDataType;
import net.maizegenetics.pal.datatype.DataType;
import net.maizegenetics.pal.ids.IdGroup;

import java.nio.ByteBuffer;

import java.util.Arrays;

/**
 *
 * @author terry
 */
public class Pack2Alignment extends AbstractAlignment {

    public static enum NUM_BITS_PER_VALUE {

        one, two, four, eight
    };
    private byte[][] myData;
    private String[] myStates;
    private NUM_BITS_PER_VALUE myNumBits;
    private int[] myVariableSites;  //positions of the variable sites
    private byte[] mySNPid; //SNP IDs
    private byte[] myVariableSitesBLOBHeader;  //header is need to clone this object or output
    private byte[] mySNPidBLOBHeader; //header for SNP id BLOB
    private byte[] myMajorAlleles,  myMinorAlleles;
    private Locus myLocus;
    private int myNumSites;
    private int mySNPidLength; //length of each SNPid in the BLOB

    public Pack2Alignment(IdGroup ids, Locus locus, int numSites, byte[][] data, NUM_BITS_PER_VALUE numBits, String[] states, byte[] variableSitesBLOB, byte[] SNPidBLOB) {
        super(ids, new TextDataType());
        myLocus = locus;
        myNumSites = numSites;
        myData = data;
        myStates = states;
        myNumBits = numBits;
        initVariableSites(variableSitesBLOB);
        initSNPids(SNPidBLOB);
        initMajorMinorAlleles();
    }

    private void initVariableSites(byte[] variableSitesBLOB) {
        ByteBuffer bb = ByteBuffer.wrap(variableSitesBLOB);
        bb.position(GdpdmBLOBUtils.numSitesField[0]);
        myNumSites = bb.getInt();
        myVariableSites = new int[getSiteCount()];
        bb.position(GdpdmBLOBUtils.totalHeaderPadding);
        for (int i = 0; i < getSiteCount(); i++) {
            myVariableSites[i] = bb.getInt();
        }
        myVariableSitesBLOBHeader = new byte[GdpdmBLOBUtils.totalHeaderPadding];
        bb.position(0);
        bb.get(myVariableSitesBLOBHeader, 0, GdpdmBLOBUtils.totalHeaderPadding);
    }

    private void initSNPids(byte[] SNPidBLOB) {
        ByteBuffer bb = ByteBuffer.wrap(SNPidBLOB);
        mySNPidBLOBHeader = new byte[GdpdmBLOBUtils.totalHeaderPadding];
        bb.get(mySNPidBLOBHeader, 0, GdpdmBLOBUtils.totalHeaderPadding);
        bb.position(GdpdmBLOBUtils.idLengthField);
        mySNPidLength = bb.getInt();
        mySNPid = new byte[SNPidBLOB.length - GdpdmBLOBUtils.totalHeaderPadding];
        bb.position(GdpdmBLOBUtils.totalHeaderPadding);
        bb.get(mySNPid);
    }

    private void initMajorMinorAlleles() {
        myMajorAlleles = new byte[getSiteCount()];
        myMinorAlleles = new byte[getSiteCount()];
        for (int i = 0; i < getSiteCount(); i++) {
            byte[] alleles = getAlleles(i);
            myMajorAlleles[i] = (alleles.length > 0) ? alleles[0] : DataType.UNKNOWN_BYTE;
            myMinorAlleles[i] = (alleles.length > 1) ? alleles[1] : DataType.UNKNOWN_BYTE;
        }
    }

    public byte[] getVariableSitesBLOB() {
        byte[] vsBLOB = new byte[GdpdmBLOBUtils.totalHeaderPadding + (4 * myVariableSites.length)];
        ByteBuffer bb = ByteBuffer.wrap(vsBLOB);
        bb.put(myVariableSitesBLOBHeader);
        for (int v : myVariableSites) {
            bb.putInt(v);
        }
        return vsBLOB;
    }

    public byte[] getSNPidBLOB() {
        byte[] SNPidBLOB = new byte[mySNPidBLOBHeader.length + mySNPid.length];
        ByteBuffer bb = ByteBuffer.wrap(SNPidBLOB);
        bb.put(mySNPidBLOBHeader);
        bb.put(mySNPid);
        return SNPidBLOB;
    }

    public char getBaseChar(int seq, int site) {
        return (char) getBase(seq, site);
    }

    public byte getBase(int seq, int site) {
        return AlignmentUtils.getAlleleValue(myData[seq], site, myNumBits);
    }

    public String getBaseString(int taxon, int site) {
        return myStates[getBase(taxon, site)];
    }

    @Override
    public byte getMajorAllele(int site) {
        return myMajorAlleles[site];
    }

    @Override
    public byte getMinorAllele(int site) {
        return myMinorAlleles[site];
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

    public int getSiteCount() {
        return myNumSites;
    }

    public int getSNPidLength() {
        return mySNPidLength;
    }

    public byte[] getSNPid() {
        return mySNPid;
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
            return "rsXXXXXXX";
        }
        return snpID;
    }

    public int getLocusSiteCount(Locus locus) {
        return getSiteCount();
    }

    public int getPositionInLocus(int site) {
        try {
            return myVariableSites[site];
        } catch (Exception e) {
            return site;
        }
    }

    public int getSiteOfPhysicalPosition(int physicalPosition, Locus locus) {
        try {
            return Arrays.binarySearch(myVariableSites, physicalPosition);
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
        return myLocus;
    }

    public Locus[] getLoci() {
        return new Locus[]{myLocus};
    }

    public int getNumLoci() {
        return 1;
    }

}