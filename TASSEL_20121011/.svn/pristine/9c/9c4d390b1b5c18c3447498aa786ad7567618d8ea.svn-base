/*
SimpleAlignment.java
 */
package net.maizegenetics.pal.alignment;

import java.util.Arrays;
import java.util.TreeSet;

import net.maizegenetics.pal.datatype.DataType;
import net.maizegenetics.pal.datatype.TextDataType;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.SimpleIdGroup;

/**
 */
public class SimpleAlignment extends AbstractAlignment {

    private static final long serialVersionUID = 4303224913340358191L;
    /** The Sequences */
    private final String[] mySequences;
    /** Physical Positions */
    private final int[] myPhysicalPositions;
    /** Position Types */
    private final byte[] myPositionTypes;
    /** Locus */
    private final Locus myLocus;
    /** Number of sites */
    private final int myNumSites;
    /** Quality Scores */
    private final float[][] mySiteScores;
    /** Site Names */
    private final String[] mySiteNames;

    /**
     * Constructor
     */
    public SimpleAlignment(IdGroup group, String[] sequences, DataType dt, byte[] reference, int[] physicalPositions, byte[] positionTypes, Locus locus, float[][] siteScores, String[] siteNames, boolean clean) {

        super(group, dt);

        if ((sequences == null) || (sequences.length == 0) || (sequences[0].length() == 0)) {
            throw new IllegalArgumentException("SimpleAlignment: init: must have sequence data.");
        }
        mySequences = sequences;
        myNumSites = sequences[0].length();
        for (int i = 0; i < sequences.length; i++) {
            if ((clean) && !(dt instanceof TextDataType)) mySequences[i] = AlignmentUtils.cleanSequenceString(mySequences[i]);
            if (myNumSites != sequences[i].length()) {
                throw new IllegalArgumentException("SimpleAlignment: init: sequences differ in length.");
            }
        }

        setReference(reference);

        if ((physicalPositions != null) && (physicalPositions.length != 0) && (myNumSites != physicalPositions.length)) {
            throw new IllegalArgumentException("SimpleAlignment: init: physical positions differs in length with sequences.");
        }
        myPhysicalPositions = physicalPositions;

        if ((positionTypes != null) && (positionTypes.length != 0) && (myNumSites != positionTypes.length)) {
            throw new IllegalArgumentException("SimpleAlignment: init: position types differs in length with sequences.");
        }
        myPositionTypes = positionTypes;

        myLocus = locus;

        if ((siteScores != null) && (siteScores.length != 0) && (getSequenceCount() != siteScores.length) &&
                (siteScores[0].length != 0) && (myNumSites != siteScores[0].length)) {
            throw new IllegalArgumentException("SimpleAlignment: init: quality score matrix differs in size of sequences.");
        }
        mySiteScores = siteScores;

        if ((siteNames == null) || (siteNames.length == 0)) {
            mySiteNames = null;
        } else if (siteNames.length != myNumSites) {
            throw new IllegalArgumentException("SimpleAlignment: init: number of site names should equal number of sites.");
        } else {
            mySiteNames = siteNames;
        }

        validateInstance();

    }

    /**
     * Constructor
     */
    public SimpleAlignment(IdGroup group, String[] sequences, DataType dt, byte[] reference, int[] physicalPositions, byte[] positionTypes, Locus locus, float[][] siteScores, String[] siteNames) {
    	this(group, sequences, dt, reference, physicalPositions, positionTypes, locus, siteScores, siteNames, true);
    }

    /**
     * Constructor
     */
    public SimpleAlignment(IdGroup group, String[] sequences, DataType dt, byte[] reference, int[] physicalPositions, byte[] positionTypes, String locus, float[][] siteScores, String[] siteNames) {
        this(group, sequences, dt, reference, physicalPositions, positionTypes, new Locus(locus, null, -1, -1, null, null), siteScores, siteNames);
    }

    /**
     * Clone constructor for Alignment
     */
    public static SimpleAlignment getInstance(Alignment a) {

        IdGroup group = a.getIdGroup();
        DataType dt = a.getDataType();
        String[] sequences = SimpleAlignment.getSequences(a);
        byte[] reference = a.getReference();
        int[] physicalPositions = getPhysicalPositions(a);
        byte[] positionTypes = a.getPositionTypes();
        Locus locus = a.getLocus(0);
        float[][] siteScores = a.getSiteScores();
        String[] siteNames = a.getSNPIDs();
        return new SimpleAlignment(group, sequences, dt, reference, physicalPositions, positionTypes, locus, siteScores, siteNames);

    }

    public static SimpleAlignment getInstance(IdGroup group, String[] sequences, DataType dt) {
        return new SimpleAlignment(group, sequences, dt, null, new int[0], new byte[0], "", new float[0][0], null);
    }

    /**
     * If you have known positions (or position offsets compensating for gaps) that you want to use (added by Jeff)
     */
    public static SimpleAlignment getInstance(IdGroup group, String[] sequences, DataType dt, int[] physicalPositions) {
        return new SimpleAlignment(group, sequences, dt, null, physicalPositions, new byte[0], "", new float[0][0], null);
    }

    /**
     * Constructor for defined ids and sequences, but using an existing annotation
     * TODO - I think this should be eliminated, this should be a filter class
     */
    public static SimpleAlignment getInstance(IdGroup group, String[] sequences, Alignment a) {

        DataType dt = a.getDataType();
        byte[] reference = a.getReference();
        int[] physicalPositions = getPhysicalPositions(a);
        byte[] positionTypes = a.getPositionTypes();
        Locus locus = a.getLocus(0);
        float[][] siteScores = a.getSiteScores();
        String[] siteNames = a.getSNPIDs();
        return new SimpleAlignment(group, sequences, dt, reference, physicalPositions, positionTypes, locus, siteScores, siteNames);

    }

    /**
     * This constructor will subset the alignment based on the taxa in IdGroup.
     * It does the intersection just to make sure there isn't something missing.
     * TODO - I think this should be eliminated, this should be a filter class
     */
    public static SimpleAlignment getInstance(Alignment a, IdGroup subGroup) {

        TreeSet<String> a1IDs = new TreeSet<String>();
        TreeSet<String> a2IDs = new TreeSet<String>();
        IdGroup oldIds = a.getIdGroup();
        for (int i = 0; i < oldIds.getIdCount(); i++) {
            a1IDs.add(oldIds.getIdentifier(i).getName());
        }
        for (int i = 0; i < subGroup.getIdCount(); i++) {
            a2IDs.add(subGroup.getIdentifier(i).getName());
        }
        TreeSet<String> finalIDs = new TreeSet<String>(a1IDs);
        finalIDs.retainAll(a2IDs); //intersection of the two sets
        String[] IDs = new String[finalIDs.size()];
        finalIDs.toArray(IDs);
        IdGroup newGroup = new SimpleIdGroup(IDs);

        String[] s = new String[newGroup.getIdCount()];
        for (int i = 0; i < newGroup.getIdCount(); i++) {
            int oldI = oldIds.whichIdNumber(newGroup.getIdentifier(i).getName());
            s[i] = a.getAlignedSequenceString(oldI);
        }
        DataType dt = a.getDataType();
        byte[] reference = a.getReference();
        int[] physicalPositions = getPhysicalPositions(a);
        byte[] positionTypes = a.getPositionTypes();
        Locus locus = a.getLocus(0);
        float[][] siteScores = a.getSiteScores();
        String[] siteNames = a.getSNPIDs();
        return new SimpleAlignment(newGroup, s, dt, reference, physicalPositions, positionTypes, locus, siteScores, siteNames);

    }

    public static SimpleAlignment getInstance(IdGroup group, char[][] cSequences, DataType dt) {

        String[] sequences = new String[cSequences.length];
        for (int i = 0; i < sequences.length; i++) {
            sequences[i] = new String(cSequences[i]);
        }
        return getInstance(group, sequences, dt);

    }

    private void validateInstance() {
        if (getPositionInLocus(0) > getPositionInLocus(getSiteCount() - 1)) {
            throw new IllegalStateException("SimpleAlignment: validateInstance: Start position is larger than end position.");
        }
    }

    public static int[] getPhysicalPositions(Alignment a) {
        int[] result = new int[a.getSiteCount()];
        for (int i = 0, n = a.getSiteCount(); i < n; i++) {
            result[i] = a.getPositionInLocus(i);
        }
        return result;
    }

    public static String[] getSequences(Alignment a) {

        String[] s = new String[a.getSequenceCount()];
        for (int i = 0; i < a.getSequenceCount(); i++) {
            s[i] = a.getAlignedSequenceString(i);
        }

        return s;

    }

    public int getPositionInLocus(int site) {
        try {
            return myPhysicalPositions[site];
        } catch (Exception e) {
            return site;
        }
    }

    public byte getPositionType(int site) {
        if ((myPositionTypes == null) || (myPositionTypes.length == 0)) {
            return PositionType.ALL_GROUP;
        }
        return myPositionTypes[site];
    }

    public String getLocusName(int site) {
        return myLocus.getName();
    }

    public float getSiteScore(int seq, int site) {

        if ((mySiteScores == null) || mySiteScores.length == 0) {
            return -9;
        } else {
            return mySiteScores[seq][site];
        }

    }

    public boolean hasSiteScores() {
        return ((mySiteScores != null) && (mySiteScores.length != 0));
    }

    public String getFormattedSequenceString(int seq) {

        StringBuffer sb = new StringBuffer();
        DataType dt = getDataType();
        for (int i = 0; i < getSiteCount(); i++) {
            int dataLen = dt.getMaxFormatedStringLength();
            String s2 = dt.getFormattedString(getBaseChar(seq, i));
            String s1 = "    ".substring(0, dataLen - s2.length()) + s2;
            sb.append(s1);
        }

        return sb.toString();

    }

    /**
     * get the number of the columns
     *
     * @return columns names
     */
    public int getColumnCount() {
        return getSiteCount() + 1;
    }

    /**
     * get the number of columns
     *
     * @return columns names
     */
    public int getRowCount() {
        return getSequenceCount();
    }

    /**
     * get the total number of elements in the dataset. Elements=rowCount * columnCount;
     *
     * @return columns names
     */
    public int getElementCount() {
        return getRowCount() * getColumnCount();
    }

    public char getBaseChar(int seq, int site) {
        return mySequences[seq].charAt(site);
    }

    /**
     * Returns a string representing a single sequence (including gaps)
     * from this alignment.
     */
    public String getAlignedSequenceString(int seq) {
        return mySequences[seq];
    }

    public byte[] getPositionTypes() {
        return myPositionTypes;
    }

    public float[][] getSiteScores() {
        return mySiteScores;
    }

    public byte getBase(int taxon, int site) {
        return (byte) getBaseChar(taxon, site);
    }

    /**
     * Returns value for given taxon , site.
     *
     * @param taxon taxon
     * @param site site
     * @param allele Not used. Only one allele per site.
     *
     * @return value
     */
    public byte getBase(int taxon, int site, int allele) {
        return (byte) getBaseChar(taxon, site);
    }

    /**
     * Returns value for given taxon, locus, physical
     * position, and allele.
     *
     * @param taxon taxon
     * @param locus locus
     * @param physicalPosition physical position
     * @param allele Not used. Only one allele per site.
     *
     * @return value
     */
    public byte getBase(int taxon, Locus locus, int physicalPosition, int allele) {

        if (locus != myLocus) {
            throw new IllegalArgumentException("SimpleAlignment: getBase: locus not defined.");
        }

        int site = Arrays.binarySearch(myPhysicalPositions, physicalPosition);

        return (byte) getBaseChar(taxon, site);

    }

    public int getSiteCount() {
        return myNumSites;
    }

    public int getLocusSiteCount(Locus locus) {

        if (locus == myLocus) {
            return myNumSites;
        } else {
            return 0;
        }

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

    public String[] getSNPIDs() {
        return mySiteNames;
    }

    public String getSNPID(int site) {
        if ((mySiteNames == null) || (mySiteNames.length == 0)) {
            return null;
        }
        return mySiteNames[site];
    }

    public int getSiteOfPhysicalPosition(int physicalPosition, Locus locus) {
        if ((myPhysicalPositions == null) || (myPhysicalPositions.length == 0)) {
            return physicalPosition;
        }
        return Arrays.binarySearch(myPhysicalPositions, physicalPosition);
    }
}
