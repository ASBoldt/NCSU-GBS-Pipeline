// SitePattern.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)
// Known bugs and limitations:
// - computational complexity O(numSeqs*numSites)
package net.maizegenetics.pal.alignment;

import net.maizegenetics.pal.datatype.DataType;
import net.maizegenetics.pal.datatype.TextDataType;
import net.maizegenetics.pal.ids.IdGroup;

/**
 * takes an Alignment and determines its site patterns
 *
 * @version $Id: SitePattern.java,v 1.2 2009/07/24 06:24:41 tcasstevens Exp $
 *
 * @author Korbinian Strimmer
 * @author Alexei Drummond
 */
public class SitePattern extends AbstractAlignment {

    /** number of site patterns */
    private int myNumPatterns;
    /** site -> site pattern */
    private int[] myAlias;
    /** weights of each site pattern */
    private int[] myWeight;
    /** site patterns [sequence][site pattern] */
    private byte[][] myPattern;
    /** Number of sites */
    private final int myNumSites;

    /**
     * infer site patterns for a given alignment
     *
     * @param a alignment
     */
    public SitePattern(Alignment a) {
        super(a);
        DataType dt = a.getDataType();
        if (dt == null) {
            System.out.println("Warning: Input alignment for SitePattern has null datatype");
            dt = new TextDataType();
        }

        myNumSites = a.getSiteCount();
        makeSitePattern(a);

    }

    /**
     * construct SitePattern from scratch
     *
     * @param dataType data type
     * @param numSites number of sites
     * @param numSeqs number of sequences
     * @param idGroup sequence identifiers
     * @param numPatterns number of site patterns
     * @param alias link site -> site pattern
     * @param weight frequency of a site pattern
     * @param pattern site patterns
     */
    public SitePattern(DataType dataType, int numSites, int numSeqs, IdGroup idGroup,
            int numPatterns, int[] alias, int[] weight, byte[][] pattern) {
        super(idGroup, dataType);
        myNumSites = numSites;

        myNumPatterns = numPatterns;
        myAlias = alias;
        myWeight = weight;
        myPattern = pattern;

        AlignmentUtils.estimateFrequencies(this);  // Bernard Suh <bsuh@tigr.org>
    }

    /**
     *
     * @param a An alignment
     * @return alignment as a site pattern if it isn't already one (other wise just returns alighnment)
     */
    public static final SitePattern getSitePattern(Alignment a) {
        if (a instanceof SitePattern) {
            return (SitePattern) a;
        }
        return new SitePattern(a);
    }

    public char getBaseChar(int seq, int site) {
        return getDataType().getChar(myPattern[seq][myAlias[site]]);
    }

    public final char getPatternData(int seq, int patternSite) {
        return getDataType().getChar(myPattern[seq][patternSite]);
    }

    public final int getPatternState(int seq, int patternSite) {
        return myPattern[seq][patternSite];
    }

    /**
     * Accessor method for weight
     */
    public int[] getSiteWeights() {
        return myWeight;
    }

    /**
     * Accessor method for numPatterns
     */
    public int getNumberOfPatterns() {
        return myNumPatterns;
    }

    private int[] patSort;

    private void makeSitePattern(Alignment alignment) {
        myAlias = new int[getSiteCount()];

        if (getSequenceCount() > 0 && getSiteCount() > 0) {
            patSort = new int[getSiteCount()];

            getNumpattern(alignment);
            myPattern = new byte[getSequenceCount()][myNumPatterns];
            myWeight = new int[myNumPatterns];
            copypattern(alignment);

            patSort = null;
        } else {
            myNumPatterns = 0;
            myPattern = null;
            myWeight = null;
        }
    }

    private int stateData(Alignment al, int seq, int site) {
        int state = getDataType().getState(al.getBaseChar(seq, site));
        if (getDataType().isUnknownState(state)) {
            return getDataType().getNumStates();
        }
        return state;
    }

    private void getNumpattern(Alignment alignment) {
        int tpmradix = getDataType().getNumStates();

        int[] awork = new int[getSiteCount()];
        int[] count = new int[tpmradix + 1];

        for (int j = 0; j < getSiteCount(); j++) {
            patSort[j] = j;
        }
        for (int i = getSequenceCount() - 1; i >= 0; i--) {
            for (int k = 0; k < tpmradix + 1; k++) {
                count[k] = 0;
            }
            for (int j = 0; j < getSiteCount(); j++) {
                count[stateData(alignment, i, patSort[j])]++;
            }
            for (int k = 1; k < tpmradix + 1; k++) {
                count[k] += count[k - 1];
            }
            for (int j = getSiteCount() - 1; j >= 0; j--) {
              awork  [ --count[stateData(alignment, i, patSort[j])] ] = patSort[j];
			}
			for (int j = 0; j < getSiteCount(); j++) {
                patSort[j] = awork[j];
            }
        }
        awork = null;
        count = null;

        myNumPatterns = 1;
        for (int j = 1; j < getSiteCount(); j++) {
            int s = patSort[j - 1];
            int t = patSort[j];
            for (int i = 0; i < getSequenceCount(); i++) {
                if (stateData(alignment, i, t) !=
                        stateData(alignment, i, s)) {
                    myNumPatterns++;
                    break;
                }
            }
        }
    }

    void copypattern(Alignment alignment) {
        int k, n;
        boolean isSame;

        n = 0;
        k = patSort[n];
        for (int i = 0; i < getSequenceCount(); i++) {
            myPattern[i][n] = (byte) stateData(alignment, i, k);
        }
        myWeight[n] = 1;
        myAlias[k] = 0;

        for (int j = 1; j < getSiteCount(); j++) {
            k = patSort[j];

            isSame = true;
            for (int i = 0; i < getSequenceCount(); i++) {
                if (myPattern[i][n] != (byte) stateData(alignment, i, k)) {
                    isSame = false;
                    break;
                }
            }

            if (isSame) {
                myWeight[n]++;
                myAlias[k] = n;
            } else {
                n++;
                for (int i = 0; i < getSequenceCount(); i++) {
                    myPattern[i][n] = (byte) stateData(alignment, i, k);
                }
                myWeight[n] = 1;
                myAlias[k] = n;
            }
        }
    }

    public byte getBase(int taxon, int site) {
        return (byte) getBaseChar(taxon, site);
    }

    public byte getBase(int taxon, int site, int allele) {
        return (byte) getBaseChar(taxon, site);
    }

    public byte getBase(int taxon, Locus locus, int physicalPosition, int allele) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public boolean hasReference() {
        return false;
    }

    public String[] getSNPIDs() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public String getSNPID(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public int getSiteCount() {
        return myNumSites;
    }

    public int getLocusSiteCount(Locus locus) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public int getPositionInLocus(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public int getSiteOfPhysicalPosition(int physicalPosition, Locus locus) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public byte getPositionType(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public byte[] getPositionTypes() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public Locus getLocus(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public Locus[] getLoci() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public int getNumLoci() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public int getIndelSize(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public boolean isIndel(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public byte getReferenceAllele(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public int getNumPatterns() {
        return myNumPatterns;
    }

    public int[] getWeight() {
        return myWeight;
    }

    public void setWeight(int[] weight) {
        myWeight = weight;
    }

    public byte[][] getPattern() {
        return myPattern;
    }

    public int[] getAlias() {
        return myAlias;
    }
}
