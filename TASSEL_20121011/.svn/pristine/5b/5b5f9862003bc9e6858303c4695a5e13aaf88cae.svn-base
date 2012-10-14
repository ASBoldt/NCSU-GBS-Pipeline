/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.pal.alignment;

import java.nio.ByteBuffer;
import java.nio.IntBuffer;
import net.maizegenetics.pal.ids.IdGroup;

/**
 *
 * @author ed
 */
public class Pack1AlignmentWithCoverage extends Pack1Alignment {

    protected int[][] coveredSites = null; //null unless covered sites loaded; arrays of start stop positions for each taxon
    private boolean isOnlyPolymorphic = false;
    private int minPosition,  maxPosition;

    //One thought is to have a inflated cache of 1M elements.  We can inflate 300M elements in 2s, so this would not
    //be very expensive and provide high performance for analyses in smaller ranges.
    public Pack1AlignmentWithCoverage(byte[][] alleleBLOB, byte[] variableSitesBLOB, byte[] SNPidBLOB, byte[][] coverageBLOB, byte[] refSeqBLOB) {
        super(alleleBLOB, variableSitesBLOB, SNPidBLOB);
        initCoveredSites(coverageBLOB);
        //todo use coverage to calc numSites, minPosition, maxPosition
    }

    private void initCoveredSites(byte[][] coverageBLOB) {
        minPosition = Integer.MAX_VALUE;
        maxPosition = Integer.MIN_VALUE;
        IdGroup coverageTaxa = AllelePositionBLOBUtils.getTaxa(coverageBLOB);
        if (this.getSequenceCount() != coverageTaxa.getIdCount()) {
            throw new IllegalArgumentException("Pack1AlignmentWithCoverage: initCoveredSites: number of taxa do not agree");
        }
        coveredSites = new int[coverageBLOB.length][];
        for (int t = 0; t < coverageTaxa.getIdCount(); t++) {
            int currTaxon = this.getIdGroup().whichIdNumber(coverageTaxa.getIdentifier(t).getFullName());
            if (currTaxon < 0) {
                throw new IllegalArgumentException("Pack1AlignmentWithCoverage: initCoveredSites: taxa names do not agree");
            }
            ByteBuffer bb = ByteBuffer.wrap(coverageBLOB[t]);
            bb.position(GdpdmBLOBUtils.numSitesField[0]);
            int segments = bb.getInt();
            IntBuffer ib = ByteBuffer.wrap(coverageBLOB[t], GdpdmBLOBUtils.totalHeaderPadding,
                    coverageBLOB[t].length - GdpdmBLOBUtils.totalHeaderPadding).asIntBuffer();
            coveredSites[currTaxon] = new int[segments * 2];
            for (int i = 0; i < coveredSites[currTaxon].length; i++) {
                coveredSites[currTaxon][i] = ib.get();
            }
            if (coveredSites[currTaxon][0] < minPosition) {
                minPosition = coveredSites[currTaxon][0];
            }
            if (coveredSites[currTaxon][coveredSites[currTaxon].length - 1] > maxPosition) {
                maxPosition = coveredSites[currTaxon][coveredSites[currTaxon].length - 1];
            }
        }
        numSites = maxPosition - minPosition + 1;
    }

    public char getDataChar(int seq, int site) {
        //todo use coverage
        return (char) getData(seq, site);
    }

    public byte getData(int seq, int site) {
        int position = site + minPosition;
        if (CoverageBLOBUtils.isCovered(coveredSites[seq], position)) {
            int varSite = this.getSiteOfPhysicalPosition(position, null);
            if (varSite > -1) {
                return AllelePositionBLOBUtils.getBaseFromAlleleBLOB(alleleBLOB[seq], varSite);
            }
            return 'R';
        } else {
            return GdpdmBLOBUtils.UNKNOWN_CHARACTER;
        }

        // if (isOnlyPolymorphic) {
        //    return AllelePositionBLOBUtils.getBaseFromAlleleBLOB(alleleBLOB[seq], site);
        // }
        // throw new UnsupportedOperationException("Not supported yet.");
    }

    public byte[] getDataRange(int seq, int startSite, int endSite) {
        //todo use coverage
        if (isOnlyPolymorphic) {
            return AllelePositionBLOBUtils.getBaseRangeFromAlleleBLOB(alleleBLOB[seq], startSite, endSite);
        }
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public int getSiteCount() {
        //todo use coverage
        return numSites;
    }

    public byte getData(int taxon, int site, int allele) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public byte getData(int taxon, Locus locus, int physicalPosition, int allele) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public boolean hasReference() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public String[] getSNPIDs() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public String getSNPID(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public int getLocusSiteCount(Locus locus) {
        return getSiteCount();
    }

    public int getPositionInLocus(int site) {
        return variableSites[site];
    }

    public int[] getPhysicalPositions() {
        return variableSites;
    }

    public byte getPositionType(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public byte[] getPositionTypes() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public Locus getLocus(int site) {
        return locus;
    }

    public byte getReferenceAllele(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
