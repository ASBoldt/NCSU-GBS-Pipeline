/*
 * CombineAlignment
 */
package net.maizegenetics.pal.alignment;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import net.maizegenetics.pal.datatype.DataType;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.IdGroupUtils;

/**
 *
 * @author terry
 */
public class CombineAlignment extends AbstractAlignment {

    private static final long serialVersionUID = -5197800047652332969L;
    private final Alignment[] myAlignments;
    private final int[] mySiteOffsets;
    private final Map myLoci = new HashMap();

    private CombineAlignment(IdGroup subIdGroup, DataType dt, Alignment[] alignments) {
        super(subIdGroup, dt);
        myAlignments = alignments;
        mySiteOffsets = new int[alignments.length + 1];

        mySiteOffsets[0] = 0;
        int count = 0;
        for (int i = 0; i < alignments.length; i++) {
            count = alignments[i].getSiteCount() + count;
            mySiteOffsets[i + 1] = count;

            Locus[] loci = alignments[i].getLoci();
            for (int j = 0; j < loci.length; j++) {
                myLoci.put(loci[j], alignments[i]);
            }
        }

    }

    /**
     * This factory method combines given alignments.
     * It only one alignment, then it is returned unchanged.
     * Otherwise, this requires that each alignment
     * has the same Identifiers in the same order.
     *
     * @param alignments
     * @return
     */
    public static Alignment getInstance(Alignment[] alignments) {

        if ((alignments == null) || (alignments.length == 0)) {
            throw new IllegalArgumentException("CombineAlignment: getInstance: must provide alignments.");
        }

        if (alignments.length == 1) {
            return alignments[0];
        }

        IdGroup firstGroup = alignments[0].getIdGroup();
        for (int i = 1; i < alignments.length; i++) {
            if (!areIdGroupsEqual(firstGroup, alignments[i].getIdGroup())) {
                throw new IllegalArgumentException("CombineAlignment: getInstance: IdGroups do not match.");
            }
        }

        // Maybe need more logic here... -terry
        DataType dt = alignments[0].getDataType();

        return new CombineAlignment(firstGroup, dt, alignments);

    }

    /**
     * This factory method combines given alignments.
     * It only one alignment, then it is returned unchanged.
     * If isUnion equals true, a union join of the Identifiers
     * will be used to construct the combination.  Any alignment
     * not containing one of the Identifiers will return
     * unknown value for those locations.  If isUnion equals
     * false, a intersect join of the Identifiers will
     * be used.
     *
     * @param alignments alignments to combine
     * @param isUnion whether to union or intersect join
     * @return
     */
    public static Alignment getInstance(Alignment[] alignments, boolean isUnion) {

        if ((alignments == null) || (alignments.length == 0)) {
            throw new IllegalArgumentException("CombineAlignment: getInstance: must provide alignments.");
        }

        if (alignments.length == 1) {
            return alignments[0];
        }

        IdGroup[] groups = new IdGroup[alignments.length];
        for (int i = 0; i < alignments.length; i++) {
            groups[i] = alignments[i].getIdGroup();
        }
        IdGroup newTaxa = null;
        if (isUnion) {
            newTaxa = IdGroupUtils.getAllIds(groups);
        } else {
            newTaxa = IdGroupUtils.getCommonIds(groups);
        }

        Alignment[] newAlignments = new Alignment[alignments.length];
        for (int i = 0; i < alignments.length; i++) {
            newAlignments[i] = FilterAlignment.getInstance(alignments[i], newTaxa);
        }
        // Maybe need more logic here... -terry
        DataType dt = alignments[0].getDataType();

        return new CombineAlignment(newTaxa, dt, newAlignments);

    }

    private static boolean areIdGroupsEqual(IdGroup first, IdGroup second) {

        if (first.getIdCount() != second.getIdCount()) {
            return false;
        }

        for (int i = 0, n = first.getIdCount(); i < n; i++) {
            if (!first.getIdentifier(i).equals(second.getIdentifier(i))) {
                return false;
            }
        }

        return true;

    }

    public char getBaseChar(int taxon, int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getBaseChar(taxon, site - mySiteOffsets[translate]);
    }

    public byte getBase(int taxon, int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getBase(taxon, site - mySiteOffsets[translate]);
    }

    public byte[] getBaseRange(int taxon, int startSite, int endSite) {

        int siteCount = getSiteCount();
        if ((startSite == 0) && (endSite == siteCount - 1)) {
            byte[] result = new byte[siteCount];
            int count = 0;
            for (int i = 0; i < myAlignments.length; i++) {
                int currentNumSites = myAlignments[i].getSiteCount();
                for (int j = 0; j < currentNumSites; j++) {
                    result[count++] = myAlignments[i].getBase(taxon, j);
                }
            }
            return result;
        } else {
            return super.getBaseRange(taxon, startSite, endSite);
        }

    }

    public byte getBase(int taxon, int site, int allele) {
        int translate = translateSite(site);
        return myAlignments[translate].getBase(taxon, site - mySiteOffsets[translate], allele);
    }

    public byte getBase(int taxon, Locus locus, int physicalPosition, int allele) {
        int site = getSiteOfPhysicalPosition(physicalPosition, locus);
        int translate = translateSite(site);
        return myAlignments[translate].getBase(taxon, site - mySiteOffsets[translate], allele);
    }

    /**
     * Returns which alignment to use.
     *
     * @param site
     * @return alignment index.
     */
    public int translateSite(int site) {

        for (int i = 1; i < mySiteOffsets.length; i++) {
            if (mySiteOffsets[i] > site) {
                return i - 1;
            }
        }
        throw new IndexOutOfBoundsException("CombineAlignment: translateSite: index out of range: " + site);

    }

    public boolean hasReference() {

        for (int i = 0; i < myAlignments.length; i++) {
            if (!myAlignments[i].hasReference()) {
                return false;
            }
        }

        return true;
    }

    public String[] getSNPIDs() {

        int numSites = getSiteCount();
        String[] result = new String[numSites];
        int count = 0;
        for (int i = 0; i < myAlignments.length; i++) {
            for (int j = 0, n = myAlignments[i].getSiteCount(); j < n; j++) {
                result[count++] = myAlignments[i].getSNPID(j);
            }
        }

        return result;

    }

    public String getSNPID(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getSNPID(site - mySiteOffsets[translate]);
    }

    public int getSiteCount() {
        return mySiteOffsets[mySiteOffsets.length - 1];
    }

    public int getLocusSiteCount(Locus locus) {
        return ((Alignment) myLoci.get(locus)).getLocusSiteCount(locus);
    }

    public int getPositionInLocus(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getPositionInLocus(site - mySiteOffsets[translate]);
    }

    public int getSiteOfPhysicalPosition(int physicalPosition, Locus locus) {
        Alignment align = ((Alignment) myLoci.get(locus));
        int i = -1;
        for (int j = 0; j < myAlignments.length; j++) {
            if (myAlignments[j] == align) {
                i = j;
                break;
            }
        }
        if (i == -1) {
            return -1;
        }
        return mySiteOffsets[i] + align.getSiteOfPhysicalPosition(physicalPosition, locus);
    }

    public byte getPositionType(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getPositionType(site - mySiteOffsets[translate]);
    }

    public byte[] getPositionTypes() {

        int numSites = getSiteCount();
        byte[] result = new byte[numSites];
        for (int i = 0; i < numSites; i++) {
            int translate = translateSite(i);
            result[i] = myAlignments[translate].getPositionType(i - mySiteOffsets[translate]);
        }

        return result;

    }

    public Locus getLocus(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getLocus(site - mySiteOffsets[translate]);
    }

    public Locus[] getLoci() {

        List<Locus> loci = new ArrayList();
        for (int i = 0; i < myAlignments.length; i++) {
            loci.addAll(Arrays.asList(myAlignments[i].getLoci()));
        }

        Locus[] result = new Locus[loci.size()];
        return loci.toArray(result);

    }

    public int getNumLoci() {
        return getLoci().length;
    }

    public float[][] getSiteScores() {

        if (!hasSiteScores()) {
            return null;
        }

        int numSeqs = getSequenceCount();
        float[][] result = new float[numSeqs][getSiteCount()];
        for (int a = 0, n = myAlignments.length; a < n; a++) {
            if (myAlignments[a].hasSiteScores()) {
                for (int s = 0, m = myAlignments[a].getSiteCount(); s < m; s++) {
                    for (int t = 0; t < numSeqs; t++) {
                        result[t][mySiteOffsets[a] + s] = myAlignments[a].getSiteScore(t, s);
                    }
                }
            }
        }

        return result;

    }

    public float getSiteScore(int seq, int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getSiteScore(seq, site - mySiteOffsets[translate]);
    }

    public boolean hasSiteScores() {
        for (Alignment align : myAlignments) {
            if (align.hasSiteScores()) {
                return true;
            }
        }
        return false;
    }

    public int getIndelSize(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getIndelSize(site - mySiteOffsets[translate]);
    }

    public boolean isIndel(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].isIndel(site - mySiteOffsets[translate]);
    }

    public byte getReferenceAllele(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getReferenceAllele(site - mySiteOffsets[translate]);
    }

    public Alignment[] getAlignments() {
        return myAlignments;
    }

    public byte getMajorAllele(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getMajorAllele(site - mySiteOffsets[translate]);
    }

    public byte getMinorAllele(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getMinorAllele(site - mySiteOffsets[translate]);
    }

    public byte[] getMinorAlleles(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getMinorAlleles(site - mySiteOffsets[translate]);
    }

    public byte[] getAlleles(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getAlleles(site - mySiteOffsets[translate]);
    }

    public double getMinorAlleleFrequency(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getMinorAlleleFrequency(site - mySiteOffsets[translate]);
    }

    public int[][] getAllelesSortedByFrequency(int site, boolean includeGAP) {
        int translate = translateSite(site);
        return myAlignments[translate].getAllelesSortedByFrequency(site - mySiteOffsets[translate], includeGAP);
    }
}
