/*
 * FilterAlignment
 */
package net.maizegenetics.pal.alignment;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import net.maizegenetics.pal.datatype.DataType;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;

/**
 * All taxa and site filtering should be controlled through this class.
 * It essentially creates views of the baseAlignment
 * @author terry
 */
public class FilterAlignment extends AbstractAlignment {

    private static final long serialVersionUID = -5197800047652332969L;
    private final boolean myIsTaxaFilter;
    private final boolean myIsSiteFilter;
    private final boolean myIsSiteFilterByRange;
    private final Alignment myBaseAlignment;
    private final int[] myTaxaRedirect;
    private final int[] mySiteRedirect;
    private final int myRangeStart;
    private final int myRangeEnd;
    private final Locus[] myLoci;

    private FilterAlignment(Alignment a, IdGroup subIdGroup, int[] taxaRedirect, FilterAlignment original) {

        super(subIdGroup, a.getDataType());

        if (subIdGroup.getIdCount() != taxaRedirect.length) {
            throw new IllegalArgumentException("FilterAlignment: init: subIdGroup should be same size as taxaRedirect.");
        }

        myIsTaxaFilter = true;
        myBaseAlignment = a;
        myTaxaRedirect = taxaRedirect;

        if (original == null) {
            myIsSiteFilter = false;
            myIsSiteFilterByRange = false;
            mySiteRedirect = null;
            myRangeStart = -1;
            myRangeEnd = -1;
            myLoci = myBaseAlignment.getLoci();
        } else {
            myIsSiteFilter = original.isSiteFilter();
            myIsSiteFilterByRange = original.isSiteFilterByRange();
            mySiteRedirect = original.getSiteRedirect();
            myRangeStart = original.getRangeStart();
            myRangeEnd = original.getRangeEnd();
            myLoci = original.getLoci();
        }

    }

    /**
     * This returns FilterAlignment with only specified subIdGroup.  Defaults
     * to retain unknown taxa.
     * 
     * @param a alignment
     * @param subIdGroup subset id group
     * 
     * @return filter alignment 
     */
    public static FilterAlignment getInstance(Alignment a, IdGroup subIdGroup) {
        return getInstance(a, subIdGroup, true);
    }

    /**
     * This returns FilterAlignment with only specified subIdGroup.  If
     * retainUnknownTaxa is true then Alignment will return unknown values
     * for missing taxa.
     * 
     * @param a alignment
     * @param subIdGroup subset id group
     * @param retainUnknownTaxa whether to retain unknown taxa
     * 
     * @return filter alignment 
     */
    public static FilterAlignment getInstance(Alignment a, IdGroup subIdGroup, boolean retainUnknownTaxa) {

        Alignment baseAlignment = a;
        FilterAlignment original = null;
        if (baseAlignment instanceof FilterAlignment) {
            original = (FilterAlignment) a;
            baseAlignment = ((FilterAlignment) a).getBaseAlignment();
        }

        List<Integer> taxaRedirectList = new ArrayList<Integer>();
        List<Identifier> idList = new ArrayList<Identifier>();
        for (int i = 0; i < subIdGroup.getIdCount(); i++) {
            int ion = a.getIdGroup().whichIdNumber(subIdGroup.getIdentifier(i));
            if ((retainUnknownTaxa) && (ion < 0)) {
                taxaRedirectList.add(-1);
                idList.add(subIdGroup.getIdentifier(i));
            } else {
                int idn = baseAlignment.getIdGroup().whichIdNumber(subIdGroup.getIdentifier(i));
                if (idn > -1) {
                    taxaRedirectList.add(idn);
                    idList.add(baseAlignment.getIdGroup().getIdentifier(idn));
                } else if (retainUnknownTaxa) {
                    taxaRedirectList.add(-1);
                    idList.add(subIdGroup.getIdentifier(i));
                }
            }
        }

        int[] taxaRedirect = new int[taxaRedirectList.size()];
        for (int j = 0, n = taxaRedirectList.size(); j < n; j++) {
            taxaRedirect[j] = (int) taxaRedirectList.get(j);
        }

        Identifier[] ids = new Identifier[idList.size()];
        idList.toArray(ids);
        IdGroup resultIdGroup = new SimpleIdGroup(ids);

        return new FilterAlignment(baseAlignment, resultIdGroup, taxaRedirect, original);

    }

    /**
     * Removes specified IDs.
     *
     * @param a alignment to filter
     * @param subIdGroup specified IDs
     *
     * @return Filtered Alignment
     */
    public static FilterAlignment getInstanceRemoveIDs(Alignment a, IdGroup subIdGroup) {

        List result = new ArrayList();
        IdGroup current = a.getIdGroup();
        for (int i = 0, n = current.getIdCount(); i < n; i++) {
            if (subIdGroup.whichIdNumber(current.getIdentifier(i)) == -1) {
                result.add(current.getIdentifier(i));
            }
        }
        Identifier[] ids = new Identifier[result.size()];
        result.toArray(ids);
        return FilterAlignment.getInstance(a, new SimpleIdGroup(ids));

    }

    /**
     * Constructor
     *
     * @param a base alignment
     * @param startSite start site (included)
     * @param endSite end site (included)
     */
    private FilterAlignment(Alignment a, int startSite, int endSite, FilterAlignment original) {

        super(original == null ? a.getIdGroup() : original.getIdGroup(), a.getDataType());

        if (startSite > endSite) {
            throw new IllegalArgumentException("FilterAlignment: init: start site: " + startSite + " is larger than end site: " + endSite);
        }

        if ((startSite < 0) || (startSite > a.getSiteCount() - 1)) {
            throw new IllegalArgumentException("FilterAlignment: init: start site: " + startSite + " is out of range.");
        }

        if ((endSite < 0) || (endSite > a.getSiteCount() - 1)) {
            throw new IllegalArgumentException("FilterAlignment: init: end site: " + endSite + " is out of range.");
        }

        myBaseAlignment = a;
        myIsSiteFilterByRange = true;
        myIsSiteFilter = false;
        myRangeStart = startSite;
        myRangeEnd = endSite;
        mySiteRedirect = null;
        myLoci = getLociFromBase();

        if (original == null) {
            myIsTaxaFilter = false;
            myTaxaRedirect = null;
        } else {
            myIsTaxaFilter = original.isTaxaFilter();
            myTaxaRedirect = original.getTaxaRedirect();
        }

    }

    /**
     * Constructor
     *
     * @param a base alignment
     * @param subSites site to include
     */
    private FilterAlignment(Alignment a, int[] subSites, FilterAlignment original) {

        super(original == null ? a.getIdGroup() : original.getIdGroup(), a.getDataType());

        myBaseAlignment = a;
        myIsSiteFilter = true;
        myIsSiteFilterByRange = false;
        if ((subSites == null) || (subSites.length == 0)) {
            mySiteRedirect = new int[0];
        } else {
            mySiteRedirect = new int[subSites.length];
            Arrays.sort(subSites);
            System.arraycopy(subSites, 0, mySiteRedirect, 0, subSites.length);
        }
        myRangeStart = -1;
        myRangeEnd = -1;
        myLoci = getLociFromBase();

        if (original == null) {
            myIsTaxaFilter = false;
            myTaxaRedirect = null;
        } else {
            myIsTaxaFilter = original.isTaxaFilter();
            myTaxaRedirect = original.getTaxaRedirect();
        }

    }

    public static FilterAlignment getInstance(Alignment a, int[] subSites) {

        if (a instanceof FilterAlignment) {
            FilterAlignment original = (FilterAlignment) a;
            Alignment baseAlignment = ((FilterAlignment) a).getBaseAlignment();
            if (original.isSiteFilter()) {
                int[] newSubSites = new int[subSites.length];
                for (int i = 0; i < subSites.length; i++) {
                    newSubSites[i] = original.translateSite(subSites[i]);
                }
                return new FilterAlignment(baseAlignment, newSubSites, original);
            } else if (original.isSiteFilterByRange()) {
                int[] newSubSites = new int[subSites.length];
                for (int i = 0; i < subSites.length; i++) {
                    newSubSites[i] = original.translateSite(subSites[i]);
                }
                return new FilterAlignment(baseAlignment, newSubSites, original);
            } else if (original.isTaxaFilter()) {
                return new FilterAlignment(baseAlignment, subSites, original);
            } else {
                throw new IllegalStateException("FilterAlignment: getInstance: original not in known state.");
            }
        } else {
            return new FilterAlignment(a, subSites, null);
        }

    }

    public static FilterAlignment getInstance(Alignment a, String[] siteNamesToKeep) {
        
        Arrays.sort(siteNamesToKeep);
        int[] temp = new int[siteNamesToKeep.length];
        int count = 0;
        for (int i = 0, n = a.getSiteCount(); i < n; i++) {
            if (Arrays.binarySearch(siteNamesToKeep, a.getSNPID(i)) >= 0) {
                temp[count++] = i;
                if (count == siteNamesToKeep.length) {
                    break;
                }
            }
        }
        
        int[] result = null;
        if (count == siteNamesToKeep.length) {
            result = temp;
        } else {
            result = new int[count];
            System.arraycopy(temp, 0, result, 0, count);
        }
        return getInstance(a, result);
        
    }
    
    public static FilterAlignment getInstanceRemoveSiteNames(Alignment a, String[] siteNamesToRemove) {

        Arrays.sort(siteNamesToRemove);
        int[] temp = new int[a.getSiteCount()];
        int count = 0;
        for (int i = 0, n = a.getSiteCount(); i < n; i++) {
            if (Arrays.binarySearch(siteNamesToRemove, a.getSNPID(i)) < 0) {
                temp[count++] = i;
            }
        }

        int[] result = null;
        if (count == temp.length) {
            result = temp;
        } else {
            result = new int[count];
            System.arraycopy(temp, 0, result, 0, count);
        }
        return getInstance(a, result);

    }

    public static FilterAlignment getInstance(Alignment a, int startSite, int endSite) {

        if (a instanceof FilterAlignment) {
            FilterAlignment original = (FilterAlignment) a;
            Alignment baseAlignment = ((FilterAlignment) a).getBaseAlignment();
            if (original.isSiteFilter()) {
                int[] subSites = new int[endSite - startSite + 1];
                int[] originalSites = original.getSiteRedirect();
                for (int i = startSite; i <= endSite; i++) {
                    subSites[i - startSite] = originalSites[i];
                }
                return new FilterAlignment(baseAlignment, subSites, original);
            } else if (original.isSiteFilterByRange()) {
                return new FilterAlignment(baseAlignment, original.translateSite(startSite), original.translateSite(endSite), original);
            } else if (original.isTaxaFilter()) {
                return new FilterAlignment(baseAlignment, startSite, endSite, original);
            } else {
                throw new IllegalStateException("FilterAlignment: getInstance: original not in known state.");
            }
        } else {
            return new FilterAlignment(a, startSite, endSite, null);
        }

    }

    public byte getBase(int taxon, int site) {
        int taxaIndex = translateTaxon(taxon);
        if (taxaIndex == -1) {
            return DataType.UNKNOWN_BYTE;
        } else {
            return myBaseAlignment.getBase(taxaIndex, translateSite(site));
        }
    }

    public byte[] getBaseRange(int taxon, int startSite, int endSite) {

        int siteCount = endSite - startSite + 1;
        byte[] result = new byte[siteCount];
        int taxaIndex = translateTaxon(taxon);
        if (taxaIndex == -1) {
            for (int i = 0; i < siteCount; i++) {
                result[i] = DataType.UNKNOWN_BYTE;
            }
            return result;
        } else {
            if (myIsSiteFilterByRange) {
                for (int i = 0; i < siteCount; i++) {
                    result[i] = myBaseAlignment.getBase(taxaIndex, startSite + myRangeStart + i);
                }
                return result;
            } else if (myIsSiteFilter) {
                for (int i = 0; i < siteCount; i++) {
                    result[i] = myBaseAlignment.getBase(taxaIndex, mySiteRedirect[startSite + i]);
                }
                return result;
            } else {
                for (int i = 0; i < siteCount; i++) {
                    result[i] = myBaseAlignment.getBase(taxaIndex, startSite + i);
                }
                return result;
            }
        }

    }

    public byte getBase(int taxon, int site, int allele) {
        int taxaIndex = translateTaxon(taxon);
        if (taxaIndex == -1) {
            return DataType.UNKNOWN_BYTE;
        } else {
            return myBaseAlignment.getBase(taxaIndex, translateSite(site), allele);
        }
    }

    public byte getBase(int taxon, Locus locus, int physicalPosition, int allele) {
        int site = getSiteOfPhysicalPosition(physicalPosition, locus);
        int taxaIndex = translateTaxon(taxon);
        if (taxaIndex == -1) {
            return DataType.UNKNOWN_BYTE;
        } else {
            return myBaseAlignment.getBase(taxaIndex, site, allele);
        }
    }

    public char getBaseChar(int taxon, int site) {
        int taxaIndex = translateTaxon(taxon);
        if (taxaIndex == -1) {
            return DataType.UNKNOWN_CHARACTER;
        } else {
            return myBaseAlignment.getBaseChar(taxaIndex, translateSite(site));
        }
    }

    public int translateSite(int site) {

        if (myIsSiteFilterByRange) {
            return site + myRangeStart;
        } else if (myIsSiteFilter) {
            return mySiteRedirect[site];
        } else {
            return site;
        }

    }

    /**
     * Returns site of this FilterAlignment based on give site from
     * embedded Alignment.
     * 
     * @param site site
     * @return site in this alignment
     */
    public int reverseTranslateSite(int site) {

        if (myIsSiteFilterByRange) {
            return site - myRangeStart;
        } else if (myIsSiteFilter) {
            return Arrays.binarySearch(mySiteRedirect, site);
        } else {
            return site;
        }

    }

    public int translateTaxon(int taxon) {

        if (myIsTaxaFilter) {
            return myTaxaRedirect[taxon];
        } else {
            return taxon;
        }

    }

    private Locus[] getLociFromBase() {

        if ((!myIsSiteFilter) && (!myIsSiteFilterByRange)) {
            return myBaseAlignment.getLoci();
        }

        int numSites = getSiteCount();
        List loci = new ArrayList();
        for (int i = 0; i < numSites; i++) {
            Locus current = getLocus(i);
            if (!loci.contains(current)) {
                loci.add(current);
            }
        }

        Locus[] result = new Locus[loci.size()];
        loci.toArray(result);
        return result;

    }

    public int getIndelSize(int site) {
        return myBaseAlignment.getIndelSize(translateSite(site));
    }

    public Locus[] getLoci() {
        return myLoci;
    }

    public Locus getLocus(int site) {
        return myBaseAlignment.getLocus(translateSite(site));
    }

    public int getLocusSiteCount(Locus locus) {

        int numSites = getSiteCount();
        int result = 0;
        for (int i = 0; i < numSites; i++) {
            if (getLocus(translateSite(i)) == locus) {
                result++;
            }
        }

        return result;

    }

    public int getNumLoci() {
        return myLoci.length;
    }

    public int getPositionInLocus(int site) {
        return myBaseAlignment.getPositionInLocus(translateSite(site));
    }

    public byte getPositionType(int site) {
        return myBaseAlignment.getPositionType(translateSite(site));
    }

    public byte[] getPositionTypes() {

        int numSites = getSiteCount();
        byte[] result = new byte[numSites];
        for (int i = 0; i < numSites; i++) {
            result[i] = myBaseAlignment.getPositionType(translateSite(i));
        }

        return result;

    }

    public float[][] getSiteScores() {

        if (!myBaseAlignment.hasSiteScores()) {
            return null;
        }

        int numSites = getSiteCount();
        int numSeqs = getSequenceCount();
        float[][] result = new float[numSeqs][numSites];
        for (int i = 0; i < numSites; i++) {
            for (int j = 0; j < numSeqs; j++) {
                int taxaIndex = translateTaxon(j);
                if (taxaIndex == -1) {
                    result[j][i] = -9;
                } else {
                    result[j][i] = myBaseAlignment.getSiteScore(taxaIndex, translateSite(i));
                }
            }
        }

        return result;

    }

    public byte getReferenceAllele(int site) {
        return myBaseAlignment.getReferenceAllele(translateSite(site));
    }

    public int getSiteCount() {

        if (myIsSiteFilterByRange) {
            return myRangeEnd - myRangeStart + 1;
        } else if (myIsSiteFilter) {
            return mySiteRedirect.length;
        } else {
            return myBaseAlignment.getSiteCount();
        }

    }

    public String getSNPID(int site) {
        return myBaseAlignment.getSNPID(translateSite(site));
    }

    public String[] getSNPIDs() {

        int numSites = getSiteCount();
        String[] result = new String[numSites];
        for (int i = 0; i < numSites; i++) {
            result[i] = myBaseAlignment.getSNPID(translateSite(i));
        }

        return result;

    }

    public boolean hasReference() {
        return myBaseAlignment.hasReference();
    }

    public boolean isIndel(int site) {
        return myBaseAlignment.isIndel(translateSite(site));
    }

    public int getSiteOfPhysicalPosition(int physicalPosition, Locus locus) {
        int temp = myBaseAlignment.getSiteOfPhysicalPosition(physicalPosition, locus);
        if (temp < 0) {
            temp = -temp;
        }
        return reverseTranslateSite(temp);
    }

    public Alignment getBaseAlignment() {
        return myBaseAlignment;
    }

    public boolean isTaxaFilter() {
        return myIsTaxaFilter;
    }

    public boolean isSiteFilter() {
        return myIsSiteFilter;
    }

    public boolean isSiteFilterByRange() {
        return myIsSiteFilterByRange;
    }

    protected int[] getTaxaRedirect() {
        return myTaxaRedirect;
    }

    protected int[] getSiteRedirect() {
        return mySiteRedirect;
    }

    protected int getRangeStart() {
        return myRangeStart;
    }

    protected int getRangeEnd() {
        return myRangeEnd;
    }

    public byte getMajorAllele(int site) {
        return myBaseAlignment.getMajorAllele(translateSite(site));
    }

    public byte getMinorAllele(int site) {
        return myBaseAlignment.getMinorAllele(translateSite(site));
    }

    public byte[] getMinorAlleles(int site) {
        return myBaseAlignment.getMinorAlleles(translateSite(site));
    }

    public byte[] getAlleles(int site) {
        return myBaseAlignment.getAlleles(translateSite(site));
    }

    public double getMinorAlleleFrequency(int site) {
        return myBaseAlignment.getMinorAlleleFrequency(translateSite(site));
    }

    public int[][] getAllelesSortedByFrequency(int site, boolean includeGAP) {
        return myBaseAlignment.getAllelesSortedByFrequency(translateSite(site), includeGAP);
    }
}
