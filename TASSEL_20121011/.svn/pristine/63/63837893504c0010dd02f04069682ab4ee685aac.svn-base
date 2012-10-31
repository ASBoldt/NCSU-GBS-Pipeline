package net.maizegenetics.pal.alignment;

import net.maizegenetics.pal.datatype.DataType;
import net.maizegenetics.pal.datatype.TextDataType;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Vector;

/**
 * Created by IntelliJ IDEA.
 * User: ed
 * Date: Nov 12, 2006
 * Time: 8:19:09 PM
 * To change this template use File | Settings | File Templates.
 */
public class AnnotatedAlignmentUtils {

    /**
     * Basic constructor.  All annotation is based off the first site in the AnnotationAlignment.
     * This Alignment should not span multiple loci.
     * @param saa the SimpleAnnotatedAlignment to filter
     * @param anchored sets to score anchored indels as same position
     */
    public static Alignment extractIndels(Alignment saa, boolean anchored) {
        Alignment maa = null;
        String[] sequences = new String[saa.getSequenceCount()];
        Vector indel = new Vector();
        StringBuffer[] tempSeq = new StringBuffer[saa.getSequenceCount()];
        findIndels(saa, anchored, tempSeq, indel);
        for (int i = 0; i < sequences.length; i++) {
            sequences[i] = tempSeq[i].toString();
        }
        String[] states = new String[50];
        for (int i = 0; i < 50; i++) {
            states[i] = Integer.toString(i);
        }
        TextDataType dt = new TextDataType(states);
        maa = SimpleAlignment.getInstance(saa.getIdGroup(), sequences, dt);
        IndelPosition ip;
        for (int i = 0; i < maa.getSiteCount(); i++) {
            ip = (IndelPosition) indel.elementAt(i);
            //maa.setChromosome(saa.getChromosome(0), i);
            //maa.setChromosomePosition(saa.getChromosomePosition(0), i);
            //maa.setLocusName(saa.getLocusName(0), i);
            //maa.setLocusPosition(ip.start, i);  //the start of the indel is used for the position
            //maa.setWeightedLocusPosition(i, saa.getWeightedLocusPosition(ip.start));
            //maa.setPositionType(i, saa.getPositionType(ip.start));
        }
        return maa;
    }

    private static void findIndels(Alignment rawAlignment, boolean anchored, StringBuffer[] tempSeq, Vector indel) {
        DataType theRawDataType = rawAlignment.getDataType();
        int rawNumSites = rawAlignment.getSiteCount();
        for (int i = 0; i < rawAlignment.getSequenceCount(); i++) {
            tempSeq[i] = new StringBuffer();
        }
        for (int j = 1; j < rawNumSites - 1; j++) {
            for (int i = 0; i < rawAlignment.getSequenceCount(); i++) {
                if (rawAlignment.getBaseChar(i, j) == DataType.GAP_CHARACTER) {
                    if (rawAlignment.getBaseChar(i, j - 1) != DataType.GAP_CHARACTER) //this is the beginning of a gap
                    {
                        int p = j + 1;
                        while ((rawAlignment.getBaseChar(i, p) == DataType.GAP_CHARACTER) && (p < (rawNumSites - 1))) {
                            p++;
                        }
                        //This tries to prevent indels at the beginning of end of sequence from being scored.
                        if ((theRawDataType.isUnknownChar(rawAlignment.getBaseChar(i, p)) == false) || (rawAlignment.getBaseChar(i, p) == DataType.GAP_CHARACTER)) {
                            IndelPosition currIndel = new IndelPosition(j, p - 1, anchored);
                            if (!indel.contains(currIndel)) {
                                indel.addElement(currIndel);
                                scoreIndelsInAllSequence(currIndel, tempSeq, rawAlignment, anchored);
                            }
                        }
                    }
                }//end c1==gap
            }//end of i
        }//end of j
    }

    private static void scoreIndelsInAllSequence(IndelPosition currIndel, StringBuffer[] tempSeq, Alignment rawAlignment, boolean anchored) {
        int firstchar = 65;
        int j, forwardSize, backwardSize, size;
        int nSites = rawAlignment.getSiteCount() - 1;
        //NumericDataType theNumericDataType = new NumericDataType();
        DataType theRawDataType = rawAlignment.getDataType();
        char nextChar = '?';
        if (anchored) //this finds anchored indels, indels must only share a flanking end
        {
            for (int i = 0; i < rawAlignment.getSequenceCount(); i++) {
                forwardSize = backwardSize = 0;
                nextChar = DataType.UNKNOWN_CHARACTER;
                //System.out.println("tempSeq["+i+"]="+tempSeq[i].toString());
                if ((rawAlignment.getBaseChar(i, currIndel.start - 1) != DataType.GAP_CHARACTER) && (rawAlignment.getBaseChar(i, currIndel.start) == DataType.GAP_CHARACTER)) {
                    j = currIndel.start;
                    while ((rawAlignment.getBaseChar(i, j) == DataType.GAP_CHARACTER) && (j < nSites)) {
                        forwardSize++;
                        j++;
                    }
                    if (theRawDataType.isUnknownChar(rawAlignment.getBaseChar(i, j))) {
                        nextChar = DataType.UNKNOWN_CHARACTER;
                    } else {
                        //nextChar = theNumericDataType.getNumericCharFromNumericIndex(forwardSize);
                        nextChar = (char) (forwardSize + firstchar);
                    }
                } else if ((rawAlignment.getBaseChar(i, currIndel.end + 1) != DataType.GAP_CHARACTER) && (rawAlignment.getBaseChar(i, currIndel.end) == DataType.GAP_CHARACTER)) {
                    j = currIndel.end;
                    while ((rawAlignment.getBaseChar(i, j) == DataType.GAP_CHARACTER) && (j > 0)) {
                        backwardSize++;
                        j--;
                    }
                    if (theRawDataType.isUnknownChar(rawAlignment.getBaseChar(i, j))) {
                        nextChar = DataType.UNKNOWN_CHARACTER;
                    } else {
                        //nextChar = theNumericDataType.getNumericCharFromNumericIndex(backwardSize);
                        nextChar = (char) (backwardSize + firstchar);
                    }
                } else if ((rawAlignment.getBaseChar(i, currIndel.start - 1) == DataType.GAP_CHARACTER) || (rawAlignment.getBaseChar(i, currIndel.end + 1) == DataType.GAP_CHARACTER)
                        || (theRawDataType.isUnknownChar(rawAlignment.getBaseChar(i, currIndel.start))) || (theRawDataType.isUnknownChar(rawAlignment.getBaseChar(i, currIndel.end)))) {
                    nextChar = DataType.UNKNOWN_CHARACTER;
                } else {
                    //nextChar = theNumericDataType.getNumericCharFromNumericIndex(0);
                    nextChar = (char) (firstchar);
                }
                tempSeq[i].append(nextChar);
                //    System.out.println("tempSeq["+i+"]="+tempSeq[i].toString());
            }
        } else //This finding perfectly matched indels
        {
            for (int i = 0; i < rawAlignment.getSequenceCount(); i++) {
                forwardSize = backwardSize = 0;
                //nextChar = theNumericDataType.getNumericCharFromNumericIndex(0);
                nextChar = (char) (firstchar);
                if ((rawAlignment.getBaseChar(i, currIndel.start - 1) == DataType.GAP_CHARACTER) && (rawAlignment.getBaseChar(i, currIndel.end + 1) == DataType.GAP_CHARACTER)) //if the GAP extend beyond both then this is not the same gap
                {
                    nextChar = DataType.UNKNOWN_CHARACTER;
                }
                j = currIndel.start;
                while ((rawAlignment.getBaseChar(i, j) == DataType.GAP_CHARACTER) && (j < nSites)) {
                    forwardSize++;
                    j++;
                }
                if (theRawDataType.isUnknownChar(rawAlignment.getBaseChar(i, j))) {
                    nextChar = DataType.UNKNOWN_CHARACTER;
                }  //if missing data within set to missing
                j = currIndel.end;
                while ((rawAlignment.getBaseChar(i, j) == DataType.GAP_CHARACTER) && (j > 0)) {
                    backwardSize++;
                    j--;
                }
                if (theRawDataType.isUnknownChar(rawAlignment.getBaseChar(i, j))) {
                    nextChar = DataType.UNKNOWN_CHARACTER;
                }
                if (forwardSize == backwardSize) {
                    //nextChar = theNumericDataType.getNumericCharFromNumericIndex(backwardSize);
                    nextChar = (char) (backwardSize + firstchar);
                } else {
                    nextChar = DataType.UNKNOWN_CHARACTER;
                }
                tempSeq[i].append(nextChar);
            }
        }
    }

    /**
     * Remove sites based on site position (excluded sites are <firstSite and >lastSite)
     * This not effect any prior exclusions.
     * @param aa the AnnotatedAlignment to filter
     * @param firstSite first site to keep in the range
     * @param lastSite  last site to keep in the range
     */
    public static Alignment removeSitesOutsideRange(Alignment aa, int firstSite, int lastSite) {
        if ((firstSite < 0) || (firstSite > lastSite)) {
            return null;
        }
        if (lastSite > aa.getSiteCount() - 1) {
            return null;
        }
        return FilterAlignment.getInstance(aa, firstSite, lastSite);
    }

    /**
     * Remove sites based on site types.  Only the included types are kept.
     * This does not effect any prior exclusions.  If types are not set (char=0) they will be included.
     *
     * @param incTypes array of types keep in the alignment
     */
    public static Alignment includeSitesByType(Alignment aa, char[] incTypes) {
        ArrayList<Integer> includeAL = new ArrayList<Integer>();
        boolean keep;
        byte currType;
        for (int i = 0; i < aa.getSiteCount(); i++) {
            keep = false;
            currType = aa.getPositionType(i);
            if (currType == 0) {
                keep = true;
                //  continue;
            }
            for (int j = 0; j < incTypes.length; j++) {
                if (currType == incTypes[j]) {
                    keep = true;
                    break;
                }
            }
            if (keep == true) {
                includeAL.add(i);  //this was set to keep = false, which seems backward Ed
            }
        }
        int[] includeSites = new int[includeAL.size()];
        for (int i = 0; i < includeAL.size(); i++) {
            includeSites[i] = includeAL.get(i);
        }
        Alignment mlaa = FilterAlignment.getInstance(aa, includeSites);
        return mlaa;
    }

    /**
     * remove constant sites but ignore gaps and missing data (- and ?)
     * @param aa the AnnotatedAlignment to filter
     */
    public static Alignment removeConstantSitesIgnoreGapsMissing(Alignment aa) {
        return removeSitesBasedOnFreqIgnoreGapsMissing(aa, 0.000001, 1.0, 2);
    }

    /**
     * remove constant sites, ignoring missing data (N or ?) BUT NOT GAPS (gaps
     * are considered legitimate states)
     * @param aa the AnnotatedAlignment to filter
     */
    public static Alignment removeConstantSitesIgnoreMissing(Alignment aa) {
        return removeSitesBasedOnFreqIgnoreMissing(aa, 0.000001, 1.0, 2);
    }

    public static Alignment removeSitesBasedOnFreqIgnoreGapsMissing(Alignment aa, double minimumProportion, int minimumCount) {
        return removeSitesBasedOnFreqIgnoreGapsMissing(aa, minimumProportion, 1.0, minimumCount);
    }

    /**
     * remove sites based on minimum frequency (the count of good bases)
     * and based on the proportion of good sites different from consensus
     *
     * @param aa the AnnotatedAlignment to filter
     * @param minimumProportion minimum proportion of sites different from the consensus
     * @param minimumCount      minimum number of sequences with a good bases (not - or ?)
     */
    public static Alignment removeSitesBasedOnFreqIgnoreGapsMissing(Alignment aa, double minimumProportion, double maximumProportion, int minimumCount) {
        int[] includeSites = getIncludedSitesBasedOnFreqIgnoreGapsMissing(aa, minimumProportion, maximumProportion, minimumCount);
        Alignment mlaa = FilterAlignment.getInstance(aa, includeSites);
        return mlaa;
    }

    public static Alignment removeSitesBasedOnFreqIgnoreMissing(Alignment aa, double minimumProportion, int minimumCount) {
        return removeSitesBasedOnFreqIgnoreGapsMissing(aa, minimumProportion, 1.0, minimumCount);
    }

    /**
     * remove sites based on minimum frequency (the count of good bases, INCLUDING GAPS)
     * and based on the proportion of good alleles (including gaps) different from consensus
     *
     * @param aa the AnnotatedAlignment to filter
     * @param minimumProportion minimum proportion of sites different from the consensus
     * @param minimumCount      minimum number of sequences with a good bases (not N or ?), where GAP IS CONSIDERED A GOOD BASE
     */
    public static Alignment removeSitesBasedOnFreqIgnoreMissing(Alignment aa, double minimumProportion, double maximumProportion, int minimumCount) {
        int[] includeSites = getIncludedSitesBasedOnFreqIgnoreMissing(aa, minimumProportion, maximumProportion, minimumCount);
        Alignment mlaa = FilterAlignment.getInstance(aa, includeSites);
        return mlaa;
    }

    /**
     * @deprecated use new getIncludedSitesBasedOnFreqIgnoreGapsMissing()
     *
     * get sites to be included based on minimum frequency (the count of good bases)
     * and based on the proportion of good sites different from consensus
     *
     * @param aa the AnnotatedAlignment to filter
     * @param minimumProportion minimum proportion of sites different from the consensus
     * @param minimumCount      minimum number of sequences with a good bases (not - or ?)
     */
    public static int[] getIncludedSitesBasedOnFreqIgnoreGapsMissingDeprecated(Alignment aa, double minimumProportion, int minimumCount) {
        ArrayList<Integer> includeAL = new ArrayList<Integer>();
        char firstGood;
        int countGoodLikeFirst, countGoodNotLikeFirst;
        for (int i = 0; i < aa.getSiteCount(); i++) {
            countGoodLikeFirst = countGoodNotLikeFirst = 0;
            firstGood = DataType.UNKNOWN_CHARACTER;
            for (int j = 0; j < aa.getSequenceCount(); j++) {
                char c = aa.getBaseChar(j, i);
                if ((c != DataType.GAP_CHARACTER) && (c != DataType.UNKNOWN_CHARACTER)) {
                    if (firstGood == DataType.UNKNOWN_CHARACTER) {
                        firstGood = c;
                    }
                    if (c == firstGood) {
                        countGoodLikeFirst++;
                    } else {
                        countGoodNotLikeFirst++;
                    }
                }
            }
            int totalGood = countGoodNotLikeFirst + countGoodLikeFirst;
            double obsMinProp;
            if (countGoodLikeFirst > countGoodNotLikeFirst) {
                obsMinProp = (double) countGoodNotLikeFirst / (double) totalGood;
            } else {
                obsMinProp = (double) countGoodLikeFirst / (double) totalGood;
            }
            if ((totalGood > 0) && (totalGood >= minimumCount) && (obsMinProp >= minimumProportion)) {
                includeAL.add(i);
            }
        }
        int[] includeSites = new int[includeAL.size()];
        for (int i = 0; i < includeAL.size(); i++) {
            includeSites[i] = includeAL.get(i);
        }
        return includeSites;
    }

    public static int[] getIncludedSitesBasedOnFreqIgnoreGapsMissing(Alignment aa, double minimumProportion, int minimumCount) {
        return getIncludedSitesBasedOnFreqIgnoreGapsMissing(aa, minimumProportion, 1.0, minimumCount);
    }

    /**
     * get sites to be included based on minimum frequency (the count of good bases)
     * and based on the proportion of good sites different from consensus
     *
     * @param aa the AnnotatedAlignment to filter
     * @param minimumProportion minimum proportion of sites different from the consensus
     * @param minimumCount      minimum number of sequences with a good bases (not - or ?)
     */
    public static int[] getIncludedSitesBasedOnFreqIgnoreGapsMissing(Alignment aa, double minimumProportion, double maximumProportion, int minimumCount) {

        ArrayList<Integer> includeAL = new ArrayList<Integer>();

        for (int i = 0; i < aa.getSiteCount(); i++) {
            //Map charCount = new HashMap();
            //char[] sortChars = sortAllelesByFrequency(aa, i, charCount);

            int[][] sortAlleles = aa.getAllelesSortedByFrequency(i, false);

            //System.out.println("site: " + i);
            //for (int z = 0; z < sortChars.length; z++) {
            //    System.out.println("char: " + sortChars[z] + "  count: " + charCount.get(sortChars[z]));
            //}
            //System.out.println("");

            int totalNonMissing = 0;
            int numAlleles = 0;
            if (sortAlleles != null) {
                numAlleles = sortAlleles[1].length;
            }
            for (int j = 0; j < numAlleles; j++) {
                totalNonMissing = totalNonMissing + sortAlleles[1][j];
            }
            
            //Iterator itr = charCount.values().iterator();
            //while (itr.hasNext()) {
            //    totalNonMissing = totalNonMissing + ((Integer) itr.next()).intValue();
            //}

            double obsMinProp;
            if ((sortAlleles != null) && (sortAlleles[0].length >= 2)) {
                obsMinProp = (double) sortAlleles[1][1] / (double) totalNonMissing;
            } else {
                obsMinProp = 0.0;
            }
            
            //if (sortChars.length >= 2) {
            //    obsMinProp = ((Integer) charCount.get(sortChars[1])).doubleValue() / totalNonMissing;
            //} else {
            //    obsMinProp = 0.0;
            //}

            if ((totalNonMissing > 0) && (totalNonMissing >= (minimumCount * 2)) && (obsMinProp >= minimumProportion) && (obsMinProp <= maximumProportion)) {
                includeAL.add(i);
            }

        }

        int[] includeSites = new int[includeAL.size()];
        for (int i = 0; i < includeAL.size(); i++) {
            includeSites[i] = includeAL.get(i);
        }
        return includeSites;

    }

    public static int[] getIncludedSitesBasedOnFreqIgnoreMissing(Alignment aa, double minimumProportion, int minimumCount) {
        return getIncludedSitesBasedOnFreqIgnoreMissing(aa, minimumProportion, 1.0, minimumCount);
    }

    /**
     * get sites to be included based on minimum frequency (the count of good
     * bases, INCLUDING GAPS) and based on the proportion of good sites (INCLUDING
     * GAPS) different from consensus
     *
     * @param aa the AnnotatedAlignment to filter
     * @param minimumProportion minimum proportion of sites different from the consensus
     * @param minimumCount      minimum number of sequences with a good base or a gap (but not N or ?)
     */
    public static int[] getIncludedSitesBasedOnFreqIgnoreMissing(Alignment aa, double minimumProportion, double maximumProportion, int minimumCount) {
        ArrayList<Integer> includeAL = new ArrayList<Integer>();
        for (int i = 0; i < aa.getSiteCount(); i++) {
            //Map charCount = new HashMap();
            //char[] sortChars = sortAllelesByFrequencyIncludingGaps(aa, i, charCount);
            
            int[][] sortAlleles = aa.getAllelesSortedByFrequency(i, true);
            
            int totalNonMissing = 0;
            int numAlleles = 0;
            if (sortAlleles != null) {
                numAlleles = sortAlleles[1].length;
            }
            for (int j = 0; j < numAlleles; j++) {
                totalNonMissing = totalNonMissing + sortAlleles[1][j];
            }
            
            //Iterator itr = charCount.values().iterator();
            //while (itr.hasNext()) {
            //    totalNonMissing = totalNonMissing + ((Integer) itr.next()).intValue();
            //}
            
            double obsMinProp;
            if ((sortAlleles != null) && (sortAlleles[0].length >= 2)) {
                obsMinProp = (double) sortAlleles[1][1] / (double) totalNonMissing;
            } else {
                obsMinProp = 0.0;
            }
            
            //if (sortChars.length >= 2) {
            //    obsMinProp = ((Integer) charCount.get(sortChars[1])).doubleValue() / totalNonMissing;
            //} else {
            //    obsMinProp = 0.0;
            //}
            
            if ((totalNonMissing > 0) && (totalNonMissing >= (minimumCount * 2)) && (obsMinProp >= minimumProportion) && (obsMinProp <= maximumProportion)) {
                includeAL.add(i);
            }
        }
        int[] includeSites = new int[includeAL.size()];
        for (int i = 0; i < includeAL.size(); i++) {
            includeSites[i] = includeAL.get(i);
        }
        return includeSites;
    }

    /**
     * Return sorted list of alleles from highest frequency to lowest
     * at given site.
     *
     * @param site site
     * @param charCount map used to store counts
     *
     * @return sorted list of alleles
     */
    /*
    public static char[] sortAllelesByFrequency(Alignment a, int site, Map charCount) {
        if (charCount == null) {
            charCount = new HashMap();
        }
        for (int j = 0; j < a.getSequenceCount(); j++) {
            char current = a.getBaseChar(j, site);
            if ((current != 'N') && (current != '?') && (current != DataType.GAP_CHARACTER)) {
                Integer count = (Integer) charCount.get(current);
                if (count != null) {
                    charCount.put(current, new Integer(count.intValue() + 1));
                } else {
                    charCount.put(current, 1);
                }
            }
        }
        int currentSize = charCount.size();
        char[] sortChars = new char[currentSize];
        Iterator itr = charCount.keySet().iterator();
        int count = 0;
        while (itr.hasNext()) {
            sortChars[count++] = ((Character) itr.next()).charValue();
        }
        boolean change = true;
        while (change) {
            change = false;
            for (int k = 0; k < currentSize - 1; k++) {
                Integer first = (Integer) charCount.get(sortChars[k]);
                Integer second = (Integer) charCount.get(sortChars[k + 1]);
                if (first.compareTo(second) < 0) {
                    Character temp = sortChars[k];
                    sortChars[k] = sortChars[k + 1];
                    sortChars[k + 1] = temp;
                    change = true;
                }
            }
        }
        return sortChars;
    }
     */
    

    /**
     * Return sorted list of alleles (INCLUDING GAPS) from highest frequency
     * to lowest at given site.
     *
     * @param site site
     * @param charCount map used to store counts
     *
     * @return sorted list of alleles (including gaps, if present)
     */
    /*
    public static char[] sortAllelesByFrequencyIncludingGaps(Alignment a, int site, Map charCount) {
        if (charCount == null) {
            charCount = new HashMap();
        }
        for (int j = 0; j < a.getSequenceCount(); j++) {
            char current = a.getBaseChar(j, site);
            if ((current != 'N') && (current != '?')) {
                Integer count = (Integer) charCount.get(current);
                if (count != null) {
                    charCount.put(current, new Integer(count.intValue() + 1));
                } else {
                    charCount.put(current, 1);
                }
            }
        }
        int currentSize = charCount.size();
        char[] sortChars = new char[currentSize];
        Iterator itr = charCount.keySet().iterator();
        int count = 0;
        while (itr.hasNext()) {
            sortChars[count++] = ((Character) itr.next()).charValue();
        }
        boolean change = true;
        while (change) {
            change = false;
            for (int k = 0; k < currentSize - 1; k++) {
                Integer first = (Integer) charCount.get(sortChars[k]);
                Integer second = (Integer) charCount.get(sortChars[k + 1]);
                if (first.compareTo(second) < 0) {
                    Character temp = sortChars[k];
                    sortChars[k] = sortChars[k + 1];
                    sortChars[k + 1] = temp;
                    change = true;
                }
            }
        }
        return sortChars;
    }
     */

    /**
     * This sets all heterozygous sites to missing for imputing purposes
     *
     */
    public static Alignment setHeteroSitesToMissing(Alignment aa) {
        StringBuffer[] tempSeq = new StringBuffer[aa.getSequenceCount()];
        for (int i = 0; i < aa.getSequenceCount(); i++) {
            tempSeq[i] = new StringBuffer();
        }

        for (int i = 0; i < aa.getSiteCount(); i++) {
            for (int j = 0; j < aa.getSequenceCount(); j++) {
                char c = aa.getBaseChar(j, i);
                if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != DataType.GAP_CHARACTER && c != DataType.UNKNOWN_CHARACTER) {
                    tempSeq[j].append(DataType.UNKNOWN_CHARACTER);
                } else {
                    tempSeq[j].append(c);
                }
            }
        }

        String[] sequences = new String[aa.getSequenceCount()];
        for (int i = 0; i < sequences.length; i++) {
            sequences[i] = tempSeq[i].toString();
        }
        Alignment maa = SimpleAlignment.getInstance(aa.getIdGroup(), sequences, aa);
        return maa;
    }

    /**
     * This will sets the third and rarer states at a locus to missing
     *
     */
    public static Alignment setRareMultiAllelicStatesToMissing(Alignment aa) {
        //Vector siteVector = new Vector();
        StringBuffer[] tempSeq = new StringBuffer[aa.getSequenceCount()];
        for (int j = 0; j < aa.getSequenceCount(); j++) {
            tempSeq[j] = new StringBuffer();
        }
        int[] charCount = new int[aa.getSequenceCount()];

        for (int i = 0; i < aa.getSiteCount(); i++) {
            for (int j = 0; j < aa.getSequenceCount(); j++) {
                charCount[j] = 0;
            }
            // count how often each character appears in this column
            for (int j = 0; j < aa.getSequenceCount(); j++) {
                char c = aa.getBaseChar(j, i);
                if ((charCount[j] == 0) && (c != DataType.GAP_CHARACTER) && (c != DataType.UNKNOWN_CHARACTER)) {
                    charCount[j] = 1;
                    //char c = rawAlignment.getBase(j, i);
                    for (int k = j + 1; k < aa.getSequenceCount(); k++) {
                        if (c == aa.getBaseChar(k, i)) {
                            charCount[j]++;
                            charCount[k] = -1;
                        }
                    }
                }
            }

            int maxProportion = 0, secondProportion = 0;
            int maxJ = -1, secondJ = -1;
            for (int j = 0; j < aa.getSequenceCount(); j++) {
                if (charCount[j] > 0) {
                    if (charCount[j] > maxProportion) {
                        secondProportion = maxProportion;
                        secondJ = maxJ;
                        maxProportion = charCount[j];
                        maxJ = j;
                    } else if (charCount[j] > secondProportion) {
                        secondJ = j;
                        secondProportion = charCount[j];
                    }
                }
            }

            // This handles the monomorphic situation.  And when
            // all sites are GAP or UNKNOWN.  GAPs are changed
            // to UNKNOWN in all cases.
            if (maxJ < 0) {
                for (int j = 0; j < aa.getSequenceCount(); j++) {
                    tempSeq[j].append(DataType.UNKNOWN_CHARACTER);
                }
            } else {

                if (secondJ < 0) {
                    secondJ = maxJ;
                }

                char primaryState = aa.getBaseChar(maxJ, i);
                char secondaryState = aa.getBaseChar(secondJ, i);
                for (int j = 0; j < aa.getSequenceCount(); j++) {
                    char c = aa.getBaseChar(j, i);
                    if ((c != primaryState) && (c != secondaryState)) {
                        tempSeq[j].append(DataType.UNKNOWN_CHARACTER);
                    } else {
                        tempSeq[j].append(c);
                    }
                }
            }

        }

        String[] sequences = new String[aa.getSequenceCount()];
        for (int i = 0; i < sequences.length; i++) {
            sequences[i] = tempSeq[i].toString();
        }
        Alignment maa = SimpleAlignment.getInstance(aa.getIdGroup(), sequences, aa);
        return maa;

    }

    /**
     * Basic constructor.  All annotation is based off the first site in the AnnotationAlignment.
     * This Alignment should not span multiple loci.
     * @param aa the AnnotatedAlignment to create haplotypes from.  It should be filtered
     * @param window the numbers of polymorphisms to be combined
     * @param step the numbers of polymorphisms to step when creating sets of polymorphisms (this is the overlap of the haplotypes)
     */
    public static Alignment extractSlidingHaplotypes(Alignment aa, int window, int step) {
        Alignment maa = null, tempMaa = null;
        String[] sequences = new String[aa.getSequenceCount()];
        //TextDataType[] tdts = new TextDataType[(aa.getSiteCount() / step) + 1];
        TextDataType tdts = new TextDataType();
        StringBuffer[] tempSeq = new StringBuffer[aa.getSequenceCount()];
        for (int i = 0; i < sequences.length; i++) {
            tempSeq[i] = new StringBuffer();
        }
        int siteCount = 0;
        DataType dt = aa.getDataType();
        for (int i = 0; i < aa.getSiteCount() - window + 1; i += step) {
            //tdts[siteCount] = new TextDataType();
            for (int t = 0; t < sequences.length; t++) {
                String s = "";
                for (int h = i; h < i + window; h++) {
                    s = s + dt.getFormattedString(aa.getBaseChar(t, h));
                }
                //char base = tdts[siteCount].getCharFromTextRepresentation(s);
                char base = tdts.getCharFromTextRepresentation(s);
                tempSeq[t].append(base);
            }
            siteCount++;
        }
        for (int i = 0; i < sequences.length; i++) {
            sequences[i] = tempSeq[i].toString();
        }
        maa = SimpleAlignment.getInstance(aa.getIdGroup(), sequences, tdts);
        //siteCount = 0;
        //for (int i = 0; i < aa.getSiteCount() - window + 1; i += step) {
        //int max = i + window - 1;
        //maa.setChromosome((aa.getChromosome(i) + aa.getChromosome(max)) / 2, siteCount);
        //maa.setChromosomePosition((aa.getChromosomePosition(0) + aa.getChromosomePosition(0)) / 2f, siteCount);
        //maa.setLocusName(aa.getLocusName(0), siteCount);
        //maa.setLocusPosition((i + max) / 2, siteCount);  //the start of the indel is used for the position
        //maa.setWeightedLocusPosition(siteCount, (aa.getWeightedLocusPosition(i) + aa.getWeightedLocusPosition(max)) / 2);
        //maa.setPositionType(siteCount, aa.getPositionType(i));
        //siteCount++;
        //}
        return maa;
    }

    static class IndelPosition implements java.io.Serializable {

        int start, end;
        boolean anchored;

        IndelPosition(int s, int e, boolean a) {
            start = s;
            end = e;
            anchored = a;
        }

        public boolean equals(Object obj) {
            if (obj instanceof IndelPosition) {
                IndelPosition pt = (IndelPosition) obj;
                if (anchored) {
                    return (start == pt.start) || (end == pt.end);
                } else {
                    return (start == pt.start) && (end == pt.end);
                }
            }
            return super.equals(obj);
        }
    }
}
