/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import net.maizegenetics.genome.BaseEncoder;
import net.maizegenetics.pal.alignment.AbstractAlignment;
import net.maizegenetics.pal.alignment.GdpdmBLOBUtils;
import net.maizegenetics.pal.alignment.Locus;
import net.maizegenetics.pal.datatype.DataType;
import net.maizegenetics.pal.datatype.IUPACNucleotides;
import net.maizegenetics.pal.ids.SimpleIdGroup;

/**
 *
 * @author edbuckler
 */
public class ConvertGBSToAlignment extends AbstractAlignment {
    public static final byte refAllele=GdpdmBLOBUtils.bases[0]; // A
    public static final byte altAllele=GdpdmBLOBUtils.bases[1]; // C
    public static final byte hetAllele=GdpdmBLOBUtils.bases[9]; // M = AC
    public static final byte gAllele=GdpdmBLOBUtils.bases[2];  //  G
    public static final byte tAllele=GdpdmBLOBUtils.bases[3];  //  T
    private byte[][] seq;

    private Locus myLocus;
    private int[] variableSites;
    private byte[] strand;
    private int[] alleleCount;
    private int siteNumber;

    public ConvertGBSToAlignment(ReadsByTaxa rbt, ReadsWPhysicalMap rwpm, double maxDCOP, int minFreq, byte specificChr) {
        super(new SimpleIdGroup(rbt.getTaxaNames()),new IUPACNucleotides());
        rwpm.sortTable(false);  //sort by chromosome order
        myLocus=new Locus(""+specificChr, ""+specificChr, 0,0,null,null);
        prepareArrays(rwpm, maxDCOP, specificChr);
        populateAlleles(rwpm, rbt,maxDCOP,specificChr);
        //populateAllelesWithAllelism(rwpm, rbt,maxDCOP);
    }

    public ConvertGBSToAlignment(ReadsByTaxa rbt, ReadsWPhysicalMap rwpm, double maxRecombAlleleA, double maxRecombAlleleC, byte specificChr) {
        // this is for handling segregation that is NOT 1:1, as in John Doebley's maize-teo BC2S3
        super(new SimpleIdGroup(rbt.getTaxaNames()),new IUPACNucleotides());
        rwpm.sortTable(false);  //sort by chromosome order
        myLocus=new Locus(""+specificChr, ""+specificChr, 0,0,null,null);
//        prepareArraysBC2S3old(rwpm, maxRecomb, specificChr);
//        populateTwoAllelesBC2S3(rwpm, rbt, maxDCOP, specificChr);
        prepareArraysBC2S3(rwpm, maxRecombAlleleA, specificChr);
        populateAllelesBC2S3(rwpm, rbt, maxRecombAlleleA, maxRecombAlleleC, specificChr);
    }

    public ConvertGBSToAlignment(ReadsWPhysicalMap rwpm, ReadsByTaxa rbt, double maxRecombAlleleA, double maxRecombAlleleC, byte specificChr) {
        // this is for handling segregation that is NOT 1:1, as in John Doebley's maize-teo BC2S3
        super(new SimpleIdGroup(rbt.getTaxaNames()),new IUPACNucleotides());
        rwpm.sortTable(false);  //sort by chromosome order
        myLocus=new Locus(""+specificChr, ""+specificChr, 0,0,null,null);
        prepareArraysPhased(rwpm, maxRecombAlleleA, maxRecombAlleleC, specificChr);
        populateAllelesPhased(rwpm, rbt, maxRecombAlleleA, maxRecombAlleleC, specificChr);
    }

    public ConvertGBSToAlignment(ReadsByTaxa rbt, ReadsWPhysicalMap refGenome, ReadsWPhysicalMap rwpm, double maxRecombMinor, double maxRecombMajor, byte specificChr, String AlleleInfoFileName) {
        // this is for providing John Doebley's group with the 64 base teosinte and W22 alleles in the framework maize-teo BC2S3 map
        // assumes that the ref genome is already sorted by position
        super(new SimpleIdGroup(rbt.getTaxaNames()),new IUPACNucleotides());
        rwpm.sortTable(false);  //sort by chromosome order
        myLocus=new Locus(""+specificChr, ""+specificChr, 0,0,null,null);
        writeAlleleInfo(rwpm, rbt, refGenome, maxRecombMinor, maxRecombMajor, specificChr, AlleleInfoFileName);
    }

    public ConvertGBSToAlignment(ReadsByTaxa rbt, ReadsWPhysicalMap rwpm, byte specificChr) {
        // written for the Sorghum stem borer RIL pop, with no framework map, and where neither of the parents are the reference genome
        super(new SimpleIdGroup(rbt.getTaxaNames()),new IUPACNucleotides());
        rwpm.sortTable(false);  //sort by chromosome order
        myLocus=new Locus(""+specificChr, ""+specificChr, 0,0,null,null);
        prepareArrays2alleles(rwpm, rbt, specificChr);
        populateTwoAlleles(rwpm, rbt, specificChr);
    }

    private boolean isInChromosomeSet(byte specificChr, byte currChr) {
        if(specificChr==currChr) return true;
        if((specificChr==99)&&(currChr>0)) return true;
        return false;
    }

    private void prepareArrays(ReadsWPhysicalMap workingMap, double maxDCOP, byte specificChr) {
        ArrayList<ReadWithPositionInfo> posList=new ArrayList<ReadWithPositionInfo>();
        PositionComparator thePC=new PositionComparator();
        ReadWithPositionInfo theLastPos=new ReadWithPositionInfo(new long[2]);
        for (int i=0; i<workingMap.getReadTotal(); i++) {
            ReadWithPositionInfo workReadPos=workingMap.getReadWithPosition(i);
            if(workReadPos.isChromosomeKnown()==false) continue;
            if(workReadPos.getDcoP()>maxDCOP || Float.isNaN(workReadPos.getDcoP())) continue;
            if(!isInChromosomeSet(specificChr,workReadPos.getChromosome())) continue;
            if(thePC.compare(theLastPos, workReadPos)!=0) {
                posList.add(workReadPos);
//                 if(workReadPos.positionMin==theLastPos.positionMin) System.out.println(i+":"+workReadPos.toString());
                 theLastPos=workReadPos;
            }
         }
        siteNumber=posList.size();
        seq=new byte[this.getSequenceCount()][siteNumber];
        variableSites=new int[siteNumber];
        strand=new byte[siteNumber];
        alleleCount=new int[siteNumber];
        for (int i = 0; i < siteNumber; i++) {
            variableSites[i]=posList.get(i).positionMin;
            strand[i]=posList.get(i).strand;
            for(int j=0; j<this.getSequenceCount(); j++) seq[j][i]=DataType.UNKNOWN_BYTE;
        }
        System.out.println("Reads/Alleles:"+workingMap.getReadTotal()+" SitePassing:"+siteNumber);
    }

    private void prepareArraysBC2S3old(ReadsWPhysicalMap workingMap, double maxDCOP, byte specificChr) {
        double bestDCOP;
        ArrayList<ReadWithPositionInfo> alleleList=new ArrayList<ReadWithPositionInfo>();
        ArrayList<ReadWithPositionInfo> posList=new ArrayList<ReadWithPositionInfo>();
        PositionComparator thePC=new PositionComparator();

        ReadWithPositionInfo currentPos=workingMap.getReadWithPosition(0);
        alleleList.add(currentPos);
        bestDCOP = currentPos.dcoP;
        for (int i=1; i<workingMap.getReadTotal(); i++) {
            ReadWithPositionInfo nextPos=workingMap.getReadWithPosition(i);

            // collect alleles with the same position
            while ( currentPos.isChromosomeKnown()
                    && isInChromosomeSet(specificChr,currentPos.getChromosome())
                    && thePC.compare(currentPos, nextPos) == 0
                    && i<workingMap.getReadTotal()  ) {
                alleleList.add(nextPos);
                if (nextPos.dcoP < bestDCOP) {bestDCOP = nextPos.dcoP;}
                ++i;
                if ( i<workingMap.getReadTotal() ) {nextPos=workingMap.getReadWithPosition(i);}
            }
            if (alleleList.size() == 2 && bestDCOP <= maxDCOP) {
                posList.add(alleleList.toArray(new ReadWithPositionInfo[alleleList.size()])[0]);
            }
            alleleList=new ArrayList<ReadWithPositionInfo>();
            currentPos = nextPos;
            alleleList.add(currentPos);
            bestDCOP = currentPos.dcoP;
        }
        siteNumber=posList.size();
        seq=new byte[this.getSequenceCount()][siteNumber];
        variableSites=new int[siteNumber];
        strand=new byte[siteNumber];
        alleleCount=new int[siteNumber];
        for (int i = 0; i < siteNumber; i++) {
            variableSites[i]=posList.get(i).positionMin;
            strand[i]=posList.get(i).strand;
            for(int j=0; j<this.getSequenceCount(); j++) seq[j][i]=DataType.UNKNOWN_BYTE;
        }
        System.out.println("Reads/Alleles:"+workingMap.getReadTotal()+" SitePassing:"+siteNumber);
    }

    private void populateAlleles(ReadsWPhysicalMap workingMap, ReadsByTaxa rbt, double maxDCOP, byte specificChr) {
        int currSite=0;
        for (int i=0; i<workingMap.getReadTotal(); i++) {
            ReadWithPositionInfo workReadPos=workingMap.getReadWithPosition(i);
            if(workReadPos.isChromosomeKnown()==false) continue;
            if(workReadPos.getDcoP()>maxDCOP || Float.isNaN(workReadPos.getDcoP())) continue;
            if(!isInChromosomeSet(specificChr,workReadPos.getChromosome())) continue;
            while((workReadPos.positionMin!=variableSites[currSite])
                ||(workReadPos.strand!=strand[currSite])) {currSite++;}
            alleleCount[currSite]++;
            int currReadIndexInRbt=rbt.getReadIndex(workReadPos.getRead());
            byte currBase=(workReadPos.divergence==0)?refAllele:altAllele;
 //           if(currBase>refAllele) currBase='C';
            for(int j=0; j<this.getSequenceCount(); j++) {
                if(rbt.getReadCountForTaxa(currReadIndexInRbt, j)>0) {
                    if(seq[j][currSite]!=DataType.UNKNOWN_BYTE) {seq[j][currSite]=hetAllele;}
                    else {seq[j][currSite]=currBase;}
                }
            }
        }
    }

    private void populateTwoAllelesBC2S3(ReadsWPhysicalMap workingMap, ReadsByTaxa rbt, double maxDCOP, byte specificChr) {
        int currSite=0;
        double bestDCOP;
        ArrayList<ReadWithPositionInfo> alleleList=new ArrayList<ReadWithPositionInfo>();
        PositionComparator thePC=new PositionComparator();
        ReadWithPositionInfo[] AlleleArray;

        ReadWithPositionInfo currentPos=workingMap.getReadWithPosition(0);
        alleleList.add(currentPos);
        bestDCOP = currentPos.dcoP;
        for (int i=1; i<workingMap.getReadTotal(); i++) {
            ReadWithPositionInfo nextPos=workingMap.getReadWithPosition(i);

            // collect alleles with the same position
            while ( currentPos.isChromosomeKnown()
                    && isInChromosomeSet(specificChr,currentPos.getChromosome())
                    && thePC.compare(currentPos, nextPos) == 0
                    && i<workingMap.getReadTotal()  ) {
                alleleList.add(nextPos);
                if (nextPos.dcoP < bestDCOP) {bestDCOP = nextPos.dcoP;}
                ++i;
                if ( i<workingMap.getReadTotal() ) {nextPos=workingMap.getReadWithPosition(i);}
            }
            if (alleleList.size() == 2 && bestDCOP <= maxDCOP) {
                // sort the alleles by DCOP
                AlleleArray = new ReadWithPositionInfo[alleleList.size()];
                AlleleArray = alleleList.toArray(AlleleArray);
                Arrays.sort(AlleleArray, new DCOPComparator());

                // record the alleles
                for (int allele=0; allele<AlleleArray.length; ++allele) {
                    while((AlleleArray[allele].positionMin!=variableSites[currSite])
                        ||(AlleleArray[allele].strand!=strand[currSite])) {currSite++;}
                    alleleCount[currSite]++;
                    int currReadIndexInRbt=rbt.getReadIndex(AlleleArray[allele].getRead());
                    byte currBase = DataType.UNKNOWN_BYTE;
                    if      (allele==0) {currBase = refAllele;}  // A
                    else if (allele==1) {currBase = altAllele;}  // C
                    for(int j=0; j<this.getSequenceCount(); j++) {
                        if(rbt.getReadCountForTaxa(currReadIndexInRbt, j)>0) {
                            byte recordedBase = seq[j][currSite];
                            if(recordedBase==DataType.UNKNOWN_BYTE) {seq[j][currSite]=currBase;}
                            else {
                                seq[j][currSite] = getIUPAC(recordedBase, currBase);
                            }
                        }
                    }
                }
            }
            alleleList=new ArrayList<ReadWithPositionInfo>();
            currentPos = nextPos;
            alleleList.add(currentPos);
            bestDCOP = currentPos.dcoP;
        }
    }

    private void prepareArraysBC2S3(ReadsWPhysicalMap workingMap, double recombThreshAlleleA, byte specificChr) {
        double bestRecomb;
//        int largestCount;
        ArrayList<ReadWithPositionInfo> alleleList=new ArrayList<ReadWithPositionInfo>();
        ArrayList<ReadWithPositionInfo> posList=new ArrayList<ReadWithPositionInfo>();
        PositionComparator thePC=new PositionComparator();

        ReadWithPositionInfo currentPos=workingMap.getReadWithPosition(0);
        alleleList.add(currentPos);
        bestRecomb = currentPos.dcoP;
//        largestCount = rbt.getReadCount(rbt.getReadIndex(currentPos.getRead()));
        for (int i=1; i<workingMap.getReadTotal(); i++) {
            ReadWithPositionInfo nextPos=workingMap.getReadWithPosition(i);

            // collect alleles with the same position
            while ( currentPos.isChromosomeKnown()
                    && isInChromosomeSet(specificChr,currentPos.getChromosome())
                    && thePC.compare(currentPos, nextPos) == 0
                    && i<workingMap.getReadTotal()  ) {
                alleleList.add(nextPos);
                if (nextPos.dcoP < bestRecomb) {bestRecomb = nextPos.dcoP;}
//                int nextCount = rbt.getReadCount(rbt.getReadIndex(nextPos.getRead()));
//                if (nextCount > largestCount) {largestCount = nextCount;}
                ++i;
                if ( i<workingMap.getReadTotal() ) {nextPos=workingMap.getReadWithPosition(i);}
            }
            if (alleleList.size() > 1 && bestRecomb <= recombThreshAlleleA) {
                posList.add(alleleList.toArray(new ReadWithPositionInfo[alleleList.size()])[0]);
            }
            alleleList=new ArrayList<ReadWithPositionInfo>();
            currentPos = nextPos;
            alleleList.add(currentPos);
            bestRecomb = currentPos.dcoP;
//            largestCount = rbt.getReadCount(rbt.getReadIndex(currentPos.getRead()));
        }
        siteNumber=posList.size();    // there will be more sites in the arrays than are actually used, b/c of additional conditions in the populate step (unused sites will be all missing)
        seq=new byte[this.getSequenceCount()][siteNumber];
        variableSites=new int[siteNumber];
        strand=new byte[siteNumber];
        alleleCount=new int[siteNumber];
        for (int i = 0; i < siteNumber; i++) {
            variableSites[i]=posList.get(i).positionMin;
            strand[i]=posList.get(i).strand;
            for(int j=0; j<this.getSequenceCount(); j++) seq[j][i]=DataType.UNKNOWN_BYTE;
        }
        System.out.println("Reads/Alleles:"+workingMap.getReadTotal()+" SitePassing:"+siteNumber);
    }

    private void populateAllelesBC2S3(ReadsWPhysicalMap workingMap, ReadsByTaxa rbt, double recombThreshAlleleA, double recombThreshAlleleC, byte specificChr) {
        // Uses the allele with the lowest recomb (relative to teo allele [A] in external framework map; n >= 5)
        //   as the teo allele and the allele with the largest number of reads as the W22 allele
        int currSite=0;
        double bestRecomb;
        int largestCount;
        ArrayList<ReadWithPositionInfo> alleleList=new ArrayList<ReadWithPositionInfo>();
        PositionComparator thePC=new PositionComparator();
        ReadWithPositionInfo[] AlleleArray;

        ReadWithPositionInfo currentPos=workingMap.getReadWithPosition(0);
        alleleList.add(currentPos);
        bestRecomb = currentPos.dcoP;
        int currReadIndexInRbt = rbt.getReadIndex(currentPos.getRead());
        largestCount = (currReadIndexInRbt < 0) ? 0 : rbt.getReadCount(currReadIndexInRbt);
        for (int i=1; i<workingMap.getReadTotal(); i++) {
            ReadWithPositionInfo nextPos=workingMap.getReadWithPosition(i);

            // collect alleles with the same position
            while ( currentPos.isChromosomeKnown()
                    && isInChromosomeSet(specificChr,currentPos.getChromosome())
                    && thePC.compare(currentPos, nextPos) == 0
                    && i<workingMap.getReadTotal()  ) {
                alleleList.add(nextPos);
                if (nextPos.dcoP < bestRecomb) {bestRecomb = nextPos.dcoP;}
                int nextReadIndexInRbt = rbt.getReadIndex(nextPos.getRead());
                int nextCount = (nextReadIndexInRbt<0)? 0 : rbt.getReadCount(nextReadIndexInRbt);
                if (nextCount > largestCount) {largestCount = nextCount;}
                ++i;
                if ( i<workingMap.getReadTotal() ) {nextPos=workingMap.getReadWithPosition(i);}
            }
            if (alleleList.size() > 1 && bestRecomb <= recombThreshAlleleA) {
                AlleleArray = new ReadWithPositionInfo[alleleList.size()];
                AlleleArray = alleleList.toArray(AlleleArray);
                Arrays.sort(AlleleArray, new DCOPComparator());  // What will this sort do when DCO = NAN?  Make sure that DCO gets set for all reads.

                // record the alleles
                for (int allele=0; allele<AlleleArray.length; ++allele) {
                    while((variableSites[currSite]!=AlleleArray[allele].positionMin)
                        ||(strand[currSite]!=AlleleArray[allele].strand)) {currSite++;}
                    
                    currReadIndexInRbt=rbt.getReadIndex(AlleleArray[allele].getRead());
                    int currTagReadCount = (currReadIndexInRbt < 0)? 0 : rbt.getReadCount(currReadIndexInRbt);
                    byte currBase = DataType.UNKNOWN_BYTE;
                    if (allele==0 && currTagReadCount <= largestCount) {
                        currBase = refAllele;  // A = teosinte allele
                        alleleCount[currSite]++;
                    }
                    else if (allele>0 && currTagReadCount == largestCount && AlleleArray[allele].mapP <= recombThreshAlleleC) {
                        currBase = altAllele;  // C = W22 allele
                        alleleCount[currSite]++;
                    }
                    if (currBase != DataType.UNKNOWN_BYTE) {
                        for(int j=0; j<this.getSequenceCount(); j++) {
                            if(currReadIndexInRbt>-1 && rbt.getReadCountForTaxa(currReadIndexInRbt, j)>0) {
                                byte recordedBase = seq[j][currSite];
                                if(recordedBase==DataType.UNKNOWN_BYTE) {seq[j][currSite]=currBase;}
                                else {
                                    seq[j][currSite] = getIUPAC(recordedBase, currBase);
                                }
                            }
                        }
                    }
                }
            }
            alleleList=new ArrayList<ReadWithPositionInfo>();
            currentPos = nextPos;
            alleleList.add(currentPos);
            bestRecomb = currentPos.dcoP;
            currReadIndexInRbt = rbt.getReadIndex(currentPos.getRead());
            largestCount = (currReadIndexInRbt < 0) ? 0 : rbt.getReadCount(currReadIndexInRbt);
        }
    }

    private void prepareArraysPhased(ReadsWPhysicalMap workingMap, double recombThreshAlleleA, double recombThreshAlleleC, byte specificChr) {
        double bestRecombA, bestRecombC;
        ArrayList<ReadWithPositionInfo> posList=new ArrayList<ReadWithPositionInfo>();
        PositionComparator thePC=new PositionComparator();

        ArrayList<ReadWithPositionInfo> alleleList=new ArrayList<ReadWithPositionInfo>();
        ReadWithPositionInfo currentPos=workingMap.getReadWithPosition(0);
        alleleList.add(currentPos);
        bestRecombA = currentPos.dcoP;
        bestRecombC = currentPos.mapP;
        for (int i=1; i<workingMap.getReadTotal(); i++) {
            ReadWithPositionInfo nextPos=workingMap.getReadWithPosition(i);

            // collect alleles with the same position
            while ( currentPos.isChromosomeKnown()
                    && isInChromosomeSet(specificChr,currentPos.getChromosome())
                    && thePC.compare(currentPos, nextPos) == 0
                    && i<workingMap.getReadTotal()  ) {
                alleleList.add(nextPos);
                if (nextPos.dcoP < bestRecombA) {bestRecombA = nextPos.dcoP;}
                if (nextPos.mapP < bestRecombC) {bestRecombC = nextPos.mapP;}
                ++i;
                if ( i<workingMap.getReadTotal() ) {nextPos=workingMap.getReadWithPosition(i);}
            }
            if (alleleList.size() > 1 && bestRecombA <= recombThreshAlleleA && bestRecombC <= recombThreshAlleleC) {
                posList.add(alleleList.toArray(new ReadWithPositionInfo[alleleList.size()])[0]);
            }
            alleleList=new ArrayList<ReadWithPositionInfo>();
            currentPos = nextPos;
            alleleList.add(currentPos);
            bestRecombA = currentPos.dcoP;
            bestRecombC = currentPos.mapP;
        }
        siteNumber=posList.size();
        seq=new byte[this.getSequenceCount()][siteNumber];
        variableSites=new int[siteNumber];
        strand=new byte[siteNumber];
        alleleCount=new int[siteNumber];
        for (int i = 0; i < siteNumber; i++) {
            variableSites[i]=posList.get(i).positionMin;
            strand[i]=posList.get(i).strand;
            for(int j=0; j<this.getSequenceCount(); j++) seq[j][i]=DataType.UNKNOWN_BYTE;
        }
        System.out.println("Reads/Alleles:"+workingMap.getReadTotal()+" SitePassing:"+siteNumber);
    }

    private void populateAllelesPhased(ReadsWPhysicalMap workingMap, ReadsByTaxa rbt, double recombThreshAlleleA, double recombThreshAlleleC, byte specificChr) {
        int currSite=0;
        double bestRecombA, bestRecombC;
        PositionComparator thePC=new PositionComparator();
        ReadWithPositionInfo[] AlleleArray;

        ArrayList<ReadWithPositionInfo> alleleList=new ArrayList<ReadWithPositionInfo>();
        ReadWithPositionInfo currentPos=workingMap.getReadWithPosition(0);
        alleleList.add(currentPos);
        bestRecombA = currentPos.dcoP;
        bestRecombC = currentPos.mapP;
        for (int i=1; i<workingMap.getReadTotal(); i++) {
            ReadWithPositionInfo nextPos=workingMap.getReadWithPosition(i);

            // collect alleles with the same position
            while ( currentPos.isChromosomeKnown()
                    && isInChromosomeSet(specificChr,currentPos.getChromosome())
                    && thePC.compare(currentPos, nextPos) == 0
                    && i<workingMap.getReadTotal()  ) {
                alleleList.add(nextPos);
                if (nextPos.dcoP < bestRecombA) {bestRecombA = nextPos.dcoP;}
                if (nextPos.mapP < bestRecombC) {bestRecombC = nextPos.mapP;}
                ++i;
                if ( i<workingMap.getReadTotal() ) {nextPos=workingMap.getReadWithPosition(i);}
            }
            if (alleleList.size() > 1 && bestRecombA <= recombThreshAlleleA && bestRecombC <= recombThreshAlleleC) {
                AlleleArray = new ReadWithPositionInfo[alleleList.size()];
                AlleleArray = alleleList.toArray(AlleleArray);
                Arrays.sort(AlleleArray, new DCOPComparator());  // What will this sort do when DCO = NAN?  Make sure that DCO gets set for all reads.

                // record the alleles
                for (int allele=0; allele<AlleleArray.length; ++allele) {
                    while((variableSites[currSite]!=AlleleArray[allele].positionMin)
                        ||(strand[currSite]!=AlleleArray[allele].strand)) {currSite++;}

                    int currReadIndexInRbt=rbt.getReadIndex(AlleleArray[allele].getRead());
                    byte currBase = DataType.UNKNOWN_BYTE;
                    if (allele==0) {
                        currBase = refAllele;  // A = allele with minimum dcoP
                        alleleCount[currSite]++;
                    }
                    else if (allele>0 && AlleleArray[allele].mapP <= recombThreshAlleleC) {  // what about ties for bestRecombC?
                        currBase = altAllele;  // C = allele with minimum mapP
                        alleleCount[currSite]++;
                    }
                    if (currBase != DataType.UNKNOWN_BYTE) {
                        for(int j=0; j<this.getSequenceCount(); j++) {
                            if(currReadIndexInRbt>-1 && rbt.getReadCountForTaxa(currReadIndexInRbt, j)>0) {
                                byte recordedBase = seq[j][currSite];
                                if(recordedBase==DataType.UNKNOWN_BYTE) {seq[j][currSite]=currBase;}
                                else {
                                    seq[j][currSite] = getIUPAC(recordedBase, currBase);
                                }
                            }
                        }
                    }
                }
            }
            alleleList=new ArrayList<ReadWithPositionInfo>();
            currentPos = nextPos;
            alleleList.add(currentPos);
            bestRecombA = currentPos.dcoP;
            bestRecombC = currentPos.mapP;
        }
    }

    private void prepareArrays2alleles(ReadsWPhysicalMap workingMap, ReadsByTaxa rbt, byte specificChr) {
        ArrayList<ReadWithPositionInfo> alleleList=new ArrayList<ReadWithPositionInfo>();
        ArrayList<ReadWithPositionInfo> posList=new ArrayList<ReadWithPositionInfo>();
        PositionComparator thePC=new PositionComparator();

        ReadWithPositionInfo currentPos=workingMap.getReadWithPosition(0);
        alleleList.add(currentPos);
        for (int i=1; i<workingMap.getReadTotal(); i++) {
            ReadWithPositionInfo nextPos=workingMap.getReadWithPosition(i);

            // collect bi-alleles with the same position
            while ( currentPos.isChromosomeKnown()
                    && isInChromosomeSet(specificChr,currentPos.getChromosome())
                    && thePC.compare(currentPos, nextPos) == 0
                    && i<workingMap.getReadTotal()  ) {
                alleleList.add(nextPos);
                ++i;
                if ( i<workingMap.getReadTotal() ) {nextPos=workingMap.getReadWithPosition(i);}
            }
            if (alleleList.size() == 2) {
                posList.add(alleleList.toArray(new ReadWithPositionInfo[alleleList.size()])[0]);
            }
            alleleList=new ArrayList<ReadWithPositionInfo>();
            currentPos = nextPos;
            alleleList.add(currentPos);
        }
        siteNumber=posList.size();
        seq=new byte[this.getSequenceCount()][siteNumber];
        variableSites=new int[siteNumber];
        strand=new byte[siteNumber];
        alleleCount=new int[siteNumber];
        for (int i = 0; i < siteNumber; i++) {
            variableSites[i]=posList.get(i).positionMin;
            strand[i]=posList.get(i).strand;
            for(int j=0; j<this.getSequenceCount(); j++) seq[j][i]=DataType.UNKNOWN_BYTE;
        }
        System.out.println("Reads/Alleles:"+workingMap.getReadTotal()+" SitePassing:"+siteNumber);
    }

    private void populateTwoAlleles(ReadsWPhysicalMap workingMap, ReadsByTaxa rbt, byte specificChr) {
        int currSite=0;
        ArrayList<ReadWithPositionInfo> alleleList=new ArrayList<ReadWithPositionInfo>();
        PositionComparator thePC=new PositionComparator();
        ReadWithPositionInfo[] AlleleArray;

        ReadWithPositionInfo currentPos=workingMap.getReadWithPosition(0);
        alleleList.add(currentPos);
        for (int i=1; i<workingMap.getReadTotal(); i++) {
            ReadWithPositionInfo nextPos=workingMap.getReadWithPosition(i);

            // collect alleles with the same position
            while ( currentPos.isChromosomeKnown()
                    && isInChromosomeSet(specificChr,currentPos.getChromosome())
                    && thePC.compare(currentPos, nextPos) == 0
                    && i<workingMap.getReadTotal()  ) {
                alleleList.add(nextPos);
                ++i;
                if ( i<workingMap.getReadTotal() ) {nextPos=workingMap.getReadWithPosition(i);}
            }
            if (alleleList.size() == 2) {
                // sort the alleles by divergence from the reference
                AlleleArray = new ReadWithPositionInfo[alleleList.size()];
                AlleleArray = alleleList.toArray(AlleleArray);
                Arrays.sort(AlleleArray, new DivergenceComparator());
                int[] readIndexByAllele = new int[2];
                for (int allele = 0; allele < 2; ++allele) {
                    readIndexByAllele[allele] = rbt.getReadIndex(AlleleArray[allele].getRead());
                }
                if (     readIndexByAllele[0] > -1 && readIndexByAllele[0] < rbt.getReadTotal()
                      && readIndexByAllele[1] > -1 && readIndexByAllele[1] < rbt.getReadTotal() ) {
                    for (int allele=0; allele<AlleleArray.length; ++allele) {  // record the alleles
                        while((AlleleArray[allele].positionMin!=variableSites[currSite])
                            ||(AlleleArray[allele].strand!=strand[currSite])) {currSite++;}
                        alleleCount[currSite]++;
                        byte currBase = DataType.UNKNOWN_BYTE;
                        if      (allele==0) {currBase = refAllele;}  // A
                        else if (allele==1) {currBase = altAllele;}  // C
                        for(int j=0; j<this.getSequenceCount(); j++) {
                            if(rbt.getReadCountForTaxa(rbt.getReadIndex(AlleleArray[allele].getRead()), j)>0) {
                                byte recordedBase = seq[j][currSite];
                                if(recordedBase==DataType.UNKNOWN_BYTE) {seq[j][currSite]=currBase;}
                                else {
                                    seq[j][currSite] = getIUPAC(recordedBase, currBase);
                                }
                            }
                        }
                    }
                }
            }
            alleleList=new ArrayList<ReadWithPositionInfo>();
            currentPos = nextPos;
            alleleList.add(currentPos);
        }
    }

    private void writeAlleleInfo(ReadsWPhysicalMap workingMap, ReadsByTaxa rbt, ReadsWPhysicalMap refGenome, double maxRecombMinor, double maxRecombMajor, byte specificChr, String AlleleInfoFileName) {
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(AlleleInfoFileName), 1000000);
            bw.write("Marker\tAGPv2Chr\tAGPv2StartPos\tAGPv2EndPos\tStrand\tAllele_A\tDivergence_A\tRecomb_A\tAllele_C\tDivergence_C\tRecomb_C\tB73_allele\n");
            double bestRecomb;
            int largestCount;
            PositionComparator thePC=new PositionComparator();
            ReadWithPositionInfo[] AlleleArray;

            ArrayList<ReadWithPositionInfo> alleleList=new ArrayList<ReadWithPositionInfo>();
            ReadWithPositionInfo currentPos=workingMap.getReadWithPosition(0);
            alleleList.add(currentPos);
            bestRecomb = currentPos.dcoP;
            int currReadIndexInRbt = rbt.getReadIndex(currentPos.getRead());
            largestCount = (currReadIndexInRbt < 0) ? 0 : rbt.getReadCount(currReadIndexInRbt);
            for (int i=1; i<workingMap.getReadTotal(); i++) {
                ReadWithPositionInfo nextPos=workingMap.getReadWithPosition(i);

                // collect alleles with the same position
                while ( currentPos.isChromosomeKnown()
                        && isInChromosomeSet(specificChr,currentPos.getChromosome())
                        && thePC.compare(currentPos, nextPos) == 0
                        && i<workingMap.getReadTotal()  ) {
                    alleleList.add(nextPos);
                    if (nextPos.dcoP < bestRecomb) {bestRecomb = nextPos.dcoP;}
                    int nextReadIndexInRbt = rbt.getReadIndex(nextPos.getRead());
                    int nextCount = (nextReadIndexInRbt<0)? 0 : rbt.getReadCount(nextReadIndexInRbt);
                    if (nextCount > largestCount) {largestCount = nextCount;}
                    ++i;
                    if ( i<workingMap.getReadTotal() ) {nextPos=workingMap.getReadWithPosition(i);}
                }
                if (alleleList.size() > 1 && bestRecomb <= maxRecombMinor) {
                    // sort the alleles by DCOP
                    AlleleArray = new ReadWithPositionInfo[alleleList.size()];
                    AlleleArray = alleleList.toArray(AlleleArray);
                    Arrays.sort(AlleleArray, new DCOPComparator());

                    int minAllele = -1, majAllele = -1;
                    for (int allele=0; allele<AlleleArray.length; ++allele) {
                        currReadIndexInRbt=rbt.getReadIndex(AlleleArray[allele].getRead());
                        int currTagReadCount = (currReadIndexInRbt < 0)? 0 : rbt.getReadCount(currReadIndexInRbt);
                        if (allele==0 && currTagReadCount <= largestCount) {
                            minAllele = 0;
                        }
                        else if (allele>0 && currTagReadCount == largestCount && AlleleArray[allele].mapP <= maxRecombMajor) {
                            majAllele = allele;
                        }
                    }

                    if (minAllele == 0 && majAllele > 0) {
                        StringBuilder sb = new StringBuilder(currentPos.getChromosome() + "_" + currentPos.positionMin + "\t");
                        sb.append(currentPos.getChromosome() + "\t");
                        sb.append(currentPos.positionMin + "\t");
                        sb.append(currentPos.positionMax + "\t");
                        sb.append((char)currentPos.strand + "\t");
                        sb.append(BaseEncoder.getSequenceFromLong(AlleleArray[0].getRead()) + "\t");
                        sb.append(AlleleArray[0].divergence + "\t");
                        sb.append(AlleleArray[0].getDcoP() + "\t");
                        sb.append(BaseEncoder.getSequenceFromLong(AlleleArray[majAllele].getRead()) + "\t");
                        sb.append(AlleleArray[majAllele].divergence + "\t");
                        sb.append(AlleleArray[majAllele].getMapP() + "\t");
                        sb.append(BaseEncoder.getSequenceFromLong(refGenome.getRead(
                                     refGenome.getReadIndex(currentPos.getChromosome(),currentPos.strand,currentPos.positionMin)))
                                + "\n");
                        bw.write(sb.toString());
                    }
                }
                alleleList=new ArrayList<ReadWithPositionInfo>();
                currentPos = nextPos;
                alleleList.add(currentPos);
                bestRecomb = currentPos.dcoP;
                currReadIndexInRbt = rbt.getReadIndex(currentPos.getRead());
                largestCount = (currReadIndexInRbt < 0) ? 0 : rbt.getReadCount(currReadIndexInRbt);
            }
            bw.close();
        }
        catch (Exception e) {
            System.out.println("Error writing writeAlleleInfo: "+e);
        }
    }

    private void writeAlleleInfoOld(ReadsWPhysicalMap workingMap, ReadsWPhysicalMap refGenome, double maxDCOP, byte specificChr, String AlleleInfoFileName) {
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(AlleleInfoFileName), 1000000);
            bw.write("Marker\tAGPv2Chr\tAGPv2StartPos\tAGPv2EndPos\tStrand\tAllele_A\tDivergence_A\tDCO_p_A\tAllele_C\tDivergence_C\tDCO_p_C\tB73_allele\n");
            double bestDCOP;
            PositionComparator thePC=new PositionComparator();
            ReadWithPositionInfo[] AlleleArray;

            ArrayList<ReadWithPositionInfo> alleleList=new ArrayList<ReadWithPositionInfo>();
            ReadWithPositionInfo currentPos=workingMap.getReadWithPosition(0);
            alleleList.add(currentPos);
            bestDCOP = currentPos.dcoP;
            for (int i=1; i<workingMap.getReadTotal(); i++) {
                ReadWithPositionInfo nextPos=workingMap.getReadWithPosition(i);

                // collect alleles with the same position
                while ( currentPos.isChromosomeKnown()
                        && isInChromosomeSet(specificChr,currentPos.getChromosome())
                        && thePC.compare(currentPos, nextPos) == 0
                        && i<workingMap.getReadTotal()  ) {
                    alleleList.add(nextPos);
                    if (nextPos.dcoP < bestDCOP) {bestDCOP = nextPos.dcoP;}
                    ++i;
                    if ( i<workingMap.getReadTotal() ) {nextPos=workingMap.getReadWithPosition(i);}
                }
                if (alleleList.size() == 2 && bestDCOP <= maxDCOP) {
                    StringBuilder sb = new StringBuilder(currentPos.getChromosome() + "_" + currentPos.positionMin + "\t");
                    sb.append(currentPos.getChromosome() + "\t");
                    sb.append(currentPos.positionMin + "\t");
                    sb.append(currentPos.positionMax + "\t");
                    sb.append((char)currentPos.strand + "\t");

                    // sort the alleles by DCOP
                    AlleleArray = new ReadWithPositionInfo[alleleList.size()];
                    AlleleArray = alleleList.toArray(AlleleArray);
                    Arrays.sort(AlleleArray, new DCOPComparator());

                    sb.append(BaseEncoder.getSequenceFromLong(AlleleArray[0].getRead()) + "\t");
                    sb.append(AlleleArray[0].divergence + "\t");
                    sb.append(AlleleArray[0].getDcoP() + "\t");
                    sb.append(BaseEncoder.getSequenceFromLong(AlleleArray[1].getRead()) + "\t");
                    sb.append(AlleleArray[1].divergence + "\t");
                    sb.append(AlleleArray[1].getDcoP() + "\t");
                    sb.append(BaseEncoder.getSequenceFromLong(refGenome.getRead(
                                 refGenome.getReadIndex(currentPos.getChromosome(),currentPos.strand,currentPos.positionMin)))
                            + "\n");
                    bw.write(sb.toString());
                }
                alleleList=new ArrayList<ReadWithPositionInfo>();
                currentPos = nextPos;
                alleleList.add(currentPos);
                bestDCOP = currentPos.dcoP;
            }
            bw.close();
        }
        catch (Exception e) {
            System.out.println("Error writing writeAlleleInfo: "+e);
        }
    }

    private void populateAllelesWithAllelism(ReadsWPhysicalMap workingMap, ReadsByTaxa rbt, double maxDCOP) {
        int currSite=0;
        int currPos=-1;
        AlleleStructure as=new AlleleStructure();
        for (int i=0; i<workingMap.getReadTotal(); i++) {
            ReadWithPositionInfo workReadPos=workingMap.getReadWithPosition(i);
            if(workReadPos.isChromosomeKnown()!=true) continue;
            if(workReadPos.getChromosome()<0) continue;
            System.out.println(workReadPos.toString());
            if(workReadPos.positionMin!=currPos) {
                if(as.getNumberOfAlleles()>1) {
                //run analysis
                //add alleles to alignment
                    System.out.println(currPos);
                    System.out.println(as.toString());
                }
                as=new AlleleStructure();
                currPos=workReadPos.positionMin;
            }
            int currReadIndexInRbt=rbt.getReadIndex(workReadPos.getRead());
            as.addAllele(workReadPos.getRead(), rbt.getReadCountsForTaxa(currReadIndexInRbt));
            }
//                    (Arrays.binarySearch(variableSites, workReadPos.positionMin)>=0)) {
//                while(workReadPos.positionMin!=variableSites[currSite]) {currSite++;}
//                alleleCount[currSite]++;
//                int currReadIndexInRbt=rbt.getReadIndex(workReadPos.readSet);
//                byte currBase=GdpdmBLOBUtils.bases[workReadPos.divergence];
//           if(currBase>'A') currBase='C';
//                for(int j=0; j<this.getSequenceCount(); j++) {
//                    if(rbt.getReadCountForTaxa(currReadIndexInRbt, j)>0) {
//                        if(seq[j][currSite]!=DataType.UNKNOWN_BYTE) {seq[j][currSite]=GdpdmBLOBUtils.bases[4];}
//                        else {seq[j][currSite]=currBase;}
//                    }
//                }
//            }
//         }
    }

    private byte getIUPAC(byte allele1, byte allele2) {
        byte A  = GdpdmBLOBUtils.bases[0]; // A
        byte C  = GdpdmBLOBUtils.bases[1]; // C
        byte G  = GdpdmBLOBUtils.bases[2]; // G
        byte T  = GdpdmBLOBUtils.bases[3]; // T
        byte AC = GdpdmBLOBUtils.bases[9]; // M
        byte AG = GdpdmBLOBUtils.bases[4]; // R
        byte AT = GdpdmBLOBUtils.bases[7]; // W
        byte CG = GdpdmBLOBUtils.bases[6]; // S
        byte CT = GdpdmBLOBUtils.bases[5]; // Y
        byte GT = GdpdmBLOBUtils.bases[8]; // K
        byte V =  GdpdmBLOBUtils.bases[13]; // use for unknown geno

        byte temp;
        if (allele1 > allele2) {
            temp = allele1;
            allele1 = allele2;
            allele2 = temp;
        }
        if (allele1 == A && allele2 == A) {return A;}
        if (allele1 == C && allele2 == C) {return C;}
        if (allele1 == G && allele2 == G) {return G;}
        if (allele1 == T && allele2 == T) {return T;}
        if (allele1 == A && allele2 == C) {return AC;}
        if (allele1 == A && allele2 == G) {return AG;}
        if (allele1 == A && allele2 == T) {return AT;}
        if (allele1 == C && allele2 == G) {return CG;}
        if (allele1 == C && allele2 == T) {return CT;}
        if (allele1 == G && allele2 == T) {return GT;}
        return V;
    }



    @Override
    public char getBaseChar(int taxon, int site) {
        return (char)getBase(taxon,site);
    }

    @Override
    public byte getBase(int taxon, int site) {
        return seq[taxon][site];
    }

    @Override
    public byte getBase(int taxon, int site, int allele) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte getBase(int taxon, Locus locus, int physicalPosition, int allele) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public String[] getSNPIDs() {
        int sites = getSiteCount();
        String[] SNPids = new String[sites];
        for (int i = 0; i < sites; i++) {
            SNPids[i] = getSNPID(i);
        }
        return SNPids;
    }

    @Override
    public String getSNPID(int site) {
        return getLocus(site).getChromosomeName()+"_"+variableSites[site];
    }

    @Override
    public int getSiteCount() {
        return siteNumber;
    }

    @Override
    public int getLocusSiteCount(Locus locus) {
        return getSiteCount();
    }

    @Override
    public int getPositionInLocus(int site) {
        return variableSites[site];
    }

    @Override
    public int getSiteOfPhysicalPosition(int physicalPosition, Locus locus) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte getPositionType(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte[] getPositionTypes() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public Locus getLocus(int site) {
        return myLocus;
    }

    @Override
    public Locus[] getLoci() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getNumLoci() {
        return 1;
    }

}
