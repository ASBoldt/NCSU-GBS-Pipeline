/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS;

import java.util.Arrays;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

/**
 *
 * @author ed
 */
public class CompareReadDistribution {
    public static enum Analysis {EvalVirtualDigest, EvalCountFile, CompVD_Count};
    private static Analysis mode=Analysis.EvalVirtualDigest;


    public CompareReadDistribution(Analysis mode, String refGenomeReadFile, String compReadFile) {
    //    theHaplotypeDistMap=new HaplotypeDistMap(inDistFileS, true);

       // mapAllMarkers(outMapFileS);
        switch (mode) {
            case EvalVirtualDigest: {evaluatePhysicalMapFile(refGenomeReadFile,compReadFile);}

        }

    }

    private static void evaluatePhysicalMapFile(String refPhysMapName, String hapFreqFile) {
        ReadsWPhysicalMap refPhysMap=new ReadsWPhysicalMap(refPhysMapName,true);
        refPhysMap.sortTable(true);
        ReadCounts shf=null;
        if(hapFreqFile!=null) {
            shf=new ReadCounts(hapFreqFile,true);
        }
        final int maxDist=5000;
//        ArrayList<int[]> superCopyNumber=new ArrayList<int[]>();  //holds the superrepetive elements
        TreeMap<Integer,Integer> superCopyNumber = new TreeMap<Integer,Integer>();
        int[] allDist=new int[maxDist];
        int[] singleDist=new int[maxDist];
        int[] repDist=new int[maxDist];
        final int maxCopyNum=1000;
        int[] copyNumber=new int[maxCopyNum];

        int[][] mappingProp=null;
        if(shf==null) {
            cycleThroughPhysMap(refPhysMap, allDist, singleDist, repDist, copyNumber, superCopyNumber);
        } else {
            mappingProp=cycleThroughHapFreqList(shf, refPhysMap, allDist, singleDist, repDist, copyNumber, superCopyNumber);
        }
        StringBuilder sb=new StringBuilder(100000);
        if(mappingProp!=null) sb.append(Arrays.deepToString(mappingProp));
        sb.append("CopyNumber\tFrequency\n");
        for (int i = 0; i < maxCopyNum; i++) {if(copyNumber[i]>0) sb.append(i+"\t"+copyNumber[i]+"\n");}
        Set superCopyNumberSet = superCopyNumber.entrySet();
        Iterator SCNsetIter = superCopyNumberSet.iterator();
        while(SCNsetIter.hasNext()) {
            Map.Entry me = (Map.Entry)SCNsetIter.next();
            sb.append(me.getKey()+"\t"+me.getValue()+"\n");
        }
        sb.append("CutSiteDistance\tAllFrequency\tUniqueFrequency\tRepetitiveFrequency\n");
        for (int i = 0; i < maxDist; i++) {sb.append(i+"\t"+allDist[i]+"\t"+singleDist[i]+"\t"+repDist[i]+"\n");}
        System.out.println(sb.toString());
    }

    private static void cycleThroughPhysMap(ReadsWPhysicalMap refPhysMap, int[] allDist, int[] singleDist,
            int[] repDist, int[] copyNumber, TreeMap<Integer,Integer> superCopyNumber) {
        final int maxDist=allDist.length;
        final int maxCopyNum=copyNumber.length;
        int currentCopy=0;
        for (int i = 0; i < refPhysMap.getSize(); i++) {
            int currDist=refPhysMap.getReadWithPosition(i).nextCutDistance;
//            if(refPhysMap.shp[i].readSet[1]==0) System.out.println(BaseEncoder.getSequenceFromLong(refPhysMap.shp[i].readSet));
            if(currDist>=maxDist) currDist=maxDist-1;
            if(currDist<0) currDist=maxDist-2;
            allDist[currDist]++;
            if((i < refPhysMap.getSize()-1)&&(Arrays.equals(refPhysMap.getRead(i), refPhysMap.getRead(i+1))==true)) {
                currentCopy++;
                repDist[currDist]++;
            } else {
                currentCopy++;
                if(currentCopy>maxCopyNum) {
                    if ( superCopyNumber.containsKey(currentCopy) )   {
                        int hits = superCopyNumber.get(currentCopy);
                        superCopyNumber.put(currentCopy,++hits);
                    }
                    else {
                        superCopyNumber.put(currentCopy,1);
                    }
                } else {copyNumber[currentCopy]++;}
                if(currentCopy>1) {repDist[currDist]++;} else {singleDist[currDist]++;}
                currentCopy=0;
            }

        }
    }
    /**
     * This cycles through a readSet frequency list and compares it to the reference genome
     * @param refPhysMap
     * @param allDist
     * @param singleDist
     * @param repDist
     * @param copyNumber
     * @param superCopyNumber
     */
    private static int[][] cycleThroughHapFreqList(ReadCounts shf, ReadsWPhysicalMap refPhysMap, int[] allDist, int[] singleDist,
            int[] repDist, int[] copyNumber, TreeMap<Integer,Integer> superCopyNumber) {
        final int maxDist=allDist.length;
        final int maxCopyNum=copyNumber.length;
        int[][] mappingProp=new int[2][3];  //readSet class in row 0, sum of counts in row 1
        //columns all (0), perfectSingle (1), perfectMultiple (2)
        for (int i = 0; i < shf.getSize(); i++) {
            long[] currHap=shf.getRead(i);
//            if(currHap[1]==0) System.out.println(BaseEncoder.getSequenceFromLong(currHap));
            int[] phyHits=refPhysMap.getReadIndexSet(currHap);
            mappingProp[0][0]++;
            mappingProp[1][0]+=shf.getReadCount(i);
            if(phyHits==null) continue; //no hit
            int avgDist=0;
            for(int p: phyHits) avgDist+=refPhysMap.getReadWithPosition(p).nextCutDistance;
            avgDist/=phyHits.length;  //phyHits.length is the number of hits in the physical genome
            if(avgDist>=maxDist) avgDist=maxDist-1;
            if(avgDist<0) avgDist=maxDist-2;
            allDist[avgDist]+=shf.getReadCount(i);
            if(phyHits.length>maxCopyNum) {
                if ( superCopyNumber.containsKey(phyHits.length) )   {
                    int count = superCopyNumber.get(phyHits.length) + shf.getReadCount(i);
                    superCopyNumber.put(phyHits.length,count);
                }
                else {
                    superCopyNumber.put(phyHits.length,shf.getReadCount(i));
                }
            } else {copyNumber[phyHits.length]+=shf.getReadCount(i);}
            if(phyHits.length==1) {
                mappingProp[0][1]++;
                mappingProp[1][1]+=shf.getReadCount(i);
                singleDist[avgDist]+=shf.getReadCount(i);
            } else {
                mappingProp[0][2]++;
                mappingProp[1][2]+=shf.getReadCount(i);
                repDist[avgDist]+=shf.getReadCount(i);
            }
        }
        return mappingProp;
    }


//    public static void main(String[] args) {
//       if((args.length!=1)) {
//           System.out.println("Require args: mode virtualDigestFile readCountFile");
//           System.exit(1);
//       } else if (args.length==2) {
//           CompareReadDistribution be=new CompareReadDistribution(mode, args[0], args[1]);
//       }
//  }


}
