package net.maizegenetics.genome.GBS;

import cern.colt.GenericSorting;
import java.io.File;
import java.util.Arrays;

/**
 * Reads and collapses ReadCount files from multiple taxa (entire directory of read
 * count files.  It tries to hold as
 * many haplotypes as memory will allow, and then it eliminates the rare ones.
 *
 * User: ed
 * Date: Jan 26, 2008
 * Time: 8:12:12 AM
 * To change this template use File | Settings | File Templates.
 */
public class CombineReadCounts extends ReadCounts {
    private int maxRows=50000000;  // 50000000 works for Jeff on Aztec

    public CombineReadCounts(String inDir, String outFile, int minCount, boolean binary) {
        super();
        init();
        processDirectoryForFusion(new File(inDir), new File(outFile), minCount, binary);
    }
    
    public CombineReadCounts(Reads baseReads, String inDir, String outFile, int minCount, boolean binary) {
        super();
        init();
        addBaseReads(baseReads, minCount);
        baseReads=null;
        System.gc();
        GenericSorting.quickSort(0, haplotype[0].length, comp, swapper);
        processDirectoryForFusion(new File(inDir), new File(outFile), minCount, binary);
    }

    public CombineReadCounts(Reads baseReads1, Reads baseReads2, String inDir, String outFile, int minCount, boolean binary) {
        super();
        init();
        addBaseReads(baseReads1, minCount);
        baseReads1=null;
        System.gc();
        addBaseReads(baseReads2, minCount);
        baseReads2=null;
        System.gc();
        GenericSorting.quickSort(0, haplotype[0].length, comp, swapper);
        processDirectoryForFusion(new File(inDir), new File(outFile), minCount, binary);
    }

    public CombineReadCounts(String inFile1, String inFile2, String outFile, int minCount, boolean binary) {
        super();
        init();
        ReadCounts rc1 = new ReadCounts(inFile1, true);
        addReadCounts(rc1);
        rc1=null;
        System.gc();
        ReadCounts rc2 = new ReadCounts(inFile2, true);
        addReadCounts(rc2);
        rc2=null;
        System.gc();
        GenericSorting.quickSort(0, haplotype[0].length, comp, swapper);
        collapseCounts();
        writeCountFile(new File(outFile), minCount, binary);
    }


    private void init() {
      haplotype=new long[2][maxRows];
      hapcount=new int[maxRows];
        for (int i = 0; i < getSize(); i++) {
            haplotype[0][i]=Long.MAX_VALUE;
            haplotype[1][i]=Long.MAX_VALUE;
        }
    }

    /**
     * This adds a series of reference reads to the list, and by setting them to 
     * the minCount they will always remain
     * @param baseReads
     * @param minCount
     */
    private void addBaseReads(Reads baseReads, int minCount) {
        long[] lastRead=new long[2];
        int duplicates=0;
        for (int i = 0; i < baseReads.getReadTotal(); i++) {
           long[] ls = baseReads.getRead(i);
           if(Arrays.equals(ls, lastRead)) {duplicates++; continue;}
           haplotype[0][currentRows]=ls[0];
           haplotype[1][currentRows]=ls[1];
           hapcount[currentRows]=minCount;
           currentRows++;
           lastRead=ls;
        }
        System.out.println("Reads:"+baseReads.getReadTotal()+"Duplicates Found:"+duplicates);
    }


    /**
     * This adds a series of ReadCounts to the list, preserving the count that each read has.
     * Assumes no duplicates, so Duplicates Found should be 0.
     * Designed to make it easy to combine ReadCounts from multiple flowcells.
     * @param ReadCountsToAdd
     */
    private void addReadCounts(ReadCounts ReadCountsToAdd) {
        long[] lastRead=new long[2];
        int duplicates=0;
        for (int i = 0; i < ReadCountsToAdd.getReadTotal(); i++) {
           long[] ls = ReadCountsToAdd.getRead(i);
           if(Arrays.equals(ls, lastRead)) {duplicates++; continue;}
           haplotype[0][currentRows]=ls[0];
           haplotype[1][currentRows]=ls[1];
           hapcount[currentRows]=ReadCountsToAdd.getReadCount(i);
           currentRows++;
           lastRead=ls;
        }
        System.out.println("Reads:"+ReadCountsToAdd.getReadTotal()+"Duplicates Found:"+duplicates);
    }


    public void processDirectoryForFusion(File inDirectory, File outFile, int minCount, boolean binary) {
        System.out.println("Current Rows:"+currentRows);
        int totalHapsRead=0;
        int fileCnt=0;
        for(File fn: inDirectory.listFiles())  {
            if(fn.getPath().contains("DS_Store")) continue;  //a mac problem
       //    printRows(50);
            fileCnt++;
           ReadCounts currSHF=new ReadCounts(fn.getPath(),true);
           if(currSHF.getSize()+currentRows>maxRows) {
               //eliminate the infrequent rows
               resetRareHaplotypes(2);
               GenericSorting.quickSort(0, haplotype[0].length, comp, swapper);
           }
           addHaplotypes(currSHF);
           GenericSorting.quickSort(0, haplotype[0].length, comp, swapper);
           totalHapsRead+=currSHF.getReadTotal();
           System.out.println("Current Rows: "+currentRows+"   totalHapsRead: "+totalHapsRead+"   from: "+fileCnt+" files");
 //          System.out.println("Rows Sorted");printRows(100);
 //          collapseCounts();
 //          System.out.println("Rows collapsed");printRows(100);
       }
       writeCountFile(outFile, minCount, binary);
    }

   protected void addHaplotypes(ReadCounts currSHF) {
       //the is just binary add not mergesort.  I could try merge sort
       int start=0, end=currentRows;
       for (int i = 0; i < currSHF.getSize(); i++) {
           long[] ls = {currSHF.haplotype[0][i],currSHF.haplotype[1][i]};
           int hit=this.binarySearchOfHaps(ls, start, end);
           if(hit<0) {
               //add a new read the end of the list
               haplotype[0][currentRows]=currSHF.haplotype[0][i];
               haplotype[1][currentRows]=currSHF.haplotype[1][i];
               hapcount[currentRows]=currSHF.hapcount[i];
               currentRows++;
           } else {
               //add to the current list
               hapcount[hit]+=currSHF.hapcount[i];
           }
       }
   }

   private void resetRareHaplotypes(int minCount) {
       for (int i = 0; i < haplotype[0].length; i++) {
           if((hapcount[i]>0)&&(hapcount[i]<minCount)) {
                haplotype[0][i]=Long.MAX_VALUE;
                haplotype[1][i]=Long.MAX_VALUE;
                hapcount[i]=0;
                currentRows--;
           }
       }
       GenericSorting.quickSort(0, haplotype[0].length, comp, swapper);
   }


   public static void main(String[] args) {
       if((args.length<2)||(args[0].equals("-h"))) {
           System.out.println("Require args: inDirectory outFile (-c minCount) (-b txt/binary)");
           System.exit(1);
       } else if (args.length==2) {
           CombineReadCounts be=new CombineReadCounts(args[0], args[1], 2, true);
       }  else if (args.length>2) {
           System.out.println("Need to implement the parser");
           System.exit(1);
       }
  }
}
