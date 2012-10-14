package net.maizegenetics.genome.GBS;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.util.Arrays;
import net.maizegenetics.genome.BaseEncoder;

/**
 * Holds haplotypes data compressed in long and a counts.  Has basic filtering methods.
 *
 * User: ed
 * Date: Jan 26, 2008
 * Time: 8:12:12 AM
 * To change this template use File | Settings | File Templates.
 */
public class ReadCounts implements Reads {

    long[][] haplotype;  //for memory efficiency the rows first and second half of the read
    //columns are the index of the reads.
    int[] hapcount;
    int currentRows;

    public ReadCounts() {
        currentRows = 0;
    }

    public ReadCounts(String inFile, String outFile, int minCount, boolean binary, boolean simpleFilter, boolean combineIdenticalTaxa) {
        currentRows = 0;
        if(simpleFilter) {
            simpleFilter(new File( inFile), new File(outFile), minCount, binary);
        } 
        else if (combineIdenticalTaxa) {
            combineIdenticalTaxa(new File(inFile), binary);
            processDirectoryForCounting(new File(inFile), new File(outFile), minCount, binary);
        }
        else {
            processDirectoryForCounting(new File(inFile), new File(outFile), minCount, binary);
        }
    }

    public ReadCounts(String inFile, String rowOut) {
        currentRows = 0;
        File fn=new File(inFile);
        int rows=(int)(fn.length()/20);  //20 is the number of bytes 2 longs (2*8b) + 1 int (4b)
        haplotype=new long[2][rows];
        hapcount=new int[rows];
        readSolexaFastaToCountList(new File(inFile),rows);
        this.printRows(Integer.parseInt(rowOut));
    }

    public ReadCounts(String inFile, boolean binary) {
        currentRows = 0;
        File fn=new File(inFile);
        int rows=(int)(fn.length()/20);  //20 is the number of bytes 2 longs (2*8b) + 1 int (4b)
        haplotype=new long[2][rows];
        hapcount=new int[rows];
        readSolexaFastaToCountList(new File(inFile),rows);
    }

    public ReadCounts(long[][] haps, int[] hapcnt, int maxCount) {
        currentRows = 0;
        if((haps.length!=2)||(haps[0].length!=hapcnt.length)) return;
        if((maxCount>haps[0].length)||(maxCount<0)) maxCount=haps[0].length;
        this.haplotype=new long[haps.length][maxCount];
        this.hapcount=new int[maxCount];
        for (int i = 0; i < maxCount; i++) {
            haplotype[0][i]=haps[0][i];
            haplotype[1][i]=haps[1][i];
            hapcount[i] = hapcnt[i];
            ++currentRows;
        }
        printRows(100);
        GenericSorting.quickSort(0, haplotype[0].length, comp, swapper);
        printRows(100);
    }

    public int getReadCount(int index) {
        if(index>=hapcount.length) return -1;
        return hapcount[index];
    }

    public long[] getRead(int index) {
        if(index>=haplotype[0].length) return null;
        long[] result={haplotype[0][index],haplotype[1][index]};
        return result;
    }

    @Override
    public int getReadIndex(long[] read) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    /**
     * This is the number of different reads in the list (NOT THE SUM OF THE COUNTS)
     * The index will vary from 0 to (ReadTotal-1)
     * This is the number of distinct reads if readUnique is true
     * @return total number of reads
     */
    public int getReadTotal() {
        return haplotype[0].length;
    }

    private void processDirectoryForCounting(File inDirectory, File outDirectory, int minCount, boolean binary) {
        for(File fn: inDirectory.listFiles())  {
           int rows=(int)(fn.length()/20);  //20 is the number of bytes 2 longs (2*8b) + 1 int (4b)
           currentRows = 0;
           haplotype=new long[2][rows];
           hapcount=new int[rows];
           readSolexaFastaToCountList(fn, rows);
 //           System.out.println("Rows Read");printRows(100);
           GenericSorting.quickSort(0, haplotype[0].length, comp, swapper);
 //          System.out.println("Rows Sorted");printRows(100);
           collapseCounts();
 //          System.out.println("Rows collapsed");printRows(100);
           File outFile=new File(outDirectory.getPath()+"/"+fn.getName());
           writeCountFile(outFile, minCount, binary);
       }
    }

    private void combineIdenticalTaxa(File inDirectory, boolean binary) {
        File[] FileList = inDirectory.listFiles();
        Arrays.sort(FileList);
        System.out.println("\n\nSorted readCount files in input directory:");
        for(File fn: FileList)  {
           System.out.println(fn.getName());
        }
        System.out.println("\n");
        File prevFile = null;
        String prevTaxon = "";
        String prevFlowcell = "";
        for(File currFile: FileList)  {
            String[] fileParts = currFile.getName().split("_");
            String currTaxon = fileParts[0];
            String currFlowcell = fileParts[1];
            fileParts = currFile.getName().split("\\.");
            String ext = fileParts[fileParts.length-1];
            if (currTaxon.equals(prevTaxon) && ext.equals("cnt")) {  // the second condition was added to prevent unintentional deletion of files
                System.out.println(prevFile.getName() + " and " + currFile.getName() + " are being merged");
                boolean sameFlowcell;
                if (currFlowcell.equals(prevFlowcell)) {
                    sameFlowcell = true;
                }
                else {
                    sameFlowcell = false;
                    currFlowcell = "merge";
                }
                currFile = mergeIdenticalTaxa(prevFile, currFile, currTaxon, sameFlowcell, binary);
            }
            prevFile = currFile;
            prevTaxon = currTaxon;
            prevFlowcell = currFlowcell;
        }
    }

    private File mergeIdenticalTaxa(File prevFile, File currFile, String Taxon, boolean sameFlowcell, boolean binary) {
        String mergedFileName;
        if (sameFlowcell) {
           mergedFileName = currFile.getParent() + "/" + Taxon + "_" + currFile.getName().split("_")[1] + "_merged.cnt";
        }
        else {
           mergedFileName = currFile.getParent() + "/" + Taxon + "_merged.cnt";
        }
        ReadCounts mergedRC = new ReadCounts();
        ReadCounts rc1 = new ReadCounts(prevFile.getPath(), true);
        ReadCounts rc2 = new ReadCounts(currFile.getPath(), true);
        int rows = rc1.currentRows + rc2.currentRows;
        mergedRC.haplotype=new long[2][rows];
        mergedRC.hapcount=new int[rows];
        mergedRC.addReadCounts(rc1);
        rc1=null;
        prevFile.delete();
        System.gc();
        mergedRC.addReadCounts(rc2);
        rc2=null;
        currFile.delete();
        System.gc();
        File mergedFile = new File(mergedFileName);
        mergedRC.writeCountFile(mergedFile, 0, binary);
        return mergedFile;
    }


    /**
     * This adds a series of ReadCounts to the list, preserving the count that each read has.
     * Designed to make it easy to combine ReadCounts from multiple barcodes, lanes or flowcells.
     * @param ReadCountsToAdd
     */
   private void addReadCounts(ReadCounts ReadCountsToAdd) {
        for (int i = 0; i < ReadCountsToAdd.getReadTotal(); i++) {
           long[] ls = ReadCountsToAdd.getRead(i);
           this.haplotype[0][currentRows]=ls[0];
           this.haplotype[1][currentRows]=ls[1];
           this.hapcount[currentRows]=ReadCountsToAdd.getReadCount(i);
           currentRows++;
        }
    }

    protected void printRows(int numRows) {
        for(int i=0; i<numRows; i++) {
            System.out.println(BaseEncoder.getSequenceFromLong(haplotype[0][i])+
                    BaseEncoder.getSequenceFromLong(haplotype[1][i])+
                    " "+
                   hapcount[i]
                    );
        }
    }

    protected void readSolexaFastaToCountList(File currentFile, int hapsToRead) {
        int totalReads=0;
        try{
         DataInputStream dis = new DataInputStream(new BufferedInputStream (new FileInputStream(currentFile),65536));
        System.out.println("File = " + currentFile);
            String temp, sl;
            for(int i=0; i<hapsToRead; i++) {
                haplotype[0][i]=dis.readLong();
                haplotype[1][i]=dis.readLong();
                hapcount[i]=dis.readInt();
                totalReads++;
                ++currentRows;
            }
        dis.close();
        }
        catch(Exception e) {
            System.out.println("Error c="+totalReads+" e="+e);
        }
       System.out.println("Count of Tags="+totalReads);
    }

   protected void collapseCounts() {
       //requires that the reads are sorted
       int collapsedRows=0;
       for(int i=1; i<currentRows; i++) {  // formerly: i<haplotype[0].length  (changed by Jeff)
           if((haplotype[0][i-1]==haplotype[0][i])&&(haplotype[1][i-1]==haplotype[1][i])) {
               hapcount[i]+=hapcount[i-1];
               hapcount[i-1]=0;
               collapsedRows++;
           }
       }
       System.out.println("Rows collapsed:"+collapsedRows);
   }

   protected void writeCountFile(File outFile, int minCount, boolean binary) {
      // int c=1;
       int hapsOutput=0;
       try{
       DataOutputStream fw=new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outFile),4000000));       
       for(int i=0; i<haplotype[0].length; i++) {
//         if((haplotype[0][i]==haplotype[0][i+1])&&(haplotype[1][i]==haplotype[1][i+1])) {
//             c++;
//         }
//         else {
         if(hapcount[i]>=minCount) {
             if(!binary) {fw.writeBytes(
             BaseEncoder.getSequenceFromLong(haplotype[0][i])+
             BaseEncoder.getSequenceFromLong(haplotype[1][i])+" "+hapcount[i]+"\n");}
             else {fw.writeLong(haplotype[0][i]);
                fw.writeLong(haplotype[1][i]);
                fw.writeInt(hapcount[i]);
             }
            hapsOutput++;
//          c=1;
         }
       }
       fw.flush();
       fw.close();
       System.out.println("Haplotypes written to:"+outFile.toString());
       System.out.println("Number of Haplotypes in file:"+hapsOutput);
       }
       catch(Exception e) {
             System.out.println("Catch in writing output file e="+e);
       }
   }

    private void simpleFilter(File infile, File outfile, int minCount, boolean binary) {
       int hapsOutput=0;
       try{
       DataInputStream dis = new DataInputStream(new BufferedInputStream (new FileInputStream(infile),65536));
       DataOutputStream dos=new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfile),65536));
       System.out.println(infile+" being subsetted to:"+ outfile);
       int rows=(int)(infile.length()/20);  //20 is the number of bytes 2 longs (2*8b) + 1 int (4b)
       System.out.println("Number of Haplotypes in file:"+rows);
       long[] hap=new long[2];
       int hcnt;
       for(int i=0; i<rows; i++) {
         hap[0]=dis.readLong();
         hap[1]=dis.readLong();
         hcnt=dis.readInt();
         if(hcnt>=minCount) {
             if(!binary) {dos.writeBytes(
             BaseEncoder.getSequenceFromLong(hap[0])+
             BaseEncoder.getSequenceFromLong(hap[1])+" "+hcnt+"\n");}
             else {dos.writeLong(hap[0]);
                dos.writeLong(hap[1]);
                dos.writeInt(hcnt);
             }
            hapsOutput++;
         }
       }
       dos.flush();
       dos.close();
       dis.close();
       System.out.println("Number of Haplotypes written to outfile:"+hapsOutput);
       }
       catch(Exception e) {
             System.out.println("Catch in writing output file e="+e);
       }
   }

    public int getSize() {
        return haplotype[0].length;
    }

    public int getTotalCount() {
        int totalCount=0;
        for (int i = 0; i < hapcount.length; i++) {
            totalCount+=hapcount[i];
        }
        return totalCount;
    }


   public int binarySearchOfHaps(long[] key, int start, int end) {
      int hit1=Arrays.binarySearch(haplotype[0], start, end, key[0]);
      if(hit1<0) return -1;
      int hit2=hit1;
      while((hit2>0)&&(haplotype[0][hit2-1]==key[0])) {
          hit2--;
      }
      while((haplotype[0][hit2]==key[0])&&(hit2<haplotype[0].length)) {
          if(haplotype[1][hit2]==key[1]) {
              return hit2;
          }
          hit2++;
      }
      return -1;
  }


 Swapper swapper = new Swapper() {
   public void swap(int a, int b) {
      long t1, t2;
      t1 = haplotype[0][a]; haplotype[0][a] = haplotype[0][b];	haplotype[0][b] = t1;
      t2 = haplotype[1][a]; haplotype[1][a] = haplotype[1][b]; haplotype[1][b] = t2;
      int t3;
      t3=hapcount[a]; hapcount[a]=hapcount[b]; hapcount[b]=t3;
   }
};

IntComparator comp = new IntComparator() {
   public int compare(int a, int b) {
      if (haplotype[0][a]==haplotype[0][b]) return haplotype[1][a]==haplotype[1][b] ? 0 : (haplotype[1][a]<haplotype[1][b] ? -1 : 1);
      return haplotype[0][a]<haplotype[0][b] ? -1 : 1;
   }
};

   public static void main(String[] args) {
       if((args.length<1)||(args[0].equals("-h"))) {
           System.out.println("Require args: inDirectory (outDirectory) (-c minCount) (-b txt/binary)");
           System.exit(1);
       } else if (args[0].equals("-f")) {
           ReadCounts be=new ReadCounts(args[1],args[2], 10, true, true, false);
       } else if (args[0].equals("-p")) {
           ReadCounts be=new ReadCounts(args[1],args[2]);
       } else if (args.length==1) {
           ReadCounts be=new ReadCounts(args[0],args[0], 1, true, false, false);
       } else if (args.length==2) {
           ReadCounts be=new ReadCounts(args[0], args[1], 1, true, false, false);
       }  else if (args.length>2) {
           System.out.println("Need to implement the parser");
           System.exit(1);
       }
  }

    @Override
    public int[] getReadIndexSet(long[] read) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean areReadsUnique() {
        throw new UnsupportedOperationException("Not supported yet.");
    }


}
