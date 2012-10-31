package net.maizegenetics.genome.GBS;

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
import java.util.TreeMap;
import net.maizegenetics.genome.BaseEncoder;

/**
 * Created by IntelliJ IDEA.
 * User: ed
 * Date: Jan 26, 2008
 * Time: 8:12:12 AM
 * To change this template use File | Settings | File Templates.
 */
public class CreateReadsByTaxa {

    private long[][] haplotype;
    private byte[][] hapDist;
    private int taxa=0;
    private String[] taxaNames;
    private TagMatchFinder theTagMatchFinder;

    public CreateReadsByTaxa(String inCountFileString, String inDirString, String outDistFileString, boolean binary) {
        File inCountFile=new File(inCountFileString);
        int rows=(int)(inCountFile.length()/20);  //20 is the number of bytes 2 longs (2*8b) + 1 int (4b)
        haplotype=new long[2][rows];
        readReferenceCountFile(inCountFile, rows);
        theTagMatchFinder=new TagMatchFinder(haplotype);
        processDirectoryForCounting(new File(inDirString));
        this.writeDistFile(new File(outDistFileString), binary);
        System.out.println("");
        for (int i = 0; i < 100; i++) {
          System.out.print(BaseEncoder.getSequenceFromLong(haplotype[0][i])+BaseEncoder.getSequenceFromLong(haplotype[1][i])+"\t");
          for(int t=0; t<taxa; t++)  {
              System.out.print(hapDist[i][t]+"\t");
            }
            System.out.println("");
        }
    }

    // new constructor with minNTaxa added so that only reads with >= minNTaxa are written to the file
    public CreateReadsByTaxa(String inCountFileString, String inDirString, String outDistFileString, int minNTaxa, boolean binary) {
        File inCountFile=new File(inCountFileString);
        int rows=(int)(inCountFile.length()/20);  //20 is the number of bytes 2 longs (2*8b) + 1 int (4b)
        haplotype=new long[2][rows];
        readReferenceCountFile(inCountFile, rows);
        theTagMatchFinder=new TagMatchFinder(haplotype);
        processDirectoryForCounting(new File(inDirString));
        this.writeDistFile(new File(outDistFileString), minNTaxa, binary);
        System.out.println("");
        for (int i = 0; i < 100; i++) {
          System.out.print(BaseEncoder.getSequenceFromLong(haplotype[0][i])+BaseEncoder.getSequenceFromLong(haplotype[1][i])+"\t");
          for(int t=0; t<taxa; t++)  {
              System.out.print(hapDist[i][t]+"\t");
            }
            System.out.println("");
        }
    }


    public void processDirectoryForCounting(File inDirectory) {
        File[] fn=inDirectory.listFiles();
        taxa=fn.length;
      //  taxa=30;
        hapDist=new byte[haplotype[0].length][taxa];
        taxaNames=new String[taxa];
        for(int i=0; i<taxa; i++)  {
            taxaNames[i]=fn[i].getName().split("_")[0];  //taxa name has to be first part of file
            readTaxaCountFile(fn[i], i);
       }
    }

    private void readTaxaCountFile(File inCountFile, int taxon) {
        int totalReads=0, totalReadsMatched=0, totalCntMatched=0;
        int totalReadsUnmatched=0, totalCntUnmatched=0;
        int onlyHitOld=0, onlyHitNew=0;
        try{
        int rows=(int)(inCountFile.length()/20);
        DataInputStream dis = new DataInputStream(new BufferedInputStream (new FileInputStream(inCountFile),65536));
        System.out.println("File = " + inCountFile);
        long[] hap=new long[2];
        int hcnt;
        for(int i=0; i<rows; i++) {
            hap[0]=dis.readLong();
            hap[1]=dis.readLong();
            hcnt=dis.readInt();
            totalReads++;
            int hit=this.getReadIndex(hap);
//            int hit=binarySearchOfHaps(hap);
//            int hit1=binarySearchOfHapsV2(hap,true,true);
//            if((hit>-1)&&(hit1==-1)) onlyHitOld++;
//            if((hit==-1)&&(hit1>-1)) onlyHitNew++;
//            hit=hit1;
            if(hit>-1) {
                totalReadsMatched++;
                totalCntMatched+=hcnt;
                int sum=hapDist[hit][taxon]+hcnt;
                hapDist[hit][taxon]=(sum>127)?127:(byte)sum;
            } else {
                totalReadsUnmatched++;
                totalCntUnmatched+=hcnt;
//                System.out.println("NoMatch:"+BaseEncoder.getSequenceFromLong(hap[0])+BaseEncoder.getSequenceFromLong(hap[1]));
            }
        }
        dis.close();
        }
        catch(Exception e) {
            System.out.println("Error c="+totalReads+" e="+e);
        }
       System.out.print("ReadInFile:"+totalReads+" ReadsMatched:"+totalReadsMatched+" ReadCntMatched:"+totalCntMatched);
       System.out.println(" ReadsUnMatched:"+totalReadsUnmatched+" ReadCntUnMatched:"+totalCntUnmatched);
     //  System.out.println("onlyHitOld:"+onlyHitOld+" onlyHitNew:"+onlyHitNew);
    }

    private void readReferenceCountFile(File currentFile, int hapsToRead) {
        int totalReads=0;
        try{
         DataInputStream dis = new DataInputStream(new BufferedInputStream (new FileInputStream(currentFile),65536));
        System.out.println("File = " + currentFile);
            for(int i=0; i<hapsToRead; i++) {
                haplotype[0][i]=dis.readLong();
                haplotype[1][i]=dis.readLong();
                dis.readInt();
                totalReads++;
//                if(totalReads%100000==0)  System.out.println("Current Reads="+totalReads);
            }
        dis.close();
        }
        catch(Exception e) {
            System.out.println("Error c="+totalReads+" e="+e);
        }
       System.out.println("Count of Tags="+totalReads);
    }

     /**
     * Gets the first index of a read (the only one if a unique list).
     * If the read is not found then it return
     * a negative value indicating its insertion point.
     * @param read as a compressed long array
     * @return the index of the read in the array
     */
  
    public int getReadIndex(long[] read) {
        int hit=Arrays.binarySearch(haplotype[0], read[0]);
        if(hit<1) return hit;
        while((hit > 0) && (haplotype[0][hit-1]==read[0])) {
			hit--;
		}
        while((haplotype[0][hit]==read[0])&&(hit<haplotype[0].length-1)&&(haplotype[1][hit]<read[1])) {hit++;}
        if(((haplotype[0][hit]==read[0])&&(haplotype[1][hit]==read[1]))) return hit;
        return -hit;
    }

  private int binarySearchOfHaps(long[] key) {
      int hit1=Arrays.binarySearch(haplotype[0], key[0]);
      if(hit1<0) return -1;
      int hit2=hit1;
      while((hit2>0)&&(haplotype[0][hit2-1]==key[0])) {
          hit2--;
      }
      while((haplotype[0][hit2]==key[0])&&(hit2<haplotype[0].length)) {
          if(haplotype[1][hit2]==key[1]) {
              return hit2;
          }
          if(BaseEncoder.seqDifferences(haplotype[1][hit2], key[1], 2)<3){
              return hit2;
          }
          hit2++;
      }
      return -1;
  }

  private int binarySearchOfHapsV2(long[] query, boolean onlyOneBest, boolean breakTiesToCommon) {
      TreeMap<Integer,Integer> result=theTagMatchFinder.findMatchesWithIntLengthWords(query, 4, true);
      if((result.size()==0)||(onlyOneBest&&(result.size()>1))) return -1;
      if(result.size()==1) {
          return result.firstEntry().getKey();
      }
      if(breakTiesToCommon) {
        //TODO implement this
      }
      return -1;
  }


   void writeDistFile(File outFile, boolean binary) {
       int hapsOutput=0;
       try{
       DataOutputStream fw=new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outFile),65536));
       if(binary) {
           fw.writeInt(taxa);
           fw.writeInt(haplotype[0].length);
           for(int t=0; t<taxa; t++)  {fw.writeUTF(taxaNames[t]);}
       }
       else {
           fw.writeBytes(taxa+"\t"+haplotype[0].length+"\n");
           for(int t=0; t<taxa; t++)  {fw.writeBytes("\t"+taxaNames[t]);}
       }
       if(!binary) {fw.writeBytes("\n");}
       for(int i=0; i<haplotype[0].length; i++) {
          if(!binary) {
             fw.writeBytes(
             BaseEncoder.getSequenceFromLong(haplotype[0][i])+
             BaseEncoder.getSequenceFromLong(haplotype[1][i])+"\t");
             for(int t=0; t<taxa; t++)  {fw.writeBytes(hapDist[i][t]+"\t");}
             fw.writeBytes("\n");
            }
         else {
            fw.writeLong(haplotype[0][i]);
            fw.writeLong(haplotype[1][i]);
            for(int t=0; t<taxa; t++)  {fw.writeByte(hapDist[i][t]);}
            //fw.writeInt(c);
         }
         hapsOutput++;
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

    void writeDistFile(File outFile, int minNTaxa, boolean binary) {
        int hapsOutput=0;
        int nReadsWEnoughTaxa = 0;
       try{
        boolean[] enoughTaxa = new boolean[haplotype[0].length];  // does the indexed haplotype have enough taxa to output?
        for(int i=0; i<haplotype[0].length; i++) {
            int nTaxaWData = 0;
            for(int t=0; t<taxa; t++)  {
               if (hapDist[i][t] > 0) ++nTaxaWData;
            }
            enoughTaxa[i] = nTaxaWData >= minNTaxa ? true : false;
            if (enoughTaxa[i]) ++nReadsWEnoughTaxa;
        }
        DataOutputStream fw=new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outFile),65536));
        if(binary) {
           fw.writeInt(taxa);
           fw.writeInt(nReadsWEnoughTaxa);
           for(int t=0; t<taxa; t++)  {fw.writeUTF(taxaNames[t]);}
        }
        else {
           fw.writeBytes(taxa+"\t"+nReadsWEnoughTaxa+"\n");
           for(int t=0; t<taxa; t++)  {fw.writeBytes("\t"+taxaNames[t]);}
        }
        if(!binary) {fw.writeBytes("\n");}
        for(int i=0; i<haplotype[0].length; i++) {
            if (enoughTaxa[i]) {
                if(!binary) {
                    fw.writeBytes(
                    BaseEncoder.getSequenceFromLong(haplotype[0][i])+
                    BaseEncoder.getSequenceFromLong(haplotype[1][i])+"\t");
                    for(int t=0; t<taxa; t++)  {fw.writeBytes(hapDist[i][t]+"\t");}
                    fw.writeBytes("\n");
                }
                else {
                    fw.writeLong(haplotype[0][i]);
                    fw.writeLong(haplotype[1][i]);
                    for(int t=0; t<taxa; t++)  {fw.writeByte(hapDist[i][t]);}
                }
                hapsOutput++;
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

 Swapper swapper = new Swapper() {
   public void swap(int a, int b) {
      long t1, t2;
      t1 = haplotype[0][a]; haplotype[0][a] = haplotype[0][b];	haplotype[0][b] = t1;
      t2 = haplotype[1][a]; haplotype[1][a] = haplotype[1][b]; haplotype[1][b] = t2;
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
           System.out.println("Require args: inRefCountFile inCountDirectory outDistFile  (-b binaryT/F)");
           System.exit(1);
       } else if (args.length==3) {
           CreateReadsByTaxa be=new CreateReadsByTaxa(args[0],args[1], args[2], true);
       }  else if (args.length>3) {
           System.out.println("Need to implement the parser");
           System.exit(1);
       }
  }
}
