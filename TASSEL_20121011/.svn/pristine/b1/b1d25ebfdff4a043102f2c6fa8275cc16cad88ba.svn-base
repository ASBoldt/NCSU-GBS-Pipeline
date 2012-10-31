package net.maizegenetics.genome;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;

/**
 * Created by IntelliJ IDEA.
 * User: ed
 * Date: Jan 26, 2008
 * Time: 8:12:12 AM
 * To change this template use File | Settings | File Templates.
 */
public class SolexaFreqCreatorV2 {
    //todo count the fraction of MAGIs hit
    static int chunkSize=net.maizegenetics.genome.BaseEncoder.chunkSize;
    String baseDir="E:/SolexaAnal/";
 //   String baseDir="C:/EdStuff/Solexa/";
    //Dir of fastq, change to qseq
 //   File parentReadDir=new File(baseDir+"NGG282/42AMUAAXX-reanalyze/fastq");
//    File parentReadDir=new File(baseDir+"NGG_IBM/434GFAAXX/qseq");
    File parentReadDir=new File(baseDir+"NGG_IBM/434LFAAXX/qseq");
    File parentOutDir=new File(baseDir+"test");

    //tagdist is haplotype by taxa counts
 //    File tagDist=new File(baseDir+"test/42AMUAAXX091021Dist.txt");

    //tagcount is haplotype and freq over all samples
    File tagCount=new File(baseDir+"test/count091221d.txt");

    long[][] haplotype;

    static String nullS="AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    static String[] fullCutSites={"GCAGC","GCTGC"};
    static String[] overhang={"CAGC","CTGC"};
    static String[] cutSitesWithUniveral={"GCAGGATC","GCTGGATC"};
    static String[] allCutSites={fullCutSites[0],fullCutSites[1],cutSitesWithUniveral[0],cutSitesWithUniveral[1]};

    static int maxBarcodeLength=8, minBarcodeLength=4;

    public SolexaFreqCreatorV2() {
       //long[] allReads=readSolexaFastaToCountList();
        processDirectoryForCounting(parentReadDir,parentOutDir,".0cnt");
      // readTagCountToFASTA(tagCount,tagCountFasta);

    }

  public static void readTagCountToFASTA(File inputFile, File outputFile) {
    int c=0;
    try{
       FileReader fr=new FileReader(inputFile);
       BufferedReader br=new BufferedReader(fr,5000000);
       FileWriter fw=new FileWriter(outputFile);
       BufferedWriter bw=new BufferedWriter(fw,5000000);
       String[] s;
       String temp;
       while(br.ready()) {
            temp=br.readLine();
            s=temp.split(" ");
            if(s.length>0)
                {bw.write(">H"+c+"\n");
                 bw.write(s[0]+"\n");
                 c++;
            }
            if(c%100000<1)  System.out.println("read number = " + c);
            }
        fr.close();
        br.close();
        System.out.println("read number = " + c);
    }
    catch(Exception e) {
        System.out.println("Catch c="+c+" e="+e);
    }
  }

  public static void subSampleFastQToFASTA(File inputFile, File outputFile, int subsampleRate) {
    int c=0;
    try{
       FileReader fr=new FileReader(inputFile);
       BufferedReader br=new BufferedReader(fr,5000000);
       FileWriter fw=new FileWriter(outputFile);
       BufferedWriter bw=new BufferedWriter(fw,1000000);
       String s;
       String temp;
       while(br.ready()) {
            temp=br.readLine();
           if(temp.startsWith("@")) {                    //fastq starts with @ rather than >
                    s=br.readLine();
                    s=s.substring(8, 8+64);
            c++;
            if((c%subsampleRate==0))
                {bw.write(">T"+c+"\n");
                 bw.write(s+"\n");
                 System.out.println("write read number = " + c);
            }
            if(c%100000<1)  System.out.println("read number = " + c);
            if(c%1000000<1)  bw.flush();
            }
       }
        fr.close();
        br.close();
        System.out.println("total read number = " + c);
    }
    catch(Exception e) {
        System.out.println("Catch c="+c+" e="+e);
    }
  }


    public static String removeSeqAfterSecondCutSite(String seq) {
         //this looks for a second restriction site, and then turns the remaining sequence to AAAA
        int pos=9999;
        int startSearchBase=3;
        for(String fcs: fullCutSites) {
            int p=seq.indexOf(fcs, startSearchBase);
            if((p>startSearchBase)&&(p<pos)) {
                pos=p;
                int slen=seq.length();
                seq=seq.substring(0,pos+4)+nullS.substring(0, slen-pos-4);//this makes all bases after the second site AAAAAAA
            }
        }
        return seq;
    }

    public static String[] findBarcodeHaplotype(String seq, int seqLength) {
        for(String fcs: overhang) {
            int p=seq.indexOf(fcs,minBarcodeLength-1);
            if((p>=minBarcodeLength)&&(p<=maxBarcodeLength)) {
                String[] seqs=new String[2];
                seqs[0]=seq.substring(0,p);
                seqs[1]=seq.substring(p, p+seqLength);
                return seqs;
            }
        }
        return null;
    }

    public void processDirectoryForCounting(File inDirectory, File outDirectory, String suffix) {
       int size=40000000;
        for(File fn: inDirectory.listFiles())  {
           haplotype=new long[2][size];
//           determineFileRows(fn);
           readSolexaFastaToCountList(fn, false);
           GenericSorting.quickSort(0, haplotype[0].length, comp, swapper);
           File outFile=new File(outDirectory.getPath()+"\\"+fn.getName()+suffix);
           writeCountFile(outFile, 1, true);
       }
    }

    public int determineFileRows(File currentFile) {
        int count=0;
        String line;
        try{
        FileInputStream  fr=new FileInputStream (currentFile);
        System.out.println("File = " + currentFile);
            BufferedReader br=new BufferedReader(new InputStreamReader(fr, "ISO-8859-1"),65536);
           while ((line = br.readLine()) != null) {
                count++;
                if(count%500000==0)  System.out.println("Total Reads="+count+" linelength"+line.length());
            }
        br.close();
        }
        catch(Exception e) {
            System.out.println("Catch c="+count+" e="+e);
        }
       System.out.println("Count of Tags="+count);
       return count;
    }

    public void readSolexaFastaToCountList(File currentFile, boolean fastq) {
        int goodBarcodedReads=0, totalReads=0;;
        try{
        FileInputStream  fr=new FileInputStream (currentFile);
        System.out.println("File = " + currentFile);
        BufferedReader br=new BufferedReader(new InputStreamReader(fr),65536);
            String temp, sl;
            while (((temp = br.readLine()) != null)&&(goodBarcodedReads<haplotype[0].length)) {
                if(fastq) {  //load fastq format
                    if(!temp.startsWith("@")) continue;    //fastq starts with @ rather than >
                    sl=br.readLine();
                    if(br.readLine().startsWith("+")) {
                        br.readLine();  //skip the quality scores
                    }
                }  else { //load prefiltered data
                    String[] jj=temp.split("\t");
                    sl=jj[8];
                }
                totalReads++;
                String[] seqs=findBarcodeHaplotype(sl,2*chunkSize);
                if(seqs==null) continue;  //overhang missing so skip
                if(seqs[1].indexOf('N')>-1) continue;  //bad sequence so skip
                if(seqs[1].indexOf('.')>-1) continue;  //bad sequence so skip
                String hap=removeSeqAfterSecondCutSite(seqs[1]);
                haplotype[0][goodBarcodedReads]=BaseEncoder.getLongFromSeq(hap.substring(0,chunkSize));
                haplotype[1][goodBarcodedReads]=BaseEncoder.getLongFromSeq(hap.substring(chunkSize,2*chunkSize));
                goodBarcodedReads++;
                if(goodBarcodedReads%100000==0)  System.out.println("Total Reads="+totalReads+" good read= " + goodBarcodedReads );
            }
        fr.close();
        }
        catch(Exception e) {
            System.out.println("Catch c="+goodBarcodedReads+" e="+e);
        }
       System.out.println("Count of Tags="+goodBarcodedReads);
    }

   void writeCountFile(File outFile, int minCount, boolean binary) {
       int c=1;
       int hapsOutput=0;
       try{
       DataOutputStream fw=new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outFile),4000000));       
       for(int i=0; i<haplotype[0].length-1; i++) {
         if((haplotype[0][i]==haplotype[0][i+1])&&(haplotype[1][i]==haplotype[1][i+1])) {
             c++;
         }
         else {
         if(c>=minCount) {
             if(!binary) {fw.writeBytes(
             BaseEncoder.getSequenceFromLong(haplotype[0][i],(byte)chunkSize)+
             BaseEncoder.getSequenceFromLong(haplotype[1][i],(byte)chunkSize)+" "+c+"\n");}
             else {fw.writeLong(haplotype[0][i]);
                fw.writeLong(haplotype[1][i]);
                fw.writeInt(c);
             }
         }
         hapsOutput++;
          c=1;
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
       SolexaFreqCreatorV2 be=new SolexaFreqCreatorV2();
  }
}
