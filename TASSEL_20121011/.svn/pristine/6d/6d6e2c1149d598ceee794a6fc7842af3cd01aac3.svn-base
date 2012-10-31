package net.maizegenetics.genome.GBS;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeMap;
import net.maizegenetics.genome.BaseEncoder;


/**
 * Takes a directory of qseq files and processes the file into count files by taxa
 * based on a key file
 * Developer: ed
 * 
 */
public class QseqToReadByTaxa {
    private static int chunkSize=BaseEncoder.chunkSize;
    //String baseDir="E:/SolexaAnal/";
    private static int maximumMismatchInBarcodeAndOverhang=1;
    private static final String[] overhang={"CAGC","CTGC"};
    private static String nullS="AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
//    static String[] fullCutSites={"GCAGC","GCTGC"};
    public static final String[] likelyReadEnd={"GCTGGATC","GCAGGATC","GCTGAGAT", "GCAGAGAT","GCAGC","GCTGC"};
    private static int maxLengthOnReadEnd=9;

    //Original design CWGGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG  -- Common Adapter
    //Redesign  CWGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG  -- Common Adapter

    private static int maxBarcodeLength=8;
//    private QBarcode[] theBarcodes;
//    private long[] quickBarcodeList;



    public QseqToReadByTaxa(String inDir, String keyFile, ReadsByTaxa outReadsByTaxa, int minQual, boolean fastq) {
        processDirectoryForCounting(new File(inDir),new File(keyFile), outReadsByTaxa, minQual, fastq);
    }



    private static QBarcode[] setupBarcodeFiles(File keyFile, String flowcell, String lane, ReadsByTaxa theRBC) {
        QBarcode[] theBarcodes=null;
        try{
            BufferedReader br=new BufferedReader(new FileReader(keyFile),65536);
            ArrayList<QBarcode> theBarcodesArrayList=new ArrayList<QBarcode>();
            String temp;
            while (((temp = br.readLine()) != null)) {
                String[] s=temp.split("\\s");  //split by whitespace
                if(s[0].equals(flowcell)&&s[1].equals(lane)) {
//                    File nf=new File(outDir.getPath()+"/"+s[3]+"_"+flowcell+"_"+lane+".cnt");
//                    DataOutputStream dos=new DataOutputStream(new BufferedOutputStream(new FileOutputStream(nf),65536));
                    QBarcode theBC=new QBarcode(s[2],overhang,s[3], flowcell, lane);
                    theBC.rbtTaxaIndex=theRBC.getIndexOfTaxaName(s[3]);
             if(theBC.rbtTaxaIndex<0) System.out.println("NotRecognizedInRefList:"+temp);

//                    theBC.setTheDOS(dos);
                    theBarcodesArrayList.add(theBC);
                    System.out.println(theBC.barcodeS+" "+theBC.taxa);
                }
            }
            theBarcodes=new QBarcode[theBarcodesArrayList.size()];
            theBarcodesArrayList.toArray(theBarcodes);
            Arrays.sort(theBarcodes);
 //           quickBarcodeList=new long[theBarcodes.length];
//            for(int i=0; i<theBarcodes.length; i++) {
//                quickBarcodeList[i]=theBarcodes[i].barOverLong[0];
//                //System.out.println(quickBarcodeList[i]+" "+BaseEncoder.getSequenceFromLong(quickBarcodeList[i]));
//            }
        } catch(Exception e){
            System.out.println("Error with setupBarcodeFiles: "+e);
        }
        return theBarcodes;
    }

    private static QBarcode findBestBarcode(QBarcode[] theBarcodes, long[] quickBarcodeList,
            String queryS, int maxDivergence) {
        long query=BaseEncoder.getLongFromSeq(queryS.substring(0, chunkSize));
        int closestHit=Arrays.binarySearch(quickBarcodeList, query);
        if((closestHit<-1)&&(theBarcodes[-(closestHit+2)].compareSequence(query, 1)==0)) {
            return theBarcodes[-(closestHit+2)];
        } else if(maxDivergence==0) return null;
        int maxLength=0, minDiv=maxDivergence+1, countBest=0;
        QBarcode bestBC=null;
        for(QBarcode bc: theBarcodes) {
            int div=bc.compareSequence(query, maxDivergence+1);
 //           System.out.println(queryS+" "+bc.barcodeS+" "+div);
            if(div<=minDiv) {
                if((div<minDiv)||(bc.barOverLength>maxLength)) {
                    minDiv=div;
                    maxLength=bc.barOverLength;
                    bestBC=bc;
                    countBest=1;
                   //if(minDiv==0) break;
                } else {  //it is a tie, so return that not resolvable
                    bestBC=null;
                    countBest++;
                }
            }
        }
        return bestBC;
    }

    /**
     * The barcode libaries used for this study can include two types of extraneous sequence
     * at the end of reads.  The first are chimeras created with the free ends.  These will
     * recreate the restriction site.  The second are short regions (less than 64bp), so that will they
     * will contain a portion of site and the univeral adapter.
     * This finds the first of site in likelyReadEnd, keeps the restriction site overhang and then sets everything
     * to polyA afterwards
     * @param seq
     * @return processed sequence
     */
     public static String removeSeqAfterSecondCutSite(String seq) {
         //this looks for a second restriction site, and then turns the remaining sequence to AAAA
        int pos=9999;
        int startSearchBase=3;
        for(String fcs: likelyReadEnd) {
            int p=seq.indexOf(fcs, startSearchBase);
            if((p>startSearchBase)&&(p<pos)) {
                pos=p;
            }
        }
        if(pos<chunkSize*2) {
            int slen=seq.length();
              //  System.out.print("seq:"+seq);
            seq=seq.substring(0,pos+4)+nullS.substring(0, slen-pos-4);//this makes all bases after the second site AAAAAAA
                //we cut after base 4 as it could ligated to an adapter or another sequence
                //this is the length of the overhang+1

              //   System.out.println(" rseq:"+seq);
        }
        return seq;
    }

    private static void divideSolexaFileToTaxaCountList(QBarcode[] theBarcodes, long[] quickBarcodeList,
            File inFile, ReadsByTaxa theRBT, ReadBLASTer theReadBLASTer, boolean fastq, int minQual) {
        int barcodedReads=0, highQualReads=0, totalReads=0, readsInRBT=0, noPeriodReads=0;
        int max=(int)100e6;
        try{
        FileInputStream  fr=new FileInputStream (inFile);
        System.out.println("File = " + inFile);
        BufferedReader br=new BufferedReader(new InputStreamReader(fr),65536);
            String temp, sl, qualS="";
            long[] h=new long[2];
            while (((temp = br.readLine()) != null)&&(barcodedReads<max)) {
                if(fastq) {  //load fastq format
                    if(!temp.startsWith("@")) continue;    //fastq starts with @ rather than >
                    sl=br.readLine();
                    if(br.readLine().startsWith("+")) {
                        qualS = br.readLine();
                    }
                }  else { //load prefiltered data
                    String[] jj=temp.split("\t");
                    sl=jj[8];
                    qualS=jj[9];
                }
                totalReads++;
                int firstBadBase=BaseEncoder.getFirstLowQualityPos(qualS, minQual);
                if(firstBadBase>(64+8)) {highQualReads++;} else {continue;}
                int miss = -1;
                if (fastq) { miss=sl.indexOf('N'); } else { miss=sl.indexOf('.'); }
                if((miss!=-1)&&(miss<(maxBarcodeLength+2*chunkSize))) continue;  //bad sequence so skip
                noPeriodReads++;
                QBarcode bestBarcode=findBestBarcode(theBarcodes, quickBarcodeList,sl,maximumMismatchInBarcodeAndOverhang);
                if(bestBarcode==null) continue;  //overhang missing so skip
                barcodedReads++;
//        System.out.println(temp);
//        System.out.println("s1.length="+sl.length());
  //              String genomicSeq=sl.substring(bestBarcode.barLength, bestBarcode.barLength+2*chunkSize+maxLengthOnReadEnd);
                 String genomicSeq=sl.substring(bestBarcode.barLength, sl.length());
                String hap=removeSeqAfterSecondCutSite(genomicSeq);  //this is slow 20% of total time
  //              String hap=genomicSeq;
                h[0]=BaseEncoder.getLongFromSeq(hap.substring(0,chunkSize));
                h[1]=BaseEncoder.getLongFromSeq(hap.substring(chunkSize,2*chunkSize));
//                bestBarcode.getTheDOS().writeBytes(sl+" 0"+"\n");
                int rbtIndex=theRBT.getReadIndex(h);

                TreeMap<Integer, Integer> al=theReadBLASTer.findMatchesWithIntLengthWords(h, 36, true);
                if(al.size()!=1) continue;
                if(bestBarcode.rbtTaxaIndex<0) continue;
//                if(rbtIndex<0) continue;
//                theRBT.addToReadCountForTaxa(rbtIndex, 1, 1);
                theRBT.addToReadCountForTaxa(al.firstKey(), bestBarcode.rbtTaxaIndex, 1);
                readsInRBT++;
                if(barcodedReads%100000==0)  System.out.println("Total Reads="+totalReads+
                         "  HighQuality="+highQualReads+" noPeriodReads="+noPeriodReads+
                         " barcoded= " + barcodedReads+
                         " MatchInRBT="+readsInRBT);
            }
        fr.close();
        }
        catch(Exception e) {
            System.out.println("Catch c="+barcodedReads+" e="+e.getMessage());
            e.printStackTrace();
        }
       System.out.println("Count of Tags="+barcodedReads);
    }

    public static void processDirectoryForCounting(File inDirectory, File keyFile,
        ReadsByTaxa outRBT, int minQual, boolean fastq) {
        QBarcode[] theBarcodes;
        long[] quickBarcodeList=new long[48];
        ReadBLASTer theReadBLASTer = new ReadBLASTer(outRBT);
        try{
            for(File fn: inDirectory.listFiles())  {
                System.out.println("Lane qseq: "+fn.toString());
                String[] np=fn.getName().split("_");
                if(np.length!=3) {
                    System.out.println("Error in parsing file name (e.g. flowcell_lane_qseq.txt):"+fn.toString());
                    continue; // skip over files that don't match the naming convention
                }
                theBarcodes=setupBarcodeFiles(keyFile, np[0], np[1], outRBT);
                quickBarcodeList=new long[theBarcodes.length];
                for(int i=0; i<theBarcodes.length; i++) {
                    quickBarcodeList[i]=theBarcodes[i].barOverLong[0];
                    //System.out.println(quickBarcodeList[i]+" "+BaseEncoder.getSequenceFromLong(quickBarcodeList[i]));
                    }
                int totalBarcodes=theBarcodes.length;
                System.out.println("Total barcodes found in flowcell:"+totalBarcodes);
                divideSolexaFileToTaxaCountList(theBarcodes, quickBarcodeList,fn, outRBT, 
                        theReadBLASTer, fastq, minQual);
            //    for(QBarcode bc: theBarcodes) bc.getTheDOS().close();
            }
       } catch (Exception e) {
           System.out.println("Error processDirectoryForCounting:"+e);
       }
    }


}

class QBarcode implements Comparable<QBarcode> {
    String barcodeS;
    String[] overhangS;
    String taxa, flowcell, lane;
    long[] barOverLong;
    int barOverLength, barLength;
    int rbtTaxaIndex=-1;



    public QBarcode(String barcodeS, String[] overhangSunsort, String taxa, String flowcell, String lane) {
        this.barcodeS=barcodeS;
        Arrays.sort(overhangSunsort);
        this.overhangS=overhangSunsort;
        this.flowcell=flowcell;
        this.lane=lane;
        this.taxa=taxa;
        barOverLong=new long[overhangS.length];
        for(int i=0; i<overhangS.length; i++) {
            barOverLong[i]=BaseEncoder.getLongFromSeq(barcodeS+overhangS[i]);
        }
        barOverLength=barcodeS.length()+overhangS[0].length();
        barLength=barcodeS.length();
    }


    public int compareSequence(long queryLong, int maxDivCheck) {
        int div=barOverLength;
        for(long targetLong: barOverLong) {
            int c=BaseEncoder.seqDifferencesForSubset(targetLong,  queryLong, barOverLength, maxDivCheck);
            if(c<div) div=c;
        }
        return div;
    }

    @Override
    public int compareTo(QBarcode anotherBarcode) {
        if(this.barOverLong[0]<anotherBarcode.barOverLong[0]) return -1;
        if(this.barOverLong[0]>anotherBarcode.barOverLong[0]) return 1;
        return 0;
    }
}
