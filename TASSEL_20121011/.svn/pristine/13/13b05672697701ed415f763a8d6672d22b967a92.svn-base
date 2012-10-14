package net.maizegenetics.genome.GBS.switchgrass;

import net.maizegenetics.genome.GBS.*;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.util.zip.GZIPInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import net.maizegenetics.genome.BaseEncoder;


/**
 * Takes a directory of qseq files and processes the file into count files by taxa
 * based on a key file
 * Developer: ed
 * 
 */
public class ParseBarcodeFiles {
    private static int chunkSize=BaseEncoder.chunkSize;
    //String baseDir="E:/SolexaAnal/";
    private int maximumMismatchInBarcodeAndOverhang=0;
	public static final String[] overhang={"CAGC","CTGC"};  // for ApeKI
//  public static final String[] overhang={"TGCAG"}; // for PstI
//	public static final String[] overhang = {"TGCAT"}; // for Ecot22I
    static String nullS="AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";

	public static final String[] likelyReadEnd={"GCTGAGAT", "GCAGAGAT","GCAGC","GCTGC"};  // ApeKI
//  public static final String[] likelyReadEnd={"CTGCAAGATC", "CTGCAG"};  // PstI 
//	public static final String[] likelyReadEnd={"ATGCAAGAT","ATGCAT"};  // Ecot22I

    //Original design CWGGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG  -- Common Adapter
    //Redesign  CWGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG  -- Common Adapter

    static int maxBarcodeLength=8;
    private Barcode[] theBarcodes;
    private long[] quickBarcodeList;
	private HashMap<Long,Integer> quickMap;

    private int nShortSeqs = 0;


    public ParseBarcodeFiles(String inDir, String keyFile, String outDir, int minQual, boolean fastq) {
        processDirectoryForCounting(new File(inDir),new File(keyFile), new File(outDir), minQual, fastq);
    }



    private int setupBarcodeFiles(File keyFile, File outDir, String flowcell, String lane) {
        try{
            BufferedReader br=new BufferedReader(new FileReader(keyFile),65536);
            ArrayList<Barcode> theBarcodesArrayList=new ArrayList<Barcode>();
            String temp;
            while (((temp = br.readLine()) != null)) {
                String[] s=temp.split("\\s");  //split by whitespace
                if(s[0].equals(flowcell)&&s[1].equals(lane)) {
                    File nf=new File(outDir.getPath()+"/"+s[3]+"_"+flowcell+"_"+lane+"_"+s[2]+".cnt");
                    DataOutputStream dos=new DataOutputStream(new BufferedOutputStream(new FileOutputStream(nf),65536));
                    Barcode theBC=new Barcode(s[2],overhang,s[3], flowcell, lane);
                    theBC.setTheDOS(dos);
                    theBarcodesArrayList.add(theBC);
                    System.out.println(theBC.barcodeS+" "+theBC.taxa);
                }
            }
            theBarcodes=new Barcode[theBarcodesArrayList.size()];
            theBarcodesArrayList.toArray(theBarcodes);
			Arrays.sort(theBarcodes);
            int nBL = theBarcodes[0].barOverLong.length;
            quickBarcodeList=new long[theBarcodes.length*nBL];
            quickMap = new HashMap();
            for (int i = 0; i < theBarcodes.length; i++) {
            	for (int j=0; j<nBL; j++){
                    quickBarcodeList[i*nBL+j]=theBarcodes[i].barOverLong[j];
                    quickMap.put(theBarcodes[i].barOverLong[j], i);
            	}

            }
/*
            quickBarcodeList=new long[theBarcodes.length];
            for(int i=0; i<theBarcodes.length; i++) {
                quickBarcodeList[i]=theBarcodes[i].barOverLong[0];
                //System.out.println(quickBarcodeList[i]+" "+BaseEncoder.getSequenceFromLong(quickBarcodeList[i]));
            }
 */
        } catch(Exception e){
            System.out.println("Error with setupBarcodeFiles: "+e);
        }
        return theBarcodes.length;
    }
    private Barcode findBestBarcode(String queryS, int maxDivergence) {
        long query=BaseEncoder.getLongFromSeq(queryS.substring(0, chunkSize));
        int closestHit=Arrays.binarySearch(quickBarcodeList, query);
		if (closestHit < -1) {
			int index = quickMap.get(quickBarcodeList[-(closestHit+2)]);
			if (theBarcodes[index].compareSequence(query, 1)==0) {
				return theBarcodes[index];
			}
		}
		else if(maxDivergence==0) return null;
/*
        if((closestHit<-1)&&(theBarcodes[-(closestHit+2)].compareSequence(query, 1)==0)) {
            return theBarcodes[-(closestHit+2)];
        } else if(maxDivergence==0) return null; 
 */
        int maxLength=0, minDiv=maxDivergence+1, countBest=0;
        Barcode bestBC=null;
        for(Barcode bc: theBarcodes) {
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
     public String removeSeqAfterSecondCutSite(String seq) {
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
            seq=seq.substring(0,pos+overhang[0].length())+nullS.substring(0, slen-pos-overhang[0].length());//this makes all bases after the second site AAAAAAA
                //we cut after base 4 as it could ligated to an adapter or another sequence
                //this is the length of the overhang+1

              //   System.out.println(" rseq:"+seq);
            ++nShortSeqs;
        }
        return seq;
    }

    private void divideSolexaFileToTaxaCountList(File currentFile, boolean fastq, int minQual) {
        int goodBarcodedReads=0, highQual=0, totalReads=0;
        int max=(int)100e7;
        try{
            System.out.println("File = " + currentFile);
            
            //Read in qseq file as a gzipped text stream if its name ends in ".gz", otherwise read as text
            BufferedReader br;
            if(currentFile.getName().endsWith(".gz")){
                br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(currentFile))));
            }else{
                br=new BufferedReader(new FileReader(currentFile),65536);
            }
            
            String temp, sl, qualS="";
            long h0, h1;
            while (((temp = br.readLine()) != null)&&(goodBarcodedReads<max)) {
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
                if(firstBadBase>(64+8)) {highQual++;} else {continue;}
                int miss = -1;
                if (fastq) { miss=sl.indexOf('N'); } else { miss=sl.indexOf('.'); }
                if((miss!=-1)&&(miss<(maxBarcodeLength+2*chunkSize))) continue;  //bad sequence so skip
                Barcode bestBarcode=findBestBarcode(sl,maximumMismatchInBarcodeAndOverhang);
                if(bestBarcode==null) continue;  //overhang missing so skip
                String genomicSeq=sl.substring(bestBarcode.barLength, sl.length());
                String hap=removeSeqAfterSecondCutSite(genomicSeq);  //this is slow 20% of total time
  //              String hap=genomicSeq;
                h0=BaseEncoder.getLongFromSeq(hap.substring(0,chunkSize));
                h1=BaseEncoder.getLongFromSeq(hap.substring(chunkSize,2*chunkSize));
//                bestBarcode.getTheDOS().writeBytes(sl+" 0"+"\n");
                bestBarcode.getTheDOS().writeLong(h0);
                bestBarcode.getTheDOS().writeLong(h1);
                bestBarcode.getTheDOS().writeInt(1);
                goodBarcodedReads++;
                if(goodBarcodedReads%100000==0)  System.out.println("Total Reads="+totalReads+
                         "  HighQuality="+highQual+ " barcoded= " + goodBarcodedReads);
            }
            br.close();
        }
        catch(Exception e) {
            System.out.println("Catch c="+goodBarcodedReads+" e="+e);
        }
       System.out.println("Count of Tags="+goodBarcodedReads);
    }

    public void processDirectoryForCounting(File inDirectory, File keyFile, File outDir, int minQual, boolean fastq) {
       try{
            for(File fn: inDirectory.listFiles())  {
                System.out.println("Lane qseq: "+fn.toString());
                String[] np=fn.getName().split("_");
                int totalBarcodes;
                if (np.length==3) {
                    totalBarcodes=setupBarcodeFiles(keyFile, outDir, np[0], np[1]);
                }
                else if (np.length==5) {
                    totalBarcodes=setupBarcodeFiles(keyFile, outDir, np[1], np[3]);
                }
                else {
                    System.out.println("Error in parsing file name (e.g. flowcell_lane_qseq.txt).  The following"+
                            " file doesn't contain either 3 or 5 underscore-delimited words:"+fn.toString());
                    continue; // skip over files that don't match the naming convention
                }
                System.out.println("Total barcodes found in key file for this lane:"+totalBarcodes);
                divideSolexaFileToTaxaCountList(fn, fastq, minQual);
                for(Barcode bc: theBarcodes) bc.getTheDOS().close();
            }
            System.out.println("TOTAL short seqs (<64) found on flow cell: " +nShortSeqs);
       } catch (Exception e) {
           System.out.println("Error processDirectoryForCounting:"+e);
       }
    }

   public static void main(String[] args) {
       if(args.length!=3) {
           System.out.println("Require args: inDirectory keyFile outDirectory");
           System.exit(1);
       }
       ParseBarcodeFiles be=new ParseBarcodeFiles(args[0], args[1], args[2],0,false);
  }
}

class Barcode implements Comparable<Barcode> {
    String barcodeS;
    String[] overhangS;
    String taxa, flowcell, lane;
    long[] barOverLong;
    int barOverLength, barLength;
    DataOutputStream theDOS=null;


    public Barcode(String barcodeS, String[] overhangSunsort, String taxa, String flowcell, String lane) {
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

    public DataOutputStream getTheDOS() {
        return theDOS;
    }

    public void setTheDOS(DataOutputStream theDOS) {
        this.theDOS = theDOS;
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
    public int compareTo(Barcode anotherBarcode) {
        if(this.barOverLong[0]<anotherBarcode.barOverLong[0]) return -1;
        if(this.barOverLong[0]>anotherBarcode.barOverLong[0]) return 1;
        return 0;
    }
}
