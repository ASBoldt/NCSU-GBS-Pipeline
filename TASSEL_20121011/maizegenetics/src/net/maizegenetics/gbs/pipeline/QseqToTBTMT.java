/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.gbs.pipeline;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import net.maizegenetics.gbs.homology.ParseBarcodeRead;
import net.maizegenetics.gbs.homology.ReadBarcodeResult;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.tagdist.TagsByTaxaBit;
import org.biojava3.core.util.ConcurrencyTools;
/**
 *Multithreaded pipeline to convert a series of Qseq files to TagsByTaxa file.  It requires a list of
 * existing tags, which may come from a TOPM file or TagCount file.
 *
 * This is a multiple threaded version, it is not much faster (2-3x) than the other one, but once
 * nio2 becomes available it should be much faster.
 * @author james/ed
 */
public class QseqToTBTMT {
    int goodBarcodedReads=0, allReads=0, goodMatched=0;

    private QseqToTBTMT(String keyFileS, String enzyme, String[] qseqFileS, String tagMapFileS, String outReadByTaxaS) {
        testBasicPipeline(keyFileS, enzyme, qseqFileS, tagMapFileS, outReadByTaxaS); //Feed all data into function
    }

    public static void main(String[] args) {
//	String qseqDirectory = new String("/media/CA568FCE568FB9A9/NAM");
//        String keyFileS = new String("/home/jvh39/Documents/code/gbs_pipeline/NAM_IBM_MDP_Ames_key_20110310.txt");
//	String tagMapFileS = new String("/home/jvh39/Documents/code/gbs_pipeline/14FCGBS.tg.ndup.bin");
//        String qseqDirectory = new String("/Volumes/ESB_GBS1/qseq/NAM");
//        String keyFileS = new String("/Volumes/ESB_GBS1/qseq/NAM_IBM_MDP_Ames_key_20110310.txt");
//	String tagMapFileS = new String("/Users/edbuckler/SolexaAnal/GBS/reftags/14FCGBS.tg.ndup.bin");
//        String outReadByTaxaS = new String("/Users/edbuckler/SolexaAnal/GBS/test/NAM110315.tbt.txt"); //Put output file in same directory as input, use first input name, and append ".result" to it
	String[] qseqFileS ={"/Users/edbuckler/SolexaAnal/NGG_IBM/434LFAAXX/qseq/434LFAAXX_2_qseq.txt"};
        String keyFileS ="/Users/edbuckler/SolexaAnal/GBS/NAM_IBM_MDP_key.txt";
        String tagMapFileS = "/Users/edbuckler/SolexaAnal/GBS/reftags/14FCGBS.tg.ndup.bin";
        String outReadByTaxaS="/Users/edbuckler/SolexaAnal/GBS/test/junk.tbt.txt";
        String enzyme = "ApeKI";
//        File[] testFiles = DirectoryCrawler.listFiles("_*qseq.txt", new File(qseqDirectory));
//        String[] qseqFileS = new String[testFiles.length];
//	for(int i=7; i<testFiles.length; i++){
//            qseqFileS[i] = testFiles[i].getAbsolutePath(); //Put names of found files into array
//            System.out.println(qseqFileS[i]);
//        }

        QseqToTBTMT theQHPMT=new QseqToTBTMT(keyFileS, enzyme, qseqFileS, tagMapFileS, outReadByTaxaS);
        
    }


 

    /**@param keyFileS A key file (list of populations by barcode).
     * @param qseqFileS qseq files (Illumina-created files with raw read sequence, quality score, machine name, etc.)
     * @param tagMapFileS A tag map file (reads, trimmed to cut sites, with their genome coordinates).
     * @param outReadByTaxaS A tag-by-taxa file (read sequence and presence/absence in one or more taxa).*/
    public void testBasicPipeline(String keyFileS, String enzyme, String[] qseqFileS, String tagMapFileS, String outReadByTaxaS) {
	final int cores=Runtime.getRuntime().availableProcessors();
        System.out.println("testBasicPipeline processors available:"+cores);
        ExecutorService pool = Executors.newCachedThreadPool();
        ThreadPoolExecutor tpe=(ThreadPoolExecutor)pool;
        TagsOnPhysicalMap theTOPM=new TagsOnPhysicalMap(tagMapFileS,true);  //Reads tag map file into an object
        theTOPM.sortTable(true);
        int uniqueCnt=0;
        int pauses=0;
//        for (int i = 0; i < theTOPM.getSize()-1; i++) {
//            if(!Arrays.equals(theTOPM.getTag(i), theTOPM.getTag(i+1))) uniqueCnt++;
//        }
        System.out.println("The Physical Map Files has UniqueTags:"+uniqueCnt+" TotalLocations:"+theTOPM.getSize());
        TagsByTaxa theTBT=null;
//        int goodBarcodedReads=0, allReads=0, goodMatched=0;
        for(int laneNum=0; laneNum<qseqFileS.length; laneNum++) {
            System.out.println("Lane qseq: "+qseqFileS[laneNum]);
            File qseqFile=new File(qseqFileS[laneNum]);
            String[] np=qseqFile.getName().split("_");
            ParseBarcodeRead thePBR;
            if(np.length==3) {thePBR=new ParseBarcodeRead(keyFileS, enzyme, np[0], np[1]);}
            else if(np.length==5) {thePBR=new ParseBarcodeRead(keyFileS, enzyme, np[1], np[3]);}
            else {
                System.out.println("Error in parsing file name:");
                System.out.println("   The filename does not contain either 3 or 5 underscore-delimited values.");
                System.out.println("   Expect: flowcell_lane_qseq.txt OR code_flowcell_s_lane_qseq.txt");
                System.out.println("   Filename: "+qseqFileS[laneNum]);
                return;
            }
            System.out.println("Total barcodes found in flowcell lane:"+thePBR.getBarCodeCount());
            String[] taxaNames=new String[thePBR.getBarCodeCount()];
            for (int i = 0; i < taxaNames.length; i++) {
                taxaNames[i]=thePBR.getTheBarcodes(i).getTaxaName();
            }
            if(laneNum==0) {theTBT=new TagsByTaxaBit(taxaNames,theTOPM);}
            else {theTBT.addTaxa(taxaNames);}
            int max=200000000;
            String temp="";
            goodBarcodedReads=0; allReads=0; goodMatched=0;
            long time=System.currentTimeMillis();
            try{
                System.out.println("qseq File = " + qseqFileS[laneNum]);
//                Charset charset = Charset.forName("US-ASCII");
//                Path pt=Paths.get(qseqFileS[flowNum]);
//                BufferedReader br = java.nio.file.Files.newBufferedReader(pt, charset);  //this will work in Java 7, but current Mac version does not have this method implemented
                BufferedReader br=new BufferedReader(new FileReader(qseqFileS[laneNum]));
                time=System.currentTimeMillis();
                int cacheMax=500;
                String[] cacheReadLines=new String[cacheMax];
                int cacheIndex=0;
                while (((temp = br.readLine()) != null)&&(goodBarcodedReads<max)) {
                      cacheReadLines[cacheIndex]=temp;
                      cacheIndex++;
                      if(cacheIndex==cacheReadLines.length) {
                          ProcessQseqLine pql=new ProcessQseqLine(thePBR, theTBT, cacheReadLines);
                          tpe.submit(pql);
                          cacheReadLines=new String[cacheMax];
                          cacheIndex=0;
                          while(tpe.getActiveCount()>cores) {
                               try{Thread.sleep(0,100000); pauses++;} catch (Exception e) {e.printStackTrace();}
                           }
                      }
                    allReads++;
                    if(allReads%1000000==0) {
                        double rate=(double)allReads/(double)(System.currentTimeMillis()-time);
                        System.out.println("Total Reads:"+allReads+" goodReads:"+goodBarcodedReads+" goodMatched:"+goodMatched+" Rate:"+rate+" reads/millisec thread"+tpe.getActiveCount());
//                        System.out.println(sl);
                    }
                }
                br.close();
                System.out.println("Main ThreadsCnt:"+Thread.activeCount()+" AlignmentThreadsCnt:"+ConcurrencyTools.getThreadPool().getActiveCount());
                try {// Wait a while for existing tasks to terminate
                     if (!pool.awaitTermination(5, TimeUnit.SECONDS)) {
                       pool.shutdownNow(); // Cancel currently executing tasks
                       // Wait a while for tasks to respond to being cancelled
                       if (!pool.awaitTermination(5, TimeUnit.SECONDS)) {System.err.println("Pool did not terminate");}
                       else {System.out.println("Pool did terminate");}
                     }
                   } catch (InterruptedException ie) {
                       System.err.println("Pool did not terminate");
                     // (Re-)Cancel if current thread also interrupted
                     pool.shutdownNow();
                     // Preserve interrupt status
                     Thread.currentThread().interrupt();
                }
                System.out.println("TC:"+Thread.activeCount()+" BJC"+ConcurrencyTools.getThreadPool().getActiveCount());
            }
            catch(Exception e) {
                System.err.println("Catch testBasicPipeline c="+goodBarcodedReads+" e="+e);
                System.err.println(temp);
                System.err.println(temp.length());
                e.printStackTrace();
            }
            if(laneNum%10==0) theTBT.writeDistFile(new File(outReadByTaxaS.replace(".tbt.","."+laneNum+".tbt.")),FilePacking.Bit,5);
        }
       theTBT.writeDistFile(new File(outReadByTaxaS.replace(".txt", ".bibin")),FilePacking.Bit,5);
       theTBT.writeDistFile(new File(outReadByTaxaS),FilePacking.Text,5);
       System.out.println("Count of Tags="+goodBarcodedReads);

    }

    private class ProcessQseqLine extends Thread {
     //   String readLine;
        ParseBarcodeRead thePBR;
        TagsByTaxa theTBT;
        String[] cacheReadLines;
        int seqLength;

        public ProcessQseqLine(ParseBarcodeRead thePBR, TagsByTaxa theTBT, String[] cacheReadLines) {
           // this.readLine=readLine;
            this.thePBR=thePBR;
            this.theTBT=theTBT;
            this.cacheReadLines=cacheReadLines;
            String[] jj=cacheReadLines[0].split("\t");
            seqLength=jj[8].length();
        }

        public void run() {
            for(String readLine: cacheReadLines) {
                if(readLine==null) continue;
//                String[] jj=readLine.split("\t");
//                String sl=jj[8];
//                String qualS=jj[9];
                int qualEnd=readLine.length()-2;
                int qualStart=qualEnd-seqLength;
                int seqEnd=qualStart-1;
                int seqStart=seqEnd-seqLength;
                String sl=readLine.substring(seqStart, seqEnd);
                String qualS=readLine.substring(qualStart, qualEnd);
                ReadBarcodeResult rr=thePBR.parseReadIntoTagAndTaxa(sl, qualS, false, 0);
                if(rr!=null) {
                    goodBarcodedReads++;
                    int t=theTBT.getIndexOfTaxaName(rr.getTaxonName());
                    int h=theTBT.getTagIndex(rr.getRead());
                    if(h>-1) {
                        theTBT.addReadsToTagTaxon(h, t, 1);
                        goodMatched++;
                    }
                }
            }

        }
    }
  
}

