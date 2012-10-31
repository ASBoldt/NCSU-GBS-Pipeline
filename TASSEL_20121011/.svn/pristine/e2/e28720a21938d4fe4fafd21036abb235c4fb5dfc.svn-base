/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.gbs.pipeline;

import cern.jet.random.Binomial;
import edu.cornell.lassp.houle.RngPack.RandomJava;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.Arrays;
import java.util.Random;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.gbs.util.OpenBitSet;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.tagdist.TagsByTaxaBitFileMap;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByteFileMap;
import net.maizegenetics.gbs.util.SBitAlignment;
import net.maizegenetics.genome.BaseEncoder;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.ids.IdGroup;

/**
 * Maps tags genetically against an anchor map.
 *
 * The tags that are mapped come from a
 * TagsByTaxa file, while the AnchorMap can be a series of HapMap files.
 *
 * A physical location map can be used to annotate the reporting.  The output is the
 * best genetic location for a given tag.
 * 
 * Future ideas - run the 10
 * 
 * @author edbuckler
 */
public class TagCallerAgainstAnchorMT {
    static int minFreq=1;
    static int minComp=10;
    static double ldThreshold=0.000001;
    static boolean requireNCO=false;
    static int ancChrCnt=10;
    
    TagsByTaxa theTBT=null;
    TagsOnPhysicalMap theTOPM=null;
    File outfileSiteAttributes;
    int ldPermutation=1;
    boolean isUsingContigency=true;
    int[][] tbt2anchRedirect;
    SBitAlignment[] theAnchorChrInBits;

   // IdGroup anchorIdGroup;


    public TagCallerAgainstAnchorMT(String tagsByTaxaFileName, String anchorFileName, 
            String blastMapFile, String outfileGoodSitesName, int tagStart, int tagStop) {
        outfileSiteAttributes=new File(outfileGoodSitesName);
        loadAnchorMaps(anchorFileName);
 //       anchorIdGroup=theAnchorChrInBits[0].getIdGroup();
        System.out.println("Alignment Loaded:"+anchorFileName);
        System.out.println("Read TagsByTaxa:"+tagsByTaxaFileName);
        if(tagsByTaxaFileName.contains("byte")) {theTBT=new TagsByTaxaByteFileMap(tagsByTaxaFileName);}
        else {theTBT=new TagsByTaxaBitFileMap(tagsByTaxaFileName);}
        String[] taxaNames=theTBT.getTaxaNames();
        int taxaCnt=taxaNames.length;
        if(tagStart<0) tagStart=0;
        if((tagStop<0)||(tagStop>=theTBT.getTagCount())) tagStop=theTBT.getTagCount();
        tbt2anchRedirect = new int[ancChrCnt][theTBT.getTaxaNames().length];
        for(int c=0; c<ancChrCnt; c++) {
            for (int t = 0; t < taxaCnt; t++) {
                IdGroup anchorIdGroup=theAnchorChrInBits[c].getIdGroup();
                tbt2anchRedirect[c][t] = anchorIdGroup.whichIdNumber(taxaNames[t]);
            }
        }
       // permuteTaxa();
        if(blastMapFile!=null) theTOPM=new TagsOnPhysicalMap(blastMapFile,true);
        run(tagStart, tagStop);

    }

    private void permuteTaxa() {
        Random rnd=new Random();
        for(int c=0; c<ancChrCnt; c++) {
            for (int t = 0; t < tbt2anchRedirect[c].length; t++) {
                int rt=rnd.nextInt(tbt2anchRedirect[c].length);
                int temp=tbt2anchRedirect[c][t];
                tbt2anchRedirect[c][t]=tbt2anchRedirect[c][rt];
                tbt2anchRedirect[c][rt]=temp;
            }
        }
    }

    private void loadAnchorMaps(String anchorFileName) {
        if(anchorFileName.contains("+")||anchorFileName.contains("=")) {
            theAnchorChrInBits=new SBitAlignment[ancChrCnt];
            for (int i = 0; i < ancChrCnt; i++) {
                System.out.println("anchorFileName"+anchorFileName);
                String file=anchorFileName.replace("+",""+(i+1));
                file=file.replace("=",""+(i+1));
                System.out.println("Reading:"+file);
                Alignment a;
                if(file.endsWith(".txt")) {a= ImportUtils.readFromHapmap(file,""+(i+1));}
                else {a= ImportUtils.readFromZip(file);}
 //               if(file.endsWith(".txt")) ExportUtils.writeToZip(a, file.replace("hmp.txt", ""));
                theAnchorChrInBits[i] = new SBitAlignment(a);
                System.out.printf("Chr %s Sites %d Taxa %d %n", theAnchorChrInBits[i].getLocus(0),
                        theAnchorChrInBits[i].getSiteCount(),theAnchorChrInBits[i].getSequenceCount());
            }
        } else {
            Alignment[] a = ImportUtils.readFromHapmap(anchorFileName).getAlignments();
            theAnchorChrInBits=new SBitAlignment[a.length];
            for (int i = 0; i < a.length; i++) {
                theAnchorChrInBits[i] = new SBitAlignment(a[i]);
    //            System.out.printf("Chr %d Sites %d Taxa %d %n", theAnchorChrInBits[i].getChr(),
    //                    theAnchorChrInBits[i].getNumSites(),theAnchorChrInBits[i].getNumTaxa());
                System.out.printf("Chr %s Sites %d Taxa %d %n", theAnchorChrInBits[i].getLocus(0),
                        theAnchorChrInBits[i].getSiteCount(),theAnchorChrInBits[i].getSequenceCount());
            }
        }
    }

    void run(int startTag, int stopTag) {
        if(theTBT==null) return;
        if(theAnchorChrInBits==null) return;
        System.out.println("Starting TagCallerAgainstAnchor on:");
        int newGeneticPositions=0, agreeWithBLAST=0;
        try {
            BufferedWriter fileOutSiteAttr = new BufferedWriter(new FileWriter(outfileSiteAttributes), 100000);
            String row = null;
            int rawCount = 0,
            countGreaterThanMinTaxaWithSNP = 0,
            countGreaterThanMinLogPValue = 0;           
            long starttime = System.currentTimeMillis();
            String header="TestTag\tTestTagNum\tBlastChr\tBlastPos\trefDiv\tLDChr\tLDSite\tLDPos\tBinomP\tSigTests\tTagTaxaCnt\tChrSig\tLRatioB:2\tLRatioB:M\tSiteOnBestChrThanNextChr\tMinSigPos\tMaxSigPos";
            fileOutSiteAttr.write(header+"\n");
            System.out.println(header);
            for (int i = startTag; i < stopTag; i++) {
          //      if(i%100==0) System.out.println("Test Site:"+i+" of "+theTBT.getTagCount());
                if(theTBT.getNumberOfTaxaWithTag(i)<5) continue;
                double[] bestR={-1,-1,-1, 1, -1};
                int blastChr=-1, blastPos=-1, bestAlignment=-1, refDiv=-1;
                long[] testTag=theTBT.getTag(i);
                if(theTOPM!=null) {
                    int index=theTOPM.getTagIndex(testTag);
                    if(index>-1) {
                        blastChr=theTOPM.getChromosome(index);
                        blastPos=theTOPM.getStartPosition(index);
                        refDiv=theTOPM.getDivergence(index);
                    }
                }
// if((blastChr!=10)||(blastPos<30000000)) continue;
                long[] testTagDist;
                //=getTagsInBits(theTBT,i,tbt2anchRedirect,anchorIdGroup.getIdCount());
                double[][] theResults=new double[theAnchorChrInBits.length][];
                ScanChromosomeMTR[] threads = new ScanChromosomeMTR[theAnchorChrInBits.length];
                for (int cn=0; cn<theAnchorChrInBits.length; cn++) {
                    testTagDist=getTagsInBits(theTBT,i,tbt2anchRedirect[cn],theAnchorChrInBits[cn].getSequenceCount());
                    threads[cn] = new ScanChromosomeMTR(theAnchorChrInBits[cn],testTagDist, cn, theResults,0.000001,1, blastPos);
                    threads[cn].start();
                }
                for (int t=0; t<theAnchorChrInBits.length; t++) {
                    try {threads[t].join();}
                    catch (Exception e) {e.printStackTrace();}
                }
                int countRealSig=0;
                double[] pRank=new double[theAnchorChrInBits.length];
                for (int cn = 0; cn < theAnchorChrInBits.length; cn++) {
                    double[] r=theResults[cn];
                    pRank[cn]=r[3];
                    if(r[3]<bestR[3]) {bestR=r.clone(); bestAlignment=cn;}
                    if(r[3]<0.000001) countRealSig++;
                }
                Arrays.sort(pRank);
                double[][] bestResWithNewThreshold=new double[1][];
                testTagDist=getTagsInBits(theTBT,i,tbt2anchRedirect[bestAlignment],theAnchorChrInBits[bestAlignment].getSequenceCount());
                ScanChromosomeMTR bestChrNewThres = new ScanChromosomeMTR(theAnchorChrInBits[bestAlignment],testTagDist, 0, bestResWithNewThreshold,pRank[1],1, blastPos);
                bestChrNewThres.run();
                double ratioOf1to2Chr=Math.log10(pRank[1]/pRank[0]);  //keep if greater than 2
                int countOfSitesBetterThanNextBestChr=(int)bestResWithNewThreshold[0][4];
                double rate=(double)(i-startTag+1)/(double)(System.currentTimeMillis()-starttime);
                String s=String.format("%s %d %d %d %d %d %d %d %g %d %d %d %g %g %d %d %d %n",BaseEncoder.getSequenceFromLong(testTag),i,blastChr, blastPos, refDiv,
                        (int)bestR[0],(int)bestR[1],(int)bestR[2],bestR[3], (int)bestR[4], theTBT.getNumberOfTaxaWithTag(i), 
                        countRealSig, Math.log10(pRank[1]/pRank[0]), Math.log10(pRank[theAnchorChrInBits.length/2]/pRank[0]), countOfSitesBetterThanNextBestChr,
                        bestChrNewThres.minSigPos, bestChrNewThres.maxSigPos);
                if(bestR[3]<0.000001) {/*System.out.print(rate+"\t");*/
                    System.out.print(s);
                    if(blastChr==(int)bestR[0]) {agreeWithBLAST++;}
                    else if(Math.log10(pRank[1]/pRank[0])>3.0) newGeneticPositions++;
                }
                
                fileOutSiteAttr.write(s);
                if(i%100==0) System.out.printf("Rate:%g TagNum:%d BLASTAgree:%d NewPos:%d %n",rate,i,agreeWithBLAST,newGeneticPositions);
            }
            fileOutSiteAttr.close();
            if(theTOPM!=null) theTOPM.writeBinaryFile(new File("/Users/edbuckler/SolexaAnal/GBS/test/14FCGBS_gen110404.tg.ndup.bin"));
            System.out.println("SNP number:" + rawCount +
                                " count over taxa threshold+:" + countGreaterThanMinTaxaWithSNP +
                                " count over logp threshold:" + countGreaterThanMinLogPValue);
            System.out.println( " finished!");
        } catch (Exception e) {
            System.err.println("File IO in TagCallerAgainstAnchor: " + e);
            e.printStackTrace();
        }
    }




    private class ScanChromosomeMTR extends Thread {
        SBitAlignment refAlignment;
        OpenBitSet obsymj;
        Binomial binomFunc=new Binomial(5, 0.5, new RandomJava());
        int resultIndex;
        double[][] resultReport;
        double sigThreshold;
        int minSigPos=-1, maxSigPos=-1;
        int step=1;
        double bionomialThreshold=0.2;
        int blockPosition=-1000;  //position to block from testing.
        int blockWindow=100;

        public ScanChromosomeMTR(SBitAlignment refAlignment, long[] ymj, int resultIndex, double[][] resultReport, double sigThreshold, int step, int blockPosition) {
            this.refAlignment=refAlignment;
            obsymj=new OpenBitSet(ymj,ymj.length);
            this.resultIndex=resultIndex;
            this.resultReport=resultReport;
            this.sigThreshold=sigThreshold;
            this.blockPosition=blockPosition;
            this.step=step;
        }

        public void run() {
            //in the future this may call other threads to split the effort up more
            long tests=0;
            int bestSite=-1, countSig=0;
            double bestP=2;
            for (int i = 0; i < refAlignment.getSiteCount(); i+=step) {
             //   if(i%10000==0) System.out.println("scanChromosome chr:"+refAlignment.getChr()+"+site:"+i+" with Chr:"+refAlignment.getChr()+" test:"+tests);
                if(Math.abs(refAlignment.getPositionInLocus(i)-blockPosition)<blockWindow) continue;
                OpenBitSet mj=refAlignment.getSiteBitsNoClone(i, 0);
                OpenBitSet mn=refAlignment.getSiteBitsNoClone(i, 1);
                if(mn.lastCalculatedCardinality()>4) {
                    double p=testSites(obsymj, mj, mn, binomFunc);
                    if(p<bestP) {bestP=p; bestSite=i;}
                    if(p<sigThreshold) {
                        countSig++;
                        if(minSigPos==-1) minSigPos=refAlignment.getPositionInLocus(i);
                        maxSigPos=refAlignment.getPositionInLocus(i);
                    }
                }
                tests++;
            }
            int chr=Integer.parseInt(refAlignment.getLocus(bestSite).getChromosomeName());
            double[] result={chr, bestSite, refAlignment.getPositionInLocus(bestSite),bestP, countSig};
          //  if(bestP<0.000001) System.out.println(Arrays.toString(result));
            resultReport[resultIndex]=result;
        }

    }


    public static double testSites(OpenBitSet ymj, OpenBitSet xmj, OpenBitSet xmn, Binomial binomFunc) {
        double result=1;
        int mn=0, tag_mn=0, mj=0, tag_mj=0;
        tag_mn=(int)OpenBitSet.intersectionCount(ymj, xmn);
        tag_mj=(int)OpenBitSet.intersectionCount(ymj,xmj);
        mn=(int)xmn.lastCalculatedCardinality();
        mj=(int)xmj.lastCalculatedCardinality();

        int sumAnc = mn + mj;
        int sumTag = tag_mn + tag_mj;
        if(sumTag<4) return result;
        double tag_rat=(tag_mn<tag_mj)?(double)tag_mn/(double)tag_mj:(double)tag_mj/(double)tag_mn;
 //       if(tag_rat>0.2) return result;
        double minorProb = (tag_mn<tag_mj)?(double)mn/(double)sumAnc:(double)mj/(double)sumAnc;
        binomFunc.setNandP(sumTag,minorProb);
        try {
            result = (tag_mn<tag_mj)?binomFunc.cdf(tag_mn):binomFunc.cdf(tag_mj);
  //          if(result<0.000001) System.out.printf("TS: %d %d %d %d %g %g %n",tag_mn, tag_mj, mn, mj, result, tag_rat);
        } catch (Exception e) {
            System.err.println("Error in the BinomialDistributionImpl");
        }
      //   System.out.printf("%d %d %d %d %g %n",c00,c01,c10,c11, fishersExact.getTwoTailedP(c00, c01, c10, c11));
        return result;
     }


     public long[] getTagsInBits(TagsByTaxa aTBT, int tagIndex, int[] reDirect, int anchorTaxa) {
        int lgPerSite = (anchorTaxa / 64) + 1;
        long[] seq = new long[lgPerSite];
        for (int j = 0; j < aTBT.getTaxaCount(); j++) {
            if(reDirect[j]<0) continue;
            int index=reDirect[j]/64;
            int offset=reDirect[j]%64;
            if (aTBT.getReadCountForTagTaxon(tagIndex, j)>0) {  //reference alleles
                seq[index]=seq[index]|(1L<<offset);
            }
        }
        return seq;
    }


     public static void main(String[] args) {
	System.out.println("Running main method in TagCallerAgainstAnchorMT");
//        String blastMapFile="/Users/edbuckler/SolexaAnal/GBS/reftags/14FCGBS.tg.ndup.bin";
        String blastMapFile="/Volumes/LaCie/mergedNam282Ames_072011_mappedOnly.topm.bin";
//        String anchorMapFile="/Users/edbuckler/SolexaAnal/GBS/hmp/HS55K_110215.hmp.txt";
//        String anchorMapFile="/Users/edbuckler/SolexaAnal/GBS/hmp/HiSeq_55K282_110221.imp.hmp.txt";
 //       String anchorMapFile="/Users/edbuckler/SolexaAnal/GBS/hmp/ch*_NAMwLowHetF10LD_110303imp.hmp.txt";
 //       String anchorMapFile="/Users/edbuckler/SolexaAnal/GBS/imp/allfusion_110401.lh.ld.c*.imp.hmp.txt";
  //      String anchorMapFile="/Volumes/LaCie/build20111217/bpec/allZea20111217_scv10mF8maf002_mgs_E1pLD5.c+.hmp.txt";
        String anchorMapFile="/Users/edbuckler/SolexaAnal/GBS/build111217/imp/allZea20111217_scv10mF8maf002_mgs_E1pLD5.imp95_1024.c+.hmp.txt";
 //       String anchorMapFile="/Users/edbuckler/SolexaAnal/GBS/imp/allfusion_110401.lh.ld.c*.imp.BLOB.zip";
 //       String anchorMapFile="/Users/edbuckler/SolexaAnal/GBS/hmp/allfusion_110401.lh.ld.c*.hmp.txt";
  //      String tagsByTaxaFile="/Users/edbuckler/SolexaAnal/GBS/dist/HiSeq282_110214.dist.txt";
//        String tagsByTaxaFile="/Users/edbuckler/SolexaAnal/GBS/dist/NAMsort_110303.bibin";
//        String tagsByTaxaFile="/Users/edbuckler/SolexaAnal/GBS/dist/allfusion.bibin";
        String tagsByTaxaFile="/Volumes/LaCie/allZea20111215.tbt.byte";
    //    String tagsByTaxaFile="/Users/edbuckler/SolexaAnal/GBS/dist/h10000HiSeq282_110214.dist.txt";
        String outfile="/Users/edbuckler/SolexaAnal/GBS/build111217/test/mergedNAM282Ames_mapOnly072011_chr10_30_p4_impblock100.txt";
        
        /*
         * TagCallerAgainstAnchorMT tagsByTaxaFile anchorMapFile outFile blastMapFile (-b) startTag endTag
         */
        int tagStart=Integer.MIN_VALUE, tagStop=Integer.MIN_VALUE;
        if(args.length>=3) {
            tagsByTaxaFile=args[0];
            anchorMapFile=args[1];
            outfile=args[2];
            tagStart=Integer.parseInt(args[3]);
            tagStop=Integer.parseInt(args[4]);
            blastMapFile=null;
        }
        
//        if(args.length>=4) {
//            int offset=0;
//            if(!args[3].equals("-b")) {blastMapFile=args[3]; offset=1;}
//            if(args[3].equals("-b")||args[4].equals("-b")) {
//                tagStart=Integer.parseInt(args[4+offset]);
//                tagStop=Integer.parseInt(args[5+offset]);
//            }
//
//        }
//        blastMapFile=null;
       // tagStart=0; tagStop=2000;
         System.out.printf("PA: %s %s %s %s %d %d %n", tagsByTaxaFile, anchorMapFile, blastMapFile, outfile, tagStart, tagStop);
        TagCallerAgainstAnchorMT tcaa=new TagCallerAgainstAnchorMT(tagsByTaxaFile, anchorMapFile, blastMapFile, outfile, tagStart, tagStop);

    }



}
