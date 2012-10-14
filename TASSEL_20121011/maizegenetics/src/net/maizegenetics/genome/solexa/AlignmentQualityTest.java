/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.solexa;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import net.maizegenetics.genome.HaplotypeLengthDistributionV2;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.Pack1Alignment;
import net.maizegenetics.pal.alignment.SiteSummary;
import net.maizegenetics.pal.datatype.DataType;
import net.maizegenetics.pal.ids.IdGroup;

/**
 *
 * @author edbuckler
 */
public class AlignmentQualityTest {

    public AlignmentQualityTest(String ibdFile, String testAlignmentFile, Alignment intersectAlignment, String outfile) {
        System.out.println("Alignment Loading:"+testAlignmentFile);
        Pack1Alignment p1a = (Pack1Alignment) ImportUtils.createPack1AlignmentFromFile(testAlignmentFile, "", "");
        System.out.println("Alignment Loaded:"+testAlignmentFile);
        System.out.printf("Sites %d Taxa %d %n", p1a.getSiteCount(),p1a.getSequenceCount());
        p1a=(Pack1Alignment)HaplotypeLengthDistributionV2.makeHomozygousAlignment(p1a);
//        Pack1Alignment intersectAlignment=null;
//        if(insectionAlignmentFile!=null) {
//            intersectAlignment = (Pack1Alignment) ImportUtils.createPack1AlignmentFromFile(testAlignmentFile, "", "");
//        }
      //  HaplotypeLengthDistributionV3 hld3=new HaplotypeLengthDistributionV3(p1a, 50, null);

        String[] labels="Chr StartSite EndSite StartPosition EndPosition Taxa1 Taxa2 PerfectMatchLength".split(" ");
        Object[] t={String.class,Integer.class,Integer.class,Integer.class,Integer.class,String.class,String.class,Integer.class};
        SimpleTextFile theSTF=new SimpleTextFile(ibdFile,labels,t,true);

        IdGroup idg=p1a.getIdGroup();
        int depthR=0, depthW=0, tripR=0, tripW=0;
        int totalRight=0, totalMinorRight=0, totalWrong=0;  //all alleles and then just the alternate alleles.
        int tripID=idg.whichIdNumber("TDD39103");
        int[] rightBySiteCnt=new int[p1a.getSiteCount()];
        int[] wrongBySiteCnt=new int[p1a.getSiteCount()];
        int[][] errorByQScore=new int[3][101];
        
        System.out.println("Chr StartSite EndSite StartPosition EndPosition Taxa1 Taxa2 Right Wrong OrigPerfectMatchLength");
        for (int i = 0; i < theSTF.getRows(); i++) {
            //System.out.println(p1a.getLocus(0).toString()+":"+theSTF.getStringElement(0, i));
            if(!p1a.getLocus(0).toString().equals(theSTF.getStringElement(0, i))) continue;
            int right=0, wrong=0, minorRight=0;
            int taxa1=idg.whichIdNumber(theSTF.getStringElement(5, i));
            int taxa2=idg.whichIdNumber(theSTF.getStringElement(6, i));
//            taxa2=1;
            if((taxa1<0)||(taxa2<0)) {
                continue;
            }
            int startSite=p1a.getSiteOfPhysicalPosition(theSTF.getIntElement(3, i), null);
            if(startSite<0) startSite=-(startSite+1);
            int endSite=p1a.getSiteOfPhysicalPosition(theSTF.getIntElement(4, i), null);
            if(endSite<0) endSite=-(endSite+1);
            int currentQ=0;
//            int inset=(endSite-startSite)/4;
//            startSite=startSite+inset;
//            endSite=endSite-inset;
            for (int b = startSite; b < endSite; b++) {
                if((intersectAlignment!=null)&&(intersectAlignment.getSiteOfPhysicalPosition(p1a.getPositionInLocus(b), null)<0)) {
                    System.out.printf("%d %d %s %s %n",b, p1a.getPositionInLocus(b), ""+(char)p1a.getMajorAllele(b), ""+(char)p1a.getMinorAllele(b));
                    continue;
                }
                SiteSummary theSS=p1a.getSiteSummary(b);
                if((theSS.getAlleleCounts().length>1)&&(theSS.getAlleleCounts()[1]<3)) continue;  //this is a filter to eliminate singletons
//                if(theSS.getNumberMissing()>50) continue;
//                if((p1a.getMinorAllele(b)==(byte)'-')||(p1a.getMajorAllele(b)=='-')) continue;

                if(p1a.hasSiteScores()) currentQ=Math.round(p1a.getSiteScore(0, b));
                if(currentQ<0) currentQ=0; //solves a problem with the initialization of the matrix
                if((p1a.getBase(taxa1, b)==DataType.UNKNOWN_BYTE)||(p1a.getBase(taxa2, b)==DataType.UNKNOWN_BYTE)) continue;
                if(p1a.getBase(taxa1, b)==p1a.getBase(taxa2, b)) {
                    right++;
                    errorByQScore[0][currentQ]++;
                    depthR+=p1a.getSequenceCount()-p1a.getSiteSummary(b).getNumberMissing();
                    if(p1a.getBase(tripID, b)!=DataType.UNKNOWN_BYTE) tripR++;
                    if(p1a.getBase(taxa1, b)==(byte)p1a.getMinorAllele(b)) {
                        minorRight++;
                        errorByQScore[1][currentQ]++;
                      //  System.out.println(b+" "+p1a.getPositionInLocus(b)+" "+p1a.getBase(taxa1, b));
                    }
                    rightBySiteCnt[b]++;
                } else {
                    wrong++;
                    errorByQScore[2][currentQ]++;
                    depthW+=p1a.getSequenceCount()-p1a.getSiteSummary(b).getNumberMissing();
                    if(p1a.getBase(tripID, b)!=DataType.UNKNOWN_BYTE) tripW++;
                    //System.out.println(b+" "+p1a.getPositionInLocus(b));
                    wrongBySiteCnt[b]++;
                   // System.out.println(b+" "+p1a.getPositionInLocus(b));
                }
            }
            System.out.printf("%s %d %d %d %d %s %s %d %d %d %n", p1a.getLocus(startSite).toString(), startSite, endSite, p1a.getPositionInLocus(startSite),
                    p1a.getPositionInLocus(endSite),
                    p1a.getTaxaName(taxa1),p1a.getTaxaName(taxa2),
                    right, wrong, theSTF.getIntElement(7, i));
            totalRight+=right;
            totalMinorRight+=minorRight;
            totalWrong+=wrong;
        }
        double errorRate=(double)totalWrong/(double)(totalRight+totalWrong);
        double errorMinorRate=(double)totalWrong/(double)(totalMinorRight+totalWrong);
 //       System.out.printf("TotalRight %d  TotalWrong %d ErrorRate %g %n", totalright, totalwrong, errorRate);
        System.out.println("File\tSites\tTaxa\tTotalRight\tTotalMinorRight\tTotalWrong\tErrorRate\tErrorMinorRate");
        System.out.printf("%s\t%d\t%d\t%d\t%d\t%d\t%g\t%g %n", testAlignmentFile, p1a.getSiteCount(),p1a.getSequenceCount(),
                totalRight, totalMinorRight, totalWrong, errorRate, errorMinorRate);
//        System.out.printf("Depth Right: %d  Wrong: %d %n",depthR/totalRight, depthW/totalWrong);
//        System.out.printf("Trip Right: %g  Wrong: %g %n",(double)tripR/(double)totalRight, (double)tripW/(double)totalWrong);

        for (int q= 0; q < errorByQScore[0].length; q++) {
            double errorRateQ=(double)errorByQScore[2][q]/(double)(errorByQScore[0][q]+errorByQScore[2][q]);
            System.out.printf("%d %d %d %d %g %n",q,errorByQScore[0][q],errorByQScore[1][q],errorByQScore[2][q],errorRateQ);
        }
//        System.out.println("Site\tPosition\tSNPDef\tIBDRight\tIBDwrong\tMAF\tExpErrort\tPOfError");
//        for (int i = 0; i < wrongBySiteCnt.length; i++) {
//            if(rightBySiteCnt[i]+wrongBySiteCnt[i]>20) {
//                double maf=p1a.getMinorAlleleFrequency(i);
//                double expError=maf*(1-maf);
//                double p=-1.0;
//                int totalTests=rightBySiteCnt[i]+wrongBySiteCnt[i];
//                BinomialDistributionImpl x = new BinomialDistributionImpl(totalTests, expError);
//                try {
//                    p = x.cumulativeProbability(wrongBySiteCnt[i]);
//                } catch (Exception e) {
//                    System.err.println("Error in the BinomialDistributionImpl");
//                }
//                String snpDef=(char)p1a.getMajorAllele(i)+"/"+p1a.getMinorAllele(i);
//                System.out.printf("%d %d %s %d %d %g %g %g %n", i, p1a.getPositionInLocus(i), snpDef, rightBySiteCnt[i], wrongBySiteCnt[i], maf, expError, p);
//            }
//        }

    }


    public static int[] createSiteList(String infile, int siteIndex) {
        ArrayList<Integer> sites=new ArrayList();
        try{
            int count=0;
            BufferedReader fileIn = new BufferedReader(new FileReader(infile), 100000);
            fileIn.readLine();
            String sl;
            while ((sl=fileIn.readLine())!=null) {
                sites.add(Integer.parseInt(sl.split("\\s")[siteIndex]));
                count++;
            }
            fileIn.close();
        } catch (IOException e) {
            System.out.println("Error with createSiteList:"+e);
            e.printStackTrace();
        }
        int[] result=new int[sites.size()];
        for (int i = 0; i < sites.size(); i++) {
            result[i] = sites.get(i);
        }
        System.out.printf("File:%s Sites:%d %n",infile,result.length);
        return result;
    }

    public static void compareWithFile(String infile, int siteIndex, int[] sites, boolean keepSites, String outfile) {
        Arrays.sort(sites);
        int count=0;
        try{
            BufferedReader fileIn = new BufferedReader(new FileReader(infile), 100000);
            BufferedWriter fileOut = new BufferedWriter(new FileWriter(outfile), 100000);
            String sl=fileIn.readLine();
            fileOut.write(sl+"\n");
            while ((sl=fileIn.readLine())!=null) {
                int currSite=(int)Math.round(Double.parseDouble(sl.split("\\s")[siteIndex]));
                int hit=Arrays.binarySearch(sites, currSite);
 //               if(sl.contains("-")) continue;
                if(keepSites&&(hit>-1)) {fileOut.write(sl+"\n"); count++;}
                else if((keepSites==false)&&(hit<0)) {fileOut.write(sl+"\n"); count++;}
            }
            fileIn.close();
            fileOut.close();
        } catch (IOException e) {
             System.out.println("Error with createSiteList:"+e);
            e.printStackTrace();
        }
        System.out.printf("File:%s Sites:%d %n",outfile,count);
    }
    /**
     * 
     * @param infile1
     * @param siteIndex1
     * @param infile2
     * @param siteIndex2
     * @param outfile
     */
    public static void mergeFile(String infile1, int siteIndex1, String infile2, int siteIndex2, String outfile) {
        int count1=0, count2=0, countOut=0;
        try{
            BufferedReader fileIn = new BufferedReader(new FileReader(infile1), 100000);
            BufferedReader fileIn2 = new BufferedReader(new FileReader(infile2), 100000);
            BufferedWriter fileOut = new BufferedWriter(new FileWriter(outfile), 100000);
            String sl=fileIn.readLine();
            String sl2=fileIn2.readLine();
            fileOut.write((sl+"\t"+sl2).replaceAll("\\s", ",")+"\n");
            boolean moreToGo=true;
            sl=fileIn.readLine();
            sl2=fileIn2.readLine();
            while ((moreToGo)&&(sl!=null)&&(sl2!=null)) {
                String[] sp1=sl.split("\\s");
                String[] sp2=sl2.split("\\s");
                int currSite1=(int)Math.round(Double.parseDouble(sp1[siteIndex1]));
                int currSite2=(int)Math.round(Double.parseDouble(sp2[siteIndex2]));
                if(currSite1<currSite2) {sl=fileIn.readLine(); count1++;}
                else if(currSite1>currSite2) {sl2=fileIn2.readLine(); count2++;}
                else {//they must be equal
                    //if(sp1[siteIndex1+3].equals(sp2[siteIndex2+1]))  {
                        String s=(sl+sl2).replaceAll("\\s", ",")+"\n";
                        fileOut.write(s);
                        countOut++;
                    //}
                    sl=fileIn.readLine();count1++;
                    sl2=fileIn2.readLine();count2++;
                    
                }
            }
            fileIn.close();
            fileIn2.close();
            fileOut.close();
        } catch (IOException e) {
             System.out.println("Error with createSiteList:"+e);
            e.printStackTrace();
        }
        System.out.printf("File:%s Sites:%d %n",infile1,count1);
        System.out.printf("File:%s Sites:%d %n",infile2,count2);
        System.out.printf("File:%s Sites:%d %n",outfile,countOut);

    }

    public static void main(String[] args) throws Exception {
        String ibdfile="/Users/edbuckler/SolexaAnal/HapMapV2/IBD/chrAll.ibd.txt";
    //    String ibdfile="/Users/edbuckler/SolexaAnal/HapMapV2/IBD/AgInbredchr10.ibd.txt";
     //   String infile="/Users/edbuckler/SolexaAnal/HapMapV2/anchorLDHap/chr10.genotypes.log2_h90_f50.good.good.hmp.txt";
//        String infile="/Users/edbuckler/SolexaAnal/HapMapV2/anchor/chr10.genotypes.log2_h90_f50.good.hmp.txt";
       // String infile="/Users/edbuckler/SolexaAnal/HapMapV2/hp1/chr10.genotypes.log2_h90_f50.hmp.txt";
        String infile="/Users/edbuckler/SolexaAnal/HapMapV2/test/chr10.08232010.genotypes.log0cx_h0_f12.hmp.txt";
        String outfile="/Users/edbuckler/SolexaAnal/HapMapV2/test/fusion101223.txt";
        //String parentDirectoryS="/Users/edbuckler/SolexaAnal/HapMapV2/test/comps/";
 //       String parentDirectoryS="/Users/edbuckler/SolexaAnal/HapMapV2/jermpipe/SNPS201010/fusion/";
        String parentDirectoryS="/Users/edbuckler/SolexaAnal/HapMapV2/jermpipe/SNPS201010/CSHL-ALL/HP1/";

//        mergeFile("/Users/edbuckler/SolexaAnal/HapMapV2/test/keep/chr10.08232010.genotypes.callrep.txt",1,
//                "/Users/edbuckler/SolexaAnal/HapMapV2/test/chr10.08232010.genotypes.log0cx_h0_f12.ier.txt",1,
//                "/Users/edbuckler/SolexaAnal/HapMapV2/test/merge3.txt");
//        System.exit(0);
//        mergeFile("/Users/edbuckler/SolexaAnal/HapMapV2/test/cenchr10.08232010.genotypes.callrep.txt",1,
//                "/Users/edbuckler/SolexaAnal/AGPv2/sub16mer.txt",1,
//                "/Users/edbuckler/SolexaAnal/HapMapV2/test/merge5.txt");
//        System.exit(0);


//        int[] bgiSites=createSiteList("/Users/edbuckler/SolexaAnal/HapMapV2/calls/chr10.genotypes.bgi.txt",1);
//        compareWithFile("/Users/edbuckler/SolexaAnal/HapMapV2/test/comps/chr10.genotypes.log2_h90_f50.good.good.hmp.txt",3,bgiSites,false,outfile);
//        int[] bgiSites=createSiteList("/Users/edbuckler/SolexaAnal/HapMapV2/test/anchorNIndel_NBGI.hmp.txt",3);
//        int[] bgiSites=createSiteList("/Users/edbuckler/SolexaAnal/HapMapV2/test/comps/chr10.08232010.genotypes.log2cx_h99_sing_f12.hmp.txt",3);
//        compareWithFile("/Users/edbuckler/SolexaAnal/HapMapV2/test/comps/chr10.genotypes.bgi.log2cx_h99_sing_f12.hmp.txt",3,bgiSites,true,outfile);
//        System.exit(0);
      //  int[] bgiSites=createSiteList("/Users/edbuckler/SolexaAnal/HapMapV2/test/anchorNIndel_NBGI.hmp.txt",3);
//        int[] bgiSites=createSiteList("/Users/edbuckler/SolexaAnal/HapMapV2/test/PositionFailingIBDInAnchor.txt",1);
//        compareWithFile("/Users/edbuckler/SolexaAnal/HapMapV2/test/keep/chr10.08232010.genotypes.callrep.txt",1,bgiSites,true,outfile);
//        System.exit(0);
        Alignment intersectAlignment=null;
 //       intersectAlignment=(Pack1Alignment) ImportUtils.createPack1AlignmentFromFile("/Users/edbuckler/SolexaAnal/HapMapV2/test/comps/chr10.genotypes.bgi.log0cx_h0_f12.hmp.txt", "", "");
        //intersectAlignment=(Pack1Alignment) ImportUtils.createPack1AlignmentFromFile("/Users/edbuckler/SolexaAnal/HapMapV2/test/comps/chr10.08232010.genotypes.log0cx_h0_f12.hmp.txt", "", "");
//        Alignment anchorAlign=(Pack1Alignment) ImportUtils.createPack1AlignmentFromFile("/Users/edbuckler/SolexaAnal/HapMapV2/test/comps/chr10.genotypes.log2_h90_f50.good.good.hmp.txt", "", "");
//        try{
//            int count=0;
//            BufferedWriter fileOut = new BufferedWriter(new FileWriter(outfile), 100000);
//            for (int i = 0; i < anchorAlign.getSiteCount(); i++) {
//                if(intersectAlignment.getSiteOfPhysicalPosition(anchorAlign.getPositionInLocus(i), null)<0) {
//                    fileOut.write(anchorAlign.getPositionInLocus(i)+" "+(char)anchorAlign.getMajorAllele(i)+" "+(char)anchorAlign.getMinorAllele(i)+"\n");
//                    count++;
//                }
//            }
//            fileOut.close();
//        } catch (IOException e) {
//
//        }
        File parentDirectory = new File(parentDirectoryS);
        File[] fileList = parentDirectory.listFiles();
        for(File fn: fileList) {
            if(!fn.getName().endsWith(".txt")) continue;
             AlignmentQualityTest hld3=new AlignmentQualityTest(ibdfile, fn.getPath(), intersectAlignment, outfile);
        }
     System.exit(0);
        AlignmentQualityTest hld3=new AlignmentQualityTest(ibdfile, infile, intersectAlignment, outfile);
//System.exit(0);
infile="/Users/edbuckler/SolexaAnal/HapMapV2/test/chr10.genotypes.bgi.log0cx_h0_f12.hmp.txt";
//        infile="/Users/edbuckler/SolexaAnal/HapMapV2/hp1/chr10.genotypes.log2_h90_f50.hmp.txt";
        hld3=new AlignmentQualityTest(ibdfile, infile, intersectAlignment, outfile);

        infile="/Users/edbuckler/SolexaAnal/HapMapV2/hp1/chr10.genotypes.log2ct_h90.hmp.txt";
        hld3=new AlignmentQualityTest(ibdfile, infile, intersectAlignment, outfile);

    }


}
