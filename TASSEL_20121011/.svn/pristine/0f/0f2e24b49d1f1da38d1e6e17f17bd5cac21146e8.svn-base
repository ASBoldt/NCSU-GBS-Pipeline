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
import java.util.Arrays;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.Pack1Alignment;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.statistics.ContigencyTable;
import net.maizegenetics.pal.statistics.FisherExact;

/**
 *
 * @author edbuckler
 */
public class SNPCallerAgainstAnchor {
    
    LDbyBitsV3 anchorLDBits;
    HaplotypeLengthDistributionV3 anchorHaplotypes;
    File infileSiteCalls=null, outfileGoodSites, outfileSiteAttributes, anchorFile;
    int minTaxaWithSNP=0, minHomozygousCnt=1;
    boolean extractSingleton=false;
    double minLogPValueSNP=1, minHomoProp=0.00, minAlleleQual=20;
    int ldPermutation=1;
    boolean isUsingContigency=true;
 //   int minPos=47000, maxPos=570000;
//    int minPos=48000000, maxPos=57000000;
    int[] ldHapRedirect, hapRedict;

    private static int maxReadValue = 100000;
    private ContigencyTable contingencyTable;
    private FisherExact fishersExact;
    IdGroup anchorIdGroup;
    int isBGI=0;  //0=CSHL, 1=BGI

    public SNPCallerAgainstAnchor(String infileSiteCallsName, String anchorFileName, String outfileGoodSitesName,
            String outfileSiteAttributesName) {
        if(!infileSiteCallsName.endsWith(".txt")) return;
        if(infileSiteCallsName.toLowerCase().contains("bgi")) {isBGI=1;}
        else if(infileSiteCallsName.toLowerCase().contains("cshl")) {isBGI=0;}
        else {System.out.println("Error all input file need to be tagged with bgi or cshl.");  System.exit(0);}
        infileSiteCalls=new File(infileSiteCallsName);
        outfileGoodSites=new File(outfileGoodSitesName);
        outfileSiteAttributes=new File(outfileSiteAttributesName);
        anchorFile=new File(anchorFileName);
        contingencyTable = new ContigencyTable(maxReadValue);
        fishersExact = new FisherExact(maxReadValue);
        Pack1Alignment p1a = (Pack1Alignment) ImportUtils.createPack1AlignmentFromFile(anchorFileName, "", "");
        anchorIdGroup=p1a.getIdGroup();
        System.out.println("Alignment Loaded:"+anchorFileName);
        System.out.printf("Sites %d Taxa %d %n", p1a.getSiteCount(),p1a.getSequenceCount());
        anchorLDBits=new LDbyBitsV3(p1a,1);
       // anchorHaplotypes=new HaplotypeLengthDistributionV3(p1a);

    }

    void run() {
        if(infileSiteCalls==null) return;
        System.out.println("Starting idSNPGenotypeFileToQSNPFile on:"+infileSiteCalls.getName());
        double[] readByTaxaPerPresentSite=null;
        System.out.println(Arrays.toString(readByTaxaPerPresentSite));
        try {
            BufferedReader fileIn = new BufferedReader(new FileReader(infileSiteCalls), 100000);
            BufferedWriter fileOut = new BufferedWriter(new FileWriter(outfileGoodSites), 100000);
            BufferedWriter fileOutSiteAttr = new BufferedWriter(new FileWriter(outfileSiteAttributes), 100000);
            String taxaHeaderString=fileIn.readLine().replaceAll("_", "");
            taxaHeaderString=taxaHeaderString.split(":\\s")[1];
            String[] taxaNames=taxaHeaderString.split("[\\s]+");
//            taxaHeaderString=taxaHeaderString.split(":\t")[1];
//            String[] taxaNames=taxaHeaderString.split("\t");
            int taxaCnt=taxaNames.length;
            String row = null;
            System.out.println("Processing: " + infileSiteCalls.getName());
            int rawCount = 0,
                    countGreaterThanMinTaxaWithSNP = 0,
                    countGreaterThanMinLogPValue = 0;
            int tripIndex=Arrays.binarySearch(taxaNames, "TDD39103");
            ldHapRedirect = new int[taxaCnt];
            for (int t = 0; t < taxaCnt; t++) {
                ldHapRedirect[t] = anchorIdGroup.whichIdNumber(taxaNames[t]);
            }
            long currentTime;
            long starttime = System.currentTimeMillis();            
            fileOut.write(SNPDistV2.toStringPropHeader(taxaNames) + "\n");
            fileOutSiteAttr.write(SNPDistV2.toStringAttributesHeader());
            fileOutSiteAttr.write("HomoRefCnt\tHomoAltCnt\tPropHomo\t");
            fileOutSiteAttr.write("RefProp\tLD_P_RefAlt\tLD_N_RefAlt\tLD_PermP_RefAlt\tLD_BestPos_RefAlt\t");
         //   fileOutSiteAttr.write("LD_P_AltRest\tLD_N_AltRest\tLD_PermP_AltRest\tLD_P_MissPres\tLD_N_MissPres\tLD_PermP_MissPres\t");
         //   fileOutSiteAttr.write("ImpCorr\tImpTotal\tImpLongCorr\tImpLongTotal\tImpMinorCorr\tImpMinorTotal\tAvgLengthMinorMatch\tAvgLengthMinorMismatch\t");
         //   fileOutSiteAttr.write("Tripsacum\tComplexity\t");
            fileOutSiteAttr.write("LogRegScore\t");
            fileOutSiteAttr.write("\n");
            int errorCnt=0;
            while ((row=fileIn.readLine())!=null) {
                try {
                    if ((!row.contains("Sample"))&&(!row.contains("line"))) {
                        SNPDistV2[] theSNPDist=new SNPDistV2[2];
                        theSNPDist[0] = new SNPDistV2(taxaCnt, row, contingencyTable, fishersExact,1);
                        if(theSNPDist[0].afreq[2]>0) theSNPDist[1] = new SNPDistV2(taxaCnt, row, contingencyTable, fishersExact,2);
//                if((theSNPDist[0].startPos<minPos)) {fileIn.skip(100000000); fileIn.readLine(); System.out.println(fileIn.readLine()); continue;}
//
//                if((theSNPDist[0].startPos>maxPos)) break;
                        for(SNPDistV2 tsd:theSNPDist) {
                            if(tsd==null) continue;
                            if (tsd.taxaNumWithReads < minTaxaWithSNP) continue;
                            if(tsd.avgQual[0]<minAlleleQual) continue;
                            if(tsd.avgQual[tsd.altAlleleNumber]<minAlleleQual) continue;
                            int[] homoCnts=tsd.getHomozygousCounts();
                            double propHomozygous=((double)homoCnts[0]+(double)homoCnts[1])/(double)homoCnts[2];
                 if(extractSingleton&&((homoCnts[1]!=1)&&(homoCnts[0]!=1))) continue;
                            if((homoCnts[0]<minHomozygousCnt)||(homoCnts[1]<minHomozygousCnt)) continue;
                            if(propHomozygous<minHomoProp) continue;
                            double altProp=(double)tsd.afreq[tsd.altAlleleNumber]/(double)tsd.totalFreq;
                            countGreaterThanMinTaxaWithSNP++;
                            //B73 present
                            int refHetero=0, refPresent=0;
                            if(Double.isNaN(tsd.minorAlleleProp[0])) {refPresent=1;}
                            else if (tsd.minorAlleleProp[0] > 0) {refHetero=1;}
                            tsd.scoreSNPFastX2();tsd.scoreSNP();
//                            if(isUsingContigency) {tsd.scoreSNPX2ThenContigency();}  //this the the more rigorous approach
//                            else  {tsd.scoreSNPFastX2(); tsd.snpLogP=tsd.snpX2P;}  //fast but not as good
                          if(tsd.snpLogP<minLogPValueSNP) continue;
                            countGreaterThanMinLogPValue++;
                            

                            //LD with surrounding, can code so only test homozygous, or can just test alternate allele
                     //       long[][] siteInBits=tsd.getAllelesInBits(true);
                            long[][] siteInBits=tsd.getAllelesInBits(true, ldHapRedirect, anchorIdGroup.getIdCount());
                            double[] minLDP=anchorLDBits.findNeighbors(siteInBits, (int)tsd.startPos, 200, 500, ldPermutation, LDbyBitsV3.AlleleMode.RefvAlt);
                            double qualPred=this.logisticPred(propHomozygous, altProp, refHetero, refPresent,
                                tsd.totalFreq, tsd.afreq[tsd.altAlleleNumber], Math.log10(minLDP[0]), tsd.snpLogP,
                                homoCnts[1], homoCnts[0], isBGI);
                            tsd.maxMLScore=qualPred;
                            fileOut.write(tsd.toStringProp() + "\n");
                            fileOutSiteAttr.write(tsd.toStringAttributes());
                            fileOutSiteAttr.write(homoCnts[0]+"\t"+homoCnts[1]+"\t"+propHomozygous+"\t");
                            fileOutSiteAttr.write(tsd.minorAlleleProp[0]+"\t");
                            fileOutSiteAttr.write(minLDP[0]+"\t"+minLDP[1]+"\t"+minLDP[2]+"\t"+minLDP[3]+"\t");
                            fileOutSiteAttr.write(qualPred+"\t");

//                            siteInBits=tsd.getAllelesInBits(false);
//                            minLDP=anchorLDBits.findNeighbors(siteInBits, (int)tsd.startPos, 200, 500, ldPermutation, LDbyBitsV3.AlleleMode.AltVRest);
//                            fileOutSiteAttr.write(minLDP[0]+"\t"+minLDP[1]+"\t"+minLDP[2]+"\t");
//                            minLDP=anchorLDBits.findNeighbors(siteInBits, (int)tsd.startPos, 200, 500, ldPermutation, LDbyBitsV3.AlleleMode.MissVPres);
//                            fileOutSiteAttr.write(minLDP[0]+"\t"+minLDP[1]+"\t"+minLDP[2]+"\t");
                            //Presense in special lines - teosinte, tripsacum
                            //haplotype test
//                            byte[] siteInBytes=tsd.getAllelesInBytes(true,0.25,0.75);
//                            byte[] siteInBytes=tsd.getAllelesInBytes(true,0.25,0.75, ldHapRedirect, anchorIdGroup.getIdCount());
//                            boolean tryInc=anchorHaplotypes.incrementToThisPosition((int)tsd.startPos);
//                            int[] imputeResults;
//                            if(anchorHaplotypes.currSplitSite==(int)tsd.startPos) {imputeResults=anchorHaplotypes.estimateImputationAccuracy(siteInBytes,true);}
//                            else {imputeResults=anchorHaplotypes.estimateImputationAccuracy(siteInBytes,false);}
//                            fileOutSiteAttr.write(imputeResults[0]+"\t"+imputeResults[1]+"\t"+imputeResults[2]+"\t"+imputeResults[3]+"\t"+imputeResults[4]+"\t"+imputeResults[5]+"\t");
//                            fileOutSiteAttr.write(imputeResults[6]+"\t"+imputeResults[7]+"\t");
//                            int hasTripsacum=(Double.isNaN(tsd.minorAlleleProp[tripIndex]))?0:1;
//                            fileOutSiteAttr.write(hasTripsacum+"\t"+0+"\t");
                            fileOutSiteAttr.write("\n");
                        }
                    }
                    rawCount++;
                    if (rawCount % 1000 == 0) {
                        currentTime = System.currentTimeMillis();
                        System.out.print("File:" + infileSiteCalls.getName() + ": ");
                        System.out.println("SNP number:" + rawCount +
                                " count over taxa threshold+:" + countGreaterThanMinTaxaWithSNP +
                                " count over logp threshold:" + countGreaterThanMinLogPValue +
                                " Time:" + ((currentTime - starttime) / 1000)+
                                " Error:"+errorCnt);
                        starttime = currentTime;
                    }
                } catch (Exception e) {
                    errorCnt++;
                    System.out.println("ERROR: " + e + "\t RowNumber:"  + rawCount+" ErrorCnt:"+errorCnt);
                    e.printStackTrace();
                    System.out.println("ERROR row text: " +row );
                }
            }
            fileIn.close();
            fileOut.close();
            fileOutSiteAttr.close();
            System.out.println("SNP number:" + rawCount +
                                " count over taxa threshold+:" + countGreaterThanMinTaxaWithSNP +
                                " count over logp threshold:" + countGreaterThanMinLogPValue);
            System.out.println(infileSiteCalls.getName() + " finished!");
        } catch (Exception e) {
            System.err.println("File IO in NucleotideDiversityProcessing: " + e);
            e.printStackTrace();
        }
    }

    public int getMaxReadValue() {
        return maxReadValue;
    }

    public void setMaxReadValue(int maxReadValue) {
        this.maxReadValue = maxReadValue;
    }

    public double getMinAlleleQual() {
        return minAlleleQual;
    }

    public void setMinAlleleQual(double minAlleleQual) {
        this.minAlleleQual = minAlleleQual;
    }

    public double getMinHomoProp() {
        return minHomoProp;
    }

    public void setMinHomoProp(double minHomoProp) {
        this.minHomoProp = minHomoProp;
    }

    public double getMinLogPValueSNP() {
        return minLogPValueSNP;
    }

    public void setMinLogPValueSNP(double minLogPValueSNP) {
        this.minLogPValueSNP = minLogPValueSNP;
    }

    public int getMinTaxaWithSNP() {
        return minTaxaWithSNP;
    }

    public void setMinTaxaWithSNP(int minTaxaWithSNP) {
        this.minTaxaWithSNP = minTaxaWithSNP;
    }


    public double logisticPred(double propHomo, double altProp, int refHetero, int referPresent,
            int TotalCnt, int Alt1Cnt, double logLD, double SNP_LOG_P,
            double HomoAltCnt, double HomoRefCnt, int bgi) {
       double[][] coeff={{-10.16742052,9.780953018,3.861959977,-1.456553448,-0.127789127,0.005314515,0.001787567,-0.303015999,0.246268072,-0.154366045,-0.018506626,1.662238612},
            {-6.599434773,7.704389552,5.507571399,-1.462107886,-0.314379786,0.001770701,0.007625059,-0.470990594,0.282982691,-0.503013911,-0.006300654,0.681185121},
            {-0.490579332,2.663613433,0,-1.245669382,-0.244655892,0,0,0,0.3076597,-0.009343453,0,0.338999997}};
       int cb;
       if(HomoAltCnt>5) {cb=0;}
       else if(HomoAltCnt>1) {cb=1;}
       else{cb=2;}
       double z=coeff[cb][0];
       z+=(propHomo*coeff[cb][1])+(altProp*coeff[cb][2])+(refHetero*coeff[cb][3])+(referPresent*coeff[cb][4]);
       z+=(TotalCnt*coeff[cb][5])+(Alt1Cnt*coeff[cb][6])+(logLD*coeff[cb][7])+(SNP_LOG_P*coeff[cb][8]);
       z+=(HomoAltCnt*coeff[cb][9])+(HomoRefCnt*coeff[cb][10])+(bgi*coeff[cb][11]);
       double pred=1.0/(1.0+Math.exp(-z));
       return pred;
    }



}
