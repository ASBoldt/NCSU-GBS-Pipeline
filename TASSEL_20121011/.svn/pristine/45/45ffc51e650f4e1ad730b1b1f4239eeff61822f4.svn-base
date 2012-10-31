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
public class SNPCallerAgainstAnchorV2 {
    
    LDbyBitsV3 anchorLDBits;
    HaplotypeLengthDistributionV3 anchorHaplotypes;
    File infileSiteCalls=null, outfileGoodSites, outfileSiteAttributes, anchorFile;
    int minTaxaWithSNP=12, minHomozygousCnt=1;
    boolean extractSingleton=false;
    double minLogPValueSNP=1, minHomoProp=0.90, minAlleleQual=20, maximumHeterozygousRatio=2.001;
    int ldPermutation=1;
    boolean isUsingContigency=true;
    int[] ldHapRedirect, hapRedict;

    private static int maxReadValue = 100000;
    private ContigencyTable contingencyTable;
    private FisherExact fishersExact;
    IdGroup anchorIdGroup;
    int isBGI=0;  //0=CSHL, 1=BGI

    public SNPCallerAgainstAnchorV2(String infileSiteCallsName, String anchorFileName, String outfileGoodSitesName,
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
            int taxaCnt=taxaNames.length;
            String row = null;
            System.out.println("Processing: " + infileSiteCalls.getName());
            int rawCount = 0,
            countGreaterThanMinTaxaWithSNP = 0,
            countGreaterThanMinLogPValue = 0;
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
            fileOutSiteAttr.write("LogRegScore\t");
            fileOutSiteAttr.write("\n");
            int errorCnt=0;
            while ((row=fileIn.readLine())!=null) {
                try {
                    if ((!row.contains("Sample"))&&(!row.contains("line"))) {
                        SNPDistV2[] theSNPDist=new SNPDistV2[2];
                        theSNPDist[0] = new SNPDistV2(taxaCnt, row, contingencyTable, fishersExact,1);
                        if(theSNPDist[0].afreq[2]>0) theSNPDist[1] = new SNPDistV2(taxaCnt, row, contingencyTable, fishersExact,2);
                        for(SNPDistV2 tsd:theSNPDist) {
                            if(tsd==null) continue;
                            if (tsd.taxaNumWithReads < minTaxaWithSNP) continue;  //eliminate low freq sites, regions mostly only in reference
                            if(tsd.avgQual[0]<minAlleleQual) continue;  //removes low quality reads
                            if(tsd.avgQual[tsd.altAlleleNumber]<minAlleleQual) continue;  //removes low quality reads
                            int[] homoCnts=tsd.getHomozygousCounts();
                            double propHomozygous=((double)homoCnts[0]+(double)homoCnts[1])/(double)homoCnts[2];
                            int heterozygousLines=homoCnts[2]-homoCnts[0]-homoCnts[1];
                            int linesHomoForMinorAllele=(homoCnts[0]>homoCnts[1])?homoCnts[1]:homoCnts[0];
                            if(linesHomoForMinorAllele<minHomozygousCnt) continue;  //removes SNPs not homozygous in at least one inbred line
                            if(propHomozygous<minHomoProp) continue;  //most SNPs should be homozygous in over 90% of the lines, as these are inbred samples
                            double hetToMinorHomoRatio = (double)heterozygousLines/(double)linesHomoForMinorAllele;
                            if(hetToMinorHomoRatio>maximumHeterozygousRatio) continue;  //removes many of the real but paralogous variants
                            countGreaterThanMinTaxaWithSNP++;
                            tsd.scoreSNPFastX2();tsd.scoreSNP();
                            if(tsd.snpLogP<minLogPValueSNP) continue;
                            countGreaterThanMinLogPValue++;
                            //LD with surrounding, can code so only test homozygous, or can just test alternate allele
                            long[][] siteInBits=tsd.getAllelesInBits(true, ldHapRedirect, anchorIdGroup.getIdCount());
                            double[] minLDP=anchorLDBits.findNeighbors(siteInBits, (int)tsd.startPos, 200, 500, ldPermutation, LDbyBitsV3.AlleleMode.RefvAlt);
                            double qualPred=logisticPred(linesHomoForMinorAllele, propHomozygous, tsd.snpLogP,
                                 Math.log10(minLDP[0]), isBGI);
                            tsd.maxMLScore=qualPred;
                            fileOut.write(tsd.toStringProp() + "\n");
                            fileOutSiteAttr.write(tsd.toStringAttributes());
                            fileOutSiteAttr.write(homoCnts[0]+"\t"+homoCnts[1]+"\t"+propHomozygous+"\t");
                            fileOutSiteAttr.write(tsd.minorAlleleProp[0]+"\t");
                            fileOutSiteAttr.write(minLDP[0]+"\t"+minLDP[1]+"\t"+minLDP[2]+"\t"+minLDP[3]+"\t");
                            fileOutSiteAttr.write(qualPred+"\t");
                            fileOutSiteAttr.write("\n");
                        }
                    }
                    rawCount++;
                    if (rawCount % 10000 == 0) {
                        currentTime = System.currentTimeMillis();
                        System.out.print("File:" + infileSiteCalls.getName() + ": ");
                        System.out.println("SNP number:" + rawCount +
                                " count over taxa&homozygous ratios+:" + countGreaterThanMinTaxaWithSNP +
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

    public double getMaximumHeterozygousRatio() {
        return maximumHeterozygousRatio;
    }

    public void setMaximumHeterozygousRatio(double maximumHeterozygousRatio) {
        this.maximumHeterozygousRatio = maximumHeterozygousRatio;
    }



    public double logisticPred(int HomoAltCnt, double propHomo, double SNP_LOG_P, double logLD, int bgi) {
       double[][] coeff={{-17.0773128,17.8947227,0.4257686,-0.1174901,0.00},
            {4.8842601,-3.7770660,0.4065829,0.00, 0.3534554}};
       int cb;
       if(HomoAltCnt>1) {cb=0;}
       else if(HomoAltCnt==1) {cb=1;}
       else{return -1;}
       double z=coeff[cb][0];  //intercept
       z+=(propHomo*coeff[cb][1])+(SNP_LOG_P*coeff[cb][2])+(logLD*coeff[cb][3])+(bgi*coeff[cb][4]);
       double pred=1.0/(1.0+Math.exp(-z));
       return pred;
    }

/**
 * High freq model
 * HomoAltCnt > 1 & hetaltrat<=2.001 & PropHomo > 0.9 & Alt1Qual > 20
 * Coefficients:
              Estimate Std. Error z value Pr(>|z|)
(Intercept) -17.077313   0.234471  -72.83   <2e-16
PropHomo     17.894723   0.246806   72.50   <2e-16
SNP_LOG_P     0.425769   0.009048   47.06   <2e-16
logLD        -0.117490   0.001743  -67.40   <2e-16
 *
 * Singleton model from R
 * HomoAltCnt == 1 & hetaltrat<=2.001 & PropHomo > 0.9 & Alt1Qual > 20
 * Coefficients:
            Estimate Std. Error z value Pr(>|z|)
(Intercept)  4.88426    1.22411   3.990  6.6e-05
PropHomo    -3.77707    1.26162  -2.994  0.00275
SNP_LOG_P    0.40658    0.01398  29.078  < 2e-16
bgi          0.35346    0.02323  15.217  < 2e-16
 */

}
