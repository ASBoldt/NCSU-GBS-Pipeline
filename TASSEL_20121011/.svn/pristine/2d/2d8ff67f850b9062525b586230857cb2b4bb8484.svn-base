/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.solexa;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Arrays;
import net.maizegenetics.pal.statistics.ContigencyTable;
import net.maizegenetics.pal.statistics.FisherExact;

/**
 *
 * @author edbuckler
 */
public class SNPCaller {
    private static int maxReadValue = 100000;

    /**
     *
     * @param infileName
     * @param outfileName
     * @param minTaxaWithSNP
     * @param minLogPValueSNP
     * @param minHomoProp
     * @param minAlleleQual
     * @param isFreq
     * @param isUsingContigency
     */
    public static void idSNPGenotypeFileToQSNPFile(String infileName, String depthFileName, String outfileName,
            int minTaxaWithSNP, double minLogPValueSNP, double minHomoProp, double minAlleleQual,
            int minHomoLines, boolean isUsingContigency, boolean extractSingleton) {
        System.out.println("Starting idSNPGenotypeFileToQSNPFile on:"+infileName);
        ContigencyTable contingencyTable = new ContigencyTable(maxReadValue);
        FisherExact fishersExact = new FisherExact(maxReadValue);
        double[] readByTaxaPerPresentSite=null;
        if(depthFileName!=null) {readByTaxaPerPresentSite=getAvgCoverageFromFile(depthFileName);}
        System.out.println(Arrays.toString(readByTaxaPerPresentSite));
        try {
            if(!infileName.endsWith(".txt")) return;
            BufferedReader fileIn = new BufferedReader(new FileReader(infileName), 100000);
            BufferedWriter fileOut = new BufferedWriter(new FileWriter(outfileName), 100000);
            String taxaHeaderString=fileIn.readLine();
            taxaHeaderString=taxaHeaderString.split(":\\s")[1];
            String[] taxaNames=taxaHeaderString.split("[\\s]+");
            int taxaCnt=taxaNames.length;
            String row = null;
            System.out.println("Processing: " + infileName);
            int rawCount = 0,
                    countGreaterThanMinTaxaWithSNP = 0,
                    countGreaterThanMinLogPValue = 0;
            long currentTime;
            long starttime = System.currentTimeMillis();            
            fileOut.write(SNPDistV2.toStringPropHeader(taxaNames) + "\n");
            int errorCnt=0;
            while ((row=fileIn.readLine())!=null) {
                try {
                    if ((!row.contains("Sample"))&&(!row.contains("line"))) {
                        SNPDistV2[] theSNPDist=new SNPDistV2[2];
                        theSNPDist[0] = new SNPDistV2(taxaCnt, row, contingencyTable, fishersExact,1);
                        if(theSNPDist[0].afreq[2]>0) theSNPDist[1] = new SNPDistV2(taxaCnt, row, contingencyTable, fishersExact,2);
                        for(SNPDistV2 tsd:theSNPDist) {
                            if(tsd==null) continue;
                            if (tsd.taxaNumWithReads < minTaxaWithSNP) continue;
                            if(tsd.avgQual[0]<minAlleleQual) continue;
                            if(tsd.avgQual[tsd.altAlleleNumber]<minAlleleQual) continue;
                            int[] homoCnts=tsd.getHomozygousCounts();
                            double propHomozygous=((double)homoCnts[0]+(double)homoCnts[1])/(double)homoCnts[2];
                            if((homoCnts[0]<minHomoLines)||(homoCnts[1]<minHomoLines)) continue;
                            if(propHomozygous<minHomoProp) continue;
               if(extractSingleton&&(homoCnts[1]!=1)&&(homoCnts[0]!=1)) continue;
                            countGreaterThanMinTaxaWithSNP++;
                            //median depth
                            //this is a good idea, but the 36bp reads are very focused and increase the variance
                       //     double medianDepth=tsd.getMedianRelativeDepthOfPresent(readByTaxaPerPresentSite);
                       //     if((medianDepth>=2)||(medianDepth<=0.5)) continue;
                            //B73 present
                            if(Double.isNaN(tsd.minorAlleleProp[0])||(tsd.minorAlleleProp[0]!=0)) continue;
                            if(isUsingContigency) {tsd.scoreSNPX2ThenContigency();}  //this the the more rigorous approach
                            else  {tsd.scoreSNPFastX2(); tsd.snpLogP=tsd.snpX2P;}  //fast but not as good
                            if(tsd.snpLogP<minLogPValueSNP) continue;
                            countGreaterThanMinLogPValue++;
                            fileOut.write(tsd.toStringProp() + "\n");
                        }
                    }
                    rawCount++;
                    if (rawCount % 10000 == 0) {
                        currentTime = System.currentTimeMillis();
                        System.out.print("File:" + infileName + ": ");
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
            System.out.println("SNP number:" + rawCount +
                                " count over taxa threshold+:" + countGreaterThanMinTaxaWithSNP +
                                " count over logp threshold:" + countGreaterThanMinLogPValue);
            System.out.println(infileName + " finished!");
        } catch (Exception e) {
            System.err.println("File IO in NucleotideDiversityProcessing: " + e);
            e.printStackTrace();
        }
    }

     public static double[] idSNPGenotypeFileAvgCoverage(String infileName, String outfileName) {
         System.out.println("Starting idSNPGenotypeFileAvgCoverage on:"+infileName);
         double[] readByTaxaPerSite=null;
         double[] readByTaxaPerPresentSite=null;
         int rowNumber=SimpleTextFile.countRows(infileName);
         short[][] depthByTaxaSite;
         try {
            if(!infileName.endsWith(".txt")) return null;
            System.out.println(infileName+" rowCount:"+rowNumber);
            BufferedReader fileIn = new BufferedReader(new FileReader(infileName), 100000);
            BufferedWriter fileOut = new BufferedWriter(new FileWriter(outfileName), 100000);
            String taxaHeaderString=fileIn.readLine();
            taxaHeaderString=taxaHeaderString.split(":\t")[1];
            String[] taxaNames=taxaHeaderString.split("\t");
            int taxaCnt=taxaNames.length;
            long[] readCntByTaxa=new long[taxaCnt];
            int[] presentSiteCntByTaxa=new int[taxaCnt];
            depthByTaxaSite=new short[taxaCnt][rowNumber];
            String row = null;
            System.out.println("Processing: " + infileName);
            int rawCount = 0;
            long currentTime;
            long starttime = System.currentTimeMillis();
            for(String tn: taxaNames) {fileOut.write(tn+"\t");} fileOut.write("\n");
            int errorCnt=0;
            while ((row=fileIn.readLine())!=null) {
                try {
                    if ((!row.contains("Sample"))&&(!row.contains("line"))) {
                        SNPDistV2 theSNPDist = new SNPDistV2(taxaCnt, row, null, null,1);
                        for (int t = 0; t < taxaCnt; t++) {
                            depthByTaxaSite[t][rawCount]=(short)theSNPDist.totalReads4Taxa[t];
                            readCntByTaxa[t]+=depthByTaxaSite[t][rawCount];
                            if(depthByTaxaSite[t][rawCount]>0) presentSiteCntByTaxa[t]++;
                        }
                    }
                    rawCount++;
                    if (rawCount % 100000 == 0) {
                        currentTime = System.currentTimeMillis();
                        System.out.print("File:" + infileName + ": ");
                        System.out.println("SNP number:" + rawCount +
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
            readByTaxaPerSite=new double[taxaCnt];
            readByTaxaPerPresentSite=new double[taxaCnt];
            for (int t = 0; t < taxaCnt; t++) {
                readByTaxaPerSite[t]=(double)readCntByTaxa[t]/(double)rawCount;
                readByTaxaPerPresentSite[t]=(double)readCntByTaxa[t]/(double)presentSiteCntByTaxa[t];
            }
            for(long tn: readCntByTaxa) {fileOut.write(tn+"\t");} fileOut.write("\n");
            for(long tn: presentSiteCntByTaxa) {fileOut.write(tn+"\t");} fileOut.write("\n");
            for(double tn: readByTaxaPerSite) {fileOut.write(tn+"\t");} fileOut.write("\n");
            for(double tn: readByTaxaPerPresentSite) {fileOut.write(tn+"\t");} fileOut.write("\n");
            for (int t = 0; t < presentSiteCntByTaxa.length; t++) {
                Arrays.sort(depthByTaxaSite[t]);
                int index=depthByTaxaSite[t].length/2;
                short s = depthByTaxaSite[t][index];
                fileOut.write(s+"\t");
             }
            fileOut.write("\n");
            for (int t = 0; t < presentSiteCntByTaxa.length; t++) {
                int index=((2*depthByTaxaSite[t].length)-presentSiteCntByTaxa[t])/2;
                short s = depthByTaxaSite[t][index];
                fileOut.write(s+"\t");
             }
            fileOut.write("\n");
            fileIn.close();
            fileOut.close();
            System.out.println("SNP number:" + rawCount +
                                " count over taxa threshold+:" + 0 +
                                " count over logp threshold:" + 0 );
            System.out.println(infileName + " finished!");
        } catch (Exception e) {
            System.err.println("File IO in idSNPGenotypeFileAvgCoverage: " + e);
            e.printStackTrace();
        }
        return readByTaxaPerPresentSite;
    }

     public static double[] getAvgCoverageFromFile(String infileName) {
         double[] readByTaxaPerPresentSite=null;
         try{
            if(!infileName.endsWith(".depth.txt")) return null;
            BufferedReader fileIn = new BufferedReader(new FileReader(infileName), 100000);
            fileIn.readLine();fileIn.readLine();fileIn.readLine();fileIn.readLine();fileIn.readLine();fileIn.readLine();
            String[] data=fileIn.readLine().split("\\s");
            readByTaxaPerPresentSite=new double[data.length];
            for (int i = 0; i < data.length; i++) {
                 readByTaxaPerPresentSite[i]=Double.parseDouble(data[i]);
             }
         } catch (Exception e) {
            System.err.println("File IO in getAvgCoverageFromFile: " + e);
            e.printStackTrace();
            return null;
        }
         return readByTaxaPerPresentSite;
     }

}
