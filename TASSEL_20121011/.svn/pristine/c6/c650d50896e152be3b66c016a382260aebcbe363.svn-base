package net.maizegenetics.genome.solexa;

/**
 * User: Ed and Jason
 * Date: Feb 17, 2009
 */
import net.maizegenetics.baseplugins.FilterAlignmentPlugin;
import net.maizegenetics.baseplugins.SequenceDiversityPlugin;
import net.maizegenetics.pal.alignment.AllelePositionBLOBUtils;
import net.maizegenetics.pal.alignment.AnnotatedAlignmentUtils;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.popgen.LinkageDisequilibrium;
import net.maizegenetics.pal.report.SimpleTableReport;
import net.maizegenetics.pal.report.TableReportUtils;
import net.maizegenetics.pal.statistics.ContigencyTable;
import net.maizegenetics.pal.statistics.FisherExact;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;

import java.io.*;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.SiteSummary;

public class EdDivProfilerV2 implements Runnable {

    private static char NEW_UNKNOWN_CHARACTER = 'N';
//    private static String argAlt="-D E:/SolexaAnal/AGPSNPCover/dummy/ E:/SolexaAnal/test/"; //make coverage SNPs on desktop
//    private static String argAlt="-ERR E:/SolexaAnal/AGPSNPCover/SNPIndel2hit/ E:/SolexaAnal/test/"; //make error distributions on desktop
//    private static String argAlt="-CS C:/EdStuff/Solexa/AGPLog2ML/combined1_13/ C:/EdStuff/Solexa/test/ C:/EdStuff/Solexa/AGPCoverage/maxmap1/ C:/EdStuff/Solexa/AGPdividedByChr/"; //make coverage SNPs on laptop
//    private static String argAlt="-CS E:/SolexaAnal/AGPLog2ML/combined1_13/ E:/SolexaAnal/test/ E:/SolexaAnal/AGPCoverage/maxmap1/"; //make coverage SNPs on desktop
//    private static String argAlt="-CS E:/SolexaAnal/AGPLog2ML/combined1_13/ E:/SolexaAnal/test/ E:/SolexaAnal/AGPCoverage/maxmap1/ E:/SolexaAnal/maize_pseudo/AGPdividedByChr/"; //make coverage SNPs on desktop with reference and alternate alleles
//    private static String argAlt="-LD E:/SolexaAnal/AGPSNPCover/SNPIndel1hitThres5poly/ E:/SolexaAnal/test/"; //make coverage SNPs on desktop with reference and alternate alleles
//  private static String argAlt="-SNP E:/SolexaAnal/AGPSNPCover/SNPIndel1hitThres5poly/ E:/SolexaAnal/test/ E:/SolexaAnal/maize_pseudo/AGPdividedByChr/"; //make coverage SNPs on desktop with reference and alternate alleles
    // private static String argAlt="-HP1 C:/EdStuff/Solexa/AGPLog2ML/combined1_13/ C:/EdStuff/Solexa/test/"; //make hapmap format from machine learn file
//    private static String argAlt = "-HP1 E:/SolexaAnal/AGPLog2ML/combined1_13/ E:/SolexaAnal/test/"; //make hapmap format from machine learn file
//    private static String argAlt = "-ALB1 E:/SolexaAnal/AGPLog2ML/combined1_13HP/ E:/SolexaAnal/test/"; //convert hapmap format to GDPDM ALLELE BLOBs

//    private static String argAlt="-CS E:/SolexaAnal/AGPLog2ML/combined2_13/ E:/SolexaAnal/test/ E:/SolexaAnal/AGPCoverage/maxmap2/"; //make coverage SNPs on desktop
//    private static String argAlt="-S E:/SolexaAnal/AGPallcalls/combined1_13/ E:/SolexaAnal/test/"; //make  SNPs calls on desktop
//    private static String argAlt="-S /Users/edbuckler/SolexaAnal/HapMapV2/calls/ /Users/edbuckler/SolexaAnal/HapMapV2/test/"; //make  SNPs calls on desktop
private static String argAlt="-HP1 /Users/edbuckler/SolexaAnal/HapMapV2/log2/ /Users/edbuckler/SolexaAnal/HapMapV2/test/"; //make  SNPs calls on desktop
    //    private static String argAlt="-LD /Users/edbuckler/SolexaAnal/HapMapV2/hp1/ /Users/edbuckler/SolexaAnal/HapMapV2/test/"; //make  SNPs calls on desktop
    //private static String argAlt="-S E:/HapMapV2/calls/ E:/HapMapV2/test/"; //make  SNPs calls on desktop
//    private static String argAlt="-HP1 E:/HapMapV2/log2/ E:/HapMapV2/test/"; //make hapmap format from machine learn file

//    private static String argAlt="-D C:/EdStuff/Solexa/AGPSNPCover/SNP1hitThres5/ C:/EdStuff/Solexa/test/"; //make diversity calls on laptop
//    private static String argAlt="-D E:/SolexaAnal/AGPSNPCover/SNPIndel1hitThres5/ E:/SolexaAnal/test/"; //make diversity calls on desktop
//        private static String argAlt="-OUT E:/SolexaAnal/AGPSNPCover/SNPIndel1hitThres5Ref/ E:/SolexaAnal/test/ E:/SolexaAnal/AGPSorghum/maize2sorghum_20090601/"; //make diversity calls on desktop
//    private static String argAlt="-CON C:/EdStuff/Solexa/AGPSNPCover/SNPIndel1hitThres5/ C:/EdStuff/Solexa/test/"; //make errors calls on laptop
    //   private String parentInDirectory, parentInCoverDirectory, parentOutDirectory;
    private static int maxTaxa = 54,  maxReadValue = 1000000,  minTaxaWithSNP = 12;
//    private static double minLogPValueSNP = 2,  minMachineLearnScore = 0.85, minHomoProp=0.9;
    private static double minLogPValueSNP = 2,  minMachineLearnScore = 0.85, minHomoProp=0.95, minAlleleQual=10;
 //   private static double minLogPValueSNP = -1,  minMachineLearnScore = -1;
    private boolean codeIndelsAsSNPs = false;
//    private static String   parentInDirectory  ="D:/Solexa/Maize0902/log2/";
//    private static String   parentInDirectory  ="E:/SolexaAnal/AGPLog2ML/";
//    private static String   parentInCoverDirectory  ="E:/SolexaAnal/AGPCoverage/";
//    private static String   parentInDirectory  ="E:/SolexaAnal/AGPallcalls/combined1_13/";
//      private static String   parentInDirectory  ="E:/SolexaAnal/AGPSNPCover/SNP1hitThres5/";
//    private static String   parentInDirectory  ="E:/SolexaAnal/AGPSNPCover/SNPIndel1hitThres5/";
//    private static String   parentInDirectory  ="E:/SolexaAnal/AGPLog2ML/dummy/";
//    private static String   parentInCoverDirectory  ="E:/SolexaAnal/AGPCoverage/dummy/";
//    private static String   parentInDirectory  ="E:/SolexaAnal/AGPLog2ML/combined1_13/";
//    private static String   parentInCoverDirectory  ="E:/SolexaAnal/AGPCoverage/maxmap1/";
//    private static String   parentOutDirectory ="E:/SolexaAnal/test/";
    //   private static String   coverageDataInfile  ="E:/SolexaAnal/AGPCoverage/chr8combined.txt";
//    private static String  validDataInfile="E:/SolexaAnal/AGPValid/PANZEA_genotype090410.txt";
    private static String validDataInfile = "C:/EdStuff/Solexa/AGPValid/PANZEA_genotype090410.txt";
//    private static String   parentInDirectory  ="C:/EdStuff/Solexa/coverage/";
    //   private static String   parentInDirectory  ="D:/Solexa/Maize0902/coverage/";
//    private static String   parentInDirectory  ="C:/EdStuff/Solexa/AGPallcalls/";
//    private static String   parentInDirectory  ="C:/EdStuff/Solexa/AGPSNPCover/SNP1hitThres5/";
//    private static String   parentInDirectory  ="C:/EdStuff/Solexa/AGPLog2ML/";
//    private static String   coverageDataInfile  ="C:/EdStuff/Solexa/AGPCoverage/chr8combined.txt";
//    private static String   parentOutDirectory ="C:/EdStuff/Solexa/test/";
    private ContigencyTable contingencyTable;
    private FisherExact fishersExact;
    String[] taxaOrder = {"B73", "B97", "CML103", "CML228", "CML247", "CML277", "CML322",
        "CML333", "CML52", "CML52R", "CML69", "HP301", "IL14H", "KI11", "KI3", "KY21", "M162W",
        "M37W", "MO17", "MO18W", "MS71", "NC350", "NC358", "OH43", "OH7B", "P39", "TX303", "TZI8"};
    String[] taxaOrderName = {"B73:temperate:mays:mays:Zea", "B97:temperate:mays:mays:Zea", "CML103:tropical:mays:mays:Zea",
        "CML228:tropical:mays:mays:Zea", "CML247:tropical:mays:mays:Zea", "CML277:tropical:mays:mays:Zea", "CML322:tropical:mays:mays:Zea",
        "CML333:tropical:mays:mays:Zea", "CML52:BC1:mays:mays:Zea", "CML52R:tropical:mays:mays:Zea",
        "CML69:tropical:mays:mays:Zea", "HP301:popcorn:mays:mays:Zea", "IL14H:sweet:mays:mays:Zea",
        "KI11:tropical:mays:mays:Zea", "KI3:tropical:mays:mays:Zea", "KY21:temperate:mays:mays:Zea", "M162W:temperate:mays:mays:Zea",
        "M37W:tropical:mays:mays:Zea", "MO17:temperate:mays:mays:Zea", "MO18W:tropical:mays:mays:Zea",
        "MS71:temperate:mays:mays:Zea", "NC350:tropical:mays:mays:Zea", "NC358:tropical:mays:mays:Zea",
        "OH43:temperate:mays:mays:Zea", "OH7B:temperate:mays:mays:Zea", "P39:sweet:mays:mays:Zea",
        "TX303:tropical:mays:mays:Zea", "TZI8:tropical:mays:mays:Zea"};
    private String mode = "",  infileForThread = "",  infile2ForThread = "",  outfileForThread = "",  refGenomeFile = "";
    private int currentSNP = 100000000;

    /**=CONSTRUCTOR-instantiates a Contingency table as well as
     * the Fishers Exact Test Constructor*/
    /*    public EdDivProfilerV2() {
    this("", "", "");
    }

    public EdDivProfilerV2(String infile, String outfile) {
    this(infile, "", outfile);
    }

    public EdDivProfilerV2(String infile, String infile2, String outfile) {
    this.contingencyTable = new ContigencyTable(maxReadValue);
    this.fishersExact     = new FisherExact(maxReadValue);
    this.infileForThread=infile;
    this.infile2ForThread=infile2;
    this.outfileForThread=outfile;
    }
     */
    public EdDivProfilerV2(String[] args, String infile) {
        this.contingencyTable = new ContigencyTable(maxReadValue);
        this.fishersExact = new FisherExact(maxReadValue);
        this.mode = args[0];
        this.infileForThread = args[1] + infile;
        this.infile2ForThread = "";
        this.outfileForThread = "";
        String stem = infile;
        String[] namesp = stem.split("\\.");
        if (mode.equals("-CS")) {
            infile2ForThread = args[3] + namesp[0] + ".coverage.1.combined.txt";
            outfileForThread = args[2] + stem.replaceFirst(".txt", ".wcover.txt");
            if (args.length == 5) {
                refGenomeFile = args[4] + "chr" + namesp[0] + ".txt";
            }
        } else if (mode.equals("-D")) {
            outfileForThread = args[2] + stem.replaceFirst(".txt", ".div.txt");   //for SNP calling file
        //   outfileForThread=args[2]+"Diversity090508a.txt";   //for SNP calling file
        } else if (mode.equals("-LD")) {
            outfileForThread = args[2] + stem.replaceFirst(".txt", ".LD.txt");   //for SNP calling file
            outfileForThread = args[2] + "LD_N20_F2_DistMin500_090606.txt";   //for LD calling file
        //           outfileForThread=args[2]+"LD_cum200_200_090608.txt";   //for LD calling file
        } else if (mode.equals("-SNP")) {
            outfileForThread = args[2] + stem.replaceFirst(".txt", ".SNP.txt");   //for SNP calling file
            outfileForThread = args[2] + "IlluminaSNPs_090616.txt";   //for LD calling file
            refGenomeFile = args[3] + "chr" + namesp[0] + ".txt";
        //           outfileForThread=args[2]+"LD_cum200_200_090608.txt";   //for LD calling file
        } else if (mode.equals("-S")) {
            outfileForThread = args[2] + stem.replaceFirst(".txt", ".log2ml.txt");   //for SNP calling file
        } else if (mode.equals("-HP1")) {
            outfileForThread = args[2] + stem.replaceFirst(".txt", ".HP1.txt");   //for SNP calling file
        } else if (mode.equals("-ERR")) {
            outfileForThread = args[2] + "NoThresdiversityB97RecErrors090512.err.txt";   //for SNP calling file
        //   outfileForThread=args[2]+stem.replaceFirst(".txt",".err.txt");   //for SNP calling file
        } else if (mode.equals("-CON")) {
            outfileForThread = args[2] + "diversityConcordance090512.err.txt";   //for determining concordance
        //   outfileForThread=args[2]+stem.replaceFirst(".txt",".err.txt");   //for determining concordance
        } else if (mode.equals("-OUT")) {
            infile2ForThread = args[3] + "coverage2sorghum_chr" + namesp[0];
            outfileForThread = args[2] + stem.replaceFirst(".txt", ".wOut.txt");
        }
    }

    /**=Read in a file containing snp calls across 27 founders
     * output bac accession, position in bac, polymorphism type, allele definition,
     * snp call, reference and alt reads in each taxon, total taxa with reads,
     * tropical with reads, temperate with reads, snp log from contingency table,
     * threshold for call, pi, variance of pi, fst, logP fst*/
    private void NucleotideDiversityProcessing(String infileName, String outfileName, boolean isFreq, boolean isTraining, int trainingLine) {
        try {
//             FPCMap theFPCMap= new FPCMap(new File(FPCmap));
            BufferedReader fileIn = new BufferedReader(new FileReader(infileName), 100000);
            BufferedWriter fileOut = new BufferedWriter(new FileWriter(outfileName), 100000);
//            BufferedWriter fileTrainOut = null;
            String taxaHeaderString=fileIn.readLine();
            taxaHeaderString=taxaHeaderString.split(":\t")[1];
            String[] taxaNames=taxaHeaderString.split("\t");
//            if (isTraining) {
//                fileTrainOut = new BufferedWriter(new FileWriter(outfileName + ".TrainML" + trainingLine), 1000000);
//            }
            String row = null;
            System.out.println("Processing: " + infileName);
            int rawCount = 0,
                    countGreaterThanMinTaxaWithSNP = 0,
                    countGreaterThanMinLogPValue = 0;
            long currentTime;
            long starttime = System.currentTimeMillis();
            if (isFreq) {//fileOut.write(SNPDistV2.toStringPropQualHeader(1)+"\n");
            } else {
                fileOut.write(SNPDistV2.toStringPropHeader(taxaNames) + "\n");
            }
//            if (isTraining) {
//                fileTrainOut.write(SNPDistV2.toStringMachineHeader(trainingLine));
//                fileTrainOut.write("ExpIBD\n");
//            }
            
            while //(true) {
                ((fileIn.ready()==true)||(fileIn.ready()==false)) {
                try {
                    row = fileIn.readLine();
                    if ((!row.contains("Sample"))&&(!row.contains("line"))) {
                        if (row.contains("SNP") == false) {
                            //                           System.out.println(row);
                        }
                        SNPDistV2 theSNPDist = new SNPDistV2(maxTaxa, row, contingencyTable, fishersExact);
                        theSNPDist.taxaOrder=taxaNames;

                        if (theSNPDist.taxaNumWithReads > minTaxaWithSNP) {
                            countGreaterThanMinTaxaWithSNP++;
                          //  theSNPDist.scoreSNPX2ThenContigency();
                                theSNPDist.scoreSNPFastX2(); theSNPDist.snpLogP=theSNPDist.snpX2P;  //fast but not as good
                            //theSNPDist.scoreMaxMachineLearnScore();
                            int[] homoCnts=theSNPDist.getHomozygousCounts();
                            double propHomozygous=((double)homoCnts[0]+(double)homoCnts[1])/(double)homoCnts[2];
                            if ((theSNPDist.snpLogP >= minLogPValueSNP) &&
                                    (propHomozygous>=minHomoProp)&&
                                    (homoCnts[0]>=2)&&(homoCnts[1]>=2)&&
                                    (theSNPDist.avgQual[0]>minAlleleQual)&&
                                    (theSNPDist.avgQual[1]>minAlleleQual)
                                    )
                                //    ||  (theSNPDist.maxMLScore > minMachineLearnScore))
                                {
                               //     System.out.println(theSNPDist.snpLogP+" "+theSNPDist.snpX2P);
                                countGreaterThanMinLogPValue++;
                                theSNPDist.calcMaxThresholdPiFst();
                                if (isFreq) {
                                    fileOut.write(theSNPDist.toStringNuc() + "\n");
                                    System.out.println("SNP number:");
                                } else {
                                    fileOut.write(theSNPDist.toStringProp() + "\n");
                                }
                            //                      fileOut.write(chrinfo[0]+"\t"+chrinfo[1]+"\t"+ibd+"\n");
                            }
                        }
//                        if (isTraining) {
//                            String ls = theSNPDist.toStringMachine(trainingLine, false);
//                            if (ls != null) {
//                                fileTrainOut.write(ls);
//                                fileTrainOut.write(getCML52RILidentityWB73(theSNPDist.startPos) + "\n");
//                            }
//                        }
                    }
                    rawCount++;
                    if (rawCount % 10000 == 0) {
                        currentTime = System.currentTimeMillis();
                        System.out.print("File:" + infileName + ": ");
                        System.out.println("SNP number:" + rawCount +
                                " count over taxa threshold+:" + countGreaterThanMinTaxaWithSNP +
                                " count over logp threshold:" + countGreaterThanMinLogPValue +
                                " Time:" + ((currentTime - starttime) / 1000));
                        starttime = currentTime;
                    }
                } catch (Exception e) {
                    System.err.println("ERROR: " + e + "\t" + row + "\t" + rawCount);
                }
            }
            fileIn.close();
            fileOut.close();
          //  fileTrainOut.close();
            System.out.println("SNP number:" + rawCount +
                                " count over taxa threshold+:" + countGreaterThanMinTaxaWithSNP +
                                " count over logp threshold:" + countGreaterThanMinLogPValue);
            System.out.println(infileName + " finished!");
        } catch (Exception e) {
            System.err.println("File IO in NucleotideDiversityProcessing: " + e);
        }
    }

    /**
     * This creates the Human Genetics community HapMap format from our quantiative SNP call files
     * 
     * @param infileName 
     * @param outfileName
     */
    private void QSNPCallsToHapMapFormat(String infileName, String outfileName, boolean filterByThreshold, boolean includeIndels) {
        int totalSNPIndelSites = 0;
        DecimalFormat chrForm = new DecimalFormat("00");
        DecimalFormat siteForm = new DecimalFormat("00000000");
        int numSeqs = 54;
        int nanCnt=0, refCnt=0, altCnt=0, hetCnt=0;
        try {
            BufferedReader fileIn = new BufferedReader(new FileReader(infileName), 1000000);
            BufferedWriter fileOut = new BufferedWriter(new FileWriter(outfileName, false), 100000);
            String[] header=fileIn.readLine().split("\t");
            taxaOrderName=new String[numSeqs];
            for (int i = 0; i < numSeqs; i++) {
                taxaOrderName[i] = header[i+10];
                if(taxaOrderName[i].startsWith("TI")) {
                    taxaOrderName[i]=taxaOrderName[i]+":TEO";
                }  else {
                    taxaOrderName[i]=taxaOrderName[i]+":MZ";
                }
            }
            fileOut.write("rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode\t");
            for (String t : taxaOrderName) {
                fileOut.write(t + "\t");
            }
            fileOut.write("\n");
            System.out.println("Setting states based on " + infileName);
            String sl="";

            while ((sl=fileIn.readLine())!=null) {
                //String sl = fileIn.readLine();
                SNPDistV2 theSNPDist = new SNPDistV2(numSeqs, sl, null, null);
                if ((filterByThreshold) && (theSNPDist.maxThreshold != 0.5)) {
                    continue;
                }
                if ((includeIndels == false) && (!theSNPDist.varType.startsWith("SNP"))) {
                    continue;
                }
                StringBuilder sb = new StringBuilder();
                sb.append("PZE");
                sb.append(chrForm.format(theSNPDist.chromosome));
                sb.append(siteForm.format(theSNPDist.startPos));
                sb.append("\t");
               // String[] states = theSNPDist.alleleDef.split("/");
                sb.append(theSNPDist.alleleDef[0] + "/" + theSNPDist.alleleDef[1] + "\t");
                sb.append(theSNPDist.chromosome + "\t");
                sb.append((int) theSNPDist.startPos + "\t");
                sb.append("+\tAGPv1\tMaizeDiv\tSBS\tMHPv1\tNAMfnd\tQC+\t");
                char altBase = 'S', refBase = 'S', hetBase='S';
                refBase = theSNPDist.alleleDef[0].charAt(0);
                altBase = theSNPDist.alleleDef[1].charAt(0);
                byte[] stateByte={(byte)refBase, (byte)altBase};
                hetBase=(char)AllelePositionBLOBUtils.getBaseFromHalfByte(AllelePositionBLOBUtils.getHalfByteFromSNPValue(stateByte));
                if(theSNPDist.varType.equals("IDP")) {
                    //ref
                    altBase=theSNPDist.alleleDef[1].charAt(0);
                    refBase=(altBase=='-')?'+':'-';
                   // System.out.println(theSNPDist.alleleDef);
                    hetBase='0';
                }
  //              System.out.println(refBase+" "+altBase+" "+hetBase);
                boolean isVariant = false;  //there are rare situations where a site may have the correct threshold, but no single line passes
                //the test below.  We track this below.
                for (double map : theSNPDist.minorAlleleProp) {
                    if (Double.isNaN(map)) {
                        sb.append(NEW_UNKNOWN_CHARACTER + "\t");
                        nanCnt++;
                    } else if (map > .99) {
                        sb.append(altBase + "\t");
                        isVariant = true;
                        altCnt++;
                    } else if (map <0.01 ) {
                        sb.append(refBase + "\t");
                        refCnt++;
                    } else {
                        sb.append(hetBase + "\t");
                        hetCnt++;
                    }
                }
//                for(double map: theSNPDist.minorAlleleProp) {
//                    if(Double.isNaN(map)) {sb.append(NEW_UNKNOWN_CHARACTER+""+NEW_UNKNOWN_CHARACTER+"\t");}
//                    else if(map>(theSNPDist.maxThreshold*1.5)) {sb.append(altBase+""+altBase+"\t"); isVariant=true;}
//                    else if(map>(theSNPDist.maxThreshold*0.5)) {sb.append(refBase+""+altBase+"\t");}
//                    else {sb.append(refBase+""+refBase+"\t");}
//                }
                if (isVariant) {
                    fileOut.write(sb.toString() + "\n");
                }
                totalSNPIndelSites++;
                if (totalSNPIndelSites % 10000 == 0) {
                    System.out.println("SNP output:"+totalSNPIndelSites+" NaN:"+nanCnt+" ref:"+refCnt+" alt:"+altCnt+" het:"+hetCnt);
                }
            }
            fileOut.flush();
            fileOut.close();
            fileIn.close();
            System.out.println("Number of sites in SNPIndel file" + totalSNPIndelSites);
        } catch (Exception e) {
            System.err.println("File IO in setStatesBasedOnCoverageFile: " + e);
        }

    }

    private void genotypingSNPs(String infileName, String outfileName) {
        try {
            BufferedReader fileIn = new BufferedReader(new FileReader(infileName));
            BufferedWriter fileOut = new BufferedWriter(new FileWriter(outfileName));
            String row = null;
            System.out.println("Processing: " + infileName);
            SNPDistV2 lastSNP = null, currSNP = null, nextSNP = null;
            int rawCount = 0;
            String currBAC = "";
            long currentTime;
            long starttime = System.currentTimeMillis();
            fileOut.write("BAC_ACCESSION\tPOSITION_IN_BAC\tPOLYMORPHISM_TYPE\tALLELE_DEFINITION\t");
            for (int i = 0; i < SNPDistV2.taxaOrder.length; i++) {
                fileOut.write(SNPDistV2.taxaOrder[i] + "\t");
            }
            fileOut.write("TAXA_W/_READS\tTROPICAL_W/_READS\tTEMPERATE_W/_READS\tSNP_LOG_P\t" +
                    "MAX_THRESHOLD\tPI\tPI_VAR\tFST\tFST_LOG_P\n");

            while (fileIn.ready()) {
                try {
                    row = fileIn.readLine();

                    if (!(row.contains("Sample") || (row.contains("ACCESSION")))) {
                        SNPDistV2 theSNPDist = new SNPDistV2(maxTaxa, row, contingencyTable, fishersExact);
                        lastSNP = currSNP;
                        currSNP = nextSNP;
                        nextSNP = theSNPDist;
                        rawCount++;
                    }
                    if (lastSNP != null) {
                        double ldist = 1000, rdist = 1000;
                        //        if(lastSNP.bacName.equals(currSNP.bacName)) {ldist=Math.abs(currSNP.bacPos-lastSNP.bacPos);}
                        //      if(lastSNP.bacName.equals(currSNP.bacName)) {rdist=Math.abs(currSNP.bacPos-nextSNP.bacPos);}
                        if (//(ldist>10)&&(rdist>10)&&       //ensures away from other SNPs
                                //       (currSNP.varType.equals("snp"))&& //ensures only SNPs
                                (currSNP.maxThreshold == 0.5) && //ensures biallelic
                                //   (currSNP.pi>0.2)&&      //ensures reasonable frequency
                                //      (currSNP.minorAlleleProp[0]<0.05)//&&
                                (currSNP.taxaNumWithReads > 19)) {  //ensures agrees with reference
                            fileOut.write(currSNP.toStringProp() + "\n");
                        }
                    }
                    if (rawCount % 10000 == 0) {
                        currentTime = System.currentTimeMillis();
                        System.out.println("SNP number:" + rawCount + " Time:" + ((currentTime - starttime) / 1000));
                        starttime = currentTime;
                    }
                } catch (Exception e) {
                    System.err.println("ERROR: " + e + "\t" + row + "\t" + rawCount);
                }
            }
            fileIn.close();
            fileOut.close();
            System.out.println(infileName + " finished!");
        } catch (Exception e) {
            System.err.println("File IO: " + e);
        }
    }

    private void evaluateConcordance(String infileName, String validInfileName, String outfileName) {
   /*     codeIndelsAsSNPs = false;
        System.out.println("evaluateAgreement processing: " + infileName);
        NextGenAlignmentWithCoverage ngai2 = new NextGenAlignmentWithCoverage(numberOfTaxa, infileName, codeIndelsAsSNPs);
        NextGenAlignmentIdeas validNAGI = new NextGenAlignmentIdeas(validInfileName);
        try {
            BufferedWriter fileOut = new BufferedWriter(new FileWriter(outfileName, true), 1000000);
            String row = null;
            System.out.println("Processing: " + infileName);
            fileOut.write("chr\tsite\tline\tb73Solexachar\tb73ValidChar\tlineSolexachar\tlineValidChar\tagree\n");
            int rawCount = 0,
                    countGreaterThanMinTaxaWithSNP = 0,
                    countInTrainingDataSet = 0;
            int currChr = ngai2.getChromosome(0);
            for (int i = 0; i < validNAGI.getSiteCount(); i++)
           {
                 if(validNAGI. currChr) continue;
            for (int i = 0; i < validNAGI.getSiteCount();i++) {
                 if(validNAGI.getChromosome(i)!=currChr) continue;
                 int site=ngai2.getSite(validNAGI.getChromosome(i), validNAGI.getChromosomePosition(i));
                 if(site>0) {
                    char b73Solexachar = ngai2.getDataChar(0, site);
                    char b73ValidChar = validNAGI.getDataChar(0, i);
                    for (int j = 0; j < ngai2.getSequenceCount(); j++) {
                        char lineSolexachar = ngai2.getDataChar(j, site);
                        char lineValidChar = validNAGI.getDataChar(j, i);
                        String outVals = b73Solexachar + "\t" + b73ValidChar + "\t" + lineSolexachar + "\t" + lineValidChar + "\t";
                        if (outVals.contains(DataType.UNKNOWN_CHARACTER + "")) {
                            continue;
                        }
                        boolean isCorrect = ((b73Solexachar == lineSolexachar) && (b73ValidChar == lineValidChar)) ||
                                ((b73Solexachar != lineSolexachar) && (b73ValidChar != lineValidChar));
                        fileOut.write(currChr + "\t" + site + "\t" + ngai2.getIdentifier(j).getName() + "\t" + outVals + "\t" + isCorrect + "\n");
                    }
                }
            }
            fileOut.close();
            System.out.println(infileName + " finished!");
        } catch (Exception e) {
            System.err.println("File IO: " + e);
        }
    */
    }

    private void createPositiveTraining(String infileName, String validInfileName, String outfileName) {
/*        NextGenAlignmentIdeas ngai = new NextGenAlignmentIdeas(validInfileName);
        try {
//             FPCMap theFPCMap= new FPCMap(new File(FPCmap));
            BufferedReader fileIn = new BufferedReader(new FileReader(infileName), 1000000);
            //           BufferedWriter fileOut  = new BufferedWriter(new FileWriter(outfileName),1000000);
            BufferedWriter fileTrainOut = new BufferedWriter(new FileWriter(outfileName + ".TrainALL"), 1000000);
            String row = null;
            System.out.println("Processing: " + infileName);
            int rawCount = 0,
                    countGreaterThanMinTaxaWithSNP = 0,
                    countInTrainingDataSet = 0;
            long currentTime;
            long starttime = System.currentTimeMillis();
            fileTrainOut.write(SNPDistV2.toStringMachineHeader(0) + "PanzeaB73\tPanzeaLine\tPanzeaPred\n");
            while (fileIn.ready()) {
                try {
                    row = fileIn.readLine();

                    if (!row.contains("Sample")) {
                        SNPDistV2 theSNPDist = new SNPDistV2(numberOfTaxa, row, contingencyTable, fishersExact);
                        theSNPDist.scoreSNPFastX2();
                        int site = ngai.getSite(theSNPDist.chromosome, theSNPDist.startPos);
                        if (site > 0) {
                            countInTrainingDataSet++;
                            char b73char = ngai.getDataChar(0, site);
                            for (int i = 0; i < ngai.getSequenceCount(); i++) {
                                if ((ngai.getDataChar(i, site) != DataType.UNKNOWN_CHARACTER) && (Double.isNaN(theSNPDist.minorAlleleProp[i])) == false) {
                                    String ls = theSNPDist.toStringMachine(i, true);
                                    int pred = (ngai.getDataChar(i, site) == b73char) ? 0 : 1;
                                    if (ls != null) {
                                        fileTrainOut.write(ls + b73char + "\t" + ngai.getDataChar(i, site) + "\t" + pred + "\n");
                                    }
                                    //                            System.out.println(ls+"\t"+b73char+"\t"+ngai.getData(i,site)+"\t"+pred);
                                    //                            System.out.println(theSNPDist.toStringProp());
                                    fileTrainOut.flush();
                                }
                            }
                        }

                    }
                    rawCount++;
                    if (rawCount % 10000 == 0) {
                        currentTime = System.currentTimeMillis();
                        System.out.println("SNP number:" + rawCount +
                                " count over taxa threshold+:" + countGreaterThanMinTaxaWithSNP +
                                " count over logp threshold:" + countInTrainingDataSet +
                                " Time:" + ((currentTime - starttime) / 1000));
                        starttime = currentTime;
                    }
                } catch (Exception e) {
                    System.err.println("ERROR: " + e + "\t" + row + "\t" + rawCount);
                }
            }
            fileIn.close();
            fileTrainOut.close();
            System.out.println(infileName + " finished!");
        } catch (Exception e) {
            System.err.println("File IO: " + e);
        }
 */
    }



    private void createFileWithOutgroup(String infileName, String outgroupInfileName, String outfileName) {
        try {
            BufferedReader fileIn = new BufferedReader(new FileReader(infileName), 1000000);
            BufferedReader fileOutGroup = new BufferedReader(new FileReader(outgroupInfileName), 1000000);
            BufferedWriter fileTrainOut = new BufferedWriter(new FileWriter(outfileName), 1000000);
            String row = null;
            System.out.println("Processing: " + infileName);
            int rawCount = 0,
                    countInTrainingDataSet = 0;
            //fileTrainOut.write(SNPDistV2.toStringMachineHeader(0)+"PanzeaB73\tPanzeaLine\tPanzeaPred\n");
            String[] snpLineParts, outGroupLineParts;
            fileTrainOut.write(fileIn.readLine() + "alignBasesNoGaps\tidWithRef\tidWithAltSNP\n");
            outGroupLineParts = fileOutGroup.readLine().split("\t");
            double snpStart = -1, snpEnd = -1, outStart = -1, outEnd = -1;
            outStart = Double.parseDouble(outGroupLineParts[1]);
            outEnd = Double.parseDouble(outGroupLineParts[2]);
            while (fileIn.ready()) {
                String snpLine = fileIn.readLine();
                snpLineParts = snpLine.split("\t");
                fileTrainOut.write(snpLine);
                snpStart = Double.parseDouble(snpLineParts[1]);
                snpEnd = Double.parseDouble(snpLineParts[2]);
                if (snpStart == outStart) {
                    if (snpEnd == outEnd) {
                        if (outGroupLineParts[1].equals("863543.00")) {
                            System.out.println("Waiting:\t" + outGroupLineParts[1] + "\t" + outGroupLineParts[2]);
                        }
                        int alignLength = outGroupLineParts[9].length();
                        int idLength = outGroupLineParts[9].replaceAll(" ", "").length();
                        int totalGaps = outGroupLineParts[7].length() + outGroupLineParts[8].length() -
                                outGroupLineParts[7].replaceAll("-", "").length() - outGroupLineParts[8].replaceAll("-", "").length();
                        int alignBasesNoGaps = alignLength - totalGaps;
                        int idWithRef = idLength;
                        int idWithAltSNP = (outGroupLineParts[8].equals(snpLineParts[5])) ? 1 : 0;
                        fileTrainOut.write("\t" + outGroupLineParts[8] + "\t" + alignBasesNoGaps + "\t" + idWithRef + "\t" + idWithAltSNP);
                        snpStart = snpEnd;
                    } else { //likely a fused SNP before
                    }
                //                System.out.println(snpLineParts[1]+"\t"+outGroupLineParts[1]+"\t"+outGroupLineParts[2]);
                } else {
                    fileTrainOut.write("\tN\t0\t0\t0");
                }
                if (fileOutGroup.ready() && (snpEnd >= outEnd)) {
                    outGroupLineParts = fileOutGroup.readLine().split("\t");
                    outStart = Double.parseDouble(outGroupLineParts[1]);
                    outEnd = Double.parseDouble(outGroupLineParts[2]);
                    //          System.out.println("Waiting:\t"+outGroupLineParts[1]+"\t"+outGroupLineParts[2]);
                    if (outGroupLineParts[1].equals("3170778.00")) {
                        System.out.println("Waiting:\t" + outGroupLineParts[1] + "\t" + outGroupLineParts[2]);
                    }
                }
                fileTrainOut.write("\n");
            }
            fileIn.close();
            fileTrainOut.close();
            System.out.println(infileName + " finished!");
        } catch (Exception e) {
            System.err.println("File IO in createFileWithOutgroup: " + e);
        }
    }

  
    private void characterizeErrors(String infileName, String outfileName) {
  /*      codeIndelsAsSNPs = false;
        System.out.println("diversityPipeline processing: " + infileName);
        if (infileName.contains(".txt") == false) {
            System.out.println("File type unrecognized ");
            return;
        }
        NextGenAlignmentWithCoverage ngai = new NextGenAlignmentWithCoverage(numberOfTaxa, infileName, codeIndelsAsSNPs);
        StringBuilder output = new StringBuilder();
        output.append("infileName\tngai.getSequenceCount()\tngai.getSiteCount()\tcurrSite\t");
        output.append("snpSite\tindelSite\t");
        output.append("snpSiteScoredInRef\tindelSitesScoredInRef\t");
        output.append("sitesScoredUnknown\tsitesScoredRef\tsitesScoredNonRef\tsitesScoredGap\n");
        //count B73 problems
        int sitesScoredRef = 0, sitesScoredNonRef = 0, sitesScoredUnknown = 0, sitesScoredGap = 0, snpSite = 0, indelSite = 0, snpSiteScoredInRef = 0,
                indelSitesScoredInRef = 0, refLine = 1;
        for (int i = 0; i < ngai.getSiteCount(); i++) {
            if ((false) && (i % 500000 == 0)) {
                output.append(infileName + "\t" + ngai.getSequenceCount() + "\t" + ngai.getSiteCount() + "\t" + i + "\t");
                output.append(snpSite + "\t" + indelSite + "\t");
                output.append(snpSiteScoredInRef + "\t" + indelSitesScoredInRef + "\t");
                output.append(sitesScoredUnknown + "\t" + sitesScoredRef + "\t" + sitesScoredNonRef + "\t" + sitesScoredGap + "\t");
                output.append("\n");
                //         System.out.println(output.toString());
                sitesScoredRef = sitesScoredNonRef = sitesScoredUnknown = sitesScoredGap = snpSite = indelSite = snpSiteScoredInRef = 0;
                indelSitesScoredInRef = 0;
            }
            boolean indel = false, snp = false;
            //classify the site either SNP or indel
            for (int j = 0; j < ngai.getSequenceCount(); j++) {
                char c = ngai.getDataChar(j, i);
                if (c == 'C') {
                    snpSite++;
                    if (ngai.getDataChar(refLine, i) != DataType.UNKNOWN_CHARACTER) {
                        snpSiteScoredInRef++;
                    }
                    break;
                } else if (c == DataType.PRIMARY_SUGGESTED_GAP_CHARACTER) {
                    indelSite++;
                    if (ngai.getDataChar(refLine, i) != DataType.UNKNOWN_CHARACTER) {
                        indelSitesScoredInRef++;
                    }
                    break;
                }
            }
            // if(indel) continue;
            char b = ngai.getDataChar(refLine, i);
            if (b == DataType.UNKNOWN_CHARACTER) {
                sitesScoredUnknown++;
            } else if (b == 'A') {
                sitesScoredRef++;
            } else if (b == 'C') {
                sitesScoredNonRef++;
            } else if (b == DataType.PRIMARY_SUGGESTED_GAP_CHARACTER) {
                sitesScoredGap++;
            }
        }
        output.append(infileName + "\t" + ngai.getSequenceCount() + "\t" + ngai.getSiteCount() + "\t" + ngai.getSiteCount() + "\t");
        output.append(snpSite + "\t" + indelSite + "\t");
        output.append(snpSiteScoredInRef + "\t" + indelSitesScoredInRef + "\t");
        output.append(sitesScoredUnknown + "\t" + sitesScoredRef + "\t" + sitesScoredNonRef + "\t" + sitesScoredGap + "\t");
        output.append("\n");
        System.out.println(output.toString());
        try {
            BufferedWriter fileOut = new BufferedWriter(new FileWriter(outfileName, true), 100000);
            fileOut.write(output.toString());
            fileOut.close();
        } catch (Exception e) {
            System.out.println("File output in diversityPipeline: " + e);
        }
*/
    }

    private void diversityPipeline(String infileName, String outfileName) {
        System.out.println("diversityPipeline processing: " + infileName);
        if (infileName.contains(".txt") == false) {
            System.out.println("File type unrecognized ");
            return;
        }
        Alignment ngai = ImportUtils.createPack1AlignmentFromFile(infileName, null, null);
      //  NextGenAlignmentWithCoveragev1_5 ngai = new NextGenAlignmentWithCoveragev1_5(numberOfTaxa, infileName, codeIndelsAsSNPs);
        //    System.out.println("ngai.getChromosomePosition(0)= "+ngai.getChromosomePosition(0));
        Datum dtAll = new Datum("Alignment", ngai, "Blah");
        int[] segSites = AnnotatedAlignmentUtils.getIncludedSitesBasedOnFreqIgnoreGapsMissing(ngai, 0.01, 12);
//    System.out.println("segSites.length:"+segSites.length);
        FilterAlignmentPlugin fap = new FilterAlignmentPlugin(null, false);
        fap.setEnd(ngai.getSiteCount());
        fap.setMinFreq(-1);
        fap.setMinCount(12);
        int window = 50000, step = 25000;
        ArrayList<SimpleTableReport> combinedrs = new ArrayList<SimpleTableReport>();
        for (int i = 0; i < ngai.getSiteCount() - window; i += step) {

            fap.setStart(i);
            fap.setEnd(i + window);
            System.out.println("Current i:" + i);
            Datum dtGood = fap.processDatum(dtAll, false);
            Alignment aa = (Alignment) dtGood.getData();
            //     System.out.println("aa sites:"+aa.getSiteCount());
            SequenceDiversityPlugin sdp = new SequenceDiversityPlugin(null, false);
            sdp.setEndSite(aa.getSiteCount() - 1);
//            sdp.setEndSite(window);
            DataSet rs = sdp.processDatum(dtGood);
            SimpleTableReport tr = (SimpleTableReport) rs.getData(1).getData();
            combinedrs.add(tr);
            System.out.println(TableReportUtils.toDelimitedString(tr, "\t"));
        }
        SimpleTableReport combinedReport = (SimpleTableReport) SimpleTableReport.getInstance(combinedrs);
        String output = TableReportUtils.toDelimitedString(combinedReport, "\t");
        System.out.println(output);
        try {
            BufferedWriter fileOut = new BufferedWriter(new FileWriter(outfileName, true), 100000);
            fileOut.write(output);
            fileOut.close();
        } catch (Exception e) {
            System.out.println("File output in diversityPipeline: " + e);
        }
    }


    private void ldPipeline(String infileName, String outfileName) {
        int minDistanceApart = 200;
        boolean permute=false;
        int bins = 5000;
        int window = 100;
        int minInbredComparisons=20;
        int minMinorInbred=2;
        int maxSite=200000;
        double[] sumR2Bins = new double[bins];
        int[] cntR2Bins = new int[bins];
        int binSize=100;

        System.out.println("LDPipeline processing: " + infileName);
        if (infileName.contains(".txt") == false) {
            System.out.println("File type unrecognized ");
            return;
        }
        System.out.println("Loading Alignment");
        Alignment ngai = ImportUtils.createPack1AlignmentFromFile(infileName, null, null);
        maxSite=ngai.getSiteCount();
        fishersExact=new FisherExact(ngai.getSequenceCount()+10);
        System.out.println("Finding sites with appropriate frequencies");
        boolean[] worthySites=new boolean[ngai.getSiteCount()];
        for (int i = 0; i < ngai.getSiteCount(); i++) {
            int majorCnt=countAlleles(ngai,i,ngai.getMajorAllele(i));
            int minorCnt=countAlleles(ngai,i,(byte)ngai.getMinorAllele(i));
            if((minorCnt>=minMinorInbred)&&((majorCnt+minorCnt)>=minInbredComparisons))
                {worthySites[i]=true;} else {worthySites[i]=false;}
        }
        double[][] bestLDpSite = new double[6][ngai.getSiteCount()];  //near r2, near P, >min distance r2, >min distance P
        Arrays.fill(bestLDpSite[0], 0);
        Arrays.fill(bestLDpSite[1],-1);
        Arrays.fill(bestLDpSite[2], 2);
        Arrays.fill(bestLDpSite[3], 0);
        Arrays.fill(bestLDpSite[4], -1);
        Arrays.fill(bestLDpSite[5], 2);
        System.out.println("Starting LD calculations");
        int conLostResNull=0, conLostResNaN=0, conLostWorthy=0;
        for (int i = 0; i < maxSite - window; i++) {
            if(worthySites[i]==false) {conLostWorthy+=window; continue;}
            for (int ts = i+1; ts < i + window; ts++) {
                int j=-1;
                if(!permute) {
                    j=ts;}
                else {
                    j=(int)Math.floor(Math.random()*(ngai.getSiteCount()-1));
                    if(Math.abs(i-j)<1000) continue;
                    }
                if(worthySites[j]==false) {conLostWorthy++; continue;}
                double[] results=LinkageDisequilibrium.getLDForSitePair(ngai, i, ngai, j,
                        minInbredComparisons, minMinorInbred, fishersExact);
                if(results==null) {conLostResNull++; continue;}
                if(Double.isNaN(results[1])) {conLostResNaN++; continue;}  //why does this happen
                bestLDpSite[0][i]++;
                bestLDpSite[0][j]++;
                if(results[3]<bestLDpSite[2][i]) {bestLDpSite[2][i]=results[3];bestLDpSite[1][i]=results[1];}
                if(results[3]<bestLDpSite[2][j]) {bestLDpSite[2][j]=results[3];bestLDpSite[1][j]=results[1];}
                int dist = Math.abs(ngai.getPositionInLocus(j) - ngai.getPositionInLocus(i));
                int cbin=dist/binSize;
                if(cbin>=bins) cbin=bins-1;
                sumR2Bins[cbin]+=results[1];  cntR2Bins[cbin]++;
                //if(results[3]<0.00001) System.out.println(i+","+j+":"+Arrays.toString(results)+" dist:"+dist);
                if(dist>minDistanceApart) {
                    bestLDpSite[3][i]++;
                    bestLDpSite[3][j]++;
                    if(results[3]<bestLDpSite[5][i]) {bestLDpSite[5][i]=results[3];bestLDpSite[4][i]=results[1];}
                    if(results[3]<bestLDpSite[5][j]) {bestLDpSite[5][j]=results[3];bestLDpSite[4][j]=results[1];}
                }

            }
            if(i%10000==0) {System.out.println("LD at site:"+i);
                System.out.println("conLostResNull:"+conLostResNull+" conLostResNaN:"+conLostResNaN+
                        " conLostWorthy:"+conLostWorthy);
            }
        }

        try {
            BufferedWriter fileOut = new BufferedWriter(new FileWriter(outfileName, false), 100000);
            BufferedWriter fileAvgLDOut = new BufferedWriter(new FileWriter(outfileName + "avg.txt", false), 100000);
            DecimalFormat f2 = new DecimalFormat("#0.00");
            fileOut.write("Chr\tPosition\tSite\tN\tr2\tP\tfar_N\tfar_r2\tfar_P\tCommonCnt\tLessCommCnt\n");
            for (int i = 0; i < maxSite; i++) {
                int majorCnt=countAlleles(ngai,i,ngai.getMajorAllele(i));
                int minorCnt=countAlleles(ngai,i,(byte)ngai.getMinorAllele(i));
                fileOut.write(ngai.getLocusName(i) + "\t" + f2.format(ngai.getPositionInLocus(i)) + "\t");
                fileOut.write(i+"\t");
                fileOut.write(bestLDpSite[0][i] + "\t" + bestLDpSite[1][i] + "\t"+ bestLDpSite[2][i] + "\t"+
                        bestLDpSite[3][i] + "\t" + bestLDpSite[4][i] + "\t" + bestLDpSite[5][i] + "\t");
                fileOut.write(majorCnt+"\t" + minorCnt +"\t");
                fileOut.write("\n");
            }
            fileAvgLDOut.write("Chr\tBinMin\tBinMax\tsumR2\tcntOfComps\tAvgR2\n");
            for (int i = 0; i < bins; i++) {
                if(cntR2Bins[i]==0) continue;
                fileAvgLDOut.write(ngai.getLocusName(0) + "\t" + (i*binSize) +"\t"+ ((i+1)*binSize-1) + "\t");
                fileAvgLDOut.write(sumR2Bins[i] + "\t" + cntR2Bins[i] + "\t");
                fileAvgLDOut.write(sumR2Bins[i]/cntR2Bins[i] + "\t");
                fileAvgLDOut.write("\n");
            }
            fileOut.close();
            fileAvgLDOut.close();
        } catch (Exception e) {
            System.out.println("File output in diversityPipeline: " + e);
        }
    }

    /**
     * Thoughts permute to determine null
     *
     * high speed LD - calculate the MAF - compare the folded MAF, dependins on fold on make specific contrasts
     *
     */

    private int countAlleles(Alignment a, int site, byte allele) {
        int cnt=0;
        for (int i = 0; i < a.getSequenceCount(); i++) {
            if(a.getBase(i, site)==allele) cnt++;
        }
        return cnt;
    }

    private void oldldPipeline(String infileName, String outfileName) {
        int minDistanceApart = 500;
        int bins = 5000;
        double[] sumR2Bins = new double[bins];
        int[] cntR2Bins = new int[bins];
        System.out.println("diversityPipeline processing: " + infileName);
        if (infileName.contains(".txt") == false) {
            System.out.println("File type unrecognized ");
            return;
        }
         Alignment ngai = ImportUtils.createPack1AlignmentFromFile(infileName, null, null);
  //      NextGenAlignmentWithCoveragev1_5 ngai = new NextGenAlignmentWithCoveragev1_5(numberOfTaxa, infileName, codeIndelsAsSNPs);
        Datum dtAll = new Datum("Alignment", ngai, "Blah");
        //       int[] segSites=AnnotatedAlignmentUtils.getIncludedSitesBasedOnFreqIgnoreGapsMissing(ngai,0.08,20);
//    System.out.println("segSites.length:"+segSites.length);
        FilterAlignmentPlugin fap = new FilterAlignmentPlugin(null, false);
        fap.setEnd(ngai.getSiteCount());
        fap.setMinFreq(-1);
        fap.setMinCount(1);
        int window = 200, step = 100;
        double[][] maxLDr2PerSite = new double[2][ngai.getSiteCount()];
        Arrays.fill(maxLDr2PerSite[0], -1);
        Arrays.fill(maxLDr2PerSite[1], -1);
//        ArrayList<SimpleTableReport> combinedrs=new ArrayList<SimpleTableReport>();
        for (int i = 0; i < ngai.getSiteCount() - window; i += step) {
            fap.setStart(i);
            fap.setEnd(i + window);
//     System.out.println("Current i:"+i);
            Datum dtGood = fap.processDatum(dtAll, false);
            Alignment aa = (Alignment) dtGood.getData();
            LinkageDisequilibrium ld = new LinkageDisequilibrium(aa, 20);
            ld.run();
            for (int j = 0; j < ld.getSiteCount(); j++) {
                for (int k = 0; k < j; k++) {
                    double dist = Math.abs(ld.getAnnotatedAlignment().getPositionInLocus(j) - ld.getAnnotatedAlignment().getPositionInLocus(k));
                    if (Double.isNaN(ld.getRSqr(j, k)) == false) {
                        int cbin = (int) Math.round(dist / 100);
                        if ((ngai.getSiteSummary(j).getAlleleCounts()[1] > 1) && (ngai.getSiteSummary(k).getAlleleCounts()[1] > 1) && (cbin < bins - 1)) {
                            sumR2Bins[cbin] += ld.getRSqr(j, k);
                            cntR2Bins[cbin]++;
                            if (ld.getRSqr(j, k) > maxLDr2PerSite[0][j + i]) {
                                maxLDr2PerSite[0][j + i] = ld.getRSqr(j, k);
                            }
                            if (ld.getRSqr(j, k) > maxLDr2PerSite[0][k + i]) {
                                maxLDr2PerSite[0][k + i] = ld.getRSqr(j, k);
                            }
                            if (dist > minDistanceApart) {
                                if (ld.getRSqr(j, k) > maxLDr2PerSite[1][j + i]) {
                                    maxLDr2PerSite[1][j + i] = ld.getRSqr(j, k);
                                }
                                if (ld.getRSqr(j, k) > maxLDr2PerSite[1][k + i]) {
                                    maxLDr2PerSite[1][k + i] = ld.getRSqr(j, k);
                                }
                            }
                        }

                    }
                }
            }
//            SimpleTableReport tr=(SimpleTableReport)rs.getData(1).getData();
//            combinedrs.add(tr);
//            System.out.println(TableReportUtils.toDelimitedString(tr, "\t"));
        }
//         SimpleTableReport combinedReport = (SimpleTableReport) SimpleTableReport.getInstance(combinedrs);
//         String output=TableReportUtils.toDelimitedString(combinedReport, "\t");
        try {
            BufferedWriter fileOut = new BufferedWriter(new FileWriter(outfileName, true), 100000);
            BufferedWriter fileAvgLDOut = new BufferedWriter(new FileWriter(outfileName + "avg.txt", true), 100000);
            DecimalFormat f2 = new DecimalFormat("#0.00");
            for (int i = 0; i < ngai.getSiteCount(); i++) {
                int[] counts = ngai.getSiteSummary(i).getAlleleCounts();
                if ((maxLDr2PerSite[0][i] < 0) || (counts[1] < 2)) {
                    continue;
                }
                fileOut.write(ngai.getLocusName(i) + "\t" + f2.format(ngai.getPositionInLocus(i)) + "\t" +
                        maxLDr2PerSite[0][i] + "\t" + maxLDr2PerSite[1][i] + "\t");
                fileOut.write(counts[0] + "\t" + counts[1] + "\t" + counts[2]);
                fileOut.write("\n");
            }
            for (int i = 0; i < bins; i++) {
                fileAvgLDOut.write(ngai.getLocusName(0) + "\t" + i + "\t" + sumR2Bins[i] + "\t" + cntR2Bins[i] + "\n");
            }
            fileOut.close();
            fileAvgLDOut.close();
        } catch (Exception e) {
            System.out.println("File output in diversityPipeline: " + e);
        }
    }

    private void illuminaSNPPipeline(String infileName, String refGenomeFileName, String outfileName) {
        System.out.println("diversityPipeline processing: " + infileName);
        if (infileName.contains(".txt") == false) {
            System.out.println("File type unrecognized ");
            return;
        }

        Alignment ngai = ImportUtils.createPack1AlignmentFromFile(infileName, null, null);
        int outCount = 0;
        try {
            BufferedWriter fileOut = new BufferedWriter(new FileWriter(outfileName, true), 100000);
            RandomAccessFile refGenome = new RandomAccessFile(refGenomeFileName, "r");
            DecimalFormat f2 = new DecimalFormat("#0.00");
            byte[] refRightSeqContext = new byte[50];
            byte[] refLeftSeqContext = new byte[50];
            for (int i = 1; i < ngai.getSiteCount() - 1; i++) {
                double pos = ngai.getPositionInLocus(i);
                double lpos = ngai.getPositionInLocus(i - 1);
                double rpos = ngai.getPositionInLocus(i + 1);
                int[] counts = ngai.getSiteSummary(i).getAlleleCounts();
                if ((Math.round(pos) == pos) && ((pos - lpos) > 20) && ((rpos - pos) > 20) && (counts[1] > 2) && (counts[0] > 2)) {
                    int base = (int) pos;
                    outCount++;
                    refGenome.seek(base - 1);
                    char refChar = (char) refGenome.readByte();
                    refGenome.seek(base - 51);
                    refGenome.read(refLeftSeqContext, 0, 50);
                    refGenome.seek(base);
                    refGenome.read(refRightSeqContext, 0, 50);
                    // String s=new String(refLeftSeqContext);
                    fileOut.write("PZE-" + (currentSNP + (Integer.parseInt(ngai.getLocusName(i)) * 1000000) + i) + ",");
                    //            fileOut.write(ngai.getRefData(i)+","+ngai.getAltData(i)+","+refChar+",");
                    fileOut.write(new String(refLeftSeqContext));
                    fileOut.write("[" + ngai.getReferenceAllele(i) + "/" + ngai.getAlleles(i)[1] + "]");
                    fileOut.write(new String(refRightSeqContext) + ",");
                    fileOut.write("1,");  //build
                    fileOut.write(ngai.getLocusName(i) + ",");  //chromosome
                    fileOut.write(Math.round(pos) + ",");  //chromosome position
                    fileOut.write(",Buckler,0,diploid,Zea mays,forward,0,,,,,");
                    fileOut.write("\n");
                }
            }
            fileOut.flush();
            fileOut.close();
        } catch (Exception e) {
            System.out.println("File output in diversityPipeline: " + e);
        }
        System.out.println("infileName: " + infileName + "\toutCout" + outCount);
    }

    private double getCML52RILidentityWB73(double base) {
        double expIdentity = 0.5;
        if (base < 7072703) {
            expIdentity = 0;
        } else if ((base > 7469703) && (base < 143767753)) {
            expIdentity = 1;
        } else if ((base > 143778669) && (base < 171107143)) {
            expIdentity = 0;
        } else if (base > 171107143) {
            expIdentity = 1;
        }
        return expIdentity;
    }

    public void run() {
        if (mode.equals("-S")) {
            NucleotideDiversityProcessing(this.infileForThread, this.outfileForThread, false, false, 8);
        } else if (mode.equals("-HP1")) {
            QSNPCallsToHapMapFormat(this.infileForThread, this.outfileForThread, false, true);
        } else if (mode.equals("-D")) {
            diversityPipeline(this.infileForThread, this.outfileForThread);
        } else if (mode.equals("-LD")) {
            ldPipeline(this.infileForThread, this.outfileForThread);
        } else if (mode.equals("-SNP")) {
            illuminaSNPPipeline(this.infileForThread, this.refGenomeFile, this.outfileForThread);
        } else if (mode.equals("-ERR")) {
            characterizeErrors(this.infileForThread, this.outfileForThread);
        } else if (mode.equals("-CON")) {
            evaluateConcordance(this.infileForThread, validDataInfile, this.outfileForThread);
        } else if (mode.equals("-OUT")) {
            createFileWithOutgroup(this.infileForThread, this.infile2ForThread, this.outfileForThread);
        }
    }

    public static void main(String[] args) throws Exception {


        int setThreads = 12;
        if (args.length == 0) {
            args = argAlt.split(" ");
        }
        File parentDirectory = new File(args[1]);
        File[] fileList = parentDirectory.listFiles();
        Thread[] tmft = new Thread[setThreads];
        for (int i = 0; (i < tmft.length) && (i < fileList.length); i++) {
            tmft[i] = new Thread(new EdDivProfilerV2(args, fileList[i].getName()));
            tmft[i].setName("tmft" + i);
            tmft[i].start();
            while (tmft[i].isAlive()) {
                Thread.sleep(1000);
            }
//             while(Thread.activeCount()>3) {
//                 Thread.sleep(1000);
//             }
        }

//         theDiversityProfiler.createPositiveTraining(parentInDirectory + fileList[index].getName(), validDataInfile,
//                                      parentOutDirectory+ stem+NucDivSuffix);
//  theDiversityProfiler.CreateTrainingDataset(parentInDirectory + fileList[index].getName(),
//                                      parentOutDirectory+ stem+
//                                      NucDivSuffix);

    }
}

