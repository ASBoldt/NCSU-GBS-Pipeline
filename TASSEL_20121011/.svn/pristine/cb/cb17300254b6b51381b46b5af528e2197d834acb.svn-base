package net.maizegenetics.genome.solexa;

/**
 * User: Ed and Jason
 * Date: Feb 17, 2009
 */
import net.maizegenetics.pal.alignment.AllelePositionBLOBUtils;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.popgen.LinkageDisequilibrium;
import net.maizegenetics.pal.statistics.ContigencyTable;
import net.maizegenetics.pal.statistics.FisherExact;

import java.io.*;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.regex.Pattern;
import java.util.regex.Matcher;
import net.maizegenetics.genome.HaplotypeLengthDistributionV2;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.FilterAlignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.Pack1Alignment;
import net.maizegenetics.pal.datatype.DataType;

public class EdHapMapV2 implements Runnable {

   // static String baseDir="/Users/edbuckler/SolexaAnal/HapMapV2/";
    static String baseDir="/Users/edbuckler/SolexaAnal/HapMapV2/jermpipe/SNPS201010/";
 //   private static String argAlt="-S "+baseDir+"calls/ "+baseDir+"test/"; //make  SNPs calls on
//    private static String argAlt="-D "+baseDir+"calls/ "+baseDir+"test/"; //get average sequencing depth by line for chromosome
//private static String argAlt="-HP1 "+baseDir+"log2/ "+baseDir+"test/"; //convert to Hapmap file
 //   private static String argAlt="-HP1 "+baseDir+"jermpipe/R101022/CSHLPhase2/SWA/ "+baseDir+"test/"; //convert to Hapmap file
 //     private static String argAlt="-LD "+baseDir+"hp1/ "+baseDir+"test/"; //estimate LD within chromosome old way
    private static String argAlt="-LD "+baseDir+"fusionmaize/ "+baseDir+"test/"; //estimate LD within chromosome old way
 //     private static String argAlt="-LDR "+baseDir+"hp1/ "+baseDir+"test/"; //estimate LD within the chromosome using bit approach with distribution produced
//      private static String argAlt="-SUB "+baseDir+"hp1/ "+baseDir+"ld/ "+baseDir+"test/ -s_T"; //subset by whether in good LD or not -s_T is in LD, -s_F is bad LD
 //     private static String argAlt="-HAP "+baseDir+"anchor/ "+baseDir+"test/ "; //subset by whether in good LD or not -s_T is in LD, -s_F is bad LD
//      private static String argAlt="-SUB "+baseDir+"anchorLD/ "+baseDir+"hap/ "+baseDir+"test/ -s_T -h"; //subset by whether in good LD or not -s_T is in LD, -s_F is bad LD
//   private static String argAlt="-LDC "+baseDir+"float/ "+baseDir+"anchor/ "+baseDir+"test/"; //find LD between two alignments
   //         private static String argAlt="-HAP "+baseDir+"anchorLDHap/ "+baseDir+"test/ "; //subset by whether in good LD or not -s_T is in LD, -s_F is bad LD

     //post anchor pipeline
//     private static String argAlt="-SwA "+baseDir+"calls/ "+baseDir+"anchorLDHap/ "+baseDir+"test/"; //make SNPs calls with an anchor map



    private boolean inGoodLD=true, isLDFile=true;  //a switch to determine what to save when splitting files
    private static int maxReadValue = 1000000;

            //quick harvest
//    private static int minTaxaWithSNP = 50, minHomoLines=1;
//    private static boolean useContigencyToEvalSNPs=false, extractSingleton=true;
//    private static double minLogPValueSNP = 4, minHomoProp=1, minAlleleQual=20;

        //quick harvest
    private static int minTaxaWithSNP = 12, minHomoLines=1;
    private static boolean useContigencyToEvalSNPs=true, extractSingleton=false;
    private static double minLogPValueSNP = 1.0, minHomoProp=0.9, minAlleleQual=20;
    private static double maximumHeterozygousRatio=2.001;  //this is the ratio of heterozygous lines to homozygous minor allele lines

    //final harvest
//    private static int minTaxaWithSNP = 12, minHomoLines=1;
//    private static boolean useContigencyToEvalSNPs=true;
//    private static double minLogPValueSNP = 2, minHomoProp=0.9, minAlleleQual=20;

    //anchor map
//    private static int minTaxaWithSNP = 50, minHomoLines=4;
//    private static boolean useContigencyToEvalSNPs=true;
//    private static double minLogPValueSNP = 2.5, minHomoProp=0.9, minAlleleQual=20;

    //HP1 processing
    int minFreqInHP1=12;
    boolean includeHetsInHP1=true;
    boolean filterByLogRScore=true;
    double minLogRScore=0.00;
    double hetsThresholdInHP1=0.9;
    //LD permutation
    int ldWindow=200, ldPermutations=5000, minDistanceForLD=500;
    double ldSubsetPValue=0.00001;

    private ContigencyTable contingencyTable;
    private FisherExact fishersExact;
    private String mode = "",  infileForThread = "",  infile2ForThread = "",  outfileForThread = "", outfile2ForThread = "",
            depthfileForThread="";
    //private int currentSNP = 100000000;


    public EdHapMapV2(String[] args, String infile) {
        this.contingencyTable = new ContigencyTable(maxReadValue);
        this.fishersExact = new FisherExact(maxReadValue);
        this.mode = args[0];
        this.infileForThread = args[1] + infile;
        this.infile2ForThread = "";
        this.outfileForThread = "";
        String stem = infile;
        String[] namesp = stem.split("\\.");
        if (mode.equals("-LD")) {
            outfileForThread = args[2] + stem.replaceFirst(".hmp.txt", ".LD.txt");   //for SNP calling file
        } else if (mode.equals("-S")) {
            int homoInt=(int)(minHomoProp*100);
            int x2Int=(int)minLogPValueSNP;
            String cttext=(useContigencyToEvalSNPs)?"ct":"cx";
            outfileForThread = args[2] + stem.replaceFirst(".txt", ".log"+x2Int+cttext+"_h"+homoInt+".txt");   //for SNP calling file
        } else if (mode.equals("-SwA")) {
            int homoInt=(int)(minHomoProp*100);
            int x2Int=(int)minLogPValueSNP;
            String cttext=(useContigencyToEvalSNPs)?"ct":"cx";

			String chrsearch="chr\\d+";
			Pattern pattern = Pattern.compile(chrsearch);
			Matcher matcher = pattern.matcher(infile);

			String chrname="";
			if (matcher.find()){
				chrname=matcher.group();
			}else{
				System.out.println("Unable to locate \"" + chrsearch + "\" in \"" + infile + "\"");
				return;

                                //System.exit(0);
			}
			infile2ForThread = args[2] + chrname + ".genotypes.log2_h90_f50.good.good.hmp.txt";
            // infile2ForThread = args[2] + stem.replaceFirst(".txt", ".log2_h90_f50.good.good.hmp.txt");  //if would be better to just look in the folder and find the match

            outfileForThread = args[3] + stem.replaceFirst(".txt", ".log"+x2Int+cttext+"_h"+homoInt+".txt");   //for SNP calling file
            outfile2ForThread = args[3] + stem.replaceFirst(".txt", ".callrep.txt");   //for SNP calling file
        }else if (mode.equals("-D")) {
            outfileForThread = args[2] + stem.replaceFirst(".txt", ".depth.txt");   //for determing read depth by line file
        } else if (mode.equals("-HP1")) {
            outfileForThread = args[2] + stem.replaceFirst(".txt", "_f"+minFreqInHP1+".hmp.txt");   //for converting to hapmap file
        } else if (mode.equals("-LDR")) {
            outfileForThread = args[2] + stem.replaceFirst(".hmp.txt", ".LD.txt");   //for characterizing LD
        } else if (mode.equals("-HAP")) {
            outfileForThread = args[2] + stem.replaceFirst(".hmp.txt", ".ehh.txt");   //for characterizing haplotype structure
        } else if (mode.equals("-LDC")) {
            outfileForThread = args[3] + stem+".ULD.txt";   //for SNP calling file
            infile2ForThread=args[2];  //directory of anchor alignments
        } else if (mode.equals("-SUB")) {
            if(args.length==5&&args[4].equals("-s_F")) inGoodLD=false;
            if(args.length==6&&args[5].equals("-h")) isLDFile=false;
            if(isLDFile) {infile2ForThread = args[2] + stem.replaceFirst(".hmp.txt",".LD.txt");}
            else {infile2ForThread = args[2] + stem.replaceFirst(".hmp.txt",".ehh.txt");}
            if(inGoodLD) {outfileForThread = args[3] + stem.replaceFirst(".hmp.txt",".good.hmp.txt");}
            else {outfileForThread = args[3] + stem.replaceFirst(".hmp.txt",".bad.hmp.txt");}
        }
    }

    /**=Read in a file containing snp calls across 27 founders
     * output bac accession, position in bac, polymorphism type, allele definition,
     * snp call, reference and alt reads in each taxon, total taxa with reads,
     * tropical with reads, temperate with reads, snp log from contingency table,
     * threshold for call, pi, variance of pi, fst, logP fst*/
    private void idSNPGenotypeFileToQSNPFile(String infileName, String outfileName, boolean isFreq, boolean isUsingContigency) {
        try {
            if(!infileName.endsWith(".txt")) return;

//            System.out.println(infileName+" rowCount:"+SimpleTextFile.countRows(infileName));
            BufferedReader fileIn = new BufferedReader(new FileReader(infileName), 100000);
            BufferedWriter fileOut = new BufferedWriter(new FileWriter(outfileName), 100000);
//            BufferedWriter fileTrainOut = null;
            String taxaHeaderString=fileIn.readLine();
            taxaHeaderString=taxaHeaderString.split(":\t")[1];
            String[] taxaNames=taxaHeaderString.split("\t");
            int taxaCnt=taxaNames.length;
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
            int errorCnt=0;
            while ((row=fileIn.readLine())!=null) {
                try {
                    //row = fileIn.readLine();
                    if ((!row.contains("Sample"))&&(!row.contains("line"))) {
                        SNPDistV2[] theSNPDist=new SNPDistV2[2];
                        theSNPDist[0] = new SNPDistV2(taxaCnt, row, contingencyTable, fishersExact,1);
                        if(theSNPDist[0].afreq[2]>0) theSNPDist[1] = new SNPDistV2(taxaCnt, row, contingencyTable, fishersExact,2);
                        for(SNPDistV2 tsd:theSNPDist) {
                            if(tsd==null) continue;
                            if (tsd.taxaNumWithReads > minTaxaWithSNP) {
                                countGreaterThanMinTaxaWithSNP++;
                                if(isUsingContigency) {tsd.scoreSNPX2ThenContigency();}  //this the the more rigorous approach
                                else  {tsd.scoreSNPFastX2(); tsd.snpLogP=tsd.snpX2P;}  //fast but not as good
                                //theSNPDist.scoreMaxMachineLearnScore();
                                int[] homoCnts=tsd.getHomozygousCounts();
                                double propHomozygous=((double)homoCnts[0]+(double)homoCnts[1])/(double)homoCnts[2];
                               // System.out.println(tsd.altAlleleNumber+":"+tsd.toStringProp());
                                if ((tsd.snpLogP >= minLogPValueSNP) &&
                                        (propHomozygous>=minHomoProp)&&
                                        (homoCnts[0]>=2)&&(homoCnts[1]>=2)&&
                                        (tsd.avgQual[0]>minAlleleQual)&&
                                        (tsd.avgQual[tsd.altAlleleNumber]>minAlleleQual)
                                        )
                                    {
                                   //     System.out.println(theSNPDist.snpLogP+" "+theSNPDist.snpX2P);
                                    countGreaterThanMinLogPValue++;
    //                                theSNPDist.calcMaxThresholdPiFst();
                                    if (isFreq) {
                                        fileOut.write(tsd.toStringNuc() + "\n");
                                        System.out.println("SNP number:");
                                    } else {
                                      //  System.out.println(tsd.altAlleleNumber+":"+tsd.toStringProp());
                                        fileOut.write(tsd.toStringProp() + "\n");
                                    }
                                //                      fileOut.write(chrinfo[0]+"\t"+chrinfo[1]+"\t"+ibd+"\n");
                                }
                            }
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

    /**
     * This creates the Human Genetics community HapMap format from our quantiative SNP call files
     *
     * @param infileName
     * @param outfileName
     */
    private void QSNPFileToHapMapFormat(String infileName, String outfileName, int minPresent, boolean filterByThreshold,
            double minLogRScore, boolean includeIndels, boolean includeHets) {
        int totalSNPIndelSites = 0;
        DecimalFormat chrForm = new DecimalFormat("00");
        DecimalFormat siteForm = new DecimalFormat("00000000");
		DecimalFormat qcForm = new DecimalFormat("00");
        int numSeqs = 105;
        int nanCnt=0, refCnt=0, altCnt=0, hetCnt=0;
        try {
            BufferedReader fileIn = new BufferedReader(new FileReader(infileName), 1000000);
            BufferedWriter fileOut = new BufferedWriter(new FileWriter(outfileName, false), 100000);
            String[] header=fileIn.readLine().split("[\\s]+");
            numSeqs=header.length-21;
            String[] taxaOrderName=new String[numSeqs];
            for (int i = 0; i < numSeqs; i++) {
                taxaOrderName[i] = header[i+10];
                if(taxaOrderName[i].startsWith("TI")) {
                    taxaOrderName[i]=taxaOrderName[i]+":TEO";
                }  else if ( taxaOrderName[i].startsWith("TDD")) {
					taxaOrderName[i]=taxaOrderName[i]+":TD";
				}	else{
                    taxaOrderName[i]=taxaOrderName[i]+":MZ";
                }
            }
            fileOut.write("rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode\t");
            for (String t : taxaOrderName) {
                String tout=t.toUpperCase().replaceAll("_", "");
                fileOut.write(tout + "\t");
            }
            fileOut.write("\n");
            System.out.println("Setting states based on " + infileName);
            String sl="";

            while ((sl=fileIn.readLine())!=null) {
                int nanCntSite=0, refCntSite=0, altCntSite=0, hetCntSite=0;
                SNPDistV2 theSNPDist = new SNPDistV2(numSeqs, sl, null, null);
                if ((filterByThreshold) && (theSNPDist.maxMLScore <minLogRScore)) {
                    continue;
                }
                if ((includeIndels == false) && (!theSNPDist.varType.startsWith("SNP"))) {
                    continue;
                }
                char altBase = 'S', refBase = 'S', hetBase='S';
                refBase = theSNPDist.alleleDef[0].charAt(0);
                altBase = theSNPDist.alleleDef[1].charAt(0);
                byte[] stateByte={(byte)refBase, (byte)altBase};
                hetBase=(char)AllelePositionBLOBUtils.getBaseFromHalfByte(AllelePositionBLOBUtils.getHalfByteFromSNPValue(stateByte));
                if(theSNPDist.varType.equals("IDP")) {
                    altBase=theSNPDist.alleleDef[1].charAt(0);
                    refBase=(altBase=='-')?'+':'-';
                    hetBase='0';
                }
                if(includeHets==false) {hetBase=DataType.UNKNOWN_CHARACTER;}
                StringBuilder sb = new StringBuilder();
                sb.append("PZE");
                sb.append(chrForm.format(theSNPDist.chromosome));
                sb.append(siteForm.format(theSNPDist.startPos));
                sb.append("\t");
                //String[] states = theSNPDist.alleleDef.split("/");
                //sb.append(theSNPDist.alleleDef[0] + "/" + theSNPDist.alleleDef[1] + "\t");
                sb.append(refBase + "/" + altBase + "\t");
                sb.append(theSNPDist.chromosome + "\t");
                sb.append((int) theSNPDist.startPos + "\t");
                sb.append("+\tAGPv1\tMaizeDiv\tSBS\tMHPv1\tNAMfnd\t");
                if(filterByThreshold) {sb.append(qcForm.format(Math.round(theSNPDist.maxMLScore*100))+"\t");}
                        else {sb.append("QC+\t");}

  //              System.out.println(refBase+" "+altBase+" "+hetBase);
                boolean isVariant = false;  //there are rare situations where a site may have the correct threshold, but no single line passes
                //the test below.  We track this below.
                for (double map : theSNPDist.minorAlleleProp) {
                    if (Double.isNaN(map)) {
                        sb.append(DataType.UNKNOWN_CHARACTER + "\t");
                        nanCnt++; nanCntSite++;
                    } else if (map > hetsThresholdInHP1 ) {
                        sb.append(altBase + "\t");
                        isVariant = true;
                        altCnt++; altCntSite++;
                    } else if (map < (1-hetsThresholdInHP1 )) {
                        sb.append(refBase + "\t");
                        refCnt++; refCntSite++;
                    } else {
                        sb.append(hetBase + "\t");
                        hetCnt++; hetCntSite++;
                    }
                }
//                for(double map: theSNPDist.minorAlleleProp) {
//                    if(Double.isNaN(map)) {sb.append(NEW_UNKNOWN_CHARACTER+""+NEW_UNKNOWN_CHARACTER+"\t");}
//                    else if(map>(theSNPDist.maxThreshold*1.5)) {sb.append(altBase+""+altBase+"\t"); isVariant=true;}
//                    else if(map>(theSNPDist.maxThreshold*0.5)) {sb.append(refBase+""+altBase+"\t");}
//                    else {sb.append(refBase+""+refBase+"\t");}
//                }
                int cntSite=altCntSite+refCntSite;
                if(includeHets) cntSite+=hetCntSite;
                double propHomo=1-((double)hetCntSite/(double)cntSite);
     if(propHomo<0.9) continue;
                if ((isVariant)&&(cntSite>=minPresent)) {
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



    private void ldPipeline(String infileName, String outfileName) {
        int minDistanceApart = 200;
        boolean permute=false;
        int bins = 5000;
        int window = 300000;
        int minInbredComparisons=40;
        int minMinorInbred=2;
        int maxSite=100000;
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
        ngai = HaplotypeLengthDistributionV2.makeHomozygousAlignment((Pack1Alignment)ngai);
        
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
//        maxSite=ngai.getSiteCount();
        maxSite=ngai.getSiteOfPhysicalPosition(ngai.getPositionInLocus(ngai.getSiteCount()-1)-window,null);
        if(maxSite<0) maxSite=-(maxSite+1);
        for (int i = 0; i < maxSite; i++) {
            if(worthySites[i]==false) {conLostWorthy+=window; continue;}
            int eSite=ngai.getSiteOfPhysicalPosition(ngai.getPositionInLocus(i)+window,null);
            if(eSite<0) eSite=-(eSite+1);
            for (int ts = i+1; ts < eSite; ts++) {
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

        int[] r2NearBin=new int[101];
        int[] r2FarBin=new int[101];
        for (int i = 0; i < maxSite; i++) {
            if(bestLDpSite[1][i]<0) {r2NearBin[0]++;}
            else {r2NearBin[(int)Math.round(bestLDpSite[1][i]*100)]++;}
            if(bestLDpSite[4][i]<0) {r2FarBin[0]++;}
            else {r2FarBin[(int)Math.round(bestLDpSite[4][i]*100)]++;}
        }
        
        try {
            BufferedWriter fileOut = new BufferedWriter(new FileWriter(outfileName, false), 100000);
            BufferedWriter fileAvgLDOut = new BufferedWriter(new FileWriter(outfileName + "avg.txt", false), 100000);
            DecimalFormat f2 = new DecimalFormat("#0.00");
            fileOut.write("r2bin\tr2all\tr2far\n");
            for (int i = 0; i < r2NearBin.length; i++) {
                fileOut.write((i-1)+"\t"+r2NearBin[i]+"\t"+r2FarBin[i]+"\n");
            }
//            fileOut.write("Chr\tPosition\tSite\tN\tr2\tP\tfar_N\tfar_r2\tfar_P\tCommonCnt\tLessCommCnt\n");
//            for (int i = 0; i < maxSite; i++) {
//                int majorCnt=countAlleles(ngai,i,ngai.getMajorAllele(i));
//                int minorCnt=countAlleles(ngai,i,(byte)ngai.getMinorAllele(i));
//                fileOut.write(ngai.getLocusName(i) + "\t" + f2.format(ngai.getPositionInLocus(i)) + "\t");
//                fileOut.write(i+"\t");
//                fileOut.write(bestLDpSite[0][i] + "\t" + bestLDpSite[1][i] + "\t"+ bestLDpSite[2][i] + "\t"+
//                        bestLDpSite[3][i] + "\t" + bestLDpSite[4][i] + "\t" + bestLDpSite[5][i] + "\t");
//                fileOut.write(majorCnt+"\t" + minorCnt +"\t");
//                fileOut.write("\n");
//            }
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
     * Calculate LD within a chromosome rapidly, creates the distrubution site by site
     * @param infileName Alignment file (generally HapMap)
     * @param outfileName textfile with the LD results
     */
    private void ldRapidPipeline(String infileName, String outfileName) {
        System.out.println("LDPipeline processing: " + infileName);
        if (infileName.contains(".txt") == false) {
            System.out.println("File type unrecognized ");
            return;
        }
        System.out.println("Loading Alignment");
        Alignment ngai = ImportUtils.createPack1AlignmentFromFile(infileName, null, null);
         LDbyBitsV3.AlleleMode[] theModes={LDbyBitsV3.AlleleMode.RefvAlt};
//        LDbyBitsV2.AlleleMode[] theModes={LDbyBitsV2.AlleleMode.RefvAlt, LDbyBitsV2.AlleleMode.MissVPres};
        LDbyBitsV3 lbb=new LDbyBitsV3(ngai,theModes.length);
        ngai=null;
        lbb.runNeighbors(ldWindow, minDistanceForLD, ldPermutations, theModes);
 //       lbb.createDistForLocal(2);
        lbb.writeResultsToFile(outfileName);
    }

       /**
     * Calculate LD within a chromosome rapidly, creates the distrubution site by site
     * @param infileName Alignment file (generally HapMap)
     * @param outfileName textfile with the LD results
     */
    private void ldComparisonPipeline(String infileName, String targetFile, String outfileName) {
        System.out.println("LDPipeline processing: " + infileName);
        if (infileName.contains(".txt") == false) {
            System.out.println("File type unrecognized ");
            return;
        }
        System.out.println("Loading Alignment");

        Alignment ngai = ImportUtils.createPack1AlignmentFromFile(infileName, null, null);
        LDbyBitsV2 lbb=new LDbyBitsV2(ngai,2);

        File parentDirectory = new File(targetFile);
        File[] fileList = parentDirectory.listFiles();
        for (int i = 0; i < fileList.length; i++) {
            System.out.println("Loading Target Alignment Alignment:"+fileList[i].getName());
            if (fileList[i].getName().contains(".txt") == false) {System.out.println("Error cannot be loaded"); continue;}
            Alignment targetAlignment = ImportUtils.createPack1AlignmentFromFile(fileList[i].getPath(), null, null);
            System.out.println("Starting comparison");
            lbb.compareWithOther(targetAlignment);
            lbb.writeResultsToFile(outfileName);
        }
    }

     private void snpCallerWithAnchorPipeline(String infileSiteCallsName, String anchorFileName, String outfileGoodSitesName,
            String outfileSiteAttributesName) {
        System.out.println("snpCallerWithAnchorPipeline processing: " + infileSiteCallsName);
        if (infileSiteCallsName.contains(".txt") == false) {
            System.out.println("File type unrecognized ");
            return;
        }
        System.out.println("Loading Alignment");

        SNPCallerAgainstAnchorV2 scaa=new SNPCallerAgainstAnchorV2(infileSiteCallsName, anchorFileName, outfileGoodSitesName,
            outfileSiteAttributesName);
        scaa.setMinLogPValueSNP(minLogPValueSNP);
        scaa.setMinAlleleQual(minAlleleQual);
        scaa.setMinHomoProp(minHomoProp);
        scaa.setMinTaxaWithSNP(minTaxaWithSNP);
        scaa.setMaximumHeterozygousRatio(maximumHeterozygousRatio);
        scaa.run();
    }


     /**
     * Calculate LD within a chromosome rapidly, creates the distrubution site by site
     * @param infileName Alignment file (generally HapMap)
     * @param outfileName textfile with the LD results
     */
    private void haplotypePipeline(String infileName, String outfileName) {
        System.out.println("HaplotypePipeline processing: " + infileName);
        if (infileName.contains("hmp.txt") == false) {
            System.out.println("File type unrecognized ");
            return;
        }
        System.out.println("Loading Alignment");
        Alignment ngai = ImportUtils.createPack1AlignmentFromFile(infileName, null, null);
        System.out.println("Alignment Loaded:"+ngai.getLocusName(0));
        HaplotypeLengthDistributionV2 hld2=new HaplotypeLengthDistributionV2(ngai, outfileName);
    }

     private void subsetHP1(String infileName, String subsetFile, String outfileName, boolean outInGoodLD, boolean isLD) {
        System.out.println("Subset SNP in good and bad LD Pipeline\n processing: " + infileName);
        if (infileName.contains(".txt") == false) {
            System.out.println("File type unrecognized ");
            return;
        }
        System.out.println("Loading Alignment");
        Alignment ngai = ImportUtils.createPack1AlignmentFromFile(infileName, null, null);
        System.out.println("Alignment Loaded:"+ngai.getLocusName(0));
        String[] labels;
        Object[] o;
        if(isLD) {
            labels="Site Chr Position NbrRefvAlt_MajorCnt NbrRefvAlt_MinorCnt NbrRefvAlt_BChr NbrRefvAlt_BestPos NbrRefvAlt_BestP NbrRefvAlt_NumTests NbrRefvAlt_LEThres NbrRefvAlt_RandLTE_P".split(" ");
            Object[] t={Integer.class,Integer.class,Integer.class,Integer.class,Integer.class,Integer.class,Integer.class,
                Double.class, Integer.class,Integer.class,Double.class};
            o=t;
        } else {
            labels="Position Site MAF ObsHapLength ObsWOCurrSite RandomGreaterObs AvgRandom StDevRandom ZobsvRand ImpCorr ImpTotal ImpLongCorr ImpLongTotal ImpMinorCorr ImpMinorTotal".split(" ");
            Object[] t={Integer.class,Integer.class,Double.class,Double.class,Double.class,Integer.class,Double.class,
                Double.class, Double.class, Integer.class, Integer.class, Integer.class, Integer.class, Integer.class, Integer.class};
            o=t;
        }
        SimpleTextFile theSTF=new SimpleTextFile(subsetFile,labels,o,true);
        int[] tempSites=new int[theSTF.getRows()];
        int index=0;
        for (int i = 0; i <theSTF.getRows(); i++) {
            boolean keep;
            if(isLD) {keep=((theSTF.getDoubleElement(10, i)<ldSubsetPValue)&&(theSTF.getDoubleElement(7, i)<ldSubsetPValue));}
            else {keep=((double)(theSTF.getIntElement(14, i)-theSTF.getIntElement(13, i))/(double)theSTF.getIntElement(14, i)<0.2);}
             if(outInGoodLD==keep) {
                 if(isLD) {tempSites[index]=theSTF.getIntElement(0, i);}
                 else {tempSites[index]=theSTF.getIntElement(1, i);}
                 index++;
             }
         }
        int[] goodSites=new int[index];
        System.arraycopy(tempSites, 0, goodSites, 0, goodSites.length);
        FilterAlignment fa=FilterAlignment.getInstance(ngai,goodSites);
         System.out.println("Alignment filtered sites:"+fa.getSiteCount());
        ExportUtils.writeToHapmap(fa, false, outfileName, ' ');


    }

    /**
     * Thoughts permute to determine null
     * split alleles - alt1 and alt2
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



    public void run() {
        if(!infileForThread.endsWith(".txt")) return;
        if (mode.equals("-S")) {
//            idSNPGenotypeFileToQSNPFile(this.infileForThread, this.outfileForThread, false, useContigencyToEvalSNPs);
            SNPCaller.idSNPGenotypeFileToQSNPFile(this.infileForThread, this.depthfileForThread, this.outfileForThread, minTaxaWithSNP,
                    minLogPValueSNP, minHomoProp, minAlleleQual, minHomoLines, useContigencyToEvalSNPs, extractSingleton);
        } else if (mode.equals("-SwA")) {

        	System.out.println("EdHapMapV2 "+ this.mode + " " + infileForThread + " " + infile2ForThread);
            snpCallerWithAnchorPipeline(this.infileForThread, infile2ForThread,this.outfileForThread, this.outfile2ForThread);
        } else if (mode.equals("-D")) {
            SNPCaller.idSNPGenotypeFileAvgCoverage(this.infileForThread, this.outfileForThread);
        } else if (mode.equals("-HP1")) {
            QSNPFileToHapMapFormat(this.infileForThread, this.outfileForThread, minFreqInHP1, filterByLogRScore, this.minLogRScore, true, includeHetsInHP1);
        } else if (mode.equals("-LD")) {
            ldPipeline(this.infileForThread, this.outfileForThread);
        } else if (mode.equals("-LDR")) {
            ldRapidPipeline(this.infileForThread, this.outfileForThread);
        } else if (mode.equals("-LDC")) {
            ldComparisonPipeline(this.infileForThread, infile2ForThread,this.outfileForThread);
        } else if (mode.equals("-HAP")) {
            haplotypePipeline(this.infileForThread, this.outfileForThread);
        }
        else if (mode.equals("-SUB")) {
            subsetHP1(this.infileForThread, this.infile2ForThread, this.outfileForThread, inGoodLD, isLDFile);
        }
    }

    public static void main(String[] args) throws Exception {

        int numProcessors = Runtime.getRuntime().availableProcessors();
        int maxActiveThreads = numProcessors/2;
        if(numProcessors>4) {maxActiveThreads = numProcessors-2;}
//        if(numProcessors>4) {maxActiveThreads = 1;}
        System.out.printf("NumberProcessors %d AttemptThreads %d %n", numProcessors, maxActiveThreads);
        if (args.length == 0) {
            args = argAlt.split(" ");
        }

        System.out.println("EdHapMapV2 "+Arrays.toString(args));
		File[] fileList;

		if ( (args[0].equals("-SwA")) && (args.length == 5) ){
			fileList = new File[] { new File(args[4]) };
		}else if ( (args[0].equals("-HP1")) && (args.length == 4)){
		 	fileList = new File[] { new File(args[3]) };
		}else {
        	File parentDirectory = new File(args[1]);
        	fileList = parentDirectory.listFiles();
		}



        Thread[] tmft = new Thread[fileList.length+1];
        for (int i = 0; (i < tmft.length) && (i < fileList.length); i++) {
            tmft[i] = new Thread(new EdHapMapV2(args, fileList[i].getName()));
            tmft[i].setName("tmft" + i);
            tmft[i].start();
//           while (tmft[i].isAlive()) {
//                Thread.sleep(1000);
//            }
             while(Thread.activeCount()>maxActiveThreads) {
                 Thread.sleep(1000);
             }


        }

    }
}