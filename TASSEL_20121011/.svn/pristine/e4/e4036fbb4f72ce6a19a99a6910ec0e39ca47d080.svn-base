/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.gbs.pipeline;

import cern.jet.random.Binomial;
import edu.cornell.lassp.houle.RngPack.RandomJava;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import javax.swing.ImageIcon;
import net.maizegenetics.gbs.util.ArgsEngine;
import net.maizegenetics.gbs.util.MutableSimpleAlignment;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.FilterAlignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.datatype.DataType;
import net.maizegenetics.pal.datatype.IUPACNucleotides;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.pal.popgen.LinkageDisequilibrium;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;

/**
 * Tools for characterizing and correcting SNPs segregating in bi-parental populations.
 *
 * Error rates are bounded away from zero, but adding 0.5 error to all error rates that
 * that were observed to be zero.
 *
 * @author edbuckler
 */
public class BiParentalErrorCorrection extends AbstractPlugin {
    private double[] errorRate;
    private short[][][] popFreq;
    private double[][] popP;
    private double[][] ldByPop;
    private short[] popOfTaxa;
    private double expSegregation=0.5;
    private int minCountForDetectingError=19;
    private double maxErrorRate=0.05;
    private double minMedianLDR2=-1.0;
    private boolean myRemoveUntestedError=true;
    private boolean myRemoveUntestedLD=false;
    private String outHapMap=null;
    private String outFreqBin=null;
    private String outErrorRate=null;
    private ArrayList<String> popNames;
    private Alignment a;
    Binomial binomFunc=new Binomial(5, 0.5, new RandomJava());
    int homoweight=1;  //weight of homozygous genotypes to heterozgyous genotypes
    //1 assumes low coverage and each homozygote is only counted once, while hets get a count for each allele
    //2 assumes better coverage, and homozygotes are counted twice.
    private double minDistortionRatio=2.5;  //this is minimum distortion to be counted as an error
      //e.g. 2.0, means that the allele frequency has to be less the 25% if the expected is 50%
      


    private void removeErrorsInAlignment(double relRatio) {
        MutableSimpleAlignment msa=new MutableSimpleAlignment(a);
        for (int s = 0; s < a.getSiteCount(); s++) {
            byte mjb=a.getMajorAllele(s);
            byte mnb=a.getMinorAllele(s);
            byte hetb=IUPACNucleotides.getDegerateSNPByteFromTwoSNPs(mjb, mnb);
            byte[] pmnb=new byte[popNames.size()];
            byte[] pmjb=new byte[popNames.size()];
            double[] errVSeg=new double[popNames.size()];
            for (int pop = 0; pop < popNames.size(); pop++) {
                int pmj, pmn;  //pop major and minor

                if(popFreq[0][pop][s]<popFreq[1][pop][s]) {  //major and minor reversed in this population
                    pmj=popFreq[1][pop][s]; pmn=popFreq[0][pop][s]; pmnb[pop]=mjb; pmjb[pop]=mnb;
                } else {
                    pmj=popFreq[0][pop][s]; pmn=popFreq[1][pop][s]; pmnb[pop]=mnb; pmjb[pop]=mjb;
                }
                int present=pmj+pmn;
                if(present<1) {
                    errVSeg[pop]=0;}
                else {
                    binomFunc.setNandP(present,errorRate[s]);
                    double errorP = binomFunc.pdf(pmn);
                    binomFunc.setNandP(present,expSegregation);
                    double segP = binomFunc.pdf(pmn);
                    errVSeg[pop]=errorP/segP;
//                    System.out.println(errorRate[s]+"\t"+pmj+"\t"+pmn+"\t"+errorP+"\t"+segP+"\t"+errVSeg[pop]);
                }
            }
            for (int t = 0; t < a.getSequenceCount(); t++) {
                if(popOfTaxa[t]<0) continue;
                if(errVSeg[popOfTaxa[t]]<1) continue;
                byte b=a.getBase(t, s);
                if(b==DataType.UNKNOWN_BYTE) continue;
                if(b==pmnb[popOfTaxa[t]]) msa.setBase(t, s, DataType.UNKNOWN_BYTE);
                if(b==hetb) msa.setBase(t, s, pmjb[popOfTaxa[t]]);
            }
         }
        a=msa;
    }

    private void calcLDByPop(Alignment a) {
        int minCntForLD=10;
        ldByPop=new double[popNames.size()][a.getSiteCount()];
        IdGroup idg=a.getIdGroup();
        for (int pop = 0; pop < popNames.size(); pop++) {
            Arrays.fill(ldByPop[pop], Double.NaN);
            ArrayList<Identifier> keepNames=new ArrayList<Identifier>();
            for (int i = 0; i < popOfTaxa.length; i++) if(popOfTaxa[i]==pop) keepNames.add(idg.getIdentifier(i));
            Identifier[] ids=new Identifier[keepNames.size()];
            keepNames.toArray(ids);
            Alignment pa = FilterAlignment.getInstance(a, new SimpleIdGroup(ids));
            int[] segSites=AlignmentFilterByGBSUtils.getLowHetSNPs(pa, false, -2.0, minCntForLD, 0.15, 2);
            Alignment paf = FilterAlignment.getInstance(pa, segSites);
            int windowSize=paf.getSiteCount()/20;
            paf=new MutableSimpleAlignment(paf);
            LinkageDisequilibrium theLD=new LinkageDisequilibrium(paf,
                    minCntForLD, windowSize, LinkageDisequilibrium.testDesign.SlidingWindow);
            theLD.run();
            for (int i = 0; i < paf.getSiteCount(); i++) {
                ArrayList<Double> obsR2=new ArrayList<Double>();
                for(int j=i-(windowSize/2); j<i+(windowSize/2); j++) {
                    if(Double.isNaN(theLD.getRSqr(i, j))) continue;
                    if(Math.abs(paf.getPositionInLocus(i)-paf.getPositionInLocus(j))<100000) continue;
                    obsR2.add(theLD.getRSqr(i, j));
                }
                Collections.sort(obsR2);
                if(obsR2.size()>4) {
                    ldByPop[pop][a.getSiteOfPhysicalPosition(paf.getPositionInLocus(i), null)]=obsR2.get(obsR2.size()/2);
                }
            }
  //          System.out.println("POP:"+pop+Arrays.toString(ldByPop[pop]));
        }
    }

    private void removeHighErrorSites(double maxErrorRate, boolean removeUntested) {
        MutableSimpleAlignment msa=new MutableSimpleAlignment(a);
        int sitesWithHighError=0, untestedSites=0;
        for (int s = 0; s < a.getSiteCount(); s++) {
            if(Double.isNaN(errorRate[s])) {
                untestedSites++;
                if(removeUntested)  msa.clearSite(s);
            }
            else if(errorRate[s] > maxErrorRate) {
                msa.clearSite(s);
                sitesWithHighError++;
            }
         }
        System.out.printf("Initial Sites %d  ErrorRate %g  HighErrorSites %d UntestedSites %d %n",a.getSiteCount(),
                maxErrorRate, sitesWithHighError, untestedSites);
        msa.sortSiteByPhysicalPosition();
        a=msa;
        System.out.printf("Final Sites %d  ErrorRate %g  HighErrorSites %d %n",a.getSiteCount(), maxErrorRate, sitesWithHighError);
    }

    private void removeLowLDSites(boolean removeUntested) {
        MutableSimpleAlignment msa=new MutableSimpleAlignment(a);
        int sitesWithLowLD=0, untestedSites=0;
        for (int s = 0; s < a.getSiteCount(); s++) {
            ArrayList<Double> obsR2=new ArrayList<Double>();
            for (int pop = 0; pop < popNames.size(); pop++) {
                if(!Double.isNaN(ldByPop[pop][s])) obsR2.add(ldByPop[pop][s]);
            }
            Collections.sort(obsR2);
            double medianR2=Double.NaN;
            if(obsR2.size()>0) medianR2=obsR2.get(obsR2.size()/2);
            if(Double.isNaN(medianR2)) {
                untestedSites++;
                if(removeUntested)  msa.clearSite(s);
            }
            else if(medianR2 < minMedianLDR2) {
                msa.clearSite(s);
                sitesWithLowLD++;
            }
         }
        System.out.printf("Initial Sites %d  minMedianLDR2 %g  LowLDSites %d UntestedSites %d %n",a.getSiteCount(),
                minMedianLDR2, sitesWithLowLD, untestedSites);
        msa.sortSiteByPhysicalPosition();
        a=msa;
        System.out.printf("Final Sites %d  ErrorRate %g  LowLDSites %d %n",a.getSiteCount(), maxErrorRate, sitesWithLowLD);
    }

    private void classifyTaxaToPops(String popMask) {
        popNames=new ArrayList<String>();
        popOfTaxa=new short[a.getSequenceCount()];
        Arrays.fill(popOfTaxa, (short)-1);
        Pattern pat=Pattern.compile(popMask);
        for (int i=0; i<a.getSequenceCount(); i++) {
            Matcher m=pat.matcher(a.getTaxaName(i));
            if(m.find()) {
                String tn=a.getTaxaName(i).substring(m.start(), m.end());
                if(!popNames.contains(tn)) popNames.add(tn);
                popOfTaxa[i]=(short)popNames.indexOf(tn);
            }
        }
    }

    /**
     * Creates and reports minor allele frequency by population
     * @param a input alignment
     * @param popMask in REGEX format
     */
    private void calcSNPsFreqByPop(Alignment a, double pDevFromExp) {
        popFreq=new short[2][popNames.size()][a.getSiteCount()];
        popP=new double[popNames.size()][a.getSiteCount()];
        errorRate=new double[a.getSiteCount()];
        int[][] errorCorrCnt=new int[2][a.getSiteCount()];
        for (int s = 0; s < a.getSiteCount(); s++) {
            byte mjb=a.getMajorAllele(s);
            byte mnb=a.getMinorAllele(s);
            byte hetb=IUPACNucleotides.getDegerateSNPByteFromTwoSNPs(mjb, mnb);
//            if((mjb==DataType.UNKNOWN_BYTE)||(mnb==DataType.UNKNOWN_BYTE)) continue;
            for (int t = 0; t < a.getSequenceCount(); t++) {
                if(popOfTaxa[t]<0) continue;
                //should this include homo versus heterozygote
                byte b=a.getBase(t, s);
                if(b==DataType.UNKNOWN_BYTE) continue;
                if(b==mjb) popFreq[0][popOfTaxa[t]][s]+=homoweight;
                if(b==mnb) popFreq[1][popOfTaxa[t]][s]+=homoweight;
                if(b==hetb) {popFreq[0][popOfTaxa[t]][s]++; popFreq[1][popOfTaxa[t]][s]++;}
            }
            for (int pop = 0; pop < popNames.size(); pop++) {
                int pmj, pmn;  //pop major and minor
                popP[pop][s]=Double.NaN;
                if(popFreq[0][pop][s]<popFreq[1][pop][s]) {
                    pmj=popFreq[1][pop][s]; pmn=popFreq[0][pop][s];
                } else {
                    pmj=popFreq[0][pop][s]; pmn=popFreq[1][pop][s];
                }
                int present=pmj+pmn;
                if(present<minCountForDetectingError) continue;
                binomFunc.setNandP(present,expSegregation);
                popP[pop][s] = binomFunc.cdf(pmn);
//                System.out.println(pmj+"\t"+pmn+"\t"+popP[pop][s]);
                double minFreq=(double)pmn/(double)present;
                if((popP[pop][s]<pDevFromExp)&&(minFreq<(expSegregation/minDistortionRatio))) {
//                if((popP[pop][s]<pDevFromExp)) {
                    errorCorrCnt[0][s]+=pmn;
                    errorCorrCnt[1][s]+=present;
                }
            }
         }
        for (int s = 0; s < a.getSiteCount(); s++) {
            double error=errorCorrCnt[0][s];
            double total=errorCorrCnt[1][s];
            if(error<1) error=0.5;
            if(total<1) {errorRate[s]=Double.NaN;}  //if the SNP was never present in biparental pops then the error is undefined
            else {errorRate[s]=error/total;}
        }
    }



    public void reportSNPFrequencyByFamily() {
        //TODO make it saveable to file or stdout
        StringBuilder sb=new StringBuilder("Site\tErrorRate\t");
        for (int pop = 0; pop < popNames.size(); pop++) {
                sb.append(popNames.get(pop)); sb.append("\t");
            }
        System.out.println(sb.toString());
        DecimalFormat df=new DecimalFormat("0.00");
        for (int s = 0; s < a.getSiteCount(); s++) {
            sb=new StringBuilder(a.getSNPID(s)+"\t"+errorRate[s]+"\t");
            for (int pop = 0; pop < popNames.size(); pop++) {
                if(Double.isNaN(popP[pop][s])) {sb.append("NaN");}
                else {sb.append(df.format((double)popFreq[1][pop][s]/(double)(popFreq[0][pop][s]+popFreq[1][pop][s])));}
                sb.append("\t");
            }
            System.out.println(sb.toString());
        }
        
    }

    public void reportDistOfFrequency() {
        //TODO make it saveable to file or stdout
        int[] bins=new int[101];
        for (int pop = 0; pop < popNames.size(); pop++) {
            for (int s = 0; s < a.getSiteCount(); s++) {
                if(Double.isNaN(popP[pop][s])) continue;
                double minFreq=(double)popFreq[1][pop][s]/(double)(popFreq[0][pop][s]+popFreq[1][pop][s]);
                bins[(int)Math.round(100*minFreq)]++;
            }
        }
        System.out.println("BinFreqX100\tNumOfSites");
        for (int i = 0; i < bins.length; i++) {
            System.out.printf("%d\t%d%n",i,bins[i]);
        }
    }
    
    public void reportPercentilesOfErrorRates() {
        //TODO make it saveable to file or stdout
        double[] le=Arrays.copyOf(errorRate, errorRate.length);
        Arrays.sort(le);
        System.out.println("Percentile\tErrorRate");
        for (int i = 0; i < le.length; i+=(le.length/20)) {
            System.out.printf("%.2g\t%.3g%n",((double)i/(double)le.length),le[i]);
        }
    }

    public static void main(String[] args) {
        if(args.length==0) {
            System.out.println("Input format -hmp MapFile -o OutMapFile -popM Z[0-9]{3} -sC 1 -eC 10 -mxE 0.01 -mnD 2.0");
        //    String[] s={"-hmp","/Users/edbuckler/SolexaAnal/GBS/hmp/allfusion_gen_110413t.c+.filt.hmp.txt",
            String[] s={

//            "-hmp","/Users/edbuckler/SolexaAnal/GBS/test/RIL_Marker_data_before_imputation.hmp.txt",
//            "-o","/Users/edbuckler/SolexaAnal/GBS/test/sorghumError.poprep.txt",
//               "-hmp", "/Users/edbuckler/SolexaAnal/GBS/test/h10000mergedNAM282Ames.c9.f97m1t5s10.hmp.txt",
//               "-o", "/Users/edbuckler/SolexaAnal/GBS/test/h10000mergedNAM282Ames.c9.f97m1t5s10.popF.hmp.txt",
//          "-hmp","/Users/edbuckler/SolexaAnal/GBS/hmp/allfusion110611/mergedNAM282Ames.c+.f9m1t5s10.hmp.txt",
//                "-o","/Users/edbuckler/SolexaAnal/GBS/test/mergedNAM282Ames.c+.f9m1t5s10.popF.hmp.txt",
         "-hmp", "/Users/edbuckler/SolexaAnal/GBS/alignment_test_071811/hmp/maize071811.cov10.fT1hLDE1.mgNoHet.c10.hmp.txt",
            "-o", "/Users/edbuckler/SolexaAnal/GBS/alignment_test_071811/test/removesteak10.hmp.txt",
//               "-hmp","/Users/edbuckler/SolexaAnal/GBS/hmp/allfusion110611/DTMA.c+.filt.hmp.txt",
//                "-o","/Users/edbuckler/SolexaAnal/GBS/test/DTMA.c+.filt.popFhmp.txt",
                "-oB", "/Users/edbuckler/SolexaAnal/GBS/test/errorBin.txt",
                "-oE", "/Users/edbuckler/SolexaAnal/GBS/test/errorBySNP.txt",
//                "-popM","SgSBRIL",
                "-popM","Z[0-9]{3}",
                "-sC","10","-eC","10",
                "-mxE","0.01",
                "-mnD","2.0",
                "-mnPLD", "0.5"
                };
            args=s;
            }
        ArgsEngine engine = new ArgsEngine();
        engine.add("-hmp", "-hmpFile", true);
        engine.add("-o", "--outFile", true);
        engine.add("-oE", "--outErrorFile", true);
        engine.add("-oB", "--outBinDistFile", true);
        engine.add("-popM", "--popMask", true);
        engine.add("-popL", "--popList", true);
        engine.add("-sC", "--startChrom", true);
        engine.add("-eC", "--endChrom", true);
        engine.add("-mxE", "--maxError", true);
        engine.add("-mnD", "--minDistortion", true);
        engine.add("-mnPLD", "--minPopLD", true);
        engine.parse(args);
        int start=1, end=1;
        if(engine.getBoolean("-sC")) {start=Integer.parseInt(engine.getString("-sC"));}
        if(engine.getBoolean("-eC")) {end=Integer.parseInt(engine.getString("-eC"));}
        for (int chr = start; chr <= end; chr++) {
            ArrayList<Datum> dList=new ArrayList<Datum>();
            String infile=engine.getString("-hmp").replace("+", ""+chr);
            System.out.println("Reading: "+infile);
            Alignment a=ImportUtils.readFromHapmap(infile);
            System.out.println("Finished Reading: "+infile);
            dList.add(new Datum(infile,a,"alignment"));
            if(engine.getBoolean("-popM")) dList.add(new Datum("PopMask",engine.getString("-popM"),"popmask"));
            BiParentalErrorCorrection theBPER=new BiParentalErrorCorrection();
            if(engine.getBoolean("-o")) theBPER.setOutHapMap(engine.getString("-o").replace("+", ""+chr));
            if(engine.getBoolean("-oE")) theBPER.setOutErrorRate(engine.getString("-oE").replace("+", ""+chr));
            if(engine.getBoolean("-oB")) theBPER.setOutFreqBin(engine.getString("-oB").replace("+", ""+chr));
            if(engine.getBoolean("-mxE")) theBPER.setMaxErrorRate(Double.parseDouble(engine.getString("-mxE")));
            if(engine.getBoolean("-mnD")) theBPER.setMinDistortionRatio(Double.parseDouble(engine.getString("-mnD")));
            if(engine.getBoolean("-mnPLD")) theBPER.setMinMedianLDR2(Double.parseDouble(engine.getString("-mnPLD")));
            theBPER.performFunction(new DataSet(dList,null));
        }
    }

    public void setOutErrorRate(String outErrorRate) {
        this.outErrorRate = outErrorRate;
    }

    public void setOutHapMap(String outHapMap) {
        this.outHapMap = outHapMap;
    }

    public void setOutFreqBin(String outFreqBin) {
        this.outFreqBin = outFreqBin;
    }

    public void setMaxErrorRate(double maxErrorRate) {
        this.maxErrorRate = maxErrorRate;
    }

    public double getMinDistortionRatio() {
        return minDistortionRatio;
    }

    public void setMinDistortionRatio(double minDistortionRatio) {
        this.minDistortionRatio = minDistortionRatio;
    }

    public void setMinMedianLDR2(double minMedianLDR2) {
        this.minMedianLDR2 = minMedianLDR2;
    }




    @Override
    public DataSet performFunction(DataSet input) {
        a=(Alignment)input.getDataOfType(Alignment.class).get(0).getData();
        a=new MutableSimpleAlignment(a);
        double realDist=AlignmentFilterByGBSUtils.getErrorRateForDuplicatedTaxa(a, true,false,false);
        double randomDist=AlignmentFilterByGBSUtils.getErrorRateForDuplicatedTaxa(a, true,true,false);
        System.out.println("Ratio of RandomToReal:"+randomDist/realDist);
        String pm=(String)input.getDataOfType(String.class).get(0).getData();
        classifyTaxaToPops(pm);
        if(this.minMedianLDR2>0) {
            calcLDByPop(a);
            removeLowLDSites(myRemoveUntestedLD);
        }
        calcSNPsFreqByPop(a,0.001);
//        reportSNPFrequencyByFamily();
        reportDistOfFrequency();
        reportPercentilesOfErrorRates();
        removeErrorsInAlignment(1.0);
        removeHighErrorSites(maxErrorRate,myRemoveUntestedError);
        calcSNPsFreqByPop(a,0.001);
//        reportSNPFrequencyByFamily();
        reportDistOfFrequency();
        reportPercentilesOfErrorRates();
        realDist=AlignmentFilterByGBSUtils.getErrorRateForDuplicatedTaxa(a, true,false,true);
        randomDist=AlignmentFilterByGBSUtils.getErrorRateForDuplicatedTaxa(a, true,true,false);
        System.out.println("Ratio of RandomToReal:"+randomDist/realDist);
        if(outHapMap!=null) ExportUtils.writeToHapmap(a, false, this.outHapMap, '\t');
        return null;
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "FilterErrorForBiparental";
    }

    @Override
    public String getToolTipText() {
        return "Filters and esitmates error rates from biparental populations";
    }
}
