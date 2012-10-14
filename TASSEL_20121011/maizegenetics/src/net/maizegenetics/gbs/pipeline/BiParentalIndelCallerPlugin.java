/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.gbs.pipeline;

import cern.jet.random.Binomial;
import edu.cornell.lassp.houle.RngPack.RandomJava;
import java.awt.Frame;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.TreeMap;
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
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;
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
public class BiParentalIndelCallerPlugin extends AbstractPlugin {
    private short[] popOfTaxa;
    private double expSegregation=0.5;
    private double portionOfChromosomeToCompare=0.05;
    int minCntForLD=20;
    private String outHapMap=null;
    private ArrayList<String> popNames;
    private Alignment a;
    private Alignment indelA; //indel recoding of alignment
    private MutableSimpleAlignment outAlign;
    private int errorCntOfNonMissing=0, consistentSet=0, totalSet=0;
    private ArgsEngine engine = new ArgsEngine();
    int start=1, end=1;
    private String infile;

    Binomial binomFunc=new Binomial(5, 0.5, new RandomJava());
    private double minDistortionRatio=2.5;  //this is minimum distortion to be counted as an error
      //e.g. 2.0, means that the allele frequency has to be less the 25% if the expected is 50%
      

    public BiParentalIndelCallerPlugin(){
        super(null, false);
    }

    public BiParentalIndelCallerPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    private Alignment recodeAlignmentByIndels(Alignment snpAlign) {
        MutableSimpleAlignment msa=new MutableSimpleAlignment(snpAlign);
        for (int s = 0; s < msa.getSiteCount(); s++) {
            for (int t = 0; t < msa.getSequenceCount(); t++) {
                byte cb=msa.getBase(t, s);
                if(cb==DataType.UNKNOWN_BYTE) {msa.setBase(t, s, DataType.GAP_BYTE);}
                else {msa.setBase(t, s, (byte)'+');}
            }
        }
        return msa;
    }


   private void setAlignmentByIndels(int outSite, Alignment snpAlign, PredictorSNPOfIndel psi) {
        IdGroup theOutIdGroup=outAlign.getIdGroup();
        for (int t = 0; t < snpAlign.getSequenceCount(); t++) {
            if(snpAlign.getBase(t, psi.snpSite)==psi.snpWithIndel) {
                int outTaxa=theOutIdGroup.whichIdNumber(snpAlign.getIdGroup().getIdentifier(t).getFullName());
   //             System.out.printf("%s %s %n", snpAlign.getTaxaName(t),outAlign.getTaxaName(outTaxa));
                byte cb=outAlign.getBase(outTaxa, outSite);
                if(cb==DataType.GAP_BYTE) {consistentSet++; continue;}
                if(cb!=DataType.UNKNOWN_BYTE) {errorCntOfNonMissing++; continue;}
                outAlign.setBase(outTaxa, outSite, DataType.GAP_BYTE);
                totalSet++;
            }
        }
    }

    private void indentifyIndelsByPop(Alignment snpAlign, Alignment indelAlign) {
        IdGroup idg=snpAlign.getIdGroup();
        int indelSiteCnt=0;
        int maxPhysicalDistance=(int)Math.round(snpAlign.getPositionInLocus(snpAlign.getSiteCount()-1)*portionOfChromosomeToCompare);
        for (int pop = 0; pop < popNames.size(); pop++) {
            int popIndelSiteCnt=0;
            ArrayList<Identifier> keepNames=new ArrayList<Identifier>();
            for (int i = 0; i < popOfTaxa.length; i++) if(popOfTaxa[i]==pop) keepNames.add(idg.getIdentifier(i));
            Identifier[] ids=new Identifier[keepNames.size()];
            keepNames.toArray(ids);
            Alignment snpPA = new MutableSimpleAlignment(FilterAlignment.getInstance(snpAlign, new SimpleIdGroup(ids)));
            Alignment indelPA = new MutableSimpleAlignment(FilterAlignment.getInstance(indelAlign, new SimpleIdGroup(ids)));
     //       AlignmentFilterByGBSUtils.genotypicCountsBySite(snpPA, false, true);
            int[] segSites=AlignmentFilterByGBSUtils.getLowHetSNPs(snpPA, false, -2.0, minCntForLD, 0.25, 2);
            Alignment snpSegAlign = new MutableSimpleAlignment(FilterAlignment.getInstance(snpPA, segSites));
          //  snpSegAlign=new MutableSimpleAlignment(snpSegAlign);
     //       AlignmentFilterByGBSUtils.genotypicCountsBySite(snpSegAlign, false, true);
            int poptests=0;
            for (int indelS = 0; indelS < indelPA.getSiteCount(); indelS++) {
                int indelPhys=indelPA.getPositionInLocus(indelS);
                int tests=1;
                PredictorSNPOfIndel bestPredictor=new PredictorSNPOfIndel(-1,-1);
                TreeMap<Double, PredictorSNPOfIndel> predictorTree=new TreeMap<Double, PredictorSNPOfIndel>();
                for (int snpS = 0; snpS < snpSegAlign.getSiteCount(); snpS++) {
                    int snpPhys=snpSegAlign.getPositionInLocus(snpS);
                    if(Math.abs(snpPhys-indelPhys)>maxPhysicalDistance) continue;
                    tests++;
                    PredictorSNPOfIndel result=testIndelWSNP(indelPA, indelS, snpSegAlign, snpS);

                    if(result==null) continue;
               //     System.out.println(Arrays.toString(result));
                    if(predictorTree.size()==0) predictorTree.put(result.pObsMinor, result);
                    if(predictorTree.lastKey()>result.pObsMinor) predictorTree.put(result.pObsMinor, result);
                    if(predictorTree.size()>10) predictorTree.remove(predictorTree.lastKey());
                }
                if((predictorTree.isEmpty()==false)&&(predictorTree.firstKey()<(0.05/(double)tests))) {
                     popIndelSiteCnt++;
//                     System.out.println(predictorTree.keySet().toString());
              //       setAlignmentByIndels(indelS, snpSegAlign, predictorTree.firstEntry().getValue());
                     for(PredictorSNPOfIndel e: predictorTree.values()) {
                         setAlignmentByIndels(indelS, snpSegAlign, e);
                     }
                }
                
                poptests+=tests;
            }
            indelSiteCnt+=popIndelSiteCnt;
            System.out.printf("Pop:%s Tests:%d Indels:%d Total: %d %d %d %d %n",ids[0].getName(),poptests, popIndelSiteCnt,
                    indelSiteCnt, errorCntOfNonMissing, consistentSet, totalSet);
        }
    }



    private PredictorSNPOfIndel testIndelWSNP(Alignment indelAlign, int indelSite, Alignment snpAlign, int snpSite) {
        byte majorAllele=snpAlign.getMajorAllele(snpSite);
        byte minorAllele=snpAlign.getMinorAllele(snpSite);
        int majorPresentCnt=0, minorPresentCnt=0, majorAbsenseCnt=0, minorAbsenseCnt=0;
        PredictorSNPOfIndel psi=null;
        for (int i = 0; i < indelAlign.getSequenceCount(); i++) {
            if(indelAlign.getBase(i, indelSite)==(byte)'+') {
                byte cb=snpAlign.getBase(i, snpSite);
                if(majorAllele==cb) {majorPresentCnt++;}
                else if(minorAllele==cb) {minorPresentCnt++;}
            } else {
                byte cb=snpAlign.getBase(i, snpSite);
                if(majorAllele==cb) {majorAbsenseCnt++;}
                else if(minorAllele==cb) {minorAbsenseCnt++;}
            }
        }
        int sumMajor=majorAbsenseCnt+majorPresentCnt;
        int sumMinor=minorAbsenseCnt+minorPresentCnt;
        int sumSNP=sumMajor+sumMinor;
        int sumPresent = minorPresentCnt + majorPresentCnt;
//        System.out.printf("%d v %d: %d %d %d %d %n",indelSite, snpSite, majorPresentCnt, minorPresentCnt, majorAbsenseCnt, minorAbsenseCnt);
//        System.out.printf("%d v %d: MAF: %g %g %n", indelSite, snpSite, indelAlign.getMinorAlleleFrequency(indelSite),snpAlign.getMinorAlleleFrequency(snpSite));
        if(sumPresent<minCntForLD) return null;
        double tag_rat=(minorPresentCnt<majorPresentCnt)?(double)minorPresentCnt/(double)majorPresentCnt:(double)majorPresentCnt/(double)minorPresentCnt;
        if(tag_rat>0.01) return null;
        double minorProb = (minorPresentCnt<majorPresentCnt)?(double)sumMinor/(double)sumSNP:(double)sumMajor/(double)sumSNP;
        try {
            psi=new PredictorSNPOfIndel(indelSite, snpSite);
            binomFunc.setNandP(sumPresent,expSegregation);
            psi.sumPresent= sumPresent;
            psi.snpWithIndel= (minorPresentCnt<majorPresentCnt)?minorAllele:majorAllele;
            psi.cntPredictorSNP = (minorPresentCnt<majorPresentCnt)?minorPresentCnt:majorPresentCnt;
            binomFunc.setNandP(sumPresent,minorProb);
            psi.pObsMinor = (minorPresentCnt<majorPresentCnt)?binomFunc.cdf(minorPresentCnt):binomFunc.cdf(majorPresentCnt);
            psi.position=snpAlign.getPositionInLocus(snpSite);
  //          if(result<0.000001) System.out.printf("TS: %d %d %d %d %g %g %n",minorPresentCnt, majorPresentCnt, mn, mj, result, tag_rat);
        } catch (Exception e) {
            System.err.println("Error in the BinomialDistributionImpl");
            System.err.printf("Error: %d %d %d %d %n",majorPresentCnt, minorPresentCnt, majorAbsenseCnt, minorAbsenseCnt);
        }
        return psi;
    }

    private class PredictorSNPOfIndel {
        int indelSite, snpSite;
        int sumPresent=-1, cntPredictorSNP=-1;
        byte snpWithIndel=DataType.UNKNOWN_BYTE;
        double pObsMinor=Double.MAX_VALUE;
        int position;

        public PredictorSNPOfIndel(int indelSite, int snpSite) {
            this.indelSite=indelSite;
            this.snpSite=snpSite;
        }

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


    public void setParameters(String[] args){
         if(args.length==0) {
            System.out.println("Input format:\n -hmp MapFile\n -o OutMapFile\n -popM Z[0-9]{3}\n -sC 1\n -eC 10\n -mxE 0.01\n -mnD 2.0\n -mnPLD 0.5\n");
        //    String[] s={"-hmp","/Users/edbuckler/SolexaAnal/GBS/hmp/allfusion_gen_110413t.c+.filt.hmp.txt",
            String[] s={
          "-hmp","/Users/edbuckler/SolexaAnal/GBS/hmp/allfusion110611/mergedNAM282Ames.c+.f9m1t5s10.hmp.txt",
                "-o","/Users/edbuckler/SolexaAnal/GBS/test/mergedNAM282Ames.c+.f9m1t5s10.popF.hmp.txt",
                "-oB", "/Users/edbuckler/SolexaAnal/GBS/test/errorBin.txt",
                "-oE", "/Users/edbuckler/SolexaAnal/GBS/test/errorBySNP.txt",
//                "-popM","SgSBRIL",
                "-popM","Z[0-9]{3}",
                "-sC","9","-eC","10",
                };
            args=s;
            }

        engine.add("-hmp", "-hmpFile", true);
        engine.add("-o", "--outFile", true);
        engine.add("-popM", "--popMask", true);
        engine.add("-sC", "--startChrom", true);
        engine.add("-eC", "--endChrom", true);
        engine.parse(args);
        if(engine.getBoolean("-sC")) {start=Integer.parseInt(engine.getString("-sC"));}
        if(engine.getBoolean("-eC")) {end=Integer.parseInt(engine.getString("-eC"));}
        infile = engine.getString("-hmp");
        outHapMap = engine.getString("-o");
        performFunction(null);
    }

    public static void main(String[] args) {
        String rootIn="/Users/edbuckler/SolexaAnal/GBS/build110816/hmp/h5kmaize110816";
        String rootOut="/Users/edbuckler/SolexaAnal/GBS/build110816/test/h5kmaize110816";
        int sC=10;
        int eC=10;
        String[] nargs = new String[] {
            "-hmp", rootIn+".cov10.fT1E1pLD.mgNoHet.imp.c+.hmp.txt",
            "-o", rootOut+".cov10.fT1E1pLD.mgNoHet.gap.imp.c+.hmp.txt",
            "-popM","Z[0-9]{3}",
             "-sC", ""+sC,  // Start chromosome
            "-eC", ""+eC // End chromosome
        };
        args=nargs;
        BiParentalIndelCallerPlugin testClass = new BiParentalIndelCallerPlugin();
        testClass.setParameters(args);
        testClass.performFunction(null);
    }

    public void setOutHapMap(String outHapMap) {
        this.outHapMap = outHapMap;
    }

    public double getMinDistortionRatio() {
        return minDistortionRatio;
    }

    public void setMinDistortionRatio(double minDistortionRatio) {
        this.minDistortionRatio = minDistortionRatio;
    }


    @Override
    public DataSet performFunction(DataSet input) {
        while (start <= end) {
            int chr = start;
            ArrayList<Datum> dList=new ArrayList<Datum>();
            String currFile =infile.replace("+", ""+chr);
            System.out.println("Reading: "+currFile);
            Alignment a=ImportUtils.readFromHapmap(currFile);
            System.out.println("Finished Reading: "+currFile);
            System.out.printf("Taxa: %d Sites: %d %n", a.getSequenceCount(), a.getSiteCount());
            dList.add(new Datum(currFile,a,"alignment"));
            dList.add(new Datum("PopMask",engine.getString("-popM"),"popmask"));
            BiParentalIndelCallerPlugin theBPER=new BiParentalIndelCallerPlugin();
            if(engine.getBoolean("-o")) theBPER.setOutHapMap(engine.getString("-o").replace("+", ""+chr));
            theBPER.filter(new DataSet(dList,null));
            start++;
        }
        return null;
    }

    private void filter(DataSet input){
        List<Datum> b=input.getDataOfType(Alignment.class);
        Datum c = b.get(0);
        Object d = c.getData();
        a=(Alignment)d;
        a=new MutableSimpleAlignment(a);
        outAlign=new MutableSimpleAlignment(a);
        String pm=(String)input.getDataOfType(String.class).get(0).getData();
        classifyTaxaToPops(pm);
        this.indelA=recodeAlignmentByIndels(a);
        indentifyIndelsByPop(a,indelA);
        if(outHapMap!=null) ExportUtils.writeToHapmap(outAlign, false, this.outHapMap, '\t');
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
