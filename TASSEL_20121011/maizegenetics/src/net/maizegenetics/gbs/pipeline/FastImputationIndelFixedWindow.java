/*
 * BasicImputation
 */
package net.maizegenetics.gbs.pipeline;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.TreeMap;
import net.maizegenetics.gbs.util.MutableSimpleAlignment;
import net.maizegenetics.gbs.util.OpenBitSet;
import net.maizegenetics.gbs.util.TBitAlignmentTest;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.GdpdmBLOBUtils;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.Pack1Alignment;
import net.maizegenetics.pal.datatype.DataType;
import net.maizegenetics.pal.datatype.IUPACNucleotides;

/**
 * An extremely fast imputation approach that uses KNN, but calculates the matrix
 * with a 64 bit window, and it does all operations with bits.  This has a fixed window size
 * unlike the other BDI imputation approaches.
 *
 *All heterozygous sets are to missing.
 *
 * Definitely a work in progress.
 * Perhaps we should do ratio to major/minor
 * Perhaps we should scan for massive nearly identical regions
 * Need to decide whether P or Length/Identity approach is better
 *
 * @author ed
 */
public class FastImputationIndelFixedWindow {
   // static byte UNKNOWN=(byte)DataType.UNKNOWN_CHARACTER;
    int[][] matchInWin, diffInWin;
    float[] presentProp;  //Present proportion
    int[] presentCntForSites;
    int[] presentCntForTaxa;
    int[] hetCntForTaxa;
    boolean[] highHet;
    long totalPresent=0;
    
    TBitAlignmentTest anchorBitAlignment=null;
    Pack1Alignment anchorAlignmentReal=null;
    MutableSimpleAlignment impAlign=null;
    static int windowSize=64*64;  //left and right distance of 64bp window, so window 2 is 2+1+2=5 or 320bp
    static int minLengthOfMatch=50;
    static double minIdentity=0.99;  //std 0.99, 0.95 gave some pretty good results also
    static int minCountNNperLine=4;  //std 4
    static int minCountNNperSite=2;  //std 2
    static int segments=1;
    static double majorityRule=0.76;
    static double maxLineErrorRate=0.025; //std 0.01
    static boolean imputeGaps=false;
    static boolean maskHighHets=true;
    static double maxHetStatic=0.01; //HetStat = hetSiteCnt/covSiteCnt :  empirically derived

//    static int minimumMajorityRule=4;  //NN within  is used, and tie are not imputed


    public FastImputationIndelFixedWindow(Alignment a) {
        this.anchorAlignmentReal=(Pack1Alignment)a;
        System.out.print("Converting To TBit: ");
        this.anchorBitAlignment=new TBitAlignmentTest(a);
        System.out.println("Done ");
        double avgMissing=calcPropMissingByLine();
        System.out.println("Average Missing:"+avgMissing);
        avgMissing=calcPropMissing();
        System.out.println("Average Missing (site and taxa):"+avgMissing+" totalPresent"+totalPresent);
        System.out.println("Average MAF:"+avgMinorAlleleFrequency());
//        reportTaxaStats();
        AlignmentFilterByGBSUtils.getCoverage_MAF_F_Dist(anchorBitAlignment, false);
        double realDist=AlignmentFilterByGBSUtils.getErrorRateForDuplicatedTaxa(anchorBitAlignment, true,false,true);
        double randomDist=AlignmentFilterByGBSUtils.getErrorRateForDuplicatedTaxa(anchorBitAlignment, true,true,false);
        System.out.println("Ratio of RandomToReal:"+randomDist/realDist);
        System.out.println("Creating mutable alignment");
        impAlign=new MutableSimpleAlignment(a);
        for(int i=4096; i>=1024; i/=2) {
//            segments=i;  
            System.out.println("Starting imputation");
 //           windowSize=((a.getSiteCount()/64)/segments)*64;
            windowSize=i;
            int offset=0;
            System.out.println("Window size:"+windowSize+" offset:"+offset);
            imputeBySiteJump(windowSize,0,minLengthOfMatch,minIdentity);
//            realDist=AlignmentFilterByGBSUtils.getErrorRateForDuplicatedTaxa(impAlign, true,false,false);
//            randomDist=AlignmentFilterByGBSUtils.getErrorRateForDuplicatedTaxa(impAlign, true,true,false);
//            System.out.println("Ratio of RandomToReal:"+randomDist/realDist);
            System.out.println("Window size:"+windowSize+" offset:"+offset);
            AlignmentFilterByGBSUtils.getCoverage_MAF_F_Dist(impAlign, false);
            offset=windowSize/2;
            System.out.println("Window size:"+windowSize+" offset:"+offset);
            imputeBySiteJump(windowSize,offset,minLengthOfMatch,minIdentity);
            System.out.println("Window size:"+windowSize+" offset:"+offset);
            AlignmentFilterByGBSUtils.getCoverage_MAF_F_Dist(impAlign, false);
            System.out.println("Xratio: "+countXRateInImputation());
        }
    }

    private double avgMinorAlleleFrequency() {
        long totmj=0, totmn=0;
        for(int i=0; i<anchorBitAlignment.getSequenceCount(); i++) {
            totmj+=anchorBitAlignment.getTaxaBitsNoClone(i, 0).cardinality();
            totmn+=anchorBitAlignment.getTaxaBitsNoClone(i, 1).cardinality();
        }
        double theMAF=(double)totmn/(double)(totmj+totmn);
        return theMAF;
    }

    private double countXRateInImputation() {
        long xCnt=0, knownCnt=0;
        for(int i=0; i<impAlign.getSequenceCount(); i++) {
            for (int j = 0; j < impAlign.getSiteCount(); j++) {
                if(impAlign.getBase(i, j)!=DataType.UNKNOWN_BYTE) {
                    knownCnt++;
                    if(impAlign.getBase(i, j)==(byte)'X') xCnt++;
                }

            }
        }
        double theXRate=(double)xCnt/(double)knownCnt;
        return theXRate;
    }

    public void reportTaxaStats() {
        double sites=anchorBitAlignment.getSiteCount();
        for(int t=0; t<anchorBitAlignment.getSequenceCount(); t++) {
            System.out.printf("%s %d %.5g %d %.5g %s %n",anchorBitAlignment.getIdGroup().getIdentifier(t).getFullName(), presentCntForTaxa[t],
                    (double)presentCntForTaxa[t]/sites, hetCntForTaxa[t], (double)hetCntForTaxa[t]/sites, highHet[t]);
        }
    }

    private double calcPropMissing() {
        presentCntForTaxa=new int[anchorBitAlignment.getSequenceCount()];
        presentCntForSites=new int[anchorBitAlignment.getSiteCount()];
        hetCntForTaxa=new int[anchorBitAlignment.getSequenceCount()];
        highHet=new boolean[anchorBitAlignment.getSequenceCount()];
        totalPresent=0;
        IUPACNucleotides iupac=new IUPACNucleotides();
        double sd=anchorBitAlignment.getSiteCount();
        for(int t=0; t<anchorBitAlignment.getSequenceCount(); t++) {
            for (int s = 0; s < anchorBitAlignment.getSiteCount(); s++) {
                byte cb=anchorBitAlignment.getBase(t, s);
                if((cb!=DataType.UNKNOWN_BYTE)&&(cb!=DataType.GAP_BYTE)) {
                    presentCntForTaxa[t]++;
                    presentCntForSites[s]++;
                    if(iupac.isHeterozygote(cb)) hetCntForTaxa[t]++;
                    totalPresent++;
                }
            }
            double covProp=(double)presentCntForTaxa[t]/sd;
            double hetProp=(double)hetCntForTaxa[t]/presentCntForTaxa[t];
            highHet[t]=(hetProp>maxHetStatic)?true:false;

        }
        return (double)totalPresent/((double)anchorBitAlignment.getSequenceCount()*(double)anchorBitAlignment.getSiteCount());
    }

    private double calcPropMissingByLine() {
        presentProp=new float[anchorBitAlignment.getSequenceCount()];
        double avgMissing=0;
        for(int i=0; i<anchorBitAlignment.getSequenceCount(); i++) {
            long present=OpenBitSet.unionCount(anchorBitAlignment.getTaxaBitsNoClone(i, 0), anchorBitAlignment.getTaxaBitsNoClone(i, 1));
            presentProp[i]=(float)present/(float)anchorBitAlignment.getSiteCount();
            avgMissing+=presentProp[i];
        }
        return avgMissing/(double)anchorBitAlignment.getSequenceCount();
    }

      private void imputeBySiteJump(int window, int offset, int minLength, double minIdentity) {
        long time=System.currentTimeMillis();
        int knownSNPs = 0, unknownSNPs = 0, imputedSNPs=0;
        int numSeqs=anchorBitAlignment.getSequenceCount();
        System.out.println("Initial matrix created in "+(System.currentTimeMillis()-time));
        time=System.currentTimeMillis();
        int corrCnt=0, wrongCnt=0, gapCnt=0;
        int imputableLines=0, imputedWithHighError=0;
        TreeMap<Double,Integer> lenBestLine=new TreeMap<Double,Integer>(Collections.reverseOrder());
        for (int b = 0+offset; b < anchorBitAlignment.getSiteCount()-window; b+=window) {
     //       currWord=b>>6;
            int endBase=b+window;
            initHapLengths(b,endBase);
            double rate=(double)unknownSNPs/(double)(System.currentTimeMillis()-time);
            System.out.println("Imputed base:" + b + " known:" + knownSNPs + " unknownSNPs:" + unknownSNPs+
                    " imputed:"+imputedSNPs+" Rate:"+rate);
            for (int i = 0; i < numSeqs; i++) {
                if(maskHighHets&&highHet[i]) continue;
                lenBestLine.clear();
                for (int j = 0; j < numSeqs; j++) {
                    if(i==j) continue;
                    if(maskHighHets&&highHet[j]) continue;
                    int sum=matchInWin[i][j]+diffInWin[i][j];
                    double identity=(double)matchInWin[i][j]/(double)sum;
                    if((sum>minLength)&&(identity>minIdentity)) lenBestLine.put(identity, j);
                }
                if(lenBestLine.size()<minCountNNperLine) continue;
                imputableLines++;
                if(i==0) System.out.println(" cnt"+lenBestLine.size());
                byte[] calls=consensusCallsMulti(anchorAlignmentReal, i, b, endBase, lenBestLine, false, majorityRule, minCountNNperSite, true);
//                System.out.println(Arrays.toString(calls));
//                byte[] calls=consensusCallBit(anchorBitAlignment, i, b, endBase, lenBestLine,
  //                      false, majorityRule, minCountNNperSite, true, imputeGaps);  //determine error while ignoring known
                int[] callError=compareCallsWithActual(anchorBitAlignment, b, i, calls);
                double lineError=(double)callError[1]/(double)(callError[1]+callError[0]);
                if(lineError>maxLineErrorRate) {
                    imputedWithHighError++;
//                    System.out.println("high error line:"+anchorBitAlignment.getTaxaName(i)+" Error:"+lineError);
                    continue;}
                calls=consensusCallsMulti(anchorAlignmentReal, i, b, endBase, lenBestLine, false, majorityRule, minCountNNperSite, false);
//                calls=consensusCallBit(anchorBitAlignment, i, b, endBase, lenBestLine,
//                        false, majorityRule, minCountNNperSite, false, imputeGaps);  //recall with known; currently conflict are set to unknown
                setBaseInImputedAlignment(impAlign, b, i, calls);
  //              if(b>3000) System.out.println(Arrays.toString(callError));
                corrCnt+=callError[0];
                wrongCnt+=callError[1];
                gapCnt+=callError[3];
                double errorRate=(double)wrongCnt/(double)(corrCnt+wrongCnt);

                if(i%100==0) System.out.printf("%d-%d %d R: %d W: %d G: %d ErrorRate: %g Imputable:%d ImputedHighError:%d %n", b, endBase, i,
                        corrCnt, wrongCnt, gapCnt, errorRate, imputableLines, imputedWithHighError);
            }
        double errorRate=(double)wrongCnt/(double)(corrCnt+wrongCnt);
        System.out.printf("R: %d W: %d ErrorRate: %g %n",corrCnt, wrongCnt,errorRate);
        }
    }

   private int[] compareCallsWithActual(Alignment a, int startBase, int taxon, byte[] calls) {
       int[] result=new int[4];  //agree[0], disagree[1], noncomparison[2], gaps[3]
       byte alignB=DataType.UNKNOWN_BYTE;
       for (int s = 0; s < calls.length; s++) {
           alignB=a.getBase(taxon, s+startBase);
           if(calls[s]==DataType.GAP_BYTE) {result[3]++;}
           if(calls[s]==DataType.UNKNOWN_BYTE) {result[2]++;}
           else if(alignB==DataType.UNKNOWN_BYTE) {result[2]++;}
           else if(alignB==calls[s]) {result[0]++;}
           else {result[1]++;}
       }
       return result;
   }

   private void setBaseInImputedAlignment(MutableSimpleAlignment a, int startBase, int taxon, byte[] calls) {
       for (int s = 0; s < calls.length; s++) {
           if(calls[s]!=DataType.UNKNOWN_BYTE) {
               byte cb=a.getBase(taxon, s+startBase);
               if(cb==DataType.UNKNOWN_BYTE) {a.setBase(taxon, s+startBase, calls[s]);}
               else if(cb!=calls[s]) {a.setBase(taxon, s+startBase, (byte)'X');}
           }
       }
   }


    private byte[] consensusCalls(Alignment a, int taxon, int startBase, int endBase, TreeMap<Double,Integer> taxa,
            boolean callhets, double majority, int minCount, boolean ignoreKnownBases) {
        short[][] siteCnt=new short[2][a.getSiteCount()];
        int[] taxaIndex=new int[taxa.size()];
        ArrayList<Integer> taxaList=new ArrayList(taxa.values());
        for (int t = 0; t < taxaIndex.length; t++) {
            taxaIndex[t]=taxaList.get(t);
        }
        byte[] calls=new byte[endBase-startBase];
        Arrays.fill(calls, DataType.UNKNOWN_BYTE);
        for (int s = startBase; s < endBase; s++) {
            if(ignoreKnownBases==false) calls[s-startBase]=a.getBase(taxon,s);
            byte mj=a.getMajorAllele(s);
            byte mn=a.getMinorAllele(s);
            byte[] snpValue={mj,mn};
            byte het=IUPACNucleotides.getDegerateSNPByteFromTwoSNPs(snpValue);
            for (int t = 0; t < taxaIndex.length; t++) {
                byte ob=a.getBase(taxaIndex[t], s);
                if(ob==DataType.UNKNOWN_BYTE) continue;
                if(ob==mj) {siteCnt[0][s]++;}
                else if(ob==mn) {siteCnt[1][s]++;}
                else if(ob==het) {siteCnt[0][s]++; siteCnt[1][s]++;}
            }
            int totalCnt=siteCnt[0][s]+siteCnt[1][s];
            if(totalCnt<minCount) continue;  //no data leave missing
            if((double)siteCnt[0][s]/(double)totalCnt>majority) {calls[s-startBase]=mj;}
            else if((double)siteCnt[1][s]/(double)totalCnt>majority) {calls[s-startBase]=mn;}
            else if(callhets) {calls[s-startBase]=het;}
        }
 //       System.out.println("Byt:"+Arrays.toString(siteCnt[0]));
        return calls;
    }

     private byte[] consensusCallsMulti(Pack1Alignment a, int taxon, int startBase, int endBase, TreeMap<Double,Integer> taxa,
            boolean callhets, double majority, int minCount, boolean ignoreKnownBases) {
//        short[][] siteCnt=new short[2][a.getSiteCount()];
        int[] taxaIndex=new int[taxa.size()];
        ArrayList<Integer> taxaList=new ArrayList(taxa.values());
        for (int t = 0; t < taxaIndex.length; t++) {
            taxaIndex[t]=taxaList.get(t);
        }
        byte[] calls=new byte[endBase-startBase];
        Arrays.fill(calls, DataType.UNKNOWN_BYTE);
        for (int s = startBase; s < endBase; s++) {
            if(ignoreKnownBases==false) calls[s-startBase]=a.getBase(taxon,s);
//            byte mj=a.getMajorAllele(s);
//            byte mn=a.getMinorAllele(s);
//            byte[] snpValue={mj,mn};
//            byte het=IUPACNucleotides.getDegerateSNPByteFromTwoSNPs(snpValue);
            short[] baseCnt=new short[Byte.MAX_VALUE];
            int totalCnt=0;
            for (int t = 0; t < taxaIndex.length; t++) {
                byte ob=a.getBase(taxaIndex[t], s);
                if(ob==DataType.UNKNOWN_BYTE) continue;
                totalCnt++;
                baseCnt[ob]++;
//                if(ob==mj) {siteCnt[0][s]++;}
//                else if(ob==mn) {siteCnt[1][s]++;}
//                else if(ob==het) {siteCnt[0][s]++; siteCnt[1][s]++;}
            }
           // int totalCnt=siteCnt[0][s]+siteCnt[1][s];
            if(totalCnt<minCount) continue;  //no data leave missing
            for(byte tb: GdpdmBLOBUtils.bases) {
                if((double)baseCnt[tb]/(double)totalCnt>majority) {calls[s-startBase]=tb; break;}
            }
        }
 //       System.out.println("Byt:"+Arrays.toString(siteCnt[0]));
        return calls;
    }

    private byte[] consensusCallBit(Alignment a, int taxon, int startBase, int endBase, TreeMap<Double,Integer> taxa,
            boolean callhets, double majority, int minCount, boolean ignoreKnownBases, boolean imputeGaps) {
        int[] taxaIndex=new int[taxa.size()];
        ArrayList<Integer> taxaList=new ArrayList(taxa.values());
        for (int t = 0; t < taxaIndex.length; t++) {
            taxaIndex[t]=taxaList.get(t);
        }
        short[][] siteCnt=new short[2][endBase-startBase];
        double[] sumExpPresent=new double[endBase-startBase];
        int[] sumNNxSitePresent=new int[endBase-startBase];
        int sumNNPresent=0;
        for (int t = 0; t < taxaIndex.length; t++) sumNNPresent+=presentCntForTaxa[taxaIndex[t]];
        byte[] calls=new byte[endBase-startBase];
        Arrays.fill(calls, DataType.UNKNOWN_BYTE);
        for (int alignS = startBase; alignS < endBase; alignS+=64) {
            int currWord=alignS/64;
            int callSite=alignS-startBase;
            for (int t = 0; t < taxaIndex.length; t++) {
                long bmj=anchorBitAlignment.getTaxaBitsNoClone(taxaIndex[t], 0).getBits()[currWord];
                long bmn=anchorBitAlignment.getTaxaBitsNoClone(taxaIndex[t], 1).getBits()[currWord];
                int cs=callSite;
                for (int j = 0; j < 64; j++) {
                    boolean presentFlag=false;
                    if((bmj & 0x01)!=0) {siteCnt[0][cs]++; presentFlag=true;}
                    bmj=bmj>>1;
                    if((bmn & 0x01)!=0) {siteCnt[1][cs]++; presentFlag=true;}
                    bmn=bmn>>1;
                    sumExpPresent[cs]+=presentProp[taxaIndex[t]];
                    if(presentFlag) sumNNxSitePresent[cs]++;
                    cs++;          
                }
            }
        }
 //       System.out.println("Bit:"+Arrays.toString(siteCnt[0]));
        for (int alignS = startBase; alignS < endBase; alignS++) {
            int callS=alignS-startBase;
            byte ob=a.getBase(taxon,alignS);
            if(ignoreKnownBases) ob=DataType.UNKNOWN_BYTE;
            calls[callS]=ob;
            byte mj=a.getMajorAllele(alignS);
            byte mn=a.getMinorAllele(alignS);
            int totalCnt=siteCnt[0][callS]+siteCnt[1][callS];
            if(totalCnt<minCount) continue;  //no data leave missing
            if((double)siteCnt[0][callS]/(double)totalCnt>majority) {
                if((ob!=DataType.UNKNOWN_BYTE)&&(ob!=mj)) {calls[callS]=DataType.UNKNOWN_BYTE;}
                else {calls[callS] = mj;}
            }
            else if((double)siteCnt[1][callS]/(double)totalCnt>majority) {
                if((ob!=DataType.UNKNOWN_BYTE)&&(ob!=mn)) {calls[callS]=DataType.UNKNOWN_BYTE;}
                else {calls[callS] = mn;}
            }
            else if(callhets) {
                byte[] snpValue={mj,mn};
                byte het=IUPACNucleotides.getDegerateSNPByteFromTwoSNPs(snpValue);
                calls[callS]=het;
            }
        }
        return calls;
    }


    private void initHapLengths(int startSite, int endSite) {
        matchInWin=new int[anchorBitAlignment.getSequenceCount()][anchorBitAlignment.getSequenceCount()];
        diffInWin=new int[anchorBitAlignment.getSequenceCount()][anchorBitAlignment.getSequenceCount()];
 //       pOfMatch=new float[anchorBitAlignment.getSequenceCount()][anchorBitAlignment.getSequenceCount()];
//        int currWord=initialSite>>6;
        int maxWord=anchorBitAlignment.getTaxaBitsNoClone(0, 0).getNumWords()-1;
        int startWord=startSite>>6;
        int endWord=endSite>>6;
        if(endWord>maxWord) endWord=maxWord;
//        int startWord=(currWord-windowSize<0)?0:currWord-windowSize;
//        int endWord=(currWord+windowSize>=maxWord)?maxWord-1:currWord+windowSize;
        System.out.printf("Start Site:%d Word:%d End Site:%d Word:%d %n",startSite, startWord, endSite, endWord);
        int[] bins=new int[101];
        int[] countClose=new int[anchorBitAlignment.getSequenceCount()];
        for(int i=0; i<anchorBitAlignment.getSequenceCount(); i++) {
            long[] imj=anchorBitAlignment.getTaxaBitsNoClone(i, 0).getBits();
            long[] imn=anchorBitAlignment.getTaxaBitsNoClone(i, 1).getBits();
            for(int j=0; j<i; j++) {
                long[] jmj=anchorBitAlignment.getTaxaBitsNoClone(j, 0).getBits();
                long[] jmn=anchorBitAlignment.getTaxaBitsNoClone(j, 1).getBits();
                int same=0, diff=0, hets=0;
                int minorSame=0;
                for (int w = startWord; w <= endWord; w++) {
                    long ihetMask=~(imj[w]&imn[w]);
                    long jhetMask=~(jmj[w]&jmn[w]);
                    long imjnh=imj[w]&ihetMask;
                    long imnnh=imn[w]&ihetMask;
                    long jmjnh=jmj[w]&jhetMask;
                    long jmnnh=jmn[w]&jhetMask;
                    same+=Long.bitCount(imjnh&jmjnh)+Long.bitCount(imnnh&jmnnh);
                    diff+=Long.bitCount(imjnh&jmnnh)+Long.bitCount(imnnh&jmjnh);
//                    same+=Long.bitCount(imj[w]&jmj[w])+Long.bitCount(imn[w]&jmn[w]);
//                    diff+=Long.bitCount(imj[w]&jmn[w])+Long.bitCount(imn[w]&jmj[w]);
//                    hets+=Long.bitCount(mj[w]&mn[w])+Long.bitCount(jmn[w]&jmj[w]);
                    minorSame+=Long.bitCount(imn[w]&jmn[w]);
                }
//               matchInWin[j][i]=matchInWin[i][j]=same;
//               diffInWin[j][i]=diffInWin[i][j]=diff;
               matchInWin[j][i]=matchInWin[i][j]=same-hets;
               diffInWin[j][i]=diffInWin[i][j]=diff-hets;
               int sum=same+diff-(2*hets);
               int majorSame=same-minorSame;
               int total=same+diff;
               double div=(double)diffInWin[j][i]/(double)sum;
               double p=1.0;
               if((div<0.02)&&(sum>100)) {
                   countClose[i]++;
                   countClose[j]++;
               }
 //                if (div<0.03) System.out.printf("%s %s %d %d %d %g %n",anchorBitAlignment.getTaxaName(i),anchorBitAlignment.getTaxaName(j),matchInWin[j][i],diffInWin[j][i],sum, div);
//                   bins[(int)((diff*100)/sum)]++;
//               }
           }
        }
//        for (int i = 0; i < bins.length; i++) {
//            System.out.printf("%d %d %d%n",initialSite,i,bins[i]);
//
//        }
//        for (int i = 0; i < countClose.length; i++) {
//            System.out.printf("%d %s %d%n",initialSite, anchorBitAlignment.getTaxaName(i),countClose[i]);
//
//        }
    }


    public void writeAlignment(String outfile) {
        ExportUtils.writeToHapmap(impAlign, false, outfile, '\t');
    }

    public static void main(String[] args) {
	System.out.println("Running main method in FastImputation");

        String[] dargs={
            "-hmp", "/Users/edbuckler/SolexaAnal/GBS/build110816/hmp/h5Kmaize110816.cov10.fT1E1pLD.mgNoHet.c10.hmp.txt",
            "-o", "/Users/edbuckler/SolexaAnal/GBS/build110816/test/h5Kmaize110816.cov10.fT1E1pLD.mgNoHet.imp.c10.hmp.txt",
            "sC","10",
            "eC","10"
        };
        if(args.length==0) {
            args=dargs;
        }
        int sC=Integer.parseInt(args[5]);
        int eC=Integer.parseInt(args[7]);
        String outfile, anchorMapFile;
        for (int cNum = sC; cNum <= eC; cNum++) {
            outfile=args[3].replace("+", ""+cNum);
            anchorMapFile=args[1].replace("+", ""+cNum);
            System.out.println("Reading "+anchorMapFile);
            Alignment a=ImportUtils.readFromHapmap(anchorMapFile);
            System.out.printf("Read Alignment with %d taxa and %d sites %n",a.getSequenceCount(), a.getSiteCount());
           // System.out.println("p1a:"+AlignmentUtils.getSequenceString(a, 1));
           // a=new TBitAlignmentTest(a);
            //System.out.println("TBA:"+AlignmentUtils.getSequenceString(a, 1));
            FastImputationIndelFixedWindow fi=new FastImputationIndelFixedWindow(a);
            System.out.println("Writing "+outfile);
            fi.writeAlignment(outfile);
        }
    }
}
