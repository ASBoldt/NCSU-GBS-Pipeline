/*
 * BasicImputation
 */
package net.maizegenetics.gbs.pipeline;

import java.util.Collections;
import java.util.TreeMap;
import net.maizegenetics.gbs.util.MutableSimpleAlignment;
import net.maizegenetics.gbs.util.TBitAlignmentTest;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AlignmentUtils;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.datatype.DataType;

/**
 * An extremely fast imputation approach that uses KNN, but calculates the matrix
 * with a 64 bit window, and it does all operations with bits.  There seem to be
 * some current bugs with edge effects (i.e. the beginning of the sequence), and possibly
 * with favor common alleles too much.
 *
 * This does class does not use a fixed window.
 * Definitely a work in progress.
 *
 * @author ed
 */
public class FastImputationBit {
    static byte UNKNOWN=(byte)DataType.UNKNOWN_CHARACTER;
    int[][] leftHapLen, rightHapLen, totalLen;
    int currSplitSite; //it is included on the right side
    int[][][] currDistMat;

    static int maxMismatch=0;
    static int maxMatchCalc=300;
    TBitAlignmentTest anchorAlignment=null;
    TBitAlignmentTest impAlign2=null;
    MutableSimpleAlignment impAlign=null;

    public FastImputationBit(Alignment a) {
        this.anchorAlignment=(TBitAlignmentTest)a;
        System.out.println("Creating mutable alignment");
        impAlign2=new TBitAlignmentTest(a);
        System.out.println("Starting imputation");
        imputeBySiteJump(64,0,20);
    }


      public Alignment imputeBySiteJump(int minLength, int mismatchNum, int maxBestLines) {
        long time=System.currentTimeMillis();
        int knownSNPs = 0, unknownSNPs = 0, samePop=0, diffPop=0;
        int numSeqs=anchorAlignment.getSequenceCount();
        initHapLengths(0);
        System.out.println("Initial matrix created in "+(System.currentTimeMillis()-time));
        time=System.currentTimeMillis();
        int currWord=0;
        TreeMap<Integer,Integer> lenBestLine=new TreeMap<Integer,Integer>(Collections.reverseOrder());
        for (int b = 0; b < anchorAlignment.getSiteCount(); b+=64) {
            currWord=b>>6;
            if(b>0) {incrementHapLengths(b);}
            if (b % 64 == 0) {
                double rate=(double)unknownSNPs/(double)(System.currentTimeMillis()-time);
                System.out.println("Imputed base:" + b + " known:" + knownSNPs + " unknownSNPs:" + unknownSNPs+" Rate:"+rate);
                System.out.printf("SamePop:%d DiffPop:%d %n",samePop,diffPop);
//                System.out.println(b+AlignmentUtils.getSequenceString(anchorAlignment, 0));
//                System.out.println(b+AlignmentUtils.getSequenceString(impAlign2, 0));
            }
            for (int i = 0; i < numSeqs; i++) {
                lenBestLine.clear();
                for (int j = 0; j < numSeqs; j++) {
                    if(totalLen[i][j]>minLength) lenBestLine.put(totalLen[i][j], j);
                }
                if(i==0) System.out.println("CurrWord:"+currWord);
//                if(i==501) System.out.println(i+" tree:"+lenBestLine.toString());
                long mj=anchorAlignment.getTaxaBitsNoClone(i, 0).getBits()[currWord];
                long mn=anchorAlignment.getTaxaBitsNoClone(i, 1).getBits()[currWord];
//                if(i==0) System.out.printf("fmj:%64s%n",Long.toBinaryString(mj));
//                if(i==0) System.out.printf("fmn:%64s%n",Long.toBinaryString(mn));
                long miss=~(mj|mn);
                int missCnt=Long.bitCount(miss);
                unknownSNPs+=missCnt;
                knownSNPs+=(64-missCnt);
                int count=0;
                for(int bestLine: lenBestLine.values()) {
                    mj=mj|(miss&anchorAlignment.getTaxaBitsNoClone(bestLine, 0).getBits()[currWord]);
                    mn=mn|(miss&anchorAlignment.getTaxaBitsNoClone(bestLine, 1).getBits()[currWord]);
                    miss=~(mj|mn);
                    count++;
                    if(miss==0) break;
                    if(count>maxBestLines) break;
                }
//                if(i==0) System.out.printf("omj:%64s%n",Long.toBinaryString(mj));
//                if(i==0) System.out.printf("omn:%64s%n",Long.toBinaryString(mn));
//                if(i==0) System.out.printf("oms:%64s%n",Long.toBinaryString(miss));
//                if(i==0) System.out.println(b+"beforeSet"+AlignmentUtils.getSequenceString(impAlign2, 0));
                impAlign2.getTaxaBitsNoClone(i, 0).getBits()[currWord]=mj;
                impAlign2.getTaxaBitsNoClone(i, 1).getBits()[currWord]=mn;
//                if(i==0) System.out.println(b+" afterSet"+AlignmentUtils.getSequenceString(impAlign2, 0));
//                if(i==0) System.out.printf("gomj:%64s%n",Long.toBinaryString(impAlign2.getTaxaBitsNoClone(i, 0).getBits()[currWord]));
//                if(i==0) System.out.printf("gomn:%64s%n",Long.toBinaryString(impAlign2.getTaxaBitsNoClone(i, 1).getBits()[currWord]));
//                if(i==0) System.out.println("gtxt:"+AlignmentUtils.getSequenceString(impAlign2, 0).substring(b, b+64));
            }
        }
        return impAlign;
    }

    private void initHapLengths(int initialSite) {
        currDistMat=new int[2][anchorAlignment.getSequenceCount()][anchorAlignment.getSequenceCount()];
        leftHapLen=currDistMat[0];
        rightHapLen=currDistMat[1];
        totalLen=new int[anchorAlignment.getSequenceCount()][anchorAlignment.getSequenceCount()];
        for(int i=0; i<anchorAlignment.getSequenceCount(); i++) {
            for(int j=0; j<i; j++) {
               leftHapLen[j][i]=leftHapLen[i][j]=getLength(anchorAlignment, initialSite-1, i, j, maxMismatch, true)[0];
               rightHapLen[j][i]=rightHapLen[i][j]=getLength(anchorAlignment, initialSite, i, j, maxMismatch, false)[0];
               totalLen[j][i]=totalLen[i][j]=leftHapLen[j][i]+rightHapLen[j][i];
           }
        }
        currSplitSite=initialSite;
    }

    private void incrementHapLengths(int currentSite) {
        //change to words and split
        int splitWord=currSplitSite>>6;
        int currWord=currentSite>>6;
        if(currWord==splitWord) return;
//        System.out.println("Str"+Arrays.toString(totalLen[0]));
        if((splitWord+1)!=currWord) {
            System.out.println("Cannot increment.  Must recalculate.");
            initHapLengths(currentSite);
            currSplitSite=currentSite;
            return;
        }
        for(int i=0; i<anchorAlignment.getSequenceCount(); i++) {
            long[] t1mj=anchorAlignment.getTaxaBitsNoClone(i, 0).getBits(); //taxa1 major alleles
            long[] t1mn=anchorAlignment.getTaxaBitsNoClone(i, 1).getBits(); //taxa1 minor alleles
            for(int j=0; j<i; j++) {
                long[] t2mj=anchorAlignment.getTaxaBitsNoClone(j, 0).getBits();  //taxa2 major alleles  -- getting with clone is 10X slower
                long[] t2mn=anchorAlignment.getTaxaBitsNoClone(j, 1).getBits();  //taxa2 major alleles
                int diffSW=Long.bitCount(t1mj[splitWord]&t2mn[splitWord])+Long.bitCount(t2mj[splitWord]&t1mn[splitWord]);  //number of differences
                int sameSW=Long.bitCount(t1mj[splitWord]&t2mj[splitWord])+Long.bitCount(t2mn[splitWord]&t1mn[splitWord]);  //number of matches
                if(diffSW==0) {
                    leftHapLen[j][i]=leftHapLen[i][j]+=sameSW;
                    rightHapLen[j][i]=rightHapLen[i][j]-=sameSW;
                } else {
                    if(diffSW>maxMismatch) {leftHapLen[j][i]=leftHapLen[i][j]=0;}
                    else {leftHapLen[j][i]=leftHapLen[i][j]=getLength(anchorAlignment, currentSite, i, j, maxMismatch, true)[0];}
                    int diffCW=Long.bitCount(t1mj[currWord]&t2mn[currWord])+Long.bitCount(t2mj[currWord]&t1mn[currWord]);
                   // int sameCW=Long.bitCount(t1mj[currWord]&t2mj[currWord])+Long.bitCount(t2mn[currWord]&t1mn[currWord]);
                    if(diffCW>maxMismatch) {rightHapLen[j][i]=rightHapLen[i][j]=0;}
                    else {rightHapLen[j][i]=rightHapLen[i][j]=getLength(anchorAlignment, currentSite, i, j, maxMismatch, false)[0];}
                    totalLen[i][j]=totalLen[j][i]=leftHapLen[j][i]+rightHapLen[j][i];
                }
           }
        }
        currSplitSite=currentSite;
//        System.out.println("Inc"+Arrays.toString(totalLen[0]));
//        initHapLengths(currentSite);
//        System.out.println("Red"+Arrays.toString(totalLen[0]));
    }

    public boolean incrementToThisPosition(int targetPosition) {
        int leftPosition=(currSplitSite>0)?anchorAlignment.getPositionInLocus(currSplitSite-1):-1;
        if(leftPosition>targetPosition) return false;
        if(anchorAlignment.getPositionInLocus(currSplitSite)>=targetPosition) return true;
        while((currSplitSite<anchorAlignment.getSiteCount())&&(anchorAlignment.getPositionInLocus(currSplitSite)<targetPosition)) {
            incrementHapLengths(currSplitSite+1);
 //           System.out.printf("Incremented: %d %d %d %n", leftPosition, targetPosition, anchorAlignment.getPositionInLocus(currSplitSite));
        }
        return true;
    }


    /**
     * The provides a rapid bit based approach that should be very rapid.
     * @param align
     * @param initialSite
     * @param taxa1
     * @param taxa2
     * @param maxMismatch
     * @param isLeft
     * @return
     */
    static int[] getLength(TBitAlignmentTest align, int initialSite, int taxa1, int taxa2, int maxMismatch, boolean isLeft) {
        int countSame=0, countDiff=0;
        boolean ibd=true;
        long[] t1mj=align.getTaxaBitsNoClone(taxa1, 0).getBits();
        long[] t1mn=align.getTaxaBitsNoClone(taxa1, 1).getBits();
        long[] t2mj=align.getTaxaBitsNoClone(taxa2, 0).getBits();
        long[] t2mn=align.getTaxaBitsNoClone(taxa2, 1).getBits();
        int numWords=t1mj.length;
        int stop=isLeft?-1:numWords;
        int inc=isLeft?-1:1;
        int lastMatch=initialSite;
        int initLongIndex=initialSite >> 6;
        for(int b=initLongIndex; (b!=stop)&&ibd; b+=inc){
            int diff=Long.bitCount(t1mj[b]&t2mn[b])+Long.bitCount(t2mj[b]&t1mn[b]);
          //  int diff=BitUtil.pop(t1mj[b]&t2mn[b])+BitUtil.pop(t2mj[b]&t1mn[b]); //slower
            if(diff>0) {
                countDiff+=diff;
                if(countDiff>maxMismatch) {ibd=false;}
            } else {
                int same=Long.bitCount(t1mj[b]&t2mj[b])+Long.bitCount(t2mn[b]&t1mn[b]);
        //        int same=BitUtil.pop(t1mj[b]&t2mj[b])+BitUtil.pop(t2mn[b]&t1mn[b]);
                countSame+=same;
                lastMatch=b;
                if(countSame>maxMatchCalc) {ibd=false;}
            }

        }
        int[] result={countSame,lastMatch<<6};
        return result;
    }

    public void writeAlignment(String outfile) {
        ExportUtils.writeToHapmap(impAlign2, false, outfile, '\t');
    }

    public static void main(String[] args) {
	System.out.println("Running main method in FastImputation");
//        String anchorMapFile="/Users/edbuckler/SolexaAnal/GBS/hmp/HS55K_110215.hmp.txt";
//        String outfile="/Users/edbuckler/SolexaAnal/GBS/test/impAnchor110221.txt";
        for (int cNum = 10; cNum <= 10; cNum++) {
//            String anchorMapFile="/Users/edbuckler/SolexaAnal/GBS/hmp/ch"+cNum+"_NAMwLowHetF10LD_110303.hmp.txt";
//            String outfile="/Users/edbuckler/SolexaAnal/GBS/hmp/ch"+cNum+"_NAMwLowHetF10LD_110303imp.hmp.txt";
//                    String anchorMapFile="/Users/edbuckler/SolexaAnal/HapMapV2/jermpipe/SNPS201010/fusion/chr10.CSHLALLBGI.h90_f12.Q87Q87Union.hmp.txt";
        String outfile="/Users/edbuckler/SolexaAnal/GBS/test/allfusion_110504.lh.ld.c+.imp.hmp.txt";
            String anchorMapFile="/Users/edbuckler/SolexaAnal/GBS/hmp/allfusion110401/allfusion_110401.lh.ld.c+.hmp.txt";
    //        String outfile="/Users/edbuckler/SolexaAnal/GBS/test/c1_NAMimpAnchorBit110307b.txt";
        //    Alignment a=ImportUtils.readFromGZIP(anchorMapFile);
            outfile=outfile.replace("+", ""+cNum);
            anchorMapFile=anchorMapFile.replace("+", ""+cNum);
            System.out.println("Reading "+anchorMapFile);
            Alignment a=ImportUtils.readFromHapmap(anchorMapFile);
            System.out.println("p1a:"+AlignmentUtils.getSequenceString(a, 1));
            a=new TBitAlignmentTest(a);
            System.out.println("TBA:"+AlignmentUtils.getSequenceString(a, 1));
            FastImputationBit fi=new FastImputationBit(a);
            System.out.println("Writing "+outfile);
            fi.writeAlignment(outfile);
        }
    }
}
