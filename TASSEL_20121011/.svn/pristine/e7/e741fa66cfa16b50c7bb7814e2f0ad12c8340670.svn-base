/*
 * BasicImputation
 */
package net.maizegenetics.gbs.pipeline;

import net.maizegenetics.gbs.util.MutableSimpleAlignment;
import net.maizegenetics.gbs.util.TBitAlignmentTest;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.datatype.DataType;

/**
 * A site by site imputation approach.  Similar to KNN, but it maintains a matrix
 * of distances between taxa and updates as it walks through the alignment.
 * 
 * @author ed
 */
public class FastImputation {
    static byte UNKNOWN=(byte)DataType.UNKNOWN_CHARACTER;
    int[][] leftHapLen, rightHapLen, totalLen;
    int currSplitSite; //it is included on the right side
    int[][][] currDistMat;

    static int maxMismatch=0;
    static int maxMatchCalc=300;
    //Alignment anchorAlignment=null;
    TBitAlignmentTest anchorAlignment=null;
    MutableSimpleAlignment impAlign=null;

    public FastImputation(Alignment a) {
        this.anchorAlignment=(TBitAlignmentTest)a;
        System.out.println("Creating mutable alignment");
        impAlign=new MutableSimpleAlignment(a);
        System.out.println("Starting imputation");
        imputeBySite(20,0);
    }

    /**
     * This imputation approach looks for the longest matching haplotype surrounding
     * every missing SNP. It finds the greatest length by only counting sites where the two haplotypes
     * have non-missing data, and then will find the longest matching stretch in either direction with
     * the number of mismatches.
     *
     * The a taxa with a missing data has no match of the minLength, then the major allele is used to replace the SNP
     *
     * Currently this requires a Pack1Alignment as there not a duplication method for simple alignment.
     * @param anchorAlignment Alignment used
     * @param minLength threshold below which majority base is used.
     * @param mismatchNum the number of misMatching bases in the length search
     * @return Imputed Pack1Alignment
     */
    public Alignment imputeBySite(int minLength, int mismatchNum) {
        long time=System.currentTimeMillis();
        int knownSNPs = 0, unknownSNPs = 0;
        int numSeqs=anchorAlignment.getSequenceCount();
        initHapLengths(0);
        System.out.println("Initial matrix created in "+(System.currentTimeMillis()-time));
        time=System.currentTimeMillis();
        int currWord=0;
        for (int b = 0; b < anchorAlignment.getSiteCount(); b++) {
            if(b>0) {
                
                incrementHapLengths(b);
            }

            char[] alleles = null;
            boolean firstTime = true;
            if (b % 10 == 0) {
                double rate=(double)unknownSNPs/(double)(System.currentTimeMillis()-time);
                System.out.println("Imputed base:" + b + " known:" + knownSNPs + " unknownSNPs:" + unknownSNPs+" Rate:"+rate);
            }
            for (int i = 0; i < numSeqs; i++) {
                byte focusBase = anchorAlignment.getBase(i, b);
                if (focusBase == DataType.UNKNOWN_BYTE) {
                    int maxMatch = minLength, bestLine = -1;
                    for (int j = 0; j < numSeqs; j++) {
                        if(anchorAlignment.getBase(j, b)==DataType.UNKNOWN_BYTE) continue;
                        int lengthij = leftHapLen[i][j] + rightHapLen[i][j];
                  //      int lengthij = totalLen[i][j];
                        if (lengthij > maxMatch) {
                            maxMatch = lengthij;
                            bestLine = j;
                        }
                    }
//                    System.out.println("maxMatch:"+maxMatch);
                    unknownSNPs++;
                    if (bestLine != -1) {
                        impAlign.setBase(i, b, anchorAlignment.getBase(bestLine, b));
 //                       System.out.println(unknownSNPs+"\t"+bestLine+"\t"+maxMatch+"\t"+anchorAlignment.getBase(bestLine, b));
                    } else {
                        impAlign.setBase(i, b, anchorAlignment.getMajorAllele(b));
//                        if (firstTime) {
//                            alleles = anchorAlignment.getSiteSummary(b).getAlleles();
//                            firstTime = false;
//                        }
//                        if ((alleles == null) || (alleles.length == 0)) {
//                            impAlign.setBase(i, b, DataType.UNKNOWN_BYTE);
//                        } else {
//                            impAlign.setBase(i, b, (byte)alleles[0]);  //use majority base
//                        }
                    }
                } else {
                    knownSNPs++;
                }
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
        if((currSplitSite+1)!=currentSite) {
            System.out.println("Error in incremeting sites in incrementHapLengths");
            return;
        }
      //  if((anchorAlignment instanceof TBitAlignmentTest)&&(currentSite%64!=0)) {currSplitSite=currentSite; return;}
        for(int i=0; i<anchorAlignment.getSequenceCount(); i++) {
            for(int j=0; j<i; j++) {
                byte s1b=anchorAlignment.getBase(i,currentSite-1);
                byte s2b=anchorAlignment.getBase(j,currentSite-1);
                if(s1b==UNKNOWN) continue;
                if(s2b==UNKNOWN) continue;
                if(s1b==s2b) {
                    leftHapLen[j][i]++; leftHapLen[i][j]++;
                    rightHapLen[j][i]--; rightHapLen[i][j]--;
                } else {
                    leftHapLen[j][i]=leftHapLen[i][j]=0;
                    rightHapLen[j][i]=rightHapLen[i][j]=getLength(anchorAlignment, currentSite, i, j, maxMismatch, false)[0];
                    totalLen[j][i]=totalLen[i][j]=leftHapLen[j][i]+rightHapLen[j][i];
                }
           }
        }
        currSplitSite=currentSite;
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

    static int[] getLength(Alignment align, int initialSite, int taxa1, int taxa2, int maxMismatch, boolean isLeft) {
        int countSame=0, countDiff=0;
        boolean ibd=true;
        int stop=isLeft?-1:align.getSiteCount();
        int inc=isLeft?-1:1;
        int lastMatch=initialSite;
        for(int b=initialSite; (b!=stop)&&ibd; b+=inc){
            byte s1b=align.getBase(taxa1,b);
            byte s2b=align.getBase(taxa2,b);
            if(s1b==UNKNOWN) continue;
            if(s2b==UNKNOWN) continue;
            //homozygous filter
//            if((s1b!=align.getMajorAllele(b))&&(s1b!=align.getMinorAllele(b))) continue;
//            if((s2b!=align.getMajorAllele(b))&&(s2b!=align.getMinorAllele(b))) continue;
            if(s1b==s2b) {
                countSame++;
                lastMatch=b;
                if(countSame>maxMatchCalc) {ibd=false;}
            } else {
                countDiff++;
                if(countDiff>maxMismatch) {ibd=false;}
            }
        }
        int[] result={countSame,lastMatch};
        return result;
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
        ExportUtils.writeToHapmap(impAlign, false, outfile, '\t');
    }

    public static void main(String[] args) {
	System.out.println("Running main method in FastImputation");
//        String anchorMapFile="/Users/edbuckler/SolexaAnal/GBS/hmp/HS55K_110215.hmp.txt";
//        String outfile="/Users/edbuckler/SolexaAnal/GBS/test/impAnchor110221.txt";
//        String anchorMapFile="/Users/edbuckler/SolexaAnal/GBS/hmp/ch1_NAMwLDLowHetF10_110303.hmp.BLOB.gz";
        String anchorMapFile="/Users/edbuckler/SolexaAnal/GBS/hmp/h2000_ch1_NAMwLDLowHetF10_110303.hmp.txt";
        String outfile="/Users/edbuckler/SolexaAnal/GBS/test/c1_NAMimpAnchor110307ok.txt";
    //    Alignment a=ImportUtils.readFromGZIP(anchorMapFile);
        Alignment a=ImportUtils.readFromHapmap(anchorMapFile);
        a=new TBitAlignmentTest(a);
        FastImputation fi=new FastImputation(a);
        fi.writeAlignment(outfile);
    }
}
