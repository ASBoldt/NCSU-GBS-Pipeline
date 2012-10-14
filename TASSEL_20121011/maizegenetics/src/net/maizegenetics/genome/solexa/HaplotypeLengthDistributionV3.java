/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.solexa;

//import cern.colt.Arrays;
import net.maizegenetics.genome.*;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import net.maizegenetics.pal.datatype.DataType;
import net.maizegenetics.pal.alignment.*;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;

/**
 *
 * @author ed
 */
public class HaplotypeLengthDistributionV3 {
    static byte UNKNOWN=(byte)DataType.UNKNOWN_CHARACTER;
    int[][] leftHapLen, rightHapLen;
    int currSplitSite; //it is included on the right side
    Random randGen;
    int[][][] currDistMat;
    int[][][] minorAlleleFreqDist=null;
 //   int[][] taxaWrong=new int[2][105];

    static int maxMismatch=0;
    Alignment anchorAlignment=null;

    //do these same stats for the lines with the minor alleles
    //remove the hets


    public HaplotypeLengthDistributionV3(Alignment theAlign) {
        this.anchorAlignment=theAlign;
        this.randGen=new Random();
        //this distribution is used for the random permutations, so the same population structure is used
        this.minorAlleleFreqDist=indexOfMinorAlleleFreq(theAlign);      
        initHapLengths(0);
        leftHapLen=currDistMat[0];
        rightHapLen=currDistMat[1];
    }

    public HaplotypeLengthDistributionV3(Alignment theAlign, String outfile) {
        randGen=new Random();
        this.anchorAlignment=theAlign;
        //this distribution is used for the random permutations, so the same population structure is used
        minorAlleleFreqDist=indexOfMinorAlleleFreq(theAlign);
        //create a simple missing site, so that the same methods can be used for everything
        byte[] missingSite=new byte[theAlign.getSequenceCount()];
        Arrays.fill(missingSite,UNKNOWN);
        StringBuffer sb=new StringBuffer();
        FileWriter fw=null;

        int permutations=0;
        try{
         if(outfile!=null) {fw=new FileWriter(outfile);}
         sb.append("Position\tSite\tMAF\tObsHapLength\tObsWOCurrSite\tRandomGreaterObs\tAvgRandom\tStDevRandom\tZobsvRand\t"+
                 "ImpCorr\tImpTotal\tImpLongCorr\tImpLongTotal\tImpMinorCorr\tImpMinorTotal\n");
         if(fw!=null) {fw.write(sb.toString());}
            else System.out.print(sb.toString());
         for (int i = 0; i < theAlign.getSiteCount(); i++) {
            if(i%10000==0) {
                System.out.println(outfile+" site:"+i);
           //     for(int x=0; x<105; x++) System.out.println(theAlign.getTaxaName(x)+" "+taxaWrong[0][x]+" "+taxaWrong[1][x]);

            }
            sb=new StringBuffer();
            if(i==0) {initHapLengths(0);
                leftHapLen=currDistMat[0];
                rightHapLen=currDistMat[1];
            }
            else {incrementHapLengths(i);}
            double obsLength=this.getLengthWithReplacementSite(i, extractSite(theAlign, i))[0];
            sb.append(theAlign.getPositionInLocus(i)+"\t"+i+"\t"+theAlign.getMinorAlleleFrequency(i)+"\t"+obsLength);
            double obsSiteToMissingLength=this.getLengthWithReplacementSite(i, missingSite)[0];
            sb.append("\t"+obsSiteToMissingLength);
            int countBeatObs=0;
            double sumRandom=0, sumRandSqr=0;
            if(getRandomSiteWithPairedFrequency(theAlign,i)!=null) {  //tests whether the null distribution is too small to test accurately
                for(int n=0; n<permutations; n++) {
              //      long rand=this.getLengthWithReplacementSite(theAlign, i, currDistMat, extractSite(theAlign, i));
              //      long rand=this.getLengthWithReplacementSite(theAlign, i, currDistMat, missingSite);
              //      double rand=this.getLengthWithReplacementSite(theAlign, i, currDistMat, permuteSite(extractSite(theAlign, i)))[0];
                    double rand=this.getLengthWithReplacementSite(i,getRandomSiteWithPairedFrequency(theAlign,i))[0];

                    if(rand>obsLength) countBeatObs++;
                    sumRandom+=rand;
                    sumRandSqr+=(rand*rand);
                }
                double avgRandom=sumRandom/(double)permutations;
                double stdev=Math.sqrt((sumRandSqr/(double)permutations)-(avgRandom*avgRandom));
                double z=((double)obsLength-avgRandom)/stdev;
                sb.append("\t"+countBeatObs+"\t"+avgRandom+"\t"+stdev+"\t"+z+"\t");
            } else {
                sb.append("\tNaN\tNaN\tNaN\tNaN\t");
            }
            int[] imputeResults=estimateImputationAccuracy(extractSite(theAlign, i),true);
//            int[] imputeResults=estimateImputationAccuracy(theAlign, i, currDistMat, getRandomSiteWithPairedFrequency(theAlign,i));
            sb.append(imputeResults[0]+"\t"+imputeResults[1]+"\t"+imputeResults[2]+"\t"+imputeResults[3]+"\t"+imputeResults[4]+"\t"+imputeResults[5]+"\n");
            if(fw!=null) {fw.write(sb.toString());}
                else System.out.print(sb.toString());
            }
         if(fw!=null) {fw.flush(); fw.close();}
        } catch (IOException e) {
            System.out.println("Error in HaplotypeLengthDistributionV2(Alignment theAlign, String outfile):"+e);

        }
    }

    public HaplotypeLengthDistributionV3(Alignment theAlign, int minLengthToOutput, String outfile) {
        randGen=new Random();
        this.anchorAlignment=theAlign;
        StringBuffer sb=new StringBuffer();
        FileWriter fw=null;
        try{
         if(outfile!=null) {fw=new FileWriter(outfile);}
         sb.append("Chr\tStartSite\tEndSite\tStartPosition\tEndPosition\tTaxa1\tTaxa2\tPerfectMatchLength\n");
         if(fw!=null) {fw.write(sb.toString());}
            else System.out.print(sb.toString());
         sb=new StringBuffer();
         for (int i = 0; i < theAlign.getSequenceCount(); i++) {
             for (int j = 0; j < i; j++) {
                 int currSite=0;
                 while(currSite<theAlign.getSiteCount()-1) {
                     int[] result=this.getLength(theAlign, currSite, i, j, maxMismatch, false);
                     if(result[0]>minLengthToOutput) {
                         sb.append(String.format("%s %d %d %d %d %s %s %d %n",theAlign.getLocus(currSite),currSite, result[1], anchorAlignment.getPositionInLocus(currSite),anchorAlignment.getPositionInLocus(result[1]),
                                 anchorAlignment.getTaxaName(i), anchorAlignment.getTaxaName(j), result[0]));
                     }
                     currSite=result[1]+1;
                 }

             }
            if(fw!=null) {fw.write(sb.toString());}
                else System.out.print(sb.toString());
            sb=new StringBuffer();
            }
         if(fw!=null) {fw.flush(); fw.close();}
        } catch (IOException e) {
            System.out.println("Error in HaplotypeLengthDistributionV2(Alignment theAlign, String outfile):"+e);

        }
    }

  
    /**
     *
     * @param testSite byte array of the tested site
     * @param testAtCurrentSite if true then equal replacing the current position, if false like inserting before the current position
     * @return array of accuracy and tests{correctAll, testableAll, correctInLongHaps, testableInLongHaps, correctMinor, testableMinor}
     */
    public int[] estimateImputationAccuracy(byte[] testSite, boolean testAtCurrentSite) {
        //result correctAll, testableAll, correctInLongHaps, testableInLongHaps, correctMinor, testableMinor,correctMinorAvgLength, wrongMinorAvgLength
        int[] result=new int[8];
        int minLength=50;
        byte minorAllele=getMinorAllele(testSite);
        for (int i = 0; i < anchorAlignment.getSequenceCount(); i++) {
            int maxMatch=-1, bestLine=-1;
            if(testSite[i]==UNKNOWN) continue;
            result[1]++;
            for(int j=0; j<anchorAlignment.getSequenceCount(); j++) {
                if((i==j)||(testSite[j]==UNKNOWN)) continue;
                int lengthij=-1;
                if((testAtCurrentSite==false)||(anchorAlignment.getBase(i, currSplitSite)==UNKNOWN)||(anchorAlignment.getBase(j, currSplitSite)==UNKNOWN)) {
                    lengthij = leftHapLen[i][j] + rightHapLen[i][j];
                } else if(anchorAlignment.getBase(i, currSplitSite)==anchorAlignment.getBase(j, currSplitSite)) {lengthij = leftHapLen[i][j] + rightHapLen[i][j]-1;}
                else {lengthij = leftHapLen[i][j] + getLength(anchorAlignment, currSplitSite+1, i, j, maxMismatch, false)[0];}
                if(lengthij>maxMatch) {maxMatch=lengthij; bestLine=j;}
            }
            if(bestLine<0)   {System.out.println(Arrays.toString(testSite)); continue;}
            if(maxMatch>minLength) {
                result[3]++;
                if(testSite[i]==testSite[bestLine]) result[2]++;
            }
            if(testSite[i]==minorAllele) {
                result[5]++;
                if(testSite[i]==testSite[bestLine]) {
                    result[4]++;
                    result[6]+=maxMatch;
                } else {
                    result[7]+=maxMatch;
                }
            }
            
            if(testSite[i]==testSite[bestLine]) result[0]++;
//            if((maxMatch>50)&&(!anchorAlignment.getTaxaName(i).startsWith("W"))) {
//                System.out.printf("%d %s %s %d %n",this.currSplitSite, anchorAlignment.getTaxaName(i), anchorAlignment.getTaxaName(bestLine), maxMatch);
//            }
         
        }
        if(result[4]>0) {result[6]=result[6]/result[4];} else {result[6]=-1;}
        int wrong=result[5]-result[4];
        if(wrong>0) {result[7]=result[7]/wrong;} else {result[7]=-1;}
        return result;
    }



    private byte getMinorAllele(byte[] sites) {
        int[] counts=new int[Byte.MAX_VALUE];
        for(byte b: sites) counts[b]++;
        byte majorAllele=-1, minorAllele=-1;
        int majorCnt=0, minorCnt=0;
        for (byte i = 0; i < counts.length; i++) {
            if(i==UNKNOWN) continue;
            if(counts[i]>minorCnt) {
                if(counts[i]>majorCnt) {
                    minorAllele=majorAllele; minorCnt=majorCnt;
                    majorAllele=i; majorCnt=counts[i];
                } else {
                    minorAllele=i; minorCnt=counts[i];
                }
            }
        }
        return minorAllele;
    }

    //TODO Find SNPs that break long haplotypes
    //TODO Impute missing data based on longest shared hit
    //TODO Make haplotype groups

    //the problem with this is that the missing data rate is not the same
    //the missing rate will place a large variance on the lengths
    //I could match by frequency and presense
    public static int[][][] indexOfMinorAlleleFreq(Alignment align) {
        int[][] counts=new int[align.getSequenceCount()+1][align.getSequenceCount()/2+1];
        for (int s = 0; s < align.getSiteCount(); s++) {
            SiteSummary ss=align.getSiteSummary(s);
            if(ss.getAlleleCounts().length>1) {
                int majorFreq=ss.getAlleleCounts()[0]/2;
                int minorFreq=ss.getAlleleCounts()[1]/2;
                counts[majorFreq+minorFreq][minorFreq]++;
            } else counts[0][0]++;
        }
        int[][][] alleleFreqSite=new int[align.getSequenceCount()+1][align.getSequenceCount()/2+1][];
        for (int i = 0; i < alleleFreqSite.length; i++) {
            for (int j = 0; j < alleleFreqSite[0].length; j++) {
                alleleFreqSite[i][j]=new int[counts[i][j]];
            }
        }
        counts=new int[align.getSequenceCount()+1][align.getSequenceCount()/2+1];
        for (int s = 0; s < align.getSiteCount(); s++) {
            SiteSummary ss=align.getSiteSummary(s);
            if(ss.getAlleleCounts().length>1) {
                int majorFreq=ss.getAlleleCounts()[0]/2;
                int minorFreq=ss.getAlleleCounts()[1]/2;
                int index=counts[majorFreq+minorFreq][minorFreq];
                alleleFreqSite[majorFreq+minorFreq][minorFreq][index]=s;
                counts[majorFreq+minorFreq][minorFreq]++;
            } else counts[0][0]++;
        }
        return alleleFreqSite;
    }


    public void initHapLengths(int initialSite) {
        currDistMat=new int[2][anchorAlignment.getSequenceCount()][anchorAlignment.getSequenceCount()];
        leftHapLen=currDistMat[0];
        rightHapLen=currDistMat[1];
        for(int i=0; i<anchorAlignment.getSequenceCount(); i++) {
            for(int j=0; j<i; j++) {
               leftHapLen[j][i]=leftHapLen[i][j]=getLength(anchorAlignment, initialSite-1, i, j, maxMismatch, true)[0];
               rightHapLen[j][i]=rightHapLen[i][j]=getLength(anchorAlignment, initialSite, i, j, maxMismatch, false)[0];

           }
        }
        currSplitSite=initialSite;
    }

    public void incrementHapLengths(int currentSite) {
        if((currSplitSite+1)!=currentSite) {
            System.out.println("Error in incremeting sites in incrementHapLengths");
            return;
        }
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

    public byte[] extractSite(Alignment align, int site) {
        byte[] result=new byte[align.getSequenceCount()];
        for (int t = 0; t < result.length; t++) {
            result[t]=align.getBase(t, site);
        }
        return result;
    }

    public byte[] permuteSite(byte[] inSite) {
        byte[] randMarkers=Arrays.copyOf(inSite, inSite.length);
        for(int i=0; i<inSite.length; i++) {
                int r=randGen.nextInt(inSite.length);
                byte temp=randMarkers[i];
                randMarkers[i]=randMarkers[r];
                randMarkers[r]=temp;
            }
        return randMarkers;
    }

    public byte[] getRandomSiteWithPairedFrequency(Alignment align, int site) {
        SiteSummary ss=align.getSiteSummary(site);
        if(ss.getAlleleCounts().length<2) return null;
        int majorCnt=ss.getAlleleCounts()[0]/2;
        int minorCnt=ss.getAlleleCounts()[1]/2;
        int index=majorCnt+minorCnt;
        if(minorAlleleFreqDist[index][minorCnt].length<10) {
            System.out.printf("Site: %d MajorCnt: %d MinorCnt: %d  InsufficentCnts%n",site, majorCnt, minorCnt);
            return null;
        }  //insufficient data to produce a permutation
        int r=randGen.nextInt(minorAlleleFreqDist[index][minorCnt].length);
//        int m2=align.getSiteSummary(minorAlleleFreqDist[index][minorCnt][r]).getAlleleCounts()[1]/2;
//        System.out.println(minorCnt+" "+r+" "+minorAlleleFreqDist[index][minorCnt][r]+" "+m2);
//        System.out.println(ss.getNumberMissing()+" "+align.getSiteSummary(minorAlleleFreqDist[index][minorCnt][r]).getNumberMissing());
        return extractSite(align, minorAlleleFreqDist[index][minorCnt][r]);
    }



     public double[] getLengthWithReplacementSite(int currentSite, byte[] testSite) {
        long totalLength=0;
        double sumSqrLength=0;
        int length;
        int tests=0;
        int[] maxLength=new int[anchorAlignment.getSequenceCount()];
        for(int i=0; i<anchorAlignment.getSequenceCount(); i++) {
            byte s1b=anchorAlignment.getBase(i,currentSite);
            byte r1b=testSite[i];
            for(int j=0; j<i; j++) {
                length=0;
                byte s2b=anchorAlignment.getBase(j,currentSite);
                byte r2b=testSite[j];
                int newRightLen=((s1b!=UNKNOWN)&&(s1b==s2b))?rightHapLen[i][j]-1:rightHapLen[i][j];
                if((r1b==UNKNOWN)||(r2b==UNKNOWN)) {
                    length=newRightLen+leftHapLen[i][j];
                    if(newRightLen==0) {length+=getLength(anchorAlignment, currentSite+1, i, j, maxMismatch, false)[0];}
                }
                else if(r1b==r2b) {
                    if(rightHapLen[i][j]==0) {
                        length=getLength(anchorAlignment, currentSite+1, i, j, maxMismatch, false)[0]+leftHapLen[i][j]+1;
                    } else {
                        length=newRightLen+leftHapLen[i][j]+1;
                    }
                }
                else {//unequal
                    length=leftHapLen[i][j];
                }
 //               System.out.printf("Rep: %d %d %d %d %d %d %d %d %n",currentSite, i,j,s1b,s2b,length,r1b,r2b);
                totalLength+=length;
                tests++;
                if(length>maxLength[i]) maxLength[i]=length;
                if(length>maxLength[j]) maxLength[j]=length;
                sumSqrLength+=(double)length*(double)length;
           }
        }
        double[] result=new double[4];
        result[0]=totalLength;
        result[1]=totalLength/tests;
        result[2]=Math.sqrt((sumSqrLength/(double)tests)-(result[1]*result[1]));
        int maxSum=0;
        for(int m: maxLength) maxSum+=m;
        result[3]=(double)maxSum/(double)maxLength.length;
        return result;
    }

    public double[] getLengthWithInsertionOfSite(int currentSite, byte[] testSite) {
        long totalLength=0;
        double sumSqrLength=0;
        int length;
        int tests=0;
        int[] maxLength=new int[anchorAlignment.getSequenceCount()];
        for(int i=0; i<anchorAlignment.getSequenceCount(); i++) {
            byte r1b=testSite[i];
            for(int j=0; j<i; j++) {
                length=0;
                byte r2b=testSite[j];
                if((r1b==UNKNOWN)||(r2b==UNKNOWN)) {
                    length=rightHapLen[i][j]+leftHapLen[i][j];
                }
                else if(r1b==r2b) {
                    length=rightHapLen[i][j]+leftHapLen[i][j]+1;
                }
                else {//unequal
                    length=leftHapLen[i][j];
                }
//                System.out.printf("Rep: %d %d %d %d %d %d %d %d %n",currentSite, i,j,s1b,s2b,length,r1b,r2b);
                totalLength+=length;
                tests++;
                if(length>maxLength[i]) maxLength[i]=length;
                if(length>maxLength[j]) maxLength[j]=length;
                sumSqrLength+=(double)length*(double)length;
           }
        }
        double[] result=new double[4];
        result[0]=totalLength;
        result[1]=totalLength/tests;
        result[2]=Math.sqrt((sumSqrLength/(double)tests)-(result[1]*result[1]));
        int maxSum=0;
        for(int m: maxLength) maxSum+=m;
        result[3]=(double)maxSum/(double)maxLength.length;
        return result;
    }


    public int getLengthRemoveCurrSite(Alignment anchorAlignment, int currentSite) {
 //       int maxMismatch=0;
        int totalLength=0;
        int length;
        for(int i=0; i<anchorAlignment.getSequenceCount(); i++) {
            byte s1b=anchorAlignment.getBase(i,currentSite);
            for(int j=0; j<i; j++) {
                length=0;
                byte s2b=anchorAlignment.getBase(j,currentSite);
                if((s1b==UNKNOWN)||(s2b==UNKNOWN)) {
                    length=rightHapLen[i][j]+leftHapLen[i][j];
                }
                else if(s1b==s2b) {
                    length=rightHapLen[i][j]+leftHapLen[i][j]-1;
                }
                else {
                    length=getLength(anchorAlignment, currentSite+1, i, j, maxMismatch, false)[0]+leftHapLen[i][j];
                }
               // length=getLength(align, currentSite, i, j, maxMismatch, false)+leftHapLen[i][j];
  //              System.out.printf("Rem: %d %d %d %d %d %d %n",currentSite, i,j,s1b,s2b,length);
                totalLength+=length;
           }
        }
        return totalLength;
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
            } else {
                countDiff++;
                if(countDiff>maxMismatch) {ibd=false;}
            }
        }
        int[] result={countSame,lastMatch};
        return result;
    }

 public static void main(String[] args) throws Exception {
      //  String infile="/Users/edbuckler/SolexaAnal/SNP55K/IBM_55K_hapmap.txt";
  //      String infile="/Users/edbuckler/SolexaAnal/SNP55K/chr10From55K.hmp.txt";
        String infile="/Users/edbuckler/SolexaAnal/SNP55K/SNP55K_hapmapV2Samples_AGPv1_20100823.txt";
     //   String infile="/Users/edbuckler/SolexaAnal/HapMapV2/hp1/h100kchr10.genotypes.log2ct_h90.hmp.txt";
      //  String infile="/Users/edbuckler/SolexaAnal/HapMapV2/anchor/chr10.genotypes.log2ct_h90_f50.good.hmp.txt";
//        String infile="/Users/edbuckler/SolexaAnal/HapMapV2/anchor/chr10.genotypes.log2ct_h90_f50.good.good.hmp.txt";
//        String outfile="/Users/edbuckler/SolexaAnal/HapMapV2/test/chr10.genotypes.log2ct_h90_f50.good.good.hmp.ehh.txt";
        String outfile="/Users/edbuckler/SolexaAnal/HapMapV2/test/crap.txt";

        /**
         * Change the hets to missing as all sort of this class are having problems with this.
         *
         * Also need to examine the LD of the missing data
         */
/*
        Pack1Alignment p1a = (Pack1Alignment) ImportUtils.createPack1AlignmentFromFile(infile, "", "");


        System.out.println("Alignment Loaded:"+infile);
        System.out.printf("Sites %d Taxa %d %n", p1a.getSiteCount(),p1a.getSequenceCount());
  //      characterizeHomozygous(p1a,50000);
        p1a=(Pack1Alignment)HaplotypeLengthDistributionV2.makeHomozygousAlignment(p1a);
 //       characterizeHomozygous(p1a,50000);
        HaplotypeLengthDistributionV3 hld3=new HaplotypeLengthDistributionV3(p1a, 50, null);
   */
        Pack1Alignment p1a;
        for (int i = 1; i <= 10; i++) {
           p1a=(Pack1Alignment)ImportUtils.readFromHapmap(infile, ""+i);
           p1a=(Pack1Alignment)HaplotypeLengthDistributionV2.makeHomozygousAlignment(p1a);
           HaplotypeLengthDistributionV3 hld3=new HaplotypeLengthDistributionV3(p1a, 150, null);

     }
    }

}
