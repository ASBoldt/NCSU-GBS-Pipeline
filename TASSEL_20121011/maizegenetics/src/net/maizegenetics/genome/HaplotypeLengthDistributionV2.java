/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome;

//import cern.colt.Arrays;
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
public class HaplotypeLengthDistributionV2 {
    static byte UNKNOWN=(byte)DataType.UNKNOWN_CHARACTER;
    //int[][] leftHapLen, rightHapLen;
    int currSplitSite; //it is included on the right side
    Random randGen;
    int[][][] currDistMat;
    int[][][] minorAlleleFreqDist=null;
 //   int[][] taxaWrong=new int[2][105];

    static int maxMismatch=0;
    Alignment anchorAlignment=null;

    //do these same stats for the lines with the minor alleles
    //remove the hets


    public HaplotypeLengthDistributionV2(Alignment theAlign) {
        randGen=new Random();
        //this distribution is used for the random permutations, so the same population structure is used
        minorAlleleFreqDist=indexOfMinorAlleleFreq(theAlign);
        anchorAlignment=theAlign;
        currDistMat=initHapLengths(theAlign,0);
    }

    public HaplotypeLengthDistributionV2(Alignment theAlign, String outfile) {
        randGen=new Random();
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
            if(i==0) {currDistMat=initHapLengths(theAlign,0);}
            else {incrementHapLengths(theAlign, i, currDistMat);}
            double obsLength=this.getLengthWithReplacementSite(theAlign, i, currDistMat, extractSite(theAlign, i))[0];
            sb.append(theAlign.getPositionInLocus(i)+"\t"+i+"\t"+theAlign.getMinorAlleleFrequency(i)+"\t"+obsLength);
            double obsSiteToMissingLength=this.getLengthWithReplacementSite(theAlign, i, currDistMat, missingSite)[0];
            sb.append("\t"+obsSiteToMissingLength);
            int countBeatObs=0;
            double sumRandom=0, sumRandSqr=0;
            if(getRandomSiteWithPairedFrequency(theAlign,i)!=null) {  //tests whether the null distribution is too small to test accurately
                for(int n=0; n<permutations; n++) {
              //      long rand=this.getLengthWithReplacementSite(theAlign, i, currDistMat, extractSite(theAlign, i));
              //      long rand=this.getLengthWithReplacementSite(theAlign, i, currDistMat, missingSite);
              //      double rand=this.getLengthWithReplacementSite(theAlign, i, currDistMat, permuteSite(extractSite(theAlign, i)))[0];
                    double rand=this.getLengthWithReplacementSite(theAlign, i, currDistMat,getRandomSiteWithPairedFrequency(theAlign,i))[0];

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
            int[] imputeResults=estimateImputationAccuracy(theAlign, i, currDistMat, extractSite(theAlign, i));
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

    public HaplotypeLengthDistributionV2(Alignment theAlign, int minCount, double minFreq) {
         //System.out.println("Position\tSite\tMAF\tObsHapLength\tObsWOCurrSite\tRandomGreaterObs\tAvgRandom\tStDevRandom\tZobsvRand");
        System.out.printf("Sites %d Taxa %d %n", theAlign.getSiteCount(),theAlign.getSequenceCount());
       // IdGroup maizeIDGroup=getMaizeIDS(theAlign,true,false,false,false);
        IdGroup maizeIDGroup=getMaizeIDS(theAlign,true,true,true,false);
        theAlign=FilterAlignment.getInstance(theAlign,maizeIDGroup);
        System.out.printf("Sites %d Taxa %d %n", theAlign.getSiteCount(),theAlign.getSequenceCount());
        for(int mc=10; mc<theAlign.getSequenceCount(); mc+=10) {
            for(double mf=0.00; mf<0.21; mf+=0.05) {
             Alignment filterAlign=AnnotatedAlignmentUtils.removeSitesBasedOnFreqIgnoreGapsMissing(theAlign,mf,mc);
             int totalSites=filterAlign.getSiteCount();
             double sumLength=0;
             double sumStdDev=0;
             double sumAvgMax=0;
             for (int i = 0; i < totalSites; i++) {
                if(i==0) {currDistMat=initHapLengths(theAlign,0);}
                else {incrementHapLengths(theAlign, i, currDistMat);}
                double[] result=this.getLengthWithReplacementSite(theAlign, i, currDistMat, extractSite(theAlign, i));
                double obsLength=result[0];
                sumLength+=result[0];
                sumStdDev+=result[2];
                sumAvgMax+=result[3];
            }
            double avgLength=sumLength/(double)totalSites;
            double avgLpTaxa=(2*avgLength)/(theAlign.getSequenceCount()*(theAlign.getSequenceCount()-1));
            double avgStdDev=sumStdDev/(double)totalSites;
            double avgAvgMax=sumAvgMax/(double)totalSites;
            System.out.printf("%d %g %d %g %g %g %g %n", mc, mf, totalSites, avgLength, avgLpTaxa, avgStdDev, avgAvgMax);
            }
        }
    }

    private int[] estimateImputationAccuracy(Alignment align, int currentSite, int[][][] hapLen, byte[] testSite) {
        //result correctAll, testableAll, correctInLongHaps, testableInLongHaps, correctMinor, testableMinor,
        //int maxMismatch=0;
        int[] result=new int[6];
        int minLength=100;
        int[][] leftHapLen=hapLen[0];
        int[][] rightHapLen=hapLen[1];
        byte minorAllele=getMinorAllele(testSite);
        for (int i = 0; i < align.getSequenceCount(); i++) {
            int maxMatch=-1, bestLine=-1;
            if(testSite[i]==UNKNOWN) continue;
//            if(align.getTaxaName(i).startsWith("TI")) continue;
//            if(align.getTaxaName(i).startsWith("TDD")) continue;
            result[1]++;
            for(int j=0; j<align.getSequenceCount(); j++) {
                if((i==j)||(testSite[j]==UNKNOWN)) continue;
                int lengthij=-1;
                if(align.getBase(i, currentSite)==align.getBase(j, currentSite)) {lengthij = leftHapLen[i][j] + rightHapLen[i][j];}
                else {lengthij = leftHapLen[i][j] + getLength(align, currentSite+1, i, j, maxMismatch, false);}
//                if(rightHapLen[i][j]==0) {lengthij = leftHapLen[i][j] + getLength(align, currentSite+1, i, j, maxMismatch, false);}
//                    else {lengthij = leftHapLen[i][j] + rightHapLen[i][j];}
                if(lengthij>maxMatch) {maxMatch=lengthij; bestLine=j;}
            }
//            taxaWrong[0][i]++;
//            if(testSite[i]!=testSite[bestLine]) taxaWrong[1][i]++;

            if(maxMatch>minLength) {
                result[3]++;
                if(testSite[i]==testSite[bestLine]) result[2]++;
            }
            if(testSite[i]==minorAllele) {
                result[5]++;
                if(testSite[i]==testSite[bestLine]) result[4]++;
            }
            if(testSite[i]==testSite[bestLine]) result[0]++;
         //   System.out.printf("%s %s %d %d %d %n", (char)testSite[i], (char)testSite[bestLine],maxMatch, result[0],result[1]);
        }
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


    public int[][][] initHapLengths(Alignment align, int initialSite) {
 //       int maxMismatch=0;
        int[][][] hapLen=new int[2][align.getSequenceCount()][align.getSequenceCount()];
        int[][] leftHapLen=hapLen[0];
        int[][] rightHapLen=hapLen[1];
        for(int i=0; i<align.getSequenceCount(); i++) {
            for(int j=0; j<i; j++) {
               leftHapLen[j][i]=leftHapLen[i][j]=getLength(align, initialSite-1, i, j, maxMismatch, true);
               rightHapLen[j][i]=rightHapLen[i][j]=getLength(align, initialSite, i, j, maxMismatch, false);

           }
        }
        currSplitSite=initialSite;
        return hapLen;
    }

    public int[][][] incrementHapLengths(Alignment align, int currentSite, int[][][] hapLen) {
 //       int maxMismatch=0;
        if((currSplitSite+1)!=currentSite) {
            System.out.println("Error in incremeting sites in incrementHapLengths");
            return null;
        }
        int[][] leftHapLen=hapLen[0];
        int[][] rightHapLen=hapLen[1];
        for(int i=0; i<align.getSequenceCount(); i++) {
            for(int j=0; j<i; j++) {
                byte s1b=align.getBase(i,currentSite-1);
                byte s2b=align.getBase(j,currentSite-1);
                if(s1b==UNKNOWN) continue;
                if(s2b==UNKNOWN) continue;
                if(s1b==s2b) {
                    leftHapLen[j][i]++; leftHapLen[i][j]++;
                    rightHapLen[j][i]--; rightHapLen[i][j]--;
                } else {
                    leftHapLen[j][i]=leftHapLen[i][j]=0;
                    rightHapLen[j][i]=rightHapLen[i][j]=getLength(align, currentSite, i, j, maxMismatch, false);
                }
           }
        }
        currSplitSite=currentSite;
        return hapLen;
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



     public double[] getLengthWithReplacementSite(Alignment align, int currentSite, int[][][] hapLen, byte[] testSite) {
     //   int maxMismatch=0;
        long totalLength=0;
        double sumSqrLength=0;
        int length;
        int tests=0;
        int[][] leftHapLen=hapLen[0];
        int[][] rightHapLen=hapLen[1];
        int[] maxLength=new int[align.getSequenceCount()];
        for(int i=0; i<align.getSequenceCount(); i++) {
            byte s1b=align.getBase(i,currentSite);
            byte r1b=testSite[i];
            for(int j=0; j<i; j++) {
                length=0;
                byte s2b=align.getBase(j,currentSite);
                byte r2b=testSite[j];
                int newRightLen=((s1b!=UNKNOWN)&&(s1b==s2b))?rightHapLen[i][j]-1:rightHapLen[i][j];
                if((r1b==UNKNOWN)||(r2b==UNKNOWN)) {
                    length=newRightLen+leftHapLen[i][j];
                    if(newRightLen==0) {length+=getLength(align, currentSite+1, i, j, maxMismatch, false);}
                }
                else if(r1b==r2b) {
                    if(rightHapLen[i][j]==0) {
                        length=getLength(align, currentSite+1, i, j, maxMismatch, false)+leftHapLen[i][j]+1;
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

    public double[] getLengthWithInsertionOfSite(Alignment align, int currentSite, int[][][] hapLen, byte[] testSite) {
        //insertion is before the currentSite, i.e. between (currentSite-1) and currentSite
 //       int maxMismatch=0;
        long totalLength=0;
        double sumSqrLength=0;
        int length;
        int tests=0;
        int[][] leftHapLen=hapLen[0];
        int[][] rightHapLen=hapLen[1];
        int[] maxLength=new int[align.getSequenceCount()];
        for(int i=0; i<align.getSequenceCount(); i++) {
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


    public static int getLengthRemoveCurrSite(Alignment align, int currentSite, int[][][] hapLen) {
 //       int maxMismatch=0;
        int totalLength=0;
        int length;
        int[][] leftHapLen=hapLen[0];
        int[][] rightHapLen=hapLen[1];
        for(int i=0; i<align.getSequenceCount(); i++) {
            byte s1b=align.getBase(i,currentSite);
            for(int j=0; j<i; j++) {
                length=0;
                byte s2b=align.getBase(j,currentSite);
                if((s1b==UNKNOWN)||(s2b==UNKNOWN)) {
                    length=rightHapLen[i][j]+leftHapLen[i][j];
                }
                else if(s1b==s2b) {
                    length=rightHapLen[i][j]+leftHapLen[i][j]-1;
                }
                else {
                    length=getLength(align, currentSite+1, i, j, maxMismatch, false)+leftHapLen[i][j];
                }
               // length=getLength(align, currentSite, i, j, maxMismatch, false)+leftHapLen[i][j];
  //              System.out.printf("Rem: %d %d %d %d %d %d %n",currentSite, i,j,s1b,s2b,length);
                totalLength+=length;
           }
        }
        return totalLength;
    }

    static int sumArrays(int[][] a) {
        int sum=0;
        for(int[] i: a) {
            for(int j: i) sum+=j;
        }
        return sum;
    }

    static int sumArrays(int[][][] a) {
        int sum=0;
        for(int[][] i: a) {
            sum+=sumArrays(i);
        }
        return sum;
    }

    static int sumCountLongRuns(int[][][] a, int minLongRun) {
        int count=0;
        for (int i = 0; i < a[0].length; i++) {
            for(int j=0; j<i; j++) {
                if((a[0][i][j]+a[1][i][j])>minLongRun) count++;
            }
        }
        return count;
    }

    static void matToString(int[][] arr, int rows) {
        for (int i = 0; i < rows; i++) {
            for (int j=0; j<rows; j++) System.out.print(arr[i][j]+"\t");
            System.out.println("");

        }
    }

    static int getLength(Alignment align, int initialSite, int taxa1, int taxa2, int maxMismatch, boolean isLeft) {
        int countSame=0, countDiff=0;
        boolean ibd=true;
        int stop=isLeft?-1:align.getSiteCount();
        int inc=isLeft?-1:1;
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
            } else {
                countDiff++;
                if(countDiff>maxMismatch) {ibd=false;}
            }
        }
        return countSame;
    }

    static int[][] characterizeHomozygous(Alignment theAlign, int bin) {
        int[][] counts=new int[2][theAlign.getSequenceCount()];
        for (int i = 0; i < theAlign.getSiteCount(); i++) {
            for(int t=0; t< theAlign.getSequenceCount(); t++) {
                byte b=theAlign.getBase(t, i);
                if(b==UNKNOWN) continue;
                if(AllelePositionBLOBUtils.isBaseHomozygous(b)) {counts[0][t]++;}
                else {counts[1][t]++;}
            }
            if((i%bin==0)&&(i>1)) {
               for(int t=0; t< theAlign.getSequenceCount(); t++) {
                System.out.printf("%d %s %d %d %g %n",i,theAlign.getIdGroup().getIdentifier(t),counts[0][t],
                        counts[1][t],(double)counts[1][t]/(double)(counts[0][t]+counts[1][t]));
                }
               counts=new int[2][theAlign.getSequenceCount()];
            }
        }
        return counts;
    }

    public static Alignment makeHomozygousAlignment(Pack1Alignment inAlign) {
        for(int t=0; t< inAlign.getSequenceCount(); t++) {
                byte[] alleleBLOB=inAlign.getAlleleBLOBs(t);
                for(int s=0; s<inAlign.getSiteCount(); s++) {
                    if(!AllelePositionBLOBUtils.isBaseHomozygous(inAlign.getBase(t, s))) {
                        AllelePositionBLOBUtils.setHalfByteInAlleleBLOB(alleleBLOB, s, DataType.UNKNOWN_CHARACTER);
                    }
                }
        }
        return inAlign;
    }

    public IdGroup getMaizeIDS(Alignment a, boolean inbred, boolean landrace, boolean teo, boolean trips) {
        ArrayList<Identifier> subIDAL=new ArrayList<Identifier>();
        int set;
        for (int i = 0; i < a.getSequenceCount(); i++) {
            boolean keep=true;
            Identifier cID=a.getIdGroup().getIdentifier(i);
            if(cID.getName().startsWith("TI")) {if(teo==false) continue;}
            else if(cID.getName().startsWith("TDD")) {if(trips==false) continue;}
            else if(cID.getName().startsWith("BK")) {if (landrace==false) continue;}
            else if(inbred==false) {continue;}
            subIDAL.add(cID);
            System.out.println(cID.getName()+":"+keep);
        }
        Identifier[] subID=new Identifier[subIDAL.size()];
        IdGroup subIdGroup=new SimpleIdGroup(subIDAL.toArray(subID));
        return subIdGroup;
    }

    /**
     *
     * @param args
     * @throws java.lang.Exception
     */
    public static void main(String[] args) throws Exception {
      //  String infile="/Users/edbuckler/SolexaAnal/HapMapV2/anchor/h100kchr10.genotypes.log2ct_h90.good.hmp.txt";
     //   String infile="/Users/edbuckler/SolexaAnal/HapMapV2/hp1/h100kchr10.genotypes.log2ct_h90.hmp.txt";
      //  String infile="/Users/edbuckler/SolexaAnal/HapMapV2/anchor/chr10.genotypes.log2ct_h90_f50.good.hmp.txt";
        String infile="/Users/edbuckler/SolexaAnal/HapMapV2/anchor/chr10.genotypes.log2ct_h90_f50.good.good.hmp.txt";
        String outfile="/Users/edbuckler/SolexaAnal/HapMapV2/test/chr10.genotypes.log2ct_h90_f50.good.good.hmp.ehh.txt";
        //String outfile=null;
      //  String infile="/Users/edbuckler/SolexaAnal/HapMapV2/hp1/chr10.genotypes.log2ct_h90.hmp.txt";
  //      String infile="/Users/edbuckler/SolexaAnal/HapMapV2/hp1/done/chr8.genotypes.log2ct_h90.hmp.txt";
 
        /**
         * Change the hets to missing as all sort of this class are having problems with this.
         *
         * Also need to examine the LD of the missing data
         */

        Pack1Alignment p1a = (Pack1Alignment) ImportUtils.createPack1AlignmentFromFile(infile, "", "");


        System.out.println("Alignment Loaded:"+infile);
        System.out.printf("Sites %d Taxa %d %n", p1a.getSiteCount(),p1a.getSequenceCount());
        characterizeHomozygous(p1a,50000);
 //       makeHomozygousAlignment(p1a);
 //       characterizeHomozygous(p1a,50000);
 //       HaplotypeLengthDistributionV2 hld2=new HaplotypeLengthDistributionV2(AnnotatedAlignmentUtils.removeSitesBasedOnFreqIgnoreGapsMissing(p1a,0.00,50),outfile);
        HaplotypeLengthDistributionV2 hld2=new HaplotypeLengthDistributionV2(p1a, outfile);
       // HaplotypeLengthDistributionV2 hld2=new HaplotypeLengthDistributionV2(p1a,1,0.01);

       // System.out.println("Site\tBaseArraySum\tNoCurrArraySum\tBaseLongRuns\tNoCurrLongRuns");
     /*   for (i = 0; i < 10000; i++) {
  //          hapLen=initHapLengths(p1a,i);
//            System.out.println("Left:"+i);
//            matToString(hapLen[0],10);
//            System.out.println("Right:"+i);
//            matToString(hapLen[1],10);
//            System.out.println("Old Haplo:"+i);
//            matToString(HaplotypeLengthDistribution.maxHaplotypeLengthMatrix(p1a, i, 0),10);
            if(i==0) {hapLenInc=initHapLengths(p1a,0);}
            else {hapLenInc=incrementHapLengths(p1a, i, hapLenInc);}
        //    System.out.println("Sum:"+sumArrays(hapLenInc));
            hapLenRemCurr=removeCurrSiteHapLengths(p1a,i, initHapLengths(p1a,i));
            System.out.print(i+ "\t"+sumArrays(initHapLengths(p1a,i))+"\t"+sumArrays(hapLenRemCurr));
            System.out.print("\t"+getLengthRemoveCurrSite(p1a,i,initHapLengths(p1a,i)));
            System.out.println("\t"+sumCountLongRuns(hapLenInc,40)+"\t"+sumCountLongRuns(hapLenRemCurr,40));
//            System.out.println("Inc Left:"+i);
//            matToString(hapLenInc[0],10);
//            System.out.println("Inc Right:"+i);
//            matToString(hapLenInc[1],10);
//            if(i>997) {
//                System.out.println("Inc Left:"+i);
//                matToString(hapLenInc[0],105);
//                System.out.println("Inc Right:"+i);
//                matToString(hapLenInc[1],105);
//            }
        }
 */
      
 //       System.out.println(Arrays.deepToString(initHapLengths(p1a,i)));
//        GdpdmBLOBUtils.writePack1AlignmentToZip(p1aV2, infile2);
//        lengthDistribution(p1a);
      //  p1a=AllelePositionBLOBUtils.permuteAlignment(p1a);
      //  lengthDistPerSite(p1a);
//        maxHaplotypeLengthMatrix(p1a,0,1);

    }
}
