    /*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.solexa;

import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Random;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.statistics.FisherExact;

/**
 *
 * @author edbuckler
 */
public class LDbyBitsV2 {

    /**
     * Ideas make it so that number of tests and comparison is settable
     *
     * AltVRest is likely to map elsewhere
     * RefVRest is likely to capture the indels well
     *
     * Ref V Missing may also be useful, but if we don't believe the alt then refvRest should be the best
     * A setting should allow the hets to score presense for both categories, this should increase the power of mapping the alt in other locations
     * --in the same like we should not put any filter on hets so that we can increase these frequencies.
     *
     * Output teosinte count, only run SNPs with substantial numbers in maize
     * Mask taxa with low coverage
     * 
     */
    public enum AlleleMode {
        RefvAlt, MissVPres, AltVRest, RefVRest
    }
    AlleleMode currAlleleMode=AlleleMode.MissVPres;
    static int minFreq=1;
    static int minComp=10;
    static double ldThreshold=0.000001;
    static boolean requireNCO=false;
    int numberOfAnalyses, currentAnalysis=-1;
    String[] analysisName;
    int[] bitCounts;
    byte[][] bestChr;
    byte[][] majorAF;
    byte[][] minorAF;
    int[][] bestPosition;
    double[][] bestP;
    int[][] numTests;
    int[][] numOverThreshold;
    float[][] randLETPValue;
    int longBitLength=64;
    long sampleMask;
    AlignmentInBits refAlignment;
    FisherExact fishersExact;


    public LDbyBitsV2(Alignment a, int numberOfAnalyses) {
        this.numberOfAnalyses=numberOfAnalyses;
        analysisName=new String[numberOfAnalyses];
        fishersExact=new FisherExact(1000);
        System.out.println("Converting Alignment to bits");
        long time=System.currentTimeMillis();
        AlignmentInBits aib=new AlignmentInBits(a);
        for (int i = 0; i < aib.getNumTaxa(); i++) {
            System.out.println(a.getIdGroup().getIdentifier(i).getName()+"\t"+aib.getMajorAllelesCnt(i)+"\t"+aib.getMinorAllelesCnt(i));

        }
        refAlignment=aib;
        System.out.println("Conversion done - Total time (milli sec):"+(System.currentTimeMillis()-time));
        initMatrices(aib, numberOfAnalyses);
        setSampleMask(aib.getNumTaxa());
        setSampleMask(a, aib, 0.5, false);
  //      findNeighbors(aib);
  //      writeResultsToScreen(100);
        
    }

    public void setSampleMask(int sampleSize) {
        sampleMask=0;
        for (int i = 0; i < sampleSize; i++) {
            sampleMask=sampleMask|(1L<<i);
        }
    }

    public void setSampleMask(Alignment a, AlignmentInBits aib, double minPresent, boolean removeTeosinte) {
        sampleMask=0;
        for (int i = 0; i < aib.getNumTaxa(); i++) {
            boolean keep=true;
            if(removeTeosinte&&(a.getIdGroup().getIdentifier(i).getName().startsWith("TI"))) keep=false;
            double propPresent=(double)(aib.getMajorAllelesCnt(i)+aib.getMinorAllelesCnt(i))/(double)aib.getNumSites();
            if(propPresent<minPresent) keep=false;
            if(keep) sampleMask=sampleMask|(1L<<i);
            System.out.println(a.getIdGroup().getIdentifier(i).getName()+":"+keep);
        }
    }

    public void runNeighborsInAllModes(int minDistance, int permutations) {
        for(AlleleMode am: AlleleMode.values()) {
            currAlleleMode=am;
            findNeighbors(minDistance, permutations);
        }
    }

     public void runNeighbors(int minDistance, int permutations, AlleleMode[] theModes) {
        for(AlleleMode am: theModes) {
            currAlleleMode=am;
            findNeighbors(minDistance, permutations);
        }
    }

    void initMatrices(AlignmentInBits a, int numberOfAnalyses) {
        System.out.println("Converting Alignment to bits");
        bestChr=new byte[numberOfAnalyses][a.getNumSites()];
        majorAF=new byte[numberOfAnalyses][a.getNumSites()];
        minorAF=new byte[numberOfAnalyses][a.getNumSites()];
        bestPosition=new int[numberOfAnalyses][a.getNumSites()];
        bestP=new double[numberOfAnalyses][a.getNumSites()];
        numTests=new int[numberOfAnalyses][a.getNumSites()];
        numOverThreshold=new int[numberOfAnalyses][a.getNumSites()];
        randLETPValue=new float[numberOfAnalyses][a.getNumSites()];
        for(double[] pa: bestP) Arrays.fill(pa, 1);
    }

   private int countPatterns(long[] theData) {
       Arrays.sort(theData);
       int count=1;
       for (int i = 1; i < theData.length; i++) {
           if(theData[i]!=theData[i-1]) count++;
       }
       return count;
   }

   private long getModifiedMajor(long major, long minor) {
       switch (currAlleleMode) {
           case RefvAlt: {return (major&sampleMask);}
           case MissVPres: {return (~(major|minor))&sampleMask;}
           case AltVRest: {return minor&sampleMask;}
           case RefVRest: {return major&sampleMask;}
           default: {return major&sampleMask;}
       }
           
   }
   
   private long getModifiedMinor(long major, long minor) {
       switch (currAlleleMode) {
           case RefvAlt: {return minor&sampleMask;}
           case MissVPres: {return (major|minor)&sampleMask;}
           case AltVRest: {return (~minor)&sampleMask;}
           case RefVRest: {return (~major)&sampleMask;}
           default: {return minor&sampleMask;}
       }
           
   }


    private boolean runRecordBasicTest(int site1, long site1major, long site1minor,
            int site2, long site2major, long site2minor, int recordPosition, AlignmentInBits site2Alignment) {
        double p=testSites(site1major, site1minor, site2major, site2minor, fishersExact);
        if(!Double.isNaN(p)) {
            numTests[recordPosition][site1]++;
            if(p<ldThreshold) {
                numOverThreshold[recordPosition][site1]++;
            }
            if(p<bestP[recordPosition][site1]) {
                bestP[recordPosition][site1]=p;
                bestChr[recordPosition][site1]=(byte)site2Alignment.getChr();
                bestPosition[recordPosition][site1]=site2Alignment.getPosition(site2);
            }
            return true;
        }
        return false;
    }

    public String resultHeader() {
        String[] headSuffix={"MajorCnt","MinorCnt", "BChr", "BestPos", "BestP", "NumTests", "LEThres","RandLTE_P"};
        StringBuilder sb=new StringBuilder();
        sb.append("Site Chr Position ");
        for (int j = 0; j < this.numberOfAnalyses; j++) {
           for(String hs: headSuffix) sb.append(analysisName[j]+"_"+hs+" ");
       }
       return sb.toString();
   }

   public String resultToString(int i) {
        StringBuilder sb=new StringBuilder();
        sb.append(String.format("%d %d %d ",i,refAlignment.getChr(),refAlignment.getPosition(i)));
        for (int j = 0; j < this.numberOfAnalyses; j++) {
           sb.append(String.format("%d %d %d %d %g %d %d %g ",majorAF[j][i],minorAF[j][i],bestChr[j][i],bestPosition[j][i],
                   bestP[j][i],numTests[j][i],numOverThreshold[j][i], randLETPValue[j][i]));

       }
       return sb.toString();
   }

   public void writeResultsToFile(String outfile) {
       try {
            PrintWriter fileOut = new PrintWriter(outfile);
            fileOut.write(resultHeader()); fileOut.write("\n");
            for (int i = 0; i < bestChr[0].length; i++) {
                fileOut.write(resultToString(i)); fileOut.write("\n");
            }
            fileOut.close();
        } catch (Exception e) {
            System.out.println("File output in diversityPipeline: " + e);
        }
   }

   public void writeResultsToScreen(int maxRows) {
       System.out.println(resultHeader());
       for (int i = 0; i < maxRows; i++) {
            System.out.println(resultToString(i));
        }
   }

     void findNeighbors(int minDistance, int numPermutations) {
        long time=System.currentTimeMillis();
        this.currentAnalysis++;
        this.analysisName[currentAnalysis]="Nbr"+currAlleleMode.toString();
        long hits=0, tests=0;
        int halfWindow=100;
        for (int i = 0; i < refAlignment.getNumSites(); i++) {
//            long ymj=a.getMajorAlleles(i)[0];
//            long ymn=a.getMinorAlleles(i)[0];
            long ymj=getModifiedMajor(refAlignment.getMajorAlleles(i)[0],refAlignment.getMinorAlleles(i)[0]);
            long ymn=getModifiedMinor(refAlignment.getMajorAlleles(i)[0],refAlignment.getMinorAlleles(i)[0]);
            majorAF[currentAnalysis][i]=(byte)Long.bitCount(ymj);
            minorAF[currentAnalysis][i]=(byte)Long.bitCount(ymn);
            if(minorAF[currentAnalysis][i]<=minFreq) continue;
            if(majorAF[currentAnalysis][i]<=minFreq) continue;
            int min=(i-halfWindow<0)?0:i-halfWindow;
            int max=(i+halfWindow>=refAlignment.getNumSites())?refAlignment.getNumSites()-1:i+halfWindow;
            refAlignment.setCurrentSite(min);
            for(int c=min; c<=max; c++){
                int dist=Math.abs(refAlignment.getPosition(i)-refAlignment.getPosition(c));
                if((dist<minDistance)||(c==i)) {refAlignment.nextSite(); continue;}
                tests++;
                runRecordBasicTest(i, ymj, ymn, c, refAlignment.getCurrentMajorAlleles(0),
                        refAlignment.getCurrentMinorAlleles(0),currentAnalysis,refAlignment);
                refAlignment.nextSite();
            }
        }
        long totalTime=System.currentTimeMillis()-time;
        System.out.println("Find Neighbors Total time:"+totalTime);
        System.out.printf("Tests: %d Hits: %d Tests/millis: %d \n",tests, hits, (tests/totalTime));
        if(numPermutations>0) addRandomSampleToPValue(refAlignment,currentAnalysis,numPermutations);
    }

   void addRandomSampleToPValue(AlignmentInBits a, int currentAnalysis, int maxTests) {
       long time=System.currentTimeMillis();
        long hits=0, tests=0;
        Random rg=new Random();
        int minSitesAway=refAlignment.getNumSites()/4;
        for (int i = 0; i < refAlignment.getNumSites(); i++) {
            if(i%10000==0) System.out.println("permutedPForNeighbors chr:"+refAlignment.getChr()+"+site:"+i+" with Chr:"+refAlignment.getChr()+" test:"+tests);
//            long ymj=getModifiedMajor(a.getMajorAlleles(i)[0],a.getMinorAlleles(i)[0]);
//            long ymn=getModifiedMinor(a.getMajorAlleles(i)[0],a.getMinorAlleles(i)[0]);
            long ymj=getModifiedMajor(refAlignment.getMajorAlleles(i)[0],refAlignment.getMinorAlleles(i)[0]);
            long ymn=getModifiedMinor(refAlignment.getMajorAlleles(i)[0],refAlignment.getMinorAlleles(i)[0]);
            if(bestP[currentAnalysis][i]==1) {randLETPValue[currentAnalysis][i]=1; continue;}
            int totalTests=0;
            int testLTEObsP=0;  //tests that are less than or equal to observed P
            while((totalTests < maxTests)&&(testLTEObsP<10)) {
                int r=rg.nextInt(refAlignment.getNumSites());
                if(Math.abs(i-r)<minSitesAway) continue;
                if(refAlignment.getMinorFreq(r)<=minFreq) continue;
                refAlignment.setCurrentSite(r);
                tests++;
                double p=testSites(ymj, ymn, refAlignment.getCurrentMajorAlleles(0),
                        refAlignment.getCurrentMinorAlleles(0), fishersExact);
                if(!Double.isNaN(p)) {
                    totalTests++;
                    if(p<=bestP[currentAnalysis][i]) testLTEObsP++;
                }
            }
            if(totalTests>0) {
                randLETPValue[currentAnalysis][i]=(float)testLTEObsP/(float)totalTests;
            } else {
                randLETPValue[currentAnalysis][i]=1;
            }
        }
        long totalTime=System.currentTimeMillis()-time;
        System.out.println("Find Neighbors Total time:"+totalTime);
        System.out.printf("Tests: %d Hits: %d Tests/millis: %d \n",tests, hits, (tests/totalTime));
   }

   public void createDistForLocal(int compP) {
        this.currentAnalysis++;
        this.analysisName[currentAnalysis]="RndSChr"+currAlleleMode.toString();
        long time=System.currentTimeMillis();
        long hits=0, tests=0;
        Random rg=new Random();
        int minSitesAway=refAlignment.getNumSites()/4;
        for (int i = 0; i < refAlignment.getNumSites(); i++) {
            if(i%1000==0) System.out.println("createLocalDist chr:"+refAlignment.getChr()+"+site:"+i+" with Chr:"+refAlignment.getChr()+" test:"+tests);
            long ymj=refAlignment.getMajorAlleles(i)[0];
            long ymn=refAlignment.getMinorAlleles(i)[0];
            if(refAlignment.getMinorFreq(i)<=minFreq) continue;
            int testToDo=numTests[0][i]*10;
            while(numTests[1][i] < testToDo) {
                int r=rg.nextInt(refAlignment.getNumSites());
                if(Math.abs(i-r)<minSitesAway) continue;
                if(refAlignment.getMinorFreq(r)<=minFreq) continue;
                refAlignment.setCurrentSite(r);
                tests++;
                runRecordBasicTest(i, ymj, ymn, r, refAlignment.getCurrentMajorAlleles(0), 
                        refAlignment.getCurrentMinorAlleles(0),currentAnalysis,refAlignment);
            }
        }
        long totalTime=System.currentTimeMillis()-time;
        System.out.println("Find Neighbors Total time:"+totalTime);
        System.out.printf("Tests: %d Hits: %d Tests/millis: %d \n",tests, hits, (tests/totalTime));
    }

    public void compareWithOther(Alignment targetA) {
        this.currentAnalysis++;
        this.analysisName[currentAnalysis]="OChr"+currAlleleMode.toString();
        AlignmentInBits a=new AlignmentInBits(targetA);
        targetA=null;
        long time=System.currentTimeMillis();
        long hits=0, tests=0;
        for (int i = 0; i < 10000 /*refAlignment.getNumSites()*/; i++) {
            if(i%1000==0) System.out.println("checkRest chr:"+refAlignment.chr+"+site:"+i+" with Chr:"+a.getChr()+" test:"+tests);
            if(refAlignment.getMinorFreq(i)<=minFreq) continue;
            long ymj=refAlignment.getMajorAlleles(i)[0];
            long ymn=refAlignment.getMinorAlleles(i)[0];
//            long ymj=(-ymn);
            a.setCurrentSite(0);
            for(int j=0; j<a.getNumSites(); j++) {
                if(a.getMinorFreq(j)<=minFreq) {a.nextSite();continue;}
                tests++;
                runRecordBasicTest(i, ymj, ymn, j, a.getCurrentMajorAlleles(0), 
                        a.getCurrentMinorAlleles(0),currentAnalysis,a);
                a.nextSite();
            }
        }
        long totalTime=System.currentTimeMillis()-time;
        System.out.println("Find Neighbors Total time:"+totalTime);
        System.out.printf("Tests: %d Hits: %d Tests/millis: %d \n",tests, hits, (tests/totalTime));
    }

     private static double testSites(long ymj, long ymn, long xmj, long xmn, FisherExact fishersExact) {
        double result=1;
        int c11=Long.bitCount(ymn&xmn);
        int c01=Long.bitCount(ymj&xmn);
        int c10=Long.bitCount(ymn&xmj);
        int c00=Long.bitCount(ymj&xmj);
      //   System.out.printf("%d %d %d %d %g %n",c00,c01,c10,c11, fishersExact.getTwoTailedP(c00, c01, c10, c11));
        int maj=c11+c00;
        int min=c01+c10;
        if(maj+min<minComp) return Double.NaN;
        if(requireNCO&&(c11>0)&&(c01>0)&&(c10>0)&&(c00>0)) return result;  //ensures d' of 1, evidence of no recombination
        if((maj>min)&&(maj>minComp)) {
            if((c11>=minFreq)&&(c00>=minFreq)) {
      //          double p=fishersExact.getRightTailedPQuick(c00, c01, c10, c11,.1);
                double p=fishersExact.getRightTailedP(c00, c01, c10, c11);
                return p;
            }
        } else if((maj<min)&&(min>minComp)){
            if((c10>=minFreq)&&(c01>=minFreq)) {
       //        double p=fishersExact.getRightTailedPQuick(c01, c00, c11, c10,.1);
               double p=fishersExact.getLeftTailedP(c00, c01, c10, c11);
               return p;
            }
        }
        return result;
     }

/*
     private boolean testSites(int si, int sj)  {
        if(minoraf[si]!=minoraf[sj]) return false;
        long x=sitesInBits[0][0][si]&sitesInBits[0][0][sj];
        int ld=Long.bitCount(x);
        if(ld==minoraf[si]) return true;
        return false;
     }
*/
    void initShortToBitCnt() {
        bitCounts=new int[1 << 16];
        for (int i = 0; i < 1 << 16; i++) {
            bitCounts[i]=Integer.bitCount(i);
            
        }
    }





}
