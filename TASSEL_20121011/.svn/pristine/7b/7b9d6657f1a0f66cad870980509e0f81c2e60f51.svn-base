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
public class LDbyBits {
    static int minFreq=1;
    static int minComp=10;
    static boolean requireNCO=false;
    int[] position;
    long[][][] sitesInBits;
    byte[] minoraf;
    byte[] majoraf;
    int[] bitCounts;
    byte[][] bestChr;
    int[][] bestPosition;
    double[][] bestP;
    int[][] numTests;
    int longBitLength=64;

    public LDbyBits(Alignment a) {
        initMatrices(a);
//        long[] siteCopy=sitesInBits[0].clone();
//        System.out.printf("Sites %d Patterns %d \n",siteCopy.length, countPatterns(siteCopy));
//       learnProb();

        findNeighbors(a);
        createDistForLocal(a);
//        checkRest(a);


    }

    public LDbyBits(Alignment a, Alignment a2) {
        AlignmentInBits aib=new AlignmentInBits(a);
        AlignmentInBits aib2=new AlignmentInBits(a2);
        System.exit(0);
        initMatrices(a);
//        long[] siteCopy=sitesInBits[0].clone();
//        System.out.printf("Sites %d Patterns %d \n",siteCopy.length, countPatterns(siteCopy));
//       learnProb();

        findNeighbors(a);
        createDistForLocal(a);
//        checkRest(a);


    }

     public LDbyBits(Alignment a, boolean space) {
        AlignmentInBits aib=new AlignmentInBits(a);
        System.exit(0);
        initMatrices(a);
        findNeighbors(a);
        createDistForLocal(a);
//        checkRest(a);


    }



    void initMatrices(Alignment a) {
        System.out.println("Converting Alignment to bits");
        long time=System.currentTimeMillis();
        int arraySize=(a.getSequenceCount()/longBitLength)+1;
        position=new int[a.getSiteCount()];
        sitesInBits=new long[2][arraySize][a.getSiteCount()];
        minoraf=new byte[a.getSiteCount()];
        majoraf=new byte[a.getSiteCount()];
        bestChr=new byte[2][a.getSiteCount()];
        bestPosition=new int[2][a.getSiteCount()];
        bestP=new double[2][a.getSiteCount()];
        numTests=new int[2][a.getSiteCount()];
        Arrays.fill(bestP[0], 1); Arrays.fill(bestP[1], 1);
        int minorCnt=0, majorCnt=0;
        for (int i = 0; i < a.getSiteCount(); i++) {
            position[i]=a.getPositionInLocus(i);
            byte minorAllele=(byte)a.getMinorAllele(i);
            byte majorAllele=(byte)a.getMajorAllele(i);
            for(int j=0; j<a.getSequenceCount(); j++) {
                if(minorAllele==a.getBase(j, i)) {
                    minoraf[i]++;
                    sitesInBits[1][0][i]=setBit(sitesInBits[1][0][i],j,true);
                    minorCnt++;
                }
                if(majorAllele==a.getBase(j, i)) {
                    majoraf[i]++;
                    sitesInBits[0][0][i]=setBit(sitesInBits[0][0][i],j,true);
                    majorCnt++;
                }
            }
            if(minoraf[i]>30) {
                System.out.print(a.getSNPID(i)+" ");
                for(int j=0; j<a.getSequenceCount(); j++) {System.out.print(a.getBaseChar(j, i));}
                System.out.println("Minor Allele:"+(char)minorAllele);
            }
        }
        System.out.printf("Sites: %d Taxa: %d MajorAlleles: %d MinorAlleles: %d \n",
                    sitesInBits[0][0].length, a.getSequenceCount(), majorCnt, minorCnt);
        System.out.println("Conversion done - Total time (milli sec):"+(System.currentTimeMillis()-time));
    }

   private int countPatterns(long[] theData) {
       Arrays.sort(theData);
       int count=1;
       for (int i = 1; i < theData.length; i++) {
           if(theData[i]!=theData[i-1]) count++;
       }
       return count;
   }


     void findNeighbors(Alignment a) {
        long time=System.currentTimeMillis();
        long hits=0, tests=0;
        int sij;
        FisherExact fishersExact=new FisherExact(1000);
        byte thisChr=Byte.parseByte(a.getLocusName(0));
        for (int i = 0; i < sitesInBits[0][0].length; i++) {
            long ymj=sitesInBits[0][0][i];
            long ymn=sitesInBits[1][0][i];
            if(Long.bitCount(ymn)<=minFreq) continue;
            if(Long.bitCount(ymj)<=minFreq) continue;
         //   if(minoraf[i]<4) {tests+=(i-1); continue;}
            for(int j = 1; j < 100; j++) {
                sij=i-j;
                if(sij>-1) {
                    tests++;
                    float p=(float)testSites(ymj, ymn, sitesInBits[0][0][sij], sitesInBits[1][0][sij], fishersExact);
                    if(!Double.isNaN(p)) {
                        numTests[0][i]++;
                        if(p<0.0001) hits++;
                        if(p<bestP[0][i]) {
                            bestP[0][i]=p;
                            bestChr[0][i]=thisChr;
                            bestPosition[0][i]=a.getPositionInLocus(sij);
                        }
                    }
                }
                sij=i+j;
                if(sij<sitesInBits[0][0].length) {
                    tests++;
                    float p=(float)testSites(ymj, ymn, sitesInBits[0][0][sij], sitesInBits[1][0][sij], fishersExact);
                    if(!Double.isNaN(p)) {
                        numTests[0][i]++;
                        if(p<0.0001) hits++;
                        if(p<bestP[0][i]) {
                            bestP[0][i]=p;
                            bestChr[0][i]=thisChr;
                            bestPosition[0][i]=a.getPositionInLocus(sij);
                        }
                    }
                }
            }
        }
        long totalTime=System.currentTimeMillis()-time;
        System.out.println("Find Neighbors Total time:"+totalTime);
        System.out.printf("Tests: %d Hits: %d Tests/millis: %d \n",tests, hits, (tests/totalTime));
         for (int i = 0; i < 1000; i++) {
             if(minoraf[i]>minFreq)
             System.out.printf("%d %d %d %d %d %g %d %n",i,a.getPositionInLocus(i),minoraf[i],
                     bestChr[0][i],bestPosition[0][i],bestP[0][i], numTests[0][i]);

         }
    }

   public void writeResultsToFile(String outfile) {
       try {
            PrintWriter fileOut = new PrintWriter(outfile);
            fileOut.write("Site\tPosition\tMajorCnt\tMinorCnt\tClose_Chr\tClose_Pos\tClose_P\tClose_N\t"+
                    "Far_Chr\tFar_Pos\tFar_P\tFar_N\n");
            for (int i = 0; i < bestChr[0].length; i++) {
                fileOut.printf("%d %d %d %d %d %d %g %d %d %d %g %d %n",i,position[i],majoraf[i], minoraf[i],
                     bestChr[0][i],bestPosition[0][i],bestP[0][i],numTests[0][i],
                     bestChr[1][i],bestPosition[1][i],bestP[1][i], numTests[1][i]);
            }
            fileOut.close();
        } catch (Exception e) {
            System.out.println("File output in diversityPipeline: " + e);
        }
   }

   void createDistForLocal(Alignment a) {
        long time=System.currentTimeMillis();
        long hits=0, tests=0;
        FisherExact fishersExact=new FisherExact(1000);
        byte thisChr=Byte.parseByte(a.getLocusName(0));
        Random rg=new Random();
        for (int i = 0; i < sitesInBits[0][0].length; i++) {
            if(i%10000==0) System.out.println("checkRest i:"+i+" test:"+tests);
            long ymj=sitesInBits[0][0][i];
            long ymn=sitesInBits[1][0][i];
            if(Long.bitCount(ymn)<=minFreq) continue;
            if(Long.bitCount(ymj)<=minFreq) continue;
            while(numTests[1][i] < (numTests[0][i]*10)) {
                int r=rg.nextInt(a.getSiteCount());
                if(Math.abs(i-r)<10000) continue;
                if(minoraf[r]<=minFreq) continue;
                if(majoraf[r]<=minFreq) continue;
                tests++;
                float p=(float)testSites(ymj, ymn, sitesInBits[0][0][r], sitesInBits[1][0][r], fishersExact);
                if(!Double.isNaN(p)) {
                    numTests[1][i]++;
                    if(p<0.0001) hits++;
                    if(p<bestP[1][i]) {
                        bestP[1][i]=p;
                        bestChr[1][i]=thisChr;
                        bestPosition[1][i]=a.getPositionInLocus(r);
                    }
                }
            }
        }
        long totalTime=System.currentTimeMillis()-time;
        System.out.println("Find Neighbors Total time:"+totalTime);
        System.out.printf("Tests: %d Hits: %d Tests/millis: %d \n",tests, hits, (tests/totalTime));
         for (int i = 0; i < 10000; i++) {
             if(minoraf[i]>minFreq)
             System.out.printf("%d %d %d %d %d %d %g %d %d %d %g %d %n",i,a.getPositionInLocus(i),majoraf[i], minoraf[i],
                     bestChr[0][i],bestPosition[0][i],bestP[0][i],numTests[0][i],
                     bestChr[1][i],bestPosition[1][i],bestP[1][i], numTests[1][i]);

         }
    }

      void checkRest(Alignment a) {
        long time=System.currentTimeMillis();
        long hits=0, tests=0;
        FisherExact fishersExact=new FisherExact(1000);
        byte thisChr=Byte.parseByte(a.getLocusName(0));
        for (int i = 0; i < 1000000; i++) {
            if(i%1000==0) System.out.println("checkRest i:"+i+" test:"+tests);
            long ymj=sitesInBits[0][0][i];
            long ymn=sitesInBits[1][0][i];
            if(Long.bitCount(ymn)<=minFreq) continue;
            if(Long.bitCount(ymj)<=minFreq) continue;
            if(bestP[0][i]<1e-4) continue;
            for(int j = 1; j < sitesInBits[0][0].length; j++) {
                if(Math.abs(i-j)<1000) continue;
                if(minoraf[j]<=minFreq) continue;
                if(majoraf[j]<=minFreq) continue;
                tests++;
                float p=(float)testSites(ymj, ymn, sitesInBits[0][0][j], sitesInBits[1][0][j], fishersExact);
                if(p<0.0001) hits++;
                if(p<bestP[1][i]) {
                    bestP[1][i]=p;
                    bestChr[1][i]=thisChr;
                    bestPosition[1][i]=a.getPositionInLocus(j);
                }
            } 
        }
        long totalTime=System.currentTimeMillis()-time;
        System.out.println("Find Neighbors Total time:"+totalTime);
        System.out.printf("Tests: %d Hits: %d Tests/millis: %d \n",tests, hits, (tests/totalTime));
         for (int i = 0; i < 1000000; i++) {
             if(minoraf[i]>minFreq)
             System.out.printf("%d %d %d %d %d %d %g %d %d %g %n",i,a.getPositionInLocus(i),majoraf[i], minoraf[i],
                     bestChr[0][i],bestPosition[0][i],bestP[0][i],
                     bestChr[1][i],bestPosition[1][i],bestP[1][i]);

         }
    }

     private static double testSites(long ymj, long ymn, long xmj, long xmn, FisherExact fishersExact) {
        double result=1;
        int c11=Long.bitCount(ymn&xmn);
        int c01=Long.bitCount(ymj&xmn);
        int c10=Long.bitCount(ymn&xmj);
        int c00=Long.bitCount(ymj&xmj);
        int maj=c11+c00;
        int min=c01+c10;
        if(maj+min<minComp) return Double.NaN;
        if(requireNCO&&(c11>0)&&(c01>0)&&(c01>0)&&(c00>0)) return result;  //ensures d' of 1, evidence of no recombination
        if((maj>min)&&(maj>minComp)) {
            if((c11>=minFreq)&&(c00>=minFreq)) {
                double p=fishersExact.getRightTailedP(c00, c01, c10, c11);
                return p;
            }
        } else if((maj<min)&&(min>minComp)){
            if((c10>=minFreq)&&(c01>=minFreq)) {
                double p=fishersExact.getLeftTailedP(c00, c01, c10, c11);
                return p;
            }
        }
        return result;
     }

     private boolean testSites(int si, int sj)  {
        if(minoraf[si]!=minoraf[sj]) return false;
        long x=sitesInBits[0][0][si]&sitesInBits[0][0][sj];
        int ld=Long.bitCount(x);
        if(ld==minoraf[si]) return true;
        return false;
     }

    void initShortToBitCnt() {
        bitCounts=new int[1 << 16];
        for (int i = 0; i < 1 << 16; i++) {
            bitCounts[i]=Integer.bitCount(i);
            
        }
    }


    boolean getBit(long set, int bit) { // get bit 'bit' from the set
       return (set&(1L<<bit)) != 0;
    }

    long setBit(long set, int bit, boolean value) { // return the new set with bit 'bit' (un)set
       if (value)
          return set|(1L<<bit);
       else
          return set&~(1<<bit);
    }


}
