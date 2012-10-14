/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.solexa;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.Pack1Alignment;
import org.apache.commons.math.stat.inference.TTestImpl;

/**
 *
 * @author edbuckler
 */
public class LDofCNVs {
    String[] cnvTaxa;
    int[] cnvChr, cnvPos;
    double[][] cnvValue;
    TTestImpl tTest=new TTestImpl();
    int binSize=10000;
    int windowSize=50000;
    boolean ignoreLDInBin=false;
    boolean pullRandomSNPS=false;
    int minLines=50;
    int minMinorAlleles=5;


    public LDofCNVs(String infileCNV, String infileSNPs, String outfile) {
        //double[] x={1,2,3,4,5};
        double[] x={-1.0379, -0.2175, -1.8558, -0.1836, -0.5916, -1.1844, -0.7986, 0.5826, -0.0557, 0.1894, -0.0915, 0.2583, 0.0105, 0.5154, -0.1337, -0.7674, -0.3794, -0.1147, -0.0213, -0.0496, 0.1796, 0.2951, 0.321, 0.1444, -0.3559, 0.2475, 0.0634, -4.3826, 0.195, -0.4935, 0.1542};
    //    double[] y={3,2,4,-1};
        double[] y={0.3731, -2.8396, -3.7358, -2.558, 0.1468, -2.8904, -3.7463, -2.9185, -2.6007, -2.1279, 0.7156, -2.3687, 0.2986, 0.2257, -1.5875, -1.5107, -0.1276, -1.0568, -1.8605, -1.6963, -0.0739, -3.258};
        try{
        System.out.println("p="+tTest.homoscedasticTTest(x,y));}
        catch(Exception e) {}
        
        readCNVFile(infileCNV);
        System.out.println("Alignment Loading:"+infileSNPs);
        Pack1Alignment p1a = (Pack1Alignment) ImportUtils.createPack1AlignmentFromFile(infileSNPs, "", "");
        System.out.println("Alignment Loaded:"+infileSNPs);
        System.out.printf("Sites %d Taxa %d %n", p1a.getSiteCount(),p1a.getSequenceCount());
        int[] p1aRedirect = new int[p1a.getSequenceCount()];
        for (int t = 0; t < cnvTaxa.length; t++) {
            p1aRedirect[t] = p1a.getIdGroup().whichIdNumber(cnvTaxa[t]);
//            System.out.println(cnvTaxa[t]+":"+p1aRedirect[t]);
        }
        System.out.println("Chr cnvPos bestSite bestPosition numTests bestMinorCnt minP");
        for (int i = 0; i < cnvPos.length; i++) {
            int startPos=cnvPos[i]-windowSize;
            int endPos=cnvPos[i]+binSize+windowSize;
            int startSite=p1a.getSiteOfPhysicalPosition(startPos,null);
            if(startSite<0) startSite=-(startSite+1);
            int endSite=p1a.getSiteOfPhysicalPosition(endPos,null);
            if(endSite<0) endSite=-(endSite+1);
            if(endSite>p1a.getSiteCount()-1) endSite=p1a.getSiteCount()-1;
 //           System.out.printf("%d %d %d %d %d %n", cnvPos[i], startPos, endPos, startSite, endSite);
            double minP=1;
            int bestSite=-1, bestMinorCnt=-1;
            int tests=0;
            for (int rs = startSite; rs < endSite; rs++) {
                int s=rs;
                if(pullRandomSNPS) {
                    while(Math.abs(s-rs)<500000) {
                        s=(int)Math.floor(Math.random()*(p1a.getSiteCount()-1));
                    }
                }
                int currPos=p1a.getPositionInLocus(s);
                if(ignoreLDInBin&&(currPos>=cnvPos[i])&&(currPos<=(cnvPos[i]+binSize))) continue;
                int majorLines=0, minorLines=0;
                double[] majorValues=new double[p1a.getSequenceCount()];
                double[] minorValues=new double[p1a.getSequenceCount()];
                for (int t = 0; t < cnvTaxa.length; t++) {
                    if(p1aRedirect[t]<0) continue;  //not matching taxa
                    byte allele=p1a.getBase(p1aRedirect[t], s);
                    if(allele==p1a.getMajorAllele(s)) {
                        majorValues[majorLines]=cnvValue[i][t];
                        majorLines++;
                    } else if(allele==p1a.getMinorAllele(s)) {
                        minorValues[minorLines]=cnvValue[i][t];
                        minorLines++;
                    }
                }
                if((majorLines+minorLines)<minLines) continue;
                if(minorLines<minMinorAlleles) continue;
                double[] majorValues2=new double[majorLines];
                for (int j = 0; j < majorLines; j++) {
                    majorValues2[j] = majorValues[j];
                }
                double[] minorValues2=new double[minorLines];
                for (int j = 0; j < minorLines; j++) {
                    minorValues2[j] = minorValues[j];
                }
//                double[] majorValues2=Arrays.copyOf(majorValues, majorLines);
//                System.out.println(s+" mj:"+Arrays.toString(majorValues2));
//                double[] minorValues2=Arrays.copyOf(minorValues, minorLines);
//                System.out.println(s+" mn:"+Arrays.toString(minorValues2));
                double p=Double.NaN;

                try{
                    p=tTest.homoscedasticTTest(majorValues2, minorValues2);
                    tests++;
                  //  System.out.println("rp="+tTest.homoscedasticTTest(majorValues2, majorValues2));
              //      System.out.println(s+" p:"+p);
                  //  System.out.println("p="+tTest.homoscedasticTTest(x,y));
                }catch(Exception e) {

                }
                if((!Double.isNaN(p))&&(p<minP)) {
                    bestSite=s;
                    minP=p;
                    bestMinorCnt=minorLines;
                }
            }
            System.out.printf("%d %d %d %d %d %d %g %n",cnvChr[i], cnvPos[i], bestSite,
                    p1a.getPositionInLocus(bestSite), tests, bestMinorCnt, minP);
        }
    }

    void readCNVFile(String infile) {
        try{
            BufferedReader fileIn = new BufferedReader(new FileReader(infile), 100000);
            String h=fileIn.readLine().replace("#", "");
            cnvTaxa=h.split("\\s");
            int rows=SimpleTextFile.countRows(infile)-1;
            cnvChr=new int[rows];
            cnvPos=new int[rows];
            cnvValue=new double[rows][cnvTaxa.length];
            for (int i = 0; i < rows; i++) {
                String[] s=fileIn.readLine().split("\\s");
                cnvChr[i]=Integer.parseInt(s[0]);
                cnvPos[i]=Integer.parseInt(s[1]);
                for(int j=0; j<cnvTaxa.length; j++) {
                    cnvValue[i][j]=Double.parseDouble(s[2+j]);
                }
            }
            fileIn.close();
        } catch (IOException e) {
            System.out.println("countRows error at row");
            e.printStackTrace();
        }
    }

   public static void main(String[] args) throws Exception {
          String inDir1="/Users/edbuckler/SolexaAnal/HapMapV2/jermpipe/SNPS201010/10kb_cov_matrices/chr.#.10kb_log2_matrix.txt";
          String inDir2="/Users/edbuckler/SolexaAnal/HapMapV2/jermpipe/SNPS201010/fusion/chr#.CSHLALLBGI.h90_f12.Q87Q87Union.hmp.txt";
          String outDir1="/Users/edbuckler/SolexaAnal/HapMapV2/jermpipe/SNPS201010/null";
          for (int i = 1; i <= 10; i++) {
            String in1=inDir1.replaceFirst("#", ""+i);
            String in2=inDir2.replaceFirst("#", ""+i);
            String out1=outDir1.replaceFirst("#", ""+i);
            LDofCNVs hld3=new LDofCNVs(in1, in2, out1);

        }

    }

}
