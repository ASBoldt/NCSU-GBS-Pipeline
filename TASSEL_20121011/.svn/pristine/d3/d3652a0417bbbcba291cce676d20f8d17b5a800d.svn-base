/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.gbs.pipeline;

import java.awt.Frame;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import javax.swing.ImageIcon;
import net.maizegenetics.gbs.util.ArgsEngine;
import net.maizegenetics.gbs.util.OpenBitSet;
import net.maizegenetics.gbs.util.SBitAlignment;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.distance.DistanceMatrix;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.pal.statistics.FisherExact;
import net.maizegenetics.pal.tree.NeighborJoiningTree;
import net.maizegenetics.pal.tree.Node;
import net.maizegenetics.pal.tree.NodeUtils;
import net.maizegenetics.pal.tree.TreeUtils;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import org.apache.log4j.Logger;

/**
 * Finds high LD markers, clusters them, and identifies potential linkage groups.
 * This is evolving to just a method for using the markers that are unique to a single
 * chromosome in a population.  
 *
 * @author edbuckler
 */
public class PhaseF1TestPlugin extends AbstractPlugin{
    private static final Logger myLogger=Logger.getLogger(PhaseF1TestPlugin.class);
    static File inputDirectory=null;
    private String[] infiles=null;
    static String outfile = null;
    static ArgsEngine myArgsEngine = null;
    private SBitAlignment myDataSBit;
//    private byte[][][] myParentalChrComp;  //parent, chr, site: assigned to major or minor allele (e.g. A or C); unknown Byte.Min
    private static int numParents=2;
    private static int numChr=2; //diploid
    private static double  minPValue=0.001;  //0.001 is normal
    private byte[][][] offspring; //taxon, parent, site:  assigned to chr of parent (either 0, 1); unknown Byte.Min
    private int[] highLD;  //count of sites in high LD with this site
    private int minHighLDSites=50;  //50 is normal
    private int numParentGametes=4;
    private byte[] minorAlleleCluster;  //normally 0,1,2,3 but if there are more than 4 clusters there could be 5 or 6
    private byte[] gameteOriginOfCluster; //0=mom-g1, 1=mom-g2, 2=dad-g1, 3=dad-g2, -1 = unknown;
    private byte[] gameteOriginOfMinorAllele; //0=mom-g1, 1=mom-g2, 2=dad-g1, 3=dad-g2, -1 = unknown;



    public PhaseF1TestPlugin() {
        super(null, false);
    }

    public PhaseF1TestPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    private void init() {
        Alignment a=ImportUtils.readFromHapmap(infiles[0]);
        myDataSBit=new SBitAlignment(a);
        System.out.printf("%s Sites: %d Taxa: %d %n",infiles[0],myDataSBit.getSiteCount(), myDataSBit.getSequenceCount());
//        myParentalChrComp=new byte[numParents][numChr][myDataSBit.getSiteCount()];
//        for (int i = 0; i < numParents; i++) {
//            for (int j = 0; j < numChr; j++) {
//                Arrays.fill(myParentalChrComp[i][j], DataType.UNKNOWN_BYTE);
//            }
//        }
        offspring=new byte[myDataSBit.getSequenceCount()][numParents][myDataSBit.getSiteCount()];
        for (int i = 0; i < myDataSBit.getSequenceCount(); i++) {
            for (int j = 0; j < numParents; j++) {
                Arrays.fill(offspring[i][j], Byte.MIN_VALUE);
//                for (int k = 0; k < myDataSBit.getSiteCount(); k++) {
//                    offspring[i][j][k]=(byte)((k/3)%2);
//                    System.out.print(offspring[i][j][k]+" ");
//                    k++;
//                    System.out.print(offspring[i][j][k]+" ");
//                }
//                System.out.println("");
            }
        }
        System.out.println(countCrossOvers(false));
 //       printOffspring();
    }

    private int[][] countCrossOvers(boolean printToScreen) {
        int[] coCnt=new int[numParents];
        int[][] alleleCnt=new int[numParents][numChr];
        for (int i = 0; i < myDataSBit.getSequenceCount(); i++) {
            
            for (int j = 0; j < numParents; j++) {
                int lineCOCnt=0;
                byte lastGoodChr=Byte.MIN_VALUE;
//                for (int k = 0; k < 300; k++) {
                for (int k = 0; k < myDataSBit.getSiteCount(); k++) {
                    if(offspring[i][j][k]>Byte.MIN_VALUE) {
                        alleleCnt[j][offspring[i][j][k]]++;
                        if((offspring[i][j][k]!=lastGoodChr)) {
                            if(lastGoodChr!=Byte.MIN_VALUE) lineCOCnt++;
                            lastGoodChr=offspring[i][j][k];
                        }
                    }
                }
                coCnt[j]+=lineCOCnt;
            }
        //    if(printToScreen) System.out.printf("%d, ",lineCOCnt);
            
        }
        if(printToScreen) System.out.println();
        if(printToScreen) System.out.println(Arrays.deepToString(alleleCnt));
        //System.out.println(Arrays.deepToString(alleleCnt));
        int[][] result=new int[2][2];
        for (int j = 0; j < numParents; j++) {
            result[j][0]=coCnt[j];
            result[j][1]=Math.min(alleleCnt[j][0], alleleCnt[j][1]);
        }
        return result;
    }

    private void findBestParentalCluster() {
        int bestParentalConfig=-1;
        double bestCORatio=Double.MIN_VALUE;
        long maxConfig=(1<<((2*numParentGametes)+1))-1;
        for (int i = 0; i <= maxConfig; i++) {
            int currConfig=i;
            for (int j = 0; j < numParentGametes; j++) {
                gameteOriginOfCluster[j]=(byte)(0x3&currConfig);
                currConfig=currConfig>>2;
            }
 //           gameteOriginOfCluster={0, 2, 3, 3, 1};
            initOffspring();
            int[][] currCOs=countCrossOvers(false);
            double avgMinorRat=0;
            for (int j = 0; j < numParents; j++) {
                double minorCoRatio=(double)currCOs[j][1]/(double)currCOs[j][0];  //this statistic is not working quite right
//                System.out.println(Arrays.toString(gameteOriginOfCluster));
//                System.out.printf(" %d P%d COCnt %d %d %g %n",i,j,currCOs[j][0], currCOs[j][1], minorCoRatio);
                avgMinorRat+=minorCoRatio;
            }
 //           if(i==504) printOffspring();
            avgMinorRat/=2.0;
            if(avgMinorRat>bestCORatio) {
                bestCORatio=avgMinorRat;
                bestParentalConfig=i;
                System.out.println("BEST "+Arrays.toString(gameteOriginOfCluster)+" i:"+i);
                System.out.printf("BEST %d P0 COCnt %d %d %g %n",i,currCOs[0][0], currCOs[0][1], avgMinorRat);
                System.out.printf("BEST %d P1 C0Cnt %d %d %g %n",i,currCOs[1][0], currCOs[1][1], avgMinorRat);
                if(avgMinorRat>20.0) printOffspring();
            }
        }

    }

    private void initOffspring() {
        for (int k = 0; k < myDataSBit.getSiteCount(); k++) {
            for (int i = 0; i < myDataSBit.getSequenceCount(); i++) { //reset all values
                offspring[i][0][k]=offspring[i][1][k]=Byte.MIN_VALUE;
            }
//            if(highLD[k]<minHighLDSites) continue;
//            byte mj=myDataSBit.getMajorAllele(k);
//            byte mn=myDataSBit.getMinorAllele(k);
            OpenBitSet theMinor=myDataSBit.getSiteBitsWithClone(k, 1);
            if(minorAlleleCluster[k]<0) continue; //no hypothesis for origin yet so skip
            int parentOfMinor=gameteOriginOfCluster[minorAlleleCluster[k]]/2;
            byte gameteOfMinor=(byte)(gameteOriginOfCluster[minorAlleleCluster[k]]%2);
            if(parentOfMinor<0) continue;
            for (int i = 0; i < myDataSBit.getSequenceCount(); i++) {
                if(theMinor.fastGet(i)) offspring[i][parentOfMinor][k]=gameteOfMinor;
            }
        }
    }

    private void updateOffspring(boolean random) {
        for (int k = 0; k < myDataSBit.getSiteCount(); k++) {
            if(highLD[k]<minHighLDSites) continue;
            byte mj=myDataSBit.getMajorAllele(k);
            byte mn=myDataSBit.getMinorAllele(k);
            OpenBitSet theMinor=myDataSBit.getSiteBitsWithClone(k, 1);
            int lowCoCnt=Integer.MAX_VALUE;
            int bestHypothesis=-1;
            for (int h = 0; h < 4; h++) {
                for (int i = 0; i < myDataSBit.getSequenceCount(); i++) {
                    if(theMinor.fastGet(i)) offspring[i][h/2][k]=(byte)(h%2);
                    //if the taxon has the minor it is assignable to are
                }
                int coCnt=countCrossOvers(false)[0][0];
//                System.out.printf("%d %d %d %n",k,h,coCnt);
//                printOffspring();
                if(coCnt<lowCoCnt) {
                    lowCoCnt=coCnt;
                    bestHypothesis=h;
                }
                for (int i = 0; i < myDataSBit.getSequenceCount(); i++) offspring[i][h/2][k]=Byte.MIN_VALUE;
            }
            if(random) bestHypothesis=k%4;
            for (int i = 0; i < myDataSBit.getSequenceCount(); i++) {
                if(theMinor.fastGet(i)) offspring[i][bestHypothesis/2][k]=(byte)(bestHypothesis%2);
                //if the taxon has the minor it is assignable to are
            }
        }
    }

    private float[][] getDMatrix() {
        FisherExact fe=new FisherExact(1000);
        float[][] d=new float[myDataSBit.getSiteCount()][myDataSBit.getSiteCount()];
        highLD=new int[myDataSBit.getSiteCount()];
        for (int i = 0; i < myDataSBit.getSiteCount(); i++) {
            OpenBitSet iMJ=myDataSBit.getSiteBitsWithClone(i, 0);
            OpenBitSet iMN=myDataSBit.getSiteBitsWithClone(i, 1);
            iMJ.andNot(iMN);  //make hets in homo minor
            for (int j = 0; j <= i; j++) {
                if(Math.abs(myDataSBit.getPositionInLocus(i)-myDataSBit.getPositionInLocus(j))<0) continue;
                OpenBitSet jMJ=myDataSBit.getSiteBitsWithClone(j, 0);
                OpenBitSet jMN=myDataSBit.getSiteBitsWithClone(j, 1);
                jMJ.andNot(jMN);  //make hets in homo minor
                int AB=(int)OpenBitSet.intersectionCount(iMJ, jMJ);
                int Ab=(int)OpenBitSet.intersectionCount(iMJ, jMN);
                int aB=(int)OpenBitSet.intersectionCount(iMN, jMJ);
                int ab=(int)OpenBitSet.intersectionCount(iMN, jMN);
                d[i][j]=d[j][i]=(float)calculateDPrime(AB,Ab,aB,ab,5);
                double p=fe.getCumlativeP(AB,Ab,aB,ab);
                if(p>minPValue) {d[i][j]=d[j][i]=(float)0;}
                else if(i!=j) {
                    highLD[i]++; highLD[j]++;
//                    if((i==8)||(j==8)) System.out.printf("%d %d %g %g%n",i,j,d[i][j],p);
                }
//                System.out.printf("%d %d %d %d %d %d %g %n", i,j,AB,Ab,aB,ab, d[i][j]);
            }
        }
//        for (int i = 0; i < 800; i++) {
//            if(highLD[i]<minHighLDSites) continue;
//            System.out.print(i+" "+highLD[i]+" ");
//            for (int j = 0; j <= i; j++) {
//                if(highLD[j]<minHighLDSites) continue;
//                System.out.printf("%g ",d[i][j]);
//            }
//            System.out.println();
//        }
        return d;
    }

    private DistanceMatrix convertDtoDistanceMatrix(float[][] ldD) {
        int highLDCnt=0;
        for (int i = 0; i < ldD.length; i++) {
            if(highLD[i]>=minHighLDSites) highLDCnt++;
        }
        double[][] dm=new double[highLDCnt][highLDCnt];
        IdGroup ids=new SimpleIdGroup(highLDCnt);
        int iIndex=0;
        for (int i = 0; i < ldD.length; i++) {
            if(highLD[i]<minHighLDSites) continue;
            ids.setIdentifier(iIndex, new Identifier(""+i));
            int jhIndex=0;
            for (int j = 0; j < i; j++) {
                if(highLD[j]<minHighLDSites) continue;
             //   if(ldD[i][j]>0) {dm[iIndex][jhIndex]=dm[jhIndex][iIndex]=1-ldD[i][j];}
                if(ldD[i][j]>0) {dm[iIndex][jhIndex]=dm[jhIndex][iIndex]=0;}
                else {dm[iIndex][jhIndex]=dm[jhIndex][iIndex]=1;}
                jhIndex++;
            }
            iIndex++;
        }
        DistanceMatrix theDM=new DistanceMatrix(dm,ids);
        return theDM;
    }

    private byte[] assignToGametesByTree(DistanceMatrix theDM) {
        minorAlleleCluster=new byte[myDataSBit.getSiteCount()];
        Arrays.fill(minorAlleleCluster, Byte.MIN_VALUE);
        NeighborJoiningTree theNJT=new NeighborJoiningTree(theDM);
        Node mid=TreeUtils.findMidpointNode(theNJT);
        System.out.println("Leaf Nodes:"+theNJT.getExternalNodeCount());
        System.out.printf("Mid %d %s %g %s %n",0,mid.getIdentifier().getName(),mid.getBranchLength(),
                    Arrays.toString(NodeUtils.getPathLengthInfo(mid)));
        theNJT.reroot(mid);
        ArrayList<Node[]> gametes=new ArrayList<Node[]>();
//        int nNodes=theNJT.getInternalNodeCount();
//        double[] sortedHeight = new double[nNodes];
//        for (int n = 0; n < nNodes; n++) {
//                sortedHeight[n] = theNJT.getInternalNode(n).getNodeHeight();
//        }
//        Arrays.sort(sortedHeight);
//        double maxHeight = sortedHeight[nNodes - 5];
        for (int i = 0; i < theNJT.getInternalNodeCount(); i++) {
            Node inode=theNJT.getInternalNode(i);
            if(inode.getBranchLength()>0.25) {
                System.out.println(NodeUtils.getLeafCount(inode));
                gametes.add(NodeUtils.getExternalNodes(inode));
//                System.out.printf("%d %s %g %s %n",i,inode.getIdentifier().getName(),inode.getBranchLength(),
//                    Arrays.toString(NodeUtils.getPathLengthInfo(inode)));
            }
        }
        this.numParentGametes=gametes.size();
        gameteOriginOfCluster=new byte[numParentGametes];
        Arrays.fill(gameteOriginOfCluster, (byte)-1);
        for (byte i = 0; i < numParentGametes; i++) {
            for (Node j:  gametes.get(i)) {
                int num=Integer.parseInt(j.getIdentifier().getName());
                if(minorAlleleCluster[num]>-1) System.out.println("Error node in two branches");
                minorAlleleCluster[num]=i;
            }
        }
        System.out.println(theNJT.toString());
        System.out.println(Arrays.toString(minorAlleleCluster));
        return minorAlleleCluster;
    }

    private void assignParent(float[][] d) {
        ArrayList<Integer>[] clusters=new ArrayList[30];
       // for(ArrayList al: clusters) al=new ArrayList();
        int currentClusters=0;
        for (int i = 0; i < myDataSBit.getSiteCount(); i++) {
             if(highLD[i]<minHighLDSites) continue;
             //find best cluster
             int bCluster=-1;
             float bD=0.5f;
             for (int j = 0; j < currentClusters; j++) {
               // int ls=clusters[j].get(clusters[j].size()-1);
                 for (int k = 0; k < clusters[j].size(); k++) {
                     int ls=clusters[j].get(k);
                     if(d[i][ls]>bD) {
                        bD=d[i][ls];
                        bCluster=j;
                        }
                 }
 //               if(currentClusters==8) System.out.printf("S%d C%d j%d ls%d d%g %n",i,currentClusters, j, ls, d[i][ls]);               
            }
            if(bCluster==-1) {
                clusters[currentClusters]=new ArrayList<Integer>();
                clusters[currentClusters].add(i);
//                System.out.printf("S%d C%d B%g %n",i,currentClusters, bD);
                currentClusters++;
            } else {
                clusters[bCluster].add(i);
 //               System.out.printf("S%d C%d %n",i,bCluster);
            }
        }
        for (int c1 = 0; c1 < currentClusters; c1++) {
            System.out.println(c1+":"+clusters[c1].toString());
        }
        for (int c1 = 0; c1 < currentClusters; c1++) {
            for (int c2 = 0; c2 < c1; c2++) {
                int cnt=0;
                double avg=0;
                for(int a1: clusters[c1]) {
                    for(int a2: clusters[c2]) {
                        if(Math.abs(a1-a2)>100) continue;
                        cnt++;
                        avg+=d[a1][a2];
                    }
                }
                avg=avg/(double)cnt;
                System.out.printf("C%d C%d Cnt:%d AvgD:%g %n", c1, c2, cnt, avg);
                if(avg>0.5) {
                    clusters[c2].addAll(clusters[c1]);
                    clusters[c1].clear();
                }
            }
        }
        for (int c1 = 0; c1 < currentClusters; c1++) {
            System.out.println(c1+":"+clusters[c1].toString());
        }
    }

    private void reportLDOfClusters(float[][] d) {
       double[][] clusterCorr=new double[numParentGametes][numParentGametes];
       int[][] clusterContrastCnt=new int[numParentGametes][numParentGametes];
        for (int i = 0; i < minorAlleleCluster.length; i++) {
            for (int j = 0; j < minorAlleleCluster.length; j++) {
                if(Math.abs(i-j)>200) continue;
                if((minorAlleleCluster[i]<0)||(minorAlleleCluster[j]<0)) continue;
                clusterCorr[minorAlleleCluster[i]][minorAlleleCluster[j]]+=d[i][j];
                clusterContrastCnt[minorAlleleCluster[i]][minorAlleleCluster[j]]++;
            }

        }
        for (int i = 0; i < clusterCorr.length; i++) {System.out.print(i+"\t");}
        System.out.println("");
        for (int i = 0; i < clusterCorr.length; i++) {
            System.out.print(i+"\t");
            for (int j = 0; j <=i; j++) {
                clusterCorr[j][i]=clusterCorr[i][j]=clusterCorr[i][j]/(double)clusterContrastCnt[i][j];
                System.out.printf("%.2g\t",clusterCorr[i][j]);
  //              System.out.print(clusterCorr[i][j]+"\t");
            }
            System.out.println("");
        }
 //       System.out.println(Arrays.deepToString(clusterCorr));

    }

//    private void reportClusterDCorrelation(float[][] d) {
//        for (int c1 = 0; c1 < numParentGametes; c1++) {
//            for (int c2 = 0; c2 < c1; c2++) {
//                int cnt=0;
//                double avg=0;
//                for(int a1: clusters[c1]) {
//                    for(int a2: clusters[c2]) {
//                        if(Math.abs(a1-a2)>100) continue;
//                        cnt++;
//                        avg+=d[a1][a2];
//                    }
//                }
//                avg=avg/(double)cnt;
//                System.out.printf("C%d C%d Cnt:%d AvgD:%g %n", c1, c2, cnt, avg);
//                if(avg>0.5) {
//                    clusters[c2].addAll(clusters[c1]);
//                    clusters[c1].clear();
//                }
//            }
//        }
//    }

    static double calculateDPrime(int countAB, int countAb, int countaB, int countab, int minTaxaForEstimate) {
        //this is the normalized D' is Weir Genetic Data Analysis II 1986 p120
        double freqR, freqC, freq, countR, countC, nonmissingSampleSize;
        nonmissingSampleSize = countAB + countAb + countaB + countab;
        if (nonmissingSampleSize < minTaxaForEstimate) {
            return Double.NaN;
        }
        countR = countab + countAb;
        countC = countab + countaB;
        freqR = (nonmissingSampleSize - countR) / nonmissingSampleSize;
        freqC = (nonmissingSampleSize - countC) / nonmissingSampleSize;
        if ((freqR == 0) || (freqC == 0) || (freqR == 1) || (freqC == 1)) {
            return Double.NaN;
        }
        freq = ((double) countAB / nonmissingSampleSize) - (freqR * freqC);
        if (freq < 0) {
            return -(freq / Math.max(-freqR * freqC, -(1 - freqR) * (1 - freqC)));
        } else {
            return freq / Math.min((1 - freqR) * freqC, (1 - freqC) * freqR);
        }  //check these equations
    }

    private void printOffspring() {
         for (int s = 0; s < highLD.length; s++) {
             if(highLD[s]<minHighLDSites) continue;
             System.out.printf("%4d %d: ",s,minorAlleleCluster[s]);
             for (int t = 0; t < myDataSBit.getSequenceCount(); t++) {
                 System.out.printf("%s:%s ",nicePrint(offspring[t][0][s]),nicePrint(offspring[t][1][s]));
             }
             System.out.println();
         }
    }

    private String nicePrint(byte b) {
        if(b>Byte.MIN_VALUE) return ""+b;
        return "*";
    }


    @Override
    public DataSet performFunction(DataSet input) {
        init();
        float[][] d=getDMatrix();

        DistanceMatrix dm=convertDtoDistanceMatrix(d);
        assignToGametesByTree(dm);
        reportLDOfClusters(d);
        findBestParentalCluster();
  //      printOffspring();
  //      assignParent(d);
  //      updateOffspring(true);
//        for (int i = 0; i < 4; i++) {
//            updateOffspring(false);
//            System.out.println("Round"+i+" CO:"+countCrossOvers());
//            printOffspring();
//        }
        return null;
    }

    public void setParameters(String[] args){
        if (myArgsEngine == null) {
            myArgsEngine = new ArgsEngine();
            myArgsEngine.add("-i", "--input-file", true);
//            myArgsEngine.add("-o", "--output_file", true);
        }
        myArgsEngine.parse(args);
        if (myArgsEngine.getBoolean("-i")) {
            infiles=new String[1];
            infiles[0]=myArgsEngine.getString("-i");
        }

    }
    private void printUsage() {
            myLogger.info(
                "\nUsage is as follows:\n"
                + "-i  Input directory containing .tbt.bin files\n"
                + "-o  Output file name\n"
            );
    }

    public static void main(String[] args) {
        String[] tArgs={"-i","/Users/edbuckler/SolexaAnal/grapes/test/Vitis_062111.c+.hmp.txt"};
        args=tArgs.clone();
        PhaseF1TestPlugin thePFTP=new PhaseF1TestPlugin();
        for (int i = 1; i < 5; i++) {
            String x=tArgs[1].replace("+", ""+i);
            System.out.println(i+":"+x);
            args[1]=x;
            thePFTP.setParameters(args);
            thePFTP.performFunction(null);

        }
        
        

        
    }
    
    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "PhaseF1";
    }

   @Override
    public String getToolTipText() {
        return "PhaseF1";
    }

}
