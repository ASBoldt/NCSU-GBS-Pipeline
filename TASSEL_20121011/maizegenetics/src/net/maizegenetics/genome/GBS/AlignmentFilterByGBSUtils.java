/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS;

import cern.colt.list.IntArrayList;
import java.util.Arrays;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.GdpdmBLOBUtils;
import net.maizegenetics.pal.datatype.DataType;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.IdGroupUtils;
import net.maizegenetics.pal.popgen.LinkageDisequilibrium;
import org.apache.poi.util.IntList;
//import org.apache.poi.util.IntList;

/**
 *
 * @author edbuckler
 */
public class AlignmentFilterByGBSUtils {
    public static final byte refAllele=GdpdmBLOBUtils.bases[0]; //A
    public static final byte altAllele=GdpdmBLOBUtils.bases[1]; //C
    public static final byte hetAllele=GdpdmBLOBUtils.bases[9]; //M

    private AlignmentFilterByGBSUtils() {
    }


    public static int[][] hetsByLine(Alignment a, boolean printToScreen) {
        int[][] counts=new int[2][a.getSequenceCount()];
        int totalScored=0, totalHets=0;
        for(int j=0; j<a.getSequenceCount(); j++) {
            for (int i = 0; i < a.getSiteCount(); i++) {
                if(a.getBase(j, i)!=DataType.UNKNOWN_BYTE) {
                    counts[0][j]++;
                    if(a.getBase(j, i)==hetAllele) counts[1][j]++;
                }
            }
//            if(printToScreen) System.out.println(a.getIdGroup().getIdentifier(j).getName()+"\t"+
//                    counts[0][j]+"\t"+counts[1][j]);
        }
        for(int c: counts[0]) totalScored+=c;
        for(int c: counts[1]) totalHets+=c;
        if(printToScreen) System.out.println("Total Alleles:"+totalScored+" TotalHets:"+totalHets);
        return counts;
    }

    public static IdGroup getLowHetIdGroup(Alignment a, double maxHets, int minCount) {
        int[][] hetCnt=hetsByLine(a,false);
        boolean[] include=new boolean[a.getSequenceCount()];
        for (int i = 0; i < hetCnt[0].length; i++) {
            if(((double)hetCnt[1][i]/(double)hetCnt[0][i]>maxHets)||(hetCnt[0][i]<minCount)) {include[i]=false;}
            else {include[i]=true;}
        }
        return IdGroupUtils.idGroupSubset(a.getIdGroup(), include);
    }

     public static int[][] hetsBySite(Alignment a, boolean printToScreen) {
        int[][] counts=new int[2][a.getSiteCount()];  //total count in row 0, hets count in row 1
        if(printToScreen) System.out.println("Locus\tMAF\tLineScored\tHetNum\tHetRate");
        for (int i = 0; i < a.getSiteCount(); i++) {
            for(int j=0; j<a.getSequenceCount(); j++) {
                if(a.getBase(j, i)!=DataType.UNKNOWN_BYTE) {
                    counts[0][i]++;
                    if(a.getBase(j, i)==hetAllele) counts[1][i]++;
                }
            }
            if(printToScreen) System.out.println(a.getSNPID(i)+"\t"+a.getMinorAlleleFrequency(i)+"\t"+
                    counts[0][i]+"\t"+counts[1][i]+"\t"+(double)counts[1][i]/counts[0][i]);
        }
       //     System.out.println("counts"+Arrays.deepToString(counts));
        return counts;
    }

     public static int[] getLowHetSNPs(Alignment a, double maxHets, int minCount) {
        int[][] hetCnt=hetsBySite(a,false);
        IntArrayList goodSites=new IntArrayList();
        for (int i = 0; i < hetCnt[0].length; i++) {
            if(((double)hetCnt[1][i]/(double)hetCnt[0][i]<maxHets)&&
                    (hetCnt[0][i]>=minCount)) {goodSites.add(i);}
        }
        goodSites.trimToSize();
        return goodSites.elements();
    }

    public static int[] getGoodSitesByLD(Alignment a, double minR2, int windowSize,
            int minCnt, boolean keepUnproven) {
        IntArrayList goodSites=new IntArrayList();
        LinkageDisequilibrium theLD=new LinkageDisequilibrium(a,
            minCnt, windowSize, LinkageDisequilibrium.testDesign.SlidingWindow);
        theLD.run();
        for (int i = 0; i < a.getSiteCount(); i++) {
            int cntInformative=0;
            double obsMaxR2=-1;
            double obsMinP=1;
            for(int j=i-(windowSize/2); j<i+(windowSize/2); j++) {
            //for(int j=i-1; j<i; j++) {
    //got a real problem r2(i=j) is not NAN or 1.0

    //       if(i==j)     System.out.println(i+","+j+" r2:"+theLD.getRSqr(i, j));
                if(Double.isNaN(theLD.getRSqr(i, j))==false) {
                    cntInformative++;
                    if(theLD.getRSqr(i, j)>obsMaxR2) obsMaxR2=theLD.getRSqr(i, j);
                    if(theLD.getPVal(i, j)<obsMinP) obsMinP=theLD.getPVal(i, j);
                }
            }
 //           if(obsMaxR2<minR2) System.out.println(i+" "+cntInformative+" "+obsMaxR2+" "+obsMinP);
            if((obsMaxR2>minR2)||(keepUnproven&&(cntInformative==0))) goodSites.add(i);
        }
        goodSites.trimToSize();
        return goodSites.elements();
    }

    public static int[][] countCrossoversByLine(Alignment a) {
        int[][] cos=new int[2][a.getSequenceCount()];  //total count of useful markers in
        //row 0, crossovers in row 1
        int sumGood=0,sumCO=0;
        for(int t=0; t<a.getSequenceCount(); t++) {
            byte lastHomozygous=DataType.UNKNOWN_CHARACTER;
            for (int i = 0; i < a.getSiteCount(); i++) {
                byte currBase=a.getBase(t, i);
                if((currBase!=refAllele)&&(currBase!=altAllele)) continue;  //not useful
                cos[0][t]++;
                if(currBase!=lastHomozygous) {cos[1][t]++;}
                lastHomozygous=currBase;
            }
            sumCO+=cos[1][t];
            sumGood+=cos[0][t];
          //  System.out.println(a.getIdGroup().getIdentifier(t).getName()+" "+cos[0][t]+" "+cos[1][t]);
        }
        System.out.println("TotalHomoMarkers:"+sumGood+" TotalCrossover:"+sumCO);
//        System.out.println("By Line:"+Arrays.toString(cos[0]));
//        System.out.println("By Line:"+Arrays.toString(cos[1]));
        return cos;
    }


     public static int[][] countDCO(Alignment a, boolean byTaxa) {
        int[][] cos;
        if(byTaxa) {cos=new int[2][a.getSequenceCount()];}
        else {cos=new int[2][a.getSiteCount()];}
        int sumGood=0,sumDCO=0;
        for (int t = 0; t < a.getSequenceCount(); t++) {
          byte[] base={-1,-1,-1};
          int[] site={-1,-1,-1};
          for(int s = 0; s <a.getSiteCount(); s++) {
            if((a.getBase(t,s)==refAllele)||(a.getBase(t,s)==altAllele)) {
                base[0]=base[1]; base[1]=base[2]; base[2]=a.getBase(t,s);
                site[0]=site[1]; site[1]=site[2]; site[2]=s;
                if(base[0]==-1) continue;  //need to load the arrays first before checking
                sumGood++;
                if(byTaxa) {cos[0][t]++;} else {cos[0][site[1]]++;}
                if((base[0]==base[2])&&(base[0]!=base[1])) {
                    sumDCO++;
                    if(byTaxa) {cos[1][t]++;} else {cos[1][site[1]]++;}
                }
            }
          }
 //           System.out.println(s+" "+cos[0][s]+" "+cos[1][s]);
        }
        if(byTaxa) {System.out.print("ByTaxa:");} else {System.out.print("BySite:");}
        System.out.println("TotalHomoMarkers:"+sumGood+" TotalDoubleCrossover:"+sumDCO);
      //  System.out.println("By Line:"+Arrays.toString(cos[1]));
        return cos;
    }

    public static int[] getLowDCOSNPs(Alignment a, double maxDCOrate, int minCount) {
        int[][] dcoCnt=countDCO(a,false);
        IntList goodSites=new IntList();
        for (int i = 0; i < dcoCnt[0].length; i++) {
            if(((double)dcoCnt[1][i]/(double)dcoCnt[0][i]<maxDCOrate)&&
                    (dcoCnt[0][i]>=minCount)) {goodSites.add(i);}
        }
        return goodSites.toArray();
    }

     public static IdGroup getLowDCOIdGroup(Alignment a, double maxDCO, int minCount) {
        int[][] dcoCnt=hetsByLine(a,false);
        boolean[] include=new boolean[a.getSequenceCount()];
        for (int i = 0; i < dcoCnt[0].length; i++) {
            if(((double)dcoCnt[1][i]/(double)dcoCnt[0][i]>maxDCO)||(dcoCnt[0][i]<minCount)) {include[i]=false;}
            else {include[i]=true;}
        }
        return IdGroupUtils.idGroupSubset(a.getIdGroup(), include);
    }
}
