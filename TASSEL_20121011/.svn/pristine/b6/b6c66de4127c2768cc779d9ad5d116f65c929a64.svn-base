/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS;

import java.util.ArrayList;
import java.util.Arrays;
import net.maizegenetics.genome.BaseEncoder;
import net.maizegenetics.pal.statistics.FisherExact;

/**
 *
 * @author edbuckler
 */
public class AlleleStructure {
    private ArrayList<Allele> theAlleles=new ArrayList<Allele>();
    private int[][][] alleleSharing=null;

    public AlleleStructure() {
    }

    public void addAllele(long[] read, byte[] byTaxa) {
        theAlleles.add(new Allele(read, byTaxa));
    }

    public void addAllele(Allele aAllele) {
        theAlleles.add(aAllele);
    }

    public int getNumberOfAlleles() {
        return theAlleles.size();
    }

    public String toString() {
        FisherExact fe=new FisherExact(1200);
        if(alleleSharing==null) calcSharing();
        StringBuilder sb=new StringBuilder();
        for(Allele al: theAlleles) sb.append(al.toString()+"\n");
        for (int i = 0; i < alleleSharing.length; i++) {
            for(int j=0; j<i; j++) {
                double fes=fe.getCumlativeP(alleleSharing[i][j][0],alleleSharing[i][j][1],alleleSharing[i][j][2],alleleSharing[i][j][3]);

                sb.append(i+"(N="+theAlleles.get(i).sum+"),"+j+"(N="+theAlleles.get(j).sum+")="+
                        Arrays.toString(alleleSharing[i][j])+" P="+fes);
                sb.append(" D="+calcD(alleleSharing[i][j]));
                sb.append(" D'="+calculateDPrime(alleleSharing[i][j]));
                sb.append("\n");
               // double fes=fe.getCumlativeP(alleleSharing[i][j][0],alleleSharing[i][j][1],alleleSharing[i][j][2],alleleSharing[i][j][3]);
            }
        }
        return sb.toString();
    }

    private void calcSharing() {
        alleleSharing=new int[getNumberOfAlleles()][getNumberOfAlleles()][];
        for (int i = 0; i < getNumberOfAlleles(); i++) {
            for(int j=0; j<i; j++) {
                alleleSharing[i][j]=alleleSharing[j][i]=calcAlleleSharingPA(theAlleles.get(i),theAlleles.get(j));
            }
        }
    }

    int[] calcAlleleSharing(Allele a1, Allele a2) {
        int[] share=new int[4];  // i0=missing both, i1=present in 1, i2=present in 2, i3=present in both
        for (int i = 0; i < a1.byTaxa.length; i++) {
            int index=0;
            if(a1.byTaxa[i]>0) index+=1;
            if(a2.byTaxa[i]>0) index+=2;
            share[index]+=a1.byTaxa[i]+a2.byTaxa[i];
        }
        return share;
    }

    int[] calcAlleleSharingPA(Allele a1, Allele a2) {
        int[] share=new int[4];  // i0=missing both, i1=present in 1, i2=present in 2, i3=present in both
        for (int i = 0; i < a1.byTaxa.length; i++) {
            int index=0;
            if(a1.byTaxa[i]>0) index+=1;
            if(a2.byTaxa[i]>0) index+=2;
            share[index]++;
        }
        return share;
    }

    double calcD(int[] ac) {
        double sum=ac[0]+ac[1]+ac[2]+ac[3];
        double A2freq=(ac[1]+ac[3])/sum;
        double B2freq=(ac[2]+ac[3])/sum;
        double d=(ac[3]/sum)-A2freq*B2freq;
        return d;
    }

   double calculateDPrime(int[] ac) {
       int countAB=ac[0], countAb=ac[1], countaB=ac[2], countab=ac[3];
        //this is the normalized D' is Weir Genetic Data Analysis II 1986 p120
        double freqR, freqC, freq, countR, countC, nonmissingSampleSize;
        nonmissingSampleSize = countAB + countAb + countaB + countab;
        if (nonmissingSampleSize <= 10) {
            return Double.NaN;
        }
        countR = countab + countAb;
        countC = countab + countaB;
        freqR = (nonmissingSampleSize - countR) / nonmissingSampleSize;
        freqC = (nonmissingSampleSize - countC) / nonmissingSampleSize;
        // if((freqR==0)||(freqC==0)||(freqR==1)||(freqC==1)) return -999;  //changed by ed 8-13-2004
        if ((freqR == 0) || (freqC == 0) || (freqR == 1) || (freqC == 1)) {
            return Double.NaN;
        }
        freq = ((double) countAB / nonmissingSampleSize) - (freqR * freqC);
        if (freq < 0) {
            return freq / Math.max(-freqR * freqC, -(1 - freqR) * (1 - freqC));
        } else {
            return freq / Math.min((1 - freqR) * freqC, (1 - freqC) * freqR);
        }  //check these equations
    }
}



class Allele {
    long[] read;
    byte[] byTaxa;
    int sum;

    public Allele(long[] read, byte[] byTaxa) {
        this.read = read;
        this.byTaxa = byTaxa;
        sum=0;
        for(int t: byTaxa) {
            sum+=t;
        }
    }

    public String toString() {
        StringBuilder sb=new StringBuilder();
        sb.append(BaseEncoder.getSequenceFromLong(read)+"\t");
        for (int i = 0; i < byTaxa.length; i++) {
            sb.append(byTaxa[i]+" ");
        }
        return sb.toString();
    }
}