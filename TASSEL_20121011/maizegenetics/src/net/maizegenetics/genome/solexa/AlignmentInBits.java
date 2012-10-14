/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.solexa;

import java.util.Arrays;
import net.maizegenetics.pal.alignment.Alignment;

/**
 *
 * @author edbuckler
 */
public class AlignmentInBits {
    private int numTaxa;
    private int numSites;
    long[] seq;
    short[] minorFreq;
    int lgPerSite;
    int[] position;
    int chr;
    int getIndex=0; //it should also be at the start of major allele, add lgPerSite to get to the start of minor

    public AlignmentInBits(Alignment a) {
        numTaxa=a.getSequenceCount();
        numSites=a.getSiteCount();
        lgPerSite=(numTaxa/64) + 1;
        position=new int[numSites];
        minorFreq=new short[numSites];
        chr=Integer.parseInt(a.getLocusName(0));
        int seqSize=numSites*2*lgPerSite;
        seq=new long[seqSize];
        for (int i = 0; i < a.getSiteCount(); i++) {
            position[i]=a.getPositionInLocus(i);
            byte minorAllele=(byte)a.getMinorAllele(i);
            byte majorAllele=(byte)a.getMajorAllele(i);
            int minorCnt=0, majorCnt=0;
            for(int j=0; j<a.getSequenceCount(); j++) {
                if(minorAllele==a.getBase(j, i)) {
                    setSeqTrue(j,i,false);
                    minorCnt++;
                } else
                if(majorAllele==a.getBase(j, i)) {
                    setSeqTrue(j,i,true);
                    majorCnt++;
                }
            }
            minorFreq[i]=(minorCnt<majorCnt)?(short)minorCnt:(short)majorCnt;
//            if(minorFreq[i]==0) System.out.println("Minor zero:"+i+" pos:"+position[i]);
        }
    }

    public int getMajorAllelesCnt(int taxon) {
        int cnt=0;
        for (int i = 0; i < numSites; i++) {
            if(isMajorAllele(i,taxon)) cnt++;
        }
        return cnt;
    }

    public int getMinorAllelesCnt(int taxon) {
        int cnt=0;
        for (int i = 0; i < numSites; i++) {
            if(isMinorAllele(i,taxon)) cnt++;
        }
        return cnt;
    }

    private void setSeqTrue(int taxon, int site, boolean majorAllele) {
        int index=(2*site*lgPerSite)+(taxon/64);
        if(majorAllele==false) index+=lgPerSite;
        seq[index]=seq[index]|(1L<<(taxon%64));
    }

    public long[] getMajorAlleles(int site) {
        long[] result=new long[lgPerSite];
        getIndex=(2*site*lgPerSite);
        int index=getIndex;
        for (int i = 0; i < lgPerSite; i++) {
            result[i]=seq[index++];
        }
        return result;
    }

    public boolean isMajorAllele(int site, int taxon) {
        int index=(2*site*lgPerSite)+(taxon/64);
        return getBit(seq[index],taxon%64);
    }

    public boolean isMinorAllele(int site, int taxon) {
        int index=(2*site*lgPerSite)+lgPerSite+(taxon/64);
        return getBit(seq[index],taxon%64);
    }

    public long getCurrentMajorAlleles(int offset) {
        return seq[getIndex+offset];
    }

    public long[] getCurrentMajorAllelesx() {
        long[] result={seq[getIndex]};
        return result;
    }

    public long[] getCurrentMajorAlleles() {
        long[] result=new long[lgPerSite];
        for (int i = 0; i < lgPerSite; i++) {
            result[i]=seq[getIndex+i];
        }
        return result;
    }

    public long getCurrentMinorAlleles(int offset) {
        return seq[getIndex+lgPerSite+offset];
    }

    public long[] getCurrentMinorAllelesx() {
        long[] result={seq[getIndex+lgPerSite]};
        return result;
    }

    public long[] getCurrentMinorAlleles() {
        long[] result=new long[lgPerSite];
        for (int i = 0; i < lgPerSite; i++) {
            result[i]=seq[getIndex+lgPerSite+i];
        }
        return result;
    }

     public long[] getMinorAlleles(int site) {
        long[] result=new long[lgPerSite];
        getIndex=(2*site*lgPerSite);
        int index=getIndex+lgPerSite;
        for (int i = 0; i < lgPerSite; i++) {
            result[i]=seq[index++];
        }
        return result;
    }

    public long[][] getAlleles(int site) {
        long[][] result=new long[2][lgPerSite];
        getIndex=(2*site*lgPerSite);
        int index=getIndex;
        for (int i = 0; i < lgPerSite; i++) {
            result[0][i]=seq[index];
            result[1][i]=seq[index+lgPerSite];
            index++;
        }
        return result;
    }

    public long[][] getCurrentAlleles() {
        long[][] result=new long[2][lgPerSite];
        
        int index=getIndex;
        for (int i = 0; i < lgPerSite; i++) {
            result[0][i]=seq[index];
            result[1][i]=seq[index+lgPerSite];
            index++;
        }
        return result;
    }

    public void nextSite() {
        getIndex=getIndex+lgPerSite+lgPerSite;
    }

    public long[][] nextAlleles() {
        getIndex=getIndex+lgPerSite+lgPerSite;
        return getCurrentAlleles();
    }

    boolean getBit(long set, int bit) { // get bit 'bit' from the set
       return (set&(1L<<bit)) != 0;
    }

    public void setCurrentSite(int site) {
        getIndex=(2*site*lgPerSite);
    }

    public int getCurrentSite() {
        return getIndex/(2*lgPerSite);
    }

    long setBit(long set, int bit, boolean value) { // return the new set with bit 'bit' (un)set
       if (value)
          return set|(1L<<bit);
       else
          return set&~(1<<bit);
    }


    public int getNumSites() {
        return numSites;
    }

    public int getNumTaxa() {
        return numTaxa;
    }

    public int getChr() {
        return chr;
    }

    public int getPosition(int site) {
        return position[site];
    }

    public int getClosestSite(int testPosition) {
        int hit=Arrays.binarySearch(position, testPosition);
        if(hit>-1) return hit;
        int insert=-(hit+1);
        if(insert==0) return insert;
        if(insert>=position.length) return (position.length-1);
        if((position[insert]-testPosition)<(testPosition-position[insert-1])) return insert;
        return insert-1;
    }

    public int getMinorFreq(int site) {
        return minorFreq[site];
    }


}
