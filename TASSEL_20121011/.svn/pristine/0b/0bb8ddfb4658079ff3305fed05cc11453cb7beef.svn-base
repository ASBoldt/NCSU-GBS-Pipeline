/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS;

import java.nio.ByteBuffer;
import java.util.Arrays;
import net.maizegenetics.genome.BaseEncoder;

/**
 *
 * @author edbuckler
 */
public class ReadWithPositionInfo implements Comparable {
 //   long[] readSet=new long[2];  //the read
    static final byte unknownChr=-1;
    long read0, read1;
    byte chromosome;
    //=unknownChr[0], chromosome1=unknownChr[1];//=new byte[2];  //storing a two character string in bytes
    byte strand; //'+', '-',"?"
    int positionMin, positionMax;  //reference genome position for start and end of read
            //if the read is only mapped genetically this can be much larger than the size of the read
            //to indicate the abiguity
            //unknown is Integer.MIN_VALUE
    short nextCutDistance;  //distance between known cut sites, missing Short.MIN_VALUE
    byte divergence;  //number of diverging bp from reference, unknown Byte.MIN_VALUE
    float dcoP, mapP;  //p-values for various mapping approaches, unknown NaN
            //if these disagree with the location, then set the p to negative
    byte multimaps;  //number of locations this readSet maps to, unknown Byte.MIN_VALUE
    public static final int totalByteSize=38;
    static final int read0Off=0, read1Off=read0Off+8, chromosomeOff=read1Off+8,
            strandOff=chromosomeOff+1, positionMinOff=strandOff+1, positionMaxOff=positionMinOff+4,
            nextCutDistanceOff=positionMaxOff+4, divergenceOff=nextCutDistanceOff+2, dcoPOff=divergenceOff+1,
            mapPOff=dcoPOff+4, multimapsOff=mapPOff+4;
  //  Enum test {Copy1PerfMGenetic+, Copy1PerfMGenetic0, Copy1PerfMGenetic-};



    ReadWithPositionInfo (long[] haplotype, byte chromosome, byte strand, int positionMin,
            int positionMax, short nextCutDistance, byte divergence, float dcoP, float mapP,
            byte multimaps) {
//        this.read0=haplotype[0];
//        this.read1=haplotype[1];
        this.read0=haplotype[0];
        this.read1=haplotype[1];
        this.chromosome=chromosome;
//        if(chromosome.length()==0) {this.chromosome="  ".getBytes();}
//        else if(chromosome.length()==1) {this.chromosome=(chromosome.charAt(0)+" ").getBytes();}
//        else {this.chromosome=chromosome.substring(0, 2).getBytes();}
        this.strand=strand;
        this.positionMin=positionMin;
        this.positionMax=positionMax;
        this.nextCutDistance=nextCutDistance;
        this.divergence=divergence;
        this.dcoP=dcoP;
        this.mapP=mapP;
        this.multimaps=multimaps;
    }

    ReadWithPositionInfo (ByteBuffer bb) {
        read0=bb.getLong(read0Off);
        read1=bb.getLong(read1Off);
        chromosome=bb.get(chromosomeOff);
        strand=bb.get(strandOff);
        positionMin=bb.getInt(positionMinOff);
        positionMax=bb.getInt(positionMaxOff);
        nextCutDistance=bb.getShort(nextCutDistanceOff);
        divergence=bb.get(divergenceOff);
        dcoP=bb.getFloat(dcoPOff);
        mapP=bb.getFloat(mapPOff);
        multimaps=bb.get(multimapsOff);
    }

    /**
     * Constructor for the physical map digest
     */
    ReadWithPositionInfo (long[] haplotype, byte chromosome, byte strand, int positionMin,
            int positionMax, short nextCutDistance, byte divergence) {
        this(haplotype, chromosome, strand, positionMin, positionMax, nextCutDistance,
                divergence, Float.NaN, Float.NaN, Byte.MIN_VALUE);
    }

    /**
     * Constructor when nothing is known except the sequence
     * @param haplotype
     */
    ReadWithPositionInfo (long[] haplotype) {
        this(haplotype, unknownChr, (byte)'?', Integer.MIN_VALUE, Integer.MIN_VALUE, Short.MIN_VALUE,
                Byte.MIN_VALUE, Float.NaN, Float.NaN, Byte.MIN_VALUE);
    }
    /**
     * Constructor when nothing is known except the physical position
     */
    ReadWithPositionInfo (byte chr, byte strand, int posMin) {
        this(BaseEncoder.getLongArrayFromSeq("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
                chr, strand, posMin, posMin+63, Short.MIN_VALUE, Byte.MIN_VALUE, Float.NaN, Float.NaN, Byte.MIN_VALUE);
    }

    @Override
    public String toString() {
        String s=BaseEncoder.getSequenceFromLong(read0)+
                    BaseEncoder.getSequenceFromLong(read1)+
                    "\t"+
                   chromosome+"\t"+(char)strand+"\t"+positionMin+"\t"+
                   positionMax+"\t"+nextCutDistance+"\t"+divergence+"\t"+
                   dcoP+"\t"+ mapP +"\t"+ multimaps;
        return s;
    }

    public void setPosition(ReadWithPositionInfo newPosition) {
        this.chromosome=newPosition.chromosome;
        this.positionMin=newPosition.positionMin;
        this.positionMax=newPosition.positionMax;
        this.strand=newPosition.strand;
    }

    public byte[] toByte() {
        byte[] theResult=new byte[totalByteSize];
        ByteBuffer fw=ByteBuffer.wrap(theResult);
        fw.putLong(read0);
        fw.putLong(read1);
        fw.put(chromosome);
        fw.put(strand);
        fw.putInt(positionMin);
        fw.putInt(positionMax);
        fw.putShort(nextCutDistance);
        fw.put(divergence);
        fw.putFloat(dcoP);
        fw.putFloat(mapP);
        fw.put(multimaps);
        return theResult;
    }


    public boolean equals(Object obj) {
        final ReadWithPositionInfo other = (ReadWithPositionInfo) obj;
        if (this.read0!=other.read0) {
            return false;
        }
        if (this.read1!=other.read1) {
            return false;
        }
        if (this.chromosome != other.chromosome) {
            return false;
        }
        if (this.positionMin != other.positionMin) {
            return false;
        }
        return true;
    }

    public int compareTo(Object o) {
            ReadWithPositionInfo other = (ReadWithPositionInfo)o;
            if(this.equals(other)) return 0;
            if (read0 < other.read0) {
                    return -1;
            } else if(read0 > other.read0) return 1;
            if (read1 < other.read1) {
                    return -1;
            } else if (read1 > other.read1) {
                    return 1;
            }
            if (this.chromosome < other.chromosome) {
                return -1;
            } else if (this.chromosome > other.chromosome) {
                return 1;
            }
            if (this.positionMin < other.positionMin) {
                return -1;
            } else if  (this.positionMin > other.positionMin) {
                return 1;
            }
            return 0;
	}

    public long[] getRead() {
        long[] result={read0, read1};
        return result;
    }

    public byte getChromosome() {
        return chromosome;
    }

    public int getPositionMin() {
        return positionMin;
    }

    public int getPositionMax() {
        return positionMax;
    }

    public byte getStrand() {
        return strand;
    }

    public short getNextCutDistance() {
        return nextCutDistance;
    }

    public boolean isChromosomeKnown() {
        if(chromosome<0) return false;
        return true;
    }

    public float getDcoP() {
        return dcoP;
    }

    public void setDcoP(float dcoP) {
        this.dcoP = dcoP;
    }

    public byte getDivergence() {
        return divergence;
    }

    public void setDivergence(int divergence) {
        if(divergence>Byte.MAX_VALUE) {this.divergence=Byte.MAX_VALUE;}
        else {this.divergence = (byte)divergence;}
    }

    public float getMapP() {
        return mapP;
    }

    public void setMapP(float mapP) {
        this.mapP = mapP;
    }

    public byte getMultimaps() {
        return multimaps;
    }

    public void setMultimaps(int multimaps) {
        if(multimaps>Byte.MAX_VALUE) {this.multimaps=Byte.MAX_VALUE;}
        else {this.multimaps = (byte)multimaps;}
    }

}
