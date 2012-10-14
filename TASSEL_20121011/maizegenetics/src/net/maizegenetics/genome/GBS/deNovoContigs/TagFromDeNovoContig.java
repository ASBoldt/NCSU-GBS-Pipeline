/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.genome.GBS.deNovoContigs;

import net.maizegenetics.genome.GBS.*;
import java.nio.ByteBuffer;
import java.util.Arrays;
import net.maizegenetics.genome.BaseEncoder;

/**
 *
 * @author glaubitz
 */
public class TagFromDeNovoContig implements Comparable {

    long tag0, tag1;

    // novel contig info
    int contig;  // Shiran's actual name (number) for the contig (will need to use a numeric code for Joshua's)
    byte tagNum;  // the order of the tag in the contig (from 1 to a max of 80 in CSHL contigs)
    byte nTagsInCtg; // the total number of tags in the contig (max of 80 in CSHL contigs)
    byte inNovelCtg;  // 1 = true, 0 = false, -1 = unknown

    // genetic mapping info
    byte genChr; // the chr that the tag genetically maps to
    int startGenPos;
    int endGenPos;
    short refCnt;
    short altCnt;
    float bestR;  // best recombination rate

    // physical mapping info
    byte phyChr; // chr on B73_RefGen_v2 that the tag aligns to (with at least one mismatch); unknown is -1
    byte strand; //'+', '-', '?'  // strand on B73_RefGen_v2
    int positionMin, positionMax;  //B73_RefGen_v2 position for start and end of read (via alignment); unknown is Integer.MIN_VALUE
    short nextCutDistance;  //distance between known cut sites on B73_RefGen_v2, missing Short.MIN_VALUE
    byte divergence;  //number of diverging bp from reference, unknown Byte.MIN_VALUE
    byte multimaps;  //number of locations this readSet maps to, unknown Byte.MIN_VALUE
    byte phyConsenChr; // the consensus phys chr across the tags in the ctg that match the genetic consensus; unknown is -1
    int phyConsenPosition;  // the average positionMin across the tags in the ctg that match the genetic consensus & the phyConsenChr; unknown is Integer.MIN_VALUE

    static final byte unknownChr = -1;
    public static final int totalByteSize = 59;
    static final int 
            tag0Off = 0,
            tag1Off = tag0Off + 8,
            contigOff = tag1Off + 8,
            tagNumOff = contigOff + 4,
            nTagsInCtgOff = tagNumOff + 1,
            inNovelCtgOff = nTagsInCtgOff + 1,
            genChrOff = inNovelCtgOff + 1,
            startPosOff = genChrOff + 1,
            endPosOff = startPosOff + 4,
            refCntOff = endPosOff + 4,
            altCntOff = refCntOff + 2,
            bestROff = altCntOff + 2,
            chromosomeOff = bestROff + 4,
            strandOff = chromosomeOff + 1,
            positionMinOff = strandOff + 1,
            positionMaxOff = positionMinOff + 4,
            nextCutDistanceOff = positionMaxOff + 4,
            divergenceOff = nextCutDistanceOff + 2,
            multimapsOff = divergenceOff + 1,
            phyConsenChrOff = multimapsOff + 1,
            phyConsenPositionOff = phyConsenChrOff + 1;

    
    TagFromDeNovoContig() {
        tag0 = Long.MIN_VALUE;
        tag1 = Long.MIN_VALUE;
        contig = Integer.MIN_VALUE;
        tagNum = Byte.MIN_VALUE;
        nTagsInCtg = Byte.MIN_VALUE;
        inNovelCtg = -1;
        genChr = unknownChr;
        startGenPos = Integer.MIN_VALUE;
        endGenPos = Integer.MIN_VALUE;
        refCnt = Short.MIN_VALUE;
        altCnt = Short.MIN_VALUE;
        bestR = Float.NaN;
        phyChr = unknownChr;
        strand = (byte) '?';
        positionMin = Integer.MIN_VALUE;
        positionMax = Integer.MIN_VALUE;
        nextCutDistance = Short.MIN_VALUE;
        divergence = Byte.MIN_VALUE;
        multimaps = Byte.MIN_VALUE;
        phyConsenChr = unknownChr;
        phyConsenPosition = Integer.MIN_VALUE;
    }

    TagFromDeNovoContig(long[] tag, int novelContig, byte tagNum, byte nTagsInCtg, byte inNovelCtg, byte genChr,
            int startPos, int endPos, short refCnt, short altCnt, float bestR, byte chromosome,
            byte strand, int positionMin, int positionMax, short nextCutDistance, byte divergence,
            byte multimaps, byte phyConsenChr, int phyConsenPosition) {
        this.tag0 = tag[0];
        this.tag1 = tag[1];
        this.contig = novelContig;
        this.tagNum = tagNum;
        this.nTagsInCtg = nTagsInCtg;
        this.inNovelCtg = inNovelCtg;
        this.genChr = genChr;
        this.startGenPos = startPos;
        this.endGenPos = endPos;
        this.refCnt = refCnt;
        this.altCnt = altCnt;
        this.bestR = bestR;
        this.phyChr = chromosome;
        this.strand = strand;
        this.positionMin = positionMin;
        this.positionMax = positionMax;
        this.nextCutDistance = nextCutDistance;
        this.divergence = divergence;
        this.multimaps = multimaps;
        this.phyConsenChr = phyConsenChr;
        this.phyConsenPosition = phyConsenPosition;
    }

    TagFromDeNovoContig(ByteBuffer bb) {
        tag0 = bb.getLong(tag0Off);
        tag1 = bb.getLong(tag1Off);
        contig = bb.getInt(contigOff);
        tagNum = bb.get(tagNumOff);
        nTagsInCtg = bb.get(nTagsInCtgOff);
        inNovelCtg = bb.get(inNovelCtgOff);
        genChr = bb.get(genChrOff);
        startGenPos = bb.getInt(startPosOff);
        endGenPos = bb.getInt(endPosOff);
        refCnt = bb.getShort(refCntOff);
        altCnt = bb.getShort(altCntOff);
        bestR = bb.getFloat(bestROff);
        phyChr = bb.get(chromosomeOff);
        strand = bb.get(strandOff);
        positionMin = bb.getInt(positionMinOff);
        positionMax = bb.getInt(positionMaxOff);
        nextCutDistance = bb.getShort(nextCutDistanceOff);
        divergence = bb.get(divergenceOff);
        multimaps = bb.get(multimapsOff);
        phyConsenChr = bb.get(phyConsenChrOff);
        phyConsenPosition = bb.getInt(phyConsenPositionOff);
    }

    @Override
    public String toString() {
        String s = BaseEncoder.getSequenceFromLong(tag0) + BaseEncoder.getSequenceFromLong(tag1) + "\t"
                + contig + "\t" + tagNum + "\t" + nTagsInCtg + "\t" + inNovelCtg + "\t" + genChr + "\t" + startGenPos + "\t"
                + endGenPos + "\t" + refCnt + "\t" + altCnt + "\t" + bestR + "\t"
                + phyChr + "\t" + (char) strand + "\t" + positionMin + "\t"
                + positionMax + "\t" + nextCutDistance + "\t" + divergence + "\t" + multimaps + "\t" + phyConsenChr + "\t" + phyConsenPosition;
        return s;
    }

    public static String getHeaderString () {
        String s = "SequenceTag\tContigNum\ttagNum\tnTagsInCtg\tinNovelCtg\tBestGeneticChr\tBestGeneticStart\t"
                + "BestGeneticEnd\tRefAlleleCnt\tAltAlleleCnt\tbestR\tBLATChr\tStrand\tAGPv2Start\t"
                + "AGPv2End\tFragmentSize\tMismatchesFromAGPv2\tnBLATPositions\tphyConsenChr\tphyConsenPos";
        return s;
    }

    public void setPosition(TagFromDeNovoContig newPosition) {
        this.phyChr = newPosition.phyChr;
        this.positionMin = newPosition.positionMin;
        this.positionMax = newPosition.positionMax;
        this.strand = newPosition.strand;
    }

    public byte[] toByte() {
        byte[] theResult = new byte[totalByteSize];
        ByteBuffer fw = ByteBuffer.wrap(theResult);
        fw.putLong(tag0);
        fw.putLong(tag1);
        fw.putInt(contig);
        fw.put(tagNum);
        fw.put(nTagsInCtg);
        fw.put(inNovelCtg);
        fw.put(genChr);
        fw.putInt(startGenPos);
        fw.putInt(endGenPos);
        fw.putShort(refCnt);
        fw.putShort(altCnt);
        fw.putFloat(bestR);
        fw.put(phyChr);
        fw.put(strand);
        fw.putInt(positionMin);
        fw.putInt(positionMax);
        fw.putShort(nextCutDistance);
        fw.put(divergence);
        fw.put(multimaps);
        fw.put(phyConsenChr);
        fw.putInt(phyConsenPosition);
        return theResult;
    }

    @Override
    public boolean equals(Object obj) {
        final TagFromDeNovoContig other = (TagFromDeNovoContig) obj;
        if (this.tag0 != other.tag0) {
            return false;
        }
        if (this.tag1 != other.tag1) {
            return false;
        }
        if (this.contig != other.contig) {
            return false;
        }
        if (this.tagNum != other.tagNum) {
            return false;
        }
        return true;
    }

    @Override
    public int compareTo(Object o) {
        TagFromDeNovoContig other = (TagFromDeNovoContig) o;
        if (this.equals(other)) {
            return 0;
        }
        if (this.contig < other.contig) {
            return -1;
        } else if (this.contig > other.contig) {
            return 1;
        }
        if (this.tagNum < other.tagNum) {
            return -1;
        } else if (this.tagNum > other.tagNum) {
            return 1;
        }
        if (tag0 < other.tag0) {
            return -1;
        } else if (tag0 > other.tag0) {
            return 1;
        }
        if (tag1 < other.tag1) {
            return -1;
        } else if (tag1 > other.tag1) {
            return 1;
        }
        return 0;
    }

    public long[] getRead() {
        long[] result = {tag0, tag1};
        return result;
    }

    public short getRefCnt() {
        return refCnt;
    }

    public float getBestR() {
        return bestR;
    }

    public byte getPhyChr() {
        return phyChr;
    }

    public boolean isPhyChrKnown() {
        if (phyChr < 0) {
            return false;
        }
        return true;
    }

    public byte getDivergence() {
        return divergence;
    }

    public byte getMultimaps() {
        return multimaps;
    }

    // SETTERS
    public void setTagSeq(long[] tag) {
        this.tag0 = tag[0];
        this.tag1 = tag[1];
    }

    public void setNovelContig(int novelContig) {
        this.contig = novelContig;
    }

    public void setTagNum(int tagNum) {
        if (tagNum > Byte.MAX_VALUE) {
            this.tagNum = Byte.MAX_VALUE;
        } else {
            this.tagNum = (byte) tagNum;
        }
    }

    public void setNTagsInCtg(int nTagsInCtg) {
        if (nTagsInCtg > Byte.MAX_VALUE) {
            this.nTagsInCtg = Byte.MAX_VALUE;
        } else {
            this.nTagsInCtg = (byte) nTagsInCtg;
        }
    }

    public void setGenChr(int genChr) {
        if (genChr > Byte.MAX_VALUE) {
            this.genChr = Byte.MAX_VALUE;
        } else if (genChr < 0) {
            this.genChr = unknownChr;
        } else {
            this.genChr = (byte) genChr;
        }
    }

    public void setStartGenPos(int startPos) {
        this.startGenPos = startPos;
    }

    public void setEndGenPos(int endPos) {
        this.endGenPos = endPos;
    }

    public void setRefCnt(int refCnt) {
        if (refCnt > Short.MAX_VALUE) {
            this.refCnt = Short.MAX_VALUE;
        } else {
            this.refCnt = (short) refCnt;
        }
    }

    public void setAltCnt(int altCnt) {
        if (altCnt > Short.MAX_VALUE) {
            this.altCnt = Short.MAX_VALUE;
        } else {
            this.altCnt = (short) altCnt;
        }
    }

    public void setBestR(float bestR) {
        this.bestR = bestR;
    }

    public void setPhyChr(int phyChr) {
        if (phyChr > Byte.MAX_VALUE) {
            this.phyChr = Byte.MAX_VALUE;
        } else if (phyChr < 0) {
            this.phyChr = unknownChr;
        } else {
            this.phyChr = (byte) phyChr;
        }
    }

    public void setStrand(byte strand) {
        if (strand == (byte) '+') {
            this.strand = (byte) '+';
        } else if (strand == (byte) '-') {
            this.strand = (byte) '-';
        } else {
            this.strand = (byte) '?';
        }
    }

    public void setStrand(char strand) {
        if (strand == '+') {
            this.strand = (byte) '+';
        } else if (strand == '-') {
            this.strand = (byte) '-';
        } else {
            this.strand = (byte) '?';
        }
    }

    public void setPositionMin(int positionMin) {
        if (positionMin > 0) {
            this.positionMin = positionMin;
        } else {
            this.positionMin = Integer.MIN_VALUE;
        }
    }  
            
    public void setPositionMax(int positionMax) {
        if (positionMax > 0) {
            this.positionMax = positionMax;
        } else {
            this.positionMax = Integer.MIN_VALUE;
        }
    }

    public void setNextCutDistance(int nextCutDistance) {
        if (nextCutDistance > Short.MAX_VALUE) {
            this.nextCutDistance = Short.MAX_VALUE;
        } else {
            this.nextCutDistance = (short) nextCutDistance;
        }
    }

    public void setDivergence(int divergence) {
        if (divergence > Byte.MAX_VALUE) {
            this.divergence = Byte.MAX_VALUE;
        } else {
            this.divergence = (byte) divergence;
        }
    }

    public void setMultimaps(int multimaps) {
        if (multimaps > Byte.MAX_VALUE) {
            this.multimaps = Byte.MAX_VALUE;
        } else {
            this.multimaps = (byte) multimaps;
        }
    }

    public void setPhyConsenChr(int phyConsenChr) {
        if (phyConsenChr > Byte.MAX_VALUE) {
            this.phyConsenChr = Byte.MAX_VALUE;
        } else if (phyConsenChr < 0) {
            this.phyConsenChr = unknownChr;
        } else {
            this.phyConsenChr = (byte) phyConsenChr;
        }
    }

    public void setPhyConsenPosition(int phyConsenPosition) {
        if (phyConsenPosition > 0) {
            this.phyConsenPosition = phyConsenPosition;
        } else {
            this.phyConsenPosition = Integer.MIN_VALUE;
        }
    }
}
