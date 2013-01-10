/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.gbs.homology;

import net.maizegenetics.genome.BaseEncoder;

/**
 * Container class for returning the results of parsed barcoded sequencing read.
 * 
 * @author edbuckler
 */
public class ReadBarcodeResult {
    public String unprocessedSequence= null;
    public String processedSequence=null;
    public String paddedSequence=null;
    byte length;

    long[] read;
    String taxonName;

    //TODO this instantiation should also include the orginal unprocessedSequence, processedSequence, and paddedSequence - the the object encode it
    public ReadBarcodeResult(long[] read, byte length, String taxon) {
        this.read = read;
        this.length = length;
        this.taxonName = taxon;
    }
    
    public ReadBarcodeResult(long[] read, String paddedSequence, byte length, String taxon) {
        this.read = read;
        this.paddedSequence = paddedSequence;
        this.length = length;
        this.taxonName = taxon;
    }

    public ReadBarcodeResult(String sequence){
        unprocessedSequence = sequence;
    }

    @Override
    public String toString() {
        return BaseEncoder.getSequenceFromLong(read)+":"+(int)length+":"+taxonName;
    }

    public byte getLength() {
        return length;
    }

    public long[] getRead() {
        return read;
    }

    public String getTaxonName() {
        return taxonName;
    }
}
