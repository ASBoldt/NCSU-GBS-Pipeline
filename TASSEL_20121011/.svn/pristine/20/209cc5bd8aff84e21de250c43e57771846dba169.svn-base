/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.genome.GBS.deNovoContigs;

import net.maizegenetics.genome.GBS.*;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import net.maizegenetics.genome.BaseEncoder;

/**
 *
 * @author ed
 */
public class ContigHaplotypesByTaxa {

    public int[] contigIDs;
    public byte[][] hapDist;
    public int taxaNum = 0;
    public int contigCount = 0;
    public String[] taxaNames;

    public ContigHaplotypesByTaxa() {
    }

    public ContigHaplotypesByTaxa(String infile, boolean binary) {
        readDistFile(new File(infile), binary);
    }

    public ContigHaplotypesByTaxa(int[] contigIDs, byte[][] ctgHapDist, String[] namesForTaxa) {
        this.contigIDs = contigIDs;
        hapDist = ctgHapDist;
        taxaNames = namesForTaxa;
        taxaNum = namesForTaxa.length;
        contigCount = contigIDs.length;
    }

    public int getTaxaCount() {
        return taxaNames.length;
    }

    public String getTaxaName(int taxaIndex) {
        return taxaNames[taxaIndex];
    }

    public String[] getTaxaNames() {
        return taxaNames;
    }

    public int getIndexOfTaxaName(String taxon) {
        for (int i = 0; i < taxaNames.length; i++) {
            if (taxon.equals(taxaNames[i])) {
                return i;
            }
        }
        return -1;
    }

    public int getHaplotypeCountForTaxa(int haplotypeIndex, int taxaIndex) {
        return (int) hapDist[haplotypeIndex][taxaIndex];
    }

    public void setHaplotypeCountForTaxa(int hapIndex, int taxaIndex, int value) {
        if (value > Byte.MAX_VALUE) {
            hapDist[hapIndex][taxaIndex] = Byte.MAX_VALUE;
        } else if (value < 0) {
            hapDist[hapIndex][taxaIndex] = 0;
        } else {
            hapDist[hapIndex][taxaIndex] = (byte) value;
        }
    }

    public void addToHaplotypeCountForTaxa(int hapIndex, int taxaIndex, int addValue) {
        setHaplotypeCountForTaxa(hapIndex, taxaIndex, addValue + hapDist[hapIndex][taxaIndex]);
    }

    public byte[] getHaplotypeCountsForTaxa(int readIndex) {
        return hapDist[readIndex].clone();
    }

    public int getTaxaCountForHaplotype(int hapIndex) {   // how many taxa was a given contig haplotype seen in?
        int nTaxaWData = 0;
        for (int cnt : hapDist[hapIndex]) {
            if (cnt > 0) {
                ++nTaxaWData;
            }
        }
        return nTaxaWData;
    }

    public int getContigID(int i) {
        if (i >= contigCount) {
            return -1;
        }
        return contigIDs[i];
    }

    /**
     * Gets the first index of a contigID (the only one if a unique list).
     * If the read is not found then it return
     * a negative value indicating its insertion point.
     * @param read as a compressed long array
     * @return the index of the read in the array
     */
    public int getContigHaplotypeIndex(int contigID) {
        int hit = Arrays.binarySearch(contigIDs, contigID);
        if (hit < 1) {
            return hit;
        }
        while (contigIDs[hit - 1] == contigID) {
            hit--;
        }
        return hit;
    }

    void readDistFile(File inFile, boolean binary) {
        System.out.println("Reading Contig Haplotypes by Taxa distribution from:" + inFile.toString());
        int hapsInput = 0;
        if (binary) {
            try {
                DataInputStream rw = new DataInputStream(new BufferedInputStream(new FileInputStream(inFile), 4000000));
                taxaNum = rw.readInt();
                contigCount = rw.readInt();
                taxaNames = new String[taxaNum];
                contigIDs = new int[contigCount];
                hapDist = new byte[contigCount][taxaNum];
                for (int t = 0; t < taxaNum; t++) {
                    taxaNames[t] = rw.readUTF();
                }
                for (int i = 0; i < contigCount; i++) {
                    contigIDs[i] = rw.readInt();
                    for (int t = 0; t < taxaNum; t++) {
                        hapDist[i][t] = rw.readByte();
                    }
                    hapsInput++;
                }
                rw.close();
            } catch (Exception e) {
                System.out.println("Catch in writing output file e=" + e);
            }
        }
        else {
            ArrayList<String> inputLine;
            try {
                BufferedReader br = new BufferedReader(new FileReader(inFile), 65536);
                inputLine = new ArrayList<String>( Arrays.asList(br.readLine().split("\t")) );
                taxaNum = Integer.parseInt(inputLine.get(0));
                contigCount = Integer.parseInt(inputLine.get(1));
                taxaNames = new String[taxaNum];
                contigIDs = new int[contigCount];
                hapDist = new byte[contigCount][taxaNum];
                inputLine = new ArrayList<String>( Arrays.asList(br.readLine().split("\t")) );
                for (int t = 0; t < taxaNum; t++) {
                    taxaNames[t] = inputLine.get(t+1);  // blank cell before the taxa list
                }
                for (int i = 0; i < contigCount; i++) {
                    inputLine = new ArrayList<String>( Arrays.asList(br.readLine().split("\t")) );
                    contigIDs[i] = Integer.parseInt(inputLine.get(0));
                    for (int t = 0; t < taxaNum; t++) {
                        hapDist[i][t] = Byte.valueOf(inputLine.get(t+1));
                    }
                    hapsInput++;
                }
            } catch (Exception e) {
                System.out.println("Catch in writing output file e=" + e);
            }
        }
        System.out.println("Number of Taxa in file:" + taxaNum);
        System.out.println("Number of Haplotypes in file:" + hapsInput);
    }

    public void writeDistFile(File outFile, boolean binary) {
        int hapsOutput = 0;
        try {
            DataOutputStream fw = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outFile), 65536));
            if (binary) {
                fw.writeInt(taxaNum);
                fw.writeInt(contigCount);
                for (int t = 0; t < taxaNum; t++) {
                    fw.writeUTF(taxaNames[t]);
                }
            } else {
                fw.writeBytes(taxaNum + "\t" + contigCount + "\n");
                for (int t = 0; t < taxaNum; t++) {
                    fw.writeBytes("\t" + taxaNames[t]);
                }
                fw.writeBytes("\n");
            }
            for (int i = 0; i < contigIDs.length; i++) {
                if (!binary) {
                    fw.writeBytes(contigIDs[i] + "");
                    for (int t = 0; t < taxaNum; t++) {
                        fw.writeBytes("\t" + hapDist[i][t]);
                    }
                    fw.writeBytes("\n");
                } else {
                    fw.writeInt(contigIDs[i]);
                    for (int t = 0; t < taxaNum; t++) {
                        fw.writeByte(hapDist[i][t]);
                    }
                }
                hapsOutput++;
            }
            fw.flush();
            fw.close();
            System.out.println("Contig Haplotypes by Taxa distribution written to:" + outFile.toString());
            System.out.println("Number of contig haplotypes in file:" + hapsOutput);
        } catch (Exception e) {
            System.out.println("Catch in writing output file e=" + e);
        }
    }

    public int getHaplotypeTotalCount(int index) {
        int sum = 0;
        for (int cnt : hapDist[index]) {
            sum += cnt;
        }
        return sum;
    }

    public int getTotalNumCtgs() {
        return contigIDs.length;
    }
}
