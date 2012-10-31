/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.genome.GBS;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import net.maizegenetics.genome.BaseEncoder;

/**
 *
 * @author jcg233
 * Written to obtain the information that John Doebley wants regarding the Teo (allele A)
 * and W22 (allele B) alleles to facilitate marker development for fine mapping
 *
 * The initial AlleleInfoFile was written in ConvertGBSToAlignment::writeAlleleInfo.
 * Here I just want to read in the B73 allele and then find its AGPv1 position
 *
 */
public class AlleleInfo {

    public AlleleInfo(String AlleleInfoFileName, String VDfileName, String OutFileName, int refAlleleField) {
        String refAllele;
        try {
            BufferedReader br = new BufferedReader(new FileReader(AlleleInfoFileName), 65536);
            ReadsWPhysicalMap VD = new ReadsWPhysicalMap(VDfileName,true);
            VD.sortTable(true);  // sort by sequence tag
            BufferedWriter bw = new BufferedWriter(new FileWriter(OutFileName), 65536);
            bw.write("B73_allele\tAGPv1Chr\tAGPv1StartPos\tAGPv1EndPos\tv1Strand\n");
            bw.flush();
            String inputLine = br.readLine();  // skip the header line
            while ( (inputLine = br.readLine()) != null) {
                String[] fields = inputLine.split("\t");
                refAllele = fields[refAlleleField];
                int[] phyHits = VD.getReadIndexSet(BaseEncoder.getLongArrayFromSeq(refAllele));
                if (phyHits == null) {
                    bw.write(refAllele + "\tNULL\tNULL\tNULL\tNULL\n");
                }
                else if (phyHits.length == 1) {
                    bw.write(refAllele + "\t" + VD.getReadWithPosition(phyHits[0]).chromosome
                            + "\t" + VD.getReadWithPosition(phyHits[0]).positionMin
                            + "\t" + VD.getReadWithPosition(phyHits[0]).positionMax
                            + "\t" + (char) VD.getReadWithPosition(phyHits[0]).strand + "\n");
                }
                else {
                    bw.write(refAllele + "\tmultiple\tmultiple\tmultiple\tmultiple\n");
                }
            }
            bw.close();
        } catch (Exception e) {
            System.out.println("Error writing AlleleInfo: " + e);
        }
    }
}
