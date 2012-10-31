/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS.switchgrass;
import java.io.*;
import net.maizegenetics.genome.BaseEncoder;
import net.maizegenetics.genome.GBS.ReadsByTaxa;
/**
 *
 * @author fl262
 */
public class showTBTByte {
	ReadsByTaxa inTbt;
	int outLine = 200000;
	public showTBTByte (String infileS, String outfileS) {
		inTbt = new ReadsByTaxa(infileS, true);
		try {
			BufferedWriter bw = new BufferedWriter (new FileWriter(outfileS), 65536);
			for (int i = 0; i < inTbt.taxaNum; i++) {
				bw.write(inTbt.taxaNames[i]+"\t");
			}
			bw.newLine();
			for (int i = 100000; i < outLine; i++) {
				long[] hap = new long[2];
				hap[0] = inTbt.haplotype[0][i];
				hap[1] = inTbt.haplotype[1][i];
				bw.write(BaseEncoder.getSequenceFromLong(hap)+"\t");
				for (int j = 0; j < inTbt.taxaNum; j++) {
					bw.write(String.valueOf(inTbt.hapDist[i][j])+"\t");
				}
				bw.newLine();
			}
			bw.flush();
			bw.close();
		}
		catch (Exception e) {

		}
	}
	public static void main (String[] args) {
		String infileS = "M:/HapMapGenerator/TagsByTaxaByte/TagsByTaxaByte.bin";
		String outfileS = "M:/HapMapGenerator/TagsByTaxaByte/TagsByTaxaByte.txt";
		new showTBTByte (infileS, outfileS);
	}

}
