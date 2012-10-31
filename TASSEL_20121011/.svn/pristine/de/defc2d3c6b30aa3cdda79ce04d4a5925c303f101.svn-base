/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS.switchgrass;

import java.io.*;
import net.maizegenetics.genome.GBS.ReadsByTaxaPlasmid;

/**
 *
 * @author fl262
 */
public class persent {
	public static void main (String[] args) {
		String infile = "M:/GBS/ReadsByTaxaMin/ReadsByTaxaMin.txt";
		String outfile = "M:/switch/Persent/persent.txt";
		ReadsByTaxaPlasmid rbt = new ReadsByTaxaPlasmid(infile, true);
		float[] per = new float[rbt.haplotypeNum];
		for (int i = 0; i < rbt.haplotypeNum; i++) {
			int count = 0;
			for (int j = 0; j < rbt.taxaNum; j++) {
				if (rbt.hapDist[i][j] != 0) {
					count++;
				}
			}
			per[i] = (float)count/rbt.taxaNum;

		}
		System.out.println(rbt.haplotypeNum);
		try {
			BufferedWriter bw = new BufferedWriter (new FileWriter(outfile), 65536*16);
			for (int i = 0; i < rbt.haplotypeNum; i++) {
				bw.write(String.valueOf(per[i]));
				bw.newLine();
			}
			bw.flush();
			bw.close();
		}
		catch (Exception e) {

		}
	}
}
