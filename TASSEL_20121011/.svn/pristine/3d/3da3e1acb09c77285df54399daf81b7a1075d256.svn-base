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
public class linePersent {
	ReadsByTaxaPlasmid rbt;
	int[] linePre;
	public linePersent (String infile, String outfile) {
		rbt = new ReadsByTaxaPlasmid(infile, true);
		process ();
		output (outfile);
	}
	public void process () {
		linePre = new int[rbt.taxaNum];
		for (int i = 0; i < rbt.taxaNum; i++) {
			for (int j = 0; j < rbt.haplotypeNum; j++) {
				if (rbt.hapDist[j][i] != 0) {
					linePre[i]++;
				}
			}
		}
	}
	public void output (String outfile) {
		try  {
			BufferedWriter bw = new BufferedWriter (new FileWriter(outfile), 65536);
			for (int i = 0; i < rbt.taxaNum; i++) {
				float ratio = (float)linePre[i]/rbt.haplotypeNum;
				bw.write(rbt.taxaNames[i]+"\t"+String.valueOf(ratio));
				bw.newLine();
			}
			bw.flush();
			bw.close();
		}
		catch (Exception e) {

		}
	}
}
