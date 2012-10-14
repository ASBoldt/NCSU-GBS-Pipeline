/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS.switchgrass;

import java.io.BufferedWriter;
import java.io.FileWriter;
import net.maizegenetics.genome.GBS.ReadsByTaxaPlasmid;

/**
 *
 * @author fl262
 */
public class mkTxtFile {
	ReadsByTaxaPlasmid inRbt;
	public mkTxtFile (String infile, String outfile) {
		inRbt = new ReadsByTaxaPlasmid(infile, true);
		try {
			BufferedWriter bw = new BufferedWriter (new FileWriter(outfile), 65536);
			for (int i = 0; i < inRbt.taxaNum; i++) {
				bw.write(inRbt.taxaNames[i]+"\t");
			}
			bw.newLine();
			for (int i = 0; i < inRbt.haplotypeNum; i++) {
				for (int j = 0; j < inRbt.taxaNum; j++) {
					bw.write(String.valueOf(inRbt.hapDist[i][j])+"\t");
				}
				bw.newLine();
			}
			bw.flush();
			bw.close();
		}
		catch (Exception e) {
			
		}
	}

}
