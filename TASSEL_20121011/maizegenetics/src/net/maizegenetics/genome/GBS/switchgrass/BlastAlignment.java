/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS.switchgrass;

import java.io.BufferedWriter;
import java.io.FileWriter;
import net.maizegenetics.genome.BaseEncoder;

/**
 *
 * @author fl262
 */
public class BlastAlignment {
	String database;
	public BlastAlignment (String ReadsByTaxaFileS, String outfileS) {
		this(new ReadsByTaxa(ReadsByTaxaFileS,true), outfileS);
	}
	public BlastAlignment (ReadsByTaxa rbt, String outfileS) {
		convert2database (rbt, outfileS);
		database = outfileS;
	}
	public BlastAlignment (String databaseFileS) {
		database = databaseFileS;
	}
	public void convert2database (ReadsByTaxa rbt, String outfileS) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outfileS), 65536);
			for (int i = 0; i < rbt.haplotypeNum; i++) {
				bw.write(">"+String.valueOf(i));
				bw.newLine();
				long[] temp = new long[2];
				temp[0] = rbt.haplotype[0][i];
				temp[1] = rbt.haplotype[1][i];
				bw.write(BaseEncoder.getSequenceFromLong(temp));
				bw.newLine();
			}
			bw.flush();
			bw.close();
		}
		catch (Exception e) {
			System.out.println("Error while writing "+ outfileS);
		}
	}
	public void formatDatabase () {
		try {
			Runtime run = Runtime.getRuntime();
			String cmd = "Makeblastdb -in " + database + " -dbtype nucl";
			Process p = run.exec(cmd);
			p.waitFor();
		}
		catch (Exception e) {
			System.out.println("Error while formatting database " + database);
		}
	}
	public void blastn (String blastFileS) {
		try {
			Runtime run = Runtime.getRuntime();
			String cmd = "blastn -query " + database + " -db " + database + " -out " + blastFileS +
						" -evalue 1e-20 -gapopen 5 -gapextend 2 -strand plus";
			System.out.println(cmd);
			Process p = run.exec(cmd);
			p.waitFor();
		}
		catch (Exception e) {
			System.out.println("Error while do alignment " + blastFileS);
		}
	}
}
