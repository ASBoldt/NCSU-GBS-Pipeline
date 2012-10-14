/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS;

import java.io.*;
import net.maizegenetics.genome.BaseEncoder;
/**
 *
 * @author fl262
 */
public class getBinary4Ram {
	public getBinary4Ram (String infileStr, String outfileStr, int format) {
		File infile = new File (infileStr);
		File outfile = new File (outfileStr);
		if (format == 0) {
			forReadsCounts (infile, outfile);
		}
		else if (format == 1) {
			forReadsByTaxaInt (infile, outfile);
		}
		else if (format == 2) {
			forReadsByTaxaByte (infile, outfile);
		}
	}
	private void forReadsByTaxaInt(File inFile, File outFile) {
		try {
			DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(inFile), 65536));
			PrintWriter pw = new PrintWriter(outFile);
			int total;
			int taxaNum = dis.readInt();
			int haplotypeNum = dis.readInt();
			int totalLine = haplotypeNum;
			if (haplotypeNum < totalLine) {
				total = haplotypeNum;
			} else {
				total = totalLine;
			}
			pw.println(taxaNum + "\t" + haplotypeNum);
			for (int i = 0; i < taxaNum - 1; i++) {
				pw.printf(dis.readUTF() + "\t");
			}
			pw.println(dis.readUTF());

			for (int i = 0; i < total; i++) {
				long[] haplotype = new long[2];
				haplotype[0] = dis.readLong();
				haplotype[1] = dis.readLong();
				String seq = BaseEncoder.getSequenceFromLong(haplotype);
				pw.printf(seq + "\t");
				for (int j = 0; j < taxaNum - 1; j++) {
					pw.printf(dis.readInt() + "\t");
				}
				pw.println(dis.readInt());
			}
			pw.flush();
			pw.close();
			dis.close();
		} catch (Exception e) {
			System.out.println("Reading or writing error: " + inFile + " / " + outFile);
		}
	}
	private void forReadsByTaxaByte(File inFile, File outFile) {
		try {
			DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(inFile), 65536));
			PrintWriter pw = new PrintWriter(outFile);
			int total;
			int taxaNum = dis.readInt();
			int haplotypeNum = dis.readInt();
			int totalLine = haplotypeNum;
			if (haplotypeNum < totalLine) {
				total = haplotypeNum;
			} else {
				total = totalLine;
			}
			pw.println(taxaNum + "\t" + haplotypeNum);
			for (int i = 0; i < taxaNum - 1; i++) {
				pw.printf(dis.readUTF() + "\t");
			}
			pw.println(dis.readUTF());

			for (int i = 0; i < total; i++) {
				long[] haplotype = new long[2];
				haplotype[0] = dis.readLong();
				haplotype[1] = dis.readLong();
				String seq = BaseEncoder.getSequenceFromLong(haplotype);
				pw.printf(seq + "\t");
				for (int j = 0; j < taxaNum - 1; j++) {
					pw.printf(dis.readByte() + "\t");
				}
				pw.println(dis.readByte());
			}
			pw.flush();
			pw.close();
			dis.close();
		} catch (Exception e) {
			System.out.println("Reading or writing error: " + inFile + " / " + outFile);
		}
	}
	public void forReadsCounts(File inFile, File outFile) {
		int totalLine = 3000;
		int rows = (int) (inFile.length() / 20);
		int total;
		if (rows < totalLine) {
			total = rows;
		} else {
			total = rows;
		}
		try {
			DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(inFile), 65536));
			PrintWriter pw = new PrintWriter(outFile);
			for (int i = 0; i < total; i++) {
				long[] haplotype = new long[2];
				haplotype[0] = dis.readLong();
				haplotype[1] = dis.readLong();
				String seq = BaseEncoder.getSequenceFromLong(haplotype);
				int counting = dis.readInt();
				pw.println(seq + "\t" + counting);
			}
			pw.flush();
			pw.close();
			dis.close();
		} catch (Exception e) {
			System.out.println("Reading or writing error: " + inFile.toString() + " / " + outFile.toString());
		}
	}

	public static void main (String[] args) {
		String infileStr = "E:/GBS/ReadsByTaxa/ReadsByTaxa.txt";
		String outfileStr = "E:/out.txt";
		//0 is ReadCounts, 1 is ReadsByTaxaInt, 2 is ReadsByTaxaByte
		new getBinary4Ram (infileStr, outfileStr, 2);
	}
}
