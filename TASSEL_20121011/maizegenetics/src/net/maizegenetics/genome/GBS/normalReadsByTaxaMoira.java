/*
 Reads number in ReadsByTaxaMin is devided by the total reads number in ParsedTaxaDir for each Taxa
 */
package net.maizegenetics.genome.GBS;

import java.io.*;
import java.util.HashMap;


/**
 *
 * @author fl262
 */
public class normalReadsByTaxaMoira extends ReadsByTaxaMoira {

	HashMap<String, Integer> taxaName2ReadNum = new HashMap();
	public float[][] normalHapDist;
	int[] readNum;

	public normalReadsByTaxaMoira(String parsedDir, String infile, String outfile, boolean binary) {
		super(infile);
		getTotalReads(parsedDir);
		writeNormalFile(outfile);
	}

	public normalReadsByTaxaMoira(String infile) {
		super();
		readNormalFile(infile);
	}

	public void readNormalFile(String infile) {
		try {
			DataInputStream dis =new DataInputStream(new BufferedInputStream(new FileInputStream(infile), 65536));
			taxaNum = dis.readInt();
			haplotypeNum = dis.readInt();
			taxaNames = new String[taxaNum];
			readNum = new int[taxaNum];
			haplotype = new long[2][haplotypeNum];
			normalHapDist = new float[haplotypeNum][taxaNum];
			for (int i = 0; i < taxaNum; i++) {
				taxaNames[i] = dis.readUTF();
			}
			for (int i = 0; i < taxaNum; i++) {
				readNum[i] = dis.readInt();
			}
			for (int i = 0; i < haplotypeNum; i++) {
				haplotype[0][i] = dis.readLong();
				haplotype[1][i] = dis.readLong();
				for (int j = 0; j < taxaNum; j++) {
					normalHapDist[i][j] = dis.readFloat();
				}
			}
			dis.close();
			System.out.println(infile+" is input");
		}
		catch (Exception e) {
			System.out.println("Catch in reading infile e = "+e);
		}
	}
	
	public void writeNormalFile(String outfile) {
		normalHapDist = new float[haplotypeNum][taxaNum];
		try {
			DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfile), 65536));

				dos.writeInt(taxaNum);
				dos.writeInt(haplotypeNum);
			for (int i = 0; i < taxaNum; i++) {
				dos.writeUTF(taxaNames[i]);
			}
			for (int i = 0; i < taxaNum; i++) {
				dos.writeInt(readNum[i]);
			}
			for (int i = 0; i < haplotypeNum; i++) {
				dos.writeLong(haplotype[0][i]);
				dos.writeLong(haplotype[1][i]);
				for (int j = 0; j < taxaNum; j++) {
					dos.writeFloat((float)newHapDist[i][j]/readNum[j]);
				}
			}
			dos.flush();
			dos.close();
			System.out.println("Haplotypes are written to: " + outfile.toString());
		} catch (Exception e) {
			System.out.println("Catch in writing output file e=" + e);
		}
	}

	public void getTotalReads(String parsedDir) {
		File[] parsedFiles = new File(parsedDir).listFiles();
		String[] taxaName = new String[taxaNum];
		readNum = new int[taxaNum];
		for (int i = 0; i < taxaNum; i++) {
			String temp = parsedFiles[i].toString().replace("\\", "/");
			temp = temp.replaceAll(".+/|_.+", "");
			taxaName[i] = temp;
		}
		for (int i = 0; i < taxaNum; i++) {
			readNum[i] = (int)parsedFiles[i].length() / 20;
		}
	}

	public static void main(String args[]) {
		String parsedDir = "E:/GBS/ParsedTaxaDir/";
		String infile = "E:/GBS/ReadsByTaxaMin/ReadsByTaxaMin.txt";
		String outfile = "E:/GBS/NormalReadsByTaxa/NormalReadsByTaxa.txt";
		new normalReadsByTaxaMoira(parsedDir, infile, outfile, true);
	}
}
