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
public class FileShowing {

	String parentDir = GBSPipelineFei.parentDir.replaceFirst("/$", "_Visual/");
	String[] childDir = GBSPipelineFei.childDir.clone();
	int totalLine = 3000;

	public FileShowing(boolean ifPipe) {
		if (ifPipe) {
			setupDirectory();
			for (int i = 0; i < childDir.length; i++) {
				String dir = childDir[i];
				String sourceDir = GBSPipelineFei.parentDir + childDir[i] + "/";
				String desDir = parentDir + childDir[i] + "/";
				if (dir.equals("Genome")) {
					getGenome(sourceDir, desDir);
				} else if (dir.equals("DigestedGenome")) {
					getDigestedGenome(sourceDir, desDir);
				} else if (dir.equals("Solexa")) {
					getSolexa(sourceDir, desDir);
				} else if (dir.equals("Key")) {
					getKey(sourceDir, desDir);
				} else if (dir.equals("ParsedTaxaDir")) {
					getParsedTaxaDir(sourceDir, desDir);
				} else if (dir.equals("CollapsedReads")) {
					getCollapsedReads(sourceDir, desDir);
				} else if (dir.equals("CombinedReads")) {
					getCombinedReads(sourceDir, desDir);
				} else if (dir.equals("ReadsWorthPhysicalMap")) {
					getReadsWorthPhysicalMap(sourceDir, desDir);
				} else if (dir.equals("ReadsByTaxa")) {
					getReadsByTaxa(sourceDir, desDir);
				} else if (dir.equals("ReadsByTaxaMin")) {
					getReadsByTaxaMin(sourceDir, desDir);
				} else if (dir.equals("CommonReadsOnMap")) {
					getCommonReadsOnMap(sourceDir, desDir);
				} else if (dir.equals("GeneticMap")) {
					getGeneticMap(sourceDir, desDir);
				} else if (dir.equals("CommonReadsWithGenetic")) {
					getCommonReadsWithGenetic(sourceDir, desDir);
				} else if (dir.equals("NormalReadsByTaxa")) {
					getNormalReadsByTaxa(sourceDir, desDir);
				} else {
					System.out.println("There must be missing methods for source directories");
					System.exit(1);
				}
			}
		}
	}

	private void setupDirectory() {
		File[] FileOfData = new File[childDir.length];
		for (int i = 0; i < childDir.length; i++) {
			String path = parentDir + childDir[i] + "/";
			FileOfData[i] = new File(path);
			FileOfData[i].mkdirs();
		}
	}

	private void simpleCopyPaste(String sourceFile, String desFile) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(sourceFile), 65536);
			BufferedWriter bw = new BufferedWriter(new FileWriter(desFile), 65536);
			String strLine;
			for (int i = 0; i < totalLine; i++) {
				strLine = br.readLine();
				if (strLine == null) {
					break;
				} else {
					bw.write(strLine);
					bw.newLine();
				}
			}
			bw.flush();
			bw.close();
			br.close();
		} catch (Exception e) {
			System.out.println("Reading or writing error: " + sourceFile + " / " + desFile);
		}
	}

	private void forReadsWPhysicalMap(String sourceFile, String desFile) {
		File fn = new File(sourceFile);
		int rows = (int) (fn.length() / ReadWithPositionInfo.totalByteSize);
		int total;
		if (rows < totalLine) {
			total = rows;
		} else {
			total = totalLine;
		}
		try {
			DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(sourceFile), 65536));
			PrintWriter pw = new PrintWriter(desFile);
			for (int i = 0; i < total; i++) {
				long[] haplotype = new long[2];
				haplotype[0] = dis.readLong();
				haplotype[1] = dis.readLong();
				String seq = BaseEncoder.getSequenceFromLong(haplotype);
				byte chrB = dis.readByte();
				//  chrB[1]=dis.readByte();
				byte strand = dis.readByte();
				int positionMin = dis.readInt();
				int positionMax = dis.readInt();
				short nextCutDistance = dis.readShort();
				byte divergence = dis.readByte();
				float dcoP = dis.readFloat();
				float mapP = dis.readFloat();
				byte multimaps = dis.readByte();
				pw.printf(seq + "\t" + chrB + "\t" + strand + "\t" + positionMin + "\t" + positionMax + "\t");
				pw.println(nextCutDistance + "\t" + divergence + "\t" + dcoP + "\t" + mapP + "\t" + multimaps);
			}
			pw.flush();
			pw.close();
			dis.close();
		} catch (Exception e) {
			System.out.println("Reading or writing error: " + sourceFile + " / " + desFile);
		}

	}

	private void forReadsCounts(String sourceFile, String desFile) {
		File fn = new File(sourceFile);
		int rows = (int) (fn.length() / 20);
		int total;
		if (rows < totalLine) {
			total = rows;
		} else {
			total = totalLine;
		}
		try {
			DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(sourceFile), 65536));
			PrintWriter pw = new PrintWriter(desFile);
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
			System.out.println("Reading or writing error: " + sourceFile + " / " + desFile);
		}
	}

	private void forReadsByTaxa(String sourceFile, String desFile) {
		try {
			DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(sourceFile), 65536));
			PrintWriter pw = new PrintWriter(desFile);
			int total;
			int taxaNum = dis.readInt();
			int haplotypeNum = dis.readInt();
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

			for (int i = 0; i < haplotypeNum; i++) {
				if (i < total) {
					long[] haplotype = new long[2];
					haplotype[0] = dis.readLong();
					haplotype[1] = dis.readLong();
					String seq = BaseEncoder.getSequenceFromLong(haplotype);
					pw.printf(seq + "\t");
					for (int j = 0; j < taxaNum - 1; j++) {
						pw.printf(dis.readInt() + "\t");
					}
					pw.println(dis.readInt());
				} else {
					break;
				}
			}
			pw.flush();
			pw.close();
			dis.close();
		} catch (Exception e) {
			System.out.println("Reading or writing error: " + sourceFile + " / " + desFile);
		}
	}
	private void forNormalReadsByTaxa(String sourceFile, String desFile) {
		try {
			DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(sourceFile), 65536));
			PrintWriter pw = new PrintWriter(desFile);
			int total;
			int taxaNum = dis.readInt();
			int haplotypeNum = dis.readInt();
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
			for (int i = 0; i < taxaNum - 1; i++) {
				pw.printf(dis.readInt() + "\t");
			}
			pw.println(dis.readInt());
			for (int i = 0; i < haplotypeNum; i++) {
				if (i < total) {
					long[] haplotype = new long[2];
					haplotype[0] = dis.readLong();
					haplotype[1] = dis.readLong();
					String seq = BaseEncoder.getSequenceFromLong(haplotype);
					pw.printf(seq + "\t");
					for (int j = 0; j < taxaNum - 1; j++) {
						pw.printf(dis.readFloat() + "\t");
					}
					pw.println(dis.readFloat());
				} else {
					break;
				}
			}
			pw.flush();
			pw.close();
			dis.close();
		} catch (Exception e) {
			System.out.println("Reading or writing error: " + sourceFile + " / " + desFile);
		}
	}

	private void getFilesInDir(String sourceDir, String desDir, int choice) {
		String[] files = new File(sourceDir).list();
		for (int i = 0; i < files.length; i++) {
			String sourceFile = sourceDir + files[i];
			String desFile = desDir + files[i];
			switch (choice) {
				case 0:
					simpleCopyPaste(sourceFile, desFile);
					break;
				case 1:
					forReadsWPhysicalMap(sourceFile, desFile);
					break;
				case 2:
					forReadsCounts(sourceFile, desFile);
					break;
				case 3:
					forReadsByTaxa(sourceFile, desFile);
					break;
				case 4:
					forNormalReadsByTaxa(sourceFile, desFile);
					break;
			}
		}
		System.gc();
	}

	public void getGenome(String sourceDir, String desDir) {
		getFilesInDir(sourceDir, desDir, 0);
	}

	public void getDigestedGenome(String sourceDir, String desDir) {
		getFilesInDir(sourceDir, desDir, 1);
	}

	public void getSolexa(String sourceDir, String desDir) {
		getFilesInDir(sourceDir, desDir, 0);
	}

	public void getKey(String sourceDir, String desDir) {
		getFilesInDir(sourceDir, desDir, 0);
	}

	public void getParsedTaxaDir(String sourceDir, String desDir) {
		getFilesInDir(sourceDir, desDir, 2);
	}

	public void getCollapsedReads(String sourceDir, String desDir) {
		getFilesInDir(sourceDir, desDir, 2);
	}

	public void getCombinedReads(String sourceDir, String desDir) {
		getFilesInDir(sourceDir, desDir, 2);
	}

	public void getReadsWorthPhysicalMap(String sourceDir, String desDir) {
		getFilesInDir(sourceDir, desDir, 2);
	}

	public void getReadsByTaxa(String sourceDir, String desDir) {
		getFilesInDir(sourceDir, desDir, 3);
	}

	public void getReadsByTaxaMin(String sourceDir, String desDir) {
		getFilesInDir(sourceDir, desDir, 3);
	}

	public void getCommonReadsOnMap(String sourceDir, String desDir) {
		getFilesInDir(sourceDir, desDir, 1);
	}

	public void getGeneticMap(String sourceDir, String desDir) {
		getFilesInDir(sourceDir, desDir, 0);
	}

	public void getCommonReadsWithGenetic(String sourceDir, String desDir) {
		getFilesInDir(sourceDir, desDir, 1);
	}

	public void getNormalReadsByTaxa(String sourceDir, String desDir) {
		getFilesInDir(sourceDir, desDir, 4);
	}

	public static void main(String args[]) {
		FileShowing pipeFileshowing = new FileShowing(true);
	}
}
