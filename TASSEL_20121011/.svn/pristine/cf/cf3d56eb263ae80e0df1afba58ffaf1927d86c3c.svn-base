/*
 read marker file of 644 or 55k format, set value of genotype for Glm analysis
 */
package net.maizegenetics.genome.GBS.repeat;

import java.io.*;
import java.util.ArrayList;
import java.util.Comparator;

/**
 *
 * @author fl262
 */
public class GeneticInfoArray {

	int taxaNum, markerNum;
	String title;
	String taxaName[];
	GeneticInfo[] markers;
	int genoValue4Tassel[][];
	byte genoValue[][];

	GeneticInfoArray(File sourceGenoFile) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(sourceGenoFile), 65536);
			String temp = br.readLine();
			String[] tempS = temp.split("\\t");
			StringBuilder sb = new StringBuilder();
			for (int i = 0; i < 11; i++) {
				sb.append(tempS[i]).append("\t");
			}
			title = sb.toString();
			taxaNum = tempS.length - 11;
			taxaName = new String[taxaNum];
			for (int i = 11; i < tempS.length; i++) {
				taxaName[i - 11] = tempS[i];
			}
			ArrayList<String> fileLine = new ArrayList<String>();
			temp = br.readLine();
			while (temp != null) {
				fileLine.add(temp);
				temp = br.readLine();
			}
			br.close();
			String[] lineArray = fileLine.toArray(new String[fileLine.size()]);
			markerNum = lineArray.length;
			markers = new GeneticInfo[markerNum];
			for (int i = 0; i < markerNum; i++) {
				markers[i] = new GeneticInfo(lineArray[i], taxaNum);
			}

		} catch (Exception e) {
			System.out.println(sourceGenoFile);
		}
	}

	public void writeFile (File desFile) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(desFile), 65536);
			bw.write(title);
			for (int i = 0; i < markerNum; i++) {
				bw.write(markers[i].markerName);
				bw.write("\t");
				for (int j = 0; j < markers[i].bgInfo.length; j++) {
					bw.write(markers[i].bgInfo[j]);
					bw.write("\t");
				}
				for (int j = 0; j < markerNum - 1; j++) {
					bw.write(markers[i].genotype[j]);
					bw.write("\t");
				}
				bw.write(markers[i].genotype[markerNum-1]);
				bw.newLine();
			}
			bw.flush();
			bw.close();
		}
		catch (Exception e) {
			System.out.println(desFile.toString());
		}
	}

	public void setGenoValue() {
		String[] stanGeno = new String[taxaNum];
		for (int i = 0; i < markerNum; i++) {
			boolean jump = true;
			for (int j = 0; j < taxaNum; j++) {
				if (stanGeno[j] == null) {
					if (markers[i].genotype[j].matches("NN")) {
						jump = false;
					}
					else if (markers[i].genotype[j].matches("MM")) {
						jump = false;
					}
					else {
						if (markers[i].genotype[j].matches("(\\w)\\1")) {
							stanGeno[j] = markers[i].genotype[j];
						} else {
							jump = false;
						}
					}
				}
			}
			if (jump == true) {
				break;
			}
		}
		genoValue = new byte[markerNum][taxaNum];
		for (int i = 0; i < markerNum; i++) {
			for (int j = 0; j < taxaNum; j++) {
				if (markers[i].genotype[j].matches("NN")) {
					genoValue[i][j] =  -1;
				}
				else if (markers[i].genotype[j].matches("MM")) {
					genoValue[i][j] =  -2;
				}
				else {
					if (markers[i].genotype[j].matches("(\\w)\\1")) {
						if (markers[i].genotype[j].equals(stanGeno[j])) {
							genoValue[i][j] = 0;
						}
						else {
							genoValue[i][j] = 2;
						}

					} else {
						genoValue[i][j] = 1;
					}
				}
			}
		}
	}

	public void setGenoValue4Tassel() {
		String[] stanGeno = new String[taxaNum];
		for (int i = 0; i < markerNum; i++) {
			boolean jump = true;
			for (int j = 0; j < taxaNum; j++) {
				if (stanGeno[j] == null) {
					if (markers[i].genotype[j].matches("NN")) {
						jump = false;
					} else {
						if (markers[i].genotype[j].matches("(\\w)\\1")) {
							stanGeno[j] = markers[i].genotype[j];
						} else {
							jump = false;
						}
					}
				}
			}
			if (jump == true) {
				break;
			}
		}
		genoValue4Tassel = new int[markerNum][taxaNum];
		for (int i = 0; i < markerNum; i++) {
			for (int j = 0; j < taxaNum; j++) {
				if (markers[i].genotype[j].matches("NN")) {
					genoValue4Tassel[i][j] = -999;
				} else {
					if (markers[i].genotype[j].matches("(\\w)\\1")) {
						if (markers[i].genotype[j].equals(stanGeno[j])) {
							genoValue4Tassel[i][j] = 0;
						}
						else {
							genoValue4Tassel[i][j] = 2;
						}

					} else {
						genoValue4Tassel[i][j] = 1;
					}
				}
			}
		}
	}
}

class GeneticInfo {
	String markerName;
	byte chrom;
	int position;
	String[] genotype;
	String[] bgInfo;

	GeneticInfo(String markerLine, int taxaNum) {
		String[] temp = markerLine.split("\\t");
		genotype = new String[taxaNum];
		markerName = temp[0];
		chrom = Byte.valueOf(temp[2]);
		position = Integer.valueOf(temp[3]);
		bgInfo = new String[11];
		for (int i = 0; i < 11; i++) {
			bgInfo[i] = temp[i];
		}
		for (int i = 11; i < temp.length; i++) {
			genotype[i - 11] = temp[i];
		}
	}
}

class sortByMarkerName implements Comparator <GeneticInfo> {
	public int compare (GeneticInfo gi1, GeneticInfo gi2) {
		return gi1.markerName.compareTo(gi2.markerName);
	}
}

class sortByChrom implements Comparator <GeneticInfo> {
	public int compare (GeneticInfo gi1, GeneticInfo gi2) {
		return gi1.chrom - gi2.chrom;
	}
}
