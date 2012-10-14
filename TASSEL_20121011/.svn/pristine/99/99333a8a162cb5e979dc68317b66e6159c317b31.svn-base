/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.gbs.pav;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;

/**
 *
 * @author Fei Lu
 */
public class HapMapUtils {
	BufferedReader br;
	BufferedWriter bw;
	String header;
	String[] headerArray;
	public HapMapUtils (String hapMapFileS, String outFileS) {
		this.creatStream(hapMapFileS, outFileS);
	}

	public void creatStream (String hapMapFileS, String outFileS) {
		try {
			br = new BufferedReader (new FileReader(hapMapFileS), 65536);
			bw = new BufferedWriter (new FileWriter(outFileS), 65536);
			header = br.readLine();
			headerArray = header.split("\\s+");
		}
		catch (Exception e) {
			System.out.println(e.toString());
			System.exit(1);
		}
	}

	public void checkSegregationInFamilyAndOutput (String[][] families, float minMAFinFamily) {
		int[][] indices = new int[families.length][];
		for (int i = 0; i < families.length; i++) {
			indices[i] = this.getIndicesOfSamples(families[i]);
			int cnt = 0;
			for (int j = 0; j < indices[i].length; j++) {
				if (indices[i][j] == -1) cnt++;
			}
			System.out.println("family " + i + " is indexed. There are " + families[i].length + " samples, " + cnt + " are not found in hapmap");
		}
		try {
			for (int i = 0; i < 11; i++) {
				bw.write(headerArray[i] + "\t");
			}
			for (int i = 0; i < indices.length; i++) {
				for (int j = 0; j < indices[i].length; j++) {
					if (indices[i][j] == -1) continue;
					bw.write(headerArray[indices[i][j]]+"\t");
				}
			}
			bw.newLine();
			this.writeRecord(indices, minMAFinFamily);
			bw.flush();
			bw.close();
			br.close();
		}
		catch (Exception e) {
			System.out.println(e.toString());
			System.exit(1);
		}
	}

	private void writeRecord (int[][] indices, float minMAFinFamily) {
		try {
			String temp;
			while ((temp = br.readLine()) != null) {
				StringBuilder sb = new StringBuilder();
				String[] tem = temp.split("\\s+");
				for (int i = 0; i < 11; i++) {
					sb.append(tem[i]).append("\t");
				}
				String[] allele = tem[1].split("/");
				boolean ifStay = true;
				boolean[] ifUnderMinMAF = new boolean[indices.length];
				for (int i = 0; i < indices.length; i++) {
					int[] cnt = new int[2];
					for (int j = 0; j < indices[i].length; j++) {
						if (indices[i][j] != -1) {
							for (int k = 0; k < allele.length; k++) {
								if (tem[indices[i][j]].startsWith(allele[k])) cnt[k]++;
							}
							sb.append(tem[indices[i][j]]).append("\t");
						}
					}
					float maf;
					if (cnt[0] < cnt[1]) maf = (float)cnt[0] / (float)(cnt[0] + cnt[1]);
					else maf = (float)cnt[1] / (float)(cnt[0] + cnt[1]);
					if (maf < minMAFinFamily) ifUnderMinMAF[i] = true;
					else ifUnderMinMAF[i] = false;
				}
				for (int i = 0; i < indices.length; i++) {
					ifStay = ifUnderMinMAF[i] & ifStay;
				}
				if (!ifStay) {
					bw.write(sb.toString());
					bw.newLine();
				}
			}
		
		}
		catch (Exception e) {
			System.out.println(e.toString());
		}
	}

	public int[] getIndicesOfSamples (String[] sampleNames) {
		int[] indices = new int[sampleNames.length];
		for (int i = 0; i < indices.length; i++) {
			indices[i] = this.getIndexOfSample(sampleNames[i]);
		}
		return indices;
	}

	public int getIndexOfSample (String sampleName) {
		for (int i = 11; i < headerArray.length; i++) {
			if (headerArray[i].startsWith(sampleName)) {
				return i;
			}
		}
		return -1;
	}
}
