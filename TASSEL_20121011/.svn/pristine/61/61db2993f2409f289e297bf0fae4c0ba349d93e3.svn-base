/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS.plasmid;

import java.io.*;
import java.util.Arrays;
import java.util.Comparator;

/**
 *
 * @author fl262
 */
public class tagByTaxaShort {
	String[] taxaNames;
	int[] totalCountsByTaxa;
	tagAndDist[] tagDs;
	int totalCount;
	int tagNum;
	int taxaNum;
	public tagByTaxaShort (String[] taxaNames, int[] totalCountsByTaxa, tagCountArray combinedTCA, short[][] countDist) {
		this.taxaNames = taxaNames;
		this.totalCountsByTaxa = totalCountsByTaxa;
		this.tagNum = combinedTCA.totalTag;
		this.taxaNum = taxaNames.length;
		this.totalCount = 0;
		for (int i = 0; i < taxaNum; i++) {
			totalCount += totalCountsByTaxa[i];
		}
		tagDs = new tagAndDist[tagNum];
		for (int i = 0; i < tagNum; i++) {
			tagDs[i] = new tagAndDist(combinedTCA.tagCs[i].tag, combinedTCA.tagCs[i].count, countDist[i]);
		}
	}
	public void writeTBTshortTXT (String outFileS) {
		writeTBTshortTXT (new File(outFileS));
	}
	public void writeTBTshortTXT (File outFile) {
		Arrays.sort(tagDs, new sortByCount());
		try {
			BufferedWriter br = new BufferedWriter (new FileWriter(outFile), 65536);
			StringBuilder sb = new StringBuilder();
			sb.append(String.valueOf(taxaNum)).append("\t").append(String.valueOf(tagNum));
			br.write(sb.toString());
			br.newLine();
			sb = new StringBuilder("Sum\t");
			for (int i = 0; i < taxaNum; i++) {
				String name = taxaNames[i];
				name = name.replaceFirst("_.+", "");
				sb.append(name).append("\t");
			}
			sb.deleteCharAt(sb.length()-1);
			br.write(sb.toString());
			br.newLine();
			sb = new StringBuilder(String.valueOf(totalCount)+"\t");
			for (int i = 0; i < taxaNum; i++) {
				sb.append(String.valueOf(totalCountsByTaxa[i])).append("\t");
			}
			sb.deleteCharAt(sb.length()-1);
			br.write(sb.toString());
			br.newLine();
			for (int i = 0; i < tagNum; i++) {
				sb = new StringBuilder(baseEncoder.getSeqFromByteArray(tagDs[i].tag));
				sb.append("\t").append(tagDs[i].count).append("\t");
				for (int j = 0; j < taxaNum; j++) {
					sb.append(String.valueOf(tagDs[i].dist[j])).append("\t");
				}
				sb.deleteCharAt(sb.length()-1);
				br.write(sb.toString());
				br.newLine();
			}
			br.flush();
			br.close();
		}
		catch (Exception e) {
			System.out.println("Error occurred while writing TBT " + e.toString());
		}
	}
	public void writeTBTshort (String outFileS) {
		writeTBTshort(new File(outFileS));
	}
	public void writeTBTshort (File outFile) {
		try {
			DataOutputStream dos = new DataOutputStream(new BufferedOutputStream (new FileOutputStream(outFile), 65536));
			dos.writeInt(taxaNum);
			dos.writeInt(tagNum);
			for (int i = 0; i < taxaNum; i++) {
				dos.writeUTF(taxaNames[i]);
			}
			for (int i = 0; i < taxaNum; i++) {
				dos.writeInt(totalCountsByTaxa[i]);
			}
			dos.writeInt(totalCount);
			for (int i = 0; i < tagNum; i++) {
				dos.write(tagDs[i].tag);
				dos.writeInt(tagDs[i].count);
				for (int j = 0; j < taxaNum; j++) {
					dos.writeShort(tagDs[i].dist[j]);
				}
			}
			dos.flush();
			dos.close();
		}
		catch (Exception e) {
			System.out.println("Error occurred while writing " + outFile.toString() + " " + e.toString());
		}
	}
	class tagAndDist implements Comparable <tagAndDist> {
		byte[] tag;
		int count;
		short[] dist;
		public tagAndDist (byte[] tag, int count, short[] dist) {
			this.tag = tag;
			this.count = count;
			this.dist = dist;
		}

		public int compareTo(tagAndDist o) {
			for (int i = 0; i < tag.length; i++) {
				if (tag[i] == o.tag[i]) {
					continue;
				}
				else {
					return tag[i] - o.tag[i];
				}
			}
			return 0;
		}
	}
	class sortByCount implements Comparator <tagAndDist> {
		public int compare(tagAndDist o1, tagAndDist o2) {
			return o1.count - o2.count;
		}
	}
}
