/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS.plasmid;

import java.io.*;
import java.util.Arrays;
import net.maizegenetics.genome.GBS.ReadsByTaxaPlasmid;
import net.maizegenetics.genome.BaseEncoder;

/**
 *
 * @author fl262
 */
public class getReadsCounts {
	int allCount;
	int[] totalCounts;
	ReadsByTaxaPlasmid rbtMin;
	readAndCount[] rac;
	int cutoff = 500;
	public getReadsCounts(String totalCountsDirS, String infileS, String positionFileS, String outCountFileS, String outPosiFileS) {
		getTotalCounts(totalCountsDirS);
		getCountsByTaxa(infileS);
		getPosition(positionFileS);
		writeOutfile(outCountFileS, outPosiFileS);
	}

	public void writeOutfile(String outCountFileS, String outPosiFileS) {
		int[] sum = new int[rac.length];
		for (int i = 0; i < rac.length; i++) {
			sum[i] = rac[i].sumCount;
		}
		int hit = Arrays.binarySearch(sum, cutoff);
		while (sum[hit] > cutoff) {
			hit--;
		}
		hit++;
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outCountFileS), 65536);
			bw.write(String.valueOf(rbtMin.taxaNum)+"\t"+String.valueOf(rbtMin.haplotypeNum-hit));
			bw.newLine();
			bw.write("Sum\t");
			for (int i = 0; i < rbtMin.taxaNum-1; i++) {
				bw.write(rbtMin.taxaNames[i]+"\t");
			}
			bw.write(rbtMin.taxaNames[rbtMin.taxaNum-1]);
			bw.newLine();
			bw.write(String.valueOf(allCount)+"\t");
			for (int i = 0; i < rbtMin.taxaNum-1; i++) {
				bw.write(String.valueOf(totalCounts[i])+"\t");
			}

			bw.write(String.valueOf(totalCounts[rbtMin.taxaNum-1]));
			bw.newLine();
			for (int i = hit; i < rac.length; i++) {
				bw.write(BaseEncoder.getSequenceFromLong(rac[i].haplotype)+"\t");
				bw.write(String.valueOf(rac[i].sumCount)+"\t");
				for (int j = 0; j < rac[i].countsByTaxa.length - 1; j++) {
					bw.write(String.valueOf(rac[i].countsByTaxa[j])+"\t");
				}
				bw.write(String.valueOf(rac[i].countsByTaxa[rac[i].countsByTaxa.length - 1]));
				bw.newLine();
			}
			bw.flush();
		}
		catch (Exception e) {
			System.out.println("!Error occurred in " + outCountFileS);
		}
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outPosiFileS), 65536);
			for(int i = hit; i < rac.length; i++) {
				bw.write(rac[i].positionS);
				bw.newLine();
			}
			bw.flush();
		}
		catch (Exception e) {
			System.out.println("!Error occurred in " + outCountFileS);
		}
	}

	public void getPosition(String positionFileS) {
		try {
			DataInputStream dis = new DataInputStream(new BufferedInputStream (new FileInputStream(positionFileS),65536));
			for(int i=0; i < rac.length; i++) {
				StringBuilder sb = new StringBuilder();
				long[] haplotype = new long[2];
                haplotype[0] = dis.readLong();
                haplotype[1] = dis.readLong();
				sb.append(BaseEncoder.getSequenceFromLong(haplotype)).append("\t");
                byte chrB=dis.readByte();
				sb.append(chrB).append("\t");
                byte strand=dis.readByte();
				sb.append(strand).append("\t");
                int positionMin=dis.readInt();
				sb.append(positionMin).append("\t");
                int positionMax=dis.readInt();
				sb.append(positionMax).append("\t");
                short nextCutDistance=dis.readShort();
				sb.append(nextCutDistance).append("\t");
                byte divergence=dis.readByte();
				sb.append(divergence).append("\t");
                float dcoP=dis.readFloat();
				sb.append(dcoP).append("\t");
                float mapP=dis.readFloat();
				sb.append(mapP).append("\t");
                byte multimaps=dis.readByte();
				sb.append(multimaps).append("\t");
				String positionS = sb.toString();
				rac[i].getPosition(positionS);
            }
			dis.close();
		}
		catch (Exception e) {
			System.out.println("Reading or writing error: " + positionFileS);
		}
		Arrays.sort(rac);
	}
	public void getCountsByTaxa(String infileS) {
		rbtMin = new ReadsByTaxaPlasmid(infileS, true);
		rac = new readAndCount[rbtMin.haplotypeNum];
		for (int i = 0; i < rbtMin.haplotypeNum; i++) {
			long[] hap = new long[2];
			hap[0] = rbtMin.haplotype[0][i];
			hap[1] = rbtMin.haplotype[1][i];
			rac[i] = new readAndCount(hap, rbtMin.hapDist[i]);
			rac[i].getSumCount();
		}
	}

	public void getTotalCounts (String totalCountsDirS) {
		File[] parsedFiles = new File(totalCountsDirS).listFiles();
		int taxaNum = parsedFiles.length;
		totalCounts = new int[taxaNum];
		for (int i = 0; i < taxaNum; i++) {
			totalCounts[i] = (int)parsedFiles[i].length() / 20;
			allCount = allCount + totalCounts[i];
		}
	}
	
	public static void main (String[] args) {
		String totalCountsDirS = "E:/GBS/ParsedTaxaDir/";
		String infileS = "E:/GBS/ReadsByTaxaMin/ReadsByTaxaMin.txt";
		String positionFileS = "E:/GBS/CommonReadsOnMap/CommonReadsOnMap.txt";
		String outCountFileS = "E:/temp/read_counts_pipe.txt";
		String outPosiFileS = "E:/temp/read_position_pipe.txt";
		new getReadsCounts(totalCountsDirS, infileS, positionFileS, outCountFileS, outPosiFileS);
	}
}

class readAndCount implements Comparable <readAndCount> {
	long[] haplotype;
	int[] countsByTaxa;
	int sumCount;
	String positionS;
	readAndCount (long[] haplotype, int[] countsByTaxa) {
		this.haplotype = haplotype;
		this.countsByTaxa = countsByTaxa;
	}
	void getSumCount () {
		for (int i = 0; i < countsByTaxa.length; i++) {
			sumCount = sumCount + countsByTaxa[i];
		}
	}
	void getPosition (String positionS) {
		this.positionS = positionS;
	}
	public int compareTo(readAndCount o) {
		return sumCount - o.sumCount;
	}
}
