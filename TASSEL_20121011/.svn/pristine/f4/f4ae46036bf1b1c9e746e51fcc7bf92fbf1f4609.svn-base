/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS.switchgrass;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import net.maizegenetics.genome.BaseEncoder;

/**
 *
 * @author Fei Lu
 */
public class NoRefHapUtils {
	HapRecord[] hr;
	String[] taxaNames;
	int taxaNum;
	int hapRecordNum;
	int[] taxaIndex;

	public NoRefHapUtils (ReadsByTaxa rbt, String infileS) {
		TagPair tp = new TagPair(infileS);
		this.taxaNum = rbt.taxaNum;
		this.hapRecordNum = rbt.haplotypeNum/2;
		this.taxaNames = rbt.taxaNames;
		this.hr = new HapRecord[hapRecordNum];
		taxaIndex = new int[taxaNum];
		for (int i = 0; i < taxaIndex.length; i++) {
			taxaIndex[i] = i;
		}
		this.getPairedTBT(rbt, tp);
		this.checkFile(rbt);
		this.getHapRecord(rbt);
	}

	public NoRefHapUtils (String infileS) {
		this.readMapInfo(infileS);
		this.getTaxaIndex();
	}

	public void checkFile(ReadsByTaxa rbt) {
		int count = 0;
		int n = 1000;
		if (rbt.haplotypeNum < n) n = rbt.haplotypeNum;
		for (int i = 0; i < n; i += 2) {
			long[] query = new long[2];
			query[0] = rbt.haplotype[0][i];
			query[1] = rbt.haplotype[1][i];
			long[] hit = new long[2];
			hit[0] = rbt.haplotype[0][i+1];
			hit[1] = rbt.haplotype[1][i+1];
			count += BaseEncoder.seqDifferences(query[0], hit[0]);
			count += BaseEncoder.seqDifferences(query[1], hit[1]);
		}
		if (count == n/2) {
			System.out.println("Tags in TBTByte file are paired");
		}
		else {
			System.out.println("Tags in TBTByte file are not paired");
			System.exit(0);
		}
	}

	public int getOutputTaxaNumber () {
		return this.taxaIndex.length;
	}

	public int getHapRecordNumber () {
		return this.hr.length;
	}

	public void getPairedTBT (ReadsByTaxa rbt, TagPair tp) {
		if (rbt.haplotypeNum != tp.haplotype[0].length) {
			System.out.println("TBT file and TagPair file don't match");
			System.exit(0);
		}
		byte[][] newdist = new byte[rbt.haplotypeNum][];
		long[][] newHaplotype = new long[2][rbt.haplotypeNum];
		for (int i = 0; i < rbt.haplotypeNum; i++) {
			newdist[tp.order[i]] = rbt.hapDist[i];
			newHaplotype[0][tp.order[i]] = rbt.haplotype[0][i];
			newHaplotype[1][tp.order[i]] = rbt.haplotype[1][i];
		}
		rbt.hapDist = newdist;
		rbt.haplotype = newHaplotype;
		newdist = null;
		newHaplotype = null;
	}

	public void getTaxaIndex () {
		taxaIndex = new int[taxaNum];
		for (int i = 0; i < taxaIndex.length; i++) {
			taxaIndex[i] = i;
		}
	}

	public void getTaxaIndex (String[] selectedTaxa) {
		ArrayList<Integer> indexList = new ArrayList();
		Arrays.sort(selectedTaxa);
		for (int i = 0; i < taxaNames.length; i++) {
			int hit = Arrays.binarySearch(selectedTaxa, taxaNames[i]);
			if (hit > -1) indexList.add(i);
		}
		Integer[] tIndex =  indexList.toArray(new Integer[indexList.size()]);
		this.taxaIndex = new int[tIndex.length];
		for (int i = 0; i < tIndex.length; i++) {
			taxaIndex[i] = tIndex[i];
		}
	}

	public void getHapRecord(ReadsByTaxa rbt) {
		for (int i = 0; i < rbt.haplotypeNum; i += 2) {
			long[] query = new long[2];
			query[0] = rbt.haplotype[0][i];
			query[1] = rbt.haplotype[1][i];
			long[] hit = new long[2];
			hit[0] = rbt.haplotype[0][i+1];
			hit[1] = rbt.haplotype[1][i+1];
			hr[i/2] = new HapRecord(query, hit, rbt.hapDist[i], rbt.hapDist[i+1]);
		}
	}

	public void selectTaxa (String taxaFileS) {
		ArrayList<String> taxaList = new ArrayList();
		try {
			BufferedReader br = new BufferedReader (new FileReader(taxaFileS), 65536);
			String temp;
			while ((temp = br.readLine()) != null) {
				String[] tem = temp.split("\t");
				for (int i = 0; i < tem.length; i++) {
					taxaList.add(tem[i]);
				}
			}
		}
		catch (Exception e) {
			System.out.println("Error occured while reading taxaSelection file " + taxaFileS + " " + e.toString());
		}
		if (taxaList.size() == 0) {
			System.out.println("No taxa is selected for hapmap output");
			System.exit(0);
		}
		this.getTaxaIndex(taxaList.toArray(new String[taxaList.size()]));
	}

	public void readMapInfo (String infileS) {
		try {
			DataInputStream dis =new DataInputStream(new BufferedInputStream(new FileInputStream(infileS),65536));
			this.taxaNum = dis.readInt();
			this.hapRecordNum = dis.readInt();
			this.taxaNames = new String[this.taxaNum];
			for (int i = 0; i <  this.taxaNum; i++) {
				this.taxaNames[i] = dis.readUTF();
			}
			hr = new HapRecord[this.hapRecordNum];
			for (int i = 0; i < this.hapRecordNum; i++) {
				hr[i] = new HapRecord(2, taxaNum);
				hr[i].readRecord(dis);
			}
			System.out.println("MapInfo is read");
			System.out.println(hr.length + " TagPair in total");
		}
		catch (Exception e) {
			System.out.println("Error occurred while reading MapInfo file " + infileS + " " + e.toString());
		}
	}
	
	public void writeMapInfo (String outfileS) {
		try {
			DataOutputStream dos =new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfileS),65536));
			dos.writeInt(this.taxaNum);
			dos.writeInt(this.hapRecordNum);
			for (int i = 0; i < this.taxaNum; i++) {
				dos.writeUTF(this.taxaNames[i]);
			}
			for (int i = 0; i < this.hapRecordNum; i++) {
				hr[i].writeRecord(dos);
			}
			dos.flush();
			dos.close();
		}
		catch (Exception e) {
			System.out.println("Error occurred while writing MapInfo file " + outfileS + " " + e.toString());
		}
	}

	public void writeHapMap (String outfileS, double lMAF, double hMAF, double lCoverage) {
		this.writeHapMap(outfileS, lMAF, hMAF, lCoverage, (double)1);
	}

	public void writeHapMap (String outfileS, double lMAF, double hMAF, double lCoverage, double hCoverage) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outfileS), 65536);
			bw.write("rs	alleles	chrom	pos	strand	assembly	center	protLSID	assayLSID	panelLSID	QCcode	");
			for (int i = 0; i < this.taxaIndex.length - 1; i++) {
				bw.write(this.taxaNames[this.taxaIndex[i]] + "\t");
			}
			bw.write(this.taxaNames[this.taxaIndex.length-1]);
			bw.newLine();
			int total = 0;
			for (int i = 0; i < this.hapRecordNum; i++) {
				if (hr[i].writeHapMap(bw, taxaIndex, lMAF, hMAF, lCoverage, hCoverage, i+1)) total++;
			}
			bw.flush();
			bw.close();
			System.out.println("HapMap file is written. There are " + total + " SNPs in total");
		}
		catch (Exception e) {
			System.out.println("Error occurred while writing HapMap file " + outfileS + " " + e.toString());
		}
	}

	public void writeHapCount (String outfileS, double lMAF, double hMAF, double lCoverage) {
		this.writeHapCount(outfileS, lMAF, hMAF, lCoverage, (double)1);
	}

	public void writeHapCount (String outfileS, double lMAF, double hMAF, double lCoverage, double hCoverage) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outfileS), 65536);
			bw.write("rs	");
			for (int i = 0; i < this.taxaIndex.length; i++) {
				bw.write(this.taxaNames[this.taxaIndex[i]] + "\t");
			}
			bw.write("Count_allele1\tCount_allele2\tFrequency");
			bw.newLine();
			for (int i = 0; i < this.hapRecordNum; i++) {
				hr[i].writeHapCount(bw, taxaIndex, lMAF, hMAF, lCoverage, hCoverage, i+1);
			}
			bw.flush();
			bw.close();
			System.out.println("HapMapCount file is written");
		}
		catch (Exception e) {
			System.out.println("Error occurred while writing HapMap file " + outfileS + " " + e.toString());
		}
	}

	public void writeHapSeq (String outfileS, double lMAF, double hMAF, double lCoverage) {
		this.writeHapSeq(outfileS, lMAF, hMAF, lCoverage, (double)1);
	}

	public void writeHapSeq (String outfileS, double lMAF, double hMAF, double lCoverage, double hCoverage) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outfileS), 65536);
			for (int i = 0; i < this.hapRecordNum; i++) {
				hr[i].writeHapSeq(bw, taxaIndex, lMAF, hMAF, lCoverage, hCoverage, i+1);
			}
			bw.flush();
			bw.close();
			System.out.println("HapMapSeq file is written");
		}
		catch (Exception e) {
			System.out.println("Error occurred while writing HapMap file " + outfileS + " " + e.toString());
		}
	}

	public void writeHapAll (String outDirS, double lMAF, double hMAF, double lCoverage) {
		this.writeHapAll(outDirS, lMAF, hMAF, lCoverage, 1);
	}

	public void writeHapAll (String outDirS, double lMAF, double hMAF, double lCoverage, double hCoverage) {
		File hapMapF = new File(outDirS, "HapMap.hmp.txt");
		File hapCountF = new File(outDirS, "HapMap.hmc.txt");
		File hapSeqF = new File(outDirS, "HapMap.fas.txt");
		this.writeHapMap(hapMapF.toString(), lMAF, hMAF, lCoverage, hCoverage);
		this.writeHapCount(hapCountF.toString(), lMAF, hMAF, lCoverage, hCoverage);
		this.writeHapSeq(hapSeqF.toString(), lMAF, hMAF, lCoverage, hCoverage);
	}

	public void outputFilter (double lMAF, double hMAF, double lCoverage) {
		this.outputFilter(lMAF, hMAF, lCoverage, 1);
	}

	public void outputFilter (double lMAF, double hMAF, double lCoverage, double hCoverage) {
		ArrayList<HapRecord> hapRList = new ArrayList();
		for (int i = 0; i < this.hr.length; i++) {
			if (hr[i].checkIfOutput(taxaIndex, lMAF, hMAF, lCoverage, hCoverage)) {
				hapRList.add(hr[i]);
			}
		}
		HapRecord[] hapRArray = hapRList.toArray(new HapRecord[hapRList.size()]);
		hr = hapRArray;
		this.hapRecordNum = hr.length;
	}

	public void pseudoTestCross (String homoParent, String heteroParent) {
		ReadCounts homoRc = new ReadCounts (homoParent, true);
		ReadCounts heteroRc = new ReadCounts (heteroParent, true);
		ArrayList<HapRecord> hapRList = new ArrayList();
		for (int i = 0; i < this.hapRecordNum; i++) {
			int hit1 = homoRc.binarySearchOfHaps(this.hr[i].query);
			int hit2 = homoRc.binarySearchOfHaps(this.hr[i].hit);
			int hit3 = heteroRc.binarySearchOfHaps(this.hr[i].query);
			int hit4 = heteroRc.binarySearchOfHaps(this.hr[i].hit);
			if (hit1 >= 0 && hit2 < 0) {
				if (hit3 >= 0 && hit4 >= 0) {
					if (heteroRc.hapcount[hit3] < 5 || heteroRc.hapcount[hit4] < 5) continue;
					hapRList.add(hr[i]);
				}
			}
			else if (hit1 < 0 && hit2 >= 0) {
				if (hit3 >= 0 && hit4 >= 0) {
					if (heteroRc.hapcount[hit3] < 5 || heteroRc.hapcount[hit4] < 5) continue;
					hapRList.add(hr[i]);
				}
			}
		}
		HapRecord[] hapRArray = hapRList.toArray(new HapRecord[hapRList.size()]);
		hr = hapRArray;
		this.hapRecordNum = hr.length;
	}


	class HapRecord {
		long[] query;
		long[] hit;
		byte[] queryCounts;
		byte[] hitCounts;
		char snpQ;
		char snpH;
		char hetCode;

		HapRecord (long[] query, long[] hit, byte[] a, byte[] b) {
			this.query = query;
			this.hit = hit;
			this.queryCounts = a;
			this.hitCounts = b;
			this.getSnps();
		}

		HapRecord (int nlong, int taxaNum) {
			this.query = new long[nlong];
			this.hit = new long[nlong];
			this.queryCounts = new byte[taxaNum];
			this.hitCounts = new byte[taxaNum];
		}

		void getSnps () {
			for (int i = 0; i < query.length; i++) {
				int snpIndex = this.getSnpIndex(query[i], hit[i]);
				if (snpIndex != -1) {
					String q = BaseEncoder.getSequenceFromLong(query[i]);
					snpQ = q.charAt(snpIndex);
					String h = BaseEncoder.getSequenceFromLong(hit[i]);
					snpH = h.charAt(snpIndex);
					this.getHetCode();
					return;
				}
			}
		}

		void getHetCode () {
			if (snpQ + snpH == 132) hetCode = 'M';
			else if (snpQ + snpH == 136) hetCode = 'R';
			else if (snpQ + snpH == 149) hetCode = 'W';
			else if (snpQ + snpH == 138) hetCode = 'S';
			else if (snpQ + snpH == 151) hetCode = 'Y';
			else if (snpQ + snpH == 155) hetCode = 'K';
		}

		int getSnpIndex(long seq1, long seq2) {
			long mask=3;
				byte cnt=0;
			    long diff=seq1 ^ seq2;
			    for(int x=0; x < 32; x++) {
					if((diff&mask)>0) return 31 - x;
				    diff=diff>>2;
				}
			return -1;
		}

		void readRecord (DataInputStream dis) {
			try {
				for (int i = 0; i < query.length; i++) {
					query[i] = dis.readLong();
				}
				snpQ = dis.readChar();
				for (int i = 0; i < queryCounts.length; i++) {
					queryCounts[i] = dis.readByte();
				}
				for (int i = 0; i < hit.length; i++) {
					hit[i] = dis.readLong();
				}
				snpH = dis.readChar();
				for (int i = 0; i < hitCounts.length; i++) {
					hitCounts[i] = dis.readByte();
				}
				hetCode = dis.readChar();
			}
			catch (Exception e) {
				System.out.println(e.toString());
			}
		}

		void writeRecord (DataOutputStream dos) {
			try {
				for (int i = 0; i < query.length; i++) {
					dos.writeLong(query[i]);
				}
				dos.writeChar(snpQ);
				for (int i = 0; i < queryCounts.length; i++) {
					dos.writeByte(queryCounts[i]);
				}
				for (int i = 0; i < hit.length; i++) {
					dos.writeLong(hit[i]);
				}
				dos.writeChar(snpH);
				for (int i = 0; i < hitCounts.length; i++) {
					dos.writeByte(hitCounts[i]);
				}
				dos.writeChar(this.hetCode);
			}
			catch (Exception e) {
				System.out.println(e.toString());
			}
		}

		boolean writeHapMap (BufferedWriter bw, int[] taxaIndex, double lMAF, double hMAF, double lCoverage, double hCoverage, int ID) {
			if (this.checkIfOutput(taxaIndex, lMAF, hMAF, lCoverage, hCoverage)) {
				StringBuilder sb = new StringBuilder();
				sb.append("TP").append(ID).append("\t").append(snpQ).append("/").append(snpH).append("\t").append(0).append("\t").append(ID).append("\t+\tNA\tSWGDiv\tGBS\tNoRef\tSWGDiv\tQC+\t");
				for (int i = 0; i < taxaIndex.length; i++) {
					if (this.queryCounts[taxaIndex[i]] == 0) {
						if (this.hitCounts[taxaIndex[i]] == 0) {
							sb.append("N").append("\t");
						}
						else {
							sb.append(this.snpH).append("\t");
						}
					}
					else {
						if (this.hitCounts[taxaIndex[i]] == 0) {
							sb.append(this.snpQ).append("\t");
						}
						else {
							sb.append(this.hetCode).append("\t");
						}
					}
				}
				sb.deleteCharAt(sb.length()-1);
				try {
					bw.write(sb.toString());
					bw.newLine();
				}
				catch (Exception e) {
					System.out.println(e.toString());
				}
				return true;
			}
			return false;
		}

		void writeHapCount (BufferedWriter bw, int[] taxaIndex, double lMAF, double hMAF, double lCoverage, double hCoverage, int ID) {
			if (this.checkIfOutput(taxaIndex, lMAF, hMAF, lCoverage, hCoverage)) {
				int cnt1 = 0, cnt2 = 0;
				StringBuilder sb = new StringBuilder();
				sb.append("TP").append(ID).append("\t");
				for (int i = 0; i < taxaIndex.length; i++) {
					sb.append(this.queryCounts[taxaIndex[i]]).append("|").append(this.hitCounts[taxaIndex[i]]).append("\t");
					cnt1 += this.queryCounts[taxaIndex[i]];
					cnt2 += this.hitCounts[taxaIndex[i]];
				}
				double fre = (double)cnt1/(double)(cnt1+cnt2);
				sb.append(cnt1).append("\t").append(cnt2).append("\t").append(fre);
				try {
					bw.write(sb.toString());
					bw.newLine();
				}
				catch (Exception e) {
					System.out.println(e.toString());
				}
			}
		}

		void writeHapSeq (BufferedWriter bw, int[] taxaIndex, double lMAF, double hMAF, double lCoverage, double hCoverage, int ID) {
			if (this.checkIfOutput(taxaIndex, lMAF, hMAF, lCoverage, hCoverage)) {
				try {
					bw.write(">TP"+ID+"_query");
					bw.newLine();
					bw.write(BaseEncoder.getSequenceFromLong(query));
					bw.newLine();
					bw.write(">TP"+ID+"_hit");
					bw.newLine();
					bw.write(BaseEncoder.getSequenceFromLong(hit));
					bw.newLine();
				}
				catch (Exception e) {
					System.out.println(e.toString());
				}
			}
		}

		boolean checkIfOutput (int[] taxaIndex, double lMAF, double hMAF, double lCoverage, double hCoverage) {
			double MAF;
			double coverage;
			int cntA = 0, cntC = 0, co = 0;
			for (int i = 0; i < taxaIndex.length; i++) {
				cntA += queryCounts[taxaIndex[i]];
				cntC += hitCounts[taxaIndex[i]];
				if (queryCounts[taxaIndex[i]] > 0 || hitCounts[taxaIndex[i]] > 0) {
					co++;
				}
			}
			coverage = (double)co / (double)taxaIndex.length;
			if (coverage < lCoverage || coverage > hCoverage) return false;
			int sumCount = cntA + cntC;
			int minorCount = cntA;
			if (cntA > cntC) minorCount = cntC;
			MAF = (double)minorCount/(double)sumCount;
			if (MAF < lMAF || MAF > hMAF) return false;
			return true;
		}
	}
}
