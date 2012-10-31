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
import java.util.TreeSet;
import net.maizegenetics.genome.BaseEncoder;

/**
 *
 * @author Fei Lu
 */
public class NoRefHapUtils1 {
	HapRecord[] hr;
	String[] taxaNames;
	int taxaNum;
	int hapRecordNum;
	int[] taxaIndex;

	public NoRefHapUtils1 () {

	}

	public NoRefHapUtils1 (ReadsByTaxa rbt, String infileS) {
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

	public NoRefHapUtils1 (String infileS) {
		this.readMapInfo(infileS);
		this.getTaxaIndex();
	}

	public void mergeMapInfo (File[] hapSeqFiles, File tagCountFolder) {
		TreeSet<HapRecord> hrSet = new TreeSet();
		File[] tagCountFiles = tagCountFolder.listFiles();
		taxaNum = tagCountFiles.length;
		taxaNames = new String[taxaNum];
		for (int i = 0; i < hapSeqFiles.length; i++) {
			try {
				BufferedReader br = new BufferedReader(new FileReader(hapSeqFiles[i]),65536);
				String temp;
				while (((temp = br.readLine()) != null)) {
					if (temp.startsWith(">")) {
						temp = br.readLine();
						long[] query = new long[2];
						query[0]=BaseEncoder.getLongFromSeq(temp.substring(0, 32));
						query[1]=BaseEncoder.getLongFromSeq(temp.substring(32, 64));
						br.readLine();
						temp = br.readLine();
						long[] hit = new long[2];
						hit[0]=BaseEncoder.getLongFromSeq(temp.substring(0, 32));
						hit[1]=BaseEncoder.getLongFromSeq(temp.substring(32, 64));
						HapRecord onehr = new HapRecord(query, hit, taxaNum);
						hrSet.add(onehr);
					}
				}
			}
			catch (Exception e) {
			}
		}
		hr = hrSet.toArray(new HapRecord[hrSet.size()]);
		hapRecordNum = hr.length;
		for (int i = 0; i < tagCountFiles.length; i++) {
			String temp = tagCountFiles[i].getName();
			taxaNames[i] = temp.split("_")[0];
			ReadCounts rc = new ReadCounts (tagCountFiles[i].getAbsolutePath(), true);
			for (int j = 0; j < hr.length; j++) {
				int index = rc.getReadIndex(hr[j].query);
				if (index < 0) {
					hr[j].queryCounts[i] = 0;
				}
				else {
					hr[j].queryCounts[i] = (byte)rc.getReadCount(index);
				}
				index = rc.getReadIndex(hr[j].hit);
				if (index < 0) {
					hr[j].hitCounts[i] = 0;
				}
				else {
					hr[j].hitCounts[i] = (byte)rc.getReadCount(index);
				}
			}
		}
		for (int i = 0; i < hr.length; i++) {
			hr[i].updateStats();
		}
		this.getTaxaIndex();
	}

	private void checkFile(ReadsByTaxa rbt) {
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

	private void getPairedTBT (ReadsByTaxa rbt, TagPair tp) {
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

	private void getTaxaIndex () {
		taxaIndex = new int[taxaNum];
		for (int i = 0; i < taxaIndex.length; i++) {
			taxaIndex[i] = i;
		}
	}

	private void getTaxaIndex (String[] selectedTaxa) {
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

	private void getHapRecord(ReadsByTaxa rbt) {
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

	public void heterozygoteTest (double chiSquareCutoff) {
		for (int i = 0; i < hr.length; i++) {
			hr[i].chiSquareTestHetSum(chiSquareCutoff, taxaIndex);
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

	private void readMapInfo (String infileS) {
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

	public void writeHapMap (String outfileS) {
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
				if (hr[i].writeHapMap(bw, taxaIndex, i+1)) total++;
			}
			bw.flush();
			bw.close();
			System.out.println("HapMap file is written. There are " + total + " SNPs in total");
		}
		catch (Exception e) {
			System.out.println("Error occurred while writing HapMap file " + outfileS + " " + e.toString());
		}
	}

	public void writeHapCount (String outfileS) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outfileS), 65536);
			bw.write("rs\t");
			for (int i = 0; i < this.taxaIndex.length; i++) {
				bw.write(this.taxaNames[this.taxaIndex[i]] + "\t");
			}
			bw.write("HetCount_allele1\tHetCount_allele2\tCount_allele1\tCount_allele2\tFrequency");
			bw.newLine();
			int SNPCount = 0;
			int readCount = 0;
			for (int i = 0; i < this.hapRecordNum; i++) {
				int cnt = hr[i].writeHapCount(bw, taxaIndex, i+1);
				if (cnt != -1) {
					SNPCount++;
					readCount += cnt;
				}
			}
			bw.flush();
			bw.close();
			double cover = (double)readCount/(double)SNPCount/taxaIndex.length;
			System.out.println("HapMapCount file is written. The coverage is " + String.valueOf(cover));
		}
		catch (Exception e) {
			System.out.println("Error occurred while writing HapMap file " + outfileS + " " + e.toString());
		}
	}

	public void writeHapSeq (String outfileS) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outfileS), 65536);
			for (int i = 0; i < this.hapRecordNum; i++) {
				hr[i].writeHapSeq(bw, taxaIndex, i+1);
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
		for (int i = 0; i < hr.length; i++) {
			hr[i].checkIfOutput(taxaIndex, lMAF, hMAF, lCoverage, hCoverage);
		}
		this.writeHapMap(hapMapF.toString());
		this.writeHapCount(hapCountF.toString());
		this.writeHapSeq(hapSeqF.toString());
	}

	public void homoParentCross (String parentA, String parentB) {
		ReadCounts rcA = new ReadCounts (parentA, true);
		ReadCounts rcB = new ReadCounts (parentB, true);
		ArrayList<HapRecord> hapRList = new ArrayList();
		for (int i = 0; i < this.hapRecordNum; i++) {
			int hit1 = rcA.binarySearchOfHaps(this.hr[i].query);
			int hit2 = rcA.binarySearchOfHaps(this.hr[i].hit);
			int hit3 = rcB.binarySearchOfHaps(this.hr[i].query);
			int hit4 = rcB.binarySearchOfHaps(this.hr[i].hit);
			if (hit1 >= 0 && hit2 < 0) {
				if (hit3 < 0 && hit4 >= 0) {
					hapRList.add(hr[i]);
				}
			}
			else if (hit1 < 0 && hit2 >= 0) {
				if (hit3 >= 0 && hit4 < 0) {
					hapRList.add(hr[i]);
				}
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
					int j = i + 1;
					System.out.println("TP" + j + "\t" + this.hr[i].snpQ + "\t" + this.hr[i].hetCode);
				}
			}
			else if (hit1 < 0 && hit2 >= 0) {
				if (hit3 >= 0 && hit4 >= 0) {
					if (heteroRc.hapcount[hit3] < 5 || heteroRc.hapcount[hit4] < 5) continue;
					hapRList.add(hr[i]);
					int j = i + 1;
					System.out.println("TP" + j + "\t" + this.hr[i].snpH + "\t" + this.hr[i].hetCode);
				}
			}
		}
		HapRecord[] hapRArray = hapRList.toArray(new HapRecord[hapRList.size()]);
		hr = hapRArray;
		this.hapRecordNum = hr.length;
	}

	public void shallowPseudoTestCross (String homoParent, String heteroParent) {
		ReadCounts homoRc = new ReadCounts (homoParent, true);
		ReadCounts heteroRc = new ReadCounts (heteroParent, true);
		ArrayList<HapRecord> hapRList = new ArrayList();
		for (int i = 0; i < this.hapRecordNum; i++) {
			int hit3 = heteroRc.binarySearchOfHaps(this.hr[i].query);
			int hit4 = heteroRc.binarySearchOfHaps(this.hr[i].hit);
				if (hit3 >= 0 && hit4 >= 0) {
					if (heteroRc.hapcount[hit3] < 5 || heteroRc.hapcount[hit4] < 5) continue;
					hapRList.add(hr[i]);
					int j = i + 1;
					System.out.println("TP" + j + "\t" + this.hr[i].snpQ + "\t" + this.hr[i].hetCode);
				}
		}
		HapRecord[] hapRArray = hapRList.toArray(new HapRecord[hapRList.size()]);
		hr = hapRArray;
		this.hapRecordNum = hr.length;
	}

	private class HapRecord implements Comparable <HapRecord> {
		long[] query;
		long[] hit;
		byte[] queryCounts;
		byte[] hitCounts;
		int sumQueryCount;
		int sumHitCount;
		int hetQueryCount;
		int hetHitCount;
		char snpQ;
		char snpH;
		char hetCode;
		boolean ifSNP;

		HapRecord (long[] query, long[] hit, int taxaNum) {
			this.query = query;
			this.hit = hit;
			this.getSnps();
			queryCounts = new byte[taxaNum];
			hitCounts = new byte[taxaNum];
			ifSNP = true;
		}

		HapRecord (long[] query, long[] hit, byte[] a, byte[] b) {
			this.query = query;
			this.hit = hit;
			this.queryCounts = a;
			this.hitCounts = b;
			this.sumQueryCount = this.getSumCount(a);
			this.sumHitCount = this.getSumCount(b);
			int[] hetCounts = this.getHetCounts(a, b);
			this.hetQueryCount = hetCounts[0];
			this.hetHitCount = hetCounts[1];
			this.getSnps();
			ifSNP = true;
		}

		private void updateStats () {
			this.sumQueryCount = this.getSumCount(queryCounts);
			this.sumHitCount = this.getSumCount(hitCounts);
			int[] hetCounts = this.getHetCounts(queryCounts, hitCounts);
			this.hetQueryCount = hetCounts[0];
			this.hetHitCount = hetCounts[1];
		}

		private int[] getHetCounts (byte[] a, byte[] b) {
			int[] hetCounts = new int[2];
			for (int i = 0; i < a.length; i++) {
				if (a[i] > 0 && b[i] > 0) {
					hetCounts[0] += a[i];
					hetCounts[1] += b[i];
				}
			}
			return hetCounts;
		}

		private int[] getSubHetCounts (int[] taxaIndex) {
			int[] hetCounts = new int[2];
			for (int i = 0; i < taxaIndex.length; i++) {
				if (queryCounts[taxaIndex[i]] > 0 && hitCounts[taxaIndex[i]] > 0) {
					hetCounts[0] += queryCounts[taxaIndex[i]];
					hetCounts[1] += hitCounts[taxaIndex[i]];
				}
			}
			return hetCounts;
		}

		private int getSumCount (byte[] cnt) {
			int sum = 0;
			for (int i = 0; i < cnt.length; i++) {
				sum += cnt[i];
			}
			return sum;
		}

		private int getSubSumCountQuery (int[] taxaIndex) {
			int sum = 0;
			for (int i = 0; i < taxaIndex.length; i++) {
				sum += queryCounts[taxaIndex[i]];
			}
			return sum;
		}

		private int getSubSumCountHit (int[] taxaIndex) {
			int sum = 0;
			for (int i = 0; i < taxaIndex.length; i++) {
				sum += hitCounts[taxaIndex[i]];
			}
			return sum;
		}

		HapRecord (int nlong, int taxaNum) {
			this.query = new long[nlong];
			this.hit = new long[nlong];
			this.queryCounts = new byte[taxaNum];
			this.hitCounts = new byte[taxaNum];
		}

		private void getSnps () {
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

		private boolean chiSquareTestHetSum (double chiSquareCutOff, int[] taxaIndex) {
			int[] hetCounts = this.getSubHetCounts(taxaIndex);
			double chi = this.getChiSquareD1(hetCounts[0], hetCounts[1]);
			if (chi > chiSquareCutOff) {
				this.ifSNP = false;
			}
			return ifSNP;
		}

		private double getChiSquareD1 (int a, int b) {
			double E = ((double)a + (double)b)/2;
			double chi = Math.pow(((double)a-E), 2)/E + Math.pow(((double)b-E), 2)/E;
			return chi;
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
				for (int i = 0; i < hit.length; i++) {
					hit[i] = dis.readLong();
				}
				for (int i = 0; i < queryCounts.length; i++) {
					queryCounts[i] = dis.readByte();
				}
				for (int i = 0; i < hitCounts.length; i++) {
					hitCounts[i] = dis.readByte();
				}
				sumQueryCount = dis.readInt();
				sumHitCount = dis.readInt();
				hetQueryCount = dis.readInt();
				hetHitCount = dis.readInt();
				snpQ = dis.readChar();
				snpH = dis.readChar();
				hetCode = dis.readChar();
				ifSNP = dis.readBoolean();
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
				for (int i = 0; i < hit.length; i++) {
					dos.writeLong(hit[i]);
				}
				for (int i = 0; i < queryCounts.length; i++) {
					dos.writeByte(queryCounts[i]);
				}
				for (int i = 0; i < hitCounts.length; i++) {
					dos.writeByte(hitCounts[i]);
				}
				dos.writeInt(sumQueryCount);
				dos.writeInt(sumHitCount);
				dos.writeInt(hetQueryCount);
				dos.writeInt(hetHitCount);
				dos.writeChar(snpQ);
				dos.writeChar(snpH);
				dos.writeChar(this.hetCode);
				dos.writeBoolean(ifSNP);
			}
			catch (Exception e) {
				System.out.println(e.toString());
			}
		}

		boolean writeHapMap (BufferedWriter bw, int[] taxaIndex, int ID) {
			if (!ifSNP) return false;
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

		int writeHapCount (BufferedWriter bw, int[] taxaIndex, int ID) {
			if (!ifSNP) return -1;
			StringBuilder sb = new StringBuilder();
			sb.append("TP").append(ID).append("\t");
			for (int i = 0; i < taxaIndex.length; i++) {
				sb.append(this.queryCounts[taxaIndex[i]]).append("|").append(this.hitCounts[taxaIndex[i]]).append("\t");
			}
			int subSumQueryCount = this.getSubSumCountQuery(taxaIndex);
			int subSumHitCount = this.getSubSumCountHit(taxaIndex);
			double fre = (double)subSumQueryCount/(double)(subSumQueryCount+subSumHitCount);
			int[] subHetCounts = this.getSubHetCounts(taxaIndex);
			sb.append(subHetCounts[0]).append("\t").append(subHetCounts[1]).append("\t").append(this.getSubSumCountQuery(taxaIndex)).append("\t").append(this.getSubSumCountHit(taxaIndex)).append("\t").append(fre);
			try {
				bw.write(sb.toString());
				bw.newLine();
			}
			catch (Exception e) {
				System.out.println(e.toString());
			}
			return subSumQueryCount + subSumHitCount;
		}

		void writeHapSeq (BufferedWriter bw, int[] taxaIndex, int ID) {
			if (!ifSNP) return;
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

		boolean checkIfOutput (int[] taxaIndex, double lMAF, double hMAF, double lCoverage, double hCoverage) {
			if (!ifSNP) return false;
			double MAF;
			double coverage;
			int cntA = 0, cntC = 0, co = 0;
			cntA = this.getSubSumCountQuery(taxaIndex);
			cntC = this.getSubSumCountHit(taxaIndex);
			if (cntA * cntC == 0){
				this.ifSNP = false;
				return false;
			}
			for (int i = 0; i < taxaIndex.length; i++) {
				if (queryCounts[taxaIndex[i]] > 0 || hitCounts[taxaIndex[i]] > 0) {
					co++;
				}
			}
			coverage = (double)co / (double)taxaIndex.length;
			if (coverage < lCoverage || coverage > hCoverage) {
				this.ifSNP = false;
				return false;
			}
			int sumCount = cntA + cntC;
			int minorCount = cntA;
			if (cntA > cntC) minorCount = cntC;
			MAF = (double)minorCount/(double)sumCount;
			if (MAF < lMAF || MAF > hMAF) {
				this.ifSNP = false;
				return false;
			}
			return true;
		}

		public int compareTo(HapRecord o) {
			int v = this.compareLongArray(this.query, o.query);
			if (v != 0) return this.compareLongArray(this.hit, o.hit);
			else return v;
		}

		private int compareLongArray (long[] l1, long[] l2) {
			int v = 0;
			for (int i = 0; i < l1.length; i++) {
				v = this.compareLong(l1[i], l2[i]);
				if (v != 0) return v;
			}
			return v;
		}

		private int compareLong (long l1, long l2) {
			if (l1 == l2) return 0;
			if (l1 < l2) return -1;
			return 1;
		}
	}
}
