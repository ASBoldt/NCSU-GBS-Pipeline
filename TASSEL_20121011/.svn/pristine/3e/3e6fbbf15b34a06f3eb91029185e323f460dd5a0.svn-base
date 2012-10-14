/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.gbs.pipeline;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import net.maizegenetics.gbs.tagdist.UTagPairs;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByte;
import net.maizegenetics.genome.BaseEncoder;
/**
 *
 * @author Fei Lu
 */
public class UMapInfo {
	HapRecord[] hr;
	String[] taxaNames;
	int taxaNum;
	int tagLengthInLong;
	int hapRecordNum;

	public UMapInfo (String mapInfoFileS) {
		this.readMapInfo(mapInfoFileS);
	}

	public UMapInfo (TagsByTaxaByte tbt, UTagPairs tp) {
		this.taxaNum = tbt.getTaxaCount();
		this.tagLengthInLong = tbt.getTagSizeInLong();
		this.hapRecordNum = tbt.getTagCount()/2;
		this.taxaNames = tbt.getTaxaNames();
		this.hr = new HapRecord[hapRecordNum];	
		this.getHapRecord(tbt, tp);
	}

	private void checkFile (TagsByTaxaByte tbt, UTagPairs tp) {
		if (tbt.getTagCount() != tp.getTagNum()) {
			System.out.println("TBT file and TagPair file don't match");
			System.exit(0);
		}
		System.out.println("TBT file and TagPair file match");
	}

	private void getHapRecord (TagsByTaxaByte tbt, UTagPairs tp) {
		this.checkFile(tbt, tp);
		hr = new HapRecord[tp.getTagNum()/2];
		int cnt = 0;
		for (int i = 0; i < tp.getTagNum(); i += 2) {
			long[] query = tp.getTag(i);
			long[] hit = tp.getTag(i+1);
			byte[] queryCounts = tbt.getTaxaReadCountsForTag(tbt.getTagIndex(query));
			byte[] hitCounts = tbt.getTaxaReadCountsForTag(tbt.getTagIndex(hit));
			hr[cnt] = new HapRecord (query, hit, tp.getTagLength(i), tp.getTagLength(i+1), queryCounts, hitCounts);
			cnt++;
		}
	}

	public void readMapInfo (String infileS) {
		try {
			DataInputStream dis =new DataInputStream(new BufferedInputStream(new FileInputStream(infileS),65536));
			this.taxaNum = dis.readInt();
			this.tagLengthInLong = dis.readInt();
			this.hapRecordNum = dis.readInt();
			this.taxaNames = new String[this.taxaNum];
			for (int i = 0; i <  this.taxaNum; i++) {
				this.taxaNames[i] = dis.readUTF();
			}
			hr = new HapRecord[this.hapRecordNum];
			for (int i = 0; i < this.hapRecordNum; i++) {
				hr[i] = new HapRecord(tagLengthInLong, taxaNum);
				hr[i].readRecord(dis);
			}
			System.out.println(hr.length + " HapMap records read from the MapInfo file " + infileS);

		}
		catch (Exception e) {
			System.out.println("Error occurred while reading MapInfo file " + infileS + " " + e.toString());
		}
	}

	public void writeHapAll (String outDirS, double mnMAF, double mxMAF, double mnCall, double mxCall) {
		File hapMapF = new File(outDirS, "HapMap.hmp.txt");
		File hapCountF = new File(outDirS, "HapMap.hmc.txt");
		File hapSeqF = new File(outDirS, "HapMap.fas.txt");
		for (int i = 0; i < hr.length; i++) {
			hr[i].checkIfOutput(mnMAF, mxMAF, mnCall, mxCall);
		}
		this.writeHapMap(hapMapF.toString());
		this.writeHapCount(hapCountF.toString());
		this.writeHapSeq(hapSeqF.toString());
		System.out.println("\nUniversal Network Enabled Analysis Kit(UNEAK) is an unpublished work.");
		System.out.println("Please check in our website http://www.maizegenetics.net/ and cite the paper \"Switchgrass genomic diversity, ploidy and evolution: novel insights from a network-based SNP discovery protocol\" by Fei Lu, Alexander E. Lipka, Rob Elshire, Jeff Glaubitz, Jerome H. Cherney, Michael D. Casler, Edward Buckler, Denise E. Costich. (Citation 2012).");
	}

	public void writeHapMap (String outfileS) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outfileS), 65536);
			bw.write("rs	alleles	chrom	pos	strand	assembly	center	protLSID	assayLSID	panelLSID	QCcode	");
			for (int i = 0; i < this.taxaNames.length - 1; i++) {
				bw.write(this.taxaNames[i] + "\t");
			}
			bw.write(this.taxaNames[this.taxaNames.length-1]);
			bw.newLine();
			int total = 0;
			for (int i = 0; i < this.hapRecordNum; i++) {
				if (hr[i].writeHapMap(bw, i+1)) total++;
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
			for (int i = 0; i < this.taxaNames.length; i++) {
				bw.write(this.taxaNames[i] + "\t");
			}
			bw.write("HetCount_allele1\tHetCount_allele2\tCount_allele1\tCount_allele2\tFrequency");
			bw.newLine();
			for (int i = 0; i < this.hapRecordNum; i++) {
				hr[i].writeHapCount(bw, i+1);
			}
			bw.flush();
			bw.close();
			System.out.println("HapMapCount file is written");
		}
		catch (Exception e) {
			System.out.println("Error occurred while writing HapMap file " + outfileS + " " + e.toString());
		}
	}

	public void writeHapSeq (String outfileS) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outfileS), 65536);
			for (int i = 0; i < this.hapRecordNum; i++) {
				hr[i].writeHapSeq(bw, i+1);
			}
			bw.flush();
			bw.close();
			System.out.println("HapMapSeq file is written");
		}
		catch (Exception e) {
			System.out.println("Error occurred while writing HapMap file " + outfileS + " " + e.toString());
		}
	}

	public void writeMapInfo (String outfileS) {
		try {
			DataOutputStream dos =new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfileS),65536));
			dos.writeInt(this.taxaNum);
			dos.writeInt(this.tagLengthInLong);
			dos.writeInt(this.hapRecordNum);
			for (int i = 0; i < this.taxaNum; i++) {
				dos.writeUTF(this.taxaNames[i]);
			}
			for (int i = 0; i < this.hapRecordNum; i++) {
				hr[i].writeRecord(dos);
			}
			dos.flush();
			dos.close();
			System.out.println(String.valueOf(hapRecordNum) + " HapMap records written to MapInfo file " + outfileS);
		}
		catch (Exception e) {
			System.out.println("Error occurred while writing MapInfo file " + outfileS + " " + e.toString());
		}
	}

	private class HapRecord {
		long[] query;
		long[] hit;
		byte queryLength;
		byte hitLength;
		byte[] queryCounts;
		byte[] hitCounts;
		int sumQueryCount;
		int sumHitCount;
		int hetQueryCount;
		int hetHitCount;
		char snpQuery;
		char snpHit;
		char hetCode;
		boolean ifSNP;

		HapRecord (long[] query, long[] hit, byte queryLength, byte hitLength, byte[] queryCounts, byte[] hitCounts) {
			this.query = query;
			this.hit = hit;
			this.queryLength = queryLength;
			this.hitLength = hitLength;
			this.queryCounts = queryCounts;
			this.hitCounts = hitCounts;
			this.sumQueryCount = this.getSumCount(queryCounts);
			this.sumHitCount = this.getSumCount(hitCounts);
			int[] hetCounts = this.getHetCounts(queryCounts, hitCounts);
			this.hetQueryCount = hetCounts[0];
			this.hetHitCount = hetCounts[1];
			this.getSnps();
			ifSNP = true;
		}

		HapRecord (int nlong, int taxaNum) {
			this.query = new long[nlong];
			this.hit = new long[nlong];
			this.queryCounts = new byte[taxaNum];
			this.hitCounts = new byte[taxaNum];
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

		private int getSumCount (byte[] cnt) {
			int sum = 0;
			for (int i = 0; i < cnt.length; i++) {
				sum += cnt[i];
			}
			return sum;
		}

		private void getSnps () {
			for (int i = 0; i < query.length; i++) {
				int snpIndex = this.getSnpIndex(query[i], hit[i]);
				if (snpIndex != -1) {
					String q = BaseEncoder.getSequenceFromLong(query[i]);
					snpQuery = q.charAt(snpIndex);
					String h = BaseEncoder.getSequenceFromLong(hit[i]);
					snpHit = h.charAt(snpIndex);
					this.getHetCode();
					return;
				}
			}
		}

		private void getHetCode () {
			if (snpQuery + snpHit == 132) hetCode = 'M';
			else if (snpQuery + snpHit == 136) hetCode = 'R';
			else if (snpQuery + snpHit == 149) hetCode = 'W';
			else if (snpQuery + snpHit == 138) hetCode = 'S';
			else if (snpQuery + snpHit == 151) hetCode = 'Y';
			else if (snpQuery + snpHit == 155) hetCode = 'K';
		}

		private boolean chiSquareTestHetSum (double chiSquareCutOff) {
			double chi = this.getChiSquareD1(this.hetQueryCount, this.hetHitCount);
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
				queryLength = dis.readByte();
				hitLength = dis.readByte();
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
				snpQuery = dis.readChar();
				snpHit = dis.readChar();
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
				dos.writeByte(queryLength);
				dos.writeByte(hitLength);
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
				dos.writeChar(snpQuery);
				dos.writeChar(snpHit);
				dos.writeChar(this.hetCode);
				dos.writeBoolean(ifSNP);
			}
			catch (Exception e) {
				System.out.println(e.toString());
			}
		}

		boolean writeHapMap (BufferedWriter bw, int ID) {
			if (!ifSNP) return false;
			StringBuilder sb = new StringBuilder();
			sb.append("TP").append(ID).append("\t").append(snpQuery).append("/").append(snpHit).append("\t").append(0).append("\t").append(ID).append("\t+\tNA\tUNEAK\tGBS\tNoRef\tCustom\tQC+\t");
			for (int i = 0; i < queryCounts.length; i++) {
				if (this.queryCounts[i] == 0) {
					if (this.hitCounts[i] == 0) {
						sb.append("N").append("\t");
					}
					else {
						sb.append(this.snpHit).append("\t");
					}
				}
				else {
					if (this.hitCounts[i] == 0) {
						sb.append(this.snpQuery).append("\t");
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

		void writeHapCount (BufferedWriter bw, int ID) {
			if (!ifSNP) return;
			StringBuilder sb = new StringBuilder();
			sb.append("TP").append(ID).append("\t");
			for (int i = 0; i < queryCounts.length; i++) {
				sb.append(this.queryCounts[i]).append("|").append(this.hitCounts[i]).append("\t");
			}
			double fre = (double)this.sumQueryCount/(double)(this.sumQueryCount+this.sumHitCount);
			sb.append(hetQueryCount).append("\t").append(hetHitCount).append("\t").append(sumQueryCount).append("\t").append(sumHitCount).append("\t").append("0.").append(Math.round(fre*1000));
			try {
				bw.write(sb.toString());
				bw.newLine();
			}
			catch (Exception e) {
				System.out.println(e.toString());
			}

		}

		void writeHapSeq (BufferedWriter bw, int ID) {
			if (!ifSNP) return;
			try {
				bw.write(">TP"+ID+"_query"+"_"+String.valueOf(queryLength));
				bw.newLine();
				bw.write(BaseEncoder.getSequenceFromLong(query));
				bw.newLine();
				bw.write(">TP"+ID+"_hit"+"_"+String.valueOf(hitLength));
				bw.newLine();
				bw.write(BaseEncoder.getSequenceFromLong(hit));
				bw.newLine();
			}
			catch (Exception e) {
				System.out.println(e.toString());
			}
		}

		boolean checkIfOutput (double mnMAF, double mxMAF, double mnCall, double mxCall) {
			if (!ifSNP) return false;
			double MAF;
			double callRate;
			int cntA = 0, cntC = 0, co = 0;
			for (int i = 0; i < queryCounts.length; i++) {
				cntA += queryCounts[i];
				cntC += hitCounts[i];
				if (queryCounts[i] > 0 || hitCounts[i] > 0) {
					co++;
				}
			}
			callRate = (double)co / (double)queryCounts.length;
			if (callRate < mnCall || callRate > mxCall) {
				this.ifSNP = false;
				return false;
			}
			int sumCount = cntA + cntC;
			int minorCount = cntA;
			if (cntA > cntC) minorCount = cntC;
			MAF = (double)minorCount/(double)sumCount;
			if (MAF < mnMAF || MAF > mxMAF) {
				this.ifSNP = false;
				return false;
			}
			return true;
		}
	}
}
