/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.gbs.pav;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeSet;

/**
 *
 * @author Fei Lu
 */
public class TaxaTagCountInBin {
	int binSize = 100000;
	int taxaNum;
	String[] taxaName;
	int chrNum;
	int[] chrID;
	int[] binNumOfChr;
	int[][] binPosStart;
	int[][] binPosEnd;
	short[][][] binCount;

	public TaxaTagCountInBin (String genomeInfoFileS, String cnvTBTFileS, String tagPMapFileS, int binSize) {
		this.binSize = binSize;
		this.calBins(genomeInfoFileS, binSize);
		this.calBinCount(cnvTBTFileS, tagPMapFileS);
	}

	public TaxaTagCountInBin (String binCountFileS) {
		this.readBinaryFile(binCountFileS);
	}

	public void mergeTaxa () {
		String[] names = new String[taxaNum];
		for (int i = 0; i < taxaNum; i++) {
			String[] temp = taxaName[i].split(":");
			names[i] = temp[0];
		}
		TreeSet<String> nameSet = new TreeSet();
		for (int i = 0; i < names.length; i++) nameSet.add(names[i]);
		String[] nonreName = nameSet.toArray(new String[nameSet.size()]);
		Arrays.sort(nonreName);
		ArrayList<Integer>[] indexList = new ArrayList[nonreName.length];
		int[][] index = new int[nonreName.length][];
		for (int i = 0; i < nonreName.length; i++) {
			indexList[i] = new ArrayList();
			for (int j = 0; j < names.length; j++) {
				if (nonreName[i].equals(names[j])) {
					indexList[i].add(j);
				}
			}
			index[i] = new int[indexList[i].size()];
			for (int j = 0; j < index[i].length; j++) {
				index[i][j] = indexList[i].get(j);
			}
		}
		short[][][] newBinCount = new short[chrNum][][];
		for (int i = 0; i < chrNum; i++) {
			newBinCount[i] = new short[binNumOfChr[i]][nonreName.length];
			for (int j = 0; j < binNumOfChr[i]; j++) {
				for (int k = 0; k < nonreName.length; k++) {
					newBinCount[i][j][k] = 0;
					for (int u = 0; u < index[k].length; u++) {
						newBinCount[i][j][k]+= binCount[i][j][index[k][u]];
					}
				}
			}
		}
		taxaName = nonreName;
		taxaNum = taxaName.length;
		binCount = newBinCount;
	}

	public void readBinaryFile (String infileS) {
		try {
			DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(infileS), 65536));
			binSize = dis.readInt();
			taxaNum = dis.readInt();
			taxaName = new String[taxaNum];
			for (int i = 0; i < taxaNum; i++) taxaName[i] = dis.readUTF();
			chrNum = dis.readInt();
			chrID = new int[chrNum];
			binNumOfChr = new int[chrNum];
			binPosStart = new int[chrNum][];
			binPosEnd = new int[chrNum][];
			binCount = new short[chrNum][][];
			for (int i = 0; i < chrNum; i++) {
				chrID[i] = dis.readInt();
				binNumOfChr[i] = dis.readInt();
				binPosStart[i] = new int[binNumOfChr[i]];
				binPosEnd[i] = new int[binNumOfChr[i]];
				binCount[i] = new short[binNumOfChr[i]][taxaNum];
				for (int j = 0; j < binNumOfChr[i]; j++) {
					binPosStart[i][j] = dis.readInt();
					binPosEnd[i][j] = dis.readInt();
					for (int k = 0; k < taxaNum; k++) {
						binCount[i][j][k] = dis.readShort();
					}
				}
			}
		}
		catch (Exception e) {
				System.out.printf(e.toString());
		}
	}

	public void writeBinaryFile (String outfileS) {
		try {
			DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfileS), 65536));
			dos.writeInt(binSize);
			dos.writeInt(taxaNum);
			for (int i = 0; i < taxaNum; i++) dos.writeUTF(taxaName[i]);
			dos.writeInt(chrNum);
			for (int i = 0; i < chrNum; i++) {
				dos.writeInt(chrID[i]);
				dos.writeInt(binNumOfChr[i]);
				for (int j = 0; j < binNumOfChr[i]; j++) {
					dos.writeInt(binPosStart[i][j]);
					dos.writeInt(binPosEnd[i][j]);
					for (int k = 0; k < taxaNum; k++) {
						dos.writeShort(binCount[i][j][k]);
					}
				}
			}
			dos.flush();
			dos.close();
		}
		catch (Exception e) {
			System.out.println(e.toString());
		}
	}

	public void writeTxtFile (String outfileS) {
		try {
			BufferedWriter bw = new BufferedWriter (new FileWriter (outfileS), 65536);
			bw.write("Chr\tBinStart\tBinEnd\t");
			for (int i = 0; i < taxaName.length; i++) {
				bw.write(taxaName[i]+"\t");
			}
			bw.newLine();
			for (int i = 0; i < chrID.length; i++) {
				for (int j = 0; j < binNumOfChr[i]; j++) {
					bw.write(String.valueOf(chrID[i])+"\t"+String.valueOf(binPosStart[i][j])+"\t"+String.valueOf(binPosEnd[i][j])+"\t");
					for (int k = 0; k < taxaNum; k++) {
						bw.write(binCount[i][j][k]+"\t");
					}
					bw.newLine();
				}
			}
			bw.flush();
			bw.close();
		}
		catch (Exception e) {
			System.out.println(e.toString());
		}
	}

	public void calBinCount (String cnvTBTFileS, String tagPMapFileS) {
		TagPMap pmap  = new TagPMap(tagPMapFileS);
		try {
			DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(cnvTBTFileS), 65536*1024));
			int tagNum = dis.readInt();
            int tagLengthInLong = dis.readInt();
            taxaNum = dis.readInt();
			taxaName = new String[taxaNum];
			for (int i = 0; i < taxaNum; i++) {
				taxaName[i] = dis.readUTF();
			}
			for (int i = 0; i < chrID.length; i++) {
				for (int j = 0; j < binNumOfChr[i]; j++) {
					binCount[i][j] = new short[taxaName.length];
					for (int k = 0; k < taxaName.length; k++) {
						binCount[i][j][k] = 0;
					}
				}
			}
			byte[] tagDis = new byte[taxaNum];
			for (int i = 0; i < tagNum; i++) {
				long[] tempTag = new long[tagLengthInLong];
				for (int j = 0; j < tagLengthInLong; j++) {
					tempTag[j] = dis.readLong();
				}
				dis.readByte();
				int ind = pmap.getTagIndex(tempTag);
				if (pmap.hitNum[ind] != 1) continue;
				if (pmap.pChr[ind][0] < 1 || pmap.pChr[ind][0] > 10) continue;
				
				int chr = pmap.pChr[ind][0];
				int pos = pmap.pPos[ind][0];
				int chrIndex = Arrays.binarySearch(chrID, chr);
				if (chrIndex < 0) continue;
				int hit = Arrays.binarySearch(binPosStart[chrIndex], pos);
				if (hit < 0) hit = -hit-2;
				for (int j = 0; j < taxaNum; j++) {
					tagDis[j] = dis.readByte();
					binCount[chrIndex][hit][j]+=tagDis[j];
				}
			}
		}
		catch (Exception e) {
			System.out.println(e.toString());
		}
	}

	public void calBins (String genomeInfoFileS, int binSize) {
		Table t = new Table(genomeInfoFileS);
		chrID = new int[t.getRowNumber()];
		chrNum = chrID.length;
		binNumOfChr = new int[chrID.length];
		binPosStart = new int[chrID.length][];
		binPosEnd = new int[chrID.length][];
		binCount = new short[chrID.length][][];
		for (int i = 0; i < chrID.length; i++) {
			chrID[i] = Integer.valueOf(t.content[i][0]);
			int binNum = Integer.valueOf(t.content[i][1])/binSize + 1;
			binNumOfChr[i] = binNum;
			binPosStart[i] = new int[binNum];
			binPosEnd[i] = new int[binNum];
			binCount[i] = new short[binNum][];
			for (int j = 0; j < binNum; j++) {
				binPosStart[i][j] = binSize * j;
				binPosEnd[i][j] = binSize * (j+1);
			}
		}
	}
}
