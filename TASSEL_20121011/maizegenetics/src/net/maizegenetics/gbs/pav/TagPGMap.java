/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.gbs.pav;

import cern.colt.GenericSorting;
import cern.colt.function.IntComparator;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.genome.BaseEncoder;

/**
 * Tags don't align have the -1 at pChr and pPos. Tags align to multiple position were asigned the pChr and pPos of the best position
 * @author Fei Lu
 */
public class TagPGMap extends TagGMap {
	byte[] pChr;
	int[] pPos;
	boolean[] ifMatch;
	boolean[] ifPerfectMatch;
	boolean[] ifUniqueMatch;

	public TagPGMap (String tagPMapFileS, String tagGMapFileS, String tagCountFileS) {
		super.readTxtFile(tagGMapFileS);
		TagPMap pMap = new TagPMap (tagPMapFileS);
		TagCounts tc = new TagCounts (tagCountFileS, FilePacking.Byte);
		this.mergePMap(pMap, tc);
        this.orderFamilyGroupByP();
	}

	public TagPGMap (String tagPGMapFileS, boolean ifBinary) {
		if (ifBinary) this.readBinaryFile(tagPGMapFileS);
		else this.readTxtFile(tagPGMapFileS);
		this.sortByTag();
	}

	public void mergePMap (TagPMap pMap, TagCounts tc) {
		tagLength = new byte[this.getTagCount()];
		pChr = new byte[this.getTagCount()];
		pPos = new int[this.getTagCount()];
		ifMatch = new boolean[this.getTagCount()];
		ifPerfectMatch = new boolean[this.getTagCount()];
		ifUniqueMatch = new boolean[this.getTagCount()];
		for (int i = 0; i < this.getTagCount(); i++) {
			int index = pMap.getTagIndex(this.getTag(i));
			if (index < 0) {
				pChr[i] = -1;
				pPos[i] = -1;
				ifMatch[i] = false;
				ifPerfectMatch[i] = false;
				ifUniqueMatch[i] = false;
			}
			else {
				ifMatch[i] = true;
				if (pMap.hitNum[index] == 1) {
					ifUniqueMatch[i] = true;
				}
				else {
					ifUniqueMatch[i] = false;
				}
				ifPerfectMatch[i] = pMap.ifPerfectMatch[index][0];
				pChr[i] = pMap.pChr[index][0];
				pPos[i] = pMap.pPos[index][0];
			}
			int hit = tc.getTagIndex(this.getTag(i));
			tagLength[i] = (byte)tc.getTagLength(hit);
		}
		System.out.println("TagGMap and TagPMap is merged");
	}

    public void orderFamilyGroupByP () {
        for (int i = 0; i < this.getTagCount(); i++) {
            if (this.gNum[i] == 1) continue;
            double[] p = new double[this.gNum[i]];
            int[] index = new int[this.gNum[i]];
            for (int j = 0; j < p.length; j++) {
                p[j] = this.binomP[i][j];
                index[j] = j;
            }
            for (int j = 0; j < p.length - 1; j++) {
                for (int k = j + 1; k < p.length; k++) {
                    if (p[j] > p[k]) {
                        double tempP = p[j];
                        int tempIndex = index[j];
                        p[j] = p[k]; p[k] = tempP;
                        index[j] = index[k]; index[k] = tempIndex;
                    }
                }
            }
            this.orderFamilyGroupByIndex(index, i);
        }
        System.out.println("Family groups are ordered by p value");
    }
    
    private void orderFamilyGroupByIndex (int[] familyIndex, int tagIndex) {
        byte[] tempGChr = new byte[familyIndex.length];
        int[] tempSite = new int[familyIndex.length];
        int[] tempGPos = new int[familyIndex.length];
        double[] tempP = new double[familyIndex.length];
        int[] tempChrSig = new int[familyIndex.length];
        double[] tempLikelihood = new double[familyIndex.length];
        byte[] tempFamilyNum = new byte[familyIndex.length];
        long[] tempFamilyCode = new long[familyIndex.length];
        for (int i = 0; i < familyIndex.length; i++) {
            tempGChr[i] = this.gChr[tagIndex][familyIndex[i]];
            tempSite[i] = this.site[tagIndex][familyIndex[i]];
            tempGPos[i] = this.gPos[tagIndex][familyIndex[i]];
            tempP[i] = this.binomP[tagIndex][familyIndex[i]];
            tempChrSig[i] = this.chrSig[tagIndex][familyIndex[i]];
            tempLikelihood[i] = this.likelyhood[tagIndex][familyIndex[i]];
            tempFamilyNum[i] = this.familyNum[tagIndex][familyIndex[i]];
            tempFamilyCode[i] = this.familyCode[tagIndex][familyIndex[i]];
        }
        this.gChr[tagIndex] = tempGChr; 
        this.site[tagIndex] = tempSite;
        this.gPos[tagIndex] = tempGPos;
        this.binomP[tagIndex] = tempP;
        this.chrSig[tagIndex] = tempChrSig;
        this.likelyhood[tagIndex] = tempLikelihood;
        this.familyNum[tagIndex] = tempFamilyNum;
        this.familyCode[tagIndex] = tempFamilyCode;
        
    }
    
    public boolean[] getReverseFilter (boolean[] f1) {
        boolean[] f2 = new boolean[f1.length];
        int cnt = 0;
        for (int i = 0; i < f2.length; i++) {
            if (f1[i] == true) f2[i] = false;
            else {
                f2[i] = true;
                cnt++;
            }
        }
        System.out.println(String.valueOf(cnt) + " are in the reversed filter");
        return f2;
    }
    
    public boolean[] getFilterCross (boolean[] f1, boolean[] f2) {
        boolean[] fCross = new boolean[this.getTagCount()];
        int cnt = 0;
        for (int i = 0; i < this.getTagCount(); i++) {
            if (f1[i] == true && f2[i] == true) {
                fCross[i] = true;
                cnt++;
            }
            else fCross[i] = false;
        }
        System.out.println(String.valueOf(cnt) + " passed the cross filter");
        return fCross;
    }
    
    public boolean[] getIfUnderThresholdFirstFamilyGroup (double pThresh) {
        boolean[] ifUnderThresh = new boolean[this.getTagCount()];
        int cnt = 0;
        for (int i = 0; i < this.getTagCount(); i++) {
            if (this.binomP[i][0] > pThresh) ifUnderThresh[i] = false;
            else {
                ifUnderThresh[i] = true;
                cnt++;
            }
        }
        System.out.println(String.valueOf(cnt) + " tags are under threshold in the first family group");
        return ifUnderThresh;
    }
    
    public boolean[] getIfUnderThresholdEitherFamilyGroup (double pThresh) {
        boolean[] ifUnderThresh = new boolean[this.getTagCount()];
        int cnt = 0;
        for (int i = 0; i < this.getTagCount(); i++) {
            ifUnderThresh[i] = false;
            for (int j = 0; j < this.gNum[i]; j++) {
                if (this.binomP[i][0] <= pThresh) {
                    ifUnderThresh[i] = true;
                    cnt++;
                    break;
                }
            }
        }
        System.out.println(String.valueOf(cnt) + " tags are under threshold in either family group");
        return ifUnderThresh;
    }
    
    public boolean[] getIfUnderThresholdAllFamilyGroup (double pThresh) {
        boolean[] ifUnderThresh = new boolean[this.getTagCount()];
        int cnt = 0;
        for (int i = 0; i < this.getTagCount(); i++) {
            ifUnderThresh[i] = true;
            int c = 0;
            for (int j = 0; j < this.gNum[i]; j++) {
                if (this.binomP[i][0] <= pThresh) {
                    c++;
                    
                }
            }
            if (c == this.gNum[i]) {
                ifUnderThresh[i] = true;
                cnt++;
            }
            else ifUnderThresh[i] = false;
        }
        System.out.println(String.valueOf(cnt) + " tags are under threshold in all family groups");
        return ifUnderThresh;
    }
    
    public boolean[] getIfSingleMapping (double pThresh) {
        boolean[] ifSingle = new boolean[this.getTagCount()];
        int cnt = 0;
        for (int i = 0; i < this.getTagCount(); i++) {
            if (this.gNum[i] == 1) {
                ifSingle[i] = true;
                cnt++;
            }
            else {
                int c = 0;
                for (int j = 0; j < this.gNum[i]; j++) {
                    if (this.binomP[i][j] <= pThresh) c++;
                }
                if (c == 1) {
                    ifSingle[i] = true;
                    cnt++;
                }
                else ifSingle[i] = false;
            }
        }
        System.out.println(String.valueOf(cnt) + " tags are uniquely mapped");
        return ifSingle;
    }
    
    public boolean[] getIfMultipleMapping (double pThresh) {
        boolean[] ifMultiple = new boolean[this.getTagCount()];
        int cnt = 0;
        for (int i = 0; i < this.getTagCount(); i++) {
            if (this.gNum[i] > 1) {
                int c = 0;
                for (int j = 0; j < this.gNum[i]; j++) {
                   if (this.binomP[i][j] <= pThresh) c++;
                }
                if (c > 1) {
                    ifMultiple[i] = true;
                    cnt++;
                }
                else ifMultiple[i] = false;
            }
            else ifMultiple[i] = false;
        }
        System.out.println(String.valueOf(cnt) + " tags are multiplly mapped");
        return ifMultiple;
    }
    
    public boolean[] getIfRefTag() {
        boolean[] ifRef = new boolean[this.getTagCount()];
        int cnt = 0;
        for (int i = 0; i < this.getTagCount(); i++) {
            if (this.ifMatch[i] == true && this.ifPerfectMatch[i]) {
                ifRef[i] = true;
                cnt++;
            }
            else ifRef[i] = false;
        }
        System.out.println(String.valueOf(cnt) + " uniquely aligned ref tags");
        return ifRef;
    }
    
    public boolean[] getIfUniqueRefTag () {
        boolean[] ifRef = new boolean[this.getTagCount()];
        int cnt = 0;
        for (int i = 0; i < this.getTagCount(); i++) {
            if (this.ifMatch[i] == true && this.ifPerfectMatch[i] && this.ifUniqueMatch[i] == true) {
                ifRef[i] = true;
                cnt++;
            }
            else ifRef[i] = false;
        }
        System.out.println(String.valueOf(cnt) + " uniquely aligned ref tags");
        return ifRef;
    }
    
	public boolean[] getIfPAV (String samFileS) {
		boolean[] ifPAV = new boolean[this.getTagCount()];
		System.out.println("Start reading sam file");
		try {
			BufferedReader br = new BufferedReader (new FileReader(samFileS), 65536);
			String temp;
			ArrayList<String> tagAlnList = new ArrayList();
			ArrayList<Integer> tagIDList = new ArrayList();
			while ((temp = br.readLine()) != null) {
				if (temp.startsWith("@")) continue;
				if (tagAlnList.isEmpty()) {
					tagAlnList.add(temp);
					String[] tem = temp.split("\t");
					tagIDList.add(Integer.valueOf(tem[0]));
				}
				else {
					String[] tem = temp.split("\t");
					int tempID = Integer.valueOf(tem[0]);
					if (tempID == tagIDList.get(0)) {
						tagAlnList.add(temp);
						tagIDList.add(Integer.valueOf(tem[0]));
					}
					else {
						checkIfPAV(tagAlnList, tagIDList.get(0), ifPAV);
						tagAlnList.clear();
						tagIDList.clear();
						tagAlnList.add(temp);
						tagIDList.add(Integer.valueOf(tem[0]));
					}
				}
			}
			checkIfPAV(tagAlnList, tagIDList.get(0), ifPAV);
		}
		catch (Exception e) {System.out.println(e.toString());}
		System.out.println("Sam file read. Physical position of tags imported");
		return ifPAV;
	}

	private void checkIfPAV (ArrayList<String> tagAlnList, int tagId, boolean[] ifPAV) {
		int distance = 500000;
		int minMatch = 40;
		int index = tagId - 1;
		String[] temp = tagAlnList.get(0).split("\t");
		if (temp[2].startsWith("*")) {
			ifPAV[index] = true;
			return;
		}
		ifPAV[index] = true;
		for (int i = 0; i < tagAlnList.size(); i++) {
			temp = tagAlnList.get(0).split("\t");
			int chr = Integer.valueOf(temp[2]);
			int posi = Integer.valueOf(temp[3]);

			String c = temp[5]+"?";
			String[] tem = c.split("M");
			int match = 0;
			for (int j = 0; j < tem.length-1; j++) {
				tem[j] = tem[j].replaceFirst(".+\\D", "");
				match += Integer.valueOf(tem[j]);
			}
			
			if (chr == this.gChr[index][0] && Math.abs(posi - this.gPos[index][0]) < distance && match > minMatch) {
				ifPAV[index] = false;
				return;
			}
		}
	}
    
    public void checkMappingQualityMultiGroup (double pThresh, int familyGroupIndex) {
        int cnt = 0, B73Cnt = 0;
        int sameChrCnt = 0;
        int same100k = 0, same200k = 0, same1m = 0, same2m = 0, same10m = 0;
        for (int i = 0; i < this.getTagCount(); i++) {
            if (this.gNum[i] == 1) continue;
            boolean flag = false;
            for (int j = 0; j < this.gNum[i]; j++) {
                if (this.binomP[i][j] > pThresh) {
                    flag = true;
                    continue;
                }
            }
            if (flag) continue;
            cnt++;
            if (this.ifMatch[i] == true && this.ifPerfectMatch[i] && this.ifUniqueMatch[i] == true) {
                B73Cnt++;
                if (!(this.pChr[i] == this.gChr[i][familyGroupIndex])) continue;
                sameChrCnt++;
                if (Math.abs(this.pPos[i] - this.gPos[i][familyGroupIndex]) < 50000) same100k++;
                if (Math.abs(this.pPos[i] - this.gPos[i][familyGroupIndex]) < 100000) same200k++;
                if (Math.abs(this.pPos[i] - this.gPos[i][familyGroupIndex]) < 500000) same1m++;
                if (Math.abs(this.pPos[i] - this.gPos[i][familyGroupIndex]) < 1000000) same2m++;
                if (Math.abs(this.pPos[i] - this.gPos[i][familyGroupIndex]) < 5000000) same10m++;
            }
        }
        System.out.println(String.valueOf(cnt) + " multiple mapped tags are under the threshhold of " + String.valueOf(pThresh) + " at each family group");
        System.out.println(sameChrCnt + "/"+ B73Cnt+ " tags are checked");
        System.out.println((double)sameChrCnt/B73Cnt + " (persentage) tags has agreement on genetic Chr and physical Chr");
        System.out.println((double)same10m/B73Cnt + " are in 10M region");
        System.out.println((double)same2m/B73Cnt + " are in 2M region");
        System.out.println((double)same1m/B73Cnt + " are in 1M region");
        System.out.println((double)same200k/B73Cnt + " are in 200k region");
        System.out.println((double)same100k/B73Cnt + " are in 100k region");
    }
    
    public void checkMappingQuality (double pThresh, boolean ifFirstFamilyGroup, boolean ifSingleMapTag) {
        int B73Cnt = 0, cnt = 0;
        int sameChrCnt = 0;
        int same100k = 0, same200k = 0, same1m = 0, same2m = 0, same10m = 0;
        for (int i = 0; i < this.getTagCount(); i++) {
            if (ifSingleMapTag) {
                if (this.gNum[i] != 1) continue;
            }
            if (ifFirstFamilyGroup) {
                if (this.binomP[i][0] > pThresh) continue;
                cnt++;
                if (this.ifMatch[i] == true && this.ifPerfectMatch[i] && this.ifUniqueMatch[i] == true) {
                    B73Cnt++;
                    if (!(this.pChr[i] == this.gChr[i][0])) continue;
                    sameChrCnt++;
                    if (Math.abs(this.pPos[i] - this.gPos[i][0]) < 50000) same100k++;
                    if (Math.abs(this.pPos[i] - this.gPos[i][0]) < 100000) same200k++;
                    if (Math.abs(this.pPos[i] - this.gPos[i][0]) < 500000) same1m++;
                    if (Math.abs(this.pPos[i] - this.gPos[i][0]) < 1000000) same2m++;
                    if (Math.abs(this.pPos[i] - this.gPos[i][0]) < 5000000) same10m++;
                }
            }
            else {
                boolean flag = true;
                for (int j = 0; j < this.gNum[i]; j++) {
                    if (this.binomP[i][j] > pThresh) continue;
                    if (flag) {
                        cnt++;
                        flag = false;
                    }
                    if (this.ifMatch[i] == true && this.ifPerfectMatch[i] && this.ifUniqueMatch[i] == true) B73Cnt++;
                    break;
                }
                for (int j = 0; j < this.gNum[i]; j++) {
                    if (this.binomP[i][j] > pThresh) continue;
                    if (this.ifMatch[i] == true && this.ifPerfectMatch[i] && this.ifUniqueMatch[i] == true) {
                        if (!(this.pChr[i] == this.gChr[i][j])) continue;
                        sameChrCnt++;
                        if (Math.abs(this.pPos[i] - this.gPos[i][0]) < 50000) same100k++;
                        if (Math.abs(this.pPos[i] - this.gPos[i][0]) < 100000) same200k++;
                        if (Math.abs(this.pPos[i] - this.gPos[i][0]) < 500000) same1m++;
                        if (Math.abs(this.pPos[i] - this.gPos[i][0]) < 1000000) same2m++;
                        if (Math.abs(this.pPos[i] - this.gPos[i][0]) < 5000000) same10m++;
                        break;
                    }
                }
            }
        }
        System.out.println(String.valueOf(cnt) + " tags are under the threshhold of " + String.valueOf(pThresh));
        System.out.println(sameChrCnt + "/"+ B73Cnt+ " tags are checked");
        System.out.println((double)sameChrCnt/B73Cnt + " (persentage) tags has agreement on genetic Chr and physical Chr");
        System.out.println((double)same10m/B73Cnt + " are in 10M region");
        System.out.println((double)same2m/B73Cnt + " are in 2M region");
        System.out.println((double)same1m/B73Cnt + " are in 1M region");
        System.out.println((double)same200k/B73Cnt + " are in 200k region");
        System.out.println((double)same100k/B73Cnt + " are in 100k region");
    }

	@Override
	public void readBinaryFile (String infileS) {
		try {
			DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(infileS), 65536));
			int tagNum = dis.readInt();
			this.tagLengthInLong = dis.readInt();
			this.iniMatrix(tagNum);
			for (int i = 0; i < tagNum; i++) {
				for (int j = 0; j < tagLengthInLong; j++) {
					tags[j][i] = dis.readLong();
				}
				tagLength[i] = dis.readByte();
				pChr[i] = dis.readByte();
				pPos[i] = dis.readInt();
				ifMatch[i] = dis.readBoolean();
				ifPerfectMatch[i] = dis.readBoolean();
				ifUniqueMatch[i] = dis.readBoolean();
				tagTaxaCount[i] = dis.readInt();
				gNum[i] = dis.readByte();
				this.iniSubMatrix(gNum[i], i);
				for (int j = 0; j < gNum[i]; j++) {
					gChr[i][j] = dis.readByte();
					site[i][j] = dis.readInt();
					gPos[i][j] = dis.readInt();
					binomP[i][j] = dis.readDouble();
					chrSig[i][j] = dis.readInt();
					likelyhood[i][j] = dis.readDouble();
					familyNum[i][j] = dis.readByte();
					familyCode[i][j] = dis.readLong();
				}
			}
			System.out.println("Binary tagPGMap file is read from " + infileS);
		}
		catch (Exception e) {
			System.out.println(e.toString());
			System.exit(1);
		}
	}

	public void writeFastaFile (String outfileS) {
		try {
			BufferedWriter bw = new BufferedWriter (new FileWriter(outfileS), 65536);
			for (int i = 0; i < this.getTagCount(); i++) {
				String name = ">"+String.valueOf(i+1);
				bw.write(name);
				bw.newLine();
				long[] temp = new long[this.tagLengthInLong];
				for (int j = 0; j < temp.length; j++) {
					temp[j] = tags[j][i];
				}
				bw.write(BaseEncoder.getSequenceFromLong(temp).substring(0, this.getTagLength(i)));
				bw.newLine();
			}
			bw.flush();
			bw.close();
			System.out.println("Fasta file of tagPGMap output");
		}
		catch (Exception e) {
			System.out.println(e.toString());
			System.exit(1);
		}
	}

	@Override
	public void writeBinaryFile (String outfileS) {
		try {
			DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfileS), 65536));
			dos.writeInt(this.getTagCount());
			dos.writeInt(this.tagLengthInLong);
			for (int i = 0; i < this.getTagCount(); i++) {
				for (int j = 0; j < this.tagLengthInLong; j++){
					dos.writeLong(tags[j][i]);
				}
				dos.writeByte(tagLength[i]);
				dos.writeByte(pChr[i]);
				dos.writeInt(pPos[i]);
				dos.writeBoolean(ifMatch[i]);
				dos.writeBoolean(ifPerfectMatch[i]);
				dos.writeBoolean(ifUniqueMatch[i]);
				dos.writeInt(tagTaxaCount[i]);
				dos.writeByte(gNum[i]);
				for (int j = 0; j < gNum[i]; j++) {
					dos.writeByte(gChr[i][j]);
					dos.writeInt(site[i][j]);
					dos.writeInt(gPos[i][j]);
					dos.writeDouble(binomP[i][j]);
					dos.writeInt(chrSig[i][j]);
					dos.writeDouble(likelyhood[i][j]);
					dos.writeByte(familyNum[i][j]);
					dos.writeLong(familyCode[i][j]);
				}
			}
			dos.flush();
			dos.close();
			System.out.println("Binary tagPGMap file is written in " + outfileS);
		}
		catch (Exception e) {
			System.out.println(e.toString());
			System.exit(1);
		}
	}

	private void iniMatrix (int tagNum) {
		tags = new long[tagLengthInLong][tagNum];
		tagLength = new byte[tagNum];
		pChr = new byte[tagNum];
		pPos = new int[tagNum];
		ifMatch = new boolean[tagNum];
		ifPerfectMatch = new boolean[tagNum];
		ifUniqueMatch = new boolean[tagNum];
		tagTaxaCount = new int[tagNum];
		gNum = new byte[tagNum];
		gChr = new byte[tagNum][];
		site = new int[tagNum][];
		gPos = new int[tagNum][];
		binomP = new double[tagNum][];
		chrSig = new int[tagNum][];
		likelyhood = new double[tagNum][];
		familyNum = new byte[tagNum][];
		familyCode = new long[tagNum][];
		System.out.println("Matrix of tagPGMap is initialized with " + tagNum + " tags");
	}

	private void iniSubMatrix (int groupNum, int index) {
		gChr[index] = new byte[groupNum];
		site[index] = new int[groupNum];
		gPos[index] = new int[groupNum];
		binomP[index] = new double[groupNum];
		chrSig[index] = new int[groupNum];
		likelyhood[index] = new double[groupNum];
		familyNum[index] = new byte[groupNum];
		familyCode[index] = new long[groupNum];
	}

	@Override
	public void readTxtFile (String infileS) {
		int tagNum = this.getSize(infileS);
		this.iniMatrix(tagNum);
		try {
			BufferedReader br = new BufferedReader (new FileReader(infileS), 65536);
			String temp = br.readLine();
			for (int i = 0; i < tagNum; i++) {
				temp = br.readLine();
				this.readTxtRecord(temp, i);
			}
			System.out.println("Txt tagPGMap file is read from " + infileS);
		}
		catch (Exception e) {
			System.out.println(e.toString());
			System.exit(1);
		}
	}

	private void readTxtRecord (String line, int index) {
		String[] temp = line.split("\\s+");
		long[] temTag = BaseEncoder.getLongArrayFromSeq(temp[0]);
		for (int i = 0; i < this.tagLengthInLong; i++) {
			tags[i][index] = temTag[i];
		}
		tagLength[index] = Byte.valueOf(temp[1]);
		pChr[index] = Byte.valueOf(temp[2]);
		pPos[index] = Integer.valueOf(temp[3]);
		ifMatch[index] = Boolean.valueOf(temp[4]);
		ifPerfectMatch[index] = Boolean.valueOf(temp[5]);
		ifUniqueMatch[index] = Boolean.valueOf(temp[6]);
		tagTaxaCount[index] = Integer.valueOf(temp[7]);
		gNum[index] = Byte.valueOf(temp[8]);
		this.iniSubMatrix(gNum[index], index);
		for (int i = 0; i < gNum[index]; i++) {
			int j = i*8+9;
			gChr[index][i] = Byte.valueOf(temp[j]);
			site[index][i] = Integer.valueOf(temp[j+1]);
			gPos[index][i] = Integer.valueOf(temp[j+2]);
			binomP[index][i] = Double.valueOf(temp[j+3]);
			chrSig[index][i] = Integer.valueOf(temp[j+4]);
			likelyhood[index][i] = Double.valueOf(temp[j+5]);
			familyNum[index][i] = Byte.valueOf(temp[j+6]);
			familyCode[index][i] = Long.valueOf(temp[j+7]);
		}
	}

	private int getSize(String infileS) {
		int cnt = 0;
		try {
			BufferedReader br = new BufferedReader (new FileReader(infileS), 65536);
			String temp = br.readLine();
			temp = br.readLine();
			cnt++;
			String[] tem = temp.split("\\s+");
			long[] temTag = BaseEncoder.getLongArrayFromSeq(tem[0]);
			tagLengthInLong = temTag.length;
			while ((temp = br.readLine()) != null) cnt++;
		}
		catch (Exception e) {System.out.println(e.toString());}
		System.out.println(cnt + " tags int tagPGMap file " + infileS);
		return cnt;
	}

	public void writeFilteredTxtFile (String outfileS, boolean[] ifPAV) {
		try {
			BufferedWriter bw = new BufferedWriter (new FileWriter(outfileS), 65536);
			bw.write(this.getHeader());
			bw.newLine();
			for (int i = 0; i < this.getTagCount(); i++) {
				if (!ifPAV[i]) continue;
				bw.write(this.getOutputString(i));
				bw.newLine();
			}
			bw.flush();
			bw.close();
			System.out.println("Txt tagPGMap file is written in " + outfileS);
		}
		catch (Exception e) {
			System.out.println(e.toString());
			System.exit(1);
		}
	}

	@Override
	public void writeTxtFile (String outfileS) {
		try {
			BufferedWriter bw = new BufferedWriter (new FileWriter(outfileS), 65536);
			bw.write(this.getHeader());
			bw.newLine();
			for (int i = 0; i < this.getTagCount(); i++) {
				bw.write(this.getOutputString(i));
				bw.newLine();
			}
			bw.flush();
			bw.close();
			System.out.println("Txt tagPGMap file is written in " + outfileS);
		}
		catch (Exception e) {
			System.out.println(e.toString());
			System.exit(1);
		}
	}

	private String getOutputString (int index) {
		StringBuilder sb = new StringBuilder();	
		for (int i = 0; i < this.tagLengthInLong; i++) {
			sb.append(BaseEncoder.getSequenceFromLong(tags[i][index]));
		}
		sb.append("\t").append(String.valueOf(tagLength[index])).append("\t").append(String.valueOf(pChr[index])).append("\t").append(String.valueOf(pPos[index])).append("\t");
		sb.append(String.valueOf(ifMatch[index])).append("\t").append(String.valueOf(ifPerfectMatch[index])).append("\t").append(String.valueOf(ifUniqueMatch[index])).append("\t");
		sb.append(String.valueOf(this.tagTaxaCount[index])).append("\t").append(String.valueOf(this.gNum[index])).append("\t");
		for (int i = 0; i < gNum[index]; i++) {
			sb.append(String.valueOf(gChr[index][i])).append("\t");
			sb.append(String.valueOf(site[index][i])).append("\t");
			sb.append(String.valueOf(gPos[index][i])).append("\t");
			sb.append(String.valueOf(binomP[index][i])).append("\t");
			sb.append(String.valueOf(chrSig[index][i])).append("\t");
			sb.append(String.valueOf(likelyhood[index][i])).append("\t");
			sb.append(String.valueOf(familyNum[index][i])).append("\t");
			sb.append(String.valueOf(familyCode[index][i])).append("\t");
		}
		return sb.toString();
	}

	private String getHeader() {
		String header = "TestTag\tTagLength\tPChr\tPPos\tIfMatch\tIfPerfectMatch\tIfUniqueMatch\tTagTaxaCount\tFamilyGroupNum\t";
		for (int i = 0; i < familyGroupNum; i++) {
			String sufix = "-G"+String.valueOf(i+1);
			header = header + "LDChr" + sufix + "\t";
			header = header + "LDSite" + sufix + "\t";
			header = header + "LDPos" + sufix + "\t";
			header = header + "BinomP" + sufix + "\t";
			header = header + "ChrSig" + sufix + "\t";
			header = header + "Likelyhood" + sufix + "\t";
			header = header + "FamilyNum" + sufix + "\t";
			header = header + "FamilyCode" + sufix + "\t";
		}
		return header;
	}

	public void sortByGeneticPosition () {
		GenericSorting.quickSort(0, this.getTagCount(), compByGeneticPosition, this);
		System.out.println("TagPGMap is sorted by genetic position");
	}

	IntComparator compByGeneticPosition = new IntComparator() {
		public int compare(int a, int b) {
			if (gChr[a][0] == gChr[b][0]) {
				return gPos[a][0]-gPos[b][0];
			}
			else {
				return gChr[a][0]-gChr[b][0];
			}
		}
	};

	@Override
	public void sortByTag() {
        System.out.println("Position index sort begin.");
        GenericSorting.quickSort(0, this.getTagCount(), this, this);
        System.out.println("Position index sort end.");
		System.out.println("TagPGMap file is sorted by tags");
    }

	@Override
	public void swap(int index1, int index2) {
        long temp;
        for (int i = 0; i < tagLengthInLong; i++) {
            temp=tags[i][index1];
            tags[i][index1]=tags[i][index2];
            tags[i][index2]=temp;
        }
        byte temByte;
        temByte = tagLength[index1]; tagLength[index1] = tagLength[index2]; tagLength[index2] = temByte;
		int temInt;
		temInt = tagTaxaCount[index1]; tagTaxaCount[index1] = tagTaxaCount[index2]; tagTaxaCount[index2] = temInt;
		temByte = gNum[index1]; gNum[index1] = gNum[index2]; gNum[index2] = temByte;
		byte[] tempByte;
		tempByte = gChr[index1]; gChr[index1] = gChr[index2]; gChr[index2] = tempByte;
		int[] tempInt;
		tempInt = site[index1]; site[index1] = site[index2]; site[index2] = tempInt;
		tempInt = gPos[index1]; gPos[index1] = gPos[index2]; gPos[index2] = tempInt;
		double[] tempD;
		tempD = binomP[index1]; binomP[index1] = binomP[index2]; binomP[index2] = tempD;
		tempInt = chrSig[index1]; chrSig[index1] = chrSig[index2]; chrSig[index2] = tempInt;
		tempD = likelyhood[index1]; likelyhood[index1] = likelyhood[index2]; likelyhood[index2] = tempD;
		tempByte = familyNum[index1]; familyNum[index1] = familyNum[index2]; familyNum[index2] = tempByte;
		long[] tempLong;
		tempLong = familyCode[index1]; familyCode[index1] = familyCode[index2]; familyCode[index2] = tempLong;
		temByte = pChr[index1]; pChr[index1] = pChr[index2]; pChr[index2] = temByte;
		temInt = pPos[index1]; pPos[index1] = pPos[index2]; pPos[index2] = temInt;
		boolean tempB;
		tempB = ifMatch[index1]; ifMatch[index1] = ifMatch[index2]; ifMatch[index2] = tempB;
		tempB = ifPerfectMatch[index1]; ifPerfectMatch[index1] = ifPerfectMatch[index2]; ifPerfectMatch[index2] = tempB;
		tempB = ifUniqueMatch[index1]; ifUniqueMatch[index1] = ifUniqueMatch[index2]; ifUniqueMatch[index2] = tempB;
		
    }
}
