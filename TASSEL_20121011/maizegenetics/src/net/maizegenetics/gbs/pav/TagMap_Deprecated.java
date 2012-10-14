/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.gbs.pav;

import cern.colt.GenericSorting;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import net.maizegenetics.gbs.tagdist.AbstractTags;
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.genome.BaseEncoder;

/**
 *
 * @author Fei Lu
 */
public class TagMap_Deprecated extends AbstractTags {
	byte[] ldChr;
	int[] ldSiteIndex;
	int[] ldPos;
	double[] ldSig;
	byte[] blastSize;
	byte[][] blastChr;
	int[][] blastPos;
	boolean[][] ifPerfectMatch;

	public TagMap_Deprecated (String gMappingFileS, String samFileS, String tagCountFileS) {
		TagCounts tc = new TagCounts(tagCountFileS, FilePacking.Byte);
		this.tagLengthInLong = tc.getTagSizeInLong();
		this.readMappingFile(gMappingFileS, tc);
		this.readSamFile(samFileS, tc);
	}

	public void updateMappingFile (String gMappingFileS, String gMappingUpdateFileS, boolean ifOnlyPerfect) {
		try {
			int cnt = 0;
			BufferedReader br = new BufferedReader (new FileReader(gMappingFileS), 65536);
			BufferedWriter bw = new BufferedWriter (new FileWriter(gMappingUpdateFileS), 65536);
			String temp = br.readLine();
			bw.write(temp);
			bw.newLine();
			for (int i = 0; i < this.getTagCount(); i++) {
				temp = br.readLine();
				String[] tem = temp.split("\\s+");
				if (blastSize[i] == 1) {
					if (ifOnlyPerfect) {
						if (ifPerfectMatch[i][0]) {
							tem[2] = String.valueOf(blastChr[i][0]);
							tem[3] = String.valueOf(blastPos[i][0]);
						}
					}
					else {
						tem[2] = String.valueOf(blastChr[i][0]);
						tem[3] = String.valueOf(blastPos[i][0]);
					}
				}
				for (int j = 0; j < tem.length; j++) {
					bw.write(tem[j] + "\t");
				}
				bw.newLine();
			}
			bw.flush();
			bw.close();
			br.close();
		}
		catch (Exception e) {
			System.out.println(e.toString());
		}
	}

	public void checkLdAndBlast () {
		int nSingleBlast = 0;
		int nAgree = 0;
		for (int i = 0; i < this.getTagCount(); i++) {
			if (blastSize[i] == 1) {
				if (ifPerfectMatch[i][0] == false) continue;
				nSingleBlast++;
				if (ldChr[i] == blastChr[i][0]) {
					nAgree++;
				}
			}
		}
		System.out.println("The singleBlast tag num is " + nSingleBlast);
		System.out.println("Among them " + nAgree + " agree with genetic mapping on chromosome level");
		System.out.println("The rate is " + ((double)(nAgree)/(double)nSingleBlast));
	}

	private void readSamFile (String samFileS, TagCounts tc) {
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
						this.importBlast(tagAlnList, tagIDList.get(0), tc);
						tagAlnList.clear();
						tagIDList.clear();
						tagAlnList.add(temp);
						tagIDList.add(Integer.valueOf(tem[0]));
					}
				}
			}
		}
		catch (Exception e) {System.out.println(e.toString());}
	}

	private void importBlast (ArrayList<String> tagAlnList, int tagID, TagCounts tc) {
		String[] temp = tagAlnList.get(0).split("\t");
		if (temp[2].startsWith("*")) return;

		int index = this.getTagIndex(tc.getTag(tagID-1));
		if (index < 0) return;
		String[] tagAln = tagAlnList.toArray(new String[tagAlnList.size()]);
		blastSize[index] = (byte)tagAln.length;
		blastChr[index] = new byte[tagAln.length];
		blastPos[index] = new int[tagAln.length];
		ifPerfectMatch[index] = new boolean[tagAln.length];
		for (int i = 0; i < tagAln.length; i++) {
			temp = tagAln[i].split("\t");
			blastChr[index][i] = Byte.valueOf(temp[2]);
			blastPos[index][i] = Integer.valueOf(temp[3]);
			if (temp[5].startsWith("64M") && tagAln[i].contains("NM:i:0")) {
				ifPerfectMatch[index][i] = true;
			}
			else {
				ifPerfectMatch[index][i] = false;
			}
		}
	}

	private void readMappingFile (String gMappingFileS, TagCounts tc) {
		int cnt = 0;
		try {
			BufferedReader br = new BufferedReader (new FileReader(gMappingFileS), 65536);
			String temp = br.readLine();
			while ((temp = br.readLine()) != null) {cnt++;}
			br.close();
		}
		catch (Exception e) {System.out.println(e.toString());}
		tags = new long[tagLengthInLong][cnt];
		tagLength = new byte[cnt];
		ldChr = new byte[cnt];
		ldSiteIndex = new int[cnt];
		ldPos = new int[cnt];
		ldSig = new double[cnt];
		blastSize = new byte[cnt];
		blastChr = new byte[cnt][];
		blastPos = new int[cnt][];
		ifPerfectMatch = new boolean[cnt][];
		try {
			BufferedReader br = new BufferedReader (new FileReader(gMappingFileS), 65536);
			String temp = br.readLine();
			cnt = 0;
			while ((temp = br.readLine()) != null) {
				String[] tem = temp.split("\\s+");
				long[] temTag = BaseEncoder.getLongArrayFromSeq(tem[0]);
				int index = tc.getTagIndex(temTag);
				tagLength[cnt] = (byte)tc.getTagLength(index);
				for (int i = 0; i < temTag.length; i++) {
					tags[i][cnt] = temTag[i];
				}
				blastSize[cnt] = 0;
				ldChr[cnt] = Byte.valueOf(tem[5]);
				ldSiteIndex[cnt] = Integer.valueOf(tem[6]);
				ldPos[cnt] = Integer.valueOf(tem[7]);
				ldSig[cnt] = Double.valueOf(tem[8]);
				blastChr[cnt] = null;
				blastPos[cnt] = null;
				ifPerfectMatch[cnt] = null;
				cnt++;
			}
			br.close();
		}
		catch (Exception e) {System.out.println(e.toString());}
		//this.sortByTag();
	}

	@Override
    public void swap(int index1, int index2) {
        long temp;
        for (int i = 0; i < tagLengthInLong; i++) {
            temp=tags[i][index1];
            tags[i][index1]=tags[i][index2];
            tags[i][index2]=temp;
        }
        byte tl;
        tl=tagLength[index1];tagLength[index1]=tagLength[index2]; tagLength[index2]=tl;
		int tempSite;
		tempSite = ldSiteIndex[index1]; ldSiteIndex[index1] = ldSiteIndex[index2]; ldSiteIndex[index2] = tempSite;
		int tempPos;
		tempPos = ldPos[index1]; ldPos[index1] = ldPos[index2]; ldPos[index2] = tempPos;
		double tempSig;
		tempSig = ldSig[index1]; ldSig[index1] = ldSig[index2]; ldSig[index2] = tempSig;
		byte tempBlastSize;
		tempBlastSize = blastSize[index1]; blastSize[index1] = blastSize[index2]; blastSize[index2] = tempBlastSize;
		byte[] tempBlastChr;
		tempBlastChr = blastChr[index1]; blastChr[index1] = blastChr[index2]; blastChr[index2] = tempBlastChr;
		int[] tempBlastPos;
		tempBlastPos = blastPos[index1]; blastPos[index1] = blastPos[index2]; blastPos[index2] = tempBlastPos;
		boolean[] tempIfPerfectMatch;
		tempIfPerfectMatch = ifPerfectMatch[index1]; ifPerfectMatch[index1] = ifPerfectMatch[index2]; ifPerfectMatch[index2] = tempIfPerfectMatch;
    }

	public void sortByTag() {
        System.out.println("Position index sort begin.");
        GenericSorting.quickSort(0, this.getTagCount(), this, this);
        System.out.println("Position index sort end.");
    }
}
