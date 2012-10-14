/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.gbs.pav;

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
import java.io.FilenameFilter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeSet;
import net.maizegenetics.gbs.tagdist.TagCountMutable;
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByte;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByteFileMapSeq;
import net.maizegenetics.genome.BaseEncoder;

/**
 *
 * @author Fei Lu
 */
public class PAVUtils {
	public PAVUtils () {}

	//**Make small pieces of TBT from a large TBT file*/
	public TagsByTaxaByte[] sliceTBT (String mergedTBTFileS, int sliceNum) {
		try {
			DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(mergedTBTFileS), 65536*1024));
			int tagNum = dis.readInt();
			int tagLengthInLong = dis.readInt();
            int taxaNum = dis.readInt();
			String[] taxaNames = new String[taxaNum];
			for (int i = 0; i < taxaNum; i++) {
				taxaNames[i] = dis.readUTF();
			}
			int[] sliceSizes = new int[sliceNum];
			if (sliceNum > 1) {
				int leftOver = tagNum % sliceNum;
				int oneSliceSize = (tagNum - leftOver) / sliceNum;
				for (int i = 0; i < sliceNum; i++) {
					sliceSizes[i] = oneSliceSize;
				}
				for (int i = 0; i < leftOver; i++) {
					sliceSizes[i]++;
				}
			}
			else {
				sliceSizes[sliceNum-1] = tagNum;
			}
			TagsByTaxaByte[] tbts = new TagsByTaxaByte[sliceNum];
			int start = 0;
			for (int i = 0; i < sliceNum; i++) {
				long[][] tags = new long[tagLengthInLong][sliceSizes[i]];
				byte[] tagLength = new byte[sliceSizes[i]];
				byte[][] tagDist = new byte[taxaNum][sliceSizes[i]];
				for (int j = start; j < start+sliceSizes[i]; j++) {
					for (int k = 0; k < tagLengthInLong; k++) {
						tags[k][j-start] = dis.readLong();
					}
					tagLength[j-start] = dis.readByte();
					for (int k = 0; k < taxaNum; k++) {
						tagDist[k][j-start] = dis.readByte();
					}
				}
				tbts[i] = new TagsByTaxaByte(tags, tagLength, tagDist, taxaNames);
				start += sliceSizes[0];
				System.out.println("TBT piece " + (i + 1) + " with " + sliceSizes[i] + " tags is cut.");
				System.out.println("maxMemory: " + Runtime.getRuntime().maxMemory()/1024/1024 + "Mb; freeMemory: " + Runtime.getRuntime().freeMemory()/1024/1024 + "Mb");
			}
			return tbts;
		}
		catch (Exception e) {
			System.out.println(e.toString());
		}
		return null;
	}

	//**Make small pieces of TBT from a large TBT file*/
	public void sliceTBT (String sliceTBTDirS, String mergedTBTFileS, int sliceNum) {
		try {
			DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(mergedTBTFileS), 65536*1024));
			int tagNum = dis.readInt();
			int tagLengthInLong = dis.readInt();
            int taxaNum = dis.readInt();
			String[] taxaNames = new String[taxaNum];
			for (int i = 0; i < taxaNum; i++) {
				taxaNames[i] = dis.readUTF();
			}

			int leftOver = tagNum % sliceNum;
			int oneSliceSize = (tagNum - leftOver) / sliceNum;
			int[] sliceSizes = new int[sliceNum];
			for (int i = 0; i < sliceNum; i++) {
				sliceSizes[i] = oneSliceSize;
			}
			for (int i = 0; i < leftOver; i++) {
				sliceSizes[i]++;
			}

			sliceSizes[sliceNum-1] = leftOver + oneSliceSize;
			for (int i = 0; i < sliceNum; i++) {
				String subName = "sliceTBT-" + String.valueOf(i+1)+".tbt.byte";
				File outfile = new File (sliceTBTDirS, subName);
				DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfile), 65536*1024));
				dos.writeInt(sliceSizes[i]);
				dos.writeInt(tagLengthInLong);
				dos.writeInt(taxaNum);
				for (int j = 0; j < taxaNum; j++) {
					dos.writeUTF(taxaNames[j]);
				}
				for (int j = 0; j < sliceSizes[i]; j++) {
					for (int k = 0; k < tagLengthInLong; k++) {
						dos.writeLong(dis.readLong());
					}
					dos.writeByte(dis.readByte());
					for (int k = 0; k < taxaNum; k++) {
						dos.writeByte(dis.readByte());
					}
				}
				dos.flush();
				dos.close();
			}
		}
		catch (Exception e) {
			System.out.println(e.toString());
		}

	}

	//**Return chr and position of B73(Unique and perfect aligned) tag*/
	public void getB73TagPosCount (String tagPMapFileS, String tbtFileS, String outfileS) {
        TagPMap tp = new TagPMap(tagPMapFileS);
        TagsByTaxaByte tbt = new TagsByTaxaByte (tbtFileS, FilePacking.Byte);
        ArrayList<Integer> B73IndexList = new ArrayList();
        for (int i = 0; i < tbt.getTaxaCount(); i++) {
            if (tbt.getTaxaName(i).startsWith("B73")) {
                B73IndexList.add(i);
            }
        }
        int[] B73Index = new int[B73IndexList.size()];
        for (int i = 0; i < B73Index.length; i++) B73Index[i] = B73IndexList.get(i);
        ArrayList<Integer> tpIndexList = new ArrayList();
        for (int i = 0; i < tp.getTagCount(); i++) {
            if (tp.hitNum[i] > 1) continue;
			if (tp.ifPerfectMatch[i][0] == false) continue;
			if (tp.pChr[i][0] == 0) continue;
			if (tp.pChr[i][0] > 10) continue;
            tpIndexList.add(i);
        }
        Integer[] tpIndex = tpIndexList.toArray(new Integer[tpIndexList.size()]);
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(outfileS), 65536);
            bw.write("Tag\tChr\tPos\tCount");
            bw.newLine();
            for (int i = 0; i < tpIndex.length; i++) {
                bw.write(BaseEncoder.getSequenceFromLong(tp.getTag(tpIndex[i]))+"\t");
                bw.write(String.valueOf(tp.pChr[tpIndex[i]][0])+"\t");
                bw.write(String.valueOf(tp.pPos[tpIndex[i]][0])+"\t");
                int hit = tbt.getTagIndex(tp.getTag(tpIndex[i]));
                int sum = 0;
                for (int j = 0; j < B73Index.length; j++) {
                    sum+=tbt.getReadCountForTagTaxon(hit, B73Index[j]);
                }
                bw.write(String.valueOf(sum));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            System.out.println(e.toString());
        }
	}
	//**convert first n row of TBT to TagCount*/
	public void convertNrowTBTToTagCount (String inTBTFileS, String outTagCountsFileS, int nRow) {
		try {
			DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(inTBTFileS), 65536*1024));
			DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outTagCountsFileS), 65536*1024));
			int inTagNum = dis.readInt();
            int tagLengthInLong = dis.readInt();
            int taxaNum = dis.readInt();
			dos.writeInt(nRow);
			dos.writeInt(tagLengthInLong);
			for (int i = 0; i < taxaNum; i++) {
				dis.readUTF();
			}
			for (int i = 0; i < nRow; i++) {
				for (int j = 0; j < tagLengthInLong; j++) {
					dos.writeLong(dis.readLong());
                }
				dos.writeByte(dis.readByte());
				int cnt = 0;
				for (int j = 0; j < taxaNum; j++) {
					cnt+= dis.readByte();
				}
				dos.writeInt(cnt);
			}
			dos.flush();
			dos.close();
			dis.close();
		}
		catch (Exception e) {
			System.out.println(e.toString());
		}
	}

	//**convert TBT to TagCount*/
	public void convertTBTToTagCount (String inTBTFileS, String outTagCountsFileS) {
		try {
			DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(inTBTFileS), 65536*1024));
			DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outTagCountsFileS), 65536*1024));
			int inTagNum = dis.readInt();
            int tagLengthInLong = dis.readInt();
            int taxaNum = dis.readInt();
			dos.writeInt(inTagNum);
			dos.writeInt(tagLengthInLong);
			for (int i = 0; i < taxaNum; i++) {
				dis.readUTF();
			}
			for (int i = 0; i < inTagNum; i++) {
				for (int j = 0; j < tagLengthInLong; j++) {
					dos.writeLong(dis.readLong());
                }
				dos.writeByte(dis.readByte());
				int cnt = 0;
				for (int j = 0; j < taxaNum; j++) {
					cnt+= dis.readByte();
				}
				dos.writeInt(cnt);
			}
			dos.flush();
			dos.close();
			dis.close();
		}
		catch (Exception e) {
			System.out.println(e.toString());
		}
	}

	//**select tags with minimum count cutoff from a tagcount file*/
	public void shrinkTagCountFileByMinCount (File infile, File outfile, int minCount) {
		TagCounts tc = new TagCounts (infile.getAbsolutePath(), FilePacking.Byte);
		tc.writeTagCountFile(outfile.getAbsolutePath(), FilePacking.Byte, minCount);
	}

	//**In a tbt with all tags, filter out the tags which are mapped to the reference, the left is unmapped tags for genetic mapping*/
	public void makeUnmappedTBTFile (String allTBTFileS, String mappedTBTDir, String unMappedTBTFileS) {
		TagCountMutable allTC = getPseudoTagCountMutableFromTBT (allTBTFileS);
		boolean[] ifOutput = new boolean[allTC.getTagCount()];
		File[] mappedTBTFiles = new File (mappedTBTDir).listFiles();
		TagCountMutable[] mappedTCs = new TagCountMutable[mappedTBTFiles.length];
		for (int i = 0; i < mappedTCs.length; i++) {
			mappedTCs[i] = getPseudoTagCountMutableFromTBT (mappedTBTFiles[i].getAbsolutePath());
			mappedTCs[i].sort();
		}
		for (int i = 0; i < allTC.getTagCount(); i++) {
			ifOutput[i] = true;
			for (int j = 0; j < mappedTCs.length; j++) {
				int hit = mappedTCs[j].getTagIndex(allTC.getTag(i));
				if (hit > -1) {
					ifOutput[i] = false;
					break;
				}
			}
		}
		this.writeSubTBTByList(allTBTFileS, unMappedTBTFileS, ifOutput);
	}

	//**Merge genetic mapping result of multiple TBTs*/
	public void mergeGeneticMappingFiles (String mappingDirS) {
		File dir = new File(mappingDirS);
		File[] slices = dir.listFiles(new NameFilter("slice", 1));
		int[] ID = new int[slices.length];
		for (int i = 0; i < slices.length; i++) {
			String temp = slices[i].getName().replaceFirst(".*-", "").replaceFirst(".mapping.txt", "");
			ID[i] = Integer.valueOf(temp);
		}
		for (int i = 0; i < ID.length - 1; i++) {
			for (int j = i + 1; j < ID.length; j++) {
				if (ID[i] > ID[j]) {
					File tempFile = slices[i];
					int tempID = ID[i];
					slices[i] = slices[j];
					ID[i] = ID[j];
					slices[j] = tempFile;
					ID[j] = tempID;
				}
			}
		}
		try {
			String outfileS = slices[0].getParent() + "/merged.mapping.txt";
			BufferedWriter bw = new BufferedWriter (new FileWriter(outfileS), 65536);
			for (int i = 0; i < slices.length; i++) {
				BufferedReader br = new BufferedReader (new FileReader(slices[i]), 65536);
				String temp = br.readLine();
				if (i == 0) {
					bw.write(temp);
					bw.newLine();
				}
				while ((temp = br.readLine()) != null) {
					bw.write(temp);
					bw.newLine();
				}
				br.close();
			}
			bw.flush();
			bw.close();
		}
		catch (Exception e) {
			System.out.println(e.toString());
		}
		

	}
	
	private class NameFilter implements FilenameFilter {
		String filterS;
		int mode;
		NameFilter (String filterS, int mode) {
			this.filterS = filterS;
			this.mode = mode;
		}
		public boolean accept(File dir, String name) {
			if (mode == 0) {
				return name.startsWith(filterS);
			}
			else if (mode == 1) {
				return name.contains(filterS);
			}
			else if (mode == 2) {
				return name.endsWith(filterS);
			}
			return false;
		}

	}
	//**To get a tag list, generate a tagCountMutable from TBT, the count of each tag is 1*/
	private TagCountMutable getPseudoTagCountMutableFromTBT (String TBTFileS) {
		TagCountMutable tcm = null;
		try {
			DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(TBTFileS), 65536*1024));
			int tagNum = dis.readInt();
			int tagLengthInLong = dis.readInt();
            int taxaNum = dis.readInt();
			for (int i = 0; i < taxaNum; i++) {
				dis.readUTF();
			}
			tcm = new TagCountMutable (tagLengthInLong, tagNum);
			for (int i = 0; i < tagNum; i++) {
				long[] tag = new long[tagLengthInLong];
				for (int j = 0; j < tagLengthInLong; j++) {
					tag[j] = dis.readLong();
                }
				byte tagLength = dis.readByte();
				for (int j = 0; j < taxaNum; j++) {
					dis.readByte();
				}
				tcm.addReadCount(tag, tagLength, 1);
			}
		}
		catch (Exception e) {
			System.out.println(e.toString());
		}
		System.out.println("PseudoTagCountMutable is produced from " + TBTFileS);
		return tcm;
	}

	//**output subset of input TBT file according to the boolean array of each tag*/
	private void writeSubTBTByList (String inputTBTFileS, String outputTBTFileS, boolean[] ifOutput) {
		int outTagNum = 0;
		for (int i = 0; i < ifOutput.length; i++) {
			if (ifOutput[i]) outTagNum++;
		}
		try {
			DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(inputTBTFileS), 65536*1024));
			DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputTBTFileS), 65536*1024));
			int inTagNum = dis.readInt();
            int tagLengthInLong = dis.readInt();
            int taxaNum = dis.readInt();
			dos.writeInt(outTagNum);
			dos.writeInt(tagLengthInLong);
			dos.writeInt(taxaNum);
			for (int i = 0; i < taxaNum; i++) {
				dos.writeUTF(dis.readUTF());
			}
			for (int i = 0; i < inTagNum; i++) {
				if (!ifOutput[i]) continue;
				for (int j = 0; j < tagLengthInLong; j++) {
					dos.writeLong(dis.readLong());
                }
				dos.writeByte(dis.readByte());
				for (int j = 0; j < taxaNum; j++) {
					dos.writeByte(dis.readByte());
				}
			}
			dos.flush();
			dos.close();
			dis.close();
		}
		catch (Exception e) {

		}
	}

	//**output subset of first n row of input TBT file*/
	public void writeSubTBTByNumber (String inputTBTFileS, String outputTBTFileS, int outTagNum) {
		try {
			DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(inputTBTFileS), 65536*1024));
			DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputTBTFileS), 65536*1024));
			int inTagNum = dis.readInt();
            int tagLengthInLong = dis.readInt();
            int taxaNum = dis.readInt();
			dos.writeInt(outTagNum);
			dos.writeInt(tagLengthInLong);
			dos.writeInt(taxaNum);
			for (int i = 0; i < taxaNum; i++) {
				dos.writeUTF(dis.readUTF());
			}
			for (int i = 0; i < outTagNum; i++) {
				for (int j = 0; j < tagLengthInLong; j++) {
					dos.writeLong(dis.readLong());
                }
				dos.writeByte(dis.readByte());
				for (int j = 0; j < taxaNum; j++) {
					dos.writeByte(dis.readByte());
				}
			}
			dos.flush();
			dos.close();
			dis.close();
		}
		catch (Exception e) {
			System.out.println(e.toString());
		}
	}

	//Output subset of TBT with txt format, this is for Ali//
	public void writeSubTBTTxtByTagCount (String inputTBTFileS, String outputTBTFileS, TagCounts tc) {
		try {
			DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(inputTBTFileS), 65536*1024));
			BufferedWriter bw = new BufferedWriter (new FileWriter(outputTBTFileS), 65536);
			int inTagNum = dis.readInt();
            int tagLengthInLong = dis.readInt();
            int taxaNum = dis.readInt();
			bw.write(String.valueOf(tc.getTagCount())+"\t"+String.valueOf(tagLengthInLong)+"\t"+String.valueOf(taxaNum));
			bw.newLine();
			for (int i = 0; i < taxaNum; i++) {
				bw.write(dis.readUTF()+"\t");
			}
			bw.newLine();
			int cnt = 0;
			for (int i = 0; i < inTagNum; i++) {
				long[] tempTag = new long[tagLengthInLong];
				for (int j = 0; j < tagLengthInLong; j++) {
					tempTag[j] = dis.readLong();
				}
				int index = tc.getTagIndex(tempTag);
				if (index < 0) {
					dis.readByte();
					for (int j = 0; j < taxaNum; j++) {
						dis.readByte();
					}
				}
				else {
					bw.write(BaseEncoder.getSequenceFromLong(tempTag)+"\t"+String.valueOf(dis.readByte())+"\t");

					for (int j = 0; j < taxaNum; j++) {
						bw.write(String.valueOf(dis.readByte())+"\t");
					}
					bw.newLine();
					cnt++;
					if (cnt == tc.getTagCount()) break;
					if (cnt % 100 == 0) System.out.println(cnt);
				}
			}
			bw.flush();
			bw.close();
			dis.close();
		}
		catch (Exception e) {
			System.out.println(e.toString());
		}
	}

	//**I have no idea why I made this*/
	public void writeSubTBTTxtByTagPGMap (String inputTBTFileS, String outputTBTFileS, TagPGMap tpg) {
		try {
			DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(inputTBTFileS), 65536*1024));
			BufferedWriter bw = new BufferedWriter (new FileWriter(outputTBTFileS), 65536);
			int inTagNum = dis.readInt();
            int tagLengthInLong = dis.readInt();
            int taxaNum = dis.readInt();
			bw.write(String.valueOf(tpg.getTagCount())+"\t"+String.valueOf(tagLengthInLong)+"\t"+String.valueOf(taxaNum));
			bw.newLine();
			for (int i = 0; i < taxaNum; i++) {
				bw.write(dis.readUTF()+"\t");
			}
			bw.newLine();
			int cnt = 0;
			for (int i = 0; i < inTagNum; i++) {
				long[] tempTag = new long[tagLengthInLong];
				for (int j = 0; j < tagLengthInLong; j++) {
					tempTag[j] = dis.readLong();
				}
				int index = tpg.getTagIndex(tempTag);
				if (index < 0) {
					dis.readByte();
					for (int j = 0; j < taxaNum; j++) {
						dis.readByte();
					}
				}
				else {
					bw.write(BaseEncoder.getSequenceFromLong(tempTag)+"\t"+String.valueOf(dis.readByte())+"\t");

					for (int j = 0; j < taxaNum; j++) {
						bw.write(String.valueOf(dis.readByte())+"\t");
					}
					bw.newLine();
					cnt++;
					if (cnt == tpg.getTagCount()) break;
					if (cnt % 100 == 0) System.out.println(cnt);
				}
			}
			bw.flush();
			bw.close();
			dis.close();
		}
		catch (Exception e) {
			System.out.println(e.toString());
		}
	}

	//**make a subset of TBT by tags in tagCount file*/
	public void writeSubTBTByTagCount (String inputTBTFileS, String outputTBTFileS, TagCounts tc) {
		try {
			DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(inputTBTFileS), 65536*1024));
			DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputTBTFileS), 65536*1024));
			int inTagNum = dis.readInt();
            int tagLengthInLong = dis.readInt();
            int taxaNum = dis.readInt();
			dos.writeInt(tc.getTagCount());
			dos.writeInt(tagLengthInLong);
			dos.writeInt(taxaNum);
			for (int i = 0; i < taxaNum; i++) {
				dos.writeUTF(dis.readUTF());
			}
			int cnt = 0;
			for (int i = 0; i < inTagNum; i++) {
				long[] tempTag = new long[tagLengthInLong];
				for (int j = 0; j < tagLengthInLong; j++) {
					tempTag[j] = dis.readLong();
				}
				int index = tc.getTagIndex(tempTag);
				if (index < 0) {
					dis.readByte();
					for (int j = 0; j < taxaNum; j++) {
						dis.readByte();
					}
				}
				else {
					for (int j = 0; j < tagLengthInLong; j++) {
						dos.writeLong(tempTag[j]);
					}
					dos.writeByte(dis.readByte());
					for (int j = 0; j < taxaNum; j++) {
						dos.writeByte(dis.readByte());
					}
					cnt++;
					if (cnt == tc.getTagCount()) break;
				}
			}
			dos.flush();
			dos.close();
			dis.close();
		}
		catch (Exception e) {
			System.out.println(e.toString());
		}
	}

	//**Make fasta file of TagCounts, PolyA are removed*/
	public void writeFastaFromTagCount (String tagCountFileS, String fastaFileS) {
		TagCounts tc = new TagCounts (tagCountFileS, FilePacking.Byte);
		try {
			BufferedWriter bw = new BufferedWriter (new FileWriter(fastaFileS), 65536);
			for (int i = 0; i < tc.getTagCount(); i++) {
				String seq = BaseEncoder.getSequenceFromLong(tc.getTag(i)).substring(0, tc.getTagLength(i));
				bw.write(">"+String.valueOf(i+1));
				bw.newLine();
				bw.write(seq);
				bw.newLine();
			}
			bw.flush();
			bw.close();
		}
		catch (Exception e) {
			System.out.println(e.toString());
		}
	}

	public void makeSubTBTByTaxaTagCount (String inputTBTFileS, String tagCountFileS, String taxaListFileS, String outputTBTFileS) {
		String[] interestName = null;
		int[] totalCount = null;
		try {
			ArrayList<String> nameList = new ArrayList();
			BufferedReader br = new BufferedReader (new FileReader(taxaListFileS), 65536);
			String temp;
			while ((temp = br.readLine()) != null) {
				nameList.add(temp);
			}
			interestName = nameList.toArray(new String[nameList.size()]);
		}
		catch (Exception e) {
			System.out.println(e.toString());
		}
		Arrays.sort(interestName);
		TagCounts tc = new TagCounts(tagCountFileS, FilePacking.Byte);
		try {
			DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(inputTBTFileS), 65536*1024));
			DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputTBTFileS), 65536*1024));
			int inTagNum = dis.readInt();
            int tagLengthInLong = dis.readInt();
            int taxaNum = dis.readInt();
			String[] taxaName = new String[taxaNum];
			for (int i = 0; i < taxaNum; i++) {
				taxaName[i] = dis.readUTF();
			}
			int[] index = this.redirect(interestName, taxaName);
			ArrayList<String> nameList = new ArrayList();
			ArrayList<Integer> indexList = new ArrayList();
			for (int i = 0; i < index.length; i++) {
				if (index[i] == -1) continue;
				nameList.add(interestName[i]);
				indexList.add(index[i]);
			}
			interestName = new String[nameList.size()];
			index = new int[nameList.size()];
			for (int i = 0; i <  index.length; i++) {
				interestName[i] = nameList.get(i);
				index[i] = indexList.get(i);
			}
			dos.writeInt(tc.getTagCount());
			dos.writeInt(tagLengthInLong);
			dos.writeInt(interestName.length);
			for (int i = 0; i < interestName.length; i++) {
				dos.writeUTF(interestName[i]);
			}
			byte[] tagDis = new byte[taxaNum];
			totalCount = new int[interestName.length];
			for (int i = 0; i < totalCount.length; i++) totalCount[i] = 0;
			int cnt = 0;
			for (int i = 0; i < inTagNum; i++) {
				long[] tempTag = new long[tagLengthInLong];
				for (int j = 0; j <  tagLengthInLong; j++) {
					tempTag[j] = dis.readLong();
				}
				int ind = tc.getTagIndex(tempTag);
				if (ind < 0) {
					dis.readByte();
					for (int j = 0; j < taxaNum; j++) {
						tagDis[j] = dis.readByte();
					}
				}
				else {
					for (int j = 0; j < tagLengthInLong; j++) {
						dos.writeLong(tempTag[j]);
					}
					dos.writeByte(dis.readByte());
					for (int j = 0; j < taxaNum; j++) {
						tagDis[j] = dis.readByte();
					}
					for (int j = 0; j < index.length; j++) {
						dos.writeByte(tagDis[index[j]]);
					}
					cnt++;
				}
				for (int j = 0; j < index.length; j++) {
					totalCount[j] += tagDis[index[j]];
				}
			}
			dos.flush();
			dos.close();
			dis.close();
			System.out.println("Subset TBT of "+String.valueOf(cnt) + " tags " + interestName.length + " taxa output");

		}
		catch (Exception e) {
			System.out.println(e.toString());
		}
		for (int i = 0; i < interestName.length; i++) {
			System.out.println(interestName[i]+"\t"+totalCount[i]);
		}
	}

	private int[] redirect (String[] source, String[] target) {
		int[] index = new int[source.length];
		int cnt = 0;
		for (int i = 0; i < index.length; i++) {
			index[i] = -1;
			for (int j = 0; j < target.length; j++) {
				if (source[i].equals(target[j])) {
					index[i] = j;
					break;
				}
			}
			if (index[i] == -1) cnt++;
		}
		System.out.println(String.valueOf(cnt) + " items are not redirected to the target");
		return index;
	}

	public void makeCNVTagCountFile (String tagPMapFileS, String tagCountFileS, String cnvTagCountFileS) {
		TagPMap tp  = new TagPMap(tagPMapFileS);
		TagCounts tc = new TagCounts(tagCountFileS, FilePacking.Byte);
		int cnt = 0;
		for (int i = 0; i < tp.getTagCount(); i++) {
			if (tp.hitNum[i] != 1) continue;
			if (tp.pChr[i][0] < 1 || tp.pChr[i][0] > 10) continue;
			cnt++;
		}
		TagCounts cnvTc = new TagCounts(tc.getTagSizeInLong(), cnt);
		cnt = 0;
		for (int i = 0; i < tp.getTagCount(); i++) {
			if (tp.hitNum[i] != 1) continue;
			if (tp.pChr[i][0] < 1 || tp.pChr[i][0] > 10) continue;
			long[] tempTag = tp.getTag(i);
			int index = tc.getTagIndex(tempTag);
			cnvTc.setTag(tempTag, (byte)tc.getTagLength(index), tc.getReadCount(index), cnt);
			cnt++;
		}
		cnvTc.writeTagCountFile(cnvTagCountFileS, FilePacking.Byte, 0);
		System.out.println(String.valueOf(cnvTc.getTagCount()) + " unique tags output in CNV tagCount File");
	}

	public void makeCNVPosFile (String inputTBTFileS, String cnvPosFileS, String tagPMapFileS, String totalCountFileS) {
		TagPMap pmap = new TagPMap(tagPMapFileS);
		String[] name = null;
		int[] total = null;
		try {
			ArrayList<String> lineList = new ArrayList();
			BufferedReader br = new BufferedReader (new FileReader(totalCountFileS), 65536);
			br.readLine();
			String temp;
			while ((temp = br.readLine()) != null) {
				lineList.add(temp);
			}
			name = new String[lineList.size()];
			total = new int[lineList.size()];
			for (int i = 0; i < lineList.size(); i++) {
				String[] tem = lineList.get(i).split("\t");
				name[i] = tem[0];
				total[i] = Integer.valueOf(tem[1]);
			}
		}
		catch (Exception e) {
			System.out.println(e.toString());
		}
		try {
			DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(inputTBTFileS), 65536*1024));
			BufferedWriter bw = new BufferedWriter (new FileWriter(cnvPosFileS), 65536);
			int tagNum = dis.readInt();
            int tagLengthInLong = dis.readInt();
            int sampleNum = dis.readInt();
			String[] sampleName = new String[sampleNum];
			for (int i = 0; i < sampleNum; i++) {
				sampleName[i] = dis.readUTF();
			}
/*
			TreeSet<String> nameSet = new TreeSet();
			for (int i = 0; i < sampleName.length; i++) {
				String[] temp = sampleName[i].split(":");
				nameSet.add(temp[0]);
			}
			String[] name = nameSet.toArray(new String[nameSet.size()]);

 *
 */
			ArrayList<Integer>[] indexList = new ArrayList[name.length];
			Integer[][] index = new Integer[name.length][];
			for (int i = 0; i < name.length; i++) {
				indexList[i] = new ArrayList();
				for (int j = 0; j < sampleName.length; j++) {
					String[] temp = sampleName[j].split(":");
					if (temp[0].equals(name[i])) indexList[i].add(j);
				}
				index[i] = indexList[i].toArray(new Integer[indexList[i].size()]);
			}
			bw.write("Tags\tLength\tChr\tPos\t");
			for (int i = 0; i < name.length; i++) {
				bw.write(name[i]+"\t");
			}
			bw.newLine();
			for (int i = 0; i < tagNum; i++) {
				long[] tempTag = new long[tagLengthInLong];
				for (int j = 0; j < tagLengthInLong; j++) {
					tempTag[j] = dis.readLong();
				}
				byte length = dis.readByte();
				int ind = pmap.getTagIndex(tempTag);
				if (pmap.pChr[ind][0] < 1 || pmap.pChr[ind][0] > 10) continue;
				byte[] dist = new byte[sampleName.length];
				int cnt = 0;
				for (int j = 0; j < sampleName.length; j++) {
					dist[j] = dis.readByte();
					cnt+=dist[j];
				}
				if (cnt==0) continue;
				bw.write(BaseEncoder.getSequenceFromLong(tempTag)+"\t"+String.valueOf(length)+"\t");
				bw.write(String.valueOf(pmap.pChr[ind][0])+"\t"+String.valueOf(pmap.pPos[ind][0])+"\t");
				for (int j = 0; j < name.length; j++) {
					cnt = 0;
					for (int k = 0; k < index[j].length; k++) {
						cnt+=dist[index[j][k]];
					}
					bw.write(String.valueOf((double)cnt/total[j])+"\t");
				}
				bw.newLine();
				bw.flush();
			}
			bw.flush();
			bw.close();
			dis.close();
			System.out.println("CNV tagPos file output");
		}
		catch (Exception e) {
			System.out.println(e.toString());
		}
	}

	

	public void correct (String inputTBTFileS, String outputTBTFileS) {
		try {
			DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(inputTBTFileS), 65536*1024));
			DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputTBTFileS), 65536*1024));
			int tagNum = dis.readInt();
			tagNum = 10664067;
            int tagLengthInLong = dis.readInt();
            int taxaNum = dis.readInt();
			dos.writeInt(tagNum);
			dos.writeInt(tagLengthInLong);
			dos.writeInt(taxaNum);
			String[] taxaName = new String[taxaNum];
			for (int i = 0; i < taxaNum; i++) {
				taxaName[i] = dis.readUTF();
				dos.writeUTF(taxaName[i]);
			}
			for (int i = 0; i < tagNum; i++) {
				long[] tempTag = new long[2];
				for (int j = 0; j < tagLengthInLong; j++) {
					tempTag[j] = dis.readLong();
					dos.writeLong(tempTag[j]);
				}
				byte temp = dis.readByte();
				dos.writeByte(temp);
				for (int j = 0; j < taxaNum; j++) {
					temp = dis.readByte();
					dos.writeByte(temp);
				}
				if (i%1000 ==0) System.out.println(i);
			}
			dos.flush();
			dos.close();
			dis.close();
		}
		catch (Exception e) {

		}
	}

}
