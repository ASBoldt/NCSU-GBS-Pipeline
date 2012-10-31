/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS.plasmid;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
/**
 *
 * @author fl262
 */
public class tagCountArray {
	tagCount[] tagCs;
	int totalTag;

	public tagCountArray (String inFileS) {
		readTagCount(inFileS);
	}
	public tagCountArray (File inFile) {
		readTagCount(inFile);
	}
	public tagCountArray (int totalTag) {
		this.totalTag = totalTag;
		tagCs = new tagCount[totalTag];
	}
	public void writeCombinedTagCount (String outFileS, int currentRow) {
		File outFile = new File(outFileS);
		try {
			DataOutputStream dos = new DataOutputStream(new BufferedOutputStream (new FileOutputStream(outFile), 65536));
			for (int i = 0; i < currentRow; i++) {
				dos.write(tagCs[i].tag);
				dos.writeInt(tagCs[i].count);
			}
			dos.flush();
			dos.close();
		}
		catch (Exception e) {
			System.out.println("Error occurred while writing " + outFileS + " " + e.toString());
		}
	}
	public void writeTagCount (String outFileS, int minCount) {
		File outFile = new File(outFileS);
		try {
			DataOutputStream dos = new DataOutputStream(new BufferedOutputStream (new FileOutputStream(outFile), 65536));
			for (int i = 0; i < totalTag; i++) {
				if (tagCs[i].count >= minCount) {
					dos.write(tagCs[i].tag);
					dos.writeInt(tagCs[i].count);
				}
			}
			dos.flush();
			dos.close();
		}
		catch (Exception e) {
			System.out.println("Error occurred while writing " + outFileS + " " + e.toString());
		}
	}
	public void writeTagCount (String outFileS) {
		File outFile = new File(outFileS);
		try {
			DataOutputStream dos = new DataOutputStream(new BufferedOutputStream (new FileOutputStream(outFile), 65536));
			for (int i = 0; i < totalTag; i++) {
				dos.write(tagCs[i].tag);
				dos.writeInt(tagCs[i].count);
			}
			dos.flush();
			dos.close();
		}
		catch (Exception e) {
			System.out.println("Error occurred while writing " + outFileS + " " + e.toString());
		}
	}
	public void readTagCount (String inFileS) {
		readTagCount(new File(inFileS));
	}
	public void readTagCount (File inFile) {
		totalTag = (int) (inFile.length() / (baseEncoder.byteLen + 4));
		tagCs = new tagCount[totalTag];
		try {
			DataInputStream dis = new DataInputStream(new BufferedInputStream (new FileInputStream(inFile),65536));
			for (int i = 0; i < totalTag; i++) {
				byte[] tag = new byte[baseEncoder.byteLen];
				dis.read(tag, 0, baseEncoder.byteLen);
				tagCs[i] = new tagCount(tag, dis.readInt());
			}
			dis.close();
		}
		catch (Exception e) {
			System.out.println("Error occurred while reading " + inFile.toString() + " " +e.toString());
		}
		Arrays.sort(tagCs);
	}
	public int collapseCounts () {
		//tagCs need to be sorted by tag
		int sameCount = 0;
		ArrayList<Integer> index = new ArrayList();
		index.add(0);
		for (int i = 1; i < totalTag; i++) {
			if (tagCs[i].byteSeqEqualwith(tagCs[i-1])) {
				tagCs[i].count += tagCs[i-1].count;
				tagCs[i-1].count = 0;
				index.set(index.size()-1,i);
				sameCount++;
			}
			else {
				index.add(i);
			}
		}
		Integer[] indexArray = index.toArray(new Integer[index.size()]);
		tagCount[] tempTagCs = new tagCount[indexArray.length];
		for (int i = 0; i < tempTagCs.length; i++) {
			tempTagCs[i] = this.tagCs[indexArray[i]];
		}
		this.tagCs = tempTagCs;
		this.totalTag = tempTagCs.length;
		return sameCount;
	}

	public int mergeTagCount (String inFileS, int currentRow) {
		return mergeTagCount(new File (inFileS), currentRow);
	}
	public int mergeTagCount (File inFile, int currentRow) {
		return mergeTagCount(new tagCountArray(inFile), currentRow);
	}
	public int mergeTagCount (tagCountArray o, int currentRow) {
		int end = currentRow;
		for (int i = 0; i < o.totalTag; i++) {
			int hit = Arrays.binarySearch(tagCs, 0, end, o.tagCs[i]);
			if (hit < 0) {
				tagCs[currentRow] = (o.tagCs[i]);
				currentRow++;
			}
			else {
				tagCs[hit].count += o.tagCs[i].count;
			}
		}
		o = null;
		System.out.println("currentRow is " + currentRow);
		System.out.println("maxRow     is " + this.totalTag);
		this.sortTag(0, currentRow);
		System.gc();
		return currentRow;
	}
	public void sortCount () {
		Arrays.sort(tagCs, new sortByCount());
	}
	public void sortTag (int start, int end) {
		Arrays.sort(tagCs, start, end);
	}
	public void sortTag () {
		Arrays.sort(tagCs);
	}
	public class tagCount implements Comparable <tagCount> {
		byte[] tag;
		int count;
		public tagCount (byte[] tag, int count) {
			this.tag = tag;
			this.count = count;
		}
		public boolean byteSeqEqualwith (tagCount hTagC) {
			for (int i = 0; i < tag.length; i++) {
				if (tag[i] != hTagC.tag[i]) {
					return false;
				}
			}
			return true;
		}
		public void deepEquals(tagCount o) {
			for (int i = 0; i < baseEncoder.byteLen; i++) {
				tag[i] = o.tag[i];
			}
			count = o.count;
		}
		public int compareTo(tagCount o) {
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
	class sortByCount implements Comparator <tagCount> {
		public int compare(tagCount o1, tagCount o2) {
			return o1.count - o2.count;
		}
	}
}
