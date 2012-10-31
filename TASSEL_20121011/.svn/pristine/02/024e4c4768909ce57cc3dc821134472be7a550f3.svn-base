/*
convert BLAST+ blastn result format to table format
 */

package net.maizegenetics.genome.GBS.switchgrass;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeMap;



/**
 *
 * @author fl262
 */
public class Blastn2Table {
	String[] header = {"Query_ID", "Query_Lenth", "Subject_ID", "Subject_Lenth", "Score",
						"E-Value", "Identities", "Identities_P", "Gaps", "Gap_P", "Frame", "Gap_opens",
						"Gap_postion", "Gap_seq", "Query_Origin", "Query_End", "Subject_Origin", "Subject_End"};
	
	public Blastn2Table (String infile, String outfile) {
		readAndWriteFile(new File(infile), new File(outfile));
	}

	public void readAndWriteFile(File sourceFile, File desFile) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(sourceFile), 65536);
			BufferedWriter bw = new BufferedWriter(new FileWriter(desFile), 65536);
			bw.write(getLine(header));
			bw.newLine();
			String temp = br.readLine();
			String[] item = new String[header.length];
			while (temp != null) {
				if (temp.startsWith("Query=")) {
					item[0] = temp.substring(8);
					temp = br.readLine();
					while (!temp.startsWith("Length")) {
						item[0] = item[0] + temp;
						temp = br.readLine();
					}
					item[1] = temp.substring(7);
				}
				else if (temp.startsWith(">")) {
					item[2] = temp.substring(2);
					temp = br.readLine();
					while (!temp.startsWith("Length")) {
						item[2] = item[2] + temp;
						temp = br.readLine();
					}
					item[3] = temp.substring(7);
				}
				else if (temp.startsWith(" Score")) {
					String[] temps = temp.split(",");
					item[4] = temps[0].split(" +")[3];
					item[5] = temps[1].split(" +")[3];
					temp = br.readLine();
					temps = temp.split(",");
					String[] tem = temps[0].split(" +");
					item[6] = tem[3];
					item[7] = tem[4].replaceAll("\\(|\\)|%", "");
					tem = temps[1].split(" +");
					item[8] = tem[3];
					item[9] = tem[4].replaceAll("\\(|\\)|%", "");
					temp = br.readLine();
					item[10] = temp.substring(13);
				}
				else if (temp.startsWith("Query ")) {
					StringBuilder querySB = new StringBuilder();
					StringBuilder subSB = new StringBuilder();
					String[] temps = temp.split(" +");
					if (!item[9].startsWith("0")) {
						querySB.append(temps[2]);
					}
					item[14] = temps[1];
					item[15] = temps[3];
					temp = br.readLine();
					temp = br.readLine();
					temps = temp.split(" +");
					if (!item[9].startsWith("0")) {
						subSB.append(temps[2]);
					}
					item[16] = temps[1];
					item[17] = temps[3];
					temp = br.readLine();
					temp = br.readLine();
					while (temp.startsWith("Query ")) {
						temps = temp.split(" +");
						if (!item[9].startsWith("0")) {
							querySB.append(temps[2]);
						}
						item[15] = temps[3];
						temp = br.readLine();
						temp = br.readLine();
						temps = temp.split(" +");
						if (!item[9].startsWith("0")) {
							subSB.append(temps[2]);
						}
						item[17] = temps[3];
						temp = br.readLine();
						temp = br.readLine();
					}
					String queryAln = querySB.toString();
					String subAln = subSB.toString();
					boolean ifMinus = false;
					if (item[10].startsWith("M")) {
						String mid = item[16];
						item[16] = item[17];
						item[17] = mid;
						ifMinus = true;
					}
					if (!item[9].startsWith("0")) {
						String[] gapInfo = getGapInfo(queryAln, subAln, ifMinus, Integer.valueOf(item[16]));
						item[11] = gapInfo[0];
						item[12] = gapInfo[1];
						item[13] = gapInfo[2];
					}
					bw.write(getLine(item));
					bw.newLine();
				}
				temp = br.readLine();
			}
			bw.flush();
			bw.close();
			br.close();
		}
		catch (Exception e) {
			System.out.println("Exception occurred in " + sourceFile.toString() + " " + desFile.toString());
		}	
	}
	public String[] getGapInfo (String queryAln, String subAln, boolean ifMinus, int subOrigin) {
		String[] gapInfo = new String[3];
		ArrayList<Integer> bound = new ArrayList();
		ArrayList<Integer> tempBound = new ArrayList();
		tempBound = getGapPosi(queryAln);
		if (tempBound != null) bound.addAll(tempBound);
		tempBound = getGapPosi(subAln);
		if (tempBound != null) bound.addAll(tempBound);
		Integer[] boundArray = bound.toArray(new Integer[bound.size()]);
		int[] start = new int[boundArray.length/2];
		TreeMap<Integer, Integer> boundPair = new TreeMap();
		for (int i = 0; i < start.length; i++) {
			start[i] = boundArray[2*i];
			boundPair.put(start[i], boundArray[2*i+1]);
		}
		Arrays.sort(start);
		int gapLenth = 0;
		gapInfo[0] = String.valueOf(start.length);
		StringBuilder gapSeqSB = new StringBuilder();
		StringBuilder gapPosiSB = new StringBuilder();
		for (int i = 0; i < start.length; i++) {
			String s1 = queryAln.substring(start[i], boundPair.get(start[i]));
			String s2 = subAln.substring(start[i], boundPair.get(start[i]));
			int posi;
			if (ifMinus) {
				posi = subOrigin - start[i] + gapLenth;
			}
			else {
				posi = subOrigin + start[i] - gapLenth;
			}
			gapLenth += boundPair.get(start[i]) - start[i] + 1;
			gapSeqSB.append(s1).append("/").append(s2).append(";");
			gapPosiSB.append(String.valueOf(posi)).append(";");
		}
		gapSeqSB.deleteCharAt(gapSeqSB.length()-1);
		gapPosiSB.deleteCharAt(gapPosiSB.length()-1);
		gapInfo[1] = gapPosiSB.toString();
		gapInfo[2] = gapSeqSB.toString();
		return gapInfo;
	}
	public ArrayList<Integer> getGapPosi (String str) {
		ArrayList<Integer> posi = new ArrayList();
		ArrayList<Integer> boundPosi = new ArrayList();
		for (int i = 0; i < str.length(); i++) {
			int index = str.indexOf("-", i);
			if (index != -1) {
				posi.add(index);
				i = index;
			}
		}
		if (posi.isEmpty()) return null;
		Integer[] posiArray = posi.toArray(new Integer[posi.size()]);
		boundPosi.add(posiArray[0]);
		for (int i = 1; i < posiArray.length; i++) {
			if ((posiArray[i]-posiArray[i-1]) > 1) {
				boundPosi.add(posiArray[i-1]);
				boundPosi.add(posiArray[i]);
			}
		}
		boundPosi.add(posiArray[posiArray.length-1]);
		return boundPosi;
	}
	public String getLine(String[] item) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < item.length - 1; i++) {
			sb.append(item[i]).append("\t");
		}
		sb.append(item[item.length-1]);
		return sb.toString();
	}
}
