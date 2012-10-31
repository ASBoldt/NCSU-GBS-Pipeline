/*
 * Integrate 2 genetic maps into 1, because some line present in one genetic map may absent in other one, for which "MM" is used to denote genotypes in absent lines
 * interval is used to reduce the marker number, as excess markers are not essential for IBM
 */

package net.maizegenetics.genome.GBS.repeat;

import java.io.*;
import java.util.Arrays;
import java.util.HashSet;
import java.util.TreeMap;
/**
 *
 * @author fl262
 */
public class mAndSGeneticMap {
	static File sourceGenoFile55k = new File("E:/Database/InfoFile/SNP55hapmap05162010_del191.txt");
	static File sourceGenoFile644 = new File("E:/Database/InfoFile/IBM_644_NAM_SNP_genos_283_RILs_AGPv2_hapmap.txt");
	static String mergeFileD = setQuery.forSASDir + "marker_selection/";
	static File mergeFile = new File (mergeFileD, "644plus55k.txt");
	int interval = 10;
	public mAndSGeneticMap () {
		GeneticInfoArray gif644 = new GeneticInfoArray(sourceGenoFile644);
		GeneticInfoArray gif55k = new GeneticInfoArray(sourceGenoFile55k);
		new File(mergeFileD).mkdir();
		writeMergeFile (gif644, gif55k, mergeFile, interval);
	}

	public void writeMergeFile (GeneticInfoArray gif1, GeneticInfoArray gif2, File mergeFile, int interval) {
		String[] tName1 = gif1.taxaName.clone();
		String[] tName2 = gif2.taxaName.clone();
		TreeMap<String, Integer> treeIndex1 = new TreeMap();
		TreeMap<String, Integer> treeIndex2 = new TreeMap();
		for  (int i = 0; i < tName1.length; i++) {
			treeIndex1.put (tName1[i], i);
		}
		for (int i = 0; i < tName2.length; i++) {
			treeIndex2.put (tName2[i], i);
		}
		HashSet<String> tName =  new HashSet<String> (Arrays.asList(tName1));
		for (int i = 0; i < gif2.taxaNum; i++) {
			tName.add(tName2[i]);
		}
		String[] taxaName = tName.toArray(new String[tName.size()]);
		Arrays.sort(tName1);
		Arrays.sort(tName2);
		Arrays.sort(taxaName);
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(mergeFile), 65536);
			StringBuilder sb = new StringBuilder(gif1.title);
			for (int i = 0; i < taxaName.length; i++) {
				sb.append(taxaName[i]).append("\t");
			}
			bw.write(sb.toString());
			bw.newLine();
			for (int i = 0; i < gif1.markerNum; i++) {
				sb = new StringBuilder();
				for (int j = 0; j < 11; j++) {
					sb = sb.append(gif1.markers[i].bgInfo[j]).append("\t");
				}
				for (int j = 0; j < taxaName.length; j++ ) {
					int hit = Arrays.binarySearch(tName1, taxaName[j]);
					if (hit < 0) {
						sb.append("MM").append("\t");
					}
					else {
						sb.append(gif1.markers[i].genotype[treeIndex1.get(tName1[hit])]).append("\t");
					}
				}
				bw.write(sb.toString());
				bw.newLine();
			}
			for (int i = 0; i < gif2.markerNum; i+= interval) {
				sb = new StringBuilder();
				for (int j = 0; j < 11; j++) {
					sb = sb.append(gif2.markers[i].bgInfo[j]).append("\t");
				}
				for (int j = 0; j < taxaName.length; j++ ) {
					int hit = Arrays.binarySearch(tName2, taxaName[j]);
					if (hit < 0) {
						sb.append("MM").append("\t");
					}
					else {
						sb.append(gif2.markers[i].genotype[treeIndex2.get(tName2[hit])]).append("\t");
					}
				}
				bw.write(sb.toString());
				bw.newLine();
			}
			bw.flush();
			bw.close();
		}
		catch (Exception e) {
			System.out.println (mergeFile.toString());
		}

	}

	public static void main (String args[]) {
		mAndSGeneticMap mgm = new mAndSGeneticMap();
	}

}
