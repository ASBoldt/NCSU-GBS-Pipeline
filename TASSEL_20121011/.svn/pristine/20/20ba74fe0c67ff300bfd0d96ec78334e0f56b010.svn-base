/*
 process SAS result table, set a filter (by p value), add marker information for the table
 According to the filtered table, exact sequences and set up a query library, and run BLAST search
 This requirs BLAST+ applications
 */

package net.maizegenetics.genome.GBS.repeat;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;

/**
 *
 * @author fl262
 */
public class processGlmSAS {
	static File genoFile = mAndSGeneticMap.mergeFile;
	static File SASResultDir = new File(setQuery.forSASDir, "SAS_result/");
	static File SASResultF = new File (SASResultDir, "GLMselect.txt");
	static File desFile = new File (SASResultDir, "filteredGLMselect.txt");
	static String blastDirS = setQuery.forSASDir + "Blast_result/";
	static File queryLib = new File(blastDirS, "query.lib");
	static String blastReS = blastDirS + "out.txt";
	float eFilter = (float) 0.0001;

	public processGlmSAS () {
		GeneticInfoArray gif = new GeneticInfoArray(genoFile);
		GlmRecSASArray grsa = new GlmRecSASArray(SASResultF, eFilter);
		setMarkerAndTrait(grsa, gif);
		grsa.writeFile(desFile);
		runBlast(grsa);
	}

	public void runBlast(GlmRecSASArray grsa) {
		setQueryLib(grsa);
		String cmd = "blastn -query " + queryLib.toString() + " -db " + GenBank2Fasta.desFile.toString() + " -out " + blastReS + " -evalue 1e-5";
		try {
			Runtime rt = Runtime.getRuntime();
			Process p = rt.exec(cmd);
			p.waitFor();
		}
		catch (Exception e) {
			System.out.println ("run blast");
		}
	}

	public void setQueryLib (GlmRecSASArray grsa) {
		File[] fa = new File(setQuery.desDir).listFiles();
		int recNum = 0;
		String[] titleArray, seqArray;
		ArrayList<String> title = new ArrayList();
		ArrayList<String> seq = new ArrayList();
		try {
			BufferedReader br = new BufferedReader(new FileReader(fa[0]), 65536);
			String temp = br.readLine();
			while (temp != null) {
				recNum++;
				title.add(temp);
				seq.add(br.readLine());
				temp = br.readLine();
			}
			br.close();
		}
		catch (Exception e) {
			System.out.println (fa[0].toString());
		}
		titleArray = title.toArray(new String[recNum]);
		seqArray = seq.toArray(new String[recNum]);
		HashSet<String> hs = new HashSet();
		for (int i = 0; i < grsa.sigMarkerNum; i++) {
			hs.add(grsa.sigMarkers[i].traitName);
		}
		String[] printList = hs.toArray(new String[hs.size()]);
		Arrays.sort(printList);
		new File(blastDirS).mkdir();
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(queryLib), 65536);
			for (int i = 0; i < recNum; i++) {
				String[] temps = titleArray[i].split("\\|");
				String temp = temps[0].replaceFirst(">", "");
				int hit  = Arrays.binarySearch(printList, temp);
				if (hit >= 0) {
					bw.write(titleArray[i]);
					bw.newLine();
					bw.write(seqArray[i]);
					bw.newLine();
				}
			}
			bw.flush();
			bw.close();
		}
		catch (Exception e) {
			System.out.println(queryLib.toString());
		}
	}

	public void setMarkerAndTrait (GlmRecSASArray grsa, GeneticInfoArray gif) {
		//Arrays.sort (gif.markers, new sortByMarkerName());
		for (int i = 0; i < grsa.sigMarkerNum; i++) {
			String temp = grsa.sigMarkers[i].markerName;
			temp = temp.replaceAll("COL", "");
			int index = Integer.valueOf(temp)-1;
			grsa.sigMarkers[i].setMarker(gif.markers[index].markerName, gif.markers[index].chrom, gif.markers[index].position);
			grsa.sigMarkers[i].setTraitName(gif.markerNum);
		}
		Arrays.sort(grsa.sigMarkers);
	}

	public static void main (String args[]) {
		processGlmSAS pgs = new processGlmSAS ();
	}

}
