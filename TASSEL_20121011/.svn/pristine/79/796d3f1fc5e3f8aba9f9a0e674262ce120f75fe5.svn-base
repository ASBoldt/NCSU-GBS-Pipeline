/*
 This is for tassal Glm result, read Tassel output table, write filtered result with marker information
 Because of its low throughput and inflexibility, Tassel is not recommanded for Glm analysis
 */

package net.maizegenetics.genome.GBS.repeat;

import java.io.*;
import java.util.ArrayList;

/**
 *
 * @author fl262
 */
public class GlmRecArrayTassel {
	glmRecTassel[] sigMarkers;
	int sigMarkerNum;
	int traitNum;
	public GlmRecArrayTassel (File sourceGlmFile, float eFilter) {
		readFile (sourceGlmFile, eFilter);
	}

	public void writeFile (File desGlmFile) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(desGlmFile), 65536);
			StringBuilder title = new StringBuilder("TraitName").append("\tMarkerName\tChromosome\tPosition\tP_Marker\n");
			bw.write(title.toString());
			for (int i = 0; i < sigMarkerNum; i++) {
				StringBuilder rec = new StringBuilder(sigMarkers[i].traitName).append("\t");
				rec.append(sigMarkers[i].markerName).append("\t");
				rec.append(sigMarkers[i].chrom).append("\t");
				rec.append(sigMarkers[i].position).append("\t");
				rec.append(sigMarkers[i].pMarker).append("\n");
				bw.write(rec.toString());
			}
			bw.flush();
			bw.close();
		}
		catch (Exception e) {
			System.out.println(desGlmFile.toString());
		}
	}
	public void writeFile (String desGlmFileS) {
		writeFile (new File(desGlmFileS));
	}

	public void readFile (File sourceGlmFile, float eFilter) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(sourceGlmFile), 65536);
			String fileTitle = br.readLine();
			ArrayList<String> filterLine = new ArrayList<String>();
			String temp = br.readLine();
			while (temp != null) {
				String[] temps = temp.split("\\t");
				if (Float.valueOf(temps[7]) < eFilter) {
					filterLine.add(temp);
				}
				temp = br.readLine();
			}
			br.close();
			String[] filterLineArray = filterLine.toArray(new String[filterLine.size()]);
			sigMarkerNum = filterLineArray.length;
			sigMarkers = new glmRecTassel[sigMarkerNum];
			String tempTraitName = "";
			for (int i = 0; i < sigMarkerNum; i++) {
				String[] temps = filterLineArray[i].split("\\t");
				sigMarkers[i] = new glmRecTassel(temps[0], temps[1], Float.valueOf(temps[7]));
				if (!tempTraitName.equals(temps[0])) {
					traitNum++;
				}
				tempTraitName = temps[0];
			}
		} catch (Exception e) {
			System.out.println(sourceGlmFile);
		}
	}

}
class glmRecTassel implements Comparable <glmRecTassel> {
	String traitName;
	String markerName;
	byte chrom;
	int position;
	float pMarker;
	glmRecTassel (String traitName, String markerName, byte chrom, int position, float pMarker) {
		this.traitName = traitName;
		this.markerName = markerName;
		this.chrom = chrom;
		this.position = position;
		this.pMarker = pMarker;
	}
	glmRecTassel (String traitName, String markerName, float pMarker) {
		this.traitName = traitName;
		this.markerName = markerName;
		this.pMarker = pMarker;
	}
	public void setChAndPosi (byte chrom, int position) {
		this.chrom = chrom;
		this.position = position;
	}

	public int compareTo(glmRecTassel o) {
		if (traitName.equals(o.traitName)){
			if (pMarker < o.pMarker) return -1;
			else if (pMarker > o.pMarker) return 1;
			else {
				if (chrom < o.chrom) return -1;
				else if (chrom > o.chrom) return 1;
				else return position - o.position;
			}
		}
		else {
			return traitName.compareTo(o.traitName);
		}
	}
}
