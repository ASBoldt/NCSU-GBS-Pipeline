/*
 read SAS output table, write filtered (by p value) result with marker information
 */

package net.maizegenetics.genome.GBS.repeat;

import java.io.*;
import java.util.ArrayList;
import java.util.Comparator;
/**
 *
 * @author fl262
 */
public class GlmRecSASArray implements Cloneable {
	int sigMarkerNum;
	glmRecSAS[] sigMarkers;
	String title;

	public GlmRecSASArray (glmRecSAS[] grsArray, String title) {
		sigMarkerNum = grsArray.length;
		sigMarkers = grsArray;
		this.title = title;
	}
	public GlmRecSASArray (File SASResultF, float eFilter) {
		readFile(SASResultF, eFilter);
	}

	public GlmRecSASArray (File SASResultF, boolean ifFinal) {
		if (ifFinal) {
			readFile(SASResultF);
		}
		else {
			readFile(SASResultF, 100);
		}
	}

	@Override
	public Object clone () {
		try {
			return super.clone();
		}
		catch (CloneNotSupportedException e) {
			return null;
		}
	}

	public void writeFile (File desFile) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(desFile), 65536);
			bw.write(title);
			bw.newLine();
			for (int i = 0; i < sigMarkerNum; i++) {
				StringBuilder sb = new StringBuilder();
				sb.append(sigMarkers[i].traitName).append("\t");
				sb.append(sigMarkers[i].markerName).append("\t");
				sb.append(sigMarkers[i].chrom).append("\t");
				sb.append(sigMarkers[i].position).append("\t");
				sb.append(sigMarkers[i].DF).append("\t");
				sb.append(sigMarkers[i].Estimate).append("\t");
				sb.append(sigMarkers[i].StandardizedEst).append("\t");
				sb.append(sigMarkers[i].StdErr).append("\t");
				sb.append(sigMarkers[i].tValue).append("\t");
				sb.append(sigMarkers[i].Probt);
				bw.write(sb.toString());
				bw.newLine();
			}
			bw.flush();
			bw.close();
		}
		catch (Exception e) {
			System.out.println(desFile.toString());
		}
	}
	public void readFile (File SASResultF) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(SASResultF), 65536);
			title = br.readLine();
			ArrayList<String> line = new ArrayList<String>();
			String temp = br.readLine();
			while (temp != null) {
				line.add(temp);
				temp = br.readLine();
			}
			br.close();
			String[] lineArray = line.toArray(new String[line.size()]);
			sigMarkerNum = lineArray.length;
			sigMarkers = new glmRecSAS[sigMarkerNum];
			for (int i = 0; i < sigMarkerNum; i++) {
				sigMarkers[i] = new glmRecSAS(lineArray[i], true);
			}
		}
		catch (Exception e) {
			System.out.println("reading " + SASResultF.toString());
		}
	}
	public void readFile (File SASResultF, float eFilter) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(SASResultF), 65536);
			String temp = br.readLine();
			title = temp.replaceFirst("Effect.Parameter", "Marker\tChrom\tPosition");
			ArrayList<String> filterLine = new ArrayList<String>();
			temp = br.readLine();
			while (temp != null) {
				String[] temps = temp.split("\t");
				if (!temps[1].startsWith("Inter")) {
					if (Float.valueOf(temps[8]) < eFilter) {
						filterLine.add(temp);
					}
				}
				temp = br.readLine();
			}
			br.close();
			String[] filterLineArray = filterLine.toArray(new String[filterLine.size()]);
			sigMarkerNum = filterLineArray.length;
			sigMarkers = new glmRecSAS[sigMarkerNum];
			for (int i = 0; i < sigMarkerNum; i++ ) {
				sigMarkers[i] = new glmRecSAS(filterLineArray[i]);
			}
		}
		catch (Exception e) {

		}
	}
}
class glmRecSAS implements Comparable<glmRecSAS> {
	String traitName;
	String markerName;
	byte chrom;
	int position;
	float DF, Estimate, StandardizedEst, StdErr, tValue, Probt;

	public glmRecSAS (String line, boolean ifFinal) {
		String[] temps = line.split("\t");
		traitName = temps[0];
		markerName = temps[1];
		chrom = Byte.valueOf(temps[2]);
		position = Integer.valueOf(temps[3]);
		DF = Float.valueOf(temps[4]);
		Estimate = Float.valueOf(temps[5]);
		StandardizedEst = Float.valueOf(temps[6]);
		StdErr = Float.valueOf(temps[7]);
		tValue = Float.valueOf(temps[8]);
		Probt = Float.valueOf(temps[9]);
	}

	public glmRecSAS (String filterLine) {
		String[] temps = filterLine.split("\t");
		traitName = temps[0];
		markerName = temps[1];
		DF = Float.valueOf(temps[3]);
		Estimate = Float.valueOf(temps[4]);
		StandardizedEst = Float.valueOf(temps[5]);
		StdErr = Float.valueOf(temps[6]);
		tValue = Float.valueOf(temps[7]);
		Probt = Float.valueOf(temps[8]);
	}

	public void setMarker (String markerName, byte chrom, int position) {
		this.markerName = markerName;
		this.chrom = chrom;
		this.position  = position;
	}
	public void setTraitName (int markerNum) {
		traitName = String.valueOf(Integer.valueOf(traitName) - markerNum);
	}
	public int compareTo(glmRecSAS o) {
		if (traitName.equals(o.traitName)){
			if (chrom < o.chrom) return -1;
			else if (chrom > o.chrom) return 1;
			else return position - o.position;
		}
		else {
			return traitName.compareTo(o.traitName);
		}
	}
}
class sortByMarker implements Comparator <glmRecSAS> {
	public int compare (glmRecSAS gi1, glmRecSAS gi2) {
		if (gi1.chrom < gi2.chrom) {
			return -1;
		}
		else if (gi1.chrom > gi2.chrom) {
			return 1;
		}
		else {return gi1.markerName.compareTo(gi2.markerName);}
	}
}
