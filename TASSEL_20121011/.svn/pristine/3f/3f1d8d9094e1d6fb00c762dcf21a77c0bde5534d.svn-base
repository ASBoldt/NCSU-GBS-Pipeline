/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS.repeat;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;

/**
 *
 * @author fl262
 */
public class annotateByKeyword {
	static String sourceFileS = filterBlastRec.desFilesS;
	static String annotationDirS = setQuery.forSASDir + "annotation/";
	static File SASResultF = processGlmSAS.desFile;
	static String[] keywords = {"chloroplast", "mitochondri"};

	blastRecArray bra;
	public annotateByKeyword () {
		getData ();
		outputAnnotation ();
	}
	public void visualizeAnnotation(String sourceFileS, String keyword) {
		ArrayList<String> annoAl = new ArrayList();
		String[] annoArray;
		try {
			BufferedReader br = new BufferedReader(new FileReader(sourceFileS), 65536);
			String temp = br.readLine();
			while (temp != null) {
				annoAl.add(temp);
				temp = br.readLine();
			}
			br.close();
		}
		catch (Exception e) {
			System.out.println (e.toString());
			System.out.println ("Error occured in " + sourceFileS);
		}
		annoArray = annoAl.toArray(new String[annoAl.size()]);
		String[] nameArray = new String[annoArray.length];
		for (int i = 0; i < annoArray.length; i++) {
			String[] temp = annoArray[i].split("\t");
			String[] tem = temp[0].split("\\|");
			nameArray[i] = tem[0];
		}
		Arrays.sort(nameArray);
		GlmRecSASArray grsaIn = new GlmRecSASArray(SASResultF, true);
		ArrayList<glmRecSAS> grsAl= new ArrayList();
		for (int i = 0; i < grsaIn.sigMarkerNum; i++) {
			int hit = Arrays.binarySearch(nameArray, grsaIn.sigMarkers[i].traitName);
			if (hit >= 0) {
				grsAl.add(grsaIn.sigMarkers[i]);
			}
		}
		glmRecSAS[] grsArray = grsAl.toArray(new glmRecSAS[grsAl.size()]);
		GlmRecSASArray grsaOut = new GlmRecSASArray(grsArray, grsaIn.title);
		String desFileS = annotationDirS + keyword + "_Glm.txt";
		grsaOut.writeFile(new File(desFileS));
		String figureFileS = annotationDirS + keyword + ".png";
		new visualizeGlm(new File(desFileS), new File(figureFileS));
		

	}
	public void writeAnnotation (String keyword) {
		String desFileS = annotationDirS + keyword + ".txt";
		ArrayList<Integer> indexAl = new ArrayList();
		indexAl.add(0);
		for (int i = 1; i < bra.recNum; i++) {
			if (!bra.recs[i].item[0].equals(bra.recs[indexAl.get(indexAl.size()-1)].item[0])) {
				indexAl.add(i);
			}
		}
		Integer[] indexArray = indexAl.toArray(new Integer[indexAl.size()]);
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(desFileS), 65536);
			for (int i = 0; i < indexArray.length; i++) {
				if (bra.recs[indexArray[i]].item[2].contains(keyword)) {
					bw.write(bra.recs[indexArray[i]].item[0] + "\t" + keyword + "\t" + bra.recs[indexArray[i]].item[2] + "\t" + bra.recs[indexArray[i]].item[5]);
					bw.newLine();
				}
			}
			bw.flush();
			bw.close();
		}
		catch (Exception e) {
			System.out.println (e.toString());
			System.out.println ("Error occured in " + desFileS);
		}
		visualizeAnnotation(desFileS, keyword);
	}
	public void outputAnnotation () {
		new File (annotationDirS).mkdir();
		for (int i = 0; i < keywords.length; i++) {
			writeAnnotation(keywords[i]);
		}
	}
	public void getData() {
		bra = new blastRecArray (new File(sourceFileS));
	}
	public static void main (String[] agrs) {
		new annotateByKeyword ();
	}
}
