/*
process Tassel result table, set a filter (by p value), add marker information
 */

package net.maizegenetics.genome.GBS.repeat;

import java.io.*;
import java.util.Arrays;

/**
 *
 * @author fl262
 */
public class processGlmTassel {
	String desDir = setQuery.forTasselDir;
	String desFile = desDir+"filteredGLM.txt";
	File sourceGenoFile = setQuery.sourceGenoFile644;
	String dirFile = "all_GLM.txt";
	float eFilter = (float) 0.0001;
	GlmRecArrayTassel gra;
	GeneticInfoArray gif;

	public processGlmTassel () {
		File sourceGlmFile = new File (desDir,dirFile);
		gra = new GlmRecArrayTassel(sourceGlmFile, eFilter);
		gif = new GeneticInfoArray(sourceGenoFile);
		setChrAndPosi();
	}
	public void setChrAndPosi () {
		Arrays.sort (gif.markers,new sortByMarkerName());
		String[] temps = new String[gif.markerNum];
		for (int i = 0; i < gif.markerNum; i++) {
			temps[i] = gif.markers[i].markerName;
		}
		for (int i = 0; i < gra.sigMarkerNum; i++) {
			int hit = Arrays.binarySearch(temps, gra.sigMarkers[i].markerName);
			gra.sigMarkers[i].setChAndPosi(gif.markers[hit].chrom, gif.markers[hit].position);
		}
		Arrays.sort (gra.sigMarkers);
		writeFile ();
		System.out.println ("Totally there are "+gra.traitNum+" traits under E-value of "+eFilter);
	}
	public void writeFile () {
		gra.writeFile(desFile);
	}

	public static void main (String args[]) {
		processGlmTassel pgt = new processGlmTassel ();
	}
}
