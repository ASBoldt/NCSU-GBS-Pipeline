/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS.switchgrass;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Arrays;

/**
 *
 * @author fl262
 */
public class Distribution {
	double boundary[];
	int[] groupCount;
	int ngroup;
	public Distribution (int[] a, int ngroup) {
		Arrays.sort(a);
		defineBoundary(a, ngroup);
		this.ngroup = ngroup;
		this.calGroupCount(a, ngroup);
	}

	public int[] getGroupCount () {
		return groupCount;
	}

	public double[][] getBounds () {
		double[][] bounds = new double[ngroup][2];
		for (int i = 0; i < ngroup; i++) {
			bounds[i][0] = this.boundary[i];
			bounds[i][1] = this.boundary[i+1];
		}
		return bounds;
	}

	public int getNumberOfGroup () {
		return this.ngroup;
	}

	public void calGroupCount(int[] a, int ngroup) {
		groupCount = new int[ngroup];
		for (int i = 0; i < ngroup; i++) {
			groupCount[i] = 0;
		}
		for (int i = 0; i < a.length; i++) {
			double d = (double)a[i];
			int hit = Arrays.binarySearch(boundary, d);
			if (hit < -1 && hit >= -boundary.length) {
				groupCount[-hit-2]++;
			}
			else if (hit >= 0 && hit < boundary.length -1) {
				groupCount[hit]++;
			}
			else if (hit == boundary.length-1) {
				groupCount[hit-1]++;
			}
		}
	}

	public void defineBoundary(int[] a, int ngroup) {
		double min = (double)a[0];
		double max = (double)a[a.length-1];
		double increment = (max-min) / ngroup;
		boundary = new double[ngroup+1];
		boundary[0] = min;
		boundary[ngroup] = max;
		for (int i = 1; i < ngroup; i++) {
			boundary[i] = boundary[i-1] + increment;
		}
	}

	public void writeDistribution (String outfileS) {
		double[][] bounds = this.getBounds();
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outfileS), 65536);
			for (int i = 0; i < this.ngroup; i++) {
				bw.write(String.valueOf(bounds[i][0]) + "\t" + String.valueOf(bounds[i][1]) + "\t" + String.valueOf(this.groupCount[i]));
				bw.newLine();
			}
			bw.flush();
			bw.close();
		}
		catch (Exception e) {
			System.out.println("Error occurred while writing distribution file " + outfileS);
		}
	}
}
