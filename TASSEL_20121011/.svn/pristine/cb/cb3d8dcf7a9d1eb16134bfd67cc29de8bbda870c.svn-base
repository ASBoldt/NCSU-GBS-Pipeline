/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS.switchgrass;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.util.Arrays;
import cern.colt.GenericPermuting;
/**
 *
 * @author fl262
 */
public class LinkageGroup implements Comparable <LinkageGroup> {
	int[] marker;

	public LinkageGroup () {};

	public LinkageGroup (int[] marker) {
		this.marker = marker;
	}

	public Matrix getR2MatrixInGroup (Matrix r2Matrix) {
		float[][] va = new float[marker.length][marker.length];
		for (int i = 1; i < va.length; i++) {
			for (int j = 0; j < i; j++) {
				va[i][j] = r2Matrix.value[marker[i]][marker[j]];
				va[j][i] = va[i][j];
			}
		}
		for (int i = 0; i < va.length; i++) {
			va[i][i] = r2Matrix.value[i][i];
		}
		return new Matrix(va);
	}

	public void orderMarkers (Matrix r2Matrix, Boolean ifseed) {
		this.iniOrder(r2Matrix);
		int[] newMarker = Arrays.copyOf(marker, marker.length);
		int size;
		if (ifseed) {
			size = 10;
			if (size > newMarker.length) size = newMarker.length;
			this.getBestSeed(r2Matrix, newMarker, size);
		}
		else {
			size = 2;
		}
		for (int i = size; i < marker.length; i++) {
			//int closePlace = findClosePlace(newMarker, i1, marker[i1], r2Matrix);
			//int bestPlace = findBestPlace(newMarker, i1, marker[i1], closePlace, r2Matrix);
			int bestPlace = this.findBestPlaceBySumOfR2(newMarker, i, marker[i], r2Matrix);
			newMarker = moveNewMarker(newMarker, i, marker[i], bestPlace);
		}
		marker = newMarker;
	}

	private int[] getBestSeed(Matrix r2Matrix, int[] newMarker, int seedLength) {
		long nPermu = 1;
		for (int i = 1; i < seedLength + 1; i++) {
			nPermu *= i;
		}
		int[] order = new int[seedLength];
		int[] bestOrder = new int[seedLength];
		int[] seedMarker = new int[seedLength];
		float highestSum = 0;
		for (long i = 0; i < nPermu; i++) {
			order = GenericPermuting.permutation(i+1, seedLength);
			float sum = this.getSumOfR2(r2Matrix, newMarker, order);
			if (sum > highestSum) {
				highestSum = sum;
				bestOrder = Arrays.copyOf(order, seedLength);
			}
		}
		for (int i = 0; i < seedMarker.length; i++) {
			seedMarker[i] = newMarker[bestOrder[i]];
		}
		System.arraycopy(seedMarker, 0, newMarker, 0, seedLength);
		return newMarker;
	}

	private void  printIntArray (int[] order) {
		for (int i = 0; i < order.length; i++) {
			System.out.print(order[i]+"\t");
		}
		System.out.println();
	}

	private void iniOrder (Matrix r2Matrix) {
		float[] infoAmount = new float[marker.length];
		for (int i = 0; i < marker.length; i++) {
			infoAmount[i] = 0;
			for (int j = 0; j < marker.length; j++) {
				infoAmount[i] += r2Matrix.value[i][j];
			}
		}
		for (int i = 0; i < marker.length - 1; i++) {
			for (int j = i + 1; j < marker.length; j++) {
				if (infoAmount[i] < infoAmount[j]) {
					int midmarker = marker[i]; marker[i] = marker[j]; marker[j] = midmarker;
					float midinfo = infoAmount[i]; infoAmount[i] = infoAmount[j]; infoAmount[j] = midinfo;
				}
			}
		}
	}

	public int getSize() {
		return this.marker.length;
	}

	private float getSumOfR2 (Matrix r2Matrix, int[] newMarker, int[] order) {
		float sum = 0;
		for (int i = 0; i < order.length - 1; i++) {
			sum += r2Matrix.value[newMarker[order[i]]][newMarker[order[i+1]]];
		}
		return sum;
	}

	private int[] moveNewMarker (int[] newMarker, int size, int another, int bestPlace) {
		for (int i = size; i > bestPlace; i--) {
			newMarker[i] = newMarker[i-1];
		}
		newMarker[bestPlace] = another;
		return newMarker;
	}

	private int findBestPlaceBySumOfR2 (int[] newMarker, int size, int query, Matrix r2Matrix) {
		float[] sums = new float[size+1];
		for (int i = 0; i < sums.length; i++) {
			sums[i] = 0;
			if (i == 0) {
				sums[i] += r2Matrix.value[newMarker[i]][query];
				for (int j = 0; j < size - 1; j++) {
					sums[i] += r2Matrix.value[newMarker[j]][newMarker[j+1]];
				}
			}
			else if (i == size) {
				for (int j = 0; j < size - 1; j++) {
					sums[i] += r2Matrix.value[newMarker[j]][newMarker[j+1]];
				}
				sums[i] += r2Matrix.value[query][newMarker[i-1]];
			}
			else {
				for (int j = 0; j < i - 1; j++) {
					sums[i] += r2Matrix.value[newMarker[j]][newMarker[j+1]];
				}
				sums[i] += r2Matrix.value[query][newMarker[i-1]];
				sums[i] += r2Matrix.value[query][newMarker[i]];
				for (int j = i; j < size - 1; j++) {
					sums[i] += r2Matrix.value[newMarker[j]][newMarker[j+1]];
				}
			}
		}
		float HighestSum = 0;
		int bestPlace = 0;
		for (int i = 0; i < sums.length; i++)  {
			if (sums[i] > HighestSum) {
				HighestSum = sums[i];
				bestPlace = i;
			}
		}
		return bestPlace;
	}

	private int findBestPlace (int[] newMarker, int size, int query, int closePlace, Matrix r2Matrix) {
		int bestPlace = closePlace;
		float leftSum = 0;
		float rightSum = 0;
		if (closePlace == 0) {
			leftSum = r2Matrix.value[query][newMarker[closePlace]] + r2Matrix.value[newMarker[closePlace]][newMarker[closePlace+1]];
			rightSum = r2Matrix.value[query][newMarker[closePlace]] + r2Matrix.value[query][newMarker[closePlace+1]];
		}
		else if (closePlace == size - 1) {
			leftSum = r2Matrix.value[query][newMarker[closePlace-1]] + r2Matrix.value[query][newMarker[closePlace]];
			rightSum = r2Matrix.value[query][newMarker[closePlace]] + r2Matrix.value[newMarker[closePlace-1]][newMarker[closePlace]];
		}
		else {
			leftSum = r2Matrix.value[query][newMarker[closePlace-1]] + r2Matrix.value[query][newMarker[closePlace]] + r2Matrix.value[newMarker[closePlace]][newMarker[closePlace+1]];
			rightSum = r2Matrix.value[newMarker[closePlace-1]][newMarker[closePlace]] + r2Matrix.value[query][newMarker[closePlace]] + r2Matrix.value[query][newMarker[closePlace+1]];
		}
		if (leftSum > rightSum) bestPlace = closePlace;
		else bestPlace = closePlace + 1;
		return bestPlace;
	}

	private int findClosePlace(int[] newMarker, int size, int query, Matrix r2Matrix) {
		int closePlace = 0;
		for (int i = 1; i < size; i++) {
			if (r2Matrix.value[newMarker[closePlace]][query] < r2Matrix.value[newMarker[i]][query]) {
				closePlace = i;
			}
		}
		return closePlace;
	}

	public void writeLinkageGroup(DataOutputStream dos) {
		try {
			dos.writeInt(this.marker.length);
			for (int i = 0; i < marker.length; i++) {
				dos.writeInt(marker[i]);
			}
		}
		catch (Exception e) {
			System.out.println("Error occurred while writing LinkageGroup");
		}
	}

	public void readLinkageGroup(DataInputStream dis) {
		try {
			marker = new int[dis.readInt()];
			for (int i = 0; i < marker.length; i++) {
				marker[i] = dis.readInt();
			}
		}
		catch (Exception e) {
			System.out.println("Error occurred while reading LinkageGroup");
		}
	}

	public int compareTo(LinkageGroup o) {
		return this.getSize()-o.getSize();
	}
}
