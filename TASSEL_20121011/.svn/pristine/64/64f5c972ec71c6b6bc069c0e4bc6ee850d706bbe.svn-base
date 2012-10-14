/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS.switchgrass;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeMap;
import org.apache.commons.math.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math.random.RandomDataImpl;

/**
 *
 * @author fl262
 */
public class RelationUtil {
	CorrRecord[] rsquareMaList; //half matrix or significant rsquare of half matrix
	CorrRecord[] doubleRsquareMaList; //switch queryIndex and hitIndex, so rsquareMaList.length * 2
	Matrix rsquareMatrix; //Pearson correlations
	Matrix rMatrix;
	Matrix ldMatrix; // Inferences about linkage disequilibrium. Biometrics, 35, 235–254. Calculation without phase
	String[][] hapmapRecord;
	byte[][] genoValue;

	public RelationUtil (String infileS, boolean ifbinary) {
		if (ifbinary) {
			try {
				DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(infileS), 65536));
				rsquareMaList = new CorrRecord[dis.readInt()];
				for (int i = 0; i < rsquareMaList.length; i++) {
					rsquareMaList[i] = new CorrRecord(dis.readInt(), dis.readInt(), dis.readFloat());
				}
			} catch (Exception e) {
				System.out.println("Error occurred while reading " + infileS + " " + e.toString());
			}
		}
		Arrays.sort(this.rsquareMaList);
		System.out.println("Correlation matrix " + infileS + " is loaded");
	}

	public RelationUtil (String hapmap) {
		this.readHapmap(hapmap);
	}

	public Matrix convert2DMatrixFromRsquareMatrix (float threshhold) {
		Matrix distanceMatrix = this.rsquareMatrix.getDMatrixFromR2Matrix(threshhold);
		System.out.println("Distance matrix is generated");
		return distanceMatrix;
	}

	public Matrix convert2DMatrixFromRsquareMatrix () {
		Matrix distanceMatrix = this.rsquareMatrix.getDMatrixFromR2Matrix();
		System.out.println("Distance matrix is generated");
		return distanceMatrix;
	}

	public Matrix conver2DMatrixFromLdsquareMatrix (float threshhold) {
		Matrix distanceMatrix = this.ldMatrix.getDMatrixFromR2Matrix(threshhold);
		System.out.println("Distance matrix is generated");
		return distanceMatrix;
	}

	public void doubleMaList () {
		doubleRsquareMaList = new CorrRecord[this.rsquareMaList.length*2];
		for (int i = 0; i < this.rsquareMaList.length; i++) {
			this.doubleRsquareMaList[i] = this.rsquareMaList[i];
			this.doubleRsquareMaList[this.rsquareMaList.length+i] = new CorrRecord(this.rsquareMaList[i].hitIndex, this.rsquareMaList[i].queryIndex, this.rsquareMaList[i].rsquare);
		}
		Arrays.sort(this.doubleRsquareMaList, new sortByIndex());
	}

	public void readHapmap (String hapmap) {
		try {
			ArrayList<String[]> hapmapRecordList = new ArrayList();
			BufferedReader br = new BufferedReader (new FileReader(hapmap), 65536);
			String temp = br.readLine();
			while ((temp = br.readLine()) != null) {
				String[] tem = temp.split(("\t"));
				hapmapRecordList.add(tem);
			}
			hapmapRecord = hapmapRecordList.toArray(new String[hapmapRecordList.size()][]);
			System.out.println("HapMap is loaded. " + String.valueOf(hapmapRecord.length) + " markers");
		} catch (Exception e) {
			System.out.println("Error while reading " + hapmap + " " + e.toString());
			System.exit(0);
		}
	}

	public void assignGenoValue () {
		genoValue = new byte[hapmapRecord.length][hapmapRecord[0].length-11];
		for (int i = 0; i < hapmapRecord.length; i++) {
			double freA = this.getAllelFrequency(hapmapRecord[i]);
			for (int j = 11; j < hapmapRecord[0].length; j++) {
				int k = j - 11;
				double r;
				if (hapmapRecord[i][j].startsWith("A")) {
					r = Math.random();
					if (r < freA) genoValue[i][k] = 0;
					else genoValue[i][k] = 1;

				}
				else if (hapmapRecord[i][j].startsWith("C")) {
					r = Math.random();
					if (r < freA) genoValue[i][k] = 1;
					else genoValue[i][k] = 2;
				}
				else if (hapmapRecord[i][j].startsWith("R")) {
					genoValue[i][k] = 1;
				}
				else {
					r = Math.random();
					double bound1 = freA*freA;
					double bound2 = (1-freA)*(1-freA) + bound1;
					if (r < bound1) genoValue[i][k] = 0;
					else if (r < bound2) genoValue[i][k] = 2;
					else genoValue[i][k] = 1;
				}
			}
		}
	}

	public void assignGenoValueSimple () {
		genoValue = new byte[hapmapRecord.length][hapmapRecord[0].length-11];
		for (int i = 0; i < hapmapRecord.length; i++) {
			for (int j = 11; j < hapmapRecord[0].length; j++) {
				int k = j - 11;
				if (hapmapRecord[i][j].startsWith("A")) {
					genoValue[i][k] = 0;
				}
				else if (hapmapRecord[i][j].startsWith("C")) {
					genoValue[i][k] = 2;
				}
				else if (hapmapRecord[i][j].startsWith("R")) {
					genoValue[i][k] = 1;
				}
				else {
					genoValue[i][k] = Byte.MIN_VALUE;
					//genoValue[i][k] = 1;
				}
			}
		}
		System.out.println("GenotypeSimple(012) value is assigned. " + String.valueOf(genoValue.length) + " markers and " + String.valueOf(genoValue[0].length) + " taxa");
	}

	public void assignGenoValueMMC () {
		genoValue = new byte[hapmapRecord.length][hapmapRecord[0].length-11];
		for (int i = 0; i < hapmapRecord.length; i++) {
			for (int j = 11; j < hapmapRecord[0].length; j++) {
				int k = j - 11;
				if (hapmapRecord[i][j].startsWith("A")) {
					genoValue[i][k] = 0;
				}
				else if (hapmapRecord[i][j].startsWith("C")) {
					genoValue[i][k] = 2;
				}
				else if (hapmapRecord[i][j].startsWith("R")) {
					genoValue[i][k] = 1;
				}
				else {
					genoValue[i][k] = 1;
				}
			}
		}
		System.out.println("GenotypeMMC(012) value is assigned. " + String.valueOf(genoValue.length) + " markers and " + String.valueOf(genoValue[0].length) + " taxa");
	}

	public void assignGenoValuePA () {
		genoValue = new byte[hapmapRecord.length][hapmapRecord[0].length-11];
		for (int i = 0; i < hapmapRecord.length; i++) {
			for (int j = 11; j < hapmapRecord[0].length; j++) {
				int k = j - 11;
				if (hapmapRecord[i][k].startsWith("N")) {
					genoValue[i][k] = 0;
				}
				else {
					genoValue[i][k] = 1;
				}
			}
		}
		System.out.println("GenotypePA(01) value is assigned. " + String.valueOf(genoValue.length) + " markers and " + String.valueOf(genoValue[0].length) + " taxa");
	}

	public void getLdSquareMatrix (int sampleSize) {
		float[][] ldValue = new float[sampleSize][sampleSize];
		for (int i = 1; i < sampleSize; i++) {
			for (int j = 0; j < i; j++) {
				ldValue[i][j] = this.getLdSquare(genoValue[i], genoValue[j]);
				ldValue[j][i] = ldValue[i][j];
			}
		}
		for (int i = 0; i < sampleSize; i++) ldValue[i][i] = 1;
		this.ldMatrix = new Matrix(ldValue);
		System.out.println("Ld matrix is generated");
	}

	public void getRMatrix (int sampleSize) {
		float[][] rValue = new float[sampleSize][sampleSize];
		for (int i = 1; i < sampleSize; i++) {
			for (int j = 0; j < i; j++) {
				rValue[i][j] = this.getPearsonR(genoValue[i], genoValue[j]);
				rValue[j][i] = rValue[i][j];
			}
		}
		for (int i = 0; i < sampleSize; i++) {rValue[i][i] = 1;}
		this.rMatrix = new Matrix(rValue);
		System.out.println("R matrix is generated");
	}

	public void getRsquareMatrix (int sampleSize) {
		float[][] r2Value = new float[sampleSize][sampleSize];
		for (int i = 1; i < sampleSize; i++) {
			for (int j = 0; j < i; j++) {
				r2Value[i][j] = this.getRsquare(genoValue[i], genoValue[j]);
				r2Value[j][i] = r2Value[i][j];
			}
		}
		for (int i = 0; i < sampleSize; i++) {r2Value[i][i] = 1;}
		this.rsquareMatrix = new Matrix(r2Value);
		System.out.println("Rsquare matrix is generated");
	}


	public void getCorrelationMaList () {//lower-left matrix
		rsquareMaList = new CorrRecord[(genoValue.length-1)*genoValue.length/2];
		int count = 0;
		for (int i = 1; i < genoValue.length; i++) {
			for (int j = 0; j < i; j++) {
				rsquareMaList[count] = new CorrRecord(i, j, getRsquare(genoValue[i], genoValue[j]));
				count++;
			}
		}
		Arrays.sort(this.rsquareMaList);
		System.out.println("Half correlation matrix is calculated");
	}

	public void getSigCorrelationMaList (float sigRsquare) {//significant lower-left matrix, sigRsquare comes from pumertation
		ArrayList<CorrRecord> sigList = new ArrayList();
		for (int i = 1; i < genoValue.length; i++) {
			for (int j = 0; j < i; j++) {
				float rsquare = getRsquare(genoValue[i], genoValue[j]);
				if (rsquare > sigRsquare) sigList.add(new CorrRecord(i, j, rsquare));
			}
		}
		this.rsquareMaList = sigList.toArray(new CorrRecord[sigList.size()]);
		Arrays.sort(this.rsquareMaList);
		System.out.println("Significant correlation matrix is calculated");
	}

	public void getRsquareDistribution (String outfileS) {
		int ngroup = 100;
		float min = this.rsquareMaList[0].rsquare;
		float max = this.rsquareMaList[this.rsquareMaList.length-1].rsquare;
		float range = max - min;
		float step = range / ngroup;
		float[] border = new float[ngroup+1];
		int[] dis = new int[ngroup];
		for (int i = 0; i < border.length; i++) {
			border[i] = i * step + min;
		}
		for (int i = 0; i < this.rsquareMaList.length; i++) {
			int hit = Arrays.binarySearch(border, this.rsquareMaList[i].rsquare);
			if (hit < -1) dis[-hit-2]++;
			else if (hit >=0 && hit < ngroup) dis[hit]++;
		}
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outfileS), 65536);
			for (int i = 0; i < ngroup; i++) {
				bw.write(String.valueOf(border[i]) + "\t" + String.valueOf(dis[i]));
				bw.newLine();
			}
			bw.flush();
			bw.close();
		}
		catch (Exception e) {
			System.out.println("Error occurred while writing rsquare distribution " + outfileS + " " + e.toString());
		}
		System.out.println("Rsquare distribution is generated in " + outfileS);
	}

	public void edgeNumberFilter (int min, int max) {
		int cnt = 1;
		TreeMap<Integer, Integer> queryIndex2Edges = new TreeMap();
		ArrayList<Integer> queryList = new ArrayList();
		for (int i = 1; i < this.doubleRsquareMaList.length; i++) {
			if (this.doubleRsquareMaList[i].queryIndex == this.doubleRsquareMaList[i-1].queryIndex) {
				cnt++;
			}
			else {
				queryList.add(this.doubleRsquareMaList[i-1].queryIndex);
				queryIndex2Edges.put(this.doubleRsquareMaList[i-1].queryIndex, cnt);
				cnt = 1;
			}
		}
		queryList.add(this.doubleRsquareMaList[this.doubleRsquareMaList.length-1].queryIndex);
		queryIndex2Edges.put(this.doubleRsquareMaList[this.doubleRsquareMaList.length-1].queryIndex, cnt);
		Integer[] queryArray = queryList.toArray(new Integer[queryList.size()]);
		ArrayList<Integer> deletedList = new ArrayList();
		for (int i = 0; i < queryArray.length; i++) {
			if (queryIndex2Edges.get(queryArray[i]) < min || queryIndex2Edges.get(queryArray[i]) > max) deletedList.add(queryArray[i]);
		}
		Integer[] deletedArray = deletedList.toArray(new Integer[deletedList.size()]);
		Arrays.sort(deletedArray);
		ArrayList<CorrRecord> newList = new ArrayList();
		for (int i = 0; i < this.doubleRsquareMaList.length; i++) {
			int hit = Arrays.binarySearch(deletedArray, this.doubleRsquareMaList[i].queryIndex);
			if (hit < 0) newList.add(this.doubleRsquareMaList[i]);
		}
		CorrRecord[] newArray = newList.toArray(new CorrRecord[newList.size()]);
		Arrays.sort(newArray, new sortByIndex());
		this.doubleRsquareMaList = newArray;
		newArray = null;
	}

	public void getEdgeDistribution (String outfileS) {
		ArrayList<Integer> countList = new ArrayList();
		int cnt = 1;
		for (int i = 1; i < this.doubleRsquareMaList.length; i++) {
			if (this.doubleRsquareMaList[i].queryIndex == this.doubleRsquareMaList[i-1].queryIndex) {
				cnt++;
			}
			else {
				countList.add(cnt);
				cnt = 1;
			}
		}
		countList.add(cnt);
		Integer[] countArray = countList.toArray(new Integer[countList.size()]);
		Arrays.sort(countArray);
		int min = countArray[0];
		int max = countArray[countArray.length-1];
		int ngroup = max - min + 1;
		int[] dis = new int[ngroup];
		int[] count = new int[ngroup];
		for (int i = 0; i < dis.length; i++) {
			count[i] = min + i;
			dis[i] = 0;
		}
		for (int i = 0; i < countArray.length; i++) {
			int hit = Arrays.binarySearch(count, countArray[i]);
			dis[hit]++;
		}
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outfileS), 65536);
			for (int i = 0; i < dis.length; i++) {
				bw.write(String.valueOf(count[i]) + "\t" + String.valueOf(dis[i]));
				bw.newLine();
			}
			bw.flush();
			bw.close();
		}
		catch (Exception e) {
			System.out.println("Error occurred while writing edge distribution " + outfileS + " " + e.toString());
		}
	}

	public double getAllelFrequency (String[] record) {
		int a = 0, c = 0;
		for (int i = 11; i < record.length; i++) {
			if (record[i].startsWith("A")) a++;
			else if (record[i].startsWith("C")) c++;
			else if (record[i].startsWith("R")) {
				a++; c++;
			}
		}
		double freA = (double)a/(double)(a+c);
		return  freA;
	}
	public float getSignificantRsquare (float significance) {//used in half matrix
		int index = (int)Math.floor((1-significance) * this.rsquareMaList.length);
		return this.rsquareMaList[index].rsquare;
	}

	public float permutation (int sampleSize, float significance) {
		ArrayList<Float> pList = new ArrayList();
		int n = genoValue[0].length;
		byte[] a = new byte[n];
		for (int i = 0; i < sampleSize; i++) {
			int index1 = (int)Math.floor(Math.random()*genoValue.length);
			int index2 = (int)Math.floor(Math.random()*genoValue.length);
			while (index1 == index2) index2 = (int)Math.floor(Math.random()*genoValue.length);
			RandomDataImpl rdi = new RandomDataImpl();
			for (int j = 0; j < 1000; j++) {
				int[] order = rdi.nextPermutation(n, n);
				for (int k = 0; k < n; k++) {
					a[k] = genoValue[index1][order[k]];
					pList.add(getRsquare(a, genoValue[index2])); //for RsquareMatrix
					//pList.add(this.getLdSquare(a, genoValue[index2])); //for LdsquareMatrix
				}
			}
		}
		Float[] pArray = pList.toArray(new Float[pList.size()]);
		Arrays.sort(pArray);
		int border = (int)Math.floor((1-significance) * pArray.length);
		return pArray[border]; 
	}

	public float getLdSquare (byte[] a, byte[] b) {
		float ld = this.getLd(a, b);
		return (float)Math.pow(ld, 2);
	}

	public float getLd (byte[] a, byte[] b) {
		float[][] genoPair = new float[3][3];
		int n = a.length;
		for (int i = 0; i < a.length; i++) {
			genoPair[a[i]][b[i]]++;
		}
		for (int i = 0; i < genoPair.length; i++) {
			for (int j = 0; j < genoPair[0].length; j++) {
				genoPair[i][j] = genoPair[i][j]/n;
			}
		}
		float delta = (2 * genoPair[0][0] + genoPair[0][1] + genoPair[1][0] + (float)0.5 *  genoPair[1][1]) / n;
		float Pa = (2 * (genoPair[0][0] + genoPair[0][1] + genoPair[0][2]) + (genoPair[1][0] + genoPair[1][1] + genoPair[1][2])) / (2 * n);
		float Pb = (2 * (genoPair[0][0] + genoPair[1][0] + genoPair[2][0]) + (genoPair[0][1] + genoPair[1][1] + genoPair[2][1])) / (2 * n);
		delta = delta - Pa * Pb;
		return delta;
	}

	public float getEuclideanD (byte[] a, byte[] b) {
		float d;
		float sum2 = 0;
		for (int i = 0; i < a.length; i++) {
			sum2 += Math.pow((a[i]-b[i]), 2);
		}
		d = (float)Math.pow(sum2, 0.5);
		return d;
	}

	public float getRsquare (byte[] a, byte[] b) {
		return (float)Math.pow(this.getPearsonR(a, b), 2);
	}

	public int getMarkerCount () {
		return this.genoValue.length;
	}

	public int getTaxaCount () {
		return this.genoValue[0].length;
	}

	public float getPearsonR (byte[] a, byte[] b) {
		ArrayList<Integer> indexList = new ArrayList();
		for (int i = 0; i < a.length; i++) {
			if (a[i] + b[i] > -1) indexList.add(i);
		}
		Integer[] index = indexList.toArray(new Integer[indexList.size()]);
		double[] x = new double[index.length];
		double[] y = new double[index.length];
		for (int i = 0; i < index.length; i++) {
			x[i] = (double)a[index[i]];
			y[i] = (double)b[index[i]];
		}
		return (float)(new PearsonsCorrelation().correlation(x, y));
	}

	public void writeCorrelationMaList (String outfileS) {
		try {
			DataOutputStream dos =new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfileS),65536));
			dos.writeInt(this.rsquareMaList.length);
			for (int i = 0; i < this.rsquareMaList.length; i++) {
				rsquareMaList[i].writeBinary(dos);
			}
			dos.flush();
			dos.close();
			System.out.println("Correlation matrix is written");
		}
		catch (Exception e) {
			System.out.println("Error occurred while writing " + outfileS + " " + e.toString());
		}
	}

	public void writeMDSMatrix (String outfileS) {//just for MDS plot
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outfileS), 65536);
			bw.write(String.valueOf(rsquareMatrix.getMatrixSize())+"\t\t");
			for (int i = 0; i < rsquareMatrix.getMatrixSize(); i++) {
				bw.write(String.valueOf(i)+"\t");
			}
			bw.newLine();
			for (int i = 0; i < rsquareMatrix.getMatrixSize(); i++) {
				bw.write(String.valueOf(i)+ "\tred\t");
				for (int j = 0; j < rsquareMatrix.getMatrixSize(); j++) {
					bw.write(String.valueOf(rsquareMatrix.value[i][j])+ "\t");
				}
				bw.newLine();
			}
			bw.flush();
			bw.close();
		}
		catch (Exception e ) {
			System.out.println("Error occurred while writing FullMatrix file " + outfileS + " " + e.toString());
		}
	}

	public void writeCytoscapeFormat (String outfileS) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outfileS), 65536);
			for (int i = 0; i < this.rsquareMaList.length; i++) {
				bw.write(String.valueOf(this.rsquareMaList[i].queryIndex + "\t" + String.valueOf(this.rsquareMaList[i].hitIndex)));
				bw.newLine();
			}
			bw.flush();
			bw.close();
		}
		catch (Exception e) {
			System.out.println("Error occurred while writing " + outfileS + " " + e.toString());
		}
	}

	public void writeGenoValue4MMC (String outfileS) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outfileS), 65536);
			bw.write("Example,");
			StringBuilder sb = new StringBuilder();
			for (int i = 0; i < this.getTaxaCount(); i++) {
				sb.append("T").append(String.valueOf(i)).append(",");
			}
			sb.deleteCharAt(sb.length()-1);
			bw.write(sb.toString());
			bw.newLine();
			for (int i = 0; i < this.getMarkerCount(); i++) {
				sb = new StringBuilder();
				sb.append("m").append(String.valueOf(i)).append(",");
				for (int j = 0; j < this.getTaxaCount(); j++) {
					sb.append(String.valueOf(genoValue[i][j])).append(",");
				}
				sb.deleteCharAt(sb.length()-1);
				bw.write(sb.toString());
				bw.newLine();
			}
			bw.flush();
			bw.close();
		}
		catch (Exception e) {
			System.out.println("Error occurred while writing MMC file" + outfileS + " " + e.toString());
		}
	}
}

class PearsonCorrelation {
	float r;

	public PearsonCorrelation (byte[] x, byte[]y) {
		this.r = this.getCorrelation(x, y);
	}

	public float getMean (byte[] a) {
		int sum = 0;
		for (int i = 0; i < a.length; i++) {
			sum += a[i];
		}
		return (float)sum / (float)a.length;
	}

	public float getVariance (byte[] a, float mean) {
		float variance = 0;
		for (int i = 0; i < a.length; i++) {
			variance += Math.pow((float)a[i]-mean, 2)/a.length;
		}
		return variance;
	}

	public float getCovariance(byte[] x, byte[] y, float xmean, float ymean) {
		float covariance = 0;
		for (int i = 0; i < x.length; i++) {
			covariance += (x[i] - xmean)* (y[i] - ymean) / x.length;
		}
		return covariance;
	}

	public float getStandardDiviation(byte[] a, float mean) {
		return (float)Math.pow(getVariance(a, mean), 0.5);
	}

	public float getCorrelation (byte[] x, byte[] y) {
		float xmean = this.getMean(x);
		float ymean = this.getMean(y);
		return this.getCovariance(x, y, xmean, ymean) / this.getStandardDiviation(x, xmean) / this.getStandardDiviation(y, ymean);
	}

	public float getRsquare() {
		return this.r * this.r;
	}
}