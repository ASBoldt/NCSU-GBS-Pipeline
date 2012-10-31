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
import java.io.FileWriter;
import java.io.FileReader;
import java.util.ArrayList;

/**
 *
 * @author fl262
 */
public class Matrix {
	float[][] value;
	public Matrix (float[][] ma) {
		this.value = ma;
	}
	public Matrix (String infileS, boolean ifBinary) {
		this.readMatrix(infileS, ifBinary);
	}

	public int getMatrixSize() {
		return value.length;
	}

	public Matrix getDMatrixFromR2Matrix (float threshold) {
		float[][] newValue = new float[value.length][value.length];
		for (int i = 1; i < value.length; i++) {
			for (int j = 0; j < i; j++) {
				if (value[i][j] < threshold) {
					newValue[i][j] = 1; // or greater, like 10
				}
				else {
					newValue[i][j] = 0; // or 1 - value[i][j]
				}
				newValue[j][i] = newValue[i][j];
			}
		}
		return new Matrix(newValue);
	}

	public Matrix getDMatrixFromR2Matrix () {
		float[][] newValue = new float[value.length][value.length];
		for (int i = 1; i < value.length; i++) {
			for (int j = 0; j < i; j++) {
				newValue[i][j] = 1 -value[i][j];
				newValue[j][i] = newValue[i][j];
			}
		}
		return new Matrix(newValue);
	}

	public void readMatrix (String infileS, boolean ifBinary) {
		if (ifBinary) {
			try {
				DataInputStream dis =new DataInputStream(new BufferedInputStream(new FileInputStream(infileS),65536));
				int size = dis.readInt();
				value = new float[size][size];
				for (int i = 0; i < size; i++) {
					for (int j = 0; j < size; j++) {
						value[i][j] = dis.readFloat();
					}
				}
			}
			catch (Exception e) {
				System.out.println("Error occurred while reading Matrix " + infileS + " " + e.toString());
			}
		}
		else {
			try {
				BufferedReader br = new BufferedReader(new FileReader(infileS), 65536);
				ArrayList<String> recordList = new ArrayList();
				String temp;
				while ((temp = br.readLine()) != null) {
					recordList.add(temp);
				}
				String[] recordArray = recordList.toArray(new String[recordList.size()]);
				value = new float[recordArray.length][recordArray.length];
				for (int i = 0; i < recordArray.length; i++) {
					String[] tem = recordArray[i].split("\t");
					for (int j = 0; j < recordArray.length; j++) {
						value[i][j] = Float.valueOf(tem[j]);
					}
				}
			}
			catch (Exception e) {
				System.out.println("Error occurred while reading Matrix " + infileS + " " + e.toString());
			}
		}
		System.out.println("Matrix is load");
	}

	public void writeMatrix4Mega (String outfileS) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outfileS), 65536);
			bw.write("#mega"); bw.newLine();
			bw.write("!TITLE Distance Matrix from " + this.getMatrixSize() + " Taxa;"); bw.newLine();
			bw.write("!Format DataType=distance;"); bw.newLine(); bw.newLine();
			for (int i = 0; i < this.getMatrixSize(); i++) {
				bw.write("#"+String.valueOf(i)); bw.newLine();
			}
			bw.newLine();
			for (int i = 1; i < this.getMatrixSize(); i++) {
				StringBuilder sb = new StringBuilder("  ");
				for (int j = 0; j < i; j++) {
					String valueStr = String.valueOf((value[i][j]))+"00000";
					sb.append(valueStr.substring(0, 5)).append("   ");
				}
				bw.write(sb.substring(0, sb.length()-3)); bw.newLine();
			}
			bw.flush();
			bw.close();
		}
		catch (Exception e) {
			System.out.println("Error occurred while writing Mega distace file " + outfileS);
		}
	}
	public void writeMatrix (String outfileS, Boolean ifBinary) {
		if (ifBinary) {
			try {
				DataOutputStream dos =new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfileS),65536));
				dos.writeInt(value.length);
				for (int i = 0; i < value.length; i++) {
					for (int j = 0; j < value.length; j++) {
						dos.writeFloat(value[i][j]);
					}
				}
				dos.flush();
				dos.close();
			}
			catch (Exception e) {
				System.out.println("Error occurred while writing Matrix " + outfileS + " " + e.toString());
			}
		}
		else {
			try {
				BufferedWriter bw = new BufferedWriter(new FileWriter(outfileS), 65536);
				for (int i = 0; i < value.length; i++) {
					for (int j = 0; j < value.length; j++) {
						bw.write(String.valueOf(value[i][j])+"\t");
					}
					bw.newLine();
				}
				bw.flush();
				bw.close();
			}
			catch (Exception e) {
				System.out.println("Error occurred while writing Matrix " + outfileS + " " + e.toString());
			}
		}
		System.out.println("Full matrix is written");
	}
}
