/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS.switchgrass;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;

/**
 *
 * @author Fei
 */
public class HapMap {
	String[] header;
	Record[] hr;

	public HapMap (String infileS) {
		this.readHapMap(infileS);
	}

	public void readHapMap (String infileS) {
		try {
			ArrayList<Record> hrList = new ArrayList();
			BufferedReader br = new BufferedReader(new FileReader(infileS), 65536);
			String temp;
			temp = br.readLine();
			header = temp.split("\t");
			while ((temp = br.readLine()) != null) {
				 hrList.add(new Record(temp.split("\t")));
			}
			hr = hrList.toArray(new Record[hrList.size()]);
		}
		catch (Exception e) {
			System.out.println("Error occurred while reading infileS " + infileS + " " + e.toString());
		}
	}

	public void writeHapMap (String outfileS) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outfileS), 65536);
			StringBuilder sb = new StringBuilder();
			for (int i = 0; i < header.length; i++) {
				sb.append(header[i]).append("\t");
			}
			sb.deleteCharAt(sb.length()-1);
			bw.write(sb.toString());
			bw.newLine();
			for (int i = 0; i < this.hr.length; i++) {
				bw.write(hr[i].getHamMapOutput());
				bw.newLine();
			}
			bw.flush();
			bw.close();
		}
		catch (Exception e) {
			System.out.println("Error occurred while writing outfileS " + outfileS + " " + e.toString());
		}
	}

	public void writeMMC (String outfileS) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outfileS), 65536);
			bw.write("Example,");
			StringBuilder sb = new StringBuilder();
			for (int i = 11; i < this.header.length; i++) {
				sb.append(header[i]).append(",");
			}
			sb.deleteCharAt(sb.length()-1);
			bw.write(sb.toString());
			bw.newLine();
			for (int i = 0; i < hr.length; i++) {
				bw.write(hr[i].getMMCOutput());
				bw.newLine();
			}
			bw.flush();
			bw.close();
		}
		catch (Exception e) {
			System.out.println("Error occurred while writing infileS " + outfileS + " " + e.toString());
		}
	}

	class Record {
		String[] item;

		Record (String[] item) {
			this.item = item;
		}

		String getHamMapOutput () {
			StringBuilder sb = new StringBuilder();
			for (int i = 0; i < item.length; i++) {
				sb.append(item[i]).append("\t");
			}
			sb.deleteCharAt(sb.length()-1);
			return sb.toString();
		}

		String getMMCOutput () {
			StringBuilder sb = new StringBuilder();
			sb.append(item[0]).append(",");
			byte[] genoValue = this.getGenoValue();
			for (int i = 0; i < genoValue.length; i++) {
				sb.append(String.valueOf(genoValue[i])).append(",");
			}
			sb.deleteCharAt(sb.length()-1);
			return sb.toString();
		}

		byte[] getGenoValue() {
			byte[] genoValue = new byte[item.length-11];
			String[] alleles = item[1].split("/");
			int cnt1 = 0, cnt2 = 0;
			for (int i = 11; i < item.length; i++) {
				if (item[i].equals(alleles[0])) {
					genoValue[i-11] = 2;
					cnt1 += 2;
				}
				else if (item[i].equals(alleles[1])) {
					genoValue[i-11] = 0;
					cnt2 += 2;
				}
				else {
					if (item[i].equals("N")) {
						genoValue[i-11] = 1;
					}
					else {
						genoValue[i-11] = 1;
						cnt1++;
						cnt2++;
					}
				}
			}
			if (cnt1 < cnt2) {
				for (int i = 0; i < genoValue.length; i++) {
					if (genoValue[i] == 2) {
						genoValue[i] = 0;
					}
					else if (genoValue[i] == 0) {
						genoValue[i] = 2;
					}
				}
			}
			return genoValue;
		}
	}
}
