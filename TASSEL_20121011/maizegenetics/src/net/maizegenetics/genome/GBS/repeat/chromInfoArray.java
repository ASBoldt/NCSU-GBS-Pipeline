/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS.repeat;

import java.io.*;
import java.util.ArrayList;

/**
 *
 * @author fl262
 */
public class chromInfoArray {
	File desFile = new File("E:/Database/InfoFile/ChrLenCentPosi.txt");
	int chromNum;
	chromInfo[] chroms;
	public chromInfoArray () {
		readFile(desFile);
	}
	public void readFile(File desFile) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(desFile), 65536);
			br.readLine();
			ArrayList<String> infoLine = new ArrayList();
			String temp = br.readLine();
			while (temp != null) {
				infoLine.add(temp);
				temp = br.readLine();
			}
			br.close();
			String[] infoLineArray = infoLine.toArray(new String[infoLine.size()]);
			chroms = new chromInfo[infoLine.size()];
			chromNum = chroms.length;
			for (int i = 0; i < chromNum; i++) {
				chroms[i] = new chromInfo(infoLineArray[i]);
			}
		}
		catch (Exception e) {
			System.out.println(desFile.toString());
		}
	}
}
class chromInfo {
	int chromID;
	int chromLength;
	int centB, centE;
	public chromInfo (String line) {
		String[] temps = line.split("\t");
		chromID = Integer.valueOf(temps[0]);
		chromLength = Integer.valueOf(temps[1]);
		String[] tems = temps[2].split("\\.\\.");
		centB = Integer.valueOf(tems[0]);
		centE = Integer.valueOf(tems[1]);
	}
}
