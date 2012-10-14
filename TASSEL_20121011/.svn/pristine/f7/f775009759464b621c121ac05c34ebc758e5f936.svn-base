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
public class blastRecArray {
	String header;
	int recNum;
	blastRec[] recs;
	public blastRecArray (File sourceFile) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(sourceFile), 65536);
			header = br.readLine();
			ArrayList<String> fileLine = new ArrayList<String>();
			String temp = br.readLine();
			while (temp != null) {
				fileLine.add(temp);
				temp = br.readLine();
			}
			br.close();
			String[] lineArray = fileLine.toArray(new String[fileLine.size()]);
			recNum = lineArray.length;
			recs = new blastRec[recNum];
			for (int i = 0; i < recNum; i++) {
				recs[i] = new blastRec(lineArray[i]);
			}

		} catch (Exception e) {
			System.out.println(sourceFile);
		}
	}
	public blastRecArray (blastRec[] br, String header) {
		this.header = header;
		this.recs = br;
		recNum = br.length;
	}
	public void writeArray (String desFileS) {
		try {
			BufferedWriter bw = new BufferedWriter (new FileWriter(desFileS), 65536);
			bw.write(header);
			for (int i = 0; i < recNum; i++ ) {
				for (int j = 0; j < recs[0].item.length - 1; j++) {
					bw.write(recs[i].item[j]);
					bw.write("\t");
				}
				bw.write(recs[i].item[recs[0].item.length - 1]);
				bw.newLine();
			}
			bw.flush();
			bw.close();
		}
		catch (Exception e) {
			System.out.println (e.toString());
			System.out.println ("Error occured in " + desFileS);
		}
	}
}
class blastRec implements Comparable<blastRec> {
	String[] item;
	blastRec (String line) {
		item = line.split("\t");
	}
	public int compareTo(blastRec o) {
		return item[0].compareTo(o.item[0]);
	}
}
