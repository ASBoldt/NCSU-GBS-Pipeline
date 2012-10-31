/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS.plasmid;

import java.io.*;

/**
 *
 * @author fl262
 */
public class updateReadCountPipe {
	File parsedTaxaDir = new File("M:/Plastid_GBS/ParsedTaxaDir");
	File tbtFile = new File ("M:/Plastid_GBS/TBTtxt/TBTtxt.txt");
	public updateReadCountPipe (String inFileDirS, String outFileDirS) {
		File[] popDirs = new File (inFileDirS).listFiles();
		for (int i = 0; i < popDirs.length; i++) {
			update (popDirs[i], outFileDirS);
		}
	}
	public void update (File popDir, String outFileS) {
		File[] parsedTaxaFiles = parsedTaxaDir.listFiles();
		for (int i = 0; i < parsedTaxaFiles.length; i++) {
			parsedTaxaFiles[i].delete();
		}
		String popName = popDir.getName();
		outFileS = outFileS.replaceFirst("Z...", popName) + "read_counts_pipe.txt";
		File[] parsedFiles = popDir.listFiles();
		for (int i = 0; i < parsedFiles.length; i++) {
			File path = new File(parsedTaxaDir, parsedFiles[i].getName());
			copyFile (parsedFiles[i], path);
		}
		new plastidPipe();
		File outFile = new File (outFileS);
		if (outFile.exists()) {
			outFile.delete();
		}
		if (tbtFile.renameTo(new File(outFileS))) {
			System.out.println(outFileS + " is updated");
		}	
	}
	public void copyFile (File sourceFile, File desFile) {
		try {
			DataInputStream dis = new DataInputStream (new BufferedInputStream(new FileInputStream(sourceFile), 65536*256));
			DataOutputStream dos = new DataOutputStream (new BufferedOutputStream (new FileOutputStream(desFile), 65536*256));
			int a;
			while ((a = dis.read()) != -1) {
				dos.write(a);
			}
			dos.flush();
			dos.close();
			dis.close();
		}
		catch (Exception e) {
			System.out.println("Error occurred while copying file " + e.toString());
		}
	}
	public static void main (String[] args) {
		String inFileDirS = "M:/Database/NAM_parsed_76bp_raw_pop/pop/";
		String outFileDirS = "E:/Research/plasmid/pipe_mapping/Z001/phenotype/count_position/";
		new updateReadCountPipe (inFileDirS, outFileDirS);
	}
}
