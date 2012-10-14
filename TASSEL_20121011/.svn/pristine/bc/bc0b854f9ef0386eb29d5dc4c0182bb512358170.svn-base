/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS.plasmid;

import java.io.*;
import java.util.Arrays;

/**
 *
 * @author fl262
 */
public class tagCountArrayOperation {
	public static void collapseFile (String inFileS, String outFileS, int minCount) {
		tagCountArray tca = new tagCountArray (inFileS);
		tca.collapseCounts();
		tca.writeTagCount(outFileS, minCount);
		tca = null;
		System.gc();
		System.out.println(inFileS + " is collapsed");
	}
	public static void collapseFilesInDir (String inFileDirS, String outFileDirS, int minCount) {
		File[] inFile = new File(inFileDirS).listFiles();
		for (int i = 0; i < inFile.length; i++) {
			String outFile = outFileDirS + "\\" + inFile[i].getName();
			collapseFile(inFile[i].getPath(), outFile, minCount);
		}
	}
	public static void getRepTag (String inFileS, String outFileS, int minCount) {
		tagCountArray tca = new tagCountArray(inFileS);
		tca.writeTagCount(outFileS, minCount);
		tca = null;
	}
	public static void combineTagCountFiles (String inFileDirS, String outFileS) {
		int maxNumTag = 50000000;
		int currentRow = 0;
		tagCountArray tca = new tagCountArray(maxNumTag);
		File[] inFile = new File(inFileDirS).listFiles();
		for (int i = 0; i < inFile.length; i++) {
			currentRow = tca.mergeTagCount(inFile[i], currentRow);
			System.out.println(inFile[i].toString() + " is combined");
		}
		tca.writeCombinedTagCount(outFileS, currentRow);
		tca = null;
	}

	public static void creatTagByTaxa (String tagCountDirS, String combinedFileS, String parsedTaxaDir, String TBTFileS) {
		File[] tagCountFile = new File(tagCountDirS).listFiles();
		tagCountArray combinedTCA = new tagCountArray(combinedFileS);
		String[] taxaNames = new String[tagCountFile.length];
		int[] totalCountsByTaxa = new int[taxaNames.length];
		short[][] countDist = new short[combinedTCA.totalTag][tagCountFile.length];
		for (int i = 0; i < tagCountFile.length; i++) {
			taxaNames[i] = tagCountFile[i].getName();
			String path = parsedTaxaDir + "/" + taxaNames[i];
			totalCountsByTaxa[i] = (int) (new File(path).length() / (baseEncoder.byteLen + 4));
		}
		for (int i = 0; i < combinedTCA.totalTag; i++) {
			for (int j = 0; j < taxaNames.length; j++) {
				countDist[i][j] = 0;
			}
		}
		for (int i = 0; i < tagCountFile.length; i++) {
			tagCountArray tca = new tagCountArray(tagCountFile[i]);
			for (int j = 0; j < tca.totalTag; j++) {
				int hit = Arrays.binarySearch(combinedTCA.tagCs, tca.tagCs[j]);
				if (hit >= 0) {
					countDist[hit][i] = (short) tca.tagCs[j].count;
				}
			}
		}
		tagByTaxaShort tbt = new tagByTaxaShort (taxaNames, totalCountsByTaxa, combinedTCA, countDist);
		tbt.writeTBTshortTXT(TBTFileS);
	}
	public static void deleteFilesInDir (String inFileDirS) {
		File[] files = new File(inFileDirS).listFiles();
		for (int i = 0; i < files.length; i++) {
			files[i].delete();
		}
	}
	public static void mergeFilesInDir (String inFileDirS) {
		int[] index;
		do {
			File[] inFiles = new File (inFileDirS).listFiles();
			String[] taxaNames = new String[inFiles.length];
			for (int i = 0; i < inFiles.length; i++) {
				taxaNames[i] = inFiles[i].getName().replaceAll("_.+", "");
			}
			index = null;
			for (int i = 0; i < inFiles.length - 1; i++) {
				boolean b = false;
				for (int j = i + 1; j < inFiles.length; j++) {
					if (taxaNames[i].equals(taxaNames[j])) {
						index = new int[2];
						index[0] = i;
						index[1] = j;
						b = true;
						break;
					}
				}
				if (b) break;
			}
			if (index != null) {
				String outFileS = inFiles[index[0]].getParent() + "/" + inFiles[index[0]].getName().replaceFirst("_.+", "_merged.cnt");
				merge2Files (inFiles[index[0]], inFiles[index[1]], new File(outFileS));
			}
		} while (index != null);
	}
	public static void merge2Files (File file1, File file2, File outFile) {
		File tempFile = new File (file1.getParent() + "/temp.cnt");
		try {
			DataInputStream dis = new DataInputStream (new BufferedInputStream(new FileInputStream(file1), 65536));
			DataOutputStream dos = new DataOutputStream (new BufferedOutputStream (new FileOutputStream(tempFile), 65536));
			int a;
			while ((a = dis.read()) != -1) {
				dos.write(a);
			}
			dos.flush();
			dis.close();
			dis =  new DataInputStream (new BufferedInputStream(new FileInputStream(file2), 65536));
			while ((a = dis.read()) != -1) {
				dos.write(a);
			}
			dos.flush();
			dos.close();
			dis.close();
		}
		catch (Exception e) {
			System.out.println("Error occurred while merging 2 files " + e.toString());
		}
		file1.delete();
		file2.delete();
		if (outFile.exists()) {
			if (outFile.equals(file1) || outFile.equals(file2)) {
				tempFile.renameTo(outFile);
			}
			else {
				outFile.renameTo(file1);
				tempFile.renameTo(outFile);
			}
		}
		else {
			tempFile.renameTo(outFile);
		}
	}
}
