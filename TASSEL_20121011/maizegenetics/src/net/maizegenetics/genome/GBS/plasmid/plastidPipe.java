/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS.plasmid;

import java.io.File;

/**
 *
 * @author fl262
 */
public class plastidPipe {
	public static String parentDir = "M:/Plastid_GBS/";
	public static String[] childDir = {"Qseq", "Key", "ParsedTaxaDir", "CollapsedTags", "CombinedTags", "CombinedRepTags",  "TBTtxt"};
	File[] FileOfData = new File[childDir.length];

	public plastidPipe () {
		ini();
		checkFile();
		pipe();
		System.gc();
	}
	void ini() {
		for (int i = 0; i < childDir.length; i++) {
			FileOfData[i] = new File(parentDir, childDir[i]);
			if (!FileOfData[i].isDirectory()) {
				FileOfData[i].mkdirs();
			}
		}
	}
	public void checkFile() {
		if (FileOfData[0].list().length == 0) {
			System.out.println("Put Sequencing File into " + FileOfData[0]);
			System.exit(1);
		}
		if (FileOfData[1].list().length == 0) {
			System.out.println("Put KeyFile into " + FileOfData[1]);
			System.exit(1);
		}
	}
	void pipe() {
		String Qseq = FileOfData[0].toString();
		String Key = FileOfData[1].toString() + "/" + FileOfData[1].list()[0];
		String ParsedTaxaDir = FileOfData[2].toString();
		String CollapsedTags = FileOfData[3].toString();
		String CombinedTags = FileOfData[4].toString() + "/" + childDir[4] + ".bin";
		String CombinedRepTags = FileOfData[5].toString() + "/" + childDir[5] + ".bin";
		
		String TBTshort = FileOfData[6].toString() + "/" + childDir[6] + ".txt";

		//new parseBarcode(Qseq, Key, ParsedTaxaDir, 0);
		System.gc();

		//tagCountArrayOperation.mergeFilesInDir(ParsedTaxaDir);

		
		tagCountArrayOperation.collapseFilesInDir(ParsedTaxaDir, CollapsedTags, 2);
		System.gc();

		tagCountArrayOperation.combineTagCountFiles(CollapsedTags, CombinedTags);
		System.gc();

		tagCountArrayOperation.getRepTag(CombinedTags, CombinedRepTags, 500); //15, select tags greater than 15 in a inbred
		System.gc();

		tagCountArrayOperation.creatTagByTaxa(CollapsedTags, CombinedRepTags, ParsedTaxaDir, TBTshort);


		tagCountArrayOperation.deleteFilesInDir(CollapsedTags);


	}
	public static void main (String[] args) {
		plastidPipe Pipe = new plastidPipe ();
	}
}
