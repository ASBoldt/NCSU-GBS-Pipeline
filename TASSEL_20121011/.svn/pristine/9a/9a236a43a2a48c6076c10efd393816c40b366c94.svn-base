
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.genome.GBS;

import java.io.File;
import net.maizegenetics.gbs.pipeline.TagHomologyPhaseNoAnchorMT;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.tagdist.TagsByTaxaBitFileMap;
import net.maizegenetics.genome.GBS.switchgrass.byteToBit;

/**
 *
 * @author fl262
 */
public class qseqToHapMapNoAnchor {
	public static String parentDir = "M:/HapMapGenerator/";
	public static String[] childDir = {"Qseq", "Key", "ParsedTaxaDir", "CollapsedTags", "CombinedTags", "TagsByTaxaByte", "TagsByTaxaBit", "HapMap"};
	File[] FileOfData = new File[childDir.length];

	public qseqToHapMapNoAnchor () {
		ini();
		pipe();
		System.gc();
	}
	public qseqToHapMapNoAnchor(String dir) {
		if (!dir.endsWith("/")) {
			dir = dir + "/";
		}
		parentDir = dir;
		ini();
		checkFile ();
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
		String ConbinedTags = FileOfData[4].toString() + "/" + childDir[4] + ".bin";
		String TagsByTaxaByte = FileOfData[5].toString() + "/" + childDir[5] + ".bin";
		String TagsByTaxaBit = FileOfData[6].toString() + "/" + childDir[6] + ".bin";
		String HapMap = FileOfData[7].toString() + "/" + childDir[7] + ".txt";

		//ParseBarcodeFiles pbf = new ParseBarcodeFiles(Qseq, Key, ParsedTaxaDir, 0, false);
		System.gc();

		//ReadCounts rc1 = new ReadCounts(ParsedTaxaDir, CollapsedTags, 1, true, false, false);
		System.gc();

		//CombineReadCounts crc = new CombineReadCounts(CollapsedTags, ConbinedTags, 10, true);
		System.gc();

		//CreateReadsByTaxa crbt = new CreateReadsByTaxa(ConbinedTags, CollapsedTags, TagsByTaxaByte, true);
		System.gc();

		new byteToBit(TagsByTaxaByte, TagsByTaxaBit);
		System.gc();

		TagsByTaxa theTBT = new TagsByTaxaBitFileMap(TagsByTaxaBit);
		new TagHomologyPhaseNoAnchorMT(theTBT, HapMap);
		System.gc();
	}

	public static void main(String[] args) {
		if (args.length != 0) {
			qseqToHapMapNoAnchor Pipe = new qseqToHapMapNoAnchor (args[0]);
		}
		else {
			qseqToHapMapNoAnchor Pipe = new qseqToHapMapNoAnchor ();
		}
	}
}
