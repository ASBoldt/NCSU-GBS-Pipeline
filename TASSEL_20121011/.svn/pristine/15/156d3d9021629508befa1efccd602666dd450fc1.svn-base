/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS.switchgrass;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import net.maizegenetics.genome.BaseEncoder;
import net.maizegenetics.genome.GBS.CombineReadCounts;

/**
 *
 * @author Fei
 */
public class Cake {
	public String parentDir = "M:/Generator/";
	public String[] childDir = {"Illumina", "Key", "ParsedTaxaDir", "CollapsedTagsDir", "CombinedTags", "TagPair", "TagsByTaxaByte", "MapInfo",
									 "HapMap", "TaxaSelection"};
	File[] FileOfData = new File[childDir.length];

	public Cake () {
		this.setupDir();
	}

	public void setupDir () {
		for (int i = 0; i < childDir.length; i++) {
			FileOfData[i] = new File(parentDir, childDir[i]);
			FileOfData[i].mkdirs();
		}
	}

	public void setParameters (String infileS) {

	}

	public void runPipe() {
		//this.MapInfoPipe();
		//this.HapMapPipe();
		//this.LinkageGroupPipe();
		//this.LinkageGroupPipeReed();
		//this.txtTagCount();
		//this.LDMapPipe();
	}

	public void MapInfoPipe() {
		String Qseq = FileOfData[0].toString();
		String Key = FileOfData[1].toString() + "/" + FileOfData[1].list()[0];
		String ParsedTaxaDir = FileOfData[2].toString();
		String CollapsedTagsDir = FileOfData[3].toString();
		String ConbinedTags = FileOfData[4].toString() + "/" + childDir[4] + ".tc";
		String TagPairFile = FileOfData[5].toString() + "/" + childDir[5] + ".tp";
		String TagPairBySeqFile = FileOfData[5].toString() + "/" + childDir[5] + ".tpbs";
		String TagsByTaxaByte = FileOfData[6].toString() + "/" + childDir[6] + ".tbt";
		String HapMapInfo = FileOfData[7].toString() + "/" + childDir[7] + ".bin";

		//ParseBarcodeFiles pbf = new ParseBarcodeFiles(Qseq, Key, ParsedTaxaDir, 0, false); //for Qseq
		//ParseBarcodeFiles pbf = new ParseBarcodeFiles(Qseq, Key, ParsedTaxaDir, -100, true);  //for fastq
		System.gc();

		//ReadCounts rc1 = new ReadCounts(ParsedTaxaDir, CollapsedTagsDir, 1, true, false, true);
		System.gc();

		//CombineReadCounts crc = new CombineReadCounts(CollapsedTagsDir, ConbinedTags, 1, true);
		System.gc();

		//ReadCounts rc = new ReadCounts ("M:/Generator/Collapsed_backup/parents/P1-U518_high.cnt", true);
		//rc.writeCountFile(new File (ConbinedTags), 5, true);

		//ReadCounts rc2 = new ReadCounts (ConbinedTags, true);
		//NetworkFilter nf = new NetworkFilter (rc2, (double)0.03, TagPairFile); //error = 0.01, P(tag(64bp) of 1bp-error) = 0.33977
		System.gc();

		//TagPair tp = new TagPair(TagPairFile, TagPairBySeqFile);


		//CreateReadsByTaxa crbt = new CreateReadsByTaxa(TagPairBySeqFile, CollapsedTagsDir, TagsByTaxaByte, true);
		System.gc();

		//ReadsByTaxa rbt = new ReadsByTaxa (TagsByTaxaByte, true);
		//NoRefHapUtils1 hmu = new NoRefHapUtils1 (rbt, TagPairBySeqFile);
		//hmu.writeMapInfo(HapMapInfo);

	}

	public void HapMapPipe() {
		String HapMapInfo = FileOfData[7].toString() + "/" + childDir[7] + ".bin";
		String HapMapDir = FileOfData[8].toString() + "/";
		
		//NoRefHapUtils1 hmu = new NoRefHapUtils1 (HapMapInfo);
		//NoRefHapUtils1 hmu = new NoRefHapUtils1 ("M:/Generator/MapInfo/MapInfo_asso_min5_e3.bin");
		NoRefHapUtils1 hmu = new NoRefHapUtils1 ("M:/Generator/MapInfo/MapInfo_NAM_pop1_min5_e3.bin");
		//NoRefHapUtils1 hmu = new NoRefHapUtils1 ("M:/Reed/MapInfo/MapInfo_LP1_min5_e3_parents.bin");
		//hmu.selectTaxa("M:/Generator/TaxaSelection/pop1_full_sib.txt");
		hmu.heterozygoteTest((double)3.841); //only if it is a heterozygous species, 10.83 p=0.0001, 6.635 p=0.01, 3.841 p=0.05
		hmu.writeHapAll(HapMapDir, 0.05, 0.5, 0.3); //MAF should always > 0.04

	}

	public void LinkageGroupPipeReed () {
		String HapMapDir = FileOfData[8].toString() + "/";
		String Parent1 = "M:/Reed/Collapsed_backup/parents/RCG-32_C05F2ACXX_merged.cnt";
		String Parent2 = "M:/Reed/Collapsed_backup/parents/RCG-42_C05F2ACXX_merged.cnt";
		NoRefHapUtils1 hmu = new NoRefHapUtils1 ("M:/Reed/MapInfo/MapInfo_LP1_min5_e3.bin");
		//hmu.heterozygoteTest((double)3.841);
		//hmu.shallowPseudoTestCross(Parent2, Parent1);
		//hmu.writeHapAll(HapMapDir, 0.2, 0.3, 0.5);


	}

	public void txtTagCount () {
		String inTagCount = "M:/RCG-45_C05F2ACXX_merged.cnt";
		String outTagCount = "M:/RCG-45_C05F2ACXX_merged.txt";
		ReadCounts tc = new ReadCounts (inTagCount, true);
		try {
			BufferedWriter bw = new BufferedWriter (new FileWriter(outTagCount), 65536);
			for (int i = 0; i < tc.getSize(); i++) {
				long[] tag = new long[2];
				for (int j = 0; j < tag.length; j++) {
					tag[j] = tc.haplotype[j][i];
				}
				bw.write(BaseEncoder.getSequenceFromLong(tag)+"\t"+String.valueOf(tc.hapcount[i]));
				bw.newLine();
			}
			bw.flush();
			bw.close();
		}
		catch (Exception e) {

		}
	}

	public void LinkageGroupPipe () {
		String HapMapDir = FileOfData[8].toString() + "/";
		String Parent1 = "M:/Generator/Collapsed_backup/parents/P1-U518_high_min5.cnt";
		String Parent2 = "M:/Generator/Collapsed_backup/parents/P2-U418_high_min5.cnt";

		NoRefHapUtils1 hmu = new NoRefHapUtils1 ("M:/Generator/MapInfo/MapInfo_pop1_fullsib_min5_e3.bin");
		hmu.heterozygoteTest((double)3.841); //10.83 p=0.0001, 6.635 p=0.01, 3.841 p=0.05
		//hmu.pseudoTestCross(Parent1, Parent2);
		hmu.homoParentCross(Parent1, Parent2);
		hmu.writeHapAll(HapMapDir, 0.45, 0.55, 0.3, 1);

		HapMap hm = new HapMap("M:/Generator/HapMap/HapMap.hmp.txt");
		hm.writeMMC("M:/Generator/MMC/mmc.txt");
	}

	public void LDMapPipe () {
		File[] hapSeqFiles = new File[3];
		hapSeqFiles[0] = new File ("M:/Generator/HapMap/pop1_fullsib_0.05_0.5_0_min5_h_e3/HapMap.fas.txt");
		hapSeqFiles[1] = new File ("M:/Generator/HapMap/pop2_halfsib_0.05_0.5_min5_h_e3/HapMap.fas.txt");
		hapSeqFiles[2] = new File ("M:/Generator/HapMap/asso_0.05_0.5_min5_e3/HapMap.fas.txt");
		File tagCountFolder = new File ("M:/Generator/CollapsedTagsDir");
		String mergedMapInfo = "M:/Generator/MapInfo/Mapinfo_merged.bin";
		String HapMapDir = "M:/Generator/HapMap";
		//NoRefHapUtils1 hmu = new NoRefHapUtils1 ();
		//hmu.mergeMapInfo(hapSeqFiles, tagCountFolder);
		//hmu.writeMapInfo(mergedMapInfo);

		NoRefHapUtils1 hmu1 = new NoRefHapUtils1 (mergedMapInfo);
		hmu1.selectTaxa("M:/Generator/TaxaSelection/4x.txt");
		//hmu1.heterozygoteTest((double)3.841);
		hmu1.writeHapAll(HapMapDir, 0.00, 0.5, 0.00);
	}

	public void runOnParameters () {
		String Qseq = FileOfData[0].toString();
		String Key = FileOfData[1].toString() + "/" + FileOfData[1].list()[0];
		String ParsedTaxaDir = FileOfData[2].toString();
		String CollapsedTagsDir = FileOfData[3].toString();
		String ConbinedTags = FileOfData[4].toString() + "/" + childDir[4] + ".bin";
		String TagsByTaxaByte = FileOfData[5].toString() + "/" + childDir[5] + ".bin";
	}
	
	public static void main (String[] args) {
		Cake c = new Cake();
		c.runPipe();
	}
}
