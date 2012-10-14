/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.genome.GBS;

import java.io.File;

/**
 *
 * @author fl262
 */
public class GBSPipelineMoira {
	public static String parentDir = "E:/Research/SwitchGrass/GBS/";

	public static String[] childDir = {"Genome", "DigestedGenome", "Solexa", "Key", "ParsedTaxaDir", "CollapsedReads", "CombinedReads", "ReadsByTaxa", "ReadsByTaxaMin",
										"ReadsWorthPhysicalMap", "CommonReadsOnMap", "NormalReadsByTaxa"};
	File[] FileOfData = new File[childDir.length];

	public GBSPipelineMoira () {
		ini();
		pipe();
		System.gc();
	}
	public GBSPipelineMoira(String dir) {
		if (!dir.endsWith("/")) {
			dir = dir + "/";
		}
		parentDir = dir;
		ini();
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

		if (FileOfData[0].list().length == 0 || FileOfData[2].list().length == 0 || FileOfData[3].list().length == 0) {
			System.out.println("Put Sorghum Genome File into" + FileOfData[0]);
			System.out.println("Put Sequencing File into" + FileOfData[1]);
			System.out.println("Put Key File into" + FileOfData[2]);
			System.exit(1);
		}
	}

	void pipe() {
		String Genome = FileOfData[0].toString() + "/" + FileOfData[0].list()[0];
		String DigestedGenome = FileOfData[1].toString() + "/" + childDir[1] + ".txt";
		String Solexa = FileOfData[2].toString();
		String Key = FileOfData[3].toString() + "/" + FileOfData[3].list()[0];
		String ParsedTaxaDir = FileOfData[4].toString();
		String CollapsedReads = FileOfData[5].toString();
		String ConbinedReads = FileOfData[6].toString() + "/" + childDir[6] + ".txt";
		String ReadsByTaxa = FileOfData[7].toString() + "/" + childDir[7] + ".txt";
		String ReadsByTaxaMin = FileOfData[8].toString() + "/" + childDir[8] + ".txt";
		String ReadsWorthPhysicalMap = FileOfData[9].toString() + "/" + childDir[9] + ".txt";
		String CommonReadsOnMap = FileOfData[10].toString() + "/" + childDir[10] + ".txt";
		String NormalReadsByTaxa = FileOfData[11].toString() + "/" + childDir[11] + ".txt";
		
		if (FileOfData[1].list().length == 0) {
			VirtualDigester vd = new VirtualDigester(new File(Genome), new File(DigestedGenome));
			System.gc();

			//sort DigestedGenome
			ReadsWPhysicalMap rwpmVD = new ReadsWPhysicalMap(DigestedGenome, true);
			rwpmVD.sortTable(true);
			rwpmVD.writeCountFile(new File(DigestedGenome));
			System.gc();
		}

		//"10" means the quanlity value for any basepair should larger(>) than 64+10 (Char Value is used to denote quanlity of sequencing for any basepare)
		//throw away the reads with "." in the length of (64+barcode length)
		ParseBarcodeFiles pbf = new ParseBarcodeFiles(Solexa, Key, ParsedTaxaDir, 10, false);
		System.gc();

		//sort reads in files in ParsedTaxaDir, counting the number of identical haplotypes. "1" means printing out the haplotypes with couting number >= 1, output is collapsedReads (Tags).
		ReadCounts rc1 = new ReadCounts(ParsedTaxaDir, CollapsedReads, 1, true, false, false);
		System.gc();


		//sort all tags after everytime a collapsed files is added, then output all tags and its couting number if counting number >= 2
		CombineReadCounts crc = new CombineReadCounts(CollapsedReads, ConbinedReads, 2, true);
		System.gc();


		//input ReadsWorthPhysicalMap and CollapasedReads, search each file of CollapsedReads in ReadsWorthPhysicalMap.
        //Then output the counting number of each reads (in ReadsWorthPhysicalMap) in each taxa
		//Note: "theTagMatchFinder=new TagMatchFinder(haplotype)" is useless here.
		CreateReadsByTaxaPlasmid crbt = new CreateReadsByTaxaPlasmid(ConbinedReads, CollapsedReads, ReadsByTaxa, true);
		System.gc();

		//4 stands for that 4 barcodes are used for 1 individual, this step merge multiple barcodes ReadsByTaxa file into 1
		//minCount = 5, output the reads with >= minCount

		mergeReadsByTaxa mRBT = new mergeReadsByTaxa(ReadsByTaxa, true, 4);
		mRBT.writeNewDistFile(new File(ReadsByTaxaMin), true, 5);
		mRBT.writeReadCountFile(new File(ReadsWorthPhysicalMap), true, 5);
		System.gc();


		//get total read number of each taxa from ParsedTaxaDir, then get ReadsByTaxa format files from ReadsByTaxaMin, then get normalized frequency of each haplotypes
		new normalReadsByTaxaMoira(ParsedTaxaDir, ReadsByTaxaMin, NormalReadsByTaxa, true);
/*
		//creat a empty table of each Tag with no mapping information (chr, positon...)
		ReadCounts rc3 = new ReadCounts(ReadsWorthPhysicalMap, true);
		ReadsWPhysicalMap rwpm1 = new ReadsWPhysicalMap(rc3);
		rwpm1.writeCountFile(new File(CommonReadsOnMap));

		//process DigestedGenome(2 longs) to lookup[2][size*4], that means storing 2 longs into 4 ints for quicksort;
		//Lookup[0] is the int form of sequence, lookup[1] is the index of tags of sorted DigestedGenome
		//Then sort the Lookup array, throw out int sequences with more than 1000 replicates by assignning a MAX_VALUE to that int sequence
		//Search CommonReadsOnMap in DigestedGenome, if perfect match (one locus or multiple loci), then assignning values to correspongding fields
		//if not (divergence neocleatide <= 3), make a TreeMap (Key = index of DigestedGenome, Value = divergence number), then assignning values to correspongding fields;
		//with regard to the Reads with divergence > 3, keep the same as they were
		//rewrite the CommonReadsOnMap (update with new mapping information)
		ReadsWPhysicalMap rwpm2 = MapHaplotypes.blastUsingPhysicalMap(new ReadsWPhysicalMap(DigestedGenome, true), new ReadsWPhysicalMap(CommonReadsOnMap, true), 3);
		rwpm2.writeCountFile(new File(CommonReadsOnMap));
 *
 */
	}

	public static void main(String[] args) {
		if (args.length != 0) {
			GBSPipelineMoira Pipe = new GBSPipelineMoira (args[0]);
		}
		else {
			GBSPipelineMoira Pipe = new GBSPipelineMoira ();
		}
	}
}

class mergeReadsByTaxa extends ReadsByTaxaMoira {
	int newTaxaNum;
    
	String[] newTaxaNames;
	public mergeReadsByTaxa (String infile, boolean binary, int group) {
		super(infile, binary);
		merge(group);
	}
	public void merge (int group) {
	   newTaxaNum = taxaNum / group;
	   newTaxaNames =new String[newTaxaNum];
       newHapDist=new int[haplotypeNum][newTaxaNum];
	   for (int i = 0; i < newTaxaNum; i++) {
			newTaxaNames[i] = taxaNames[i*group];
	   }
	   for (int i = 0; i < haplotypeNum; i++) {
			for (int j = 0; j < newTaxaNum; j++) {
				for (int k = 0; k < group; k++) {
					newHapDist[i][j] =  (newHapDist[i][j] + hapDist[i][j * group + k]);
				}
			}
	   }
	   taxaNum = newTaxaNum;
	   taxaNames = newTaxaNames;
	}
}