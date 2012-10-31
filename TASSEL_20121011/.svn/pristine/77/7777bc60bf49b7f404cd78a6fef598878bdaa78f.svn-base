
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
public class GBSPipelineFei {
	public static String parentDir = "M:/GBS/";
	public static String[] childDir = {"Genome", "DigestedGenome", "Solexa", "Key", "ParsedTaxaDir", "CollapsedReads", "CombinedReads",
		"ReadsWorthPhysicalMap", "ReadsByTaxa", "ReadsByTaxaMin", "CommonReadsOnMap", "GeneticMap", "CommonReadsWithGenetic", "NormalReadsByTaxa"};
	File[] FileOfData = new File[childDir.length];

	public GBSPipelineFei () {
		ini();
		pipe();
		System.gc();
	}
	public GBSPipelineFei(String dir) {
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

		if (FileOfData[0].list().length == 0 || FileOfData[2].list().length == 0 || FileOfData[3].list().length == 0 || FileOfData[11].list().length == 0) {
			System.out.println("Put Maize Genome File into " + FileOfData[0]);
			System.out.println("Put Sequencing File into " + FileOfData[2]);
			System.out.println("Put KeyFile into " + FileOfData[3]);
			System.out.println("Put Genetic Map File into " + FileOfData[11]);
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
		String ReadsWorthPhysicalMap = FileOfData[7].toString() + "/" + childDir[7] + ".txt";
		String ReadsByTaxa = FileOfData[8].toString() + "/" + childDir[8] + ".txt";
		String ReadsByTaxaMin = FileOfData[9].toString() + "/" + childDir[9] + ".txt";
		String CommonReadsOnMap = FileOfData[10].toString() + "/" + childDir[10] + ".txt";
		String GeneticMap = FileOfData[11].toString() + "/" + FileOfData[11].list()[0];
		String CommonReadsWithGenetic = FileOfData[12].toString() + "/" + childDir[12] + ".txt";
		String NormalReadsByTaxa = FileOfData[13].toString() + "/" + childDir[13] + ".txt";
		
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

		//take virtual degisted genome as a base to get tags and its couting number(counting number for all tags in digested genome=2),
		//then input collapsedReads files to get tags and its couting number
		//sort all tags after everytime a collapsed files is added, then output all tags and its couting number if counting number >= 2
		ReadsWPhysicalMap theRef = new ReadsWPhysicalMap(DigestedGenome, true);
		theRef.sortTable(true);
		CombineReadCounts crc = new CombineReadCounts(theRef, CollapsedReads, ConbinedReads, 2, true);
		System.gc();

		//2 means the tags with counting number >=2 are chosen
//		ReadCounts rc2 = new ReadCounts(ConbinedReads, ReadsWorthPhysicalMap, 2, true, true);
//		System.gc();

		//input ReadsWorthPhysicalMap and CollapasedReads, search each file of CollapsedReads in ReadsWorthPhysicalMap.
		//Then output the counting number of each reads (in ReadsWorthPhysicalMap) in each taxa
		//Note: "theTagMatchFinder=new TagMatchFinder(haplotype)" is useless here.
		CreateReadsByTaxa crbt = new CreateReadsByTaxa(ConbinedReads, CollapsedReads, ReadsByTaxa, true);
		System.gc();

		//input ReadsByTaxa file, if total counting number for a Tag in all Taxa >=5 (MinCount), output the ReadsByTaxaMin and
		//rewrite the ReadsWorthPhysicalMap with MinCount >= 5
		ReadsByTaxa theRBT = new ReadsByTaxa(ReadsByTaxa, true);
		theRBT.writeDistFile(new File(ReadsByTaxaMin), true, 5);
		theRBT.writeReadCountFile(new File(ReadsWorthPhysicalMap), true, 5);
		System.gc();

		//get total read number of each taxa from ParsedTaxaDir, then get ReadsByTaxa format files from ReadsByTaxaMin, then get normalized frequency of each haplotypes
		//new normalReadsByTaxa(ParsedTaxaDir, ReadsByTaxaMin, NormalReadsByTaxa, true);

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

		ReadsByTaxa theRBT2 = new ReadsByTaxa(ReadsByTaxaMin, true);
		System.gc();
		System.out.println("rwpmTest2 memory" + rwpm2);

		//Check CommonReadsOnMap with genetic map to refine reads mapping
		//Select reads in CommenReadsOnMap with chromsome information (actually means perfect match to DigestedGenome), moreover the chromsome number should be 1-10
		//According to positon information of reads, find two flanking two markers
		//In all taxa appeared in ReadsByTaxaMin, compute proportion of crossover and non-crossover recombinants to get "P" and "Q" for a Bianomial distribution
		//In taxa with counting number > 0, compute proportion of crossover and non-crossover recombinants to test if the read is mapped corretly, process:
		//made a bianomial disstribution with "P", "Q" and computation in last step to have the p value
		//signigicant p value means the read is here
		//similarly, search other loci in genetic map to get the p value, which is used for comparison
		  /*if ((flankP < sigThres) || (bestNonChrP < sigThres)) {
		mapsSomewhere++;
		}
		 *///what does is mean?
		rwpm2 = MapHaplotypes.checkWorkingMapWithGenetic(GeneticMap, rwpm2, theRBT2, true, 0.00001);
		rwpm2.writeCountFile(new File(CommonReadsWithGenetic), Integer.MAX_VALUE, true, true, 0.01f, true);
	}

	public static void main(String[] args) {
		if (args.length != 0) {
			GBSPipelineFei Pipe = new GBSPipelineFei (args[0]);
		}
		else {
			GBSPipelineFei Pipe = new GBSPipelineFei ();
		}
	}
}
