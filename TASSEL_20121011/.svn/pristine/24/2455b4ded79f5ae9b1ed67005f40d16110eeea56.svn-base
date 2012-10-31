
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.genome.GBS.switchgrass;


import java.io.File;
import net.maizegenetics.genome.GBS.CombineReadCounts;

/**
 *
 * @author fl262
 */
public class qseqToHapMapNoAnchor {
	public static String parentDir = "M:/HapMapGenerator/";
	public static String[] childDir = {"Qseq", "Key", "ParsedTaxaDir", "CollapsedTagsDir", "CombinedTags", "TagsByTaxaByte",
									 "HapMap", "SimilarityMatirx", "DistanceMatrix", "NJTree", "LinkageGroup", "CollapsedParents"};
	File[] FileOfData = new File[childDir.length];

	public qseqToHapMapNoAnchor () {
		ini();
		//this.TBTPipe();
		this.HapMapPipe();
		//this.PseudoTestCrossHapMapPipe();
		//this.LinkageGroupPipe();
		//pipe();//deprecated
	}

	void ini() {
		for (int i = 0; i < childDir.length; i++) {
			FileOfData[i] = new File(parentDir, childDir[i]);
			if (!FileOfData[i].isDirectory()) {
				FileOfData[i].mkdirs();
			}
		}
		this.checkFile();
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

	public void TBTPipe () {
		String Qseq = FileOfData[0].toString();
		String Key = FileOfData[1].toString() + "/" + FileOfData[1].list()[0];
		String ParsedTaxaDir = FileOfData[2].toString();
		String CollapsedTagsDir = FileOfData[3].toString();
		String ConbinedTags = FileOfData[4].toString() + "/" + childDir[4] + ".bin";
		String TagsByTaxaByte = FileOfData[5].toString() + "/" + childDir[5] + ".bin";
//**************************
//Make TBT
		//ParseBarcodeFiles pbf = new ParseBarcodeFiles(Qseq, Key, ParsedTaxaDir, 0, false); //for Qseq
		//ParseBarcodeFiles pbf = new ParseBarcodeFiles(Qseq, Key, ParsedTaxaDir, -100, true);  //for fastq
		System.gc();

		//ReadCounts rc1 = new ReadCounts(ParsedTaxaDir, CollapsedTagsDir, 1, true, false, true);
		System.gc();

		CombineReadCounts crc = new CombineReadCounts(CollapsedTagsDir, ConbinedTags, 10, true);
		System.gc();

		CreateReadsByTaxa crbt = new CreateReadsByTaxa(ConbinedTags, CollapsedTagsDir, TagsByTaxaByte, true);
		System.gc();
//**************************
	}

	public void HapMapPipe() {
		//String TagsByTaxaByte = FileOfData[5].toString() + "/" + childDir[5] + ".bin";
		String TagsByTaxaByte = FileOfData[5].toString() + "/" + "TagsByTaxaByte_asso.bin";
		//String TagsByTaxaByte = FileOfData[5].toString() + "/" + "TagsByTaxaByte_pop2.bin";
		String HapMap = FileOfData[6].toString() + "/" + childDir[6] + ".txt";
		

//**************************
//Make HapMap
		ReadsByTaxa rbt = new ReadsByTaxa(TagsByTaxaByte, true);
		clusters cls = new clusters(rbt);
		cls.networkFilter(); //Really powerful, which generates 85% single locus SNPs in 282 maize lines
		cls.alleleFrequencyFileter(rbt, (double)0.05, (double)0.5); //Only for linkage pop which has allele frequency peaks
		//cls.heteozygoteFilter(rbt); //Seems not useful, need to be refined or removed
		cls.writeHapMap(rbt, HapMap, (float)0.3, (float)1);
		//cls.hapDetail(rbt, "M:/a.txt", (float)0.5);
		cls.writeFastA(rbt, "M:/a.txt");

		
//**************************


	}

	public void PseudoTestCrossHapMapPipe () {
		//String TagsByTaxaByte = FileOfData[5].toString() + "/" + childDir[5] + ".bin";
		String TagsByTaxaByte = FileOfData[5].toString() + "/" + "TagsByTaxaByte_pop1_linkage.bin";
		String HapMap = FileOfData[6].toString() + "/" + childDir[6] + ".txt";
		String Parent1 = FileOfData[11].toString() + "/" + "U518_merged.cnt";
		String Parent2 = FileOfData[11].toString() + "/" + "U418_merged.cnt";
		String TestCrossHapMap = FileOfData[6].toString() + "/" +"TestCrossHapMap.txt";

		ReadsByTaxa rbt = new ReadsByTaxa (TagsByTaxaByte, true);
		clusters cls = new clusters(rbt);
		cls.networkFilter();
		cls.alleleFrequencyFileter(rbt, (double)0.2, (double)0.3);
		ReadCounts homoRc = new ReadCounts (Parent1, true);
		ReadCounts heteoRc = new ReadCounts (Parent2, true);
		cls.pseudoTestCrossFilter(rbt, homoRc, heteoRc);
		cls.writeHapMap(rbt, TestCrossHapMap, (float)0.9, (float)1);
		//cls.hapDetail(rbt, "M:/detail.txt", (float)0.9, (float)1);
	}

	public void LinkageGroupPipe () {
		String SimilarityMatrix = FileOfData[7].toString() + "/" + childDir[7] + ".bin";
		String DistanceMatrix = FileOfData[8].toString() + "/" + childDir[8] + ".bin";
		String LinkageGroup =  FileOfData[10].toString() + "/" + childDir[10] + ".csv";

//**************************
//Make linkage groups by MMC
		RelationUtil ru = new RelationUtil("M:/HapMapGenerator/HapMap/TestCrossHapMap.txt");
		ru.assignGenoValueMMC();
		ru.writeGenoValue4MMC(LinkageGroup);
//**************************

//**************************
//Make ordered linkage groups, deprecated
//	Make distance matrix for NJ tree
		//RelationUtil ru = new RelationUtil("M:/HapMapGenerator/HapMap/TestCrossHapMap_U518_homo_0.90_0.97_0.23_0.27.txt");
		//RelationUtil ru = new RelationUtil("M:/HapMapGenerator/HapMap/TestCrossHapMap_U518_homo_0.23_0.27.txt");
		//ru.assignGenoValueSimple();
		//System.out.println(ru.permutation(100, (float)0.0001));
		//ru.getRsquareMatrix(ru.getMarkerCount());
		
		//Matrix dMatrix = ru.convert2DMatrixFromRsquareMatrix(); //From permutation, when FDR = 10-6, r2 = 0.2500
		//dMatrix.writeMatrix4Mega("M:/tree.meg");
		//dMatrix.writeMatrix(DistanceMatrix, true);

//	Make linkage groups by NJ tree
		//Matrix dMatrix2 = new Matrix(DistanceMatrix, true);
		//LinkageGroupUtils theLGU = new LinkageGroupUtils(dMatrix2);
		//theLGU.orderLinkageGroup(new Matrix(SimilarityMatrix, true), true);
		//theLGU.writeLinkageGroup(LinkageGroup);

//	Assembly the linkage groups
		//LinkageGroupUtils theLGU = new LinkageGroupUtils(LinkageGroup);
		//theLGU.deleteSmallSizeGroups(10);
		//theLGU.sortLinkageGroupBySize();
		//Matrix dMatrix = new Matrix(SimilarityMatrix, true);
		//Matrix m = theLGU.lGroup[theLGU.lGroup.length-1].getR2MatrixInGroup(dMatrix);
		//m.writeMatrix("M:/m.txt", false);		
//**************************
	}

	public void pipe() { //deprecated
		String Qseq = FileOfData[0].toString();
		String Key = FileOfData[1].toString() + "/" + FileOfData[1].list()[0];
		String ParsedTaxaDir = FileOfData[2].toString();
		String CollapsedTagsDir = FileOfData[3].toString();
		String ConbinedTags = FileOfData[4].toString() + "/" + childDir[4] + ".bin";
		String TagsByTaxaByte = FileOfData[5].toString() + "/" + childDir[5] + ".bin";
		String HapMap = FileOfData[6].toString() + "/" + childDir[6] + ".txt";
		String SimilarityMatirx = FileOfData[7].toString() + "/" +childDir[7] + ".bin";
		String DistanceMatrix = FileOfData[8].toString() + "/" + childDir[8] + ".bin";
		String MegaMatrix = FileOfData[8].toString() + "/" + childDir[8] + ".meg";
		String NJTreeNH = FileOfData[9].toString() + "/" +childDir[9] + ".txt";
		String LinkageGroup =  FileOfData[10].toString() + "/" +childDir[10] + ".bin";

//**************************
//Make TBT
		//ParseBarcodeFiles pbf = new ParseBarcodeFiles(Qseq, Key, ParsedTaxaDir, 0, false);
		System.gc();

		//ReadCounts rc1 = new ReadCounts(ParsedTaxaDir, CollapsedTagsDir, 1, true, false, true);
		System.gc();

		//CombineReadCounts crc = new CombineReadCounts(CollapsedTagsDir, ConbinedTags, 10, true);
		System.gc();

		//CreateReadsByTaxa crbt = new CreateReadsByTaxa(ConbinedTags, CollapsedTagsDir, TagsByTaxaByte, true);
		System.gc();
//**************************

//**************************
//Make HapMap
		//ReadsByTaxa rbt = new ReadsByTaxa(TagsByTaxaByte, true);
		//clusters cls = new clusters(rbt);
		//cls.networkFilter(); //Really powerful, which generates 85% single locus SNPs in 282 maize lines
		//cls.alleleFrequencyFileter(rbt, (double)0.13, (double)0.36); //Only for linkage pop which has allele frequency peaks
		//cls.heteozygoteFilter(rbt); //Seems not useful, need to be refined or removed
		//cls.writeHapMap(rbt, HapMap, (float)0.9);
//**************************

//**************************
//Make phylogeny

		ReadsByTaxa rbt = new ReadsByTaxa(TagsByTaxaByte, true);
		clusters cls = new clusters(rbt);
		cls.networkFilter(); //Really powerful, which generates 85% single locus SNPs in 282 maize lines
		cls.alleleFrequencyFileter(rbt, (double)0.05, (double)0.51); //Only for linkage pop which has allele frequency peaks
		//cls.heteozygoteFilter(rbt); //Seems not useful, need to be refined or removed
		cls.writeHapMap(rbt, HapMap, (float)0.5, (float)1);
		cls.hapDetail(rbt, "M:/hapdetail.txt", (float)0.5, (float)1);


//**************************
//**************************
//Make ordered linkage groups
		//RelationUtil ru = new RelationUtil("M:/HapMapGenerator/HapMap/very_small_set.txt");
		//ru.getRsquareMatrix(ru.genoValue.length);
		//ru.rsquareMatrix.writeMatrix(SimilarityMatirx, true);
		//Matrix dMatrix = ru.convert2DMatrixFromRsquareMatrix((float)0.2267); //From permutation, when FDR = 10-6, r2 = 0.2267
		//dMatrix.writeMatrix(DistanceMatrix, true);

		//Matrix dMatrix2 = new Matrix(DistanceMatrix, true);
		//LinkageGroupUtils theLGU = new LinkageGroupUtils(dMatrix2);
		//theLGU.orderLinkageGroup(new Matrix(SimilarityMatirx, true), true);
		//theLGU.writeLinkageGroup(LinkageGroup);
//**************************


		//LinkageGroupUtils theLGU2 = new LinkageGroupUtils(LinkageGroup);
		//theLGU2.deleteSmallSizeGroups(8);
		//theLGU2.sortLinkageGroupBySize();
		//Matrix r2 = new Matrix(SimilarityMatirx, true);





		//Matrix m = theLGU2.lGroup[theLGU2.getLinkageGroupNumber()-1].getR2MatrixInGroup(r2);
		//Matrix m = theLGU2.lGroup[1].getR2MatrixInGroup(r2);
		//m.writeMatrix("M:/m.txt", false);

		//System.out.println(ru.permutation(200, (double)0.000001));
		//ru.getSigCorrelationMatrix(0.2267);
		//ru.writeCorrelationMatrix(SimilarityMatirx);
		//RelationUtil cu2 = new RelationUtil(SimilarityMatirx, true);
		//cu2.doubleMatrix();
		//cu2.getEdgeDistribution("M:/dis.txt");
		//cu2.edgeNumberFilter(5, 50);
		//LinkageGroupFinder lgf = new LinkageGroupFinder(cu2.doubleRsquareMatrix);
		//lgf.greedyFinder(3);
		//lgf.writeLinkageGroup("M:/linkage.txt");
		
	}
	

	public static void main(String[] args) {
		qseqToHapMapNoAnchor Pipe = new qseqToHapMapNoAnchor ();
	}
}
