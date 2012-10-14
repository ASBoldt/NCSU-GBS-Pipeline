/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.gbs.pav;

import java.io.File;
import net.maizegenetics.gbs.pipeline.MergeTagsByTaxaFilesByRowPlugin;
import net.maizegenetics.gbs.pipeline.MergeTagsByTaxaFilesPlugin;
import net.maizegenetics.gbs.pipeline.TagAgainstAnchor;
import net.maizegenetics.gbs.pipeline.TagCallerAgainstAnchorMT;
import net.maizegenetics.gbs.pipeline.UNetworkFilter;
import net.maizegenetics.gbs.tagdist.TagCountMutable;
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByte;
import net.maizegenetics.gbs.tagdist.UTagCountMutable;

/**
 *
 * @author Fei Lu
 */
public class PAVgo {
	String parentDir;
	String[] childDir = {"tbtByLane", "mergedTBT", "tagCount", "sliceTBT", "anchor", "mappingResult", "tagAlignment", "pedigree", "PhyGenMapping", "cnv"};
	File[] fileOfData;
	
	public PAVgo (String workingDirS) {
		this.parentDir = workingDirS;
		this.creatDir(workingDirS);
		//this.test();
		//this.assoPipe();
		this.jointLinkagePipe();
		//this.cnvPipe();
	}

	public void test () {
		Test t = new Test();
		t.makeNrowOfTBT();
	}

	public void cnvPipe () {
		this.getCNVTagCount();
		this.getTBTOfCNVTag();
		this.assignPosToCnvTBTAndMergeTaxa();
		this.getTaxaTagCountInBin();
	}

	public void jointLinkagePipe () {
//*************PAV identification section******************        
		//this.mergeTBTByRow();

		//this.sliceTBT();

		//this.creatTagPMap();
		
		//this.mapTBTNAM(); //impossible in this computer

		//this.mergeMappingResult();

		//this.creatTagPGMAP();
        
        //this.checkMappingQuality();
        
        //this.filteredTagPGMapByThresh();
        
        //this.statistics();

		//this.alignTagOfPGMap();

		//this.identifyPAV(); 
        
        
//*******************Standerdize data**********************      
		//this.getB73TagPosCount();
	}

	public void assoPipe () {//LRatio > 4, 90% on chromosome, 60% 100kb region, population structure affect the mapping
		//this.mergeTBTByRow();

		//this.sliceTBT();

		this.mapTBT(); //impossible in this computer

		//this.mergeMappingResult();

		//this.creatTagMap();
	}

	public void calCnvRatioByB73 () {
		String cnvPosFileS = new File (fileOfData[9], "/tagRatio/cnv_posi.txt").getAbsolutePath();
		String cnvRatioFileS = new File (fileOfData[9], "/tagRatio/cnv_ratio.txt").getAbsolutePath();

	}

	public void assignPosToCnvTBTAndMergeTaxa () {//just a test, this method consumes too much memory, can't sort
		String cnvTBTFileS = new File (fileOfData[9], "/tbt/cnv.tbt.byte").getAbsolutePath();
		String cnvPosFileS = new File (fileOfData[9], "/tagPos/cnv_pos.txt").getAbsolutePath();
		String tagPMapFileS = new File (fileOfData[6], "26Mtag-k2.pmap").getAbsolutePath();
		String totalCountFileS = new File (fileOfData[9], "/totalCount/taxaCount.txt").getAbsolutePath();
		PAVUtils util = new PAVUtils();
		util.makeCNVPosFile(cnvTBTFileS, cnvPosFileS, tagPMapFileS, totalCountFileS);
	
	}

	public void getTaxaTagCountInBin () {
//		String cnvTBTFileS = new File (fileOfData[9], "/tbt/cnv.tbt.byte").getAbsolutePath();
//		String tagPMapFileS = new File (fileOfData[6], "26Mtag-k2.pmap").getAbsolutePath();
//		String genomeInfoFileS = "E:/Database/InfoFile/ChrLenCentPosi.txt";
		String binaryBinFileS = new File (fileOfData[9], "/tagCountInBin/binCount.bin").getAbsolutePath();
//		TaxaTagCountInBin tagBin = new TaxaTagCountInBin(genomeInfoFileS,cnvTBTFileS, tagPMapFileS, 100000);
//		tagBin.writeBinaryFile(binaryBinFileS);

		String txtBinFileS = new File (fileOfData[9], "/tagCountInBin/binCount.txt").getAbsolutePath();
		TaxaTagCountInBin tagBin2 = new TaxaTagCountInBin(binaryBinFileS);
		tagBin2.mergeTaxa();
		tagBin2.writeTxtFile(txtBinFileS);
	}

	public void getTBTOfCNVTag () {
		String cnvTagCountFileS =  new File (fileOfData[9], "/tagCount/cnv.cnt").getAbsolutePath();
		String mergedTBTFileS = "Z:/tbt_20.byte";
		//String mergedTBTFileS = new File (fileOfData[1], "tbt_test.byte").getAbsolutePath();
		String taxaNameFileS = new File (fileOfData[9], "/taxaName/interet_sample.txt").getAbsolutePath();
		String cnvTBTFileS = new File (fileOfData[9], "/tbt/cnv.tbt.byte").getAbsolutePath();
		PAVUtils util = new PAVUtils();
		util.makeSubTBTByTaxaTagCount(mergedTBTFileS, cnvTagCountFileS, taxaNameFileS, cnvTBTFileS);
	}

	public void getCNVTagCount () {
		String tagPMapFileS = new File (fileOfData[6], "26Mtag-k2.pmap").getAbsolutePath();
		String tagCountFileS = new File (fileOfData[2], "merged_20.cnt").getAbsolutePath();
		String cnvTagCountFileS =  new File (fileOfData[9], "/tagCount/cnv.cnt").getAbsolutePath();
		PAVUtils util = new PAVUtils();
		util.makeCNVTagCountFile(tagPMapFileS, tagCountFileS, cnvTagCountFileS);
	}

	public void getB73TagPosCount () {
		String tagPMapFileS = new File (fileOfData[6], "26Mtag-k2.pmap").getAbsolutePath();
		String B73TagPosCountFileS = "E:/Research/pav/attDis/B73TagPosCount.txt";
        String cnvTBTFileS = new File (fileOfData[9], "/tbt/cnv_NAMParents.tbt.byte").getAbsolutePath();
		PAVUtils util = new PAVUtils();
		util.getB73TagPosCount(tagPMapFileS, cnvTBTFileS, B73TagPosCountFileS);
	}

	public void identifyPAV () {
		String samFileS = new File (fileOfData[6], "NAM.pgmap.fil.very-sensitive-local-k10.sam").getAbsolutePath();
		String tagPGMapFileS = new File (fileOfData[8], "NAM.pgmap.fil.txt").getAbsolutePath();
		String pavPGMapFileS = new File (fileOfData[8], "pav.pgmap.fil.txt").getAbsolutePath();
		TagPGMap pgmap = new TagPGMap (tagPGMapFileS, false);
		boolean[] ifPAV = pgmap.getIfPAV(samFileS);
		pgmap.writeFilteredTxtFile(pavPGMapFileS, ifPAV);
		TagPGMap pgmap2 = new TagPGMap (pavPGMapFileS, false);
		pgmap2.sortByGeneticPosition();
		pgmap2.writeTxtFile("M:/test.txt");
	}

	public void alignTagOfPGMap () {
		String tagPGMapFileS = new File (fileOfData[8], "NAM.pgmap.fil.txt").getAbsolutePath();
		String fastaFileS = new File (fileOfData[8], "NAM.pgmap.fil.fas").getAbsolutePath();
		TagPGMap pgmap = new TagPGMap (tagPGMapFileS, false);
		pgmap.writeFastaFile(fastaFileS);
		//align those tags using bowtie2
	}

    public void statistics () {
         String filteredTagPGMapFileS = new File (fileOfData[8], "NAM.pgmap.fil.txt").getAbsolutePath();
         TagPGMap pgmap = new TagPGMap(filteredTagPGMapFileS, false);
         boolean[] ifB73 = pgmap.getIfRefTag();
         boolean[] ifNonB73 = pgmap.getReverseFilter(ifB73);
         boolean[] ifSingle = pgmap.getIfSingleMapping(0.05);
         boolean[] ifMultiple = pgmap.getIfMultipleMapping(0.05);
         pgmap.getFilterCross(ifB73, ifMultiple);
    }
    
    public void filteredTagPGMapByThresh () {
        double pThresh = 0.0001;
        String tagPGMapFileS = new File (fileOfData[8], "NAM.pgmap.txt").getAbsolutePath();
        String filteredTagPGMapFileS = new File (fileOfData[8], "NAM.pgmap.fil.txt").getAbsolutePath();
        TagPGMap pgmap = new TagPGMap(tagPGMapFileS, false);
        boolean[] ifUnderThresh = pgmap.getIfUnderThresholdEitherFamilyGroup(pThresh);
        pgmap.writeFilteredTxtFile(filteredTagPGMapFileS, ifUnderThresh);
    }
    
    public void checkMappingQuality () {
        String tagPGMapFileS = new File (fileOfData[8], "NAM.pgmap.txt").getAbsolutePath();
        TagPGMap pgmap = new TagPGMap(tagPGMapFileS, false);
        pgmap.checkMappingQuality(0.0001, true, true);
        //pgmap.checkMappingQualityMultiGroup(0.0001, 1);
    }
    
	public void creatTagPGMAP () {
		String tagPMapFileS = new File (fileOfData[6], "26Mtag-k2.pmap").getAbsolutePath();
		String tagGMapFileS = new File (fileOfData[5], "merged.mapping.txt").getAbsolutePath();
		String tagPGMapFileS = new File (fileOfData[8], "NAM.pgmap.txt").getAbsolutePath();
		String tagCountFileS = new File (fileOfData[2], "merged_20.cnt").getAbsolutePath();
		TagPGMap pgmap = new TagPGMap (tagPMapFileS, tagGMapFileS, tagCountFileS);
		pgmap.writeTxtFile(tagPGMapFileS);	
	}

	public void mapTBTNAM () {
		String tbtFileS = new File (fileOfData[1], "tbt_test.byte").getAbsolutePath();
		String anchorMapFileS = new File (fileOfData[4], "testData/HMM/NAM282_20111217_scv10mF8maf002_mgs_E1pLD5_HMM_filtSites_imp_chr+.hmp.txt").getAbsolutePath();
		//String anchorMapFileS = new File (fileOfData[4], "impNAM/NAMc+.hmp.txt").getAbsolutePath();
		String pedigreeFileS = new File (fileOfData[7], "NAM.txt").getAbsolutePath();
		String mappingResultFileS = "M:/mapR.txt";
		TagAgainstAnchorNAM taa = new TagAgainstAnchorNAM (tbtFileS, anchorMapFileS, pedigreeFileS, mappingResultFileS);

	}

	public void mapTBTSorghum () {
		String tbtFileS = "M:/sorghum/test.tbt.byte";
		String anchorMapFileS = "M:/sorghum/chr+.hmp.txt";
		//String anchorMapFileS = new File (fileOfData[4], "impNAM/NAMc+.hmp.txt").getAbsolutePath();
		String pedigreeFileS = "M:/sorghum/SbPedigree.txt";
		String mappingResultFileS = "M:/mapR.txt";
		TagAgainstAnchorBiParental taa = new TagAgainstAnchorBiParental (tbtFileS, anchorMapFileS, pedigreeFileS, mappingResultFileS);

	}

	public void selectNamFromAllAnchor () {
		String pedigreeFileS = new File (fileOfData[7], "NAM.txt").getAbsolutePath();
		String inAnchorMapDirS = new File (fileOfData[4], "impAll").getAbsolutePath();
		String outAnchorMapDirS = new File (fileOfData[4], "impNAM").getAbsolutePath();
		Pedigree ped = new Pedigree (pedigreeFileS);

	
		String[][] samples = ped.getSampleOfFamilies(ped.getFamilyStartWith("NAM_", 100), 100);
		File[] inHapMapFiles = new File (inAnchorMapDirS).listFiles();
		for (int i = 0; i < inHapMapFiles.length; i++) {
			String[] temp = inHapMapFiles[i].getName().split("\\.");
			String tem = "NAM" + temp[temp.length-3];
			String outfileS = outAnchorMapDirS + "/" + tem + ".hmp.txt";
			HapMapUtils hmu = new HapMapUtils(inHapMapFiles[i].getAbsolutePath(), outfileS);
			hmu.checkSegregationInFamilyAndOutput(samples, (float)0.2);
			System.out.println("BiParental families selected from " + inHapMapFiles[i]);
		}
	}

	public void creatTagPMap () {
		String samFileS = new File (fileOfData[6], "26Mtag-very-sensitive-local-k2.sam").getAbsolutePath();
		String tagPMapFileS = new File (fileOfData[6], "26Mtag-k2.pmap").getAbsolutePath();
		String tagCountFileS = new File (fileOfData[2], "merged_20.cnt").getAbsolutePath();
		TagCounts tc = new TagCounts(tagCountFileS, FilePacking.Byte);
		new TagPMap(tc, samFileS).writeTagPMap(tagPMapFileS);
	}

	public void creatTagMap () {
		//String gMappingFileS = new File (fileOfData[5], "merged.mapping.txt").getAbsolutePath();
		String gMappingFileS = "M:/sliceTBT-1000.mapping.txt";
		String samFileS = new File (fileOfData[6], "26Mtag-very-sensitive-local-k3.sam").getAbsolutePath();
		String tagCountFileS = new File (fileOfData[2], "merged_20.cnt").getAbsolutePath();
		//String gMappingUpdateFileS = new File (fileOfData[5], "update.mapping.txt").getAbsolutePath();
		String gMappingUpdateFileS = "M:/updateB73mapR.txt";
		TagMap_Deprecated tm = new TagMap_Deprecated (gMappingFileS, samFileS, tagCountFileS);
		tm.checkLdAndBlast();
		tm.updateMappingFile(gMappingFileS, gMappingUpdateFileS, true);

	}

	public void mergeMappingResult () {
		new PAVUtils().mergeGeneticMappingFiles(fileOfData[5].getAbsolutePath());
	}

	public void sliceTBT () {
		String sliceTBTDirS = fileOfData[3].getAbsolutePath();
		//String mergedTBT = new File (fileOfData[1], "434GFAAXX_s_2.tbt.byte").getAbsolutePath();
		String mergedTBT = "Z:/tbt_20.byte";
		int sliceNum = 200;
		PAVUtils uti = new PAVUtils();
		uti.sliceTBT(sliceTBTDirS, mergedTBT, sliceNum);
	}


	public void mapTBT () {//this step is just a test, mapping runs on clusters,
		String tbtFileS = new File (fileOfData[1], "tbt_test_50.byte").getAbsolutePath();
		//String tbtFileS = "M:/test55againstAnchor/fake1.tbt.byte";
		String anchorMapFileS = new File (fileOfData[4], "testData/impAll/Zea20120110_scv10mF8maf002_mgs_E1pLD5kpUn.imp95_1024.c+.hmp.txt").getAbsolutePath();
		//String anchorMapFileS = "N:/Zea/build20120110/imp/Zea20120110_scv10mF8maf002_mgs_E1pLD5kpUn.imp95_1024.c+.hmp.txt";
		String mappingResultFileS = "M:/mapR.txt";
		//TagCallerAgainstAnchorMT tcaa=new TagCallerAgainstAnchorMT(tbtFileS, anchorMapFileS, null, mappingResultFileS, -1, -1);
		TagAgainstAnchor taa = new TagAgainstAnchor (tbtFileS, anchorMapFileS, mappingResultFileS);
	}

	public void mergeTBTByRow () {
		MergeTagsByTaxaFilesByRowPlugin mtbt = new MergeTagsByTaxaFilesByRowPlugin ();
		mtbt.creatMergeTBTByTagCount(fileOfData[0].getAbsolutePath(), new File (fileOfData[1], "tbt_20.byte").getAbsolutePath(), new File (fileOfData[2], "merged_20.cnt").getAbsolutePath());
	}
	
	public void getTagPairs () {
		String mergedTagCountOfAllS = fileOfData[2].getAbsolutePath() + "/merged_20.cnt";
		String tagPairS = fileOfData[3].getAbsolutePath() + "/tagPair.tps";
		double etr = 0.05;
		UNetworkFilter unf = new UNetworkFilter (mergedTagCountOfAllS, etr, tagPairS);
	}

	public void creatDir (String workingDirS) {
		fileOfData = new File[childDir.length];
		for (int i = 0; i < childDir.length; i++) {
			fileOfData[i] = new File(parentDir, childDir[i]);
			fileOfData[i].mkdirs();
		}
	}

	public static void main (String[] args) {
		String workingDirS = "M:/pav/";
		new PAVgo (workingDirS);
	}
}
