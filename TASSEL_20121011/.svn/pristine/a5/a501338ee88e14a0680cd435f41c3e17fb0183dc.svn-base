/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.gbs.pipeline;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.util.HashMap;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.gbs.tagdist.TagsByTaxaBitFileMap;
import net.maizegenetics.genome.BaseEncoder;
import net.maizegenetics.genome.GBS.MutableAlignmentForGBS;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.FilterAlignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.Locus;
import net.maizegenetics.pal.ids.IdGroup;

/**
 *
 * @author glaubitz
 */
public class GeneticallyMapB73TagsVsNAMFramework {
    static String[] ZStems = {"Z001","Z002","Z003","Z004","Z005","Z006","Z007","Z008","Z009","Z010",
                              "Z011","Z012","Z013","Z014","Z015","Z016"       ,"Z018","Z019","Z020",
                              "Z021","Z022","Z023","Z024","Z025","Z026"}; // no Z017 = IBM
    static String segFamStr;


    public static void main(String[] args) {
//        convertNAMMapToHapMap();
//        removeEmptySitesFromNAMHapMap();
//        imputeMissingDataNAMHapMap();
        mapTagsVsNAMHapMap();
    }

    public static void convertNAMMapToHapMap() {
        String HMMNamMapFileName = "/usr/local/maizediv/illumina/Zea/build20111217/HMM/NAM282_20111217_scv10mF8maf002_mgs_E1pLD5_chr+.hmp.abhv2.impute5to3stateHMM.txt";
        String NamHapMapFileName = "/usr/local/maizediv/illumina/Zea/build20111217/HMM/NAM282_20111217_scv10mF8maf002_mgs_E1pLD5_HMM_chr+.hmp.txt";
        
        for (int chr = 1; chr < 11; chr++) {
            String infile=HMMNamMapFileName.replace("+", ""+chr);
            String outfile=NamHapMapFileName.replace("+", ""+chr);
            String[] inputLine;
            int lineNum=0;
            String header, temp;
            System.out.println("Reading: "+infile);
            try {
                BufferedReader br = new BufferedReader(new FileReader(new File(infile)), 65536);
                DataOutputStream fw = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfile), 65536));
                header = br.readLine();
                ++lineNum;
                fw.writeBytes(header+"\n");
                header=null;
                while ((temp = br.readLine()) != null) {
                    ++lineNum;
                    inputLine = temp.split("\t");
                    for (int col = 0; col < 11; col++) {
                        if (col == 0) {
                            fw.writeBytes(inputLine[col]);
                        } else {
                            fw.writeBytes("\t"+inputLine[col]);
                        }
                    }
                    for (int col = 11; col < inputLine.length; col++) {
                        if (inputLine[col].equals("-1")) {
                            fw.writeBytes("\tN");
                        } else if (inputLine[col].equals("0")) {
                            fw.writeBytes("\tA");
                        } else if (inputLine[col].equals("1")) {
                            fw.writeBytes("\tM");
                        } else if (inputLine[col].equals("2")) {
                            fw.writeBytes("\tC");
                        } else {
                            System.out.println("WARNING: invalid genotype found on line "+lineNum+" in column "+col+": "+inputLine[col]);
                        }
                    }
                    fw.writeBytes("\n");
                    if (lineNum%1000==0) System.out.println("    ...finished reading/writing "+lineNum+" lines");
                }
                fw.close();
            } catch (Exception e) {
                System.out.println("Catch in reading HMMNamMapFile file e=" + e);
                e.printStackTrace();
            }
        }
    }

    public static void removeEmptySitesFromNAMHapMap() {
        String NamHapMapFileName = "/usr/local/maizediv/illumina/Zea/build20111217/HMM/NAM282_20111217_scv10mF8maf002_mgs_E1pLD5_HMM_chr+.hmp.txt";
        String FilteredNamHapMapFileName = "/usr/local/maizediv/illumina/Zea/build20111217/HMM/NAM282_20111217_scv10mF8maf002_mgs_E1pLD5_HMM_filtSites_chr+.hmp.txt";
        for (int chr = 1; chr < 11; chr++) {
            String infile=NamHapMapFileName.replace("+", ""+chr);
            String outfile=FilteredNamHapMapFileName.replace("+", ""+chr);
            System.out.println("Reading: "+infile);  
            Alignment a;
            try {
                a=ImportUtils.readFromHapmap(infile);  
            } catch(Exception e) {
                System.out.println("Could not read input hapmap file for chr"+chr+":\n\t"+infile+"\n\tSkipping...");
                continue;
            }
            System.out.println("Original Alignment  Taxa:" + a.getSequenceCount() + " Sites:" + a.getSiteCount());
            if (a.getSiteCount()==0) continue;
            int[] goodLowHetSites = AlignmentFilterByGBSUtils.getLowHetSNPs(a, true, -2.0, 45, 0.01, 1.1);
            a = FilterAlignment.getInstance(a, goodLowHetSites);
            System.out.println("SiteFiltered Alignment  Taxa:" + a.getSequenceCount() + " Sites:" + a.getSiteCount());
            if (a.getSiteCount()==0) continue;
            AlignmentFilterByGBSUtils.getCoverage_MAF_F_Dist(a,false);
            ExportUtils.writeToHapmap(a, false, outfile, '\t');
            System.out.println("File written after filtering out empty sites:"+outfile);
        }
    }

    public static void imputeMissingDataNAMHapMap() {
        String FilteredNamHapMapFileName = "/usr/local/maizediv/illumina/Zea/build20111217/HMM/NAM282_20111217_scv10mF8maf002_mgs_E1pLD5_HMM_filtSites_chr+.hmp.txt";
        String imputedNamHapMapFileName = "/usr/local/maizediv/illumina/Zea/build20111217/HMM/NAM282_20111217_scv10mF8maf002_mgs_E1pLD5_HMM_filtSites_imp_chr+.hmp.txt";
        for (int chr = 1; chr < 11; chr++) {
            String infile=FilteredNamHapMapFileName.replace("+", ""+chr);
            String outfile=imputedNamHapMapFileName.replace("+", ""+chr);
            System.out.println("Reading: "+infile);
            Alignment a;
            try {
                a=ImportUtils.readFromHapmap(infile);
            } catch(Exception e) {
                System.out.println("Could not read input hapmap file for chr"+chr+":\n\t"+infile+"\n\tSkipping...");
                continue;
            }
            System.out.println("Original Alignment  Taxa:" + a.getSequenceCount() + " Sites:" + a.getSiteCount());
            if (a.getSiteCount()==0) continue;
            MutableAlignmentForGBS mutAlign=new MutableAlignmentForGBS(a);
            mutAlign.imputeMissingDataIncludingHets();
            ExportUtils.writeToHapmap(mutAlign, false, outfile, '\t');
            System.out.println("File written after imputing missing data between consecutive, equivalent markers:"+outfile);
        }
    }

    public static void mapTagsVsNAMHapMap() {
        String workDir = "/usr/local/maizediv/illumina/";
        String TOPMFileName =                workDir+"NAM_Ames_282/mergedNAM282Ames_072011.topm.bin";
        String TBTBitFileMapFileName =       workDir+"NAM_Ames_282/mergedNAM282Ames.tbt.bin";
        String usefulNAMFullTaxonNamesFile = workDir+"NAM_Ames_282/mergedNAM282Ames_tbt_TaxaList_UseVsNAMFramework.txt";
        String imputedNamHapMapFileName = workDir+"Zea/build20111217/HMM/NAM282_20111217_scv10mF8maf002_mgs_E1pLD5_HMM_filtSites_imp_chr+.hmp.txt";
        String outFile =                  workDir+"Zea/build20111217/HMM/mappedTags/mappedB73TagsFromMergedNAM282AmesTBT_ouput6_aligned_20120314.txt";
        String outLine;
        int[] segFams = null;
        int nFams = 25;
        TagsOnPhysicalMap theTOPM = new TagsOnPhysicalMap(TOPMFileName, true);
        TagsByTaxaBitFileMap theTBT = new TagsByTaxaBitFileMap(TBTBitFileMapFileName);
        String[] theTaxa = theTBT.getTaxaNames();
        boolean[] B73Mask = getB73Mask(theTaxa);
        int[] taxaMapAlignmentToTBT = new int[theTaxa.length];
        int[][] taxaMapAlignmentToTBTByFamily = new int[nFams][theTaxa.length];
        int[] familyOfEachTaxon = new int[theTaxa.length];
        int[] taxaMapAlignmentToTBTSegFams = new int[theTaxa.length];
        Alignment[] NAMMap = loadFrameWorkMap(imputedNamHapMapFileName);
        boolean[][] familyMasks = getMaskedTaxaByFamily(theTaxa, taxaMapAlignmentToTBT, NAMMap[0].getIdGroup(), taxaMapAlignmentToTBTByFamily, familyOfEachTaxon, usefulNAMFullTaxonNamesFile);
        int nB73Tags = 0, nTagsWithConsenseChr = 0;
        try {
            DataOutputStream fw=new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outFile),65536));
            outLine = "Tag\tTagIndexInTBT\tTagLen\tPhysChr\tStrand\tStartPos\tDivergence\t"+GeneticMapResult.writeHeader()+"\tGenStartPos\tGenEndPos\tnB73SoFar\tSegFams";
            System.out.println(outLine);
            fw.writeBytes(outLine); fw.writeBytes("\n");
            int tagTBTIndex;
            for (tagTBTIndex = 2890991; tagTBTIndex < theTBT.getTagCount(); ++tagTBTIndex) {
                if (tagTBTIndex % 1000 == 0) fw.flush();
                long[] tag = theTBT.getTag(tagTBTIndex);
                byte[] tagDist = theTBT.getTaxaReadCountsForTag(tagTBTIndex);
                if (getB73Presence(tagDist, B73Mask) >= 0.1) {
                    ++nB73Tags;
                    int tagTOPMIndex=theTOPM.getTagIndex(tag);
                    if (tagTOPMIndex < 0) {
                        System.out.println("Warning: tag "+BaseEncoder.getSequenceFromLong(tag)+" not found in the TOPM");
                        continue;
                    }
                    if (theTOPM.getChromosome(tagTOPMIndex)<0) continue;  // only work with aligned tags in this run
                    segFams = getSegFams(tagDist, familyMasks, NAMMap, taxaMapAlignmentToTBTByFamily);

                    boolean hasConsensusChr = false;
                    int consensChrIndex = -1;
                    for (int fam = 0; fam < segFams.length; fam++) {
                        if (segFams[fam]>0) {
                            consensChrIndex = segFams[fam] - 1;
                            nTagsWithConsenseChr++;
                            hasConsensusChr = true;
                            break;
                        }
                    }
                    if (hasConsensusChr) {
                        updateTaxaMapAlignmentToTBTSegFams(taxaMapAlignmentToTBTSegFams, taxaMapAlignmentToTBTByFamily, segFams, familyOfEachTaxon);
                        GeneticMapResult gmr=findBestHitOnChromosome(NAMMap[consensChrIndex], tagDist, taxaMapAlignmentToTBTSegFams, 30, 1);
                        if (gmr.bestP < 0.05) {
                            outLine =     BaseEncoder.getSequenceFromLong(tag)
                                    +"\t"+tagTBTIndex
                                    +"\t"+theTOPM.getTagLength(tagTOPMIndex)
                                    +"\t"+theTOPM.getChromosome(tagTOPMIndex)
                                    +"\t"+theTOPM.getStrand(tagTOPMIndex)
                                    +"\t"+theTOPM.getStartPosition(tagTOPMIndex)
                                    +"\t"+theTOPM.getDivergence(tagTOPMIndex)
                                    +"\t"+gmr.toString()
                                    +"\t"+NAMMap[consensChrIndex].getPositionInLocus(gmr.startPos)
                                    +"\t"+NAMMap[consensChrIndex].getPositionInLocus(gmr.endPos)
                                    +"\t"+nB73Tags
                                    +"\t"+segFamStr;
                            System.out.println(outLine);
                            fw.writeBytes(outLine); fw.writeBytes("\n");
                        }
                    }
                }
            }
            System.out.println("We have reached the end of the TBT: "+tagTBTIndex);
            fw.close();
        } catch (Exception e) {
            System.out.println("Error writing to file:"+e.getMessage());
            e.printStackTrace();
        }
    }

    private static boolean[] getB73Mask(String[] theTaxa) {
        boolean[] B73Mask = new boolean[theTaxa.length];
        System.out.println("\nB73 taxa in TBT:");
        for (int taxon = 0; taxon < theTaxa.length; taxon++) {
            if (theTaxa[taxon].startsWith("B73:")) {
                B73Mask[taxon] = true;
                System.out.println(theTaxa[taxon]);
            }
        }
        return B73Mask;
    }

    private static double getB73Presence(byte[] tagDist, boolean[] B73Mask) {
        int nB73 = 0, nB73WithTag = 0;
        for (int taxon = 0; taxon < tagDist.length; taxon++) {
            if (B73Mask[taxon]) {
                nB73++;
                if (tagDist[taxon]>0) {
                    nB73WithTag++;
                }
            }
        }
        return (double) nB73WithTag / (double) nB73;
    }

    private static Alignment[] loadFrameWorkMap(String fileNameTemplate) {
        Alignment[] frameWorkMap = new Alignment[10];
        for (int chr = 1; chr < 11; chr++) {
            String infile=fileNameTemplate.replace("+", ""+chr);
            System.out.println("Reading: "+infile);
            try {
                frameWorkMap[chr-1]=ImportUtils.readFromHapmap(infile);
            } catch(Exception e) {
                System.out.println("Could not read input hapmap file for chr"+chr+":\n\t"+infile+"\n\tSkipping...");
                continue;
            }
            System.out.println("Alignment  Taxa:" + frameWorkMap[chr-1].getSequenceCount() + " Sites:" + frameWorkMap[chr-1].getSiteCount());
            AlignmentFilterByGBSUtils.getCoverage_MAF_F_Dist(frameWorkMap[chr-1],false);
        }
        return frameWorkMap;
    }

    private static boolean[][] getMaskedTaxaByFamily(
            String[] taxaInTBT, int[] taxaMapAlignmentToTBT, IdGroup frameWorkMapChr1IdGroup, int[][] taxaMapAlignmentToTBTByFamily, int[] familyOfEachTaxon, String usefulNAMFullTaxonNamesFile) {

        boolean[][] familyMasks = new boolean[ZStems.length][taxaInTBT.length];  // which taxa in the TBT belong in each family
        HashMap<String,Integer> ZStemIndices = new HashMap<String,Integer>();
        for (int familyIndex = 0; familyIndex < ZStems.length; familyIndex++) {
            ZStemIndices.put(ZStems[familyIndex], familyIndex);
        }
        String[] usefulTaxa = new String[taxaInTBT.length];
        int lineNum = 0;
        String temp;
        String[] inputLine;
        HashMap<String,String> usefulTaxaKey = new HashMap<String,String>();
        System.out.println("Reading: "+usefulNAMFullTaxonNamesFile);
        try {
            BufferedReader br = new BufferedReader(new FileReader(new File(usefulNAMFullTaxonNamesFile)), 65536);
            br.readLine();  // skip the header line
            ++lineNum;
            while ((temp = br.readLine()) != null) {
                ++lineNum;
                inputLine = temp.split("\t");
                usefulTaxaKey.put(inputLine[0], inputLine[1]);
            }
            System.out.println("Finished reading "+lineNum+" lines\n\n");
        } catch (Exception e) {
            System.out.println("Catch in reading HMMNamMapFile file e=" + e);
            e.printStackTrace();
        }
        System.out.println("\n\ntaxonInTBT\tusefulTaxon");
        for (int t = 0; t < taxaInTBT.length; ++t) {
            familyOfEachTaxon[t] = -1;
            if (usefulTaxaKey.containsKey(taxaInTBT[t])) {  // all the NAM RILs are in the usefulTaxaKey list
                usefulTaxa[t] = usefulTaxaKey.get(taxaInTBT[t]);
                String ZStem = usefulTaxa[t].substring(0, 4);
                if (ZStemIndices.containsKey(ZStem)) {  // some of the NAM RILs are duplicates and should be ignored (where usefulTaxa[t] = "ignore")
                    int z = ZStemIndices.get(ZStem);
                    familyMasks[z][t] = true;
                    familyOfEachTaxon[t] = z;
                }
            } else {
                usefulTaxa[t] = "ignore";  // non NAM RILs
            }
            System.out.println(taxaInTBT[t]+"\t"+usefulTaxa[t]);
            taxaMapAlignmentToTBT[t] = frameWorkMapChr1IdGroup.whichIdNumber(usefulTaxa[t]);  // taxa to ignore are named "ignore" and will result in a value of -1
        }
        for (int family = 0; family < familyMasks.length; family++) {
            for (int t = 0; t < familyMasks[family].length; t++) {
                if (familyMasks[family][t]) {
                    taxaMapAlignmentToTBTByFamily[family][t] = taxaMapAlignmentToTBT[t];
                } else {
                    taxaMapAlignmentToTBTByFamily[family][t] = -1;
                }
            }
        }
        return familyMasks;
    }

    private static void updateTaxaMapAlignmentToTBTSegFams(int[] taxaMapAlignmentToTBTSegFams, int[][] taxaMapAlignmentToTBTByFamily, int[] segFams, int[] familyOfEachTaxon) {
        for (int taxon = 0; taxon < taxaMapAlignmentToTBTSegFams.length; taxon++) {
            if (familyOfEachTaxon[taxon] > -1 && segFams[familyOfEachTaxon[taxon]]>0) {  // a value greater than 0 is equal to the consensus chr
                taxaMapAlignmentToTBTSegFams[taxon] = taxaMapAlignmentToTBTByFamily[familyOfEachTaxon[taxon]][taxon];
            } else {
                taxaMapAlignmentToTBTSegFams[taxon] = -1;
            }
        }
    }

    private static int[] getSegFams(byte[] tagDist, boolean[][] familyMasks, Alignment[] frameWorkMap, int[][] taxaMapAlignmentToTBTByFamily) {
        int minCountForTesting = 20, initSkimRate = 20;
        double recombThresh = 0.05;
        int[] chrCounts = new int[10];  // index = (chr# - 1);  value = # of families for which that chr was the best
        int[] segFams = new int[familyMasks.length]; // non-segregating fams will have a value of 0;  segregating families will have a value of bestChr
        for (int family = 0; family < familyMasks.length; family++) {
            int bestSkimAlign = -1;
            double bestSkimP = 1;
            GeneticMapResult gmr=null;
            for (int currAlign = 0; currAlign < frameWorkMap.length; currAlign++) {
                gmr=findBestHitOnChromosome(frameWorkMap[currAlign], tagDist, taxaMapAlignmentToTBTByFamily[family], minCountForTesting, initSkimRate);
                if(gmr.bestP<bestSkimP) {bestSkimP=gmr.bestP; bestSkimAlign=currAlign;}
            }  //end of chromosomes
            if (bestSkimP < recombThresh) {
                chrCounts[bestSkimAlign]++;
                segFams[family] = bestSkimAlign + 1;  // value of 0 means that the tag did not map to any chrs for that family
            }
        }
        int consensChr = 0, maxChrBestCount = 0;
        for (int ci = 0; ci < chrCounts.length; ci++) {  // ci means "chromosome index"
            if (chrCounts[ci] > maxChrBestCount) {
                maxChrBestCount = chrCounts[ci];
                consensChr = ci+1;
            }
        }
        segFamStr = "";
        for (int family = 0; family < segFams.length; family++) {
            if (segFams[family] != consensChr) segFams[family] = 0;  // disregard families for which the tag maps to a non-consensus chr
            else segFamStr += "_"+ZStems[family];
        }
        segFamStr.replaceFirst("_", "");
        return segFams;
    }

    private static GeneticMapResult findBestHitOnChromosome(Alignment gMap, byte[] tagDist, int[] taxaMapAlignmentToTBT, int minCountForTesting, int stepSize) {
        int startBest = -1, endBest = -1, altCntBest = -1, refCntBest = -1;
        double bestP = 1;
        for (int geneticSite = 0; geneticSite < gMap.getSiteCount(); geneticSite+=stepSize) {
            int refCnt = 0, altCnt = 0;
            for (int taxon = 0; taxon < taxaMapAlignmentToTBT.length; taxon++) {
                if (taxaMapAlignmentToTBT[taxon] < 0) {
                    continue;
                }
                if (tagDist == null || tagDist[taxon] < 1) {
                    continue;
                }
                byte currBase = gMap.getBase(taxaMapAlignmentToTBT[taxon], geneticSite);
                if (currBase == AlignmentFilterByGBSUtils.refAllele) {
                    refCnt++;
                } else if (currBase == AlignmentFilterByGBSUtils.altAllele) {
                    altCnt++;
                }
            }
            int cnt = refCnt + altCnt;
            if (cnt < minCountForTesting || (double) refCnt/altCnt > 5.0) continue;
            double currP = (double) altCnt / (double) cnt;  // recombination ratio assuming a B73 tag (one-tailed test)
            if (currP < bestP) {
                bestP = currP;
                startBest = endBest = geneticSite;
                altCntBest = altCnt;
                refCntBest = refCnt;
            } else if (currP == bestP) {
                endBest = geneticSite;
            }
        }  //end of one chromosome
        GeneticMapResult gmr=new GeneticMapResult(startBest, endBest, refCntBest, altCntBest, bestP, gMap.getLocus(startBest));
        return gmr;
    }

    private static class GeneticMapResult {
        int startPos, endPos, altCnt, refCnt;
        double bestP;
        Locus theChromosome;

        GeneticMapResult(int startPos, int endPos, int refCnt, int altCnt, double bestP, Locus theChromosome) {
            this.startPos=startPos;
            this.endPos=endPos;
            this.refCnt=refCnt;
            this.altCnt=altCnt;
            this.bestP=bestP;
            this.theChromosome=theChromosome;
        }

        @Override
        public String toString() {
            String str = theChromosome.getChromosomeName()+"\t"+bestP+"\t"+refCnt+"\t"+altCnt;
            return str;
        }

        public static String writeHeader() {
            String header = "GenChr\tRecomb\tRefCount\tAltCount";
            return header;
        }
    }
}
