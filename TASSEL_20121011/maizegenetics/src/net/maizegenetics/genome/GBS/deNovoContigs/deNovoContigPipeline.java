/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS.deNovoContigs;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.util.TreeMap;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.genome.GBS.*;
import net.maizegenetics.pal.alignment.CombineAlignment;
import net.maizegenetics.pal.alignment.ImportUtils;

/**
 *
 * @author jcg233
 */
public class deNovoContigPipeline {
    public deNovoContigPipeline() {
    }

    public static void main(String[] args) {
//        physicallyMapAllTags();
        runPipeline();
//        runPipelineMo17();

//        identifyNovelContigs();
//        printSampleOfTags();
//        lookupGenMapResults();
        filterForPhyOrGenMappedTags();
//        filterForGeneticallyMappedTags();
//        filterForGenConsensus();
//        determinePhyPositCtgs();
//        mergeIntoHaplotypes();
//        geneticallyMapContigHaplotypes();
    }

    public static void physicallyMapAllTags() {
        String AGPv2VDfile = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/VirtualDigests/B73AGPv2VD.rwpm.bin";

        String assembly = "454k96II";
        String workingPath = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/de_novo_ctgs/" + assembly + "/";

        // rwpm files for the various sets of contigs, using a fake Chr20 as the physical position
        String B73K50VD454file = workingPath + "454k50VD.rwpm.bin";
        String B73K60VD454file = workingPath + "454k60VD.rwpm.bin";
        String B73K80VD454file = workingPath + "454k80VD.rwpm.bin";
        String B73K96VD454file = workingPath + "454k96VD.rwpm.bin";
        String B73K96IIVD454file = workingPath + "454k96IIVD.rwpm.bin";
        String flcDNAvdFile    = workingPath + "flcdna.rwpm.bin";
        String Mo17BGIvdFile    = workingPath + "Mo17_CAU_VD.rwpm.bin";

        // Contig key files for the various sets of contigs, providing the extent of each contig (multiple of 64) on fake Chr20
        String B73K50VD454CtgKeyFile = workingPath + "454asm_k50_tags_ctg_key.txt";
        String B73K60VD454CtgKeyFile = workingPath + "454asm_k60_tags_ctg_key.txt";
        String B73K80VD454CtgKeyFile = workingPath + "454asm_k80_tags_ctg_key.txt";
        String B73K96VD454CtgKeyFile = workingPath + "454asm_k96_tags_ctg_key.txt";
        String B73K96IIVD454CtgKeyFile = workingPath + "k96II_tags_ctg_key.txt";
        String flcDNAvdCtgKeyFile    = workingPath + "flcdna_tags_ctg_key.txt";
        String Mo17BGIvdCtgKeyFile    = workingPath + "Mo17_virtualDigest_ctg_key.txt";

        String allTagsOutFile        = workingPath + assembly + ".dnct.bin";

        int maxDiv = 20;

        System.out.println("\nPHYSICALLY MAPPING ALL DE NOVO CONTIG TAGS");
        System.out.println("------------------------------------------");
        ReadsWPhysicalMap AGPv2VD = new ReadsWPhysicalMap(AGPv2VDfile,true);
        ReadBLASTer theReadBLASTer = new ReadBLASTer(AGPv2VD);
        ReadsWPhysicalMap ContigsChr20 = new ReadsWPhysicalMap(B73K96IIVD454file,true);
        ContigsChr20.sortTable(false);  // sort by position on fake Chr20
        int nPerfectSingleMatches = 0, nPerfectMultiMatches = 0, nSingleWMismatch = 0, nMultiWMismatch = 0, noMap = 0, nWritten = 0, nTotal = 0;

        System.out.println();
        System.out.println("Finding de novo tags (maxDivergence = " + maxDiv + ")");
        System.out.println("----------------------------------------");

        try {
            DataOutputStream fw=new DataOutputStream(new BufferedOutputStream(new FileOutputStream(allTagsOutFile),4000000));
            BufferedReader br = new BufferedReader(new FileReader(B73K96IIVD454CtgKeyFile), 65536);
            String inputLine = br.readLine();  // skip the header line
            inputLine = br.readLine();
            String[] fields = inputLine.split("\t");
            int Contig = Integer.parseInt(fields[0]);
            int chr20_PosMin = Integer.parseInt(fields[1]);
            int chr20_PosMax = Integer.parseInt(fields[2]);
            for (int i = 0; i < ContigsChr20.getSize(); ++i) {
                ReadWithPositionInfo deNovoTagRWPI = ContigsChr20.getReadWithPosition(i);
                ++nTotal;
                while (deNovoTagRWPI.getPositionMax() > chr20_PosMax) {
                    inputLine = br.readLine();
                    fields = inputLine.split("\t");
                    Contig = Integer.parseInt(fields[0]);
                    chr20_PosMin = Integer.parseInt(fields[1]);
                    chr20_PosMax = Integer.parseInt(fields[2]);
                }
                if (i % 100000 == 0) {
                    System.out.println("tagIndex: " + i + ";  PositionMax: " + deNovoTagRWPI.getPositionMax() + ";  chr20_PosMax: " + chr20_PosMax);
                }
                long[] currTag = deNovoTagRWPI.getRead();
                TagFromDeNovoContig deNovoTag = new TagFromDeNovoContig();
                deNovoTag.setTagSeq(currTag);
                deNovoTag.setNovelContig(Contig);
                deNovoTag.setTagNum( (deNovoTagRWPI.getPositionMax()-chr20_PosMin+1)/64 );
                deNovoTag.setNTagsInCtg( (chr20_PosMax-chr20_PosMin+1)/64 );
                int[] phyHits = AGPv2VD.getReadIndexSet(currTag);
                if (phyHits != null) {
                    if (phyHits.length == 1) {
                        ++nPerfectSingleMatches;
                        deNovoTag.setPhyChr(AGPv2VD.getReadWithPosition(phyHits[0]).getChromosome());
                        deNovoTag.setPositionMin(AGPv2VD.getReadWithPosition(phyHits[0]).getPositionMin());
                        deNovoTag.setPositionMax(AGPv2VD.getReadWithPosition(phyHits[0]).getPositionMax());
                        deNovoTag.setStrand(AGPv2VD.getReadWithPosition(phyHits[0]).getStrand());
                        deNovoTag.setDivergence(0);
                        deNovoTag.setNextCutDistance(AGPv2VD.getReadWithPosition(phyHits[0]).getNextCutDistance());
                        deNovoTag.setMultimaps(1);
                    } else {
                        deNovoTag.setDivergence(0);
                        deNovoTag.setMultimaps(phyHits.length);
                        ++nPerfectMultiMatches;
                    }
                } else {
                    TreeMap<Integer, Integer> align;
                    align = theReadBLASTer.findMatchesWithIntLengthWords(currTag, maxDiv, true);
                    if (align.size() == 1) {
                        ReadWithPositionInfo AGPv2Tag = AGPv2VD.getReadWithPosition(align.firstEntry().getKey());
                        deNovoTag.setPhyChr(AGPv2Tag.getChromosome());
                        deNovoTag.setPositionMin(AGPv2Tag.getPositionMin());
                        deNovoTag.setPositionMax(AGPv2Tag.getPositionMax());
                        deNovoTag.setStrand(AGPv2Tag.getStrand());
                        deNovoTag.setDivergence(align.firstEntry().getValue());
                        deNovoTag.setNextCutDistance(AGPv2Tag.getNextCutDistance());
                        deNovoTag.setMultimaps(1);
                        ++nSingleWMismatch;
                    } else if (align.size() > 1) {
                        deNovoTag.setDivergence(align.firstEntry().getValue());
                        deNovoTag.setMultimaps(align.size());
                        ++nMultiWMismatch;
                    } else {
                        ++noMap;
                    }
                }
                fw.write(deNovoTag.toByte());
                ++nWritten;
            }
            fw.flush();
            fw.close();
            System.out.println("Total_nTags  " + nTotal);
            System.out.println("nTags_written  " + nWritten);
            System.out.println("nPerfectSingleMatches  " + nPerfectSingleMatches);
            System.out.println("nPerfectMultiMatches  " + nPerfectMultiMatches);
            System.out.println("nSingleWMismatch  " + nSingleWMismatch);
            System.out.println("nMultiWMismatch  " + nMultiWMismatch);
            System.out.println("nNotMapping  " + noMap);
            System.out.println("maxDivergence  " + maxDiv);
        } catch (Exception e) {
            System.out.println("Read/write error in physicallyMapAllTags: " + e);
        }
    }

    public static void printSampleOfTags() {
        String allTagsFile = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/ShiranVDs/454k96.nct.bin";
        String novelContigTagsFile = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/ShiranVDs/454k80novel.nct.bin";

        deNovoContigTags nct454 = new deNovoContigTags(allTagsFile, true);
        nct454.printRows(100000);
    }

    public static void determinePhyPositCtgs() {
        String genConsensusTagsFile = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/ShiranVDs/454k96genConsensus.nct.bin";
        String contigPhyPositFile = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/ShiranVDs/454k96genConsensusWPhyPosit.nct.txt";

        deNovoContigTags nct454genConsensus = new deNovoContigTags(genConsensusTagsFile, true);
        nct454genConsensus.determineConsensusPhyPositions(0);  // want maxDiv of 0 (perfect matches) for sequence from B73
        nct454genConsensus.writeToFile(new File(contigPhyPositFile), false);
    }

    public static void identifyNovelContigs() {
        String allTagsFile = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/ShiranVDs/454flcDNA.nct.bin";
        String novelContigTagsFile = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/ShiranVDs/454k96novel.nct.bin";

        deNovoContigTags nct454all = new deNovoContigTags(novelContigTagsFile, true);
        nct454all.findNovelContigs(0,0.25,0.8);  // want maxDiv of 0 (perfect matches) for sequence from B73
//        nct454all.writeNovelContigsOnlyToFile(new File(novelContigTagsFile), true);
//        nct454all.writeToFile(new File(allTagsFile), true);  // overwrite the original file with the NovelContig info (inNovelCtg = 1 or 0)
    }

    public static void lookupGenMapResults() {
        String allTagsFile = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/ShiranVDs/454flcDNA.nct.bin";
        String GenMapResultFile = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/ShiranVDs/IBMgeneticResults_20100828_JG.txt"; // 485860 tags

        deNovoContigTags nct454 = new deNovoContigTags(allTagsFile, true);
        nct454.lookupGeneticMappingResult(new File(GenMapResultFile), 485860, true);  // true indicates that we are working with B73 (reference) sequence
        nct454.writeToFile(new File(allTagsFile), true);  // overwrite the original file with the Genetic Mapping Results
    }

    public static void filterForPhyOrGenMappedTags() {
        String allTagsFile    = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/de_novo_ctgs/454k96II/454k96II.dnct.bin";
        String mappedTagsFile = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/de_novo_ctgs/454k96II/454k96II.POrGMapped.dnct.txt";

        deNovoContigTags nct454all = new deNovoContigTags(allTagsFile, true);
        nct454all.writeTagsWPhyOrGenPositionToFile(new File(mappedTagsFile), false);  // true = binary
    }

    public static void filterForGeneticallyMappedTags() {
        String allTagsFile = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/ShiranVDs/454flcDNA.nct.bin";
        String genMappedTagsFile = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/ShiranVDs/454flcDNAgenMapped.nct.bin";

        deNovoContigTags nct454all = new deNovoContigTags(allTagsFile, true);
        nct454all.writeTagsWGeneticPositionToFile(new File(genMappedTagsFile), true, true);  // true = reference (B73), true = binary
    }

    public static void filterForGenConsensus() {
        String genMappedTagsFile = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/ShiranVDs/454k96genMapped.nct.bin";
        String genConsensusTagsFile = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/ShiranVDs/454k96genConsensus.nct.txt";

        deNovoContigTags nct454genMapped = new deNovoContigTags(genMappedTagsFile, true);
        nct454genMapped.determineConsensusGenPositions(true);  // true = sequence is from reference (B73)
        nct454genMapped.writeTagsWGeneticPositionToFile(new File(genConsensusTagsFile), true, false);  // true = reference (B73), false = not binary
    }

    public static void mergeIntoHaplotypes() {
        String genConsensusTagsFile = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/ShiranVDs/454k96genConsensus.nct.bin";
        String rbtFileIBM =           "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/ShiranVDs/rbtIBM_Min10_20100827.bin";
        String contigsByTaxaFile =    "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/ShiranVDs/454k96ContigsByTaxa.txt";

        deNovoContigTags nct454genConsensus = new deNovoContigTags(genConsensusTagsFile, true);
        ReadsByTaxa rbtIBM = new ReadsByTaxa(rbtFileIBM, true);
        nct454genConsensus.mergeIntoHaplotypes(rbtIBM, contigsByTaxaFile);
    }

    public static void geneticallyMapContigHaplotypes() {
        String frameWorkGenMapFile55K = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/de_novo_ctgs/IBM_55K_08192010.hmp.txt";
        String frameWorkGenMapFileGBS = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/434LFAAXX/anchor_map/IBM_GBS_framework_map_20100917.hmp.txt";
        String contigsByTaxaFile =      "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/de_novo_ctgs/454k96/454k96contigsByTaxa.txt";
        String contig454ResultsFile =   "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/de_novo_ctgs/454k96/454k96ctgHaplosGenResultsGBS.txt";

        // map the contig haplotypes (merged data across multiple, consensus tags) vs either the 55K or GBS framework map
        Alignment[] genFrameMap= ((CombineAlignment)ImportUtils.readFromHapmap(frameWorkGenMapFileGBS)).getAlignments();
        ContigHaplotypesByTaxa CtgsByTaxa = new ContigHaplotypesByTaxa(contigsByTaxaFile,false);
        mapContigHaplotypes.geneticallyMapContigHaplos(genFrameMap,CtgsByTaxa,0.1,20,10,contig454ResultsFile);
    }

    public static void runPipeline() {
        String assembly = "454k96II";
        String masterPath =  "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/de_novo_ctgs/";
        String workingPath = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/de_novo_ctgs/" + assembly + "/";

        String GenMapResultFile =       masterPath + "IBMgeneticResults3FCs_20110516_JG.txt"; // 589002 tags
//        String GenMapResultFile =       masterPath + "IBMgeneticResults_20100828_JG.txt"; // 485860 tags
        String rbtFileIBM =             masterPath + "rbt3IBM_Min10_20110315.bin";
//        String rbtFileIBM =             masterPath + "rbtIBM_Min10_20100827.bin";
        String frameWorkGenMapFile55K = masterPath + "IBM_55K_08192010.hmp.txt";
        String frameWorkGenMapFileGBS = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/434LFAAXX/anchor_map/IBM_GBS_framework_map_20100917.hmp.txt";

        String deNovoContigTagsFile =  workingPath + assembly + ".dnct.bin";
        String genMappedTagsFile =     workingPath + assembly + "genMapped.dnct.bin";
        String genConsenTagsFile =     workingPath + assembly + "genConsen.dnct.bin";
        String genConsenTagsTextFile = workingPath + assembly + "genConsen.dnct.txt";
        String contigsByTaxaFile =     workingPath + assembly + "contigsByTaxa.txt";
        String contig454ResultsFile =  workingPath + assembly + "ctgHaplosGenResults55K.txt"; //  <--- ***CHANGE "55k" to "GBS" if using GBS framework map

        deNovoContigTags dnct454all = new deNovoContigTags(deNovoContigTagsFile, true);
        dnct454all.determineConsensusPhyPositions(0);  // want maxDiv of 0 (perfect matches) for sequence from B73
        dnct454all.findNovelContigs(0,0.25,0.8);  // (maxDiv,matchProportion,consenseThresh). Determine which contigs are novel.  Want maxDiv of 0 (perfect matches) for sequence from B73

        // look up each tag in the previous IBM genetic mapping results file for 589002 mapped GBS tags
        dnct454all.lookupGeneticMappingResult(new File(GenMapResultFile), 589002, true);  // true indicates that we are working with B73 (reference) sequence

        dnct454all.writeToFile(new File(deNovoContigTagsFile), true);  // overwrite the original file with the consensus phy posits, novel ctg info, and genetic map results

        // filter for the tags that were genetically mapped
        dnct454all.writeTagsWGeneticPositionToFile(new File(genMappedTagsFile), true, true);  // true = reference (B73), true = binary
        dnct454all = null;
        System.gc();

        // indicate which tags match the consensus genetic position of their contig
        deNovoContigTags dnct454genMapped = new deNovoContigTags(genMappedTagsFile, true);
        dnct454genMapped.determineConsensusGenPositions(true);  // true = sequence is from reference (B73)

        // write the binary file, filtered for tags matching the genetic consensus
        dnct454genMapped.writeTagsWGeneticPositionToFile(new File(genConsenTagsFile), true, true);  // true = reference (B73), true = binary

        // write a text file, too (for opening in Excel)
        dnct454genMapped.writeTagsWGeneticPositionToFile(new File(genConsenTagsTextFile), true, false);  // true = reference (B73), false = not binary
        dnct454genMapped = null;
        System.gc();

        // merge the countsByTaxa for each tag into contig haplotypes
        deNovoContigTags dnct454genConsen = new deNovoContigTags(genConsenTagsFile, true);
        ReadsByTaxa rbtIBM = new ReadsByTaxa(rbtFileIBM, true);
        dnct454genConsen.mergeIntoHaplotypes(rbtIBM, contigsByTaxaFile);
        dnct454genConsen = null;
        System.gc();

        // map the contig haplotypes (merged data across multiple, consensus tags) vs either the 55K or GBS framework map
        Alignment[] genFrameMap= ((CombineAlignment)ImportUtils.readFromHapmap(frameWorkGenMapFile55K)).getAlignments();
        ContigHaplotypesByTaxa ctgsByTaxa = new ContigHaplotypesByTaxa(contigsByTaxaFile,false);
        mapContigHaplotypes.geneticallyMapContigHaplos(genFrameMap,ctgsByTaxa,0.1,20,10,contig454ResultsFile);
    }

    public static void runPipelineMo17() {
        String assembly = "Mo17CAU";
        String masterPath =  "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/de_novo_ctgs/";
        String workingPath = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/de_novo_ctgs/" + assembly + "/";

        String GenMapResultFile =       masterPath + "IBMgeneticResults_20100828_JG.txt"; // 485860 tags
        String rbtFileIBM =             masterPath + "rbtIBM_Min10_20100827.bin";
        String frameWorkGenMapFile55K = masterPath + "IBM_55K_08192010.hmp.txt";
        String frameWorkGenMapFileGBS = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/434LFAAXX/anchor_map/IBM_GBS_framework_map_20100917.hmp.txt";

        String deNovoContigTagsFile =     workingPath + assembly + ".dnct.bin";
        String deNovoContigTagsTextFile = workingPath + assembly + ".dnct.txt";
        String genMappedTagsFile =        workingPath + assembly + "genMapped.dnct.bin";
        String genConsenTagsFile =        workingPath + assembly + "genConsen.dnct.bin";
        String genConsenTagsTextFile =    workingPath + assembly + "genConsen.dnct.txt";
        String contigsByTaxaFile =        workingPath + assembly + "contigsByTaxa.txt";
        String contig454ResultsFile =     workingPath + assembly + "ctgHaplosGenResults55K.txt"; //  <--- ***CHANGE "55k" to "GBS" if using GBS framework map

        deNovoContigTags dnct454all = new deNovoContigTags(deNovoContigTagsFile, true);
        dnct454all.determineConsensusPhyPositions(20);  // want maxDiv of 20 for sequence from Mo17
        dnct454all.findNovelContigs(20, 0.1, 0.6);  // (maxDiv,matchProportion,consenseThresh). Want maxDiv of 20 for sequence from Mo17

        // look up each tag in the previous IBM genetic mapping results file for 485860 mapped GBS tags
        dnct454all.lookupGeneticMappingResult(new File(GenMapResultFile), 485860, false);  // true indicates that we are working with B73 (reference) sequence

        dnct454all.writeToFile(new File(deNovoContigTagsFile), true);  // overwrite the original file with the consensus phy posits, novel ctg info, and genetic map results
        dnct454all.writeToFile(new File(deNovoContigTagsTextFile), false);  // write a text file too (probably too big for Excel)

        // filter for the tags that were genetically mapped
        dnct454all.writeTagsWGeneticPositionToFile(new File(genMappedTagsFile), false, true);  // true = reference (B73), true = binary
        dnct454all = null;
        System.gc();

        // indicate which tags match the consensus genetic position of their contig
        deNovoContigTags dnct454genMapped = new deNovoContigTags(genMappedTagsFile, true);
        dnct454genMapped.determineConsensusGenPositions(false);  // true = sequence is from reference (B73)

        // write the binary file, filtered for tags matching the genetic consensus
        dnct454genMapped.writeTagsWGeneticPositionToFile(new File(genConsenTagsFile), false, true);  // true = reference (B73), true = binary

        // write a text file, too (for opening in Excel)
        dnct454genMapped.writeTagsWGeneticPositionToFile(new File(genConsenTagsTextFile), false, false);  // true = reference (B73), false = not binary
        dnct454genMapped = null;
        System.gc();

        // merge the countsByTaxa for each tag into contig haplotypes
        deNovoContigTags dnct454genConsen = new deNovoContigTags(genConsenTagsFile, true);
        ReadsByTaxa rbtIBM = new ReadsByTaxa(rbtFileIBM, true);
        dnct454genConsen.mergeIntoHaplotypes(rbtIBM, contigsByTaxaFile);
        dnct454genConsen = null;
        System.gc();

        // map the contig haplotypes (merged data across multiple, consensus tags) vs either the 55K or GBS framework map
        Alignment[] genFrameMap= ((CombineAlignment)ImportUtils.readFromHapmap(frameWorkGenMapFile55K)).getAlignments();
        ContigHaplotypesByTaxa ctgsByTaxa = new ContigHaplotypesByTaxa(contigsByTaxaFile,false);
        mapContigHaplotypes.geneticallyMapContigHaplos(genFrameMap,ctgsByTaxa,0.1,20,10,contig454ResultsFile);
    }
}
