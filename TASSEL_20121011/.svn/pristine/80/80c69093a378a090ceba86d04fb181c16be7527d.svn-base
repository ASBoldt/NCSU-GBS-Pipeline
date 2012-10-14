/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pipeline;

import java.io.File;
import java.util.Arrays;
import net.maizegenetics.gbs.maps.SAMConverterPlugin;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.tagdist.TagsByTaxaBit;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByte;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByteHDF5TagGroups;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByteHDF5TaxaGroups;
import net.maizegenetics.gbs.tagdist.TagsByTaxaUtils;
import net.maizegenetics.gbs.util.CompareGenosBetweenHapMapFilesPlugin;
import net.maizegenetics.genome.GBS.MutableAlignmentForGBS;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pipeline.TasselPipeline;


/**
 *
 * @author jcg233
 */
public class JeffPipelines {
    public static void main(String[] args) {
//        sayHi();
//        convertTextTagCountsToBinary();
//        convertBinaryTagCountsToText();
//        convertBinaryTBTToText();
//        tagCountsToFastQ();
//        runQseqToTagCountPlugin();
//        runFastqToTagCountPlugin();
//        runMergeMultipleTagCountPlugin();
//        runQseqToTBTPlugin();
//        runFastqToTBTPlugin();
//        runMergeTagsByTaxaFilesPlugin();
//        printSumCountsInTBTByTaxa();
//        mergeTaxaInTBT();
//        convertSAMToTOPM();
//        filterTOPMWithPositions();
//        runTagsToSNPByAlignmentMTPlugin();
//        runMergeDuplicateSNPsPlugin();
//        runGBSHapMapFiltersPlugin();
//        runBiParentalErrorCorrection();
//        runMergeIdenticalTaxaPlugin();
//        runBinaryToTextPlugin();
//        filterMergeAllZea();
//        imputeAllZea();
//        runQseqToHapMapPlugin();
//        runRawReadsToHapMapPlugin();
//        runRawReadsToHapMapQuantPlugin();
//        analyzeD09FYACXX_5_PstIMaize();
//        analyzeB08AAABXX_1_PstI_IBM();
//        analyzeB08AAABXX_1_PstI_IBM_TBTByte();
//        analyzeIBM94ApeKI();
//        analyzeD0E3PACXX_5_PstI_PalmarChico();
//        testTOPMExpandMaxVariants();
//        runTagsToSNPByAlignmentPlugin();
//        filterHapMapForTaxa();
//        analyzeNIL28_Z026();
//        analyzeTeoFineMapping2012();
//        analyzeZakChr5FineMapping();
//        analyzeETb1FineMapping();
//        testTBTHDF5TagGroups();
//        testTBTHDF5TaxaGroups();
//        runModifyTBTHDF5Plugin();
//        runSAMConverterPlugin();
//        analyzeMGP1_low_vol();
        runCompareGenosBetweenHapMapFilesPlugin();
//        runFastImputationBitFixedWindowPlugin();
//        printChrsFromTOPM();
//        getHapMapReport();
    }

    public static void sayHi() {
        System.out.println("\n\nGreetings Panzean!!\n\n");
    }

    public static void convertTextTagCountsToBinary() {
        String textTagCountsFileS =   "C:/Users/jcg233/Documents/Bioinformatics/NextGen/HapMapV2/test_RandomPairedEndToTBT/FakeTagCounts.txt";
        String binaryTagCountsFileS = "C:/Users/jcg233/Documents/Bioinformatics/NextGen/HapMapV2/test_RandomPairedEndToTBT/FakeTagCounts.bin";
        TagCounts tc = new TagCounts(textTagCountsFileS, FilePacking.Text);
        tc.sort();
        tc.writeTagCountFile(binaryTagCountsFileS, FilePacking.Bit, 1);
    }

    public static void convertBinaryTagCountsToText() {
        String binaryTagCountsFileS = "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/mergedTagCounts/NAM49_50_ApeKI_mergedTags_min10.cnt";
        String textTagCountsFileS =   "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/mergedTagCounts/NAM49_50_ApeKI_mergedTags_min10.cnt.txt";
        TagCounts tc = new TagCounts(binaryTagCountsFileS, FilePacking.Bit);
        tc.writeTagCountFile(textTagCountsFileS, FilePacking.Text, 1);
    }

    public static void convertBinaryTBTToText() {
        String workDir = "H:/NAM_ApeKI_plates49_50/mergedTBT/";
        String binaryTBTFileS = workDir + "NAM49_50_ApeKI_mergedTBT_min10_20110720.mergedTaxa.tbt.bin";
        String textTBTFileS =   workDir + "NAM49_50_ApeKI_mergedTBT_min10_20110720.mergedTaxa.tbt.txt";
        TagsByTaxa tbt = new TagsByTaxaBit(binaryTBTFileS, FilePacking.Bit);
        File textTBTFile = new File(textTBTFileS);
        tbt.writeDistFile(textTBTFile, FilePacking.Text, 0);
    }
    
    public static void tagCountsToFastQ() {
        String TagCountFileName = "H:/64GJAAAXX/newPipeline/mergedTagCounts/mergedPstIIBM_min50.cnt";
        String FastQFileName    = "H:/64GJAAAXX/newPipeline/mergedTagCounts/mergedPstIIBM_min50.fastq";
        TagCounts tc = new TagCounts();
        tc.toFASTQ(TagCountFileName, FastQFileName);
    }

    public static void runQseqToTagCountPlugin() {
//        String[] NAM49_50_ApeKIargs = new String[] {
//            "-i", "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/qseq",
//            "-k", "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/NAM49_50_ApeKI_key.txt",
//            "-e", "ApeKI", // Enzyme used to create the GBS library
//            "-s", "200000000", // Max good reads per lane. (Optional. Default is 200,000,000)
//            "-c", "1", // Minimum tag count (default is 1)
//            "-o", "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/tagCounts"
//        };

        String[] args = new String[] {
            "-i", "/cbsufsrv4/data1/maizediv/illumina/NAM_PstI_plates49_50/qseq",
            "-k", "/cbsufsrv4/data1/maizediv/illumina/NAM_PstI_plates49_50/70MU0AAXX_NAM_PstI_key.txt",
            "-e", "PstI", // Enzyme used to create the GBS library
            "-s", "200000000", // Max good reads per lane. (Optional. Default is 200,000,000)
            "-c", "1", // Minimum tag count (default is 1)
            "-o", "/cbsufsrv4/data1/maizediv/illumina/NAM_PstI_plates49_50/tagCounts"
        };

        QseqToTagCountPlugin plugin = new QseqToTagCountPlugin();
        plugin.setParameters(args);
        plugin.performFunction(null);
    }

    public static void runFastqToTagCountPlugin() {
        String testWorkdir = "H:/NAM_ApeKI_plates49_50/";
        String[] testArgs = new String[] {
            "-i", testWorkdir + "fakeFastq",
            "-k", testWorkdir + "NAM49_50_ApeKI_key2.txt",
            "-e", "ApeKI", // Enzyme used to create the GBS library
            "-s", "20000000", // Max good reads per lane. (Optional. Default is 200,000,000)
            "-c", "1", // Minimum tag count (default is 1)
            "-o", testWorkdir + "testFastq/tagCounts2",
        };

        String bigTestWorkdir = "/cbsufsrv4/data1/maizediv/illumina/";
        String[] bigTestArgs = new String[] {
            "-i", bigTestWorkdir + "B08AAABXX",
            "-k", bigTestWorkdir + "glaubitz_test/fake_B08AAABXX_key.txt",
            "-e", "ApeKI", // Enzyme used to create the GBS library
//            "-s", "20000000", // Max good reads per lane. (Optional. Default is 200,000,000)
//            "-c", "1", // Minimum tag count (default is 1)
            "-o", bigTestWorkdir + "glaubitz_test/testTagCounts",
        };
        
        String[] args = testArgs;
        FastqToTagCountPlugin plugin = new FastqToTagCountPlugin();
        plugin.setParameters(args);
        plugin.performFunction(null);
    }

    public static void runMergeMultipleTagCountPlugin() {
//        String[] NAM49_50_ApeKIargs = new String[] {
//            "-i", "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/tagCounts", // Input directory containing .cnt files
//            "-o", "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/mergedTagCounts/NAM49_50_ApeKI_mergedTags_min10.fastq", //  Output file name
//            "-c", "10" // Minimum count of reads to be output (default 1)
//          , "-t" // Specifies that reads should be output in FASTQ text format.
//        };

        String[] args = new String[] {
            "-i", "/cbsufsrv4/data1/maizediv/illumina/NAM_PstI_plates49_50/tagCounts", // Input directory containing .cnt files
            "-o", "/cbsufsrv4/data1/maizediv/illumina/NAM_PstI_plates49_50/mergedTagCounts/NAM49_50_PstI_mergedTags_min50.fastq", //  Output file name
            "-c", "50" // Minimum count of reads to be output (default 1)
          , "-t" // Specifies that reads should be output in FASTQ text format.
        };

        MergeMultipleTagCountPlugin plugin = new MergeMultipleTagCountPlugin();
        plugin.setParameters(args);
        plugin.performFunction(null);
    }

    public static void runQseqToTBTPlugin() {
        String[] NAM49_50_ApeKIargs = new String[] {
            "-i", "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/qseq",
            "-k", "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/NAM49_50_ApeKI_key.txt",
            "-e", "ApeKI", // Enzyme used to create the GBS library
            "-o", "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/tbt",
            "-c", "0", // Minimum tag count (default is 1).
            "-t", "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/mergedTagCounts/NAM49_50_ApeKI_mergedTags_min10.cnt", // master Tags file
        };

        String[] NAM49_50_PstIargs = new String[] {
            "-i", "/cbsufsrv4/data1/maizediv/illumina/NAM_PstI_plates49_50/qseq",
            "-k", "/cbsufsrv4/data1/maizediv/illumina/NAM_PstI_plates49_50/70MU0AAXX_NAM_PstI_key.txt",
            "-e", "PstI", // Enzyme used to create the GBS library
            "-o", "/cbsufsrv4/data1/maizediv/illumina/NAM_PstI_plates49_50/tbt",
            "-c", "0", // Minimum tag count (default is 1).
            "-t", "/cbsufsrv4/data1/maizediv/illumina/NAM_PstI_plates49_50/mergedTagCounts/NAM49_50_PstI_mergedTags_min50.cnt", // master Tags file
        };

        String[] KassaZea_args = new String[] {
            "-i", "/cbsufsrv4/data1/maizediv/illumina/KassaTemp/",
            "-k", "/cbsufsrv4/data1/maizediv/illumina/Kassa_Zea_key_20110810.txt",
            "-e", "ApeKI", // Enzyme used to create the GBS library
            "-o", "/cbsufsrv4/data1/maizediv/illumina/KassaMaize/tbt",
            "-c", "1", // Minimum tag count (default is 1).
            "-t", "/usr/local/maizediv/illumina/NAM_Ames_282/mergedNAM282Ames_qseq.cnt", // master Tags file
        };

        String[] KassaZeaB81RP7Args = new String[] {
            "-i", "/cbsufsrv4/data1/maizediv/illumina/KassaTemp/B81RP7ABXX/qseq/",  
            "-k", "/cbsufsrv4/data1/maizediv/illumina/Kassa_Zea_key_20110810.txt",
            "-e", "ApeKI", // Enzyme used to create the GBS library
            "-o", "/cbsufsrv4/data1/maizediv/illumina/KassaMaize/tbt",
            "-c", "1", // Minimum tag count (default is 1).
            "-t", "/usr/local/maizediv/illumina/NAM_Ames_282/mergedNAM282Ames_qseq.cnt", // master Tags file
        };

        String[] DTMATestArgs = new String[] {
            "-i", "/usr/local/maizediv/illumina/Zea/DTMA/qseq",   // contains a single HiSeq lane from DTMA (B08G4ABXX lane 1)
            "-k", "/usr/local/maizediv/illumina/Zea/DTMA/DTMA_1-3_Kassa_31-38_key.txt",
            "-e", "ApeKI", // Enzyme used to create the GBS library
            "-o", "/usr/local/maizediv/illumina/Zea/DTMA/tbtJG",
            "-c", "1", // Minimum tag count (default is 1).
            "-m", "/home/glaubitz/data/nextgen/allZea/topm/mergedNAM282Ames_072011.topm.bin", // master Tags file
        };

        String[] args = DTMATestArgs;
        QseqToTBTPlugin plugin = new QseqToTBTPlugin();
        plugin.setParameters(args);
        plugin.performFunction(null);
    }

    public static void runFastqToTBTPlugin() {

        String testFastqDir = "H:/NAM_ApeKI_plates49_50/";
        String[] testFastqArgs = new String[] {
            "-i", testFastqDir+"fakeFastq",  
            "-k", testFastqDir+"NAM49_50_ApeKI_key.txt",
            "-e", "ApeKI", // Enzyme used to create the GBS library
            "-o", testFastqDir+"testFastq/tbt",
            "-c", "1", // Minimum tag count (default is 1).
            "-t", testFastqDir+"testFastq/mergedTagCounts/testMerge.cnt", // master Tags file
        };

        String[] TeoDNA_P1_Args = new String[] {
            "-i", "/cbsufsrv4/data1/maizediv/illumina/GBSFlowCells/C08L7ACXX",  
            "-k", "/usr/local/maizediv/illumina/Zea/TeoDNA_P1/TeoDNA_P1_C08L7ACXX_key.txt",
            "-e", "ApeKI", // Enzyme used to create the GBS library
            "-o", "/usr/local/maizediv/illumina/Zea/TeoDNA_P1/tbt",
            "-c", "1", // Minimum tag count (default is 1).
            "-m", "/usr/local/maizediv/illumina/Zea/allZea_mappedonly_20120115.topm", // master Tags file (topm)
            "-y"  // output tbtByte
        };

        String[] C08L7ACXX_otherLanes_Args = new String[] {
            "-i", "/cbsufsrv4/data1/maizediv/illumina/GBSFlowCells/C08L7ACXX",  
            "-k", "/usr/local/maizediv/illumina/Zea/TeoDNA_P1/C08L7ACXX_otherLanes_key.txt",
            "-e", "ApeKI", // Enzyme used to create the GBS library
            "-o", "/usr/local/maizediv/illumina/Zea/TeoDNA_P1/tbt",
            "-c", "1", // Minimum tag count (default is 1).
            "-m", "/usr/local/maizediv/illumina/Zea/allZea_mappedonly_20120115.topm", // master Tags file (topm)
            "-y"  // output tbtByte
        };

        String[] args = C08L7ACXX_otherLanes_Args;
        FastqToTBTPlugin plugin = new FastqToTBTPlugin();
        plugin.setParameters(args);
        plugin.performFunction(null);
    }

    public static void runMergeTagsByTaxaFilesPlugin() {
        String[] NAM49_50_ApeKIargs = new String[] {
            "-i", "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/tbt/test",
            "-o", "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/mergedTBT/NAM49_50_ApeKI_mergedTBT_min10_20110720.tbt.bin",
        };

        String[] NAM49_50_PstIargs = new String[] {
            "-i", "/cbsufsrv4/data1/maizediv/illumina/NAM_PstI_plates49_50/tbt/",
            "-o", "/cbsufsrv4/data1/maizediv/illumina/NAM_PstI_plates49_50/mergedTBT/NAM49_50_PstI_mergedTBT_min50.tbt.bin",
        };

        String[] HapMap2Args = new String[] {
            "-i", "/usr/local/maizediv/illumina/NAM_Ames_282/HapMap2TBT/",
            "-o", "/usr/local/maizediv/illumina/NAM_Ames_282/HapMap2MergedTBT/HapMap2mergedTBT20110809.tbt.bin",
        };
        
        String[] testFastqArgs = new String[] {
            "-i", "H:/NAM_ApeKI_plates49_50/testFastq/tbt",
            "-o", "H:/NAM_ApeKI_plates49_50/testFastq/mergedTBT/testMerged2.tbt.bin",
        };

        String[] args = testFastqArgs;
        MergeTagsByTaxaFilesPlugin plugin = new MergeTagsByTaxaFilesPlugin();
        plugin.setParameters(args);
        plugin.performFunction(null);
    }

    public static void mergeTaxaInTBT() {
        String workDir = "/usr/local/maizediv/illumina/NAM_Ames_282/HapMap2MergedTBT/";
        String inputTBTFileS =            workDir + "HapMap2mergedTBT20110809.tbt.bin";
        String outputMergedTaxaTBTFileS = workDir + "HapMap2mergedTBT20110809.mergedTaxa.tbt.bin";
        TagsByTaxaUtils.mergeTaxaByName(inputTBTFileS, outputMergedTaxaTBTFileS, FilePacking.Bit, true);
        TagsByTaxaUtils.streamBinaryToText(outputMergedTaxaTBTFileS, 10000);
        TagsByTaxaUtils.printSumCounts(outputMergedTaxaTBTFileS, FilePacking.Bit, true);
    }

    public static void printSumCountsInTBTByTaxa() {
        String workDir = "/home/glaubitz/data/nextgen/allZea/tbt/";
        String inputTBTFileS = workDir + "allZea20110811.tbt.bin";
        TagsByTaxaUtils.printSumCounts(inputTBTFileS, FilePacking.Bit, true);
    }

    public static void convertSAMToTOPM() {
//        String SAMFile =      "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/mergedTagCounts/NAM49_50_ApeKI_mergedTags_min10.sam";
//        String TOPMFile =     "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/mergedTagCounts/NAM49_50_ApeKI_mergedTags_min10.topm.bin";
//        String TOPMTextFile = "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/mergedTagCounts/NAM49_50_ApeKI_mergedTags_min10.topm.txt";

        String SAMFile =      "H:/70MU0AAXX/newestPipeline/NAM49_50_PstI_mergedTags_min50.sam";
        String TOPMFile =     "H:/70MU0AAXX/newestPipeline/test/NAM49_50_PstI_mergedTags_min50.topm";
        String TOPMTextFile = "H:/70MU0AAXX/newestPipeline/test/NAM49_50_PstI_mergedTags_min50.topm.txt";

        TagsOnPhysicalMap topm = new TagsOnPhysicalMap();
        topm.readSAMFile(SAMFile,2);
        topm.writeBinaryFile(new File(TOPMFile));
        topm.writeTextFile(new File(TOPMTextFile));
    }

    public static void filterTOPMWithPositions() {
        String TOPMFile =     "N:/cassava/cassava.topm.bin";
        String TOPMFileWPos =     "N:/cassava/cassava_wPos.topm.bin";
        String TOPMTextFileWPos = "N:/cassava/cassava_wPos.topm.txt";

        TagsOnPhysicalMap topm = new TagsOnPhysicalMap(TOPMFile, true);
        topm.writeBinaryFile(new File(TOPMFileWPos), Integer.MAX_VALUE, true, true, (float) 1.0, true);
        topm = new TagsOnPhysicalMap(TOPMFileWPos, true);
        topm.writeTextFile(new File(TOPMTextFileWPos));
    }

    public static void runTagsToSNPByAlignmentMTPluginOld() {
//        String[] NAM49_50_ApeKIargs = new String[] {
//            "-i", "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/mergedTBT/NAM49_50_ApeKI_mergedTBT_min10_20110720.tbt.bin",
//            "-o", "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/hapmapOutput/maxAllelicTags100_rep2/",  
//            "-m", "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/mergedTagCounts/NAM49_50_ApeKI_mergedTags_min10.topm.bin",
//            "-mnLCov", "0.05", // Minimum locus coverage (proportion of Taxa)
//            "-s", "1",  // Start chromosome
//            "-e", "10" // End chromosome
//        };

        String[] args = new String[] {
            "-i", "/cbsufsrv4/data1/maizediv/illumina/NAM_PstI_plates49_50/mergedTBT/NAM49_50_PstI_mergedTBT_min50.tbt.bin",
            "-o", "/cbsufsrv4/data1/maizediv/illumina/NAM_PstI_plates49_50/hapmapOutput/",
            "-m", "/cbsufsrv4/data1/maizediv/illumina/NAM_PstI_plates49_50/mergedTagCounts/NAM49_50_PstI_mergedTags_min50.topm.bin",
            "-mnLCov", "0.1", // Minimum locus coverage (proportion of Taxa)
            "-s", "1",  // Start chromosome
            "-e", "10" // End chromosome
        };

        TagsToSNPByAlignmentMTPlugin plugin = new TagsToSNPByAlignmentMTPlugin();
        plugin.setParameters(args);
        plugin.performFunction(null);
    }

    public static void runTagsToSNPByAlignmentMTPlugin() {
        String[] allZeaArgs = new String[] {
            "-i", "/home/glaubitz/data/nextgen/allZea/tbt/allZea20110811.tbt.bin",
            "-o", "/usr/local/maizediv/illumina/allZeaHapMap",
            "-m", "/home/glaubitz/data/nextgen/allZea/topm/mergedNAM282Ames_072011.topm.bin",
            "-mnF",   "0.9",
            "-mnMAF", "0.005",
            "-mnMAC", "20",
            "-mnLCov","0.10", // Minimum locus coverage (proportion of Taxa)
//            "-inclGaps",  // Include sites where major or minor allele is a GAP
            "-s", "10",  // Start chromosome
            "-e", "10"  // End chromosome
        };

        String[] DTMAArgs = new String[] {
            "-i",    "/usr/local/maizediv/illumina/Zea/DTMA/tbtJG/B08G4ABXX_1.tbt.bin",
            "-o",    "/usr/local/maizediv/illumina/Zea/DTMA/answersJG",
            "-m",    "/home/glaubitz/data/nextgen/allZea/topm/mergedNAM282Ames_072011.topm.bin",
//            "-mUpd", "/usr/local/maizediv/illumina/Zea/DTMA/topm/mergedNAM282Ames_variants2_20111018.topm.bin",
            "-mnF",   "0.9",
            "-mnMAF", "0.03",
            "-mnMAC", "20",
            "-mnLCov","0.10", // Minimum locus coverage (proportion of Taxa)
//            "-inclGaps",  // Include sites where major or minor allele is a GAP
            "-s", "1",  // Start chromosome
            "-e", "10"  // End chromosome
        };
        
        String[] args = DTMAArgs;
        TagsToSNPByAlignmentMTPlugin plugin = new TagsToSNPByAlignmentMTPlugin();
        plugin.setParameters(args);
        plugin.performFunction(null);
    }

    public static void runMergeDuplicateSNPsPlugin() {

        String testWorkdir = "C:/Users/jcg233/Documents/Bioinformatics/NextGen/NAM_ApeKI_plates49_50/hapmapOutput/testMergeDuplicateSNPsPlugin/";
        String[] testArgs = new String[] {
            "-hmp", testWorkdir+"mergedTBT_testSNPMerge_inputTestEnd_c+.hmp.txt", // Input HapMap file; use a plus sign (+) as a wild card character to specify multiple chromosome numbers
            "-o",   testWorkdir+"mergedTBT_testSNPMerge_outputTestEnd_delUnmergDups_c+.hmp.txt", // Output HapMap file
            "-misMat", "0.15", // Threshold genotypic mismatch rate above which the duplicate SNPs won't be merged (default: 0.05)
            "-callHets",       // When two genotypes disagree at a SNP, call it a heterozygote (default: off = set to missing)
//          "-kpUnmergDups",   // When two duplicate SNPs were not merged (different alleles or too many mismatches), keep them (default: off = delete them)
            "-s", "10",        // Start chromosome (default 1)
            "-e", "10"         // End chromosome (default 10)
        };

        String workdir = "/usr/local/maizediv/illumina/Zea/hapmap/build112311/";
        String[] allZeaArgs = new String[] {
            "-hmp", workdir+    "unfilt/allZea_112311.c+.hmp.txt", // Input HapMap file; use a plus sign (+) as a wild card character to specify multiple chromosome numbers
            "-o",   workdir+"mergedSNPs/allZea_112311_SNPmerge15_c+.hmp.txt", // Output HapMap file
            "-misMat", "0.15", // Threshold genotypic mismatch rate above which the duplicate SNPs won't be merged (default: 0.05)
            "-callHets",       // When two genotypes disagree at a SNP, call it a heterozygote (default: off = set to missing)
//          "-kpUnmergDups",   // When two duplicate SNPs were not merged (different alleles or too many mismatches), keep them (default: off = delete them)
            "-s", "1",        // Start chromosome (default 1)
            "-e", "10"         // End chromosome (default 10)
        };

        String[] IBM94PstIArgs = new String[] {
            "-hmp",     "/usr/local/maizediv/illumina/Zea/PstI/hapmap/quantMnF80MinSCov80/filtLD/IBM94PstI_mnF80_filt_mnT10_mnS80_LD_c+.hmp.txt",
            "-o",       "/usr/local/maizediv/illumina/Zea/PstI/hapmap/quantMnF80MinSCov80/filtLD/mergedSNPs/IBM94PstI_mnF80_filt_mnT10_mnS80_LD_mergSNP20_c+.hmp.txt",
            "-misMat",  "0.2", // Threshold genotypic mismatch rate above which the duplicate SNPs won't be merged (default: 0.05)
            "-callHets",       // When two genotypes disagree at a SNP, call it a heterozygote (default: off = set to missing)
//            "-kpUnmergDups",   // When two duplicate SNPs were not merged (different alleles or too many mismatches), keep them (default: off = delete them)
            "-s", "1",        // Start chromosome (default 1)
            "-e", "10"         // End chromosome (default 10)
        };

        String[] IBM94PstITestPedArgs = new String[] {
            "-hmp",     "/usr/local/maizediv/illumina/Zea/PstI/hapmap/quantMnF80MinSCov80/filtLD/IBM94PstI_mnF80_filt_mnT10_mnS80_LD_c+.hmp.txt",
            "-o",       "/usr/local/maizediv/illumina/Zea/PstI/hapmap/quantMnF80MinSCov80/testFakePed/mergedSNPs/IBM94PstI_mnF80_filt_mnT10_mnS80_LD_mergSNP20_c+.hmp.txt",
            "-misMat",  "0.2", // Threshold genotypic mismatch rate above which the duplicate SNPs won't be merged (default: 0.05)
            "-p",       "/usr/local/maizediv/illumina/Zea/PstI/fakePedFile6_ped.txt",
            "-callHets",       // When two genotypes disagree at a SNP, call it a heterozygote (default: off = set to missing)
//            "-kpUnmergDups",   // When two duplicate SNPs were not merged (different alleles or too many mismatches), keep them (default: off = delete them)
            "-s", "1",        // Start chromosome (default 1)
            "-e", "10"         // End chromosome (default 10)
        };

        MergeDuplicateSNPsPlugin plugin = new MergeDuplicateSNPsPlugin();
        plugin.setParameters(IBM94PstITestPedArgs);
        plugin.performFunction(null);
    }

    public static void runGBSHapMapFiltersPlugin() {

        String workdir = "/usr/local/maizediv/illumina/allZeaHapMap/";

        String[] argsNAM49_50ApeKI = new String[] {
            "-hmp",    "H:/NAM_ApeKI_plates49_50/newestPipeline/hapmapOutput/20110726/mergedTBT.c+.hmp.txt", // Input HapMap file; use a plus sign (+) as a wild card character to specify multiple chromosome numbers
            "-o",      "H:/NAM_ApeKI_plates49_50/newestPipeline/hapmapOutput/20110726/filtered/NAM49_50_ApeKI_filtered.c+.hmp.txt", // Output HapMap file
            "-mnTCov", "0.05", // Minimum taxa coverage
            "-mnSCov", "0.05", // Minimum presence
            "-mnF",    "0.8",  // Minimum F (inbreeding coefficient) (default -2.0  - no filter)
            "-mnMAF",  "0.1",  // Minimum minor allele frequency (default 0.0)
            "-mxMAF",  "0.9",  // Maximum minor allele frequency (default 1.0 - no filter)
            "-hLD",            // Filter for high LD
            "-sC",     "1",    // Start chromosome (default 1)
            "-eC",     "10"    // End chromosome (default 10)"
        };

        String[] allZeaArgs = new String[] {
            "-hmp", workdir + "unfiltered/allZea20110812_unfiltSNPs_f9maf005mac20cov10_chr+.hmp.txt",
            "-o",   workdir +   "filtered/allZea20110812_filtF9maf001mac20LCov10_chr+.hmp.txt",
//            "-mnTCov", "0.0001",
            "-mnF",    "0.9",
            "-mnMAF",  "0.001",
//            "-hLD",
            "-mnSCov", "0.10", // Minimum locus coverage (proportion of Taxa)
            "-sC",     "8",    // Start chromosome
            "-eC",     "8"    // End chromosome
        };

        String[] BrunetArgs = new String[] {
            "-hmp", "C:/Users/jcg233/Documents/Bioinformatics/NextGen/UserSupport/Brunet/BC364mergedSNPs.c+.hmp.txt",
            "-o",   "C:/Users/jcg233/Documents/Bioinformatics/NextGen/UserSupport/Brunet/filt/BC364mergedSNPs_filt.c+.hmp.txt",
            "-mnTCov", "0.1",
//            "-mnF",    "0.9",
            "-mnMAF",  "0.0",
//            "-hLD",
            "-mnSCov", "0.1", // Minimum locus coverage (proportion of Taxa)
            "-sC",     "100",    // Start chromosome
            "-eC",     "101"    // End chromosome
        };

        String[] IBM94PstIArgs = new String[] {
            "-hmp",    "/usr/local/maizediv/illumina/Zea/PstI/hapmap/quantMnF80MinLCov10/filt/IBM94PstI_mnF80_filt_mnT10_mnS10_c+.hmp.txt",
            "-o",      "/usr/local/maizediv/illumina/Zea/PstI/hapmap/quantMnF80MinSCov80/filtLD/IBM94PstI_mnF80_filt_mnT10_mnS80_LD_c+.hmp.txt",
            "-mnTCov", "0.10", // Minimum taxa coverage (proportion of non-missing sites at a taxon)
            "-mnSCov", "0.80", // Minimum presence (proportion of non-missing taxa at a site)
            "-mnF",    "0.80",  // Minimum F (inbreeding coefficient) (default -2.0  - no filter)
            "-mnMAF",  "0.20",  // Minimum minor allele frequency (default 0.0)
//            "-mxMAF",  "0.9",  // Maximum minor allele frequency (default 1.0 - no filter)
            "-hLD",            // Filter for high LD
            "-sC",     "1",    // Start chromosome (default 1)
            "-eC",     "10"    // End chromosome (default 10)
        };

        String[] IBM94ApeKIArgs = new String[] {
            "-hmp",    "/usr/local/maizediv/illumina/Zea/IBM94ApeKI/hapmap/filt/IBM94ApeKI_mnF80_filt_mnT03_mnS05_c+.hmp.txt",
            "-o",      "/usr/local/maizediv/illumina/Zea/IBM94ApeKI/hapmap/filtLD/IBM94ApeKI_mnF80_filt_mnT03_mnS05_LD_c+.hmp.txt",
            "-mnTCov", "0.03", // Minimum taxa coverage (proportion of non-missing sites at a taxon)
            "-mnSCov", "0.05", // Minimum presence (proportion of non-missing taxa at a site)
            "-mnF",    "0.80",  // Minimum F (inbreeding coefficient) (default -2.0  - no filter)
            "-mnMAF",  "0.20",  // Minimum minor allele frequency (default 0.0)
//            "-mxMAF",  "0.9",  // Maximum minor allele frequency (default 1.0 - no filter)
            "-hLD",            // Filter for high LD
            "-sC",     "1",    // Start chromosome (default 1)
            "-eC",     "10"    // End chromosome (default 10)
        };

        String[] NIL28Args = new String[] {
            "-hmp",    "/usr/local/maizediv/illumina/Zea/NIL28/ProdHapmap/oldTOPM/NIL28_mgt80_chr+.hmp.txt",
            "-o",      "/usr/local/maizediv/illumina/Zea/NIL28/ProdHapmap/oldTOPM/NIL28_mgt80_filtMnT3mnS10mnF20mnMAF2_chr+.hmp.txt",
            "-mnTCov", "0.03", // Minimum taxa coverage (proportion of non-missing sites at a taxon)
            "-mnSCov", "0.10", // Minimum presence (proportion of non-missing taxa at a site)
            "-mnF",    "0.2",  // Minimum F (inbreeding coefficient) (default -2.0  - no filter)
            "-mnMAF",  "0.02",  // Minimum minor allele frequency (default 0.0)
//            "-mxMAF",  "0.9",  // Maximum minor allele frequency (default 1.0 - no filter)
//            "-hLD",            // Filter for high LD
            "-sC",     "3",    // Start chromosome (default 1)
            "-eC",     "3"    // End chromosome (default 10)
        };

        String[] NAM12Args = new String[] {
            "-hmp",    "/usr/local/maizediv/illumina/Zea/Production/NAM12/NAM12Build20120110_scv10mF8maf002_mgs_E1pLD5kpUn.c+.hmp.txt",
            "-o",      "/usr/local/maizediv/illumina/Zea/Production/NAM12/NAM12Build20120110_scv10mF8maf002_mgs_E1pLD5kpUn_maf20.c+.hmp.txt",
//            "-mnTCov", "0.03", // Minimum taxa coverage (proportion of non-missing sites at a taxon)
//            "-mnSCov", "0.10", // Minimum presence (proportion of non-missing taxa at a site)
//            "-mnF",    "0.2",  // Minimum F (inbreeding coefficient) (default -2.0  - no filter)
            "-mnMAF",  "0.2",  // Minimum minor allele frequency (default 0.0)
//            "-mxMAF",  "0.9",  // Maximum minor allele frequency (default 1.0 - no filter)
//            "-hLD",            // Filter for high LD
            "-sC",     "10",    // Start chromosome (default 1)
            "-eC",     "10"    // End chromosome (default 10)
        };

        String[] maize282Jan2012buildArgs = new String[] {
            "-hmp",    "/usr/local/maizediv/illumina/Zea/build20120110/bpec/maize282/maize282_Zea20120110_scv10mF8maf002_mgs_E1pLD5kpUn1.c+.hmp.txt",
            "-o",      "/usr/local/maizediv/illumina/Zea/build20120110/bpec/maize282/maize282_Zea20120110_scv10mF8maf002_mgs_E1pLD5kpUn1_tc5sc5maf1.c+.hmp.txt",
            "-mnTCov", "0.05", // Minimum taxa coverage (proportion of non-missing sites at a taxon)
            "-mnSCov", "0.05", // Minimum presence (proportion of non-missing taxa at a site)
//            "-mnF",    "0.2",  // Minimum F (inbreeding coefficient) (default -2.0  - no filter)
            "-mnMAF",  "0.01",  // Minimum minor allele frequency (default 0.0)
//            "-mxMAF",  "0.9",  // Maximum minor allele frequency (default 1.0 - no filter)
//            "-hLD",            // Filter for high LD
            "-sC",     "1",    // Start chromosome (default 1)
            "-eC",     "10"    // End chromosome (default 10)
        };

        String[] IBM94PstITestPedArgs = new String[] {
            "-hmp",    "/usr/local/maizediv/illumina/Zea/PstI/hapmap/quantMnF80MinLCov10/filt/IBM94PstI_mnF80_filt_mnT10_mnS10_c+.hmp.txt",
            "-o",      "/usr/local/maizediv/illumina/Zea/PstI/hapmap/quantMnF80MinSCov80/filtLD/fakePed/IBM94PstI_mnF80_filt_mnT10_mnS80_LD_c+.hmp.txt",
            "-mnTCov", "0.10", // Minimum taxa coverage (proportion of non-missing sites at a taxon)
            "-mnSCov", "0.80", // Minimum presence (proportion of non-missing taxa at a site)
            "-mnF",    "0.80",  // Minimum F (inbreeding coefficient) (default -2.0  - no filter)
            "-p",       "/usr/local/maizediv/illumina/Zea/PstI/fakePedFile6_ped.txt",
            "-mnMAF",  "0.20",  // Minimum minor allele frequency (default 0.0)
//            "-mxMAF",  "0.9",  // Maximum minor allele frequency (default 1.0 - no filter)
            "-hLD",            // Filter for high LD
            "-mnR2",   "0.50",  // Minimum minor allele frequency (default 0.0)
            "-mnBonP", "0.001",  // Minimum minor allele frequency (default 0.0)
            "-sC",     "1",    // Start chromosome (default 1)
            "-eC",     "10"    // End chromosome (default 10)
        };

        GBSHapMapFiltersPlugin plugin = new GBSHapMapFiltersPlugin();
        plugin.setParameters(IBM94PstITestPedArgs);
        plugin.performFunction(null);
    }

    public static void runBiParentalErrorCorrection() {

        String workdir = "/usr/local/maizediv/illumina/allZeaHapMap/CintaFM/";

        String[] args = new String[] {
            "-hmp", workdir + "CintaFM20110812_Tx303_Z025_maf20_chr+.hmp.txt",
            "-o",   workdir + "CintaFM20110812_Tx303_Z025_maf20_r50_chr+.hmp.txt",
            "-oB",  workdir + "errorBin.txt",
            "-oE",  workdir + "errorBySNP.txt",
            "-popM", "Z[0-9]{3}",
            "-sC",    "8",
            "-eC",    "8",
            "-mxE",   "0.01",
            "-mnD",   "2.0",
            "-mnPLD", "0.5",
        };

        BiParentalErrorCorrection.main(args);
    }

     public static void runMergeIdenticalTaxaPlugin() {
        String basedir = "/usr/local/maizediv/illumina/Zea/build20120110/bpec/TeoDNA_P1/";
        String[] TeoDNA_P1Args = new String[] {
            "-hmp", basedir+           "TeoDNA_P1_20120110_scv10mF8maf002_mgs_E1pLD5kpUn_chr+.hmp.txt",
            "-o",   basedir+"mergedTaxa/TeoDNA_P1_20120110_scv10mF8maf002_mgs_E1pLD5kpUn_mgt8_chr+.hmp.txt",
            "-hetFreq", "0.8", // if  ( (nMajGenos+nHetGenos)/(nMajGenos+nMinGenos+2*nHetGenos) > hetFreq ) geno = homMinor;  // likewise for homMinor (otherwise = het)
//            "-xHets",   // after application of above majority rule, hets set to missing
            "-sC", "1",  // Start chromosome
            "-eC", "10" // End chromosome
        };
        
        String[] NIL28Args = new String[] {
            "-hmp", "/usr/local/maizediv/illumina/Zea/NIL28/ProdHapmap/oldTOPM/NIL28_chr+.hmp.txt",
            "-o",   "/usr/local/maizediv/illumina/Zea/NIL28/ProdHapmap/oldTOPM/NIL28_mgt80_chr+.hmp.txt",
            "-hetFreq", "0.8", // if  ( (nMajGenos+nHetGenos)/(nMajGenos+nMinGenos+2*nHetGenos) > hetFreq ) geno = homMinor;  // likewise for homMinor (otherwise = het)
//            "-xHets",   // after application of above majority rule, hets set to missing
            "-sC", "3",  // Start chromosome
            "-eC", "3" // End chromosome
        };

        String[] maize282build20120110Args = new String[] {
            "-hmp",           "/usr/local/maizediv/illumina/Zea/build20120110/bpec/maize282/maize282_Zea20120110_scv10mF8maf002_mgs_E1pLD5kpUn1_tc5sc5maf1.c+.hmp.txt",
            "-o",   "/usr/local/maizediv/illumina/Zea/build20120110/bpec/maize282/taxaMerge/maize282_Zea20120110_scv10mF8maf002_mgs_E1pLD5kpUn1_tc5sc5maf1_mgT7.c+.hmp.txt",
            "-hetFreq", "0.7",  // if  ( (nMajGenos+nHetGenos)/(nMajGenos+nMinGenos+2*nHetGenos) > hetFreq ) geno = homMinor;  // likewise for homMinor (otherwise = het)
//            "-xHets",   // after application of above majority rule, hets set to missing
            "-sC", "1",  // Start chromosome
            "-eC", "10" // End chromosome
        };

        MergeIdenticalTaxaPlugin mitp = new MergeIdenticalTaxaPlugin();
        mitp.setParameters(maize282build20120110Args);
        mitp.performFunction(null);
    }

    public static void runBinaryToTextPlugin() {
        String NAM_ApeKI_plates49_50Dir = "H:/NAM_ApeKI_plates49_50/newPipeline/";
        String[] TagCountsArgs = new String[] {
            "-i", NAM_ApeKI_plates49_50Dir+"mergedTagCounts/mergedApeKINAM49_50_min10.cnt",
            "-o", NAM_ApeKI_plates49_50Dir+"mergedTagCounts/mergedApeKINAM49_50_min10_testB2Tplugin.cnt.txt",
            "-t", "TagCounts",  // TOPM, TagCounts, TBTBit
        };

        String[] TBTArgs = new String[] {
            "-i", NAM_ApeKI_plates49_50Dir+"mergedTagsByTaxa/ApekI_NAM49_50_merged_tbt.bin",
            "-o", NAM_ApeKI_plates49_50Dir+"mergedTagsByTaxa/ApekI_NAM49_50_merged_tbt.bin.txt",
            "-t", "TBTBit",  // TOPM, TagCounts, TBTBit
        };

        String[] PstIMin5Args = new String[] {
            "-i", "/usr/local/maizediv/illumina/Zea/PstI/topm/B08AAABXX_1_IBM94PstI_min5.topm.bin",
            "-o", "/usr/local/maizediv/illumina/Zea/PstI/topm/B08AAABXX_1_IBM94PstI_min5.topm.txt",
            "-t", "TOPM",  // TOPM, TagCounts, TBTBit
        };

        String[] TestArgs = new String[] {
            "-i", "H:/NAM_ApeKI_plates49_50/testFastq/mergedTBT/testMerged2.tbt.bin",
            "-o", "H:/NAM_ApeKI_plates49_50/testFastq/mergedTBT/testMerged2.tbt.txt",
            "-t", "TBTBit",  // TOPM, TagCounts, TBTBit
        };

        String[] ProdTOPM20120110Args = new String[] {
            "-i", "/usr/local/maizediv/illumina/Zea/build20120110/topm/zea20120110c1-4.prod1-4.topm",
            "-o", "/usr/local/maizediv/illumina/Zea/build20120110/topm/zea20120110c1-4.prod1-4.topm.txt",
            "-t", "TOPM",  // TOPM, TagCounts, TBTBit
        };

        String[] PearlMilletTOPMArgs = new String[] {
            "-i", "/usr/local/maizediv/illumina/JasonW_PM_GBS/3_AlignedTags.topm.bin",
            "-o", "/usr/local/maizediv/illumina/JasonW_PM_GBS/3_AlignedTags.topm.txt",
            "-t", "TOPM",  // TOPM, TagCounts, TBTBit
        };

        BinaryToTextPlugin plugin = new BinaryToTextPlugin(null);
        plugin.setParameters(PearlMilletTOPMArgs);
        plugin.performFunction(null);
    }

    public static void filterMergeAllZea() {
        String workDir = "/usr/local/maizediv/illumina/Zea/hapmap/build111123/";
        String filePrefix = "allZea_111123_SNPmerge15";
        int sC=1;
        int eC=9;

        String[] args = new String[] {
            "-hmp", workDir+"mergedSNPs/"+filePrefix+"_c+.hmp.txt",
            "-o",       workDir+"filt/"+filePrefix+"_cov10_fT1_c+.hmp.txt",
            "-mnTCov", "0.01",
            "-mnF",    "0.9",
            "-mnMAF",  "0.002",
//            "-hLD",
            "-mnSCov", "0.10", // Minimum locus coverage (proportion of Taxa)
            "-sC", ""+sC,  // Start chromosome
            "-eC", ""+eC // End chromosome
        };
        GBSHapMapFiltersPlugin testClass = new GBSHapMapFiltersPlugin();
        testClass.setParameters(args);
        testClass.performFunction(null);

        args = new String[] {
            "-hmp", workDir+"filt/"+filePrefix+"_cov10_fT1_c+.hmp.txt",
            "-o",     workDir+"bpec/"+filePrefix+"_cov10_fT1E1pLD_c+.hmp.txt",
            "-oB", "usr/local/maizediv/illumina/allZeaHapMap/build111003/bpec/errorBin.txt",
            "-oE", "usr/local/maizediv/illumina/allZeaHapMap/build111003/bpec/errorBySNP.txt",
            "-popM", "Z[0-9]{3}",
            "-sC",    ""+sC,
            "-eC",    ""+eC,
            "-mxE",   "0.01",
            "-mnD",   "2.0",
            "-mnPLD", "0.5"
        };
        BiParentalErrorCorrectionPlugin bpec = new BiParentalErrorCorrectionPlugin();
        bpec.setParameters(args);
        bpec.performFunction(null);

        args = new String[] {
            "-hmp",      workDir+"bpec/"+filePrefix+"_cov10_fT1E1pLD_c+.hmp.txt",
            "-o",  workDir+"mergedTaxa/"+filePrefix+"_cov10_fT1E1pLD_mergedTaxa_c+.hmp.txt",
            "-xHets",
            "-hetFreq", "0.76",
            "-sC", ""+sC,  // Start chromosome
            "-eC", ""+eC // End chromosome
        };
        MergeIdenticalTaxaPlugin mitp = new MergeIdenticalTaxaPlugin();
        mitp.setParameters(args);
        mitp.performFunction(null);
    }

     public static void imputeAllZea() {
        int sC=2;
        int eC=2;
        String stemIn =  "/usr/local/maizediv/illumina/Zea/hapmap/build111123/mergedTaxa/allZea_112311_SNPmerge15";
        String stemOut = "/usr/local/maizediv/illumina/Zea/hapmap/build111123/imputed/allZea_112311_SNPmerge15";
        String[] args = new String[] {
            "-hmp", stemIn+"_cov10_fT1E1pLD_mergedTaxa_c+.hmp.txt",
            "-o", stemOut+"_cov10_fT1E1pLD_mergedTaxa_imputed_c+.hmp.txt",
            "-sC", ""+sC,  // Start chromosome
            "-eC", ""+eC // End chromosome
        };
        FastImputationBitFixedWindow.main(args);
    }

    public static void runQseqToHapMapPlugin() {

        String[] DTMAArgs = new String[] {
            "-i", "/usr/local/maizediv/illumina/Zea/DTMA/qseq",   // contains a single HiSeq lane from DTMA (B08G4ABXX lane 1)
            "-k", "/usr/local/maizediv/illumina/Zea/DTMA/DTMA_1-3_Kassa_31-38_key.txt",
            "-e", "ApeKI", // Enzyme used to create the GBS library
            "-o", "/usr/local/maizediv/illumina/Zea/DTMA/qseq2hapmap",
            "-m", "/usr/local/maizediv/illumina/Zea/DTMA/topm/mergedNAM282Ames_variants2_20111018.topm.bin", // master TOPM file with variants recorded from discovery phase
//            "-c", "1", // Minimum tag count (default is 1)
//            "-d", "0",  // Maximum divergence between new read and previously mapped read (default = 0)
        };

        String[] args = DTMAArgs;
        QseqToHapMapPlugin plugin = new QseqToHapMapPlugin();
        plugin.setParameters(args);
        plugin.performFunction(null);
    }

    public static void runRawReadsToHapMapPlugin() {

        String[] DTMAArgs = new String[] {
            "-i", "/usr/local/maizediv/illumina/Zea/DTMA/qseq",   // contains a single HiSeq lane from DTMA (B08G4ABXX lane 1)
            "-k", "/usr/local/maizediv/illumina/Zea/DTMA/DTMA_1-3_Kassa_31-38_key.txt",
            "-e", "ApeKI", // Enzyme used to create the GBS library
            "-o", "/usr/local/maizediv/illumina/Zea/DTMA/rawReads2hapmap",
            "-m", "/usr/local/maizediv/illumina/Zea/DTMA/topm/mergedNAM282Ames_variants2_20111018.topm.bin", // master TOPM file with variants recorded from discovery phase
//            "-d", "0",  // Maximum divergence between new read and previously mapped read (default = 0)
        };

        String[] NIL28Args = new String[] {
            "-i", "/usr/local/maizediv/illumina/Zea/NIL28/fastq",   // contains symbolic links to the 4 lanes of data
            "-k", "/usr/local/maizediv/illumina/Zea/NIL28/NIL28BarcodeKey.txt",
            "-e", "ApeKI", 
            "-o", "/usr/local/maizediv/illumina/Zea/NIL28/ProdHapmap/oldTOPM",
            "-m", "/usr/local/maizediv/illumina/Zea/DTMA/topm/mergedNAM282Ames_variants2_20111018.topm.bin", // old TOPM with variants recorded from discovery phase (NOT COMPLEMENTED)
//            "-m", "/usr/local/maizediv/illumina/Zea/build20120110/topm/zea20120110c1-4.prod1-4.topm", // master TOPM file with variants recorded from discovery phase
//            "-d", "0",  // Maximum divergence between new read and previously mapped read (default = 0)
        };

        String[] NAM12Args = new String[] {
            "-i", "/usr/local/maizediv/illumina/Zea/Production/NAM12/fastq",   // contains a single HiSeq lane from test of new 96 plex barcodes on NAM12 (C0EBNACXX lane 1)
            "-k", "/usr/local/maizediv/illumina/Zea/Production/NAM12/C0EBNACXX_1_key.txt",
            "-e", "ApeKI", // Enzyme used to create the GBS library
            "-o", "/usr/local/maizediv/illumina/Zea/Production/NAM12",
            "-m", "/usr/local/maizediv/illumina/Zea/build20120110/topm/bpec_filtered_042012.topm", // master TOPM file with variants recorded from discovery phase & filtered after bpec
//            "-d", "0",  // Maximum divergence between new read and previously mapped read (default = 0)
        };

        String[] args = NAM12Args;
        RawReadsToHapMapPlugin plugin = new RawReadsToHapMapPlugin();
        plugin.setParameters(args);
        plugin.performFunction(null);
    }

//    public static void runRawReadsToHapMapQuantPlugin() {
//        String[] NAM12Args = new String[] {
//            "-i", "/usr/local/maizediv/illumina/Zea/Production/NAM12/fastq",   // contains a single HiSeq lane from test of new 96 plex barcodes on NAM12 (C0EBNACXX lane 1)
//            "-k", "/usr/local/maizediv/illumina/Zea/Production/NAM12/C0EBNACXX_1_key.txt",
//            "-e", "ApeKI", // Enzyme used to create the GBS library
//            "-o", "/usr/local/maizediv/illumina/Zea/Production/NAM12/quant",
//            "-m", "/usr/local/maizediv/illumina/Zea/build20120110/topm/bpec_filtered_042012.topm", // master TOPM file with variants recorded from discovery phase & filtered after bpec
////            "-d", "0",  // Maximum divergence between new read and previously mapped read (default = 0)
//        };
//
//        String[] args = NAM12Args;
//        RawReadsToHapMapQuantPlugin plugin = new RawReadsToHapMapQuantPlugin();
//        plugin.setParameters(args);
//        plugin.performFunction(null);
//    }

    public static void analyzeD09FYACXX_5_PstIMaize() {
        String[] args;
        String baseDir = "/usr/local/maizediv/illumina/Zea/PstI/";

        String[] FastqToTagCountArgs = new String[] {
            "-i", baseDir+"fastq",
            "-k", baseDir+"D09FYACXX_5_key.txt",
            "-e", "PstI", // Enzyme used to create the GBS library
//            "-s", "200000000", // Max good reads per lane. (Optional. Default is 200,000,000)
            "-c", "50", // Minimum tag count (default is 1).
            "-o", baseDir+"tagCounts",
        };
        args = FastqToTagCountArgs;
        FastqToTagCountPlugin pluginFastqToTagCount = new FastqToTagCountPlugin();
        pluginFastqToTagCount.setParameters(args);
        pluginFastqToTagCount.performFunction(null);

        String[] FastqToTBTArgs = new String[] {
            "-i", baseDir+"fastq",
            "-k", baseDir+"D09FYACXX_5_key.txt",
            "-e", "PstI", // Enzyme used to create the GBS library
            "-o", baseDir+"tbt",
            "-c", "1", // Minimum tag count (default is 1).
            "-t", baseDir+"tagCounts/D09FYACXX_5.cnt", // master Tags file (only one lane)
        };
        args = FastqToTBTArgs;
        FastqToTBTPlugin pluginFastqToTBT = new FastqToTBTPlugin();
        pluginFastqToTBT.setParameters(args);
        pluginFastqToTBT.performFunction(null);

        String[] MergeMultipleTagCountArgs = new String[] {
            "-i", baseDir+"tagCounts", // Input directory containing .cnt files
            "-o", baseDir+"mergedTagCounts/D09FYACXX_5_PstI_min50.fastq", //  Output file name
            "-c", "50" // Minimum count of reads to be output (default 1)
          , "-t" // Specifies that reads should be output in FASTQ text format.
        };
        args = MergeMultipleTagCountArgs;
        MergeMultipleTagCountPlugin pluginMergeMultipleTagCount = new MergeMultipleTagCountPlugin();
        pluginMergeMultipleTagCount.setParameters(args);
        pluginMergeMultipleTagCount.performFunction(null);
        
        // next is BWA...
    }

    public static void analyzeB08AAABXX_1_PstI_IBM() {
        String baseDir = "/usr/local/maizediv/illumina/Zea/PstI/";

        String[] FastqToTagCountArgs = new String[] {
            "-i", baseDir+"fastq",
            "-k", baseDir+"B08AAABXX_1_IBM94PstI_key.txt",
            "-e", "PstI", // Enzyme used to create the GBS library
//            "-s", "200000000", // Max good reads per lane. (Optional. Default is 200,000,000)
            "-c", "50", // Minimum tag count (default is 1).
            "-o", baseDir+"tagCounts",
        };
        FastqToTagCountPlugin pluginFastqToTagCount = new FastqToTagCountPlugin();
        pluginFastqToTagCount.setParameters(FastqToTagCountArgs);
        pluginFastqToTagCount.performFunction(null);

        String[] FastqToTBTArgs = new String[] {
            "-i", baseDir+"fastq",
            "-k", baseDir+"B08AAABXX_1_IBM94PstI_key.txt",
            "-e", "PstI", // Enzyme used to create the GBS library
            "-o", baseDir+"tbt",
            "-c", "1", // Minimum tag count (default is 1).
            "-t", baseDir+"tagCounts/B08AAABXX_1.cnt", // master Tags file (only one lane)
        };
        FastqToTBTPlugin pluginFastqToTBT = new FastqToTBTPlugin();
        pluginFastqToTBT.setParameters(FastqToTBTArgs);
        pluginFastqToTBT.performFunction(null);

        String TagCountFileName = baseDir+"tagCounts/B08AAABXX_1.cnt";
        String FastQFileName    = baseDir+"tagCounts/B08AAABXX_1_IBM94PstI_min50.fastq";
        TagCounts tc = new TagCounts();
        tc.toFASTQ(TagCountFileName, FastQFileName);
        
        String[] SAMConverterPluginArgs = new String[] {
            "-i", baseDir+"tagCounts/B08AAABXX_1_IBM94PstI_min50.sam",
            "-o", baseDir+     "topm/B08AAABXX_1_IBM94PstI_min50.topm.bin",
        };
        SAMConverterPlugin pluginSAMConverter = new SAMConverterPlugin();
        pluginSAMConverter.setParameters(SAMConverterPluginArgs);
        pluginSAMConverter.performFunction(null);

        String[] TagsToSNPByAlignmentArgs = new String[] {
            "-i",    baseDir+"tbt/B08AAABXX_1.tbt.bin",
            "-o",    baseDir+"hapmap/minCount50",
            "-m",    baseDir+"topm/B08AAABXX_1_IBM94PstI_min50.topm.bin",
//            "-mUpd", baseDir+"",
            "-mnF",   "0.8",  // allow some more hets in IBM
            "-mnMAF", "0.2",
            "-mnMAC", "90",  // this will never be satified: this way -mnMAF overrides it
            "-mnLCov","0.80", // Minimum locus coverage (proportion of Taxa)
//            "-inclGaps",  // Include sites where major or minor allele is a GAP
            "-s", "1",  // Start chromosome
            "-e", "10"  // End chromosome
        };
        TagsToSNPByAlignmentMTPlugin pluginTagsToSNPByAlignment = new TagsToSNPByAlignmentMTPlugin();
        pluginTagsToSNPByAlignment.setParameters(TagsToSNPByAlignmentArgs);
        pluginTagsToSNPByAlignment.performFunction(null);
    }

    public static void analyzeB08AAABXX_1_PstI_IBM_TBTByte() {
        String baseDir = "/usr/local/maizediv/illumina/Zea/PstI/";

        String[] FastqToTBTArgs = new String[] {
            "-i", baseDir+"fastq",
            "-k", baseDir+"B08AAABXX_1_IBM94PstI_key.txt",
            "-e", "PstI", // Enzyme used to create the GBS library
            "-o", baseDir+"tbtByte",
            "-c", "1", // Minimum tag count (default is 1).
            "-y", // use TagsByTaxaByte
            "-t", baseDir+"tagCounts/B08AAABXX_1.cnt", // master Tags file (only one lane)
        };
        FastqToTBTPlugin pluginFastqToTBT = new FastqToTBTPlugin();
        pluginFastqToTBT.setParameters(FastqToTBTArgs);
        pluginFastqToTBT.performFunction(null);

        // convery the binary tbt.byte file to a text version
        String tbtByteFile = baseDir+"tbtByte/B08AAABXX_1.tbt.byte";
        TagsByTaxa tbtByteIBMPstI = new TagsByTaxaByte(tbtByteFile, FilePacking.Byte);
        String tbtByteTextFile = baseDir+"tbtByte/B08AAABXX_1.tbt.txt";
        tbtByteIBMPstI.writeDistFile(new File(tbtByteTextFile), FilePacking.Text, 1);

        String[] TagsToSNPByAlignmentArgs = new String[] {
            "-i",    baseDir+"tbtByte/B08AAABXX_1.tbt.byte",
            "-y", // use TagsByTaxaByte
            "-o",    baseDir+"hapmap/quantMnF80MinLCov10",
            "-m",    baseDir+"topm/B08AAABXX_1_IBM94PstI_min50.topm.bin",
//            "-mUpd", baseDir+"",
            "-mnF",   "0.8",  // allow some more hets in IBM (tried 0.8 initially, then 0.6, 0.7, and 0.8 again)
            "-mnMAF", "0.2",
            "-mnMAC", "99999",  // this will never be satified: this way -mnMAF overrides it
            "-mnLCov","0.1", // Minimum locus coverage (proportion of Taxa)
//            "-inclGaps",  // Include sites where major or minor allele is a GAP
//            "-callBiSNPsWGap",  //call sites with a biallelic SNP plus a gap (e.g., A/C/-)
            "-s", "1",  // Start chromosome
            "-e", "10"  // End chromosome
        };
        String[] args = TagsToSNPByAlignmentArgs;
        TagsToSNPByAlignmentPlugin pluginTagsToSNPByAlignment = new TagsToSNPByAlignmentPlugin();
        pluginTagsToSNPByAlignment.setParameters(args);
        pluginTagsToSNPByAlignment.performFunction(null);
    }

    public static void analyzeIBM94ApeKI() {
        String baseDir = "/usr/local/maizediv/illumina/Zea/IBM94ApeKI/";

        String[] QseqToTagCountArgs = new String[] {
            "-i", baseDir+"qseq",
            "-k", baseDir+"C05F2ACXX_5_key.txt",
            "-e", "ApeKI", // Enzyme used to create the GBS library
//            "-s", "200000000", // Max good reads per lane. (Optional. Default is 200,000,000)
            "-c", "5", // Minimum tag count (default is 1).
            "-o", baseDir+"tagCounts",
        };
        QseqToTagCountPlugin pluginQseqToTagCount = new QseqToTagCountPlugin();
        pluginQseqToTagCount.setParameters(QseqToTagCountArgs);
        pluginQseqToTagCount.performFunction(null);

        String[] QseqToTBTArgs = new String[] {
            "-i", baseDir+"qseq",
            "-k", baseDir+"C05F2ACXX_5_key.txt",
            "-e", "ApeKI", // Enzyme used to create the GBS library
            "-o", baseDir+"tbt",
            "-c", "1", // Minimum tag count (default is 1).
            "-y", // use TagsByTaxaByte
            "-t", baseDir+"tagCounts/C05F2ACXX_5.cnt", // master Tags file (only one lane)
        };
        QseqToTBTPlugin pluginQseqToTBT = new QseqToTBTPlugin();
        pluginQseqToTBT.setParameters(QseqToTBTArgs);
        pluginQseqToTBT.performFunction(null);

        String TagCountFileName = baseDir+"tagCounts/C05F2ACXX_5.cnt";
        String FastQFileName    = baseDir+"tagCounts/C05F2ACXX_5_IBM94ApeKI_min5.fastq";
        TagCounts tc = new TagCounts();
        tc.toFASTQ(TagCountFileName, FastQFileName);

        // Comment out the next two steps until BWA has been run
        // At that point, comment out the above steps
        
        String[] SAMConverterPluginArgs = new String[] {
            "-i", baseDir+"tagCounts/C05F2ACXX_5_IBM94ApeKI_min5.sam",
            "-o", baseDir+     "topm/C05F2ACXX_5_IBM94ApeKI_min5.topm.bin",
        };
        SAMConverterPlugin pluginSAMConverter = new SAMConverterPlugin();
        pluginSAMConverter.setParameters(SAMConverterPluginArgs);
        pluginSAMConverter.performFunction(null);

        String[] TagsToSNPByAlignmentArgs = new String[] {
            "-i",    baseDir+"tbt/C05F2ACXX_5.tbt.byte",
            "-y", // use TagsByTaxaByte
            "-o",    baseDir+"hapmap/unfilt",
            "-m",    baseDir+"topm/C05F2ACXX_5_IBM94ApeKI_min5.topm.bin",
//            "-mUpd", baseDir+"",
            "-mnF",   "0.8",  // allow some more hets in IBM
            "-mnMAF", "0.2",
            "-mnMAC", "99999",  // this will never be satified: this way -mnMAF overrides it
            "-mnLCov","0.05", // Minimum locus coverage (proportion of Taxa)
//            "-inclGaps",  // Include sites where major or minor allele is a GAP
//            "-callBiSNPsWGap",  //call sites with a biallelic SNP plus a gap (e.g., A/C/-)
            "-s", "1",  // Start chromosome
            "-e", "10"  // End chromosome
        };
        TagsToSNPByAlignmentPlugin pluginTagsToSNPByAlignment = new TagsToSNPByAlignmentPlugin();
        pluginTagsToSNPByAlignment.setParameters(TagsToSNPByAlignmentArgs);
        pluginTagsToSNPByAlignment.performFunction(null);
    }

    public static void analyzeD0E3PACXX_5_PstI_PalmarChico() {
        String baseDir = "/Users/jcg233/Documents/GBS/Zea/PstI_PalmarChico/";

        String[] FastqToTagCountArgs = new String[] {
            "-i", baseDir+"fastq",
            "-k", baseDir+"D0E3PACXX_5_barcode_key.txt",
            "-e", "PstI", // Enzyme used to create the GBS library
            "-o", baseDir+"tagCounts",
            "-c", "3", // Minimum tag count (default is 1).
        };
        FastqToTagCountPlugin pluginFastqToTagCount = new FastqToTagCountPlugin();
        pluginFastqToTagCount.setParameters(FastqToTagCountArgs);
        pluginFastqToTagCount.performFunction(null);

        String[] TagCountToFastqArgs = new String[] {
            "-i", baseDir+"tagCounts/D0E3PACXX_5.cnt",
            "-o", baseDir+"tagCounts/PstI_PalmarChico_tags_min5.fq",
            "-c", "5", // Minimum tag count (default is 1).
        };
        TagCountToFastqPlugin pluginTagCountToFastq = new TagCountToFastqPlugin();
        pluginTagCountToFastq.setParameters(TagCountToFastqArgs);
        pluginTagCountToFastq.performFunction(null);


        // BWA steps:
//        bwa aln -t 8 /usr/local/maizediv/genome/maize_agp_v2.fasta PstI_PalmarChico_tags_min5.fq > PstI_PalmarChico_tags_min5.sai
//        bwa samse /usr/local/maizediv/genome/maize_agp_v2.fasta PstI_PalmarChico_tags_min5.sai PstI_PalmarChico_tags_min5.fq > PstI_PalmarChico_tags_min5.sam
        
        
        String[] SAMConverterArgs = new String[] {
            "-i", baseDir+"tagCounts/PstI_PalmarChico_tags_min5.sam",
            "-o", baseDir+     "topm/PstI_PalmarChico_min5.topm",
        };
        SAMConverterPlugin pluginSAMConverter = new SAMConverterPlugin();
        pluginSAMConverter.setParameters(SAMConverterArgs);
        pluginSAMConverter.performFunction(null);

        String[] FastqToTBTArgs = new String[] {
            "-i", baseDir+"fastq",
            "-k", baseDir+"D0E3PACXX_5_barcode_key.txt",
            "-e", "PstI", // Enzyme used to create the GBS library
            "-o", baseDir+"tbt",
            "-c", "1", // Minimum tag count (default is 1).
            "-y", // use TagsByTaxaByte
            "-m", baseDir+"topm/PstI_PalmarChico_min5.topm", // master Tags file (topm)
        };
        FastqToTBTPlugin pluginFastqToTBT = new FastqToTBTPlugin();
        pluginFastqToTBT.setParameters(FastqToTBTArgs);
        pluginFastqToTBT.performFunction(null);

//         convert the binary tbt.byte file to a text version
        String tbtByteFile = baseDir+"tbt/D0E3PACXX_5.tbt.byte";
        TagsByTaxa tbtBytePCPstI = new TagsByTaxaByte(tbtByteFile, FilePacking.Byte);
        String tbtByteTextFile = baseDir+"tbt/D0E3PACXX_5.tbt.txt";
        tbtBytePCPstI.writeDistFile(new File(tbtByteTextFile), FilePacking.Text, 1);

        String[] TagsToSNPByAlignmentArgs = new String[] {
            "-i",    baseDir+"tbt/D0E3PACXX_5.tbt.byte",
            "-y", // use TagsByTaxaByte
            "-o",    baseDir+"hapmap/noFilters",
            "-m",    baseDir+"topm/PstI_PalmarChico_min5.topm",
//            "-mUpd", baseDir+"",
//            "-mnF",   " -0.2",  // used -0.2 for the initial, filtered version
//            "-mnMAF", "0.01",   // used 0.01 for the initial, filtered version
            "-mnMAC", "1",  //mnMAC of 1 ensures that all are SNPs   // used 99999 for the initial, filtered version (=irrelevant)
            "-mnLCov","0.1", // Minimum locus coverage (proportion of Taxa)
//            "-inclGaps",  // Include sites where major or minor allele is a GAP
//            "-callBiSNPsWGap",  //call sites with a biallelic SNP plus a gap (e.g., A/C/-)
            "-s", "1",  // Start chromosome
            "-e", "10"  // End chromosome
        };
        TagsToSNPByAlignmentPlugin pluginTagsToSNPByAlignment = new TagsToSNPByAlignmentPlugin();
        pluginTagsToSNPByAlignment.setParameters(TagsToSNPByAlignmentArgs);
        pluginTagsToSNPByAlignment.performFunction(null);

        // filter out the blank & the partially failed N19-53 (this will also sort the taxa in the output)
        String taxaListFileName = baseDir+"hapmap/filt/keepTaxa/PalmarChicoPstIKeepTaxa.txt";
        String infileS =          baseDir+"hapmap/noFilters/tbt.c+.hmp.txt";
        String outfileS =         baseDir+"hapmap/noFilters/keepTaxa/PalmarChicoPstI_c5_MAC1_94taxa_chr+.hmp.txt";
        String infile, outfile;
        TasselPipeline tp = null;
        String[] keepTaxaArgs;
        int sC=1;
        int eC=10;
        for (int chr=sC; chr<=eC; ++chr) {
            infile=infileS.replace("+", ""+chr);
            outfile=outfileS.replace("+", ""+chr);
            keepTaxaArgs = new String[] {
                "-fork1",
                "-h", infile,
                "-includeTaxaInFile", taxaListFileName,
                "-taxaJoinStrict", "true",
                "-export", outfile,
                "-exportType", "Hapmap",
                "-runfork1",
            };
            tp = new TasselPipeline(keepTaxaArgs, null);
        }
        
        String[] mergeDupSNPsArgs = new String[] {
            "-hmp",     baseDir+  "hapmap/noFilters/keepTaxa/PalmarChicoPstI_c5_MAC1_94taxa_chr+.hmp.txt",
            "-o",       baseDir+"hapmap/noFilters/mergedSNPs/PalmarChicoPstI_c5_MAC1_94taxa_mgSNPs10_chr+.hmp.txt",
            "-misMat",  "0.1", // Threshold genotypic mismatch rate above which the duplicate SNPs won't be merged (default: 0.05)
            "-callHets",       // When two genotypes disagree at a SNP, call it a heterozygote (default: off = set to missing)
//            "-kpUnmergDups",   // When two duplicate SNPs were not merged (different alleles or too many mismatches), keep them (default: off = delete them)
            "-s", "1",        // Start chromosome (default 1)
            "-e", "10"         // End chromosome (default 10)
        };
        MergeDuplicateSNPsPlugin pluginMergeDuplicateSNPs = new MergeDuplicateSNPsPlugin();
        pluginMergeDuplicateSNPs.setParameters(mergeDupSNPsArgs);
        pluginMergeDuplicateSNPs.performFunction(null);

        String[] GBSHapMapFiltersArgs = new String[] {
            "-hmp",    baseDir+"hapmap/noFilters/mergedSNPs/PalmarChicoPstI_c5_MAC1_94taxa_mgSNPs10_chr+.hmp.txt",
            "-o",      baseDir+      "hapmap/noFilters/filtCovPoly/PalmarChicoPstI_c5_MAC1_94taxa_mgSNPs10_chr+.hmp.txt",
            "-mnSCov", "0.9",  // Minimum SNP coverage (proportion of Taxa)
//            "-mnTCov", "0.01",
//            "-mnF",    " -0.2",
            "-mnMAF",  "0.0001",  // must be polymorphic among the 94 remaining taxa
//            "-hLD",
            "-sC",     "1",     // Start chromosome
            "-eC",     "10"      // End chromosome
        };
        GBSHapMapFiltersPlugin pluginGBSHapMapFilters = new GBSHapMapFiltersPlugin();
        pluginGBSHapMapFilters.setParameters(GBSHapMapFiltersArgs);
        pluginGBSHapMapFilters.performFunction(null);
        
//        // filter again because taxa removed (blank and N19-53) made some SNPs monomorphic  // not necessary for 2nd, unfiltered version b/c I removed taxa right after SNP calling
//        String[] GBSHapMapFilters2Args = new String[] {
//            "-hmp",    baseDir+"hapmap/filt/keepTaxa/PalmarChicoPstI_c5_maf01_Fneg20_mgSNP10_sCov90_94taxa_chr+.hmp.txt",
//            "-o",      baseDir+     "hapmap/keepTaxa/PalmarChicoPstI_c5_maf01_Fneg20_mgSNP10_sCov90_94taxa_chr+.hmp.txt",
//            "-mnSCov", "0.9",  // Minimum SNP coverage (proportion of Taxa)
////            "-mnTCov", "0.01",
//            "-mnF",    " -0.2",
//            "-mnMAF",  "0.01",
////            "-hLD",
//            "-sC",     "1",     // Start chromosome
//            "-eC",     "10"      // End chromosome
//        };
//        GBSHapMapFiltersPlugin pluginGBSHapMapFilters2 = new GBSHapMapFiltersPlugin();
//        pluginGBSHapMapFilters2.setParameters(GBSHapMapFilters2Args);
//        pluginGBSHapMapFilters2.performFunction(null);

        String[] tbt2vcfArgs = new String[] {
            "-i",    baseDir+"tbt/D0E3PACXX_5.tbt.byte",
            "-o",    baseDir+"hapmap/noFilters/vcf",
            "-m",    baseDir+"topm/PstI_PalmarChico_min5.topm",
            "-mnMAF",  "0.0001",  // must be polymorphic
            "-mnLCov", "0.1", // Minimum locus coverage (proportion of Taxa)
            "-ak", "2", // alleles kept
            "-s", "1",  // Start chromosome
            "-e", "10"  // End chromosome
        };
        tbt2vcfPlugin pluginTBT2vcf = new tbt2vcfPlugin();
        pluginTBT2vcf.setParameters(tbt2vcfArgs);
        pluginTBT2vcf.performFunction(null);
    }

    public static void  testTOPMExpandMaxVariants() {
        String TOPMwVariantsFile = "/usr/local/maizediv/illumina/Zea/DTMA/topm/mergedNAM282Ames_variants2_20111018.topm.bin";
        TagsOnPhysicalMap topm = new TagsOnPhysicalMap(TOPMwVariantsFile, true);
        topm.printRows(1000, true, true);
        topm.expandMaxVariants(8);
        topm.printRows(1000, true, true);
    }

    public static void runTagsToSNPByAlignmentPlugin() {
        String baseDirIBM94PstI = "/usr/local/maizediv/illumina/Zea/PstI/";
        String[] IBM94PstIArgs = new String[] {
            "-i",    baseDirIBM94PstI+"tbtByte/B08AAABXX_1.tbt.byte",
            "-y", // use TagsByTaxaByte
            "-o",    baseDirIBM94PstI+"hapmap/quantMnF80MinLCov10/20121003",
            "-m",    baseDirIBM94PstI+"topm/B08AAABXX_1_IBM94PstI_min50.topm.bin",
//            "-mUpd", baseDir+"",
            "-mnF",   "0.8",  // allow some more hets in IBM (tried 0.8 initially, then 0.6, 0.7, and 0.8 again)
            "-mnMAF", "0.2",
            "-mnMAC", "99999",  // this will never be satified: this way -mnMAF overrides it
            "-mnLCov","0.1", // Minimum locus coverage (proportion of Taxa)
//            "-inclGaps",  // Include sites where major or minor allele is a GAP
//            "-callBiSNPsWGap",  //call sites with a biallelic SNP plus a gap (e.g., A/C/-)
            "-s", "1",  // Start chromosome
            "-e", "10"  // End chromosome
        };

        String baseDirCassava = "N:/cassava/";
        String[] CassavaArgs = new String[] {
            "-i",    baseDirCassava+"D09WGACXX_7.tbt.byte",
            "-y", // use TagsByTaxaByte
            "-o",    baseDirCassava+"hapmap/unfilt",
            "-m",    baseDirCassava+"cassava.topm.bin",
//            "-mUpd", baseDir+"",
//            "-mnF",   "-0.5",  // cassava is outbred
            "-mnMAF", "0.05",
            "-mnMAC", "9999",  // this will never be satified: this way -mnMAF overrides it
            "-mnLCov","0.40", // Minimum locus coverage (proportion of Taxa)
//            "-inclGaps",  // Include sites where major or minor allele is a GAP
            "-s", "1",  // Start chromosome
            "-e", "13000"  // End chromosome
        };
        
        String baseDirTestFuzzyPozPstI = "/usr/local/maizediv/illumina/Zea/PstI/";
        String[] testFuzzyPozPstIArgs = new String[] {
            "-i", baseDirTestFuzzyPozPstI+"tbtByte/B08AAABXX_1.tbt.byte",
            "-y", // use TagsByTaxaByte
            "-o", baseDirTestFuzzyPozPstI+"hapmap/testFuzzyPoz",
            "-m", baseDirTestFuzzyPozPstI+"topm/B08AAABXX_1_IBM94PstI_min50.topm.bin",
//            "-p", baseDirTestFuzzyPoz+"AllZeaPedigree20120703.txt",
//            "-mUpd", baseDir+"",
            "-mnF",   "0.8",  // allow some more hets in IBM (tried 0.8 initially, then 0.6, 0.7, and 0.8 again)
            "-mnMAF", "0.2",
            "-mnMAC", "99999",  // this will never be satified: this way -mnMAF overrides it
            "-mnLCov","0.1", // Minimum locus coverage (proportion of Taxa)
//            "-inclGaps",  // Include sites where major or minor allele is a GAP
//            "-callBiSNPsWGap",  //call sites with a biallelic SNP plus a gap (e.g., A/C/-)
            "-s", "10",  // Start chromosome
            "-e", "10"  // End chromosome
        };

        String baseDirTestFuzzyPozApeKI = "/usr/local/maizediv/illumina/Zea/IBM94ApeKI/";
        String[] testFuzzyPozApeKIArgs = new String[] {
            "-i", baseDirTestFuzzyPozApeKI+"tbt/C05F2ACXX_5.tbt.byte",
            "-y", // use TagsByTaxaByte
            "-o", baseDirTestFuzzyPozApeKI+"hapmap/testFuzzyPoz/tol60",
            "-m", baseDirTestFuzzyPozApeKI+"topm/C05F2ACXX_5_IBM94ApeKI_min5.topm.bin",
//            "-mUpd", baseDir+"",
            "-mnF",   "0.8",  // allow some more hets in IBM
            "-mnMAF", "0.2",
            "-mnMAC", "99999",  // this will never be satified: this way -mnMAF overrides it
            "-mnLCov","0.1", // Minimum locus coverage (proportion of Taxa)
//            "-inclGaps",  // Include sites where major or minor allele is a GAP
//            "-callBiSNPsWGap",  //call sites with a biallelic SNP plus a gap (e.g., A/C/-)
            "-s", "10",  // Start chromosome
            "-e", "10"  // End chromosome
        };

        String baseDirTestFuzzyPoz282 = "/usr/local/maizediv/illumina/Zea/MGP1_low_vol/";
        String[] testFuzzyPoz282Args = new String[] {
            "-i", baseDirTestFuzzyPoz282+"tbt/C08L7ACXX_6.tbt.byte",
            "-y", // use TagsByTaxaByte
            "-o", baseDirTestFuzzyPoz282+"hapmap/unfilt/tol150",
            "-m", baseDirTestFuzzyPoz282+"topm/MGP1_low_vol_min3_wPosit.topm.bin",
//            "-mUpd", baseDir+"",
            "-ref", "/usr/local/maizediv/genome/maize_agp_v2.fasta",
            "-LocusBorder", "150",
            "-mnF",   "0.8",  
            "-mnMAF", "0.2",
            "-mnMAC", "99999",  // this will never be satified: this way -mnMAF overrides it
            "-mnLCov","0.1", // Minimum locus coverage (proportion of Taxa)
//            "-inclGaps",  // Include sites where major or minor allele is a GAP
//            "-callBiSNPsWGap",  //call sites with a biallelic SNP plus a gap (e.g., A/C/-)
            "-s", "10",  // Start chromosome
            "-e", "10"  // End chromosome
        };
        
        String baseDirZeaBuild20120701 = "N:/Zea/build20120701/";
        String[] ZeaBuild20120701Args = new String[] {
            "-i", baseDirZeaBuild20120701+"pivotTBT/mergedTBTHDF5_mergedtaxa_pivot_20120628.h5",
//            "-y", // use TagsByTaxaByte
            "-o", baseDirZeaBuild20120701+"hapmap/unfilt",
            "-m", baseDirZeaBuild20120701+"topm/AllZeaMasterTags_c10_20120703.topm",
            "-p", baseDirZeaBuild20120701+"key_files/AllZeaPedigree20120703.txt",
//            "-mUpd", baseDir+"",
            "-mnF",   "0.8",  
            "-mnMAF", "0.001",
            "-mnMAC", "10",  
            "-mnLCov","0.1", // Minimum locus coverage (proportion of Taxa)
//            "-inclGaps",  // Include sites where major or minor allele is a GAP
//            "-callBiSNPsWGap",  //call sites with a biallelic SNP plus a gap (e.g., A/C/-)
            "-s", "10",  // Start chromosome
            "-e", "10"  // End chromosome
        };
        
        String baseDirPearlMillet = "/usr/local/maizediv/illumina/JasonW_PM_GBS/";
        String[] PearlMilletArgs = new String[] {
            "-i", baseDirPearlMillet+"5_MergedTBT.tbt.byte",
            "-y", // use TagsByTaxaByte
            "-o", baseDirPearlMillet+"hapmap/unfilt",
            "-m", baseDirPearlMillet+"3_AlignedTags.topm.bin",
//            "-p", baseDirPearlMillet+"key_files/AllZeaPedigree20120703.txt",
//            "-mUpd", baseDir+"",
            "-mnF",   " -0.2",
            "-mnMAF", "0.01",
            "-mnMAC", "10",
            "-mnLCov","0.1", // Minimum locus coverage (proportion of Taxa)
//            "-inclGaps",  // Include sites where major or minor allele is a GAP
//            "-callBiSNPsWGap",  //call sites with a biallelic SNP plus a gap (e.g., A/C/-)
            "-s", "1",  // Start chromosome
            "-e", "503"  // End chromosome
        };

        String baseDirWillow = "/usr/local/maizediv/willow/";
        //  -TagsToSNPByAlignmentPlugin -i mergedtbt/willow_ApeKI.tbt.byte -m topm/willow_ApeKI.topm.bin -o hapmap -y -s 1 -e 1 -mnMAF 0.01 -mnMAC 10 -mnLCov 0.1 -mxSites 1000000 -endPlugin
        String[] WillowArgs = new String[] {
            "-i", baseDirWillow+"willow_ApeKI.tbt.byte",
            "-y", // use TagsByTaxaByte
            "-o", baseDirWillow+"hapmap/unfilt",
            "-m", baseDirWillow+"willow_ApeKI.topm.bin",
            "-mnMAF", "0.01",
            "-mnMAC", "10",
            "-mnLCov","0.1", // Minimum locus coverage (proportion of Taxa)
            "-mxSites", "300000",
            "-errRate", "0.0033",
            "-s", "1",  // Start chromosome
            "-e", "1"  // End chromosome
        };

        String baseDir282LowVolMin2 = "/Users/jcg233/Documents/GBS/MDP1_low_vol/";
        String[] test282LowVolMin2Args = new String[]{
            "-i", baseDir282LowVolMin2 + "C08L7ACXX_6_min2.tbt.byte",
            "-y", // use TagsByTaxaByte
            "-o", baseDir282LowVolMin2 + "hapmap/tassel3",
            "-m", baseDir282LowVolMin2 + "MGP1_low_vol_min2_wPosit.topm.bin",
            //"-mUpd", baseDir+"",
            //"-ref", "maize_agp_v2.fasta",
            //"-LocusBorder", "150",
            "-mnF", "0.8",
            "-mnMAF", "0.005",
            "-mnMAC", "99999", // this will never be satified: this way -mnMAF overrides it
            "-mnLCov", "0.1", // Minimum locus coverage (proportion of Taxa)
            //            "-inclGaps",  // Include sites where major or minor allele is a GAP
            //            "-callBiSNPsWGap",  //call sites with a biallelic SNP plus a gap (e.g., A/C/-)
            "-s", "10", // Start chromosome
            "-e", "10" // End chromosome
        };

        TagsToSNPByAlignmentPlugin plugin = new TagsToSNPByAlignmentPlugin();
        plugin.setParameters(test282LowVolMin2Args);
        plugin.performFunction(null);
    }

    public static void runQuantPipeline() {  // this method is not finished
        String[] NAM49_50_ApeKIargs = new String[] {
            "-i", "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/qseq",
            "-k", "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/NAM49_50_ApeKI_key.txt",
            "-e", "ApeKI", // Enzyme used to create the GBS library
            "-o", "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/tbt",
            "-c", "1", // Minimum tag count (default is 1).
            "-t", "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/mergedTagCounts/NAM49_50_ApeKI_mergedTags_min10.cnt", // master Tags file
        };
        String[] args = NAM49_50_ApeKIargs;
        QseqToTBTPlugin plugin = new QseqToTBTPlugin();
        plugin.setParameters(args);
        plugin.performFunction(null);
    }

     public static void filterHapMapForTaxa() {
        String baseDir          = "/usr/local/maizediv/illumina/Zea/build20120110/";
        String taxaListFileName = baseDir+"ZeaBuild20120110_maize282_fullNames.txt";
        String infileS =          baseDir+"bpec/Zea20120110_scv10mF8maf002_mgs_E1pLD5kpUn.c+.hmp.txt";
        String outfileS =         baseDir+"bpec/maize282/maize282_Zea20120110_scv10mF8maf002_mgs_E1pLD5kpUn.c+.hmp.txt";
        String infile, outfile;
        TasselPipeline tp = null;
        String[] args;
        int sC=1;
        int eC=10;
        for (int chr=sC; chr<=eC; ++chr) {
            infile=infileS.replace("+", ""+chr);
            outfile=outfileS.replace("+", ""+chr);
            args = new String[] {
                "-fork1",
                "-h", infile,
                "-includeTaxaInFile", taxaListFileName,
                "-taxaJoinStrict", "true",
                "-export", outfile,
                "-exportType", "Hapmap",
                "-runfork1",
            };
            tp = new TasselPipeline(args, null);
        }
    }

    public static void analyzeNIL28_Z026() {
        String baseDir = "/usr/local/maizediv/illumina/Zea/NIL28/";

        String[] FastqToTBTArgs = new String[] {
            "-i", baseDir+"fastq",
            "-k", baseDir+"NIL28_Z026_barcode_key_newTBT.txt",
            "-e", "ApeKI", // Enzyme used to create the GBS library
            "-o", baseDir+"tbtByte",
            "-c", "1", // Minimum tag count (default is 1).
            "-y", // use TagsByTaxaByte
            "-t", "/usr/local/maizediv/illumina/NAM_Ames_282/mergedNAM282Ames.cnt", // master Tags file (only one lane)
        };
        FastqToTBTPlugin pluginFastqToTBT = new FastqToTBTPlugin();
        pluginFastqToTBT.setParameters(FastqToTBTArgs);
        pluginFastqToTBT.performFunction(null);

        String[] MergeTagsByTaxaArgs = new String[] {
            "-i", baseDir+"tbtByte",
            "-o", baseDir+"mergedTBT/mergedNAM282AmesNIL28.tbt.byte",
        };
        String[] args = MergeTagsByTaxaArgs;
        MergeTagsByTaxaFilesPlugin pluginMergeTags = new MergeTagsByTaxaFilesPlugin();
        pluginMergeTags.setParameters(args);
        pluginMergeTags.performFunction(null);

        String[] TagsToSNPByAlignArgs = new String[] {
            "-i",    baseDir+"mergedTBT/mergedNAM282AmesNIL28.tbt.byte",
            "-y",    // use TagsByTaxaByte
            "-o",    baseDir+"hapmap",
            "-m",    "/usr/local/maizediv/illumina/NAM_Ames_282/mergedNAM282Ames_072011.topm.bin",
//            "-mUpd", baseDir+"",
            "-mnF",   "0.7",  // allow some more hets because many of the RCNILs are heterozygous (but they are a small minority prior to taxa filtration)
            "-mnMAF", "0.01",  // will filter for segregation after filter of the hapMap to the taxa of interest
            "-mnMAC", "99999",  // this will never be satified: this way -mnMAF overrides it
            "-mnLCov","0.1", // Minimum locus coverage (proportion of Taxa)  // note that for 384 plex in the future, this should be lower
//            "-inclGaps",  // Include sites where major or minor allele is a GAP
//            "-callBiSNPsWGap",  //call sites with a biallelic SNP plus a gap (e.g., A/C/-)
            "-s", "3",  // Start chromosome
            "-e", "3"  // End chromosome
        };
        TagsToSNPByAlignmentPlugin pluginTagsToSNPByAlign = new TagsToSNPByAlignmentPlugin();
        pluginTagsToSNPByAlign.setParameters(TagsToSNPByAlignArgs);
        pluginTagsToSNPByAlign.performFunction(null);


        // getHapMapReport() was run at this point


        // OOPS - forgot to mergeDuplicateSNPs!


        // filter the chr3 hapmap file for taxa of interest
        String taxaListFileName = baseDir+"NIL28NAM282Ames.c3.hmp.KeepTaxa.txt";
        String infile =           baseDir+"hapmap/NIL28NAM282Ames.c3.hmp.txt";
        String outfile =          baseDir+"hapmap/NIL28rcNILsZ026.c3.hmp.txt";
        String[] filterTaxaArgs = new String[] {
            "-fork1",
            "-h", infile,
            "-includeTaxaInFile", taxaListFileName,
            "-taxaJoinStrict", "true",
            "-export", outfile,
            "-exportType", "Hapmap",
            "-runfork1",
        };
        TasselPipeline tp = new TasselPipeline(filterTaxaArgs, null);

        String[] MergeIdenticalTaxaArgs = new String[] {
            "-hmp", baseDir+"hapmap/NIL28rcNILsZ026.c+.hmp.txt",
            "-o",   baseDir+"hapmap/NIL28rcNILsZ026_mgT76.c+.hmp.txt",
//            "-xHets",
            "-hetFreq", "0.76",
            "-sC", "3",  // Start chromosome
            "-eC", "3" // End chromosome
        };
        MergeIdenticalTaxaPlugin mitp = new MergeIdenticalTaxaPlugin();
        mitp.setParameters(MergeIdenticalTaxaArgs);
        mitp.performFunction(null);

        String[] GBSHapMapFiltersArgs = new String[] {
            "-hmp",    baseDir+"hapmap/NIL28rcNILsZ026_mgT76.c+.hmp.txt",
            "-o",      baseDir+"hapmap/NIL28rcNILsZ026_mgT76_t1_f50_maf10_ld_s5.c+.hmp.txt",
            "-mnTCov", "0.01",
            "-mnF",    "0.50",
            "-mnMAF",  "0.10",
            "-hLD",
            "-mnSCov", "0.05",  // Minimum SNP coverage (proportion of Taxa)
            "-sC",     "3",     // Start chromosome
            "-eC",     "3"      // End chromosome
        };
        GBSHapMapFiltersPlugin pluginGBSHapMapFilters = new GBSHapMapFiltersPlugin();
        pluginGBSHapMapFilters.setParameters(GBSHapMapFiltersArgs);
        pluginGBSHapMapFilters.performFunction(null);

        // remove DCOs, call hetSegs and Impute
        String GBSHapMapFileStem = baseDir+"hapmap/NIL28rcNILs_mgT76_t1_f50_maf10_ld_s5_B73hom_1SNPperTag_";
        Alignment a = ImportUtils.readFromHapmap(GBSHapMapFileStem + "c3.hmp.txt");
        System.out.println("Alignment read:  Taxa:" + a.getSequenceCount() + " Sites:" + a.getSiteCount());
        MutableAlignmentForGBS mutAlign = new MutableAlignmentForGBS(a);
        a = null;
        mutAlign.removeDCOs(7);
        System.out.println("nDCOs on chr3" + " = " + mutAlign.getnDCOs());
        mutAlign.writeToHapmap(false, false, GBSHapMapFileStem + "noDCOs_c3", '\t', 3);
        mutAlign.callHetSegments(7);
        System.out.println("nDCOs on chr3" + " = " + mutAlign.getnDCOs());
        mutAlign.writeToHapmap(false, true, GBSHapMapFileStem + "hetSegs_c3", '\t', 3);
        a = ImportUtils.readFromHapmap(GBSHapMapFileStem + "hetSegs_c3.hmp.txt");
        mutAlign = new MutableAlignmentForGBS(a);
        mutAlign.imputeMissingDataIncludingHets();
        mutAlign.writeToHapmap(false, false, GBSHapMapFileStem + "imputed_c3", '\t', 3);

        // did some manual clean up of small recombinant segments in the QTL2 interval (121 to 139 Mb) by comparing noDCOs vs imputed
        //      --> NIL28rcNILs_mgT76_t1_f50_maf10_ld_s5_B73hom_1SNPperTag_noDCOsCleaned_c3.hmp.txt
        // The below code re-does the het calling and imputation (base on the manually cleaned up noDCOs file) (removeDCOs is also repeated)
        String GBSHapMapFileStem1 = baseDir+"hapmap/NIL28rcNILs_mgT76_t1_f50_maf10_ld_s5_B73hom_1SNPperTag_noDCOsCleaned_";
        Alignment a1 = ImportUtils.readFromHapmap(GBSHapMapFileStem1 + "c3.hmp.txt");
        System.out.println("Alignment read:  Taxa:" + a1.getSequenceCount() + " Sites:" + a1.getSiteCount());
        MutableAlignmentForGBS mutAlign1 = new MutableAlignmentForGBS(a1);
        a1 = null;
        mutAlign1.removeDCOs(7);
        System.out.println("nDCOs on chr3" + " = " + mutAlign1.getnDCOs());
        mutAlign1.writeToHapmap(false, false, GBSHapMapFileStem1 + "noDCOs_c3", '\t', 3);
        mutAlign1.callHetSegments(7);
        System.out.println("nDCOs on chr3" + " = " + mutAlign1.getnDCOs());
        mutAlign1.writeToHapmap(false, true, GBSHapMapFileStem1 + "hetSegs_c3", '\t', 3);
        a1 = ImportUtils.readFromHapmap(GBSHapMapFileStem1 + "hetSegs_c3.hmp.txt");
        mutAlign1 = new MutableAlignmentForGBS(a1);
        mutAlign1.imputeMissingDataIncludingHets();
        mutAlign1.writeToHapmap(false, false, GBSHapMapFileStem1 + "imputed_c3", '\t', 3);

        // did some manual clean up of small recombinant segments in the QTL1 interval (9 to 30 Mb) by comparing noDCOs vs imputed
        //      --> NIL28rcNILs_mgT76_t1_f50_maf10_ld_s5_B73hom_1SNPperTag_noDCOsCleaned2_c3.hmp.txt
        // The below code re-does the het calling and imputation (base on the manually cleaned up noDCOs file) (removeDCOs is also repeated)
        String GBSHapMapFileStem2 = baseDir+"hapmap/NIL28rcNILs_mgT76_t1_f50_maf10_ld_s5_B73hom_1SNPperTag_noDCOsCleaned2_";
        Alignment a2 = ImportUtils.readFromHapmap(GBSHapMapFileStem2 + "c3.hmp.txt");
        System.out.println("Alignment read:  Taxa:" + a2.getSequenceCount() + " Sites:" + a2.getSiteCount());
        MutableAlignmentForGBS mutAlign2 = new MutableAlignmentForGBS(a2);
        a2 = null;
        mutAlign2.removeDCOs(7);
        System.out.println("nDCOs on chr3" + " = " + mutAlign2.getnDCOs());
        mutAlign2.writeToHapmap(false, false, GBSHapMapFileStem2 + "noDCOs_c3", '\t', 3);
        mutAlign2.callHetSegments(7);
        System.out.println("nDCOs on chr3" + " = " + mutAlign2.getnDCOs());
        mutAlign2.writeToHapmap(false, true, GBSHapMapFileStem2 + "hetSegs_c3", '\t', 3);
        a2 = ImportUtils.readFromHapmap(GBSHapMapFileStem2 + "hetSegs_c3.hmp.txt");
        mutAlign2 = new MutableAlignmentForGBS(a2);
        mutAlign2.imputeMissingDataIncludingHets();
        mutAlign2.writeToHapmap(false, false, GBSHapMapFileStem2 + "imputed_c3", '\t', 3);

        // did some manual clean up of small recombinant segments in the large interval between QTL1 & QTL2 interval (30 to 120 Mb) by comparing noDCOs vs imputed
        //      --> NIL28rcNILs_mgT76_t1_f50_maf10_ld_s5_B73hom_1SNPperTag_noDCOsCleaned3_c3.hmp.txt
        //      --> NIL28rcNILs_mgT76_t1_f50_maf10_ld_s5_B73hom_1SNPperTag_noDCOsCleaned4_c3.hmp.txt  (2nd round in same interval)
        // The below code re-does the het calling and imputation (base on the manually cleaned up noDCOs file) (removeDCOs is also repeated)
        String GBSHapMapFileStem3 = baseDir+"hapmap/NIL28rcNILs_mgT76_t1_f50_maf10_ld_s5_B73hom_1SNPperTag_noDCOsCleaned4_";
        Alignment a3 = ImportUtils.readFromHapmap(GBSHapMapFileStem3 + "c3.hmp.txt");
        System.out.println("Alignment read:  Taxa:" + a3.getSequenceCount() + " Sites:" + a3.getSiteCount());
        MutableAlignmentForGBS mutAlign3 = new MutableAlignmentForGBS(a3);
        a3 = null;
        mutAlign3.removeDCOs(7);
        System.out.println("nDCOs on chr3" + " = " + mutAlign3.getnDCOs());
        mutAlign3.writeToHapmap(false, false, GBSHapMapFileStem3 + "noDCOs_c3", '\t', 3);
        mutAlign3.callHetSegments(7);
        System.out.println("nDCOs on chr3" + " = " + mutAlign3.getnDCOs());
        mutAlign3.writeToHapmap(false, true, GBSHapMapFileStem3 + "hetSegs_c3", '\t', 3);
        a3 = ImportUtils.readFromHapmap(GBSHapMapFileStem3 + "hetSegs_c3.hmp.txt");
        mutAlign3 = new MutableAlignmentForGBS(a3);
        mutAlign3.imputeMissingDataIncludingHets();
        mutAlign3.writeToHapmap(false, false, GBSHapMapFileStem3 + "imputed_c3", '\t', 3);
    }

    public static void analyzeTeoFineMapping2012() {
        String baseDir = "/usr/local/maizediv/illumina/teosinte/FineMapping2012/";

        String[] FastqToTagCountArgs = new String[] {
            "-i", baseDir + "fastq",
            "-k", baseDir + "teoFineMapping_barcode_key20120410fastq.txt",
            "-e", "ApeKI", // Enzyme used to create the GBS library
//            "-s", "20000000", // Max good reads per lane. (Optional. Default is 200,000,000)
//            "-c", "1", // Minimum tag count (default is 1)
            "-o", baseDir + "tagCounts",
        };
        FastqToTagCountPlugin pluginFastqToTagCount = new FastqToTagCountPlugin();
        pluginFastqToTagCount.setParameters(FastqToTagCountArgs);
        pluginFastqToTagCount.performFunction(null);

        String[] QseqToTagCountArgs = new String[] {
            "-i", baseDir + "qseq",
            "-k", baseDir + "teoFineMapping_barcode_key20120410qseq.txt",
            "-e", "ApeKI", // Enzyme used to create the GBS library
//            "-s", "20000000", // Max good reads per lane. (Optional. Default is 200,000,000)
//            "-c", "1", // Minimum tag count (default is 1)
            "-o", baseDir + "tagCounts",
        };
        QseqToTagCountPlugin pluginQseqToTagCount = new QseqToTagCountPlugin();
        pluginQseqToTagCount.setParameters(QseqToTagCountArgs);
        pluginQseqToTagCount.performFunction(null);

        String[] MergeMultipleTagCountArgs = new String[] {
            "-i", baseDir+"tagCounts", // Input directory containing .cnt files
            "-o", baseDir+"mergedTagCounts/TeoFineMapping2012_mergedTagCounts_min3", //  Output file name (*.fq got appended)
            "-c", "3", // Minimum count of reads to be output (default 1)
            "-t", // Specifies that reads should be output in FASTQ text format.
        };
        MergeMultipleTagCountPlugin pluginMergeMultipleTagCount = new MergeMultipleTagCountPlugin();
        pluginMergeMultipleTagCount.setParameters(MergeMultipleTagCountArgs);
        pluginMergeMultipleTagCount.performFunction(null);

        // BWA commands:
        // bwa aln -t 8 /usr/local/maizediv/genome/maize_agp_v2.fasta  TeoFineMapping2012_mergedTagCounts_min3.fq > TeoFineMapping2012_mergedTagCounts_min3.sai
        // bwa samse /usr/local/maizediv/genome/maize_agp_v2.fasta  TeoFineMapping2012_mergedTagCounts_min3.sai  TeoFineMapping2012_mergedTagCounts_min3.fq > TeoFineMapping2012_mergedTagCounts_min3.sam

        String[] SAMConverterPluginArgs = new String[] {
            "-i", baseDir+"mergedTagCounts/TeoFineMapping2012_mergedTagCounts_min3.sam",
            "-o", baseDir+           "topm/TeoFineMapping2012_min3.topm.bin",
        };
        SAMConverterPlugin pluginSAMConverter = new SAMConverterPlugin();
        pluginSAMConverter.setParameters(SAMConverterPluginArgs);
        pluginSAMConverter.performFunction(null);

        TagsOnPhysicalMap topm = new TagsOnPhysicalMap(baseDir+"topm/TeoFineMapping2012_min3.topm.bin", true);
        topm.writeBinaryFile(new File(                 baseDir+"topm/TeoFineMapping2012_min3_wPosit.topm.bin"), Integer.MAX_VALUE, true, false, Float.MIN_VALUE, true);

        String[] FastqToTBTArgs = new String[] {
            "-i", baseDir+"fastq",
            "-k", baseDir+"teoFineMapping_barcode_key20120410fastq.txt",
            "-e", "ApeKI", // Enzyme used to create the GBS library
            "-o", baseDir+"tbtByte",
//            "-s", "200000000", // Max good reads per lane. (Optional. Default is 200,000,000)
//            "-c", "1", // Minimum tag count (default is 1).
            "-y", // use TagsByTaxaByte
            "-m", baseDir+"topm/TeoFineMapping2012_min3_wPosit.topm.bin",  // use tags with unique physical postions for smaller TBT
        };
        FastqToTBTPlugin pluginFastqToTBT = new FastqToTBTPlugin();
        pluginFastqToTBT.setParameters(FastqToTBTArgs);
        pluginFastqToTBT.performFunction(null);

        String[] QseqToTBTArgs = new String[] {
            "-i", baseDir+"qseq",
            "-k", baseDir+"teoFineMapping_barcode_key20120410qseq.txt",
            "-e", "ApeKI", // Enzyme used to create the GBS library
            "-o", baseDir+"tbtByte",
//            "-c", "1", // Minimum tag count (default is 1).
            "-y", // use TagsByTaxaByte
            "-m", baseDir+"topm/TeoFineMapping2012_min3_wPosit.topm.bin",  // use tags with unique physical postions for smaller TBT
        };
        QseqToTBTPlugin pluginQseqToTBT = new QseqToTBTPlugin();
        pluginQseqToTBT.setParameters(QseqToTBTArgs);
        pluginQseqToTBT.performFunction(null);

        String[] MergeTagsByTaxaArgs = new String[] {
            "-i", baseDir+"tbtByte",
            "-o", baseDir+"mergedTBT/TeoFineMapping2012_min3_wPosit.tbt.byte",
        };
        MergeTagsByTaxaFilesPlugin pluginMergeTags = new MergeTagsByTaxaFilesPlugin();
        pluginMergeTags.setParameters(MergeTagsByTaxaArgs);
        pluginMergeTags.performFunction(null);

        String[] TagsToSNPByAlignmentArgs = new String[] {
            "-i",    baseDir+"mergedTBT/TeoFineMapping2012_min3_wPosit.tbt.byte",
            "-y",    // use TagsByTaxaByte
            "-o",    baseDir+"hapmap/raw",
            "-m",    baseDir+"topm/TeoFineMapping2012_min3_wPosit.topm.bin",
//            "-mUpd", baseDir+"",
            "-mnF",   "0.2",  // keep this low b/c I don't know how many hets to expect
            "-mnMAF", "0.02",
            "-mnMAC", "99999",  // this will never be satified: this way -mnMAF overrides it
            "-mnLCov","0.05", // Minimum locus coverage (proportion of Taxa)
//            "-inclGaps",  // Include sites where major or minor allele is a GAP
//            "-callBiSNPsWGap",  //call sites with a biallelic SNP plus a gap (e.g., A/C/-)
            "-s", "1",  // Start chromosome
            "-e", "10"  // End chromosome
        };
        TagsToSNPByAlignmentPlugin pluginTagsToSNPByAlignment = new TagsToSNPByAlignmentPlugin();
        pluginTagsToSNPByAlignment.setParameters(TagsToSNPByAlignmentArgs);
        pluginTagsToSNPByAlignment.performFunction(null);

        String[] mergeDupSNPsArgs = new String[] {
            "-hmp",     baseDir+"hapmap/raw/mergedTBT.c+.hmp.txt",
            "-o",       baseDir+"hapmap/mergedSNPs/TeoFineMapping2012_min3_mgSNP20_chr+.hmp.txt",
            "-misMat",  "0.3", // Threshold genotypic mismatch rate above which the duplicate SNPs won't be merged (default: 0.05)
            "-callHets",       // When two genotypes disagree at a SNP, call it a heterozygote (default: off = set to missing)
//            "-kpUnmergDups",   // When two duplicate SNPs were not merged (different alleles or too many mismatches), keep them (default: off = delete them)
            "-s", "1",        // Start chromosome (default 1)
            "-e", "10"         // End chromosome (default 10)
        };
        MergeDuplicateSNPsPlugin pluginMergeDuplicateSNPs = new MergeDuplicateSNPsPlugin();
        pluginMergeDuplicateSNPs.setParameters(mergeDupSNPsArgs);
        pluginMergeDuplicateSNPs.performFunction(null);
    }

    public static void analyzeZakChr5FineMapping() {
        String baseDir = "/usr/local/maizediv/illumina/teosinte/FineMapping2012/";

        // filter the chr5 hapmap file for taxa of interest
        String taxaListFileName = baseDir+"ZakChr5FineMappingTaxa2012.txt";
        String infile =           baseDir+"hapmap/mergedSNPs/TeoFineMapping2012_min3_mgSNP20_chr5.hmp.txt";
        String outfile =          baseDir+"hapmap/mergedSNPs/Zak_min3_mgSNP20_chr5.hmp.txt";
        String[] filterTaxaArgs = new String[] {
            "-fork1",
            "-h", infile,
            "-includeTaxaInFile", taxaListFileName,
            "-taxaJoinStrict", "true",
            "-export", outfile,
            "-exportType", "Hapmap",
            "-runfork1",
        };
        TasselPipeline tp = new TasselPipeline(filterTaxaArgs, null);

        String[] MergeIdenticalTaxaArgs = new String[] {
            "-hmp", baseDir+"hapmap/mergedSNPs/Zak_min3_mgSNP20_chr+.hmp.txt",
            "-o",   baseDir+"hapmap/mergedTaxa/Zak_min3_mgSNP20_mgT76_chr+.hmp.txt",
//            "-xHets",
            "-hetFreq", "0.76",
            "-sC", "5",  // Start chromosome
            "-eC", "5" // End chromosome
        };
        MergeIdenticalTaxaPlugin mitp = new MergeIdenticalTaxaPlugin();
        mitp.setParameters(MergeIdenticalTaxaArgs);
        mitp.performFunction(null);

        String[] GBSHapMapFiltersArgs = new String[] {
            "-hmp",    baseDir+"hapmap/mergedTaxa/Zak_min3_mgSNP20_mgT76_chr+.hmp.txt",
            "-o",      baseDir+      "hapmap/filt/Zak_min3_mgSNP20_mgT76_s10f70maf5hLDmnRsq90_chr+.hmp.txt",
            "-mnSCov", "0.1",  // Minimum SNP coverage (proportion of Taxa)
//            "-mnTCov", "0.01",
            "-mnF",    "0.70",
            "-mnMAF",  "0.05",
            "-hLD",             // I temporarily set the minR2 to 0.9 and BonP to 0.001 (they are hard coded to 0.01 & 0.01)
            "-sC",     "5",     // Start chromosome
            "-eC",     "5"      // End chromosome
        };
        GBSHapMapFiltersPlugin pluginGBSHapMapFilters = new GBSHapMapFiltersPlugin();
        pluginGBSHapMapFilters.setParameters(GBSHapMapFiltersArgs);
        pluginGBSHapMapFilters.performFunction(null);

        // remove DCOs, call hetSegs and Impute         Zak_min3_mgSNP20_mgT76_s10f70maf5hLDmnRsq90_1SNPperTag_phased_manFilt2_chr5.hmp.txt = round 2
        String GBSHapMapFileStem = baseDir+"hapmap/filt/Zak_min3_mgSNP20_mgT76_s10f70maf5hLDmnRsq90_1SNPperTag_phased_manFilt2_";
        Alignment a = ImportUtils.readFromHapmap(GBSHapMapFileStem + "chr5.hmp.txt");
        System.out.println("Alignment read:  Taxa:" + a.getSequenceCount() + " Sites:" + a.getSiteCount());
        MutableAlignmentForGBS mutAlign = new MutableAlignmentForGBS(a);
        a = null;
        mutAlign.removeDCOs(7);
        System.out.println("nDCOs on chr5" + " = " + mutAlign.getnDCOs());
        mutAlign.writeToHapmap(false, false, GBSHapMapFileStem + "noDCOs_chr5", '\t', 3);
        mutAlign.callHetSegments(7);
        System.out.println("nDCOs on chr5" + " = " + mutAlign.getnDCOs());
        mutAlign.writeToHapmap(false, true, GBSHapMapFileStem + "hetSegs_chr5", '\t', 3);
        a = ImportUtils.readFromHapmap(GBSHapMapFileStem + "hetSegs_chr5.hmp.txt");
        mutAlign = new MutableAlignmentForGBS(a);
        mutAlign.imputeMissingDataIncludingHets();
        mutAlign.writeToHapmap(false, false, GBSHapMapFileStem + "imputed_chr5", '\t', 3);
    }

    public static void analyzeETb1FineMapping() {
        String baseDir = "/usr/local/maizediv/illumina/teosinte/FineMapping2012/";

        String[] MergeIdenticalTaxaArgs = new String[] {
            "-hmp", baseDir+"hapmap/mergedSNPs/eTb1FineMapping2012_min3_mgSNP20_chr+gt220Mb.hmp.txt",
            "-o",   baseDir+"hapmap/mergedTaxa/eTb1FineMapping2012_min3_mgSNP20_mgT76_chr+gt220Mb.hmp.txt",
//            "-xHets",
            "-hetFreq", "0.76",
            "-sC", "1",  // Start chromosome
            "-eC", "1" // End chromosome
        };
        MergeIdenticalTaxaPlugin mitp = new MergeIdenticalTaxaPlugin();
        mitp.setParameters(MergeIdenticalTaxaArgs);
        mitp.performFunction(null);

        String[] GBSHapMapFiltersArgs = new String[] {
            "-hmp",    baseDir+"hapmap/mergedTaxa/eTb1FineMapping2012_min3_mgSNP20_mgT76_chr+gt220Mb.hmp.txt",
            "-o",      baseDir+      "hapmap/filt/eTb1FineMapping2012_min3_mgSNP20_mgT76_s10f50maf10hLDmnRsq80_chr+gt220Mb.hmp.txt",
            "-mnSCov", "0.1",  // Minimum SNP coverage (proportion of Taxa)
//            "-mnTCov", "0.01",
            "-mnF",    "0.50",
            "-mnMAF",  "0.1",
            "-hLD",             // I temporarily set the minR2 to 0.8 and BonP to 0.001 (they are hard coded to 0.01 & 0.01) (I also tested minR2 of 0.1, 0.4, 0.5, 0.6, 0.7)
            "-sC",     "1",     // Start chromosome
            "-eC",     "1"      // End chromosome
        };
        GBSHapMapFiltersPlugin pluginGBSHapMapFilters = new GBSHapMapFiltersPlugin();
        pluginGBSHapMapFilters.setParameters(GBSHapMapFiltersArgs);
        pluginGBSHapMapFilters.performFunction(null);

        // remove DCOs, call hetSegs and Impute         eTb1_min3_mgSNP20_mgT76_s10f50maf10hLDmnRsq80_phased_1SNPperTag_manFilt1_chr1.hmp.txt = round 1
        String GBSHapMapFileStem = baseDir+"hapmap/filt/eTb1_min3_mgSNP20_mgT76_s10f50maf10hLDmnRsq80_phased_1SNPperTag_manFilt1_";
        Alignment a = ImportUtils.readFromHapmap(GBSHapMapFileStem + "chr1.hmp.txt");
        System.out.println("Alignment read:  Taxa:" + a.getSequenceCount() + " Sites:" + a.getSiteCount());
        MutableAlignmentForGBS mutAlign = new MutableAlignmentForGBS(a);
        a = null;
        mutAlign.removeDCOs(7);
        System.out.println("nDCOs on chr5" + " = " + mutAlign.getnDCOs());
        mutAlign.writeToHapmap(false, false, GBSHapMapFileStem + "noDCOs_chr1", '\t', 3);
        mutAlign.callHetSegments(7);
        System.out.println("nDCOs on chr5" + " = " + mutAlign.getnDCOs());
        mutAlign.writeToHapmap(false, true, GBSHapMapFileStem + "hetSegs_chr1", '\t', 3);
        a = ImportUtils.readFromHapmap(GBSHapMapFileStem + "hetSegs_chr1.hmp.txt");
        mutAlign = new MutableAlignmentForGBS(a);
        mutAlign.imputeMissingDataIncludingHets();
        mutAlign.writeToHapmap(false, false, GBSHapMapFileStem + "imputed_chr1", '\t', 3);
    }

     public static void testTBTHDF5TagGroups() {
        String baseDir = "/usr/local/maizediv/illumina/Zea/build20120110/tbt/";
        String inTBTFile = baseDir+"C08FFACXX_5.tbt.byte";
        String outTBTHDF5TagGroupsFile = baseDir+"C08FFACXX_5_TagGroups.tbt.h5";
        String outTBTByteFileTagGroupsTest = baseDir+"C08FFACXX_5_TagGroupsTest.tbt.byte";

        TagsByTaxaByte inTBT=new TagsByTaxaByte(inTBTFile, FilePacking.Byte);
        TagsByTaxaByteHDF5TagGroups myTBTHDF5TagGroups = new TagsByTaxaByteHDF5TagGroups(inTBT,outTBTHDF5TagGroupsFile);
        myTBTHDF5TagGroups.getFileReadyForClosing();
        myTBTHDF5TagGroups = null;
        inTBT = null;
        System.gc();
        myTBTHDF5TagGroups = new TagsByTaxaByteHDF5TagGroups(outTBTHDF5TagGroupsFile);
        TagsByTaxaByte outTBT = myTBTHDF5TagGroups.convertToTBTByte();
        myTBTHDF5TagGroups.getFileReadyForClosing();
        myTBTHDF5TagGroups = null;
        outTBT.writeDistFile(new File(outTBTByteFileTagGroupsTest), FilePacking.Byte, 0);
    }

     public static void testTBTHDF5TaxaGroups() {
        String baseDir = "/usr/local/maizediv/illumina/Zea/build20120110/tbt/";
        String inTBTFile =                    baseDir+"C08FFACXX_5.tbt.byte";
        String outTBTHDF5TaxaGroupsFile =     baseDir+"C08FFACXX_5_TaxaGroups.tbt.h5";
        String outTBTByteFileTaxaGroupsTest = baseDir+"C08FFACXX_5_TaxaGroupsTest.tbt.byte";

        TagsByTaxaByte inTBT=new TagsByTaxaByte(inTBTFile, FilePacking.Byte);
        TagsByTaxaByteHDF5TaxaGroups myTBTHDF5TaxaGroups = new TagsByTaxaByteHDF5TaxaGroups(inTBT, outTBTHDF5TaxaGroupsFile);
        myTBTHDF5TaxaGroups.getFileReadyForClosing();
        myTBTHDF5TaxaGroups = null;
        inTBT = null;
        System.gc();
        myTBTHDF5TaxaGroups = new TagsByTaxaByteHDF5TaxaGroups(outTBTHDF5TaxaGroupsFile);
        TagsByTaxaByte outTBT = myTBTHDF5TaxaGroups.convertToTBTByte();
        myTBTHDF5TaxaGroups.getFileReadyForClosing();
        myTBTHDF5TaxaGroups = null;
        outTBT.writeDistFile(new File(outTBTByteFileTaxaGroupsTest), FilePacking.Byte, 0);

        // since the file size was correct from the above but the md5sum wasn't, try another round trip to see if the taxon order gets in sync
        String outTBTHDF5TaxaGroupsFile2 =     baseDir+"C08FFACXX_5_TaxaGroups2.tbt.h5";
        String outTBTByteFileTaxaGroupsTest2 = baseDir+"C08FFACXX_5_TaxaGroupsTest2.tbt.byte";

        TagsByTaxaByte inTBT2 = new TagsByTaxaByte(outTBTByteFileTaxaGroupsTest , FilePacking.Byte);
        TagsByTaxaByteHDF5TaxaGroups myTBTHDF5TaxaGroups2 = new TagsByTaxaByteHDF5TaxaGroups(inTBT2, outTBTHDF5TaxaGroupsFile2);
        myTBTHDF5TaxaGroups2.getFileReadyForClosing();
        myTBTHDF5TaxaGroups2 = null;
        inTBT2 = null;
        System.gc();
        myTBTHDF5TaxaGroups2 = new TagsByTaxaByteHDF5TaxaGroups(outTBTHDF5TaxaGroupsFile2);
        TagsByTaxaByte outTBT2 = myTBTHDF5TaxaGroups2.convertToTBTByte();
        myTBTHDF5TaxaGroups2.getFileReadyForClosing();
        myTBTHDF5TaxaGroups2 = null;
        outTBT2.writeDistFile(new File(outTBTByteFileTaxaGroupsTest2), FilePacking.Byte, 0);

        String[] BinToTextArgs = new String[] {
            "-i", outTBTByteFileTaxaGroupsTest,
            "-o", baseDir+"C08FFACXX_5_TaxaGroupsTest.tbt.txt",
            "-t", "TBTByte",  // TOPM, TagCounts, TBTBit, TBTByte
        };
        BinaryToTextPlugin plugin = new BinaryToTextPlugin(null);
        plugin.setParameters(BinToTextArgs);
        plugin.performFunction(null);

        String[] BinToTextArgsOrig = new String[] {
            "-i", inTBTFile,
            "-o", baseDir+"C08FFACXX_5.tbt.txt",
            "-t", "TBTByte",  // TOPM, TagCounts, TBTBit, TBTByte
        };
        BinaryToTextPlugin plugin2 = new BinaryToTextPlugin(null);
        plugin2.setParameters(BinToTextArgsOrig);
        plugin2.performFunction(null);
    }
     
    public static void runModifyTBTHDF5Plugin() {
        String baseDir = "/usr/local/maizediv/illumina/Zea/build20120701/";
        String[] testMergeTBTHDF5Args = new String[] {
            "-i", baseDir+"tbt/node2_ZM_c10_TBT_20120615.hdf",
            "-o", baseDir+"mergedTBT/mergedTBTHDF5_20120618.h5",
            "-L", baseDir+"mergedTBT/mergedTBTHDF5_20120618.log"
        };
        ModifyTBTHDF5Plugin plugin = new ModifyTBTHDF5Plugin(null);
        plugin.setParameters(testMergeTBTHDF5Args);
        plugin.performFunction(null);
    }
     
    public static void runSAMConverterPlugin() {
        String baseDir = "N:/Zea/build20120701/";
        String[] testSAMConverterArgs = new String[] {
            "-i", baseDir+"sam_file/AllZeaMasterTags_c10_20120613_1st10K.sam",
            "-o", baseDir+"topm/AllZeaMasterTags_c10_20120702_1st10K.topm.txt",
            "-t"  // text output
        };
        SAMConverterPlugin plugin = new SAMConverterPlugin(null);
        plugin.setParameters(testSAMConverterArgs);
        plugin.performFunction(null);
    }
    
    public static void runCompareGenosBetweenHapMapFilesPlugin() {
        String baseDir = "/usr/local/maizediv/illumina/Zea/";
        String[] maize282Jan2012BuildVs50KAGPv2Args = new String[] {
            "-hmp1", baseDir+"build20120110/imp/maize282/282_20120110_scv10mF8maf002_mgs_E1pLD5kpUn_imp95_1024_chr+.hmp.txt",
            "-hmp2", baseDir+"50K/SNP55K_maize282_AGPv2_20100513_chr+.hmp.txt",
            "-syn",  baseDir+"build20120110/imp/maize282/maize282FullNames_GBS20120110PanzeaImputed_to_50K_AGPv2_noW64A.txt",
            "-o",    baseDir+"50K/compareGBS20120110PanzeaImputedMaize282GenosTo50K_noW64A.txt",
            "-sC",    "1",     // Start chromosome
            "-eC",   "10"     // End chromosome
        };

        baseDir = "N:/Zea/build20120701/HMP_For_Testing/";
        String[] maize282July2012BuildVs50KAGPv2Args = new String[] {
            "-hmp1", baseDir+"Bowtie2_BPEC_55KTaxa.c10.hmp.txt",
            "-hmp2", baseDir+"SNP55K_maize282_AGPv2_chr10_20100513.hmp.txt",
            "-syn",  baseDir+"LinkFileForTestAgainst55K_switch_columns.txt",
            "-o",    baseDir+"55KVGBSver2.txt",
            "-sC",   "10",     // Start chromosome
            "-eC",   "10"     // End chromosome
        };

        baseDir = "/Users/jcg233/Documents/GBS/MDP1_low_vol/hapmap/";
        String[] tassel3vs4Args = new String[] {
            "-hmp1", baseDir+"tassel3/MDP1_low_vol.c10.hmp.txt",
            "-hmp2", baseDir+"MDP1_low_vol.c10.hmp.txt",
            "-syn",  baseDir+"MDPLowVolTaxaLinkFile.txt",
            "-o",    baseDir+"MDPLowVolTassel3vs4.txt",
            "-sC",   "10",     // Start chromosome
            "-eC",   "10"     // End chromosome
        };
        CompareGenosBetweenHapMapFilesPlugin plugin = new CompareGenosBetweenHapMapFilesPlugin(null);
        plugin.setParameters(tassel3vs4Args);
        plugin.performFunction(null);
    }
    
    public static void runFastImputationBitFixedWindowPlugin() {
        String[] maize282July2012BuildArgs = new String[] {
            "-hmp", "/Volumes/nextgen/Zea/build20120701/06_HapMap/Bowtie2/04_BPECFilteredSNPs/PivotTBT.mergedSNPs.hapmapfiltered.BPEC.c10.hmp.txt",
            "-o",   "/Users/jcg233/Documents/GBS/Zea/build20120701/hapmap/imputed/PivotTBT.mergedSNPs.hapmapfiltered.BPEC.imp.c10.hmp.txt",
            "-p",   "/Volumes/nextgen/Zea/build20120701/50_KeyFiles/AllZeaPedigree20120730.txt",
        };
        FastImputationBitFixedWindowPlugin plugin = new FastImputationBitFixedWindowPlugin(null);
        plugin.setParameters(maize282July2012BuildArgs);
        plugin.performFunction(null);
    }

    public static void analyzeMGP1_low_vol() {
        String baseDir = "/usr/local/maizediv/illumina/Zea/MGP1_low_vol/";
        
        String[] FastqToTagCountArgs = new String[] {
            "-i", baseDir + "fastq",
            "-k", baseDir + "MGP1_low_vol_key.txt",
            "-e", "ApeKI", // Enzyme used to create the GBS library
//            "-s", "20000000", // Max good reads per lane. (Optional. Default is 200,000,000)
//            "-c", "1", // Minimum tag count (default is 1)
            "-o", baseDir + "tagCounts",
        };
        FastqToTagCountPlugin pluginFastqToTagCount = new FastqToTagCountPlugin();
        pluginFastqToTagCount.setParameters(FastqToTagCountArgs);
        pluginFastqToTagCount.performFunction(null);

        String[] MergeMultipleTagCountArgs = new String[] {
            "-i", baseDir+"tagCounts", // Input directory containing .cnt files
            "-o", baseDir+"sam/MGP1_low_vol_min2", //  Output file name (*.fq gets appended)
            "-c", "2", // Minimum count of reads to be output (default 1)
            "-t", // Specifies that reads should be output in FASTQ text format.
        };
        MergeMultipleTagCountPlugin pluginMergeMultipleTagCount = new MergeMultipleTagCountPlugin();
        pluginMergeMultipleTagCount.setParameters(MergeMultipleTagCountArgs);
        pluginMergeMultipleTagCount.performFunction(null);

        // BWA commands:
        // bwa aln -t 8 /usr/local/maizediv/genome/maize_agp_v2.fasta  MGP1_low_vol_min2.fq > MGP1_low_vol_min2.sai
        // bwa samse /usr/local/maizediv/genome/maize_agp_v2.fasta  MGP1_low_vol_min2.sai  MGP1_low_vol_min2.fq > MGP1_low_vol_min2.sam

        String[] SAMConverterPluginArgs = new String[] {
            "-i", baseDir+ "sam/MGP1_low_vol_min2.sam",
            "-o", baseDir+"topm/MGP1_low_vol_min2.topm.bin",
        };
        SAMConverterPlugin pluginSAMConverter = new SAMConverterPlugin();
        pluginSAMConverter.setParameters(SAMConverterPluginArgs);
        pluginSAMConverter.performFunction(null);

        // filter TOPM for only those tags with unique positions
        TagsOnPhysicalMap topm = new TagsOnPhysicalMap(baseDir+"topm/MGP1_low_vol_min2.topm.bin", true);
        topm.writeBinaryFile(new File(baseDir+"topm/MGP1_low_vol_min2_wPosit.topm.bin"), Integer.MAX_VALUE, true, false, Float.MIN_VALUE, true);

        String[] FastqToTBTArgs = new String[] {
            "-i", baseDir+"fastq",
            "-k", baseDir+"MGP1_low_vol_key.txt",
            "-e", "ApeKI", // Enzyme used to create the GBS library
            "-o", baseDir+"tbt",
//            "-s", "200000000", // Max good reads per lane. (Optional. Default is 200,000,000)
//            "-c", "1", // Minimum tag count (default is 1).
            "-y", // use TagsByTaxaByte
            "-m", baseDir+"topm/MGP1_low_vol_min2_wPosit.topm.bin",  // use tags with unique physical postions for smaller TBT
        };
        FastqToTBTPlugin pluginFastqToTBT = new FastqToTBTPlugin();
        pluginFastqToTBT.setParameters(FastqToTBTArgs);
        pluginFastqToTBT.performFunction(null);
    }

    public static void printChrsFromTOPM() {
        String topmFile = "/usr/local/maizediv/illumina/JasonW_PM_GBS/3_AlignedTags.topm.bin";
        TagsOnPhysicalMap myTOPM = new TagsOnPhysicalMap(topmFile, true);
        int[] chrs = myTOPM.getChromosomes();
        Arrays.sort(chrs);
        System.out.println("The chromosomes in this TOPM are:");
        for (int c=0; c<chrs.length; c++) {
            System.out.println(chrs[c]+"");
        }
    }

    public static void getHapMapReport() {
        String infileName = "/Users/jcg233/Documents/GBS/Zea/PstI_PalmarChico/hapmap/noFilters/filtCovPoly/PalmarChicoPstI_c5_MAC1_94taxa_mgSNPs10_chr10.hmp.txt";

        Alignment align=ImportUtils.readFromHapmap(infileName);

        int taxonCount=align.getIdGroup().getIdCount(), nonBlankTaxaCount=0;
        for (int i = 0; i < taxonCount; i++){
            String name=align.getIdGroup().getIdentifier(i).getName();
                if(!name.equalsIgnoreCase("blank") && !name.equalsIgnoreCase("empty")) nonBlankTaxaCount++;
        }

        int siteCount=align.getSiteCount();

        //Taxon properties
        int[] taxonCalls=new int[taxonCount];
        int[] taxonHetCalls=new int[taxonCount];

        //Histogram properties
        int numBins=20;
        int binWidth=100/numBins;
        int[] bins=new int[numBins];
        int[] hetBins=new int[numBins];

        System.out.println(
            infileName+"\n"+
            "Total taxa:\t"+align.getIdGroup().getIdCount()+"\n"+
            "Total non-\"blank\" taxa:\t"+nonBlankTaxaCount+"\n"+
            "Total SNPs:\t"+align.getSiteCount()+"\n"
        );

        System.out.println();
        System.out.println("Site Properties");
        System.out.println(
                "Site\t"
                + "Position\t"
                + "Calls\t"
                + "HetCalls\t"
                + "CallRate\t"
                + "ObsHet\t"
                + "MAF\t"
                + "Fis"
        );
        
        int[][] hetCnt = AlignmentFilterByGBSUtils.genotypicCountsBySite(align, false, false);

        for (int site = 0; site < siteCount; site++){
            int calls=0, hetCalls=0;
            //Skip uncalled bases.  Increment "calls" for all called bases and "hetCalls" for hets.
            for (int taxon = 0; taxon < taxonCount; taxon++) {
                if(align.getIdGroup().getIdentifier(taxon).getName().equalsIgnoreCase("blank")) continue;    //Skip blank samples

                char base = (char)align.getBase(taxon, site);
                switch(base){
                    case 'N':
                        break;
                    case 'A':  case 'C':  case 'G':  case 'T':
                        calls++;
                        taxonCalls[taxon]++;
                        break;
                    default:
                        calls++;
                        hetCalls++;
                        taxonCalls[taxon]++;
                        taxonHetCalls[taxon]++;
                        break;
                }
            }

            //Calculate call rates (site coverage) as #calls/#taxa & other stats
            double callRate = ((double)calls/(double)nonBlankTaxaCount);
            double obsHet = (calls > 0) ? ((double)hetCalls/(double)calls) : Double.NaN;
            double maf = hetCnt[0][site]>0 ? (double)(hetCnt[3][site]+((double)hetCnt[1][site]/2.0))/(double)hetCnt[0][site] : Double.NaN;
            double fis = (obsHet == Double.NaN) ? Double.NaN : 1.0 - (obsHet/(2.0*maf*(1.0-maf)));
            
            
            //Add SNP to an appropriate chart bin based on its call rate
            double percentile=callRate*100;
            for (int bin = 0; bin < bins.length; bin++) {
                int lowerBound = (bin*binWidth);
                int upperBound = ((bin+1)*binWidth);
                if (bin > 0) {
                    if(percentile <= upperBound && percentile > lowerBound)  bins[bin]+=1;
                } else {
                    if(percentile <= upperBound)  bins[bin]+=1;
                }
            }

            System.out.println(
                site+"\t"+
                align.getPositionInLocus(site)+"\t"+
                calls+"\t"+
                hetCalls+"\t"+
                callRate+"\t"+
                obsHet+"\t"+
                maf+"\t"+
                fis
            );
        }


        //Print bar chart data for sites
        System.out.println();
        System.out.println("Site Coverage:");
        System.out.println("%CalledTaxa:\t"+"nSNPs");
        for (int bin = 0; bin < bins.length; bin++) {
            int callRatePct =((bin+1)*binWidth);
            System.out.println(callRatePct+"\t"+bins[bin]);
            bins[bin]=0; //Re-zero bins after they are printed
        }

        //Print taxon coverage stats
        System.out.println();
        System.out.println("Taxon Properties");
        System.out.println(
                "FullName\t"
                +"Calls\t"
                +"HetCalls\t"
                +"CallRate\t"
                +"ObsHet"
        );

        for (int i = 0; i < taxonCount; i++) {
            if(align.getIdGroup().getIdentifier(i).getName().equalsIgnoreCase("blank")) continue;    //Skip blank samples
            double callRate=(double)taxonCalls[i]/(double)siteCount;
            double obsHet = taxonCalls[i]>0 ? (double)taxonHetCalls[i]/(double)taxonCalls[i] : Double.NaN;
            System.out.println(
                    align.getFullTaxaName(i)+"\t"
                    +taxonCalls[i]+"\t"
                    +taxonHetCalls[i]+"\t"
                    +callRate+"\t"
                    +obsHet
            );

            //Add SNP to an appropriate chart bin based on its call rate
            double percentile=callRate*100;
            for (int bin = 0; bin < bins.length; bin++) {
                int lowerBound = (bin*binWidth);
                int upperBound = ((bin+1)*binWidth);
                if (bin > 0) {
                    if(percentile <= upperBound && percentile > lowerBound)  bins[bin]+=1;
                } else {
                    if(percentile <= upperBound)  bins[bin]+=1;
                }
            }
        }

        //Print bar chart data for taxa
        System.out.println();
        System.out.println("Taxon Coverage:");
        System.out.println("%CalledSites:\t"+"nTaxa");
        for (int bin = 0; bin < bins.length; bin++){
            int callRatePct =((bin+1)*binWidth);
            System.out.println(callRatePct+"\t"+bins[bin]);
        }
        System.gc();
    }

}
