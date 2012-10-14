/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pipeline;

import java.io.File;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.tagdist.TagsByTaxaBit;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByteHDF5TaxaGroups;


/**
 *
 * @author esb
 */
public class EdPipelines {
    public static void main(String[] args) {
//        TagsByTaxaByteHDF5TaxaGroups rHDF5=new TagsByTaxaByteHDF5TaxaGroups("/Volumes/LaCie/test7KtaxaFQ.hf5");
//        System.out.println(rHDF5.getTaxaCount());
//        String[] names=rHDF5.getTaxaNames();
//        for (String s : names) {
//            System.out.println(s); 
//        }
//        convertTextTagCountsToBinary();
//        convertBinaryTagCountsToText();
//        convertBinaryTBTToText();
//        tagCountsToFastQ();
//        runQseqToTagCountPlugin();
//        runMergeMultipleTagCountPlugin();
     //   runQseqToTBTPlugin();
 //       runSeqToTBTHDF5Plugin();
        runSeqToTBTHDF5PluginMT();
//        runMergeTagsByTaxaFilesPlugin();
//        convertSAMToTOPM();
        //createBigTBT();
        //runTagsToSNPByAlignmentMTPlugin();
//        filter10KMaizeTaxa();
//        mergeTaxaInHapMap();
//        imputeIndelsOnBiparentals();
//        filterMergeImpute10KMaizeTaxa();
 //       basicImputation();

    }
    
    public static void createBigTBT() {
        String[] args = new String[] {
            "-i","/Volumes/nextgen/Zea/build20120110/tbt/",
            "-o","/Volumes/LaCie/zea20120110c510.tbt.byte",
        };

        MergeTagsByTaxaFilesByRowPlugin testClass = new MergeTagsByTaxaFilesByRowPlugin();
        testClass.setParameters(args);
        testClass.performFunction(null);
        
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
        String binaryTBTFileS = "C:/Users/jcg233/Documents/Bioinformatics/NextGen/HapMapV2/test_RandomPairedEndToTBT/tbt/Jeff1_FAKEFLOW_1.tbt.bin";
        String textTBTFileS =   "C:/Users/jcg233/Documents/Bioinformatics/NextGen/HapMapV2/test_RandomPairedEndToTBT/tbt/Jeff1_FAKEFLOW_1.tbt.txt";
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
        String[] args = new String[] {
            "-i", "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/qseq",
            "-k", "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/NAM49_50_ApeKI_key.txt",
            "-e", "ApeKI", // Enzyme used to create the GBS library\n"
            "-s", "200000000", // Max good reads per lane. (Optional. Default is 200,000,000).\n"
            "-c", "1", // Minimum tag count (default is 1).\n"
            "-o", "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/tagCounts"
        };

        QseqToTagCountPlugin testClass = new QseqToTagCountPlugin();
        testClass.setParameters(args);
        testClass.performFunction(null);
    }

    public static void runMergeMultipleTagCountPlugin() {
        String[] args = new String[] {
            "-i", "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/tagCounts", // Input directory containing .cnt files
            "-o", "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/mergedTagCounts/NAM49_50_ApeKI_mergedTags_min10.fastq", //  Output file name
            "-c", "10" // Minimum count of reads to be output (default 1)
          , "-t" // Specifies that reads should be output in FASTQ text format.
        };

        MergeMultipleTagCountPlugin testClass = new MergeMultipleTagCountPlugin();
        testClass.setParameters(args);
        testClass.performFunction(null);
    }

    public static void runQseqToTBTPlugin() {
        String[] args = new String[] {
            "-i", "/Volumes/LaCie/70980AAXX/",
            "-o", "/Users/edbuckler/SolexaAnal/GBS/test/",
            "-k", "/Users/edbuckler/SolexaAnal/GBS/NAM_IBM_MDP_key.txt",
            "-e", "ApeKI", // Enzyme used to create the GBS library\n"
            "-m", "/Volumes/LaCie/zea20120110c510.topm", // master Tags file
            "-y",
            "-c", "0" // Minimum tag count (default is 1).
        };

        QseqToTBTPlugin testClass = new QseqToTBTPlugin();
        testClass.setParameters(args);
        testClass.performFunction(null);
    }
    
    public static void runSeqToTBTHDF5Plugin() {
        String[] args = new String[] {
            "-i", "/Volumes/LaCie/70980AAXX/",
         //   "-i","/Volumes/nextgen-1/GBSFlowCells",
            "-o", "/Volumes/LaCie/testFQSMALL.h5",
            "-L", "/Volumes/LaCie/testFQSeqTBT.log",
            "-k", "/Volumes/LaCie/AllZeaKey20120531.txt",
            "-s","5000000",
            "-e", "ApeKI", // Enzyme used to create the GBS library\n"
            "-m", "/Volumes/LaCie/zea20120110c510.topm", // master Tags file
        };

        SeqToTBTHDF5Plugin testClass = new SeqToTBTHDF5Plugin();
        testClass.setParameters(args);
        testClass.performFunction(null);
    }
    
    public static void runSeqToTBTHDF5PluginMT() {
        String[] args = new String[] {
            "-i", "/Users/edbuckler/SolexaAnal/GBS/70980AAXX",
         //   "-i","/Volumes/nextgen-1/GBSFlowCells",
            "-o", "/Users/edbuckler/SolexaAnal/GBS/test/testFQSMALL.h5",
            "-L", "/Users/edbuckler/SolexaAnal/GBS/test/testFQSeqTBT.log",
            "-k", "/Users/edbuckler/SolexaAnal/GBS/AllZeaKey20120531.txt",
            "-s","5000000",
            "-e", "ApeKI", // Enzyme used to create the GBS library\n"
            "-m", "/Users/edbuckler/SolexaAnal/GBS/build110816/tbttopm/mergedNAM282Ames_072011.topm.bin", // master Tags file
        };

        SeqToTBTHDF5MultiThreadPlugin testClass = new SeqToTBTHDF5MultiThreadPlugin();
        testClass.setParameters(args);
        testClass.performFunction(null);
    }

    public static void runMergeTagsByTaxaFilesPlugin() {
        String[] args = new String[] {
            "-i", "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/tbt",
            "-o", "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/mergedTBT/NAM49_50_ApeKI_mergedTBT_min10.tbt.bin",
        };

        MergeTagsByTaxaFilesPlugin testClass = new MergeTagsByTaxaFilesPlugin();
        testClass.setParameters(args);
        testClass.performFunction(null);
    }

    public static void convertSAMToTOPM() {
        String SAMFile =      "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/mergedTagCounts/NAM49_50_ApeKI_mergedTags_min10.sam";
        String TOPMFile =     "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/mergedTagCounts/NAM49_50_ApeKI_mergedTags_min10.topm.bin";
        String TOPMTextFile = "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/mergedTagCounts/NAM49_50_ApeKI_mergedTags_min10.topm.txt";

        TagsOnPhysicalMap topm = new TagsOnPhysicalMap();
        topm.readSAMFile(SAMFile,2);
        topm.writeBinaryFile(new File(TOPMFile));
        topm.writeTextFile(new File(TOPMTextFile));
    }

    public static void mergeTaxaInHapMap() {
         String[] args = new String[] {
            "-hmp", "/Users/edbuckler/SolexaAnal/GBS/alignment_test_071811/hmp/maize071811.cov10.fT1hLDE1.c+.hmp.txt",
            "-o", "/Users/edbuckler/SolexaAnal/GBS/alignment_test_071811/test/maize071811.cov10.fT1hLDE1.mgNoHet.c+.hmp.txt",
            "-xHets",
            "-hetFreq", "0.76",
            "-sC", "1",  // Start chromosome
            "-eC", "10" // End chromosome
        };

        MergeIdenticalTaxaPlugin testClass = new MergeIdenticalTaxaPlugin();
        testClass.setParameters(args);
        testClass.performFunction(null);
    }

    public static void runTagsToSNPByAlignmentMTPlugin() {
        String[] args = new String[] {
            "-i", "/Volumes/LaCie/zea20120110c510.tbt.byte",
            "-y",
//            "-i", "/Users/edbuckler/SolexaAnal/GBS/build110816/tbttopm/DTMA.tbt.bin",
//            "-i", "/Users/edbuckler/SolexaAnal/GBS/alignment_test_071811/merged_141_lanes_chr3.tbt.bin",

//            "-m", "/Users/edbuckler/SolexaAnal/GBS/build111217/mapping/mergedNam282Ames_072011.topm.bin",
            
            "-m", "/Volumes/LaCie/zea20120110c510.topm",
            "-mUpd", "/Volumes/LaCie/zea20120110c510.prod510.topm",
            "-o", "/Volumes/LaCie/build20120110/test/",
            
//            "-m", "/Users/edbuckler/SolexaAnal/GBS/build111217/mapping/mergedNam282Ames_072011_mappedOnly.topm.bin",
//            "-mUpd", "/Users/edbuckler/SolexaAnal/GBS/build111217/testNoMap/C10mergedNam282Ames_072011_mappedOnly.prod.topm.bin",
//            "-o", "/Users/edbuckler/SolexaAnal/GBS/build111217/testNoMap/",
 //           "-mUpd", "/Volumes/nextgen/Zea/temp/mergedNAM282Ames_072011.allprod.topm.bin",
//            "-o", "/Volumes/nextgen/Zea/temp/",
            
            "-mnF",".8",
  //          "-mnMAF","0.005",
            "-mnMAC","20",
            "-mnLCov", "0.1", // Minimum locus coverage (proportion of Taxa)
            "-s", "5",  // Start chromosome
            "-e", "10" // End chromosome
        };

        TagsToSNPByAlignmentPlugin testClass = new TagsToSNPByAlignmentPlugin();
        testClass.setParameters(args);
        testClass.performFunction(null);
    }

    public static void filter10KMaizeTaxa() {
        String[] args = new String[] {
            "-hmp", "/Users/edbuckler/SolexaAnal/GBS/alignment_test_071811/hmp/alignment_test_071811.c+.cov20.hmp.txt",
            "-o", "/Users/edbuckler/SolexaAnal/GBS/alignment_test_071811/test/alignment_test_071811.cov20.fT1hLD.c+.hmp.txt",
            "-mnTCov", "0.01",
            "-mnF",".9",
            "-mnMAF","0.01",
            "-hLD",
            "-mnSCov", "0.20", // Minimum locus coverage (proportion of Taxa)
            "-sC", "6",  // Start chromosome
            "-eC", "10" // End chromosome
        };

        GBSHapMapFiltersPlugin testClass = new GBSHapMapFiltersPlugin();
        testClass.setParameters(args);
        testClass.performFunction(null);

        args = new String[] {
            "-hmp", "/Users/edbuckler/SolexaAnal/GBS/alignment_test_071811/test/alignment_test_071811.cov20.fT1hLD.c+.hmp.txt",
            "-o", "/Users/edbuckler/SolexaAnal/GBS/alignment_test_071811/test/alignment_test_071811.cov20.fT1hLDE1.c+.hmp.txt",
            "-oB", "/Users/edbuckler/SolexaAnal/GBS/test/errorBin.txt",
            "-oE", "/Users/edbuckler/SolexaAnal/GBS/test/errorBySNP.txt",
            "-popM","Z[0-9]{3}",
            "-sC","1","-eC","10",
            "-mxE","0.01",
            "-mnD","2.0",
            };
        BiParentalErrorCorrection.main(args);
    }

        public static void imputeIndelsOnBiparentals() {
        String rootIn="/Users/edbuckler/SolexaAnal/GBS/build110813/hmp/maize110812";
        String rootOut="/Users/edbuckler/SolexaAnal/GBS/build110813/test/h2Kmaize110812";
        int sC=10;
        int eC=10;
        String[] args = new String[] {
            "-hmp", rootOut+".cov10.fT1.c+.hmp.txt",
            "-o", rootOut+".cov10.fT1Indel.c+.hmp.txt",
            "-popM","Z[0-9]{3}",
             "-sC", ""+sC,  // Start chromosome
            "-eC", ""+eC // End chromosome
        };
        BiParentalIndelCallerPlugin testClass = new BiParentalIndelCallerPlugin();
        testClass.setParameters(args);
        testClass.performFunction(null);

    }

     public static void filterMergeImpute10KMaizeTaxa() {
        String rootIn="/Volumes/LaCie/build20120110/test/Zea20120110";
        String rootOut="/Users/edbuckler/SolexaAnal/GBS/build20120110/test/Zea20120110";
        int sC=5;
        int eC=10;
        String[] args = new String[] {
            "-hmp", rootIn+".c+.hmp.txt",
            "-o", rootOut+"_scv10mF8maf002.c+.hmp.txt",
            "-mnTCov", "0.01",
            "-mnF",".8",
            "-mnMAF","0.002",
//            "-hLD",
            "-mnSCov", "0.10", // Minimum locus coverage (proportion of Taxa)
            "-sC", ""+sC,  // Start chromosome
            "-eC", ""+eC // End chromosome
        };

        GBSHapMapFiltersPlugin testClass = new GBSHapMapFiltersPlugin();
        testClass.setParameters(args);
        //testClass.performFunction(null);
        args = new String[] {
            "-hmp", rootOut+"_scv10mF8maf002.c+.hmp.txt", // Input HapMap file; use a plus sign (+) as a wild card character to specify multiple chromosome numbers
            "-o", rootOut+"_scv10mF8maf002_mgs.c+.hmp.txt", // Output HapMap file
            "-misMat", "0.15", // Threshold genotypic mismatch rate above which the duplicate SNPs won't be merged (default: 0.05)
            "-callHets",       // When two genotypes disagree at a SNP, call it a heterozygote (default: off = set to missing)
//          "-kpUnmergDups",   // When two duplicate SNPs were not merged (different alleles or too many mismatches), keep them (default: off = delete them)
            "-s", ""+sC,        // Start chromosome (default 1)
            "-e", ""+eC         // End chromosome (default 10)
        };

        MergeDuplicateSNPsPlugin plugin = new MergeDuplicateSNPsPlugin();
        plugin.setParameters(args);
        plugin.performFunction(null);
            
        args = new String[] {
            "-hmp", rootOut+"_scv10mF8maf002_mgs.c+.hmp.txt",
            "-o", rootOut+"_scv10mF8maf002_mgs_E1pLD5kpUn.c+.hmp.txt",
            "-oB", "/Users/edbuckler/SolexaAnal/GBS/test/errorBin.txt",
            "-oE", "/Users/edbuckler/SolexaAnal/GBS/test/errorBySNP.txt",
            "-popM","Z[0-9]{3}",
            "-sC",""+sC,
            "-eC",""+eC,
            "-mxE","0.01",
            "-mnD","2.0",
            "-mnPLD", "0.5",
            "-kpUT"
            };
        BiParentalErrorCorrectionPlugin bpec = new BiParentalErrorCorrectionPlugin();
        bpec.setParameters(args);
        bpec.performFunction(null);
System.exit(0);
        args = new String[] {
            "-hmp", rootOut+".cov10.fT1E1pLD.c+.hmp.txt",
            "-o", rootOut+".cov10.fT1E1pLD.mgNoHet.c+.hmp.txt",
            "-xHets",
            "-hetFreq", "0.76",
            "-sC", ""+sC,  // Start chromosome
            "-eC", ""+eC // End chromosome
        };
        MergeIdenticalTaxaPlugin mitp = new MergeIdenticalTaxaPlugin();
        mitp.setParameters(args);
        mitp.performFunction(null);

        args = new String[] {
            "-hmp", rootOut+".cov10.fT1E1pLD.mgNoHet.c+.hmp.txt",
            "-o", rootOut+".cov10.fT1E1pLD.mgNoHet.imp.c+.hmp.txt",
            "-sC", ""+sC,  // Start chromosome
            "-eC", ""+eC // End chromosome
        };
        FastImputationBitFixedWindow.main(args);
    }
     
     public static void basicImputation() {
         String rootIn="/Users/edbuckler/SolexaAnal/GBS/build20120110/test/Zea20120110_scv10mF8maf002_mgs_E1pLD5kpUn";
        String rootOut="/Users/edbuckler/SolexaAnal/GBS/build20120110/test/Zea20120110_scv10mF8maf002_mgs_E1pLD5kpUn";
        int sC=5;
        int eC=10;
          String[]    args = new String[] {
            "-hmp", rootIn+".c+.hmp.txt",
            "-o", rootOut+".imp95_1024.c+.hmp.txt",
            "-sC", ""+sC,  // Start chromosome
            "-eC", ""+eC // End chromosome
        };
        FastImputationBitFixedWindow.main(args);
    }
  
}
