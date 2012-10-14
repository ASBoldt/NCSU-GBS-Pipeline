/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.genome.GBS;

import java.io.File;
import java.io.BufferedReader;
import java.io.FileReader;
import net.maizegenetics.genome.BaseEncoder;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.CombineAlignment;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.FilterAlignment;
import net.maizegenetics.pal.alignment.GdpdmBLOBUtils;
import net.maizegenetics.pal.alignment.ImportUtils;

/**
 *
 * @author edbuckler
 */
public class GBSPipelinesJG {

    public GBSPipelinesJG() {
    }

    public static void main(String[] args) {
//        createRefMarkerMap();
//        looseMatchesToWorkingReads();
//        makeReadCountsFileForTwoIBMFlowCells();
//        performVD();
//        loadExternalVD();
//        mergeVDs();
//        makeReadCountsTeosinteTestPlate();
//        makeTeosinteTestPlateGenoFile();
//        writeTeosinteTestPlateAlleleInfoFile();
//        makeReadCountsTeosinte();
//        makeReadCountsTeosinteDuplicates();
//        makeTeosinteGenoFile();
//        writeTeosinteAlleleInfoFile();
//        alignB73Htrhm();
//        alignB73();
//        analyzeIBMforGBSmethodsPaper();
//        testDuplicateTaxon();
//        makeReadsByTaxaTwoIBMFlowCells();
//        makePhysicalMapForTwoIBMFlowCells();
//        mapTwoIBMFlowCellsVs55KFrame();
//        performPstIVD();
//        makeReadCountsFileForIBMPstI();
//        makePhysicalMapForIBMPstI();
//        mapIBMPstITestVs55KFrame();
//        imputeNAMSNPFrameworkMap();
//        makeReadCountsFileForNAMPstI();
//        makePhysicalMapForNAMPstI();
//        mapNAMPstITestVsSNPFrame();
//        makeReadCountsFilesForNAMApeKI();
//        makeTagsByTaxaForNAMApeKI();
//        makePhysicalMapForNAMApeKI();
//        mapNAMApeKIVsSNPFrame();
        makeReadCountsFilesForCintaFengFineMapping();
//        makeReadCountsSorghum();
//        physicallyMapSorghumTags();
//        makeGBSFrameworkSorghum();
//        filterGBSFrameworkSorghum();
//        mapSorghumVsImputedGBSFrame();
//        makeReadCountsTeoW22QTL10();
//        writeTeoPhyGenMappedReads();
    }

    public static void createRefMarkerMap() {
        /**
         * Virtual digest of a reference
         */
        String baseDir = "/usr/local/maizediv/";
        String refGenomeFasta = baseDir + "genome/maize_agp_v2.fasta";
        String refGenomeDigested = baseDir + "genome/virtual_digest/B73refVirtualV2.rwpm.bin";
        String qseqDir = baseDir + "illumina/434LFAAXX/qseq";
        String qseqKey = baseDir + "illumina/434LFAAXX/434LFAAXX_Multiplex_key.txt";
        String parsedTaxaDir = baseDir + "illumina/434LFAAXX/taxa";
        String collapsedTaxaDir = baseDir + "illumina/434LFAAXX/taxacollapse";
        String readsWorthMapping = baseDir + "illumina/434LFAAXX/counts/ReadCountsMin10_100421b.bin";
        String obsCommonReadsOnMap = baseDir + "illumina/434LFAAXX/counts/commonReadMap_100422.rwpm.bin";

        VirtualDigester vd = new VirtualDigester(new File(refGenomeFasta), new File(refGenomeDigested));
        System.gc();

        //Sort the physical map for quick searching and save
        ReadsWPhysicalMap rwpmVD = new ReadsWPhysicalMap(refGenomeDigested, true);
        rwpmVD.sortTable(true);
        rwpmVD.writeCountFile(new File(refGenomeDigested));
        System.gc();

        //Read and parse qseq or fastq file
       ParseBarcodeFiles pbf=new ParseBarcodeFiles(qseqDir,qseqKey,parsedTaxaDir,10,false);
       System.gc();

        //Take a GBSReads folder and collapse the folder by sorting and collapse identical reads
       ReadCounts rc=new ReadCounts(parsedTaxaDir, collapsedTaxaDir, 1, true, false, false);
       System.gc();

    //Fuse multiple taxa files into one, minimum frequency is 2, binary=true, simpleFilter=false
    //Use the physical genome as a base, and then add high quality reads to the top
       ReadsWPhysicalMap theRef=new ReadsWPhysicalMap(refGenomeDigested,true);
       theRef.sortTable(true);
//       CombineReadCounts shff=new CombineReadCounts(theRef,qualcollapsedTaxaDir,combinedFlowcellsReads, 2, true, false);
       System.gc();

       //Found that using the reference reads as all greater than 2 in frequency was the best
       //so readsWorthMapping=combinedFlowcellsReads
       //now use all the reads, but only populate with those with perfect matches
//       CreateReadsByTaxa crbt=new CreateReadsByTaxa(readsWorthMapping, collapsedTaxaDir, theReadsByTaxaAll, true);
       System.gc();

//       ReadsByTaxa theRBT=new ReadsByTaxa(theReadsByTaxaAll,true);
//       theRBT.writeDistFile(new File(theReadsByTaxaMinCount), true, 5);
//       theRBT.writeReadCountFile(new File(readsWorthMapping), true, 5);
//       System.gc();

      //Create an empty mapping list from the readList

       ReadCounts rc3=new ReadCounts(readsWorthMapping, true);
       ReadsWPhysicalMap rwpm1=new ReadsWPhysicalMap(rc3);
       rwpm1.writeCountFile(new File(obsCommonReadsOnMap));
       System.gc();

//       ReadsWPhysicalMap rwpm2=new ReadsWPhysicalMap(obsCommonReadsOnMap,true);


         ReadsWPhysicalMap rwpmTest2=MapHaplotypes.blastUsingPhysicalMap(new ReadsWPhysicalMap(refGenomeDigested,true),
                  new ReadsWPhysicalMap(obsCommonReadsOnMap,true),3);
        rwpmTest2.writeCountFile(new File(obsCommonReadsOnMap));


//         ReadsByTaxa theRBT2=new ReadsByTaxa(theReadsByTaxaMinCount,true);
         System.gc();
        System.out.println("rwpmTest2 memory"+rwpmTest2);
//         rwpmTest2=MapHaplotypes.checkWorkingMapWithGenetic(baseGeneticMap,rwpmTest2,theRBT2,true,0.00001);
//         rwpmTest2.writeCountFile(new File(obsCommonReadsWithGenetic1),Integer.MAX_VALUE, true,true,0.01f,true);


    }

    public static void looseMatchesToWorkingReads() {
        String baseDir = "/usr/local/maizediv/";
        String refGenomeDigested = baseDir + "genome/virtual_digest/B73refVirtualV2.rwpm.bin";
        String mo17CtgsDigested = baseDir + "genome/Mo17/Mo17_454AllContigs_chr0.cut.bin";
        String GBSFrameWorkMapSuffix = baseDir + "illumina/434LFAAXX/GBS_framework/commonReadMap_100508t";
        String qseqDir = baseDir + "illumina/434LFAAXX/qseq";
        String qseqKey = baseDir + "illumina/434LFAAXX/434LFAAXX_Multiplex_key.txt";
        String parsedTaxaDir = baseDir + "illumina/434LFAAXX/taxa";
        String collapsedTaxaDir = baseDir + "illumina/434LFAAXX/taxacollapse";
        String referenceReadsFile = baseDir + "illumina/434LFAAXX/counts/RefReads_Q10_B73_Mo17_20100525.bin";
        String stringentReadsByTaxaFile = baseDir + "illumina/434LFAAXX/counts/ReadByTaxa_stringent_10taxa_20100525.bin";
        String RefRwPMFile = baseDir + "illumina/434LFAAXX/counts/RefReadsWPM_Q10_B73_Mo17_20100525.bin";
        String fullGeneticOuputFile = baseDir + "illumina/434LFAAXX/counts/geneticOutputWSkim_stringentRBT_20100526.txt";


        // use IBM good run to pick initial set of reference (working) reads: minQual 10
        ParseBarcodeFiles pbf = new ParseBarcodeFiles(qseqDir, qseqKey, parsedTaxaDir, 10, false);
        ReadCounts rc = new ReadCounts(parsedTaxaDir, collapsedTaxaDir, 1, true, false, false);

        // create reference reads from the B73 & Mo17 virtual digests combined with the high quality GBS reads from above
        ReadsWPhysicalMap b73Ref = new ReadsWPhysicalMap(refGenomeDigested, true);
        b73Ref.sortTable(true);

        ReadsWPhysicalMap mo17Ref=new ReadsWPhysicalMap(mo17CtgsDigested, true);
        mo17Ref.sortTable(true);

        ReadCounts referenceReads = new CombineReadCounts(b73Ref, mo17Ref, collapsedTaxaDir, referenceReadsFile, 1, true);

        // now get all of the GBS reads, regardless of quality
        pbf = new ParseBarcodeFiles(qseqDir, qseqKey, parsedTaxaDir, 0, false);
        rc = new ReadCounts(parsedTaxaDir, collapsedTaxaDir, 1, true, false, false);

        // populate the ReadsByTaxa with reads that match perfectly, only write the reads that have data in 10 or more taxa
        CreateReadsByTaxa myCRBT = new CreateReadsByTaxa(referenceReadsFile, collapsedTaxaDir, stringentReadsByTaxaFile, 10, true);

        // create an empty mapping list from the reference readList
        ReadCounts refReads=new ReadCounts(referenceReadsFile, true);
        ReadsWPhysicalMap RefRwPM=new ReadsWPhysicalMap(refReads);
        RefRwPM.writeCountFile(new File(RefRwPMFile));
        System.gc();

        // populate the physical positions on the reference rwpm
        ReadsWPhysicalMap rwpmRef=MapHaplotypes.blastUsingPhysicalMap(new ReadsWPhysicalMap(refGenomeDigested,true),
                  new ReadsWPhysicalMap(RefRwPMFile,true),3);
        rwpmRef.writeCountFile(new File(RefRwPMFile));

        // load the framework (imputed) genetic map
        Alignment[] GBSFrameworkMap = new Alignment[10];
        for (int chromosome = 1; chromosome <= 10; chromosome++) {
            GBSFrameworkMap[chromosome - 1] = ImportUtils.readFromHapmap(GBSFrameWorkMapSuffix + ".imc" + chromosome + ".hmp.txt");
            System.out.println("Chromosome " + chromosome + " GeneticMap  Taxa:" + GBSFrameworkMap[chromosome - 1].getSequenceCount() + " Sites:" + GBSFrameworkMap[chromosome - 1].getSiteCount());
        }

        // open the reference rwpm and stringent RBT
        RefRwPM = new ReadsWPhysicalMap(RefRwPMFile, true);
        ReadsByTaxa stringentRBT = new ReadsByTaxa(stringentReadsByTaxaFile, true);

        // map the high stringency RBT
        MapHaplotypes.mapWorkingMapWithGenetic(GBSFrameworkMap, RefRwPM, stringentRBT, 0.02, 10, 10, fullGeneticOuputFile);

        // look for looser matching reads (try 48, then 32 bases): compare mapping with and wout the looser reads against the reference map

    }

    public static void performVD() {
        String Mo17JGI454Chr80Contigfile = "Q:/SolexaAnal/IBM/virtual_digests/JGI_Mo17_454/Mo17JGI454AllContigs_chr70.fasta";
        String Mo17JGI454Chr80VDfile = "Q:/SolexaAnal/IBM/virtual_digests/JGI_Mo17_454/Mo17JGI454VD_chr70.rwpm.bin";
        String AGPv1Fastafile = "/usr/local/maizediv/genome/maize_agp_v1.fasta";
        String AGPv1VDfile = "/usr/local/maizediv/genome/virtual_digest/B73refVirtualV1.rwpm.bin";
        String AGPv2Fastafile = "Q:/SolexaAnal/AGPv2/maize_agp_v2.fasta";
        String AGPv2VDfile = "Q:/SolexaAnal/IBM/virtual_digests/AGPv2/B73AGPv2VD.rwpm.bin";
        String SorghumFasta = "/usr/local/maizediv/genome/sorghum/Sbicolor_79_RM_042611_10chrs.fa";
        String SorghumVDFile = "/usr/local/maizediv/genome/sorghum/Sbicolor_79_RM_042611_10chrs_VD.rwpm.bin";

        VirtualDigester myVDigester=new VirtualDigester(new File(SorghumFasta), new File(SorghumVDFile));
        System.gc();

        //Sort the physical map for quick searching and save
        ReadsWPhysicalMap rwpmMyVD = new ReadsWPhysicalMap(SorghumVDFile, true);
        rwpmMyVD.sortTable(true);
        rwpmMyVD.writeCountFile(new File(SorghumVDFile));
        System.gc();
    }

    public static void loadExternalVD() {
        String CSHLB73454VDfile = "Q:/SolexaAnal/IBM/virtual_digests/CSHL_Shiran_B73/454_tags_processed.txt";  // 1135193 tags
        String CSHLB73454VDoutfile = "Q:/SolexaAnal/IBM/virtual_digests/CSHL_Shiran_B73/B73_454_VD.rwpm.bin";
        String BGI_Mo17VDfile = "Q:/SolexaAnal/IBM/virtual_digests/BGI_Mo17/Mo17_virtualDigest_processed.txt"; // 3761353 tags
        String BGI_Mo17VDoutfile = "Q:/SolexaAnal/IBM/virtual_digests/BGI_Mo17/Mo17_BGI_VD.rwpm.bin";
        String B73K50VD454file = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/ShiranVDs/454asm_k50_tags.txt";  //  7680513 tags
        String B73K50VD454outfile = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/ShiranVDs/454k50VD.rwpm.bin";
        String B73K60VD454file = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/ShiranVDs/454asm_k60_tags.txt";  //  8539236 tags
        String B73K60VD454outfile = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/ShiranVDs/454k60VD.rwpm.bin";
        String B73K80VD454file = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/ShiranVDs/454asm_k80_tags.txt";  //  8580948 tags
        String B73K80VD454outfile = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/ShiranVDs/454k80VD.rwpm.bin";
        String B73K96VD454file = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/ShiranVDs/454asm_k96_tags.txt";  // 10426055 tags
        String B73K96VD454outfile = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/ShiranVDs/454k96VD.rwpm.bin";
        String B73K96IIVD454file =    "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/de_novo_ctgs/454k96II/k96II_tags.txt";  // 9757598 tags
        String B73K96IIVD454outfile = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/de_novo_ctgs/454k96II/454k96IIVD.rwpm.bin";
        String flcDNAvdFile = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/ShiranVDs/flcdna_tags.txt";     //   733425 tags; 727955 with ACGT only
        String flcDNAvdOutFile = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/ShiranVDs/flcdna.rwpm.bin";

        int nRows = 9757598;  //  <-----******CHANGE THIS ACCORDING TO THE SIZE OF EACH VD
        int maxTags = 10500000 + 1;
        String temp = "";
        String[] cells = new String[5];
        String ReadStr;
        long[] parts = new long[2];
        int chunkSize = BaseEncoder.chunkSize;
        byte chromosome;
        char strand; //'+', '-', "?"
        int posMin, posMax;

        try {
            FileReader fr = new FileReader(B73K96IIVD454file); //  <-----CHANGE infile for the particular VD
            BufferedReader br = new BufferedReader(fr);
            ReadsWPhysicalMap rwpmExternalVD = new ReadsWPhysicalMap(nRows);
            int nAdded = 0;
            while ((br.ready()) && (rwpmExternalVD.getReadTotal() < (maxTags - 2))) {
                temp = br.readLine().trim();
                cells = temp.split("\t");
                ReadStr = cells[0];
                parts[0] = BaseEncoder.getLongFromSeq(ReadStr.substring(0, chunkSize));
                parts[1] = BaseEncoder.getLongFromSeq(ReadStr.substring(chunkSize, 2 * chunkSize));
                chromosome = Byte.parseByte(cells[1]);
                strand = cells[2].charAt(0);
                posMin = Integer.parseInt(cells[3]);
                posMax = Integer.parseInt(cells[4]);
                if ((parts[0] != -1) && (parts[1] != -1)) {  // Nb: virtual reads containing N will encode as null
                    rwpmExternalVD.addHaplotype(parts, chromosome, (byte) strand, posMin, posMax, Short.MIN_VALUE, (byte) 0);
                    ++nAdded;
                }
            }
            System.out.println();
            System.out.println(nAdded + " reads added to rwpmExternalVD");
            System.out.println();
            rwpmExternalVD.sortTable(true);
            rwpmExternalVD.writeCountFile(new File(B73K96IIVD454outfile));  //  <-----CHANGE outfile for the particular VD
            rwpmExternalVD.printRows(1000);

        } catch (Exception e) {
            System.out.println("Catch: e=" + e);
            e.printStackTrace();
            System.exit(1);
        }
    }

    public static void mergeVDs() {
        String AGPv2VDfile = "Q:/SolexaAnal/IBM/virtual_digests/AGPv2/B73AGPv2VD.rwpm.bin";
        String Mo17JGI454VDfile = "Q:/SolexaAnal/IBM/virtual_digests/JGI_Mo17_454/Mo17JGI454VD_chr70.rwpm.bin";
        String B73CSHL454VDfile = "Q:/SolexaAnal/IBM/virtual_digests/CSHL_Shiran_B73/B73_454_VD.rwpm.bin";
        String BGIMo17VDfile = "Q:/SolexaAnal/IBM/virtual_digests/BGI_Mo17/Mo17_BGI_VD.rwpm.bin";
        String MasterVDoutfile = "Q:/SolexaAnal/IBM/virtual_digests/MasterVD.rwpm.bin";
        
        ReadsWPhysicalMap rwpmMasterVD = new ReadsWPhysicalMap(AGPv2VDfile, true, 15620352);

        ReadsWPhysicalMap Mo17JGI454VD = new ReadsWPhysicalMap(Mo17JGI454VDfile, true);
        for (int i = 0; i < Mo17JGI454VD.getReadTotal(); i++) {
            if (rwpmMasterVD.addHaplotype(Mo17JGI454VD.getReadWithPosition(i)) < 0) {
                System.out.println("Ran out of room in rwpmMasterVD. Current size = " + rwpmMasterVD.getSize());
            }
        }
        Mo17JGI454VD = null;
        System.gc();

        ReadsWPhysicalMap B73CSHL454VD = new ReadsWPhysicalMap(B73CSHL454VDfile, true);
        for (int i = 0; i < B73CSHL454VD.getReadTotal(); i++) {
            if (rwpmMasterVD.addHaplotype(B73CSHL454VD.getReadWithPosition(i)) < 0) {
                System.out.println("Ran out of room in rwpmMasterVD. Current size = " + rwpmMasterVD.getSize());
            }
        }
        B73CSHL454VD = null;
        System.gc();

        ReadsWPhysicalMap BGIMo17VD = new ReadsWPhysicalMap(BGIMo17VDfile, true);
        for (int i = 0; i < BGIMo17VD.getReadTotal(); i++) {
            if (rwpmMasterVD.addHaplotype(BGIMo17VD.getReadWithPosition(i)) < 0) {
                System.out.println("Ran out of room in rwpmMasterVD. Current size = " + rwpmMasterVD.getSize());
            }
        }
        BGIMo17VD = null;
        System.gc();

        rwpmMasterVD.sortTable(true);
        rwpmMasterVD.writeCountFile(new File(MasterVDoutfile));
        rwpmMasterVD.printRows(10000);
    }

    public static void makeReadCountsFileForTwoIBMFlowCells() {
        String qseqDirIBM1 =           "/cbsufsrv3/data2/maizediv/illumina/434LFAAXX/qseq";
        String qseqKeyIBM1 =           "/cbsufsrv3/data2/maizediv/illumina/434LFAAXX/434LFAAXX_Multiplex_key.txt";
        String parsedTaxaDirIBM1 =     "/cbsufsrv3/data2/maizediv/illumina/434LFAAXX/taxa";
        String collapsedTaxaDirIBM1 =  "/cbsufsrv3/data2/maizediv/illumina/434LFAAXX/taxacollapse";
        String combinedReadsFileIBM1 = "/cbsufsrv3/data2/maizediv/illumina/434LFAAXX/counts/CombReads_IBM1_Min5_20100525.bin";

        String qseqDirIBM2 =           "/cbsufsrv3/data2/maizediv/illumina/434GFAAXX/qseq";
        String qseqKeyIBM2 =           "/cbsufsrv3/data2/maizediv/illumina/434GFAAXX/434GFAAXX_Multiplex_key.txt";
        String parsedTaxaDirIBM2 =     "/cbsufsrv3/data2/maizediv/illumina/434GFAAXX/taxa";
        String collapsedTaxaDirIBM2 =  "/cbsufsrv3/data2/maizediv/illumina/434GFAAXX/taxacollapse";
        String combinedReadsFileIBM2 = "/cbsufsrv3/data2/maizediv/illumina/434GFAAXX/counts/CombReads_IBM2_Min5_20100525.bin";

        String combinedReadsFileIBM =  "/cbsufsrv3/data2/maizediv/illumina/434GFAAXX/counts/CombReads_BothIBM_Min10_20100525.bin";

        // START WITH IBM2 (434GFAAXX)
        // Read and parse qseq or fastq file using low stringency Q score (maximal # of reads)
        System.out.println();
        System.out.println("Parsing the qseq files for IBM2 (434GFAAXX)");
        System.out.println("-------------------------------------------");
        ParseBarcodeFiles pbfIBM2=new ParseBarcodeFiles(qseqDirIBM2, qseqKeyIBM2, parsedTaxaDirIBM2, 0, false);
        pbfIBM2 = null;
        System.gc();

        // condense the taxa ReadCounts files by sorting and collapsing identical reads
        System.out.println();
        System.out.println("Collapsing the ReadCounts files for each taxon in IBM2 (434GFAAXX)");
        System.out.println("------------------------------------------------------------------");
        ReadCounts rcIBM2=new ReadCounts(parsedTaxaDirIBM2, collapsedTaxaDirIBM2, 1, true, false, false);
        rcIBM2 = null;
        System.gc();

        // Fuse multiple taxa files into one, minimum frequency is 5, binary=true (reads that are seen repeatedly are hopefully real)
        System.out.println();
        System.out.println("Combining the ReadCounts files across all taxa in IBM2 (434GFAAXX)");
        System.out.println("------------------------------------------------------------------");
        CombineReadCounts crcIBM2=new CombineReadCounts(collapsedTaxaDirIBM2, combinedReadsFileIBM2, 5, true);
        crcIBM2 = null;
        System.gc();

        // REPEAT FOR IBM1 (434LFAAXX)
        // Read and parse qseq or fastq file using low stringency Q score (maximal # of reads)
        System.out.println();
        System.out.println("Parsing the qseq files for IBM1 (434LFAAXX)");
        System.out.println("-------------------------------------------");
        ParseBarcodeFiles pbfIBM1=new ParseBarcodeFiles(qseqDirIBM1, qseqKeyIBM1, parsedTaxaDirIBM1, 0, false);
        System.gc();

        // condense the taxa ReadCounts files by sorting and collapsing identical reads
        System.out.println();
        System.out.println("Collapsing the ReadCounts files for each taxon in IBM1 (434LFAAXX)");
        System.out.println("------------------------------------------------------------------");
        ReadCounts rcIBM1=new ReadCounts(parsedTaxaDirIBM1, collapsedTaxaDirIBM1, 1, true, false, false);
        rcIBM1 = null;
        System.gc();

        // Fuse multiple taxa files into one, minimum frequency is 5, binary=true (reads that are seen repeatedly are hopefully real)
        System.out.println();
        System.out.println("Combining the ReadCounts files across all taxa in IBM1 (434LFAAXX)");
        System.out.println("------------------------------------------------------------------");
        CombineReadCounts crcIBM1=new CombineReadCounts(collapsedTaxaDirIBM1, combinedReadsFileIBM1, 5, true);
        crcIBM1 = null;
        System.gc();

        // COMBINE TOGETHER IBM1 & 1BM2
        System.out.println();
        System.out.println("Combining IBM1 & IBM2 (minCount of 10)");
        System.out.println("--------------------------------------");
        ReadCounts crcIBM = new CombineReadCounts(combinedReadsFileIBM1, combinedReadsFileIBM2, combinedReadsFileIBM, 10, true);
    }

    public static void makeReadCountsFileForThirdIBMFlowCell() {
        String qseqDirIBM3 =           "/cbsufsrv3/data2/maizediv/illumina/61VBRAAXX/qseq";
        String qseqKeyIBM3 =           "/cbsufsrv3/data2/maizediv/illumina/IBM/61VBRAAXX/61VBRAAXX_IBM_key.txt";
        String parsedTaxaDirIBM3 =     "/cbsufsrv3/data2/maizediv/illumina/IBM/61VBRAAXX/taxa";
        String collapsedTaxaDirIBM3 =  "/cbsufsrv3/data2/maizediv/illumina/IBM/61VBRAAXX/taxacollapse";
        String combinedReadsFileIBM3 = "/cbsufsrv3/data2/maizediv/illumina/IBM/61VBRAAXX/counts/CombReads_IBM3_Min3_20110315.bin";

        String combinedReadsFileIBM12 =  "/cbsufsrv3/data2/maizediv/illumina/434GFAAXX/counts/CombReads_BothIBM_Min10_20100525.bin";
        String combinedReadsFileIBM123 =  "/cbsufsrv3/data2/maizediv/illumina/IBM/61VBRAAXX/counts/CombReads_ThreeIBM_Min10_20110315.bin";

        // Read and parse qseq or fastq file using low stringency Q score (maximal # of reads)
        System.out.println();
        System.out.println("Parsing the qseq files for IBM3 (61VBRAAXX)");
        System.out.println("-------------------------------------------");
        ParseBarcodeFiles pbfIBM3=new ParseBarcodeFiles(qseqDirIBM3, qseqKeyIBM3, parsedTaxaDirIBM3, 0, false);
        pbfIBM3 = null;
        System.gc();

        // condense the taxa ReadCounts files by sorting and collapsing identical reads
        System.out.println();
        System.out.println("Collapsing the ReadCounts files for each taxon in IBM3 (61VBRAAXX)");
        System.out.println("------------------------------------------------------------------");
        ReadCounts rcIBM3=new ReadCounts(parsedTaxaDirIBM3, collapsedTaxaDirIBM3, 1, true, false, false);
        rcIBM3 = null;
        System.gc();

        // Fuse multiple taxa files into one, minimum frequency is 3, binary=true (reads that are seen repeatedly are hopefully real)
        System.out.println();
        System.out.println("Combining the ReadCounts files across all taxa in IBM3 (61VBRAAXX)");
        System.out.println("------------------------------------------------------------------");
        CombineReadCounts crcIBM3=new CombineReadCounts(collapsedTaxaDirIBM3, combinedReadsFileIBM3, 3, true);
        crcIBM3 = null;
        System.gc();

        // COMBINE TOGETHER IBM12 & 1BM3
        System.out.println();
        System.out.println("Combining IBM12 & IBM3 (minCount of 10)");
        System.out.println("--------------------------------------");
        ReadCounts crcIBM = new CombineReadCounts(combinedReadsFileIBM12, combinedReadsFileIBM3, combinedReadsFileIBM123, 10, true);
    }

    public static void makeReadsByTaxaThreeIBMFlowCells() {
        String parsedTaxaDirIBM =     "/cbsufsrv3/data2/maizediv/illumina/IBM/taxa3FCs";
        String collapsedTaxaDirIBM =  "/cbsufsrv3/data2/maizediv/illumina/IBM/taxacollapse3FCs";
        String combinedReadsFileIBM = "/cbsufsrv3/data2/maizediv/illumina/IBM/61VBRAAXX/counts/CombReads_ThreeIBM_Min10_20110315.bin";
        String rbtFileIBM =           "/cbsufsrv3/data2/maizediv/illumina/IBM/counts3FCs/rbt3IBM_Min10_20110315.bin";

        // condense the taxa ReadCounts files by sorting and collapsing identical reads
        System.out.println();
        System.out.println("Collapsing the ReadCounts files for each taxon in IBM (3 flow cells)");
        System.out.println("-------------------------------------------------------------------");
        ReadCounts rcIBM=new ReadCounts(parsedTaxaDirIBM, collapsedTaxaDirIBM, 1, true, false, true); //binary=true,simpleFilter=false,combineIdenticalTaxa=true
        rcIBM = null;
        System.gc();

        CreateReadsByTaxa teoRBT=new CreateReadsByTaxa(combinedReadsFileIBM, collapsedTaxaDirIBM, rbtFileIBM, true);
    }

    public static void makePhysicalMapForThreeIBMFlowCells() {
        String combinedReadsFileIBM =  "/cbsufsrv3/data2/maizediv/illumina/IBM/counts3FCs/CombReads_ThreeIBM_Min10_20110315.bin";
        String AGPv2VDfile =           "/usr/local/maizediv/genome/virtual_digest/B73refVirtualV2.rwpm.bin";
        String phyMappedReadsFileIBM = "/cbsufsrv3/data2/maizediv/illumina/IBM/counts3FCs/phyMappedReadsIBM_3FC_Div20_20110315.bin";

        ReadsWPhysicalMap rwpmIBMmin10=MapHaplotypes.blastUsingPhysicalMap(new ReadsWPhysicalMap(AGPv2VDfile,true),
                  new ReadsWPhysicalMap(new ReadCounts(combinedReadsFileIBM, true)), 20);  // allow 20 mismatches
        rwpmIBMmin10.sortTable(false); // sort by position
        rwpmIBMmin10.writeCountFile(new File(phyMappedReadsFileIBM), Integer.MAX_VALUE, false, false, Float.NaN, true);
    }

    public static void mapThreeIBMFlowCellsVs55KFrame() {
        String frameWorkGenMapFile55K =   "/usr/local/maizediv/illumina/IBM/counts/IBM_55K_08192010.hmp.txt";
//        String frameWorkGenMapFile55K =   "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/ShiranVDs/IBM_55K_08192010.hmp.txt";
//        String frameWorkGenMapFileGBS =   "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/434LFAAXX/anchor_map/IBM_GBS_framework_map_20100917.hmp.txt";
        String rbtFileIBM =            "/cbsufsrv3/data2/maizediv/illumina/IBM/counts3FCs/rbt3IBM_Min10_20110315.bin";
//        String rbtFileIBM =            "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/ShiranVDs/rbtIBM_Min10_20100827.bin";
        String phyMappedReadsFileIBM = "/cbsufsrv3/data2/maizediv/illumina/IBM/counts3FCs/phyMappedReadsIBM_3FC_Div20_20110315.bin";
//        String phyMappedReadsFileIBM = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/ShiranVDs/phyMappedReadsIBMDiv20_20100827.bin";
        String GeneticResultsFileIBM = "/cbsufsrv3/data2/maizediv/illumina/IBM/counts3FCs/IBMgeneticResults3FCs_20110516.txt";
//        String GeneticResultsFileIBM = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/ShiranVDs/IBMgeneticResults_20100828_tabs.txt";

        // genetically map every tag in IBM that shows up 10 or more times across the two flow cells
        Alignment[] gMap55K= ((CombineAlignment)ImportUtils.readFromHapmap(frameWorkGenMapFile55K)).getAlignments();
        MapHaplotypes.mapWorkingMapWithGenetic(gMap55K, new ReadsWPhysicalMap(phyMappedReadsFileIBM,true),
                new ReadsByTaxa(rbtFileIBM,true), 0.05, 10, 10, GeneticResultsFileIBM);
    }


    public static void makeReadsByTaxaTwoIBMFlowCells() {
        String parsedTaxaDirIBM =     "/cbsufsrv3/data2/maizediv/illumina/IBM/taxa";
        String collapsedTaxaDirIBM =  "/usr/local/maizediv/illumina/IBM/taxacollapse";
        String combinedReadsFileIBM = "/usr/local/maizediv/illumina/IBM/counts/CombReads_BothIBM_Min10_20100525.bin";
        String rbtFileIBM =           "/usr/local/maizediv/illumina/IBM/counts/rbtIBM_Min10_20100827.bin";

        // condense the taxa ReadCounts files by sorting and collapsing identical reads
        System.out.println();
        System.out.println("Collapsing the ReadCounts files for each taxon in IBM (2 flow cells)");
        System.out.println("-------------------------------------------------------------------");
        ReadCounts rcIBM=new ReadCounts(parsedTaxaDirIBM, collapsedTaxaDirIBM, 1, true, false, true); //binary=true,simpleFilter=false,combineIdenticalTaxa=true
        rcIBM = null;
        System.gc();

        CreateReadsByTaxa teoRBT=new CreateReadsByTaxa(combinedReadsFileIBM, collapsedTaxaDirIBM, rbtFileIBM, true);
    }

    public static void makePhysicalMapForTwoIBMFlowCells() {
        String combinedReadsFileIBM =  "/usr/local/maizediv/illumina/IBM/counts/CombReads_BothIBM_Min10_20100525.bin";
        String AGPv2VDfile =           "/usr/local/maizediv/genome/virtual_digest/B73refVirtualV2.rwpm.bin";
        String phyMappedReadsFileIBM = "/usr/local/maizediv/illumina/IBM/counts/phyMappedReadsIBMDiv20_20100827.bin";

        ReadsWPhysicalMap rwpmIBMmin10=MapHaplotypes.blastUsingPhysicalMap(new ReadsWPhysicalMap(AGPv2VDfile,true),
                  new ReadsWPhysicalMap(new ReadCounts(combinedReadsFileIBM, true)), 20);  // allow 20 mismatches
        rwpmIBMmin10.sortTable(false); // sort by position
        rwpmIBMmin10.writeCountFile(new File(phyMappedReadsFileIBM), Integer.MAX_VALUE, false, false, Float.NaN, true);
    }

    public static void mapTwoIBMFlowCellsVs55KFrame() {
//        String frameWorkGenMapFile =   "/usr/local/maizediv/illumina/IBM/counts/IBM_55K_08192010.hmp.txt";
        String frameWorkGenMapFile55K =   "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/ShiranVDs/IBM_55K_08192010.hmp.txt";
        String frameWorkGenMapFileGBS =   "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/434LFAAXX/anchor_map/IBM_GBS_framework_map_20100917.hmp.txt";
//        String rbtFileIBM =            "/usr/local/maizediv/illumina/IBM/counts/rbtIBM_Min10_20100827.bin";
        String rbtFileIBM =            "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/ShiranVDs/rbtIBM_Min10_20100827.bin";
//        String phyMappedReadsFileIBM = "/usr/local/maizediv/illumina/IBM/counts/phyMappedReadsIBMDiv20_20100827.bin";
        String phyMappedReadsFileIBM = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/ShiranVDs/phyMappedReadsIBMDiv20_20100827.bin";
//        String GeneticResultsFileIBM = "/usr/local/maizediv/illumina/IBM/counts/IBMgeneticResults_20100828.txt";
        String GeneticResultsFileIBM = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/ShiranVDs/IBMgeneticResults_20100828_tabs.txt";
        String AGPv1VDfile = "Q:/SolexaAnal/virtual_digests/AGPv1/B73AGPv1VD.rwpm.bin";
        String IBMAlleleInfoOutFile = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/ShiranVDs/IBMAlleleInfo_20100914.txt";
        String readsToKeepFile = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/Chr0/mappedReadsFromCtg697.txt";
        String filteredRBTFile = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/Chr0/rbtCtg545_haplo.txt";
        String chr0CtgHaplotypesFile = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/Chr0/Chr0CtgHaplotypes.txt";
        String chr0CtgHaplotypesResultsFile = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/Chr0/Chr0CtgHaplotypesGeneticResults20100922_GBSFrame2.txt";

        // genetically map every tag in IBM that shows up 10 or more times across the two flow cells
        Alignment[] gMap55K= ((CombineAlignment)ImportUtils.readFromHapmap(frameWorkGenMapFile55K)).getAlignments();
        MapHaplotypes.mapWorkingMapWithGenetic(gMap55K, new ReadsWPhysicalMap(phyMappedReadsFileIBM,true),
                new ReadsByTaxa(rbtFileIBM,true), 0.05, 10, 10, GeneticResultsFileIBM);

        // get the AGPv1 Chr0 position of the mapped tags.  This allows me to figure out which contig they belong to
        AlleleInfo ai = new AlleleInfo(GeneticResultsFileIBM, AGPv1VDfile, IBMAlleleInfoOutFile, 7);  // the 64 base tag of the B73 allele is in column 7 (0 based)

        // make a RBT text file for each ctg that contains only the consensus tags (these were determined manually in Excel)
        ReadsByTaxa rbtIBM = new ReadsByTaxa(rbtFileIBM, true);
        rbtIBM.filterForListOfReads(new File(readsToKeepFile), new File(filteredRBTFile), false);

        // map the Chr0 haplotypes (merged data across multiple, consensus tags) vs both the 55K and GBS framework maps
        Alignment[] gMapGBS= ((CombineAlignment)ImportUtils.readFromHapmap(frameWorkGenMapFileGBS)).getAlignments();
        ReadsWPhysicalMap phyMappedReadsIBM = new ReadsWPhysicalMap(phyMappedReadsFileIBM,true);
        ReadsByTaxa chr0CtgRBT = new ReadsByTaxa(chr0CtgHaplotypesFile,false);  // NOTE: this file must be sorted by tag for the genetic mapping to work properly
        for (int i = 0; i < chr0CtgRBT.getReadTotal(); ++i) {
            System.out.println(BaseEncoder.getSequenceFromLong(chr0CtgRBT.getRead(i)));
            System.out.println("Index in RBT:" + chr0CtgRBT.getReadIndex(chr0CtgRBT.getRead(i)));
        }
        MapHaplotypes.mapWorkingMapWithGenetic(gMapGBS,phyMappedReadsIBM,chr0CtgRBT,0.2,10,10,chr0CtgHaplotypesResultsFile);
    }

    public static void performPstIVD() {
        String AGPv2Fastafile =  "S:/SolexaAnal/AGPv2/maize_agp_v2.fasta";
        String AGPv2PstIVDfile = "H:/VirtualDigests/B73/B73AGPv2PstIVD.rwpm.bin";

        VirtualDigester myVDigester=new VirtualDigester(new File(AGPv2Fastafile), new File(AGPv2PstIVDfile));
        System.gc();

        //Sort the physical map for quick searching and save
        ReadsWPhysicalMap rwpmMyVD = new ReadsWPhysicalMap(AGPv2PstIVDfile, true);
        rwpmMyVD.sortTable(true);
        rwpmMyVD.writeCountFile(new File(AGPv2PstIVDfile));
        System.gc();
    }

    public static void makeReadCountsFileForIBMPstI() {
        String qseqDir =           "H:/64GJAAAXX/qseq";
        String qseqKey =           "H:/64GJAAAXX/64GJAAAXX_PstI_IBM94_key.txt";
        String parsedTaxaDir =     "H:/64GJAAAXX/taxa";
        String collapsedTaxaDir =  "H:/64GJAAAXX/taxacollapse";
        String combinedReadsFile = "H:/64GJAAAXX/counts/CombReads_IBM_PstI_min50_20110708.bin";
        String rbtFile =           "H:/64GJAAAXX/counts/rbtIBM_PstI_Min50_20110708.bin";

//        // Read and parse qseq or fastq file using low stringency Q score (maximal # of reads)
//        System.out.println();
//        System.out.println("Parsing the qseq files for IBM PstI (64GJAAAXX)");
//        System.out.println("-----------------------------------------------");
//        ParseBarcodeFiles pbf=new ParseBarcodeFiles(qseqDir, qseqKey, parsedTaxaDir, 0, false);
//        pbf = null;
//        System.gc();
//
//        // condense the taxa ReadCounts files by sorting and collapsing identical reads
//        System.out.println();
//        System.out.println("Collapsing the ReadCounts files for each taxon in IBM PstI (64GJAAAXX)");
//        System.out.println("----------------------------------------------------------------------");
//        ReadCounts rc=new ReadCounts(parsedTaxaDir, collapsedTaxaDir, 1, true, false, false);
//        rc = null;
//        System.gc();

        // Fuse multiple taxa files into one, minimum frequency is 50, binary=true (good PstI tags should be quite common)
        System.out.println();
        System.out.println("Combining the ReadCounts files across all taxa in IBM PstI (64GJAAAXX)");
        System.out.println("----------------------------------------------------------------------");
        CombineReadCounts crc=new CombineReadCounts(collapsedTaxaDir, combinedReadsFile, 50, true);
        crc = null;
        System.gc();

        CreateReadsByTaxa rbt=new CreateReadsByTaxa(combinedReadsFile, collapsedTaxaDir, rbtFile, true);
    }

    public static void makePhysicalMapForIBMPstI() {
        String AGPv2PstIVDfile =    "H:/VirtualDigests/B73/B73AGPv2PstIVD.rwpm.bin";
        String combinedReadsFile =  "H:/64GJAAAXX/counts/CombReads_IBM_PstI_min50_20110708.bin";
        String phyMappedReadsFile = "H:/64GJAAAXX/counts/phyMappedReadsIBMPstI_min50_Div20_20110708.bin";

        ReadsWPhysicalMap rwpm=MapHaplotypes.blastUsingPhysicalMap(new ReadsWPhysicalMap(AGPv2PstIVDfile,true),
                  new ReadsWPhysicalMap(new ReadCounts(combinedReadsFile, true)), 20);  // allow 20 mismatches
        rwpm.sortTable(false); // sort by position
        rwpm.writeCountFile(new File(phyMappedReadsFile), Integer.MAX_VALUE, false, false, Float.NaN, true);
    }

    public static void mapIBMPstITestVs55KFrame() {
//        String frameWorkGenMapFile55K =   "/usr/local/maizediv/illumina/IBM/counts/IBM_55K_08192010.hmp.txt";
        String frameWorkGenMapFile55K =   "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/de_novo_ctgs/IBM_55K_08192010.hmp.txt";
//        String frameWorkGenMapFileGBS =   "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/434LFAAXX/anchor_map/IBM_GBS_framework_map_20100917.hmp.txt";

        String rbtFile =            "H:/64GJAAAXX/counts/rbtIBM_PstI_Min50_20110708.bin";
        String phyMappedReadsFile = "H:/64GJAAAXX/counts/phyMappedReadsIBMPstI_min50_Div20_20110708.bin";
        String GeneticResultsFile = "H:/64GJAAAXX/counts/IBMPstIGeneticResults_min50_r10_t15_20110708.txt";

        // genetically map every tag in IBM PstI that showed up 100 or more times
        Alignment[] gMap55K= ((CombineAlignment)ImportUtils.readFromHapmap(frameWorkGenMapFile55K)).getAlignments();
        MapHaplotypes.mapWorkingMapWithGenetic(gMap55K, new ReadsWPhysicalMap(phyMappedReadsFile,true),
                new ReadsByTaxa(rbtFile,true), 0.1, 15, 10, GeneticResultsFile);
    }

    public static void imputeNAMSNPFrameworkMap() {
        String frameWorkGenMapFile = "C:/Users/jcg233/Documents/Bioinformatics/Panzea website/NAM/NAM genos/published_on_panzea/Update_20080703/HapMapFormat/NAMgenos_1106SNP_AGPv2_20080103.hmp.txt";
        String HapMapFileStem =      "H:/70MU0AAXX/NAMgenos_1106SNP_AGPv2_imputed";

        Alignment[] NAMgMap1106SNPs =  ((CombineAlignment)ImportUtils.readFromHapmap(frameWorkGenMapFile)).getAlignments();
        for (int i = 1; i < 11; ++i) {
            MutableAlignmentForGBS  mutAlign = new MutableAlignmentForGBS(NAMgMap1106SNPs[i-1]);
            mutAlign.imputeMissingDataIncludingHets();
            ExportUtils.writeToHapmap(mutAlign, false, HapMapFileStem + "_chr" + i, '\t');
        }
    }

    public static void makeReadCountsFileForNAMPstI() {
        String qseqDir =           "H:/70MU0AAXX/qseq";
        String qseqKey =           "H:/70MU0AAXX/70MU0AAXX_NAM_PstI_key.txt";
        String parsedTaxaDir =     "H:/70MU0AAXX/taxa";
        String collapsedTaxaDir =  "H:/70MU0AAXX/taxacollapse";
        String combinedReadsFile = "H:/70MU0AAXX/counts/CombReads_NAM_PstI_min50_118taxa_20110714.bin";
        String rbtFile =           "H:/70MU0AAXX/counts/rbtNAM_PstI_Min50_118taxa_20110714.bin";

//        // Read and parse qseq or fastq file using low stringency Q score (maximal # of reads)
//        System.out.println();
//        System.out.println("Parsing the qseq files for NAM PstI (70MU0AAXX)");
//        System.out.println("-----------------------------------------------");
//        ParseBarcodeFiles pbf=new ParseBarcodeFiles(qseqDir, qseqKey, parsedTaxaDir, 0, false);
//        pbf = null;
//        System.gc();
//
//        // condense the taxa ReadCounts files by sorting and collapsing identical reads
//        System.out.println();
//        System.out.println("Collapsing the ReadCounts files for each taxon in NAM PstI (70MU0AAXX)");
//        System.out.println("----------------------------------------------------------------------");
//        ReadCounts rc=new ReadCounts(parsedTaxaDir, collapsedTaxaDir, 1, true, false, false);
//        rc = null;
//        System.gc();

        // Fuse multiple taxa files into one, minimum frequency is 50, binary=true (good PstI tags should be quite common)
        System.out.println();
        System.out.println("Combining the ReadCounts files across all taxa in NAM PstI (70MU0AAXX)");
        System.out.println("----------------------------------------------------------------------");
        CombineReadCounts crc=new CombineReadCounts(collapsedTaxaDir, combinedReadsFile, 50, true);
        crc = null;
        System.gc();

        CreateReadsByTaxa rbt=new CreateReadsByTaxa(combinedReadsFile, collapsedTaxaDir, rbtFile, true);
    }

    public static void makePhysicalMapForNAMPstI() {
        String AGPv2PstIVDfile =    "H:/VirtualDigests/B73/B73AGPv2PstIVD.rwpm.bin";
        String combinedReadsFile =  "H:/70MU0AAXX/counts/CombReads_NAM_PstI_min50_118taxa_20110714.bin";
        String phyMappedReadsFile = "H:/70MU0AAXX/counts/phyMappedReadsNAMPstI_min50_118taxa_Div20_20110714.bin";

        ReadsWPhysicalMap rwpm=MapHaplotypes.blastUsingPhysicalMap(new ReadsWPhysicalMap(AGPv2PstIVDfile,true),
                  new ReadsWPhysicalMap(new ReadCounts(combinedReadsFile, true)), 20);  // allow 20 mismatches
        rwpm.sortTable(false); // sort by position
        rwpm.writeCountFile(new File(phyMappedReadsFile), Integer.MAX_VALUE, false, false, Float.NaN, true);
    }

    public static void mapNAMPstITestVsSNPFrame() {
        String frameWorkGenMapFile = "H:/NAM1106SNPframe/NAMgenos_1106SNP_AGPv2_imputed.hmp.txt";
        String rbtFile =             "H:/70MU0AAXX/counts/rbtNAM_PstI_Min50_118taxa_20110714.bin";
        String phyMappedReadsFile =  "H:/70MU0AAXX/counts/phyMappedReadsNAMPstI_min50_118taxa_Div20_20110714.bin";
        String GeneticResultsFile =  "H:/70MU0AAXX/counts/NAMPstIGeneticResults_min50_118taxa_r10_t15_20110714.txt";

        // genetically map every tag in NAM PstI that showed up 50 or more times
        Alignment[] gMap55K= ((CombineAlignment)ImportUtils.readFromHapmap(frameWorkGenMapFile)).getAlignments();
        MapHaplotypes.mapWorkingMapWithGenetic(gMap55K, new ReadsWPhysicalMap(phyMappedReadsFile,true),
                new ReadsByTaxa(rbtFile,true), 0.1, 15, 10, GeneticResultsFile);
    }

    public static void makeReadCountsFilesForNAMApeKI() {
        String qseqDir =           "H:/NAM_ApeKI_plates49_50/qseq";
        String qseqKey =           "H:/NAM_ApeKI_plates49_50/NAM49_50_ApeKI_key.txt";
        String parsedTaxaDir =     "H:/NAM_ApeKI_plates49_50/taxa";
        String collapsedTaxaDir =  "H:/NAM_ApeKI_plates49_50/taxacollapse";

        // Read and parse qseq or fastq file using low stringency Q score (maximal # of reads)
        System.out.println();
        System.out.println("Parsing the qseq files");
        System.out.println("----------------------");
        ParseBarcodeFiles pbf=new ParseBarcodeFiles(qseqDir, qseqKey, parsedTaxaDir, 0, false);
        pbf = null;
        System.gc();

        // condense the taxa ReadCounts files by sorting and collapsing identical reads
        System.out.println();
        System.out.println("Collapsing the ReadCounts files for each taxon");
        System.out.println("----------------------------------------------");
        ReadCounts rc=new ReadCounts(parsedTaxaDir, collapsedTaxaDir, 1, true, false, false);
        rc = null;
        System.gc();
    }

    public static void makeTagsByTaxaForNAMApeKI() {
        String collapsedTaxaDir =  "H:/NAM_ApeKI_plates49_50/taxacollapse";
        String combinedReadsFile = "H:/NAM_ApeKI_plates49_50/counts/CombReads_NAM49_50_ApeKI_min10_152taxa_20110714.bin";
        String rbtFile =           "H:/NAM_ApeKI_plates49_50/counts/rbtNAM49_50_ApeKI_Min10_152taxa_20110714.bin";

        // Fuse multiple taxa files into one, minimum frequency is 10, binary=true
        System.out.println();
        System.out.println("Combining the ReadCounts files across all taxa");
        System.out.println("----------------------------------------------");
        CombineReadCounts crc=new CombineReadCounts(collapsedTaxaDir, combinedReadsFile, 10, true);
        crc = null;
        System.gc();

        CreateReadsByTaxa rbt=new CreateReadsByTaxa(combinedReadsFile, collapsedTaxaDir, rbtFile, true);
    }

    public static void makePhysicalMapForNAMApeKI() {
        String AGPv2ApeKIVDfile =   "S:/SolexaAnal/virtual_digests/AGPv2/B73AGPv2VD.rwpm.bin";
        String combinedReadsFile =  "H:/NAM_ApeKI_plates49_50/counts/CombReads_NAM49_50_ApeKI_min10_152taxa_20110714.bin";
        String phyMappedReadsFile = "H:/NAM_ApeKI_plates49_50/counts/phyMappedReadsNAMApeKI_min10_152taxa_div20_20110714.bin";

        ReadsWPhysicalMap rwpm=MapHaplotypes.blastUsingPhysicalMap(new ReadsWPhysicalMap(AGPv2ApeKIVDfile,true),
                  new ReadsWPhysicalMap(new ReadCounts(combinedReadsFile, true)), 20);  // allow 20 mismatches
        rwpm.sortTable(false); // sort by position
        rwpm.writeCountFile(new File(phyMappedReadsFile), Integer.MAX_VALUE, false, false, Float.NaN, true);
    }

    public static void mapNAMApeKIVsSNPFrame() {
        String frameWorkGenMapFile = "H:/NAM1106SNPframe/NAMgenos_1106SNP_AGPv2_imputed.hmp.txt";
        String rbtFile =             "H:/NAM_ApeKI_plates49_50/counts/rbtNAM49_50_ApeKI_Min10_152taxa_20110714.bin";
        String phyMappedReadsFile =  "H:/NAM_ApeKI_plates49_50/counts/phyMappedReadsNAMApeKI_min10_152taxa_div20_20110714.bin";
        String GeneticResultsFile =  "H:/NAM_ApeKI_plates49_50/counts/NAMApeKIGeneticResults_min10_152taxa_r10_t15_20110714.txt";

        // genetically map every tag in NAM ApeKI that showed up in 15 or more taxa (and sample size of recomb test >= 15)
        Alignment[] gMap55K= ((CombineAlignment)ImportUtils.readFromHapmap(frameWorkGenMapFile)).getAlignments();
        MapHaplotypes.mapWorkingMapWithGenetic(gMap55K, new ReadsWPhysicalMap(phyMappedReadsFile,true),
                new ReadsByTaxa(rbtFile,true), 0.1, 15, 10, GeneticResultsFile);
    }

    public static void makeReadCountsFilesForCintaFengFineMapping() {
//        String qseqDir =           "/cbsufsrv4/data1/maizediv/illumina/708CLAAXX/qseq/";
//        String qseqKey =           "/cbsufsrv4/data1/maizediv/illumina/708CLAAXX/708CLAAXX_FineMapping_key.txt";
//        String parsedTaxaDir =     "/cbsufsrv4/data1/maizediv/illumina/qual_control/708CLAAXX/fineMapping/taxa";
//        String collapsedTaxaDir =  "/cbsufsrv4/data1/maizediv/illumina/qual_control/708CLAAXX/fineMapping/taxacollapse";

        String qseqDir =           "/cbsufsrv4/data1/maizediv/illumina/C00TDABXX/qseq/";
        String qseqKey =           "/cbsufsrv4/data1/maizediv/illumina/C00TDABXX/C00TDABXX_FineMapping_key.txt";
        String parsedTaxaDir =     "/cbsufsrv4/data1/maizediv/illumina/qual_control/708CLAAXX/fineMapping/C00TDABXX/taxa";
        String collapsedTaxaDir =  "/cbsufsrv4/data1/maizediv/illumina/qual_control/708CLAAXX/fineMapping/C00TDABXX/taxacollapse";

        // Read and parse qseq or fastq file using low stringency Q score (maximal # of reads)
        System.out.println();
        System.out.println("Parsing the qseq files");
        System.out.println("----------------------");
        ParseBarcodeFiles pbf=new ParseBarcodeFiles(qseqDir, qseqKey, parsedTaxaDir, 0, false);
        pbf = null;
        System.gc();

        // condense the taxa ReadCounts files by sorting and collapsing identical reads
        System.out.println();
        System.out.println("Collapsing the ReadCounts files for each taxon");
        System.out.println("----------------------------------------------");
        ReadCounts rc=new ReadCounts(parsedTaxaDir, collapsedTaxaDir, 1, true, false, false);
        rc = null;
        System.gc();
    }


    public static void makeReadCountsTeosinteTestPlate() {
        String qseqDirTeo =           "/usr/local/maizediv/illumina/teosinte/qseq";
        String qseqKeyTeo =           "/usr/local/maizediv/illumina/teosinte/teosinte_key.txt";
        String parsedTaxaDirTeo =     "/usr/local/maizediv/illumina/teosinte/taxa";
        String collapsedTaxaDirTeo =  "/usr/local/maizediv/illumina/teosinte/taxacollapse";
        String combinedReadsFileTeo = "/usr/local/maizediv/illumina/teosinte/counts/CombReadsTeoMin2_20100712.bin";
        String rbtFileTeo =           "/usr/local/maizediv/illumina/teosinte/counts/teoRBT_Min2_20100716.bin";
        String rbtTextFileTeo =       "/usr/local/maizediv/illumina/teosinte/counts/teoRBT_Min2_20100719.txt";

        // Read and parse qseq or fastq file using low stringency Q score (maximal # of reads)
        System.out.println();
        System.out.println("Parsing the qseq files for Teosinte");
        System.out.println("-----------------------------------");
        ParseBarcodeFiles pbfTeo=new ParseBarcodeFiles(qseqDirTeo, qseqKeyTeo, parsedTaxaDirTeo, 0, false);
        System.gc();

        // condense the taxa ReadCounts files by sorting and collapsing identical reads
        System.out.println();
        System.out.println("Collapsing the ReadCounts files for each teosinte taxon");
        System.out.println("-------------------------------------------------------");
        ReadCounts rcTeo=new ReadCounts(parsedTaxaDirTeo, collapsedTaxaDirTeo, 1, true, false, false);
        rcTeo = null;
        System.gc();

        // Fuse multiple taxa files into one, minimum frequency is 2, binary=true (low min freq b/c BC2S3)
        System.out.println();
        System.out.println("Combining the ReadCounts files across all teosinte taxa");
        System.out.println("-------------------------------------------------------");
        CombineReadCounts crcTeo=new CombineReadCounts(collapsedTaxaDirTeo, combinedReadsFileTeo, 2, true);
        crcTeo = null;
        System.gc();

        CreateReadsByTaxa teoRBT=new CreateReadsByTaxa(combinedReadsFileTeo, collapsedTaxaDirTeo, rbtFileTeo, true);
    }

    public static void makeTeosinteTestPlateGenoFile() {
        // input files
        String teoFrameworkMap = "Q:/SolexaAnal/Solexa_GBS/teosinte/Teo_BC2S3_framework_HapMap.txt";
        String AGPv2VDfile = "Q:/SolexaAnal/virtual_digests/AGPv2/B73AGPv2VD.rwpm.bin";  // although JD's framework is v1, AGPv2 should suffice
        String teoCombinedReadsFile = "Q:/SolexaAnal/Solexa_GBS/teosinte/CombReadsTeoMin2_20100712.bin";
        String teoRBTFile = "Q:/SolexaAnal/Solexa_GBS/teosinte/teoRBT_Min2_20100716.bin";

        // output files
        String teoPhyMappedReadsFile = "Q:/SolexaAnal/Solexa_GBS/teosinte/20100819/phyMappedReadsTeoDiv5_20100716.bin";
        String teoPhyMappedReadsTxtFile = "Q:/SolexaAnal/Solexa_GBS/teosinte/20100819/phyMappedReadsTeoDiv5_20100716.txt";
        String teoGenMappedReadsTxtFile = "Q:/SolexaAnal/Solexa_GBS/teosinte/20100819/phyGenMappedReadsTeoDiv5_20100719.txt";
        String teoGenMappedReadsFile = "Q:/SolexaAnal/Solexa_GBS/teosinte/20100819/phyGenMappedReadsTeoDiv5_20100720.bin";
        String teoGBSHapMapFileStem = "Q:/SolexaAnal/Solexa_GBS/teosinte/20100819/teo_GBS_biallelic_genos_20100730";

        ReadsWPhysicalMap rwpmTeoPhyMapped=MapHaplotypes.blastUsingPhysicalMap(new ReadsWPhysicalMap(AGPv2VDfile,true),
                  new ReadsWPhysicalMap(new ReadCounts(teoCombinedReadsFile, true)),5);  // allow 5 mismatches for a 1st pass
        rwpmTeoPhyMapped.writeCountFile(new File(teoPhyMappedReadsFile), Integer.MAX_VALUE, true, false, Float.NaN, true);
        rwpmTeoPhyMapped.sortTable(false); // sort by position

        ReadsByTaxa teoRBT=new ReadsByTaxa(teoRBTFile,true);
        rwpmTeoPhyMapped=MapHaplotypes.checkWorkingMapWithGeneticBC2(teoFrameworkMap, rwpmTeoPhyMapped, teoRBT, 0.05);
//        rwpmTeoPhyMapped.writeCountFile(new File(teoGenMappedReadsTxtFile), Integer.MAX_VALUE, true, false, Float.NaN, false);
        rwpmTeoPhyMapped.writeCountFile(new File(teoGenMappedReadsFile), Integer.MAX_VALUE, true, false, Float.NaN, true);
        ReadsWPhysicalMap rwpmTeoGenPhyMapped = new ReadsWPhysicalMap(teoGenMappedReadsFile,true);

        for(int i=1; i<=10; i++) {
            System.out.println("WORKING ON CHROMOSOME " + i);
            Alignment a = new ConvertGBSToAlignment(teoRBT, rwpmTeoGenPhyMapped, 0.01, 0.02, (byte)i);
            System.out.println("UnFiltered Alignment  Taxa:" + a.getSequenceCount() + " Sites:" + a.getSiteCount());
            int[] goodLowHetSites = AlignmentFilterByGBSUtils.getLowHetSNPs(a, 0.20, 5);
            a = FilterAlignment.getInstance(a, goodLowHetSites);
            MutableAlignmentForGBS mutAlign=new MutableAlignmentForGBS(a);
            a = null;
            System.gc();
            mutAlign.writeToHapmapBC2S3(false, false, teoGBSHapMapFileStem + ".chr" + i + ".CAfiltered", '\t', i);  // writes only the C/A markers to the file
            Alignment teoW22BC2S3 = ImportUtils.readFromHapmap(teoGBSHapMapFileStem + ".chr" + i + ".CAfiltered" + ".hmp.txt");
            System.out.println("Filtered Alignment  Taxa:" + teoW22BC2S3.getSequenceCount() + " Sites:" + teoW22BC2S3.getSiteCount());
            mutAlign = new MutableAlignmentForGBS(teoW22BC2S3);
            teoW22BC2S3 = null;
            mutAlign.removeDCOs(7);
            System.out.println("nDCOs on chr" + i + " = " + mutAlign.getnDCOs());
            mutAlign.callHetSegments(7);
            System.out.println("nDCOs on chr" + i + " = " + mutAlign.getnDCOs());
            mutAlign.writeToHapmapBC2S3(false, true, teoGBSHapMapFileStem + ".chr" + i + ".hetSegs", '\t', i);
            teoW22BC2S3 = ImportUtils.readFromHapmap(teoGBSHapMapFileStem + ".chr" + i + ".hetSegs" + ".hmp.txt");
            mutAlign = new MutableAlignmentForGBS(teoW22BC2S3);
            mutAlign.imputeMissingDataIncludingHets();
            mutAlign.writeToHapmapBC2S3(false, false, teoGBSHapMapFileStem + ".chr" + i + ".imputed", '\t', i);
        }
        for (int base = 0; base < 16; ++base) {
            System.out.println("GdpdmBLOBUtils.bases[" + base + "] is: " + (char) GdpdmBLOBUtils.bases[base]);
        }
    }

    public static void makeReadCountsTeosinte() {
        String qseqDirTeo =           "/usr/local/maizediv/illumina/teosinte/qseq";
        String qseqKeyTeo =           "/usr/local/maizediv/illumina/teosinte/all_teoW22BC2S3_MR_key_20110401.txt";  // the .1, .2 and .5 suffixes removed to allow merging of genetically equiv. samples
        String parsedTaxaDirTeo =     "/usr/local/maizediv/illumina/teosinte/taxa";
        String collapsedTaxaDirTeo =  "/usr/local/maizediv/illumina/teosinte/taxacollapse";
        String combinedReadsFileTeo = "/usr/local/maizediv/illumina/teosinte/counts/CombReadsTeoMin2_20110401.bin";
        String rbtFileTeo =           "/usr/local/maizediv/illumina/teosinte/counts/teoRBT_Min3_20110412.bin";
        String rbtTextFileTeo =       "/usr/local/maizediv/illumina/teosinte/counts/teoRBT_Min2_20110402.txt";

        // Read and parse qseq or fastq file using low stringency Q score (maximal # of reads)
        System.out.println();
        System.out.println("Parsing the qseq files for Teosinte");
        System.out.println("-----------------------------------");
        ParseBarcodeFiles pbfTeo=new ParseBarcodeFiles(qseqDirTeo, qseqKeyTeo, parsedTaxaDirTeo, 0, false);
        System.gc();

        // condense the taxa ReadCounts files by sorting and collapsing identical reads
        System.out.println();
        System.out.println("Collapsing the ReadCounts files for each teosinte taxon");
        System.out.println("-------------------------------------------------------");
        ReadCounts rcTeo=new ReadCounts(parsedTaxaDirTeo, collapsedTaxaDirTeo, 1, true, false, true);  // final true = combine identical taxa
        rcTeo = null;
        System.gc();

        // Fuse multiple taxa files into one, minimum frequency is 2, binary=true (low min freq b/c BC2S3)
        System.out.println();
        System.out.println("Combining the ReadCounts files across all teosinte taxa");
        System.out.println("-------------------------------------------------------");
        CombineReadCounts crcTeo=new CombineReadCounts(collapsedTaxaDirTeo, combinedReadsFileTeo, 2, true);
        crcTeo = null;
        System.gc();

        CreateReadsByTaxa teoRBT=new CreateReadsByTaxa(combinedReadsFileTeo, collapsedTaxaDirTeo, rbtFileTeo, 3, true);
    }


    /**
     * Makes ReadCounts for and readsByTaxa including duplicate teosinte samples.
     *
     * This will allow a reproducibility check.  Genetically equivalent taxa,
     * derived from the same RIL (but different generations), should yield the
     * same GBS genotypes.
     *
     * @return            null
     */
    public static void makeReadCountsTeosinteDuplicates() {
        String qseqDirTeo =           "/usr/local/maizediv/illumina/teosinte/qseq";
        String qseqKeyTeoDuplicates = "/usr/local/maizediv/illumina/teosinte/all_teoW22BC2S3_duplicates_key_20110401.txt";  // the .1, .2 and .5 suffixes included to distinguish genetically equiv. samples
        String parsedTaxaDirTeoDuplicates =    "/usr/local/maizediv/illumina/teosinte/duplicatetaxa";
        String collapsedTaxaDirTeoDuplicates = "/usr/local/maizediv/illumina/teosinte/duplicatetaxacollapse";
        String collapsedTaxaDirTeo =  "/usr/local/maizediv/illumina/teosinte/taxacollapse";
        String combinedReadsFileTeo = "/usr/local/maizediv/illumina/teosinte/counts/CombReadsTeoMin2_20110401.bin";  // use the same one as before
        String rbtFileTeoDuplicates = "/usr/local/maizediv/illumina/teosinte/counts/teoRBT_withDuplicates_Min3_20110503.bin";

        // Read and parse qseq or fastq file using low stringency Q score (maximal # of reads)
        System.out.println();
        System.out.println("Parsing the qseq files for Teosinte duplicates");
        System.out.println("----------------------------------------------");
        ParseBarcodeFiles pbfTeo=new ParseBarcodeFiles(qseqDirTeo, qseqKeyTeoDuplicates, parsedTaxaDirTeoDuplicates, 0, false);
        System.gc();

        // condense the taxa ReadCounts files by sorting and collapsing identical reads
        System.out.println();
        System.out.println("Collapsing the ReadCounts files for each duplicate teosinte taxon");
        System.out.println("-----------------------------------------------------------------");
        ReadCounts rcTeoDup=new ReadCounts(parsedTaxaDirTeoDuplicates, collapsedTaxaDirTeoDuplicates, 1, true, false, false);  // final false = do not combine identical taxa (irrelevant)
        rcTeoDup = null;
        System.gc();

        // need to (manually) move the duplicate taxa ReadCounts from collapsedTaxaDirTeoDuplicates to collapsedTaxaDirTeo before running this step
        CreateReadsByTaxa teoRBT=new CreateReadsByTaxa(combinedReadsFileTeo, collapsedTaxaDirTeo, rbtFileTeoDuplicates, 3, true);
    }


    /**
     * Makes ReadCounts for and readsByTaxa including the Doebley lab's fine 
     * mapping sample W22QTL10.
     *
     * @return            null
     */
    public static void makeReadCountsTeoW22QTL10() {
        String qseqDirW22QTL10 =          "/usr/local/maizediv/illumina/teosinte/qseq/W22QTL10";
        String qseqKeyW22QTL10 =          "/usr/local/maizediv/illumina/teosinte/C00R8ABXX_W22QTL10_key.txt";  // key contains all samples on lane (so delete them)
        String parsedTaxaDirW22QTL10 =    "/usr/local/maizediv/illumina/teosinte/W22QTL10";
        String collapsedTaxaDirW22QTL10 = "/usr/local/maizediv/illumina/teosinte/W22QTL10collapse";
        String collapsedTaxaDirTeo =      "/usr/local/maizediv/illumina/teosinte/taxacollapse";
        String combinedReadsFileTeo =     "/usr/local/maizediv/illumina/teosinte/counts/CombReadsTeoMin2_20110401.bin";  // use the same one as before
        String rbtFileTeoWithW22QTL10 =   "/usr/local/maizediv/illumina/teosinte/counts/teoRBT_withW22QTL10_Min3_20110610.bin";

        // Read and parse qseq or fastq file using low stringency Q score (maximal # of reads)
        System.out.println();
        System.out.println("Parsing the qseq file containing W22QTL10");
        System.out.println("-----------------------------------------");
        ParseBarcodeFiles pbfTeo=new ParseBarcodeFiles(qseqDirW22QTL10, qseqKeyW22QTL10, parsedTaxaDirW22QTL10, 0, false);
        System.gc();

        // condense the taxa ReadCounts files by sorting and collapsing identical reads (delete the useless samples first)
//        System.out.println();
//        System.out.println("Collapsing the ReadCounts files for W22QTL10");
//        System.out.println("--------------------------------------------");
//        ReadCounts rcTeoDup=new ReadCounts(parsedTaxaDirW22QTL10, collapsedTaxaDirW22QTL10, 1, true, false, false);  // final false = do not combine identical taxa (irrelevant)
//        rcTeoDup = null;
//        System.gc();

        // need to (manually) move the W22QTL10 collapsed ReadCounts to collapsedTaxaDirTeo before running this step
//        CreateReadsByTaxa teoRBT=new CreateReadsByTaxa(combinedReadsFileTeo, collapsedTaxaDirTeo, rbtFileTeoWithW22QTL10, 3, true);
    }

    public static void makeTeosinteGenoFile() {
        // input files
        String teoFrameworkMap = "/usr/local/maizediv/illumina/teosinte/Teo_BC2S3_885RILs_framework.hmp.txt";
        String teoFrameworkGBSMap = "/usr/local/maizediv/illumina/teosinte/teo_GBS_biallelic_genos_20110422.all.imputed.hmp.txt";
        String AGPv2VDfile = "/usr/local/maizediv/genome/virtual_digest/B73refVirtualV2.rwpm.bin";  
        String teoCombinedReadsFile = "/usr/local/maizediv/illumina/teosinte/counts/CombReadsTeoMin2_20110401.bin";
        String teoRBTFile = "/usr/local/maizediv/illumina/teosinte/counts/teoRBT_Min3_20110412.bin";
        String rbtFileTeoDuplicates = "/usr/local/maizediv/illumina/teosinte/counts/teoRBT_withDuplicates_Min3_20110503.bin";
        String rbtFileTeoWithW22QTL10 =   "/usr/local/maizediv/illumina/teosinte/counts/teoRBT_withW22QTL10_Min3_20110610.bin";
        String teoPhyMappedReadsFile = "/usr/local/maizediv/illumina/teosinte/counts/phyMappedReadsTeoDiv5_20110404.bin";

        // output files
        String teoGenMappedReadsFile = "/usr/local/maizediv/illumina/teosinte/phyGenMappedReadsTeoDiv5r10_20110423.bin";
        String teoGBSHapMapFileStem = "/usr/local/maizediv/illumina/teosinte/20110610/teo_GBS_biallelic_genos_20110610";

        ReadsWPhysicalMap rwpmTeoPhyMapped=MapHaplotypes.blastUsingPhysicalMap(new ReadsWPhysicalMap(AGPv2VDfile,true),
                  new ReadsWPhysicalMap(new ReadCounts(teoCombinedReadsFile, true)),5);  // allow 5 mismatches for a 1st pass
        rwpmTeoPhyMapped.writeCountFile(new File(teoPhyMappedReadsFile), Integer.MAX_VALUE, true, false, Float.NaN, true);
//        ReadsWPhysicalMap rwpmTeoPhyMapped = new ReadsWPhysicalMap(teoPhyMappedReadsFile,true);
        rwpmTeoPhyMapped.sortTable(false); // sort by position

        ReadsByTaxa teoRBT=new ReadsByTaxa(rbtFileTeoWithW22QTL10,true);
        rwpmTeoPhyMapped=MapHaplotypes.checkWorkingMapWithGeneticBC2(teoFrameworkGBSMap, rwpmTeoPhyMapped, teoRBT, 0.1);
        rwpmTeoPhyMapped.writeCountFile(new File(teoGenMappedReadsFile), Integer.MAX_VALUE, true, false, Float.NaN, true);
        rwpmTeoPhyMapped = null;
        System.gc();
        ReadsWPhysicalMap rwpmTeoGenPhyMapped = new ReadsWPhysicalMap(teoGenMappedReadsFile,true);

        for(int i=1; i<=10; i++) {
            System.out.println("WORKING ON CHROMOSOME " + i);
            Alignment a = new ConvertGBSToAlignment(teoRBT, rwpmTeoGenPhyMapped, 0.1, 0.02, (byte)i);
            System.out.println("UnFiltered Alignment  Taxa:" + a.getSequenceCount() + " Sites:" + a.getSiteCount());
            int[] goodLowHetSites = AlignmentFilterByGBSUtils.getLowHetSNPs(a, 0.20, 5);
            a = FilterAlignment.getInstance(a, goodLowHetSites);
            MutableAlignmentForGBS mutAlign=new MutableAlignmentForGBS(a);
//            System.out.println("Filtered Alignment for Excessive Hets:   Taxa:" + mutAlign.getSequenceCount() + " Sites:" + mutAlign.getSiteCount());
            a = null;
            System.gc();
            mutAlign.writeToHapmapBC2S3CA(false, false, teoGBSHapMapFileStem + ".chr" + i + ".CAfiltered", '\t', i);  // writes only the C/A and C/A/M markers to the file
            Alignment teoW22BC2S3 = ImportUtils.readFromHapmap(teoGBSHapMapFileStem + ".chr" + i + ".CAfiltered" + ".hmp.txt");
            System.out.println("Filtered Alignment for C/A:    Taxa:" + teoW22BC2S3.getSequenceCount() + " Sites:" + teoW22BC2S3.getSiteCount());
            mutAlign = new MutableAlignmentForGBS(teoW22BC2S3);
            teoW22BC2S3 = null;
            mutAlign.removeDCOs(7);
            System.out.println("nDCOs on chr" + i + " = " + mutAlign.getnDCOs());
            mutAlign.callHetSegments(7);
            System.out.println("nDCOs on chr" + i + " = " + mutAlign.getnDCOs());  // it's possible for callHetSegments to flag more DCOs
            mutAlign.writeToHapmapBC2S3(false, true, teoGBSHapMapFileStem + ".chr" + i + ".hetSegs", '\t', i);
            teoW22BC2S3 = ImportUtils.readFromHapmap(teoGBSHapMapFileStem + ".chr" + i + ".hetSegs" + ".hmp.txt");
            mutAlign = new MutableAlignmentForGBS(teoW22BC2S3);
            mutAlign.imputeMissingDataIncludingHets();
            mutAlign.writeToHapmapBC2S3(false, false, teoGBSHapMapFileStem + ".chr" + i + ".imputed", '\t', i);
        }
        for (int base = 0; base < GdpdmBLOBUtils.bases.length; ++base) {
            System.out.println("GdpdmBLOBUtils.bases[" + base + "] is: " + (char) GdpdmBLOBUtils.bases[base]);
        }
    }

    public static void writeTeosinteAlleleInfoFile() {
        String AGPv2VDfile = "/usr/local/maizediv/genome/virtual_digest/B73refVirtualV2.rwpm.bin";
        String teoRBTFile = "/usr/local/maizediv/illumina/teosinte/counts/teoRBT_Min3_20110412.bin";
        String teoGenMappedReadsFile = "/usr/local/maizediv/illumina/teosinte/phyGenMappedReadsTeoDiv5r10_20110423.bin";
        String teoGBSHapMapFileStem = "/usr/local/maizediv/illumina/teosinte/20110423/teo_GBS_biallelic_genos_20110423";
        String AGPv1VDfile = "/usr/local/maizediv/genome/virtual_digest/B73refVirtualV1.rwpm.bin";
        String teoAlleleInfoInFile = "/usr/local/maizediv/illumina/teosinte/20110423/teo_GBS_biallelic_genos_20110426.alleleInfo.txt";
        String teoAlleleInfoOutFile = "/usr/local/maizediv/illumina/teosinte/20110423/teo_GBS_biallelic_genos_20110426.alleleInfo.withAGPv1.txt";

        // Chr0 stuff (to help figure out v2 chr0 contig names based on v1 positions)
        String chr0AlleleInfoInFile = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/434LFAAXX/IBM_GBSwQualReads_100513b_chr0.txt";
        String chr0AlleleInfoOutFile = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/434LFAAXX/IBM_GBSwQualReads_100513b_chr0.withAGPv1.M20.txt";

//        ReadsByTaxa teoRBT=new ReadsByTaxa(teoRBTFile,true);
//        ReadsWPhysicalMap rwpmTeoGenPhyMapped = new ReadsWPhysicalMap(teoGenMappedReadsFile,true);
//        ReadsWPhysicalMap AGPv2VD = new ReadsWPhysicalMap(AGPv2VDfile,true);
//        AGPv2VD.sortTable(false);  // sort by position
//
//        for(int i=1; i<=10; i++) {
//            System.out.println("WORKING ON CHROMOSOME " + i);
//            Alignment a = new ConvertGBSToAlignment(teoRBT, AGPv2VD,
//                    rwpmTeoGenPhyMapped, 0.1, 0.02, (byte)i, teoGBSHapMapFileStem+".chr"+i+".alleleInfo.txt");
//        }

        AlleleInfo ai = new AlleleInfo(teoAlleleInfoInFile, AGPv1VDfile, teoAlleleInfoOutFile, 11);  // the 64 base tag of the B73 allele is in column 11 (0 based)
    }

    public static void writeTeoPhyGenMappedReads() {
        String teoGenMappedReadsFile =      "/usr/local/maizediv/illumina/teosinte/phyGenMappedReadsTeoDiv5r10_20110423.bin";
        String teoGenMappedReadsTextFile =  "/usr/local/maizediv/illumina/teosinte/phyGenMappedReadsTeoDiv5r10_20110423.txt";
        String rbtFileTeoWithW22QTL10 =     "/usr/local/maizediv/illumina/teosinte/counts/teoRBT_withW22QTL10_Min3_20110610.bin";
        String tagsOfInterestFile =         "/usr/local/maizediv/illumina/teosinte/tags.chr10.LauraXO.txt";
        String rbtTextFileTeoWithW22QTL10 = "/usr/local/maizediv/illumina/teosinte/teoRBT_withW22QTL10_LauraXO_20110615.txt";

        ReadsWPhysicalMap rwpmTeoGenPhyMapped = new ReadsWPhysicalMap(teoGenMappedReadsFile,true);
        rwpmTeoGenPhyMapped.sortTable(false);  // sort by position
        rwpmTeoGenPhyMapped.writeCountFile(new File(teoGenMappedReadsTextFile), false);

        ReadsByTaxa teoRBT=new ReadsByTaxa(rbtFileTeoWithW22QTL10,true);
        teoRBT.filterForListOfReads(new File(tagsOfInterestFile), new File(rbtTextFileTeoWithW22QTL10), false);
    }

    public static void makeReadCountsSorghum() {
        String qseqDir =           "/cbsufsrv3/data2/maizediv/illumina/633Y5AAXX/qseq/sorghum";
        String qseqKey =           "/usr/local/maizediv/illumina/sorghum/RILs/Sorghum_RILs_key.txt";
        String parsedTaxaDir =     "/usr/local/maizediv/illumina/sorghum/RILs/taxa";
        String collapsedTaxaDir =  "/usr/local/maizediv/illumina/sorghum/RILs/taxacollapse";
        String combinedReadsFile = "/usr/local/maizediv/illumina/sorghum/RILs/counts/CombReadsSorghumMin10_20110506.bin";
        String rbtFile =           "/usr/local/maizediv/illumina/sorghum/RILs/counts/sorghumRBT_Min10_20110506.bin";

        // Read and parse qseq or fastq file using low stringency Q score (maximal # of reads)
        System.out.println();
        System.out.println("Parsing the qseq files for the Sorgum RILs");
        System.out.println("------------------------------------------");
        ParseBarcodeFiles pbf=new ParseBarcodeFiles(qseqDir, qseqKey, parsedTaxaDir, 0, false);
        System.gc();

        // condense the taxa ReadCounts files by sorting and collapsing identical reads
        System.out.println();
        System.out.println("Collapsing the ReadCounts files for each duplicate Sorghum taxon");
        System.out.println("----------------------------------------------------------------");
        ReadCounts rc = new ReadCounts(parsedTaxaDir, collapsedTaxaDir, 1, true, false, true);  // final true = combine identical taxa
        rc = null;
        System.gc();

        // Fuse multiple taxa files into one, minimum frequency is 10
        System.out.println();
        System.out.println("Combining the ReadCounts files across all Sorghum taxa");
        System.out.println("-------------------------------------------------------");
        CombineReadCounts crc = new CombineReadCounts(collapsedTaxaDir, combinedReadsFile, 10, true);
        crc = null;
        System.gc();

        CreateReadsByTaxa rbt = new CreateReadsByTaxa(combinedReadsFile, collapsedTaxaDir, rbtFile, 10, true);
    }

    public static void physicallyMapSorghumTags() {
        // input files
        String vdFile =            "/usr/local/maizediv/genome/sorghum/Sbicolor_79_RM_042611_10chrs_VD.rwpm.bin";
        String combinedReadsFile = "/usr/local/maizediv/illumina/sorghum/RILs/counts/CombReadsSorghumMin10_20110506.bin";

        // output files
        String phyMappedReadsFile1 = "/usr/local/maizediv/illumina/sorghum/RILs/counts/phyMappedReadsSorghumDiv1_20110509.bin"; // allow only one mismatch -- I want really solid alignments for an initial framework
        String phyMappedReadsFile5 = "/usr/local/maizediv/illumina/sorghum/RILs/counts/phyMappedReadsSorghumDiv5_20110516.bin"; // up to 5 mismatches

        int nMismatches = 5;

        ReadsWPhysicalMap rwpmTeoPhyMapped=MapHaplotypes.blastUsingPhysicalMap(new ReadsWPhysicalMap(vdFile,true),
                  new ReadsWPhysicalMap(new ReadCounts(combinedReadsFile, true)),nMismatches);
        rwpmTeoPhyMapped.sortTable(false); // sort by position
        rwpmTeoPhyMapped.writeCountFile(new File(phyMappedReadsFile5), Integer.MAX_VALUE, true, false, Float.NaN, true);
    }

    public static void makeGBSFrameworkSorghum() {
        // input files
        String rbtFile =           "/usr/local/maizediv/illumina/sorghum/RILs/counts/sorghumRBT_Min10_20110506.bin";
        String phyMappedReadsFile = "/usr/local/maizediv/illumina/sorghum/RILs/counts/phyMappedReadsSorghumDiv1_20110509.bin";

        // output files
        String GBSHapMapFileStem = "/usr/local/maizediv/illumina/sorghum/RILs/20110510/sorghum_GBS_biallelic_genos_20110503";

        ReadsByTaxa rbt = new ReadsByTaxa(rbtFile,true);
        ReadsWPhysicalMap rwpmPhyMapped = new ReadsWPhysicalMap(phyMappedReadsFile,true);

        for(int i=1; i<=10; i++) {
            System.out.println("WORKING ON CHROMOSOME " + i);
            Alignment a = new ConvertGBSToAlignment(rbt, rwpmPhyMapped, (byte)i);
            System.out.println("UnFiltered Alignment  Taxa:" + a.getSequenceCount() + " Sites:" + a.getSiteCount());
            int[] goodLowHetSites = AlignmentFilterByGBSUtils.getLowHetSNPs(a, 0.10, 100);
            a = FilterAlignment.getInstance(a, goodLowHetSites);
            MutableAlignmentForGBS mutAlign=new MutableAlignmentForGBS(a);
            a = null;
            System.gc();
            mutAlign.writeToHapmap(false, false, GBSHapMapFileStem + ".chr" + i, '\t', i);  
        }
    }

    public static void filterGBSFrameworkSorghum() {
        // input files
        String phasedGenoFile = "C:/Users/jcg233/Documents/Bioinformatics/NextGen/Sorghum/20110510/sorghum_GBS_biallelic_genos_20110503.chr9.filtered5.phased.hmp.txt";

        // output files
        String GBSHapMapFileStem = "C:/Users/jcg233/Documents/Bioinformatics/NextGen/Sorghum/20110525/sorghum_GBS_biallelic_genos_20110524";

        int i = 9;
        Alignment phasedGenos = ImportUtils.readFromHapmap(phasedGenoFile);
        MutableAlignmentForGBS  mutAlign = new MutableAlignmentForGBS(phasedGenos);
        phasedGenos = null;
        mutAlign.removeDCOs(4);
        System.out.println("nDCOs: "+ mutAlign.getnDCOs());
        mutAlign.callHetSegments(4);
        System.out.println("nDCOs: "+ mutAlign.getnDCOs());  // it's possible for callHetSegments to flag more DCOs
        mutAlign.writeToHapmap(false, true, GBSHapMapFileStem + ".chr" + i + ".hetSegs", '\t', i);
        phasedGenos = ImportUtils.readFromHapmap(GBSHapMapFileStem + ".chr" + i + ".hetSegs" + ".hmp.txt");
        mutAlign = new MutableAlignmentForGBS(phasedGenos);
        mutAlign.imputeMissingDataIncludingHets();
        mutAlign.writeToHapmap(false, false, GBSHapMapFileStem + ".chr" + i + ".imputed", '\t', i);
    }

    public static void mapSorghumVsImputedGBSFrame() {
        // input files
        String frameworkMapFile =   "/usr/local/maizediv/illumina/sorghum/RILs/20110525/sorghum_GBS_biallelic_genos_20110524.chr9.imputed.hmp.txt";
        String rbtFile =            "/usr/local/maizediv/illumina/sorghum/RILs/counts/sorghumRBT_Min10_20110506.bin";
        String phyMappedReadsFile = "/usr/local/maizediv/illumina/sorghum/RILs/counts/phyMappedReadsSorghumDiv5_20110516.bin"; // up to 5 mismatches

        // output files
        String GBSHapMapFileStem =     "/usr/local/maizediv/illumina/sorghum/RILs/20110525/sorghum_GBS_biallelic_genos_20110525";
        String genPhyMappedReadsFile = "/usr/local/maizediv/illumina/sorghum/RILs/20110525/genPhyMappedReadsSorgum_20110525.bin";

        ReadsByTaxa rbt=new ReadsByTaxa(rbtFile, true);
        ReadsWPhysicalMap rwpmPhyMapped = new ReadsWPhysicalMap(phyMappedReadsFile,true);
        rwpmPhyMapped=MapHaplotypes.mapGeneticallyVsFramework(frameworkMapFile, rwpmPhyMapped, rbt, 0.10, 20, false);  // multipleChrs = false
        rwpmPhyMapped.writeCountFile(new File(genPhyMappedReadsFile), Integer.MAX_VALUE, true, false, Float.NaN, true);
        rwpmPhyMapped = null;
        System.gc();
        ReadsWPhysicalMap rwpmGenPhyMapped = new ReadsWPhysicalMap(genPhyMappedReadsFile,true);

        int chr = 9;
        Alignment a = new ConvertGBSToAlignment(rwpmGenPhyMapped, rbt, 0.05, 0.05, (byte)chr);
        System.out.println("UnFiltered Alignment  Taxa:" + a.getSequenceCount() + " Sites:" + a.getSiteCount());
        int[] goodLowHetSites = AlignmentFilterByGBSUtils.getLowHetSNPs(a, 0.20, 10);
        a = FilterAlignment.getInstance(a, goodLowHetSites);
        MutableAlignmentForGBS  mutAlign = new MutableAlignmentForGBS(a);
        a = null;
        System.gc();
        mutAlign.removeDCOs(7);
        System.out.println("nDCOs on chr" + chr + " = " + mutAlign.getnDCOs());
        mutAlign.callHetSegments(7);
        System.out.println("nDCOs on chr" + chr + " = " + mutAlign.getnDCOs());  // it's possible for callHetSegments to flag more DCOs
        mutAlign.writeToHapmap(false, true, GBSHapMapFileStem + ".chr" + chr + ".hetSegs", '\t', chr);
        Alignment filteredGenos = ImportUtils.readFromHapmap(GBSHapMapFileStem + ".chr" + chr + ".hetSegs" + ".hmp.txt");
        mutAlign = new MutableAlignmentForGBS(filteredGenos);
        System.out.println("Imputed Genos: " + mutAlign.imputeMissingDataIncludingHets());
        System.out.println("nGenos written to imputed HapMap file: " + mutAlign.writeToHapmap(false, false, GBSHapMapFileStem + ".chr" + chr + ".imputed", '\t', chr));
    }

    public static void alignB73Htrhm() {
        String flowcell = "42A87AAXX";
        String DirStem = "Q:/SolexaAnal/Solexa_NGG";
        String qseqKeyFile = DirStem + "/" + flowcell + "/" + flowcell + "_Multiplex_key.txt";
        String qseqDir = DirStem + "/" + flowcell + "/qseq";
        String parsedTaxaDir = DirStem + "/" + flowcell + "/qseq/taxa";
        String collapsedTaxaDir = DirStem + "/" + flowcell + "/qseq/taxacollapse";
        String B73HtrhmReadsFile = DirStem + "/" + flowcell + "/qseq/taxacollapse/B73Htrhm_" + flowcell + "_6.cnt";
        String B73AGPv1VDfile = "Q:/SolexaAnal/IBM/virtual_digests/AGPv1/B73AGPv1VD.rwpm.bin";
        String B73HtrhmRwpmFile = DirStem + "/" + flowcell + "/qseq/counts/B73Htrhm.Q10.M5.rwpm.txt";

        // Read and parse qseq or fastq file using Q score cutoff of 10
        System.out.println("Parsing the qseq files for " + flowcell + " lane 6");
        System.out.println("-------------------------------------------");
        ParseBarcodeFiles pbf=new ParseBarcodeFiles(qseqDir, qseqKeyFile, parsedTaxaDir, 10, false);
        pbf = null;
        System.gc();

        // condense the taxa ReadCounts files by sorting and collapsing identical reads
        System.out.println();
        System.out.println("Collapsing the ReadCounts files for each taxon in " + flowcell + " lane 6");
        System.out.println("------------------------------------------------------------------");
        ReadCounts rc=new ReadCounts(parsedTaxaDir, collapsedTaxaDir, 1, true, false, false);
        rc = null;
        System.gc();

        // open the B73Htrhm reads and align them against B73 AGPv1 with no mismatch
        System.out.println();
        System.out.println("Aligning the ReadCounts from B73Htrhm vs AGPv1 (<6 mismatches)");
        System.out.println("--------------------------------------------------------------");
        ReadCounts B73HtrhmReadCounts = new ReadCounts(B73HtrhmReadsFile,true);
        ReadsWPhysicalMap B73HtrhmRwpm=MapHaplotypes.blastUsingPhysicalMap(new ReadsWPhysicalMap(B73AGPv1VDfile,true),
                  new ReadsWPhysicalMap(B73HtrhmReadCounts),5);
        B73HtrhmRwpm.writeCountFile(new File(B73HtrhmRwpmFile), false);
    }

    public static void analyzeIBMforGBSmethodsPaper() {
        String combinedReadsFileIBM1 = "/cbsufsrv3/data2/maizediv/illumina/434LFAAXX/counts/CombReads_IBM1_Min5_20100525.bin";  // min Q score of 0
        String AGPv1VDfile = "Q:/SolexaAnal/virtual_digests/AGPv1/B73AGPv1VD.rwpm.bin";
        String IBMrwpmM0File = "/cbsufsrv3/data2/maizediv/illumina/434LFAAXX/counts/IBM1perfectMatches.rwpm.txt";
        String IBMrwpmM5File = "/cbsufsrv3/data2/maizediv/illumina/434LFAAXX/counts/IBM1lt5misMatches.rwpm.txt";
//        String AGPv2VDfile = "Q:/SolexaAnal/virtual_digests/AGPv1/B73AGPv2VD.rwpm.bin";

        CompareReadDistribution crd1=new CompareReadDistribution(CompareReadDistribution.Analysis.EvalVirtualDigest, AGPv1VDfile, null);
        CompareReadDistribution crd2=new CompareReadDistribution(CompareReadDistribution.Analysis.EvalVirtualDigest, AGPv1VDfile, combinedReadsFileIBM1);

        ReadCounts IBMReadCounts = new ReadCounts(combinedReadsFileIBM1,true);
        ReadsWPhysicalMap IBMrwpm=MapHaplotypes.blastUsingPhysicalMap(new ReadsWPhysicalMap(AGPv1VDfile,true), new ReadsWPhysicalMap(IBMReadCounts),5);
        IBMrwpm.writeCountFile(new File(IBMrwpmM5File), false);
    }

    public static void testDuplicateTaxon() {
        String qseqDirIBM =           "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/434LFAAXX/qseq";
        String origKeyIBM =           "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/434LFAAXX/duplicate_taxon_test/434LFAAXX_original_key.txt";
        String duplKeyIBM =           "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/434LFAAXX/duplicate_taxon_test/434LFAAXX_duplicates_key.txt";
        String parsedTaxaDirOrig =    "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/434LFAAXX/duplicate_taxon_test/original";
        String collapsedTaxaDirOrig = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/434LFAAXX/duplicate_taxon_test/orig_collapsed";
        String parsedTaxaDirDupl =    "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/434LFAAXX/duplicate_taxon_test/with_duplicates";
        String collapsedTaxaDirDupl = "C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/434LFAAXX/duplicate_taxon_test/dupl_collapsed";

//        System.out.println();
//        System.out.println("Parsing the qseq files for IBM (434LFAAXX): with DUPLICATES");
//        System.out.println("-----------------------------------------------------------");
//        ParseBarcodeFiles pbfIBMdupl=new ParseBarcodeFiles(qseqDirIBM, duplKeyIBM, parsedTaxaDirDupl, 10, false);
//        pbfIBMdupl = null;
//        System.gc();

        ReadCounts rcIBMdupl=new ReadCounts(parsedTaxaDirDupl, collapsedTaxaDirDupl, 1, true, false, true);
    }
}

