/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS;

/**
 *
 * @author JeffGlaubitz
 */
public class GBSPipelinesJesse {
    static String baseDir="/usr/local/maizediv/";
    static String qseqDir=baseDir+"illumina/433V4AAXX/Jesse/qseq";
    static String qseqKey=baseDir+"illumina/433V4AAXX/Jesse/433V4AAXX_key.txt";
    static String parsedTaxaDir=baseDir+"illumina/433V4AAXX/Jesse/taxa";
    static String collapsedTaxaDir=baseDir+"illumina/433V4AAXX/Jesse/taxacollapse";
    static String combinedFlowcellsReads=baseDir+"illumina/433V4AAXX/Jesse/barleyReadCounts_Q0_Min5_20100519.bin";
    static String myReadsByTaxa=baseDir+"illumina/433V4AAXX/Jesse/barleyReadsByTaxa_Q0_Min5_20100519.txt";

    public GBSPipelinesJesse() {
    }

    public static void main(String[] args) {
        JessePipelineStrict();
//        JessePipelineLax();
    }

    public static void JessePipelineStrict() {
        // Read and parse qseq or fastq file using high stringency Q score (to get a good set of reference reads)
        ParseBarcodeFiles pbfStingent=new ParseBarcodeFiles(qseqDir, qseqKey, parsedTaxaDir, 10, false);
        pbfStingent = null;
        System.gc();

        // condense the taxa ReadCounts files by sorting and collapsing identical reads: binary=true, simpleFilter=false, combineIdenticalTaxa=true
        ReadCounts rcStringent=new ReadCounts(parsedTaxaDir, collapsedTaxaDir, 1, true, false, true);
        rcStringent = null;
        System.gc();

        // Fuse multiple taxa files into one, minimum frequency is 1, binary=true (low min freq b/c may only have 1 high quality read for a given allele)
        CombineReadCounts myCRCs=new CombineReadCounts(collapsedTaxaDir, combinedFlowcellsReads, 1, true);
        myCRCs = null;
        System.gc();

        // Read and parse qseq or fastq file once again, this time with a min Qscore of 0 (low stringency for maximal data completeness)
        ParseBarcodeFiles pbfLax=new ParseBarcodeFiles(qseqDir, qseqKey, parsedTaxaDir, 0, false);

        // once again, collapse the folder by sorting and collapse identical reads: binary=true, simpleFilter=false, combineIdenticalTaxa=true
        ReadCounts rcLax=new ReadCounts(parsedTaxaDir, collapsedTaxaDir, 1, true, false, true);
        System.gc();

        CreateReadsByTaxa myCRBT=new CreateReadsByTaxa(combinedFlowcellsReads, collapsedTaxaDir, myReadsByTaxa, false);
    }

    public static void JessePipelineLax() {
        // Read and parse qseq or fastq file using low stringency Q score (maximal # of reads)
        ParseBarcodeFiles pbfLax=new ParseBarcodeFiles(qseqDir, qseqKey, parsedTaxaDir, 0, false);
        System.gc();

        // condense the taxa ReadCounts files by sorting and collapsing identical reads: binary=true, simpleFilter=false, combineIdenticalTaxa=true
        ReadCounts rcLax=new ReadCounts(parsedTaxaDir, collapsedTaxaDir, 1, true, false, true);
        System.gc();

        // Fuse multiple taxa files into one, minimum frequency is 5, binary=true (reads that are seen repeatedly are hopefully real)
        CombineReadCounts myCRCs=new CombineReadCounts(collapsedTaxaDir, combinedFlowcellsReads, 5, true);
        System.gc();

        // Read and parse qseq or fastq file once again, this time with a min Qscore of 0 (low stringency for maximal data completeness)
        pbfLax=new ParseBarcodeFiles(qseqDir, qseqKey, parsedTaxaDir, 0, false);
        System.gc();

        // once again, collapse the folder by sorting and collapse identical reads: binary=true, simpleFilter=false, combineIdenticalTaxa=true
        rcLax=new ReadCounts(parsedTaxaDir, collapsedTaxaDir, 1, true, false, true);
        System.gc();

        CreateReadsByTaxa myCRBT=new CreateReadsByTaxa(combinedFlowcellsReads, collapsedTaxaDir, myReadsByTaxa, false);
    }
}
