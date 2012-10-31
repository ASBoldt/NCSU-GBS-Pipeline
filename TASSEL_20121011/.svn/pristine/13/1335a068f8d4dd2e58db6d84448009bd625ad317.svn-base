/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.gbs.pipeline;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.zip.GZIPInputStream;
import net.maizegenetics.gbs.homology.ParseBarcodeRead;
import net.maizegenetics.gbs.homology.ReadBarcodeResult;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.gbs.tagdist.Tags;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.tagdist.TagsByTaxaBit;
import net.maizegenetics.gbs.util.ArgsEngine;
import net.maizegenetics.util.DirectoryCrawler;

/**
 * This pipeline converts a series of qseq files to TagsByTaxa files (one per qseq file).
 * It requires a list of existing tags, which may come from a TOPM file or TagCounts file.
 *
 * @author james
 */
public class QseqToTBT {
    static long timePoint1;
    private QseqToTBT() {
    }

    public static void main(String[] args) {
        String[] qseqFileS = null;
        String keyFile = null;
        String enzyme = null;
        String outputDir = null;
        int minCount = 1;
        Tags masterTags = null;
        if (args.length < 5 || args == null) {
            System.out.println(
                "\nUsage is as follows:\n"
                    + "-i  Input directory containing .qseq files\n"
                    + "-k  Barcode key file\n"
                    + "-e  Enzyme used to create the GBS library\n"
                    + "-o  Output directory\n"
                    + "-c  Minimum taxa count within a qseq file for a tag to be output (default 1)\n"  // Nb: using TagsByTaxaBit, so max count PER TAXON = 1
                    + "One of either:\n"
                    + "-m  Physical map file containing alignments, OR\n"
                    + "-t  Tag count file\n"
            );
            System.exit(0);

        } else {
            ArgsEngine engine = new ArgsEngine();
            engine.add("-i", "--input-directory", true);
            engine.add("-k", "--key-file", true);
            engine.add("-e", "--enzyme", true);
            engine.add("-o", "--output-directory", true);
            engine.add("-c", "--min-count", true);
            engine.add("-m", "--physical-map", true);
            engine.add("-t", "--tag-count", true);
            engine.parse(args);

            if (engine.getBoolean("-i")) {
                File qseqDirectory = new File(engine.getString("-i"));
                if (!qseqDirectory.isDirectory()) {
                    System.out.println("The input name you supplied (option -i) is not a directory.");
                    System.exit(0);
                }
                qseqFileS = DirectoryCrawler.listFileNames(".*qseq.txt$|.*qseq.txt.gz$", qseqDirectory.getAbsolutePath());
                if (qseqFileS.length == 0 || qseqFileS == null) {
                    System.out.println("Couldn't find any files that end with \"qseq.txt\" or \"qseq.txt.gz\" in the supplied directory.");
                } else {
                    System.out.println("Using the following .qseq files:");
                    for (String filename : qseqFileS) {
                        System.out.println(filename);
                    }
                }
            }
            if (engine.getBoolean("-k")) {
                keyFile = engine.getString("-k");
            } else {
                System.out.println("Please specify a key file (option -k).");
                System.exit(0);
            }
            if(engine.getBoolean("-e")){
                enzyme = engine.getString("-e");
            } else {
                System.out.println("Please specify the enzyme used to create the GBS library (option -e).");
                System.exit(0);
            }
            if (engine.getBoolean("-o")) {
                outputDir = engine.getString("-o");
                File outDirectory = new File(outputDir);
                if (!outDirectory.isDirectory()) {
                    System.out.println("The output name you supplied (option -o) is not a directory.");
                    System.exit(0);
                }
                outDirectory = null;
            } else {
                System.out.println("Please specify an output directory (option -o).");
                System.exit(0);
            }
            if(engine.getBoolean("-c")) {
                minCount = Integer.parseInt(engine.getString("-c"));
            } else {
                minCount =1;
            }

            //Create TOPM object from TOPM file with option -m, or from read count file with option -t
            if (engine.getBoolean("-m")) {
                if (engine.getBoolean("-t")) {
                    System.out.println("Options -m and -t are mutually exclusive.");
                    System.exit(0);
                }
                masterTags = new TagsOnPhysicalMap(engine.getString("-m"), true);
            } else if (engine.getBoolean("-t")) {
                if (engine.getBoolean("-m")) {
                    System.out.println("Options -m and -t are mutually exclusive.");
                    System.exit(0);
                }
                masterTags = new TagCounts(engine.getString("-t"), FilePacking.Bit);
            } else {
                System.out.println("Please specify a TagsOnPhysicalMap file (-m) *OR* a readCounts file (-t)");
                System.exit(0);
            }
        }
        matchTagsToTaxa(qseqFileS, keyFile, enzyme, masterTags, outputDir, minCount); //Feed all data into function
    }


    /**
     * Uses an existing TagsOnPhysicalMap (or readCounts) object to create one TagsByTaxa file for each qseq file in the input directory.
     *
     * Output TBT files written to the outputDir, using qseq file names with extension changed to .tbt.bin (or .tbt.txt)
     *
     * @param keyFileS  A key file (list of taxa by barcode, lane & flow cell, including plate maps).
     * @param enzyme  The enzyme used to make the library (currently ApeKI or PstI)
     * @param qseqFileS  Array of qseq file names (Illumina-created files with raw read sequence, quality score, machine name, etc.)
     * @param theTOPM  A TagsOnPhysicalMap object: tags (barcodes trimmed off) with (-m) or without (-t) their genome coordinates.
     * @param outputDir  String containing the path of the output directory
     * @param minCount  A tag-by-taxa file (read sequence and presence/absence in one or more taxa).
     */
    public static void matchTagsToTaxa(String[] qseqFileS, String keyFileS, String enzyme, Tags theMasterTags, String outputDir, int minCount) {
        for(int laneNum=0; laneNum<qseqFileS.length; laneNum++) {
            System.out.println("\nWorking on qseq file: " + qseqFileS[laneNum]);
            TagsByTaxa theTBT=null;
            System.gc();
            int goodBarcodedReads=0, allReads=0, goodMatched=0;
            File qseqFile=new File(qseqFileS[laneNum]);
            String[] np=qseqFile.getName().split("_");

            //Create a new object to hold barcoded tags.  The constructor can optionally process a group of fastq
            //files, given a minimum quality score.
            ParseBarcodeRead thePBR;
            if(np.length==3) {thePBR=new ParseBarcodeRead(keyFileS, enzyme, np[0], np[1]);}
            else if(np.length==5) {thePBR=new ParseBarcodeRead(keyFileS, enzyme, np[1], np[3]);}
            else {
                System.out.println("Error in parsing file name:");
                System.out.println("   The filename does not contain either 3 or 5 underscore-delimited values.");
                System.out.println("   Expect: flowcell_lane_qseq.txt OR code_flowcell_s_lane_qseq.txt");
                System.out.println("   Filename: "+qseqFileS[laneNum]);
                return;
            }
            System.out.println("Total barcodes found in lane:"+thePBR.getBarCodeCount());
            if(thePBR.getBarCodeCount() == 0){
                System.out.println("No barcodes found.  Skipping this flowcell."); continue;
            }

            //Fill an array with taxon names.
            String[] taxaNames=new String[thePBR.getBarCodeCount()];
            for (int i = 0; i < taxaNames.length; i++) {
                taxaNames[i]=thePBR.getTheBarcodes(i).getTaxaName();//+":"+np[1]+":"+np[3];
            }

            theTBT=new TagsByTaxaBit(taxaNames,theMasterTags);

            //Find the union of .qseq entries and barcoded reads.
            int max=200000000;  // maximum number of good barcoded reads expected in a qseq file
            String temp="";
            goodBarcodedReads=0; allReads=0; goodMatched=0;
            try{
                //Read in qseq file as a gzipped text stream if its name ends in ".gz", otherwise read as text
                BufferedReader br;
                if(qseqFileS[laneNum].endsWith(".gz")){
                    br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(qseqFileS[laneNum]))));
                }else{
                    br=new BufferedReader(new FileReader(qseqFileS[laneNum]),65536);
                }
                String sl, qualS="";
                while (((temp = br.readLine()) != null)&&(goodBarcodedReads<max)) {
                    String[] jj=temp.split("\\s");
                    allReads++;
                    if(allReads%1000000==0) System.out.println("Total Reads:"+allReads+" goodReads:"+goodBarcodedReads+" goodMatched:"+goodMatched);
                    sl=jj[8];
                    qualS=jj[9];
                    ReadBarcodeResult rr=thePBR.parseReadIntoTagAndTaxa(sl, qualS, false, 0);
                    if(rr!=null) {
                        goodBarcodedReads++;
                        int t=theTBT.getIndexOfTaxaName(rr.getTaxonName());
                        int h=theTBT.getTagIndex(rr.getRead());
                        if(h>-1) {
                            theTBT.addReadsToTagTaxon(h, t, 1);
                            goodMatched++;
                        }
                    }
                }
                br.close();
            }
            catch(Exception e) {
                System.out.println("Catch testBasicPipeline c="+goodBarcodedReads+" e="+e);
                System.out.println(temp);
                e.printStackTrace();
            }
            System.out.println("Timing process (writing TagsByTaxa file)..."); timePoint1 = System.currentTimeMillis();
            File outfile;
            FilePacking outFormat = FilePacking.Bit;
            String outFileS = outputDir+qseqFileS[laneNum].substring(qseqFileS[laneNum].lastIndexOf(File.separator));
            if (outFormat == FilePacking.Text) {
                outfile = new File(outFileS.replaceAll(".txt$|.txt.gz$", ".tbt.txt"));
            } else {
                outfile = new File(outFileS.replaceAll(".txt$|.txt.gz$", ".tbt.bin"));
            }
            theTBT.writeDistFile(outfile,outFormat, minCount);
            System.out.println("...process took "+(System.currentTimeMillis()-timePoint1)+" milliseconds.");
            System.out.println("Total number of reads in lane=" + allReads);
            System.out.println("Total number of good, barcoded reads="+goodBarcodedReads);
            int filesDone = laneNum + 1;
            System.out.println("Finished reading "+filesDone+" of "+qseqFileS.length+" sequence files: "+qseqFileS[laneNum]+"\n");
        }
    }
}