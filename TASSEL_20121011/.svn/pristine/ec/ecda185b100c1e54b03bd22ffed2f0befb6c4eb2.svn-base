/*
 * QseqToTBTPlugin
 */
package net.maizegenetics.gbs.pipeline;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import net.maizegenetics.util.MultiMemberGZIPInputStream;
import net.maizegenetics.gbs.homology.ParseBarcodeRead;
import net.maizegenetics.gbs.homology.ReadBarcodeResult;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.gbs.tagdist.Tags;
import net.maizegenetics.gbs.util.ArgsEngine;
import net.maizegenetics.util.DirectoryCrawler;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import java.awt.Frame;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeMap;
import javax.swing.ImageIcon;
import net.maizegenetics.gbs.homology.TagMatchFinder;
import net.maizegenetics.gbs.util.MutableSimpleAlignment;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.Locus;
import net.maizegenetics.pal.datatype.DataType;
import org.apache.log4j.ConsoleAppender;
import org.apache.log4j.Logger;
import org.apache.log4j.SimpleLayout;


/**
 * This pipeline converts a series of qseq files to TagsByTaxa files (one per qseq file).
 * It requires a list of existing tags (Tags object), which may come from a TagCounts file or TOPM file.
 *
 * @author james
 */
public class QseqToHapMapPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(QseqToHapMapPlugin.class);
    private ArgsEngine myArgsEngine = null;

    private String[] myQseqFileS = null;
    private String   myKeyFile = null;
    private String   myEnzyme = null;
    private String   myOutputDir = null;
    private int      myMinCount = 1;
    private static TagsOnPhysicalMap  topm = null;
    private static TagMatchFinder tmf = null;
    private static int maxDivergence=0;
    public static int[][] TagVarToMSASite;
    public QseqToHapMapPlugin() {
        super(null, false);
    }

    public QseqToHapMapPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    @Override
    public DataSet performFunction(DataSet input) {

        if ((myEnzyme == null) || (myEnzyme.length() == 0)) {
            printUsage();
            throw new IllegalStateException("performFunction: enzyme must be set.");
        }
        // TODO - More checks to validate parameters...

        matchTagsToTaxa(myQseqFileS, myKeyFile, myEnzyme, topm, myOutputDir, myMinCount);
        return null;
    }

    private void printUsage() {
        myLogger.info(
            "\nUsage is as follows:\n"
                + "-i  Input directory containing .qseq files\n"
                + "-k  Barcode key file\n"
                + "-e  Enzyme used to create the GBS library\n"
                + "-o  Output directory\n"
                + "-c  Minimum taxa count within a qseq file for a tag to be output (default 1)\n"  // Nb: using TagsByTaxaBit, so max count PER TAXON = 1
                + "-m  Physical map file containing alignments\n"
                + "-d  Maximum divergence (edit distance) between new read and previously mapped read\n");
    }

    @Override
    public void setParameters(String[] args) {
        if(args.length==0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        
        if (myArgsEngine == null) {
            myArgsEngine = new ArgsEngine();
            myArgsEngine.add("-i", "--input-directory", true);
            myArgsEngine.add("-k", "--key-file", true);
            myArgsEngine.add("-e", "--enzyme", true);
            myArgsEngine.add("-o", "--output-directory", true);
            myArgsEngine.add("-c", "--min-count", true);
            myArgsEngine.add("-t", "--tag-count", true);
            myArgsEngine.add("-m", "--physical-map", true);
            myArgsEngine.add("-d", "--divergence", true);
        }
        
        myArgsEngine.parse(args);
        myLogger.addAppender(new ConsoleAppender(new SimpleLayout()));
        
        String tempDirectory = myArgsEngine.getString("-i");
        if (tempDirectory != null) {
            File qseqDirectory = new File(tempDirectory);
            if (!qseqDirectory.isDirectory()) {
                printUsage();
                throw new IllegalArgumentException("setParameters: The input name you supplied is not a directory: " + tempDirectory);
            }
            myQseqFileS = DirectoryCrawler.listFileNames(".*qseq\\.txt$|.*qseq\\.txt\\.gz$", qseqDirectory.getAbsolutePath());
            if (myQseqFileS.length == 0 || myQseqFileS == null) {
                printUsage();
                throw new IllegalArgumentException("Couldn't find any files that end with \"qseq.txt\" or \"qseq.txt.gz\" in the supplied directory: " + tempDirectory);
            } else {
                myLogger.info("QseqToTBTPlugin: setParameters: Using the following .qseq files:");
                for (String filename : myQseqFileS) {
                    myLogger.info(filename);
                }
            }
        }
        
        if (myArgsEngine.getBoolean("-k")) {
            myKeyFile = myArgsEngine.getString("-k");
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify a key file (option -k).");
        }
        if (myArgsEngine.getBoolean("-e")) {
            myEnzyme = myArgsEngine.getString("-e");
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify the enzyme used to create the GBS library.");
        }
        if (myArgsEngine.getBoolean("-o")) {
            myOutputDir = myArgsEngine.getString("-o");
            File outDirectory = new File(myOutputDir);
            if (!outDirectory.isDirectory()) {
                printUsage();
                throw new IllegalArgumentException("The output name you supplied (option -o) is not a directory: " + myOutputDir);
            }
            outDirectory = null;
        }else{
            printUsage();
            throw new IllegalArgumentException("Please specify an output directory (option -o).");
        }
        if (myArgsEngine.getBoolean("-c")) {
            myMinCount = Integer.parseInt(myArgsEngine.getString("-c"));
        } else {
            myMinCount = 1;
        }

        if (myArgsEngine.getBoolean("-d")) {
            maxDivergence = Integer.parseInt(myArgsEngine.getString("-d"));
        }
        
        // Create Tags object from tag count file with option -t, or from TOPM file with option -m
        if (myArgsEngine.getBoolean("-m")){
            topm = new TagsOnPhysicalMap(myArgsEngine.getString("-m"), true);
//            topm.writeTextFile(new File("/media/Data/Paola/test/output.topm.txt"));
            tmf = new TagMatchFinder(topm);
        }else{
            printUsage();
            throw new IllegalArgumentException("Please specify a TagsOnPhysicalMap file (-m)");
        }
    }


    /**
     * Uses an existing Tags object to create one TagsByTaxa file for each qseq file in the input directory.
     *
     * Output TBT files written to the outputDir, using qseq file names with extension changed to .tbt.bin (or .tbt.txt)
     *
     * @param qseqFileS      Array of qseq file names (Illumina-created files with raw read sequence, quality score, machine name, etc.)
     * @param keyFileS       A key file (list of taxa by barcode, lane & flow cell, including plate maps)
     * @param enzyme         The enzyme used to make the library (currently ApeKI or PstI)
     * @param theMasterTags  A Tags object: list of tags to be included in the final TBT
     * @param outputDir      String contaJining the path of the output directory to contain tags-by-taxa files
     * @param minCount       The minimum number of times a tag must show up in a qseq file before it is included in the corresponding TBT file
     */
    public static void matchTagsToTaxa(String[] qseqFileS, String keyFileS, String enzyme, Tags theMasterTags, String outputDir, int minCount) {

        myLogger.info("Counting sites in TOPM file.");
        ArrayList<int[]> uniquePositions = new ArrayList<int[]>();
        int[] chromosomes = topm.getChromosomes();
        Locus[] loci = topm.getLoci();
        for (int i = 0; i < chromosomes.length; i++) {
            uniquePositions.add(topm.uniquePositions(chromosomes[i]));
        }

         for(int laneNum=0; laneNum<qseqFileS.length; laneNum++) {
            System.out.println("\nWorking on qseq file: " + qseqFileS[laneNum]);
            System.gc();
            int goodMatchedReads=0, allReads=0, goodMatched=0;
            File qseqFile=new File(qseqFileS[laneNum]);
            String[] np=qseqFile.getName().split("_");

            //Create a new object to hold barcoded tags.  The constructor can optionally process a group of fastq
            //files.  A minimum quality score for inclusion of a read can also be provided.
            ParseBarcodeRead thePBR;
            if(np.length==3) {thePBR=new ParseBarcodeRead(keyFileS, enzyme, np[0], np[1]);}
            else if(np.length==4) {thePBR=new ParseBarcodeRead(keyFileS, enzyme, np[0], np[2]);}
            else if(np.length==5) {thePBR=new ParseBarcodeRead(keyFileS, enzyme, np[1], np[3]);}
            else {
                System.out.println("Error in parsing file name:");
                System.out.println("   The filename does not contain either 3, 4, or 5 underscore-delimited values.");
                System.out.println("   Expect: flowcell_lane_qseq.txt OR flowcell_s_lane_qseq.txt OR code_flowcell_s_lane_qseq.txt");
                System.out.println("   Filename: "+qseqFileS[laneNum]);
                return;
            }
            System.out.println("Total barcodes found in lane:"+thePBR.getBarCodeCount());
            if(thePBR.getBarCodeCount() == 0){
                System.out.println("No barcodes found.  Skipping this flowcell lane."); continue;
            }

            myLogger.info("Creating alignments.");
            String[] taxaNames=thePBR.getTaxaNames();
            MutableSimpleAlignment[] outMSA =new MutableSimpleAlignment[chromosomes.length];
             for (int i = 0; i < outMSA.length; i++) {
                 outMSA[i]=new MutableSimpleAlignment(taxaNames, uniquePositions.get(i).length, new Locus[] {loci[i]} );
             }
            HashMap<String,Integer> taxaNameIndices = outMSA[0].taxonMap(); //Find the indices of taxa names within the MSA for quick lookup

            //Transfer locus, strand, & position from TOPM to appropriate MSA
            myLogger.info("Adding sites from TOPM file to alignments.");
             for (int i = 0; i < outMSA.length; i++) {
            int currSite=0;
                 for (int j = 0; j < uniquePositions.get(i).length; j++) {
                     int position=uniquePositions.get(i)[j];
                     outMSA[i].setLocusOfSite(currSite, Integer.toString(chromosomes[i]));
                     outMSA[i].setStrandOfSite(currSite, (byte)'+');
                     outMSA[i].setPositionOfSite(currSite, position);
                     currSite++;
                 }
                outMSA[i].sortSiteByPhysicalPosition();
             }
            
            // Read the qseq file and assign set base for a given site
            myLogger.info("Looking for known SNPs in sequence reads...");
            String temp="";
            goodMatchedReads=0; allReads=0; goodMatched=0; int perfectMatches=0; int imperfectMatches=0; int singleImperfectMatches=0;
            try{
                BufferedReader br;
                //Read in qseq file as a gzipped text stream if its name ends in ".gz", otherwise read as text
                if(qseqFileS[laneNum].endsWith(".gz")){
                    br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(qseqFileS[laneNum]))));
                }else{
                    br=new BufferedReader(new FileReader(qseqFileS[laneNum]),65536);
                }
                String sl, qualS="";
                while (((temp = br.readLine()) != null)) {
                    String[] jj=temp.split("\\s");
                    allReads++;
                    if(allReads%1000000==0){
                        System.out.println(
                                "Total Reads:"+allReads+
                                " goodReads:"+goodMatchedReads+
                                " goodMatched: "+goodMatched+
                                " perfect matches: "+perfectMatches+
                                " near matches: "+imperfectMatches+
                                " single imperfect matches: "+singleImperfectMatches
                        );
                    }
                    sl=jj[8];
                    qualS=jj[9];
                    ReadBarcodeResult rr=thePBR.parseReadIntoTagAndTaxa(sl, qualS, true, 0);
                    if (rr != null) {
                        goodMatchedReads++;
                        int tagIndex = topm.getTagIndex(rr.getRead());
                        if(tagIndex>=0){
                            perfectMatches++;
//                            TreeMap<Integer, Integer> bestHitsAndDiv = tmf.findMatchesWithIntLengthWords(rr.getRead(), 3, false);
//                            System.out.println(bestHitsAndDiv.toString()+" "+tagIndex);
                        }

                        //If binary search for exact match fails, try searching for near matches.  Continue if this also fails.
                        if (tagIndex < 0) {
                            TreeMap<Integer, Integer> bestHitsAndDiv = tmf.findMatchesWithIntLengthWords(rr.getRead(), 1, false);
                            if (bestHitsAndDiv.size() > 0) {
                                imperfectMatches++;
                                if(bestHitsAndDiv.size()==1)singleImperfectMatches++;
                                tagIndex = bestHitsAndDiv.firstKey();
                            }else{
                                continue;
                            }
                        }
                        
                        int taxon = taxaNameIndices.get(rr.getTaxonName());
                        
                            goodMatched++;
                            int chromosome = topm.getChromosome(tagIndex);
                            if(chromosome== Integer.MIN_VALUE) continue;
                            int startPos = topm.getStartPosition(tagIndex);
                            Locus locus = topm.getLocus(tagIndex);
                            int chrIndex=0;
                            for (int i = 0; i < chromosomes.length; i++) {
                                if(chromosomes[i]==chromosome){chrIndex=i; break;}
                            }
                            for (int variant = 0; variant < topm.maxVariants; variant++){
                                byte currBase = topm.getVariantDef(tagIndex, variant);
                                if( (currBase==topm.byteMissing) || (currBase==DataType.UNKNOWN_BYTE) ) continue;
                                int offset= topm.getVariantPosOff(tagIndex, variant);
                                int pos=startPos+offset;
                                int currSite=outMSA[chrIndex].getSiteOfPhysicalPosition(pos, locus);
                                if(currSite<0) continue;
                                byte prevBase = outMSA[chrIndex].getBase(taxon, currSite);
//
                                if(prevBase == DataType.UNKNOWN_BYTE) {
                                    outMSA[chrIndex].setBase(taxon, currSite, currBase);
                                } else if (currBase != prevBase) {
                                    outMSA[chrIndex].setBase(taxon,currSite,TagsToSNPByAlignmentMTPlugin.resolveSNPByteFromCallPair(prevBase, currBase));
                                }
                            }
                        }                    
                }
                br.close();
            } catch(Exception e) {
                System.out.println("Catch testBasicPipeline c="+goodMatchedReads+" e="+e);
                System.out.println(temp);
                e.printStackTrace();
            }
            

            for (int i = 0; i < outMSA.length; i++) {

                outMSA[i].sortSiteByPhysicalPosition();
                AlignmentFilterByGBSUtils.getCoverage_MAF_F_Dist(outMSA[i], false);
//                 for (int j = 0; j < outMSA[i].getSiteCount(); j++) {
//                     int pos = outMSA[i].getPositionInLocus(j);
//                     if(JVHPipeline.searchArray(JVHPipeline.sitesUniqueInDiscovery, pos)){
//                         for(Integer taxon: outMSA[i].taxonMap().values()){
//                             System.out.print((char)outMSA[i].getBase(taxon.intValue(), j));
//                         }
//                         System.out.println();
//                    }
//                }

                 
                 String outFileS = outputDir + qseqFileS[laneNum].substring(qseqFileS[laneNum].lastIndexOf(File.separator));
                 outFileS = outFileS.replaceAll("_qseq\\.txt$|_qseq\\.txt\\.gz$", ".c"+chromosomes[i]);  // ".hmp.txt" gets added by ExportUtils.writeToHapmap
                 ExportUtils.writeToHapmap(outMSA[i], false, outFileS, '\t');
             }
            
            System.out.println("Total number of reads in lane=" + allReads);
            System.out.println("Total number of good, barcoded reads="+goodMatchedReads);
            int filesDone = laneNum + 1;
            System.out.println("Finished reading "+filesDone+" of "+qseqFileS.length+" sequence files: "+qseqFileS[laneNum]+"\n");
        }
    }

    @Override
    public ImageIcon getIcon() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getButtonName() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getToolTipText() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
