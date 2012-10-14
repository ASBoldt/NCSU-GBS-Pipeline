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
import java.util.ArrayList;
import java.util.Arrays;
//import org.apache.tools.bzip2.CBZip2InputStream;  // See http://www.kohsuke.org/bzip2//
import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.gbs.tagdist.Tags;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.tagdist.TagsByTaxaBit;
import net.maizegenetics.gbs.util.ArgsEngine;
import net.maizegenetics.genome.BaseEncoder;
import net.maizegenetics.util.DirectoryCrawler;

/**
 * Reads sequences from HapMap2 project fastq files (random sheared paired end) 
 * and performs a virtual digest on them to extract 64 base GBS tags and fill in 
 * a TBT for the provided list of master tags.  Only one taxon in each TBT file.
 * 
 * @author jcg233
 */
public class RandomPairedEndToTBT {
    static long timePoint1;
    private RandomPairedEndToTBT() {
    }

    public static void main(String[] args) {
        String[] fastqFileNames = null;
        String enzyme = null;
        String outputDir = null;
        Tags masterTags = null;
        if (args.length < 4 || args == null) {
            System.out.println(
                "\nUsage is as follows:\n"
                    + "-i  Input directory containing HapMap fastq files\n"
                    + "-e  Enzyme used to create the GBS library\n"
                    + "-o  Output directory\n"
                    + "One of either:\n"
                    + "    -m  Physical map file containing alignments, OR\n"
                    + "    -t  Tag count file\n"
            );
            System.exit(0);
        } else {
            ArgsEngine engine = new ArgsEngine();
            engine.add("-i", "--input-directory", true);
            engine.add("-e", "--enzyme", true);
            engine.add("-o", "--output-directory", true);
            engine.add("-m", "--physical-map", true);
            engine.add("-t", "--tag-count", true);
            engine.parse(args);

            if (engine.getBoolean("-i")) {
                File qseqDirectory = new File(engine.getString("-i"));
                if (!qseqDirectory.isDirectory()) {
                    System.out.println("The input name you supplied (option -i) is not a directory.");
                    System.exit(0);
                }
                fastqFileNames = DirectoryCrawler.listFileNames(".*.txt.bz2$|.*.txt$", qseqDirectory.getAbsolutePath());
                if (fastqFileNames.length == 0 || fastqFileNames == null) {
                    System.out.println("Couldn't find any files that end with \"txt.bz2\" in the supplied directory.");
                } else {
                    Arrays.sort(fastqFileNames);
                    System.out.println("Using the following compressed fastq files:");
                    for (String filename : fastqFileNames) {
                        System.out.println(filename);
                    }
                }
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

            //Create TOPM object from TOPM file with option -m, or from tag count file with option -t
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
        parseTagsFromRandomShearedReads(fastqFileNames, masterTags, enzyme, outputDir);
    }


    /**
     * Uses an existing readCounts (or TagsOnPhysicalMap) object to create one TagsByTaxa file 
     * for each pair of fastq files in the input directory (only 1 taxon per TBT file).
     *
     * Assumes the fastq files are from random sheared libraries with paired end sequencing.  
     * Objective is to populate a TBT for an existing set of GBS tags with high coverage data
     * from random sheared libraries.  These high coverage samples will be useful for imputation.
     * Output TBT files written to the outputDir, using fastq file names minus the _1 or _2 and 
     * with extension changed to .tbt.bin (or .tbt.txt).
     *
     * @param fastQFileNames  Array of fastq file names (from paired end sequencing of random sheared DNA)
     * @param theMasterTags  A Tags object containing the GBS tags of interest
     * @param enzyme  The enzyme used to make the GBS libraries from which theMasterTags were derived (currently ApeKI only)
     * @param outputDir  String containing the path of the output directory
     */
    public static void parseTagsFromRandomShearedReads(String[] fastQFileNames, Tags theMasterTags, String enzyme, String outputDir) {
        TagsByTaxa theTBT = null;
        String prevSampFCLn = "NONE";
        int readsWithTag=0, allReads=0, tagsFound=0, tagsMatching=0;
        for(int FileIndex=0; FileIndex<fastQFileNames.length; FileIndex++) {
            boolean tooShort = false;
            System.out.println("\nWorking on fastq file: " + fastQFileNames[FileIndex]);
            File fastqFile=new File(fastQFileNames[FileIndex]);
            String[] np=fastqFile.getName().split("_");
            String sampFCLn;
            if(np.length==4) {sampFCLn = np[0] + ":" + np[1] + ":" + np[2];}  // sample:flowcell:lane
            else {
                System.out.println("Error in parsing file name:");
                System.out.println("   The filename does not contain 4 underscore-delimited values.");
                System.out.println("   Expect: sample_flowcell_lane_end.txt.bz2 where \'end\' is 1 or 2 (i.e., which paired end)");
                System.out.println("   Filename: "+fastQFileNames[FileIndex]);
                return;
            }
            if (!sampFCLn.equals(prevSampFCLn)) {
                String[] taxonName=new String[1];
                taxonName[0] = sampFCLn;
                theTBT=null;
                System.gc();
                theTBT=new TagsByTaxaBit(taxonName, theMasterTags);
                readsWithTag=0; allReads=0; tagsFound=0; tagsMatching=0;
            }

            String temp="";
            try{
                BufferedReader br;
//                if(fastQFileNames[FileIndex].endsWith(".bz2")){
//                    FileInputStream fis = new FileInputStream(fastQFileNames[FileIndex]);
//                    fis.read(); fis.read();  // extract the header "BZ"  // See http://www.kohsuke.org/bzip2//
//                    br = new BufferedReader(new InputStreamReader(new CBZip2InputStream(fis)));
//                }else{
                    br=new BufferedReader(new FileReader(fastQFileNames[FileIndex]),65536);
//                }
                String sl="", qualS="";
                while ((temp = br.readLine()) != null) {
                    if(!temp.startsWith("@")) continue;    //fastq starts with @ rather than >
                    sl=br.readLine();
                    if (allReads == 0 && sl.length() < 70) {
                        System.out.println("Reads are too short ("+sl.length()+" bases) in this lane...skipping");
                        tooShort = true;
                        break;  // no need to read the whole file
                    }
                    if(br.readLine().startsWith("+")) {
                        qualS = br.readLine();
                    }
                    allReads++;
                    TagRecord[] tagsInRead = findTags(sl, enzyme);
                    if (tagsInRead != null && tagsInRead.length > 0) {
                        readsWithTag++;
                        for (int i=0; i<tagsInRead.length; ++i) {
                            tagsFound++;
                            int h=theTBT.getTagIndex(tagsInRead[i].getTag());
                            int t=theTBT.getIndexOfTaxaName(sampFCLn);
                            if(h>-1) {
                                theTBT.addReadsToTagTaxon(h, t, 1);
                                tagsMatching++;
                            }
                        }
                    }
                    if(allReads%1000000==0) System.out.println("Total Reads:"+allReads+"   Reads containing a tag:"+readsWithTag
                            +"   Tags found:"+tagsFound+"   Tags matching:"+tagsMatching);
                }
                br.close();
            }
            catch(Exception e) {
                System.out.println("Catch testBasicPipeline c="+readsWithTag+" e="+e);
                System.out.println(temp);
                e.printStackTrace();
            }
            int filesDone = FileIndex + 1;
            System.out.println("Finished reading "+filesDone+" of "+fastQFileNames.length+" sequence files: "+fastQFileNames[FileIndex]);
            if (!tooShort && (sampFCLn.equals(prevSampFCLn) || sampFCLn.equals("P39:42UWMAAXX:8"))) {  // the 2nd read  failed for P39:42UWMAAXX:8 (only one file)
                System.out.println("Finished reading both ends of sample "+sampFCLn);
                System.out.println("Timing process (writing TagsByTaxa file)..."); timePoint1 = System.currentTimeMillis();
                File outfile;
                FilePacking outFormat = FilePacking.Bit;
                String outFileS = outputDir+fastQFileNames[FileIndex].substring(fastQFileNames[FileIndex].lastIndexOf(File.separator));
                if (outFormat == FilePacking.Text) {
                    outfile = new File(outFileS.replaceAll("_1.txt$|_1.txt.bz2$|_2.txt$|_2.txt.bz2$", ".tbt.txt"));
                } else {
                    outfile = new File(outFileS.replaceAll("_1.txt$|_1.txt.bz2$|_2.txt$|_2.txt.bz2$", ".tbt.bin"));
                }
                theTBT.writeDistFile(outfile,outFormat, 0);  // minCount of 0
                System.out.println("...process took "+(System.currentTimeMillis()-timePoint1)+" milliseconds.");
                System.out.println("Total number of reads in lane=" + allReads);
                System.out.println("Number of reads with >=1 ApeKI GBS tag="+readsWithTag);
                System.out.println("Number of tags found in the reads="+tagsFound);
                System.out.println("Number of tags matching one of the master tags="+tagsMatching);
            }
            prevSampFCLn = sampFCLn;
        }
    }
    
    
    /**
     * Returns any full length (64 base) GBS tags found in a read.  If two or more
     * cut sites are found, the corresponding GBS tags are returned.
     * 
     * @param read  A read from random-sheared, paired end Illumina sequencing
     * @param enzyme  The enzyme used to make the GBS master tags of interest
     * @return a long[][] containing the GBS tags found in the read (assumes 64 base)
     */
    public static TagRecord[] findTags(String read, String enzyme) {
        String [] recogSeqs;
        String tagAsString;
        ArrayList<TagRecord> tags = new ArrayList<TagRecord>();
        String polyA = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        
        if(enzyme.matches("(?i)apek[i1]")) {
            recogSeqs = new String[]{"GCAGC","GCTGC"};
        } else {
            System.out.println("Unrecognized cut site.  Only ApeKI is currently recognized.");
            return null;
        }
        ArrayList<Integer> cutStartsAL = new ArrayList<Integer>();
        int cutStart, lookFrom;
        for (int cs=0; cs<recogSeqs.length; cs++) {
            lookFrom = 0;
            do {
                cutStart = read.indexOf(recogSeqs[cs], lookFrom);
                if (cutStart != -1) {
                    cutStartsAL.add(cutStart);
                    lookFrom = cutStart + 1;
                }
            } while (cutStart != -1);
        }
        if (cutStartsAL.size()>0) {
            Integer[] cutStarts = cutStartsAL.toArray(new Integer[0]);
            Arrays.sort(cutStarts);
            if (cutStarts[0] >= 60) {  // leftmost tag
                tagAsString = read.substring(cutStarts[0]-60, cutStarts[0]+4);
                tags.add(new TagRecord(BaseEncoder.getLongArrayFromSeq(BaseEncoder.getReverseComplement(tagAsString))));
            }
            if (cutStartsAL.size() > 1) { // middle tags
                for (int cs=1; cs<cutStarts.length; cs++) {
                    if (cutStarts[cs]-cutStarts[cs-1]>3) {  // avoid overlapping cut sites (GCWGCWGC)
                        tagAsString = read.substring(cutStarts[cs-1]+1,cutStarts[cs]+4);
                        String pad = "";
                        if (tagAsString.length()<64) {
                            pad = polyA.substring(0, 64-tagAsString.length());
                            tags.add(new TagRecord(BaseEncoder.getLongArrayFromSeq(tagAsString+pad))); // top strand
                            tags.add(new TagRecord(BaseEncoder.getLongArrayFromSeq(BaseEncoder.getReverseComplement(tagAsString)+pad)));  // bottom strand
                        } else if (tagAsString.length()>64) {
                            String topTag = tagAsString.substring(0, 64);
                            tags.add(new TagRecord(BaseEncoder.getLongArrayFromSeq(topTag+pad))); // top strand
                            tagAsString = tagAsString.substring(tagAsString.length()-64,tagAsString.length());
                            tags.add(new TagRecord(BaseEncoder.getLongArrayFromSeq(BaseEncoder.getReverseComplement(tagAsString)+pad)));  // bottom strand
                        }
                    }
                }
            }
            if (read.length()-cutStarts[cutStarts.length-1] >= 65) { // rightmost tag
                tagAsString = read.substring(cutStarts[cutStarts.length-1]+1, cutStarts[cutStarts.length-1]+65);
                tags.add(new TagRecord(BaseEncoder.getLongArrayFromSeq(tagAsString)));
            }
            TagRecord[] tagsInRead = tags.toArray(new TagRecord[0]);
            return tagsInRead;
        }
        return null;
    }

}

    
class TagRecord {
    long[] tag = new long[2];
    public TagRecord(long[] tag) {
        this.tag = tag;
    }

    public long[] getTag() {
        return tag;
    }
}
