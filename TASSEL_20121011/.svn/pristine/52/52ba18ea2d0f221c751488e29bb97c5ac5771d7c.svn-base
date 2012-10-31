package net.maizegenetics.gbs.pipeline;

import java.awt.Frame;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.Arrays;
import net.maizegenetics.util.MultiMemberGZIPInputStream;
import javax.swing.ImageIcon;
import net.maizegenetics.gbs.homology.ParseBarcodeRead;
import net.maizegenetics.gbs.homology.ReadBarcodeResult;
import net.maizegenetics.gbs.tagdist.TagCountMutable;
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.tagdist.UTagCountMutable;
import net.maizegenetics.gbs.util.ArgsEngine;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.util.DirectoryCrawler;
import org.apache.log4j.Logger;

/** 
 * Derives a tagCount list for each fastq file in the input directory.
 *
 * Keeps only good reads having a barcode and a cut site and no N's in the
 * useful part of the sequence.  Trims off the barcodes and truncates sequences
 * that (1) have a second cut site, or (2) read into the common adapter.
 * 
 */
public class UFastqToTagCountPlugin extends AbstractPlugin {  
    static long timePoint1;
    private ArgsEngine engine = null;
    private Logger logger = Logger.getLogger(UFastqToTagCountPlugin.class);
	private String parentDir = null;
    String directoryName=null;
    String keyfile=null;
    String enzyme = null;
    int maxGoodReads = 200000000;
    int minCount =1;
    String outputDir=null;

    public UFastqToTagCountPlugin() {
        super(null, false);
    }

    public UFastqToTagCountPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    private void printUsage(){
        logger.info(
             "\n\nUsage is as follows:\n"
			+ " -w  Working directory to contain subdirectories\n"
            + " -e  Enzyme used to create the GBS library\n"
            + " -c  Maximum number of good, barcoded reads per lane. Default: 200,000,000"
        );
    }

    public DataSet performFunction(DataSet input){
        File qseqDirectory = new File(directoryName);
        if (!qseqDirectory.isDirectory()) {
            printUsage();
            throw new IllegalStateException("The input name you supplied is not a directory.");
        }
        //countTags(keyfile, enzyme, directoryName, outputDir, maxGoodReads, minCount);
		this.countTagsOutput(keyfile, enzyme, parentDir, outputDir, maxGoodReads, minCount, 2);
        return null;
    }

    @Override
    public void setParameters(String[] args) {
        if(args.length==0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
//        try{
            if(engine == null){
                engine = new ArgsEngine();
                engine.add("-w", "--working-directory", true);
                engine.add("-e", "--enzyme", true);
                engine.add("-s", "--max-reads", true);
                engine.add("-c", "--min-count", true);
                engine.parse(args);
            }

            if (engine.getBoolean("-w")) {
				parentDir = engine.getString("-w");
				File pdf = new File (parentDir);
				directoryName = new File (pdf, UCreatWorkingDirPlugin.childDir[0]).toString();
				outputDir = new File (pdf, UCreatWorkingDirPlugin.childDir[2]).toString();
				File[] potenKey = new File (pdf, UCreatWorkingDirPlugin.childDir[1]).listFiles();
				for (int i = 0; i < potenKey.length; i++) {
					if (potenKey[i].getName().startsWith(".")) continue;
					keyfile = potenKey[i].getAbsolutePath();
					break;
				}
				//keyfile = new File (pdf, UCreatWorkingDirPlugin.childDir[1]).listFiles()[0].toString();
				System.out.println(keyfile);
			}
            else{ printUsage(); throw new IllegalArgumentException("Please specify the working directory."); }

            if(engine.getBoolean("-e")){ enzyme = engine.getString("-e"); }
            else{ printUsage(); throw new IllegalArgumentException("Please specify the enzyme used to create the GBS library.");}

            if(engine.getBoolean("-s")){ maxGoodReads = Integer.parseInt(engine.getString("-s")); }

            if (engine.getBoolean("-c")) { minCount = Integer.parseInt(engine.getString("-c")); }
//        }catch (Exception e){
//            System.out.println("Caught exception while setting parameters of "+this.getClass()+": "+e);
//        }
    }

    /**
     * Derives a tagCount list for each fastq file in the fastqDirectory.
     *
     * @param keyFileS        A key file (a sample key by barcode, with a plate map included).
     * @param enzyme          The enzyme used to create the library (currently ApeKI or PstI).
     * @param fastqDirectory  Directory containing the fastq files (will be recursively searched).
     * @param outputDir       Directory to which the tagCounts files (one per fastq file) will be written.
     * @param maxGoodReads    The maximum number of barcoded reads expected in a fastq file
     * @param minCount        The minimum number of occurrences of a tag in a fastq file for it to be included in the output tagCounts file
     */
    public static void countTags(String keyFileS, String enzyme, String fastqDirectory, String outputDir, int maxGoodReads, int minCount) {
        BufferedReader br;
        String[] countFileNames = null;

        File inputDirectory = new File(fastqDirectory);
        File[] fastqFiles = DirectoryCrawler.listFiles("(?i).*\\.fq$|.*\\.fq\\.gz$|.*\\.fastq$|.*_fastq\\.txt$|.*_fastq\\.gz$|.*_fastq\\.txt\\.gz$|.*_sequence\\.txt$|.*_sequence\\.txt\\.gz$", inputDirectory.getAbsolutePath());
//                                                      (?i) denotes case insensitive;                 \\. denotes escape . so it doesn't mean 'any char' & escape the backslash
        if(fastqFiles.length==0 || fastqFiles == null){
            System.out.println("Couldn't find any files that end with \".fq\", \".fq.gz\", \".fastq\", \"_fastq.txt\", \"_fastq.gz\", \"_fastq.txt.gz\", \"_sequence.txt\", or \"_sequence.txt.gz\" in the supplied directory.");
        } else {
            System.out.println("Using the following FASTQ files:");
            countFileNames = new String[fastqFiles.length];
            for (int i=0; i<fastqFiles.length; i++) {
                countFileNames[i] = fastqFiles[i].getName().replaceAll
                    ("(?i)\\.fq$|\\.fq\\.gz$|\\.fastq$|_fastq\\.txt$|_fastq\\.gz$|_fastq\\.txt\\.gz$|_sequence\\.txt$|_sequence\\.txt\\.gz$", ".cnt");
//                        \\. escape . so it doesn't mean 'any char' & escape the backslash                
                System.out.println(fastqFiles[i].getAbsolutePath());
            }
        }
        int allReads=0, goodBarcodedReads=0;
        for(int laneNum=0; laneNum<fastqFiles.length; laneNum++) {
            System.out.println("Reading FASTQ file: "+fastqFiles[laneNum]);
            String[] filenameField=fastqFiles[laneNum].getName().split("_");
            ParseBarcodeRead thePBR;  // this reads the key file and store the expected barcodes for this lane
            if(filenameField.length==3) {thePBR=new ParseBarcodeRead(keyFileS, enzyme, filenameField[0], filenameField[1]);}
            else if(filenameField.length==4) {thePBR=new ParseBarcodeRead(keyFileS, enzyme, filenameField[0], filenameField[2]);}  // B08AAABXX_s_1_sequence.txt.gz
            else if(filenameField.length==5) {thePBR=new ParseBarcodeRead(keyFileS, enzyme, filenameField[1], filenameField[3]);}
            else {
                System.out.println("Error in parsing file name:");
                System.out.println("   The filename does not contain either 3, 4, or 5 underscore-delimited values.");
                System.out.println("   Expect: flowcell_lane_fastq.txt.gz OR flowcell_s_lane_fastq.txt.gz OR code_flowcell_s_lane_fastq.txt.gz");
                System.out.println("   Filename: "+fastqFiles[laneNum]);
                continue;
            }
            System.out.println("Total barcodes found in lane:"+thePBR.getBarCodeCount());
            if(thePBR.getBarCodeCount() == 0){
                System.out.println("No barcodes found.  Skipping this flowcell lane."); continue;
            }
            String[] taxaNames=new String[thePBR.getBarCodeCount()];
            for (int i = 0; i < taxaNames.length; i++) {
                taxaNames[i]=thePBR.getTheBarcodes(i).getTaxaName();
            }
			Arrays.sort(taxaNames);
			UTagCountMutable[] theTCs = new UTagCountMutable[taxaNames.length];
            try{
                //Read in qseq file as a gzipped text stream if its name ends in ".gz", otherwise read as text
                if(fastqFiles[laneNum].getName().endsWith(".gz")){
                    br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(fastqFiles[laneNum]))));
                }else{
                    br=new BufferedReader(new FileReader(fastqFiles[laneNum]),65536);
                }
                String sequence="", qualityScore="";
                String temp;

                try{
                    for (int i = 0; i < theTCs.length; i++) {
						theTCs[i] = new UTagCountMutable(2, 0);
					}
                }catch(OutOfMemoryError e){
                    System.out.println(
                        "Your system doesn't have enough memory to store the number of sequences"+
                        "you specified.  Try using a smaller value for the minimum number of reads."
                    );
                }
                int currLine=0;
                allReads = 0;
                goodBarcodedReads = 0;
                while (((temp = br.readLine()) != null) && goodBarcodedReads < maxGoodReads) {
                    currLine++;
                    try{
                        //The quality score is every 4th line; the sequence is every 4th line starting from the 2nd.
                        if((currLine+2)%4==0){
                            sequence = temp;
                        }else if(currLine%4==0){
                            qualityScore = temp;
                            allReads++;
                            //After quality score is read, decode barcode using the current sequence & quality  score
                            ReadBarcodeResult rr = thePBR.parseReadIntoTagAndTaxa(sequence, qualityScore, true, 0);
                            if (rr != null){
                                goodBarcodedReads++;
                                int hit = Arrays.binarySearch(taxaNames, rr.getTaxonName());
								theTCs[hit].addReadCount(rr.getRead(), rr.getLength(), 1);
                            }
                            if (allReads % 1000000 == 0) {
                                System.out.println("Total Reads:" + allReads + " Reads with barcode and cut site overhang:" + goodBarcodedReads);
                            }
                        }
                    }catch(NullPointerException e){
                        System.out.println("Unable to correctly parse the sequence and "
                        + "quality score from fastq file.  Your fastq file may have been corrupted.");
                        System.exit(0);
                    }
                }
                System.out.println("Total number of reads in lane=" + allReads);
                System.out.println("Total number of good barcoded reads=" + goodBarcodedReads);
                System.out.println("Timing process (sorting, collapsing, and writing TagCount to file).");
                timePoint1 = System.currentTimeMillis();
               for (int i = 0; i < taxaNames.length; i++) {
					theTCs[i].toArray();
					theTCs[i].collapseCounts();
					taxaNames[i] = taxaNames[i].replaceAll(":", "_");
					String[] tem = taxaNames[i].split(":");
					String outfileS = new File(outputDir, taxaNames[i]).getAbsolutePath() + ".cnt";
					theTCs[i].writeTagCountFile(outfileS, FilePacking.Bit, minCount);
				}
                System.out.println("Process took " + (System.currentTimeMillis() - timePoint1) + " milliseconds.");
                br.close();

            } catch(Exception e) {
                System.out.println("Catch testBasicPipeline c="+goodBarcodedReads+" e="+e);
                e.printStackTrace();
            }
            System.out.println("Finished reading "+(laneNum+1)+" of "+fastqFiles.length+" sequence files.");
        }
    }

	 private void countTagsOutput(String keyFileS, String enzyme, String fastqDirectory, String outputDir, int maxGoodReads, int minCount, int tagLengthInLong) {
		File pd = new File (parentDir);
		String tagParseDirS = new File (pd, "parse").getAbsolutePath();

        BufferedReader br;
        String[] countFileNames = null;
        File inputDirectory = new File(fastqDirectory);
        File[] fastqFiles = DirectoryCrawler.listFiles("(?i).*\\.fq$|.*\\.fq\\.gz$|.*\\.fastq$|.*_fastq\\.txt$|.*_fastq\\.gz$|.*_fastq\\.txt\\.gz$|.*_sequence\\.txt$|.*_sequence\\.txt\\.gz$", inputDirectory.getAbsolutePath());
//                                                      (?i) denotes case insensitive;                 \\. denotes escape . so it doesn't mean 'any char' & escape the backslash
        if(fastqFiles.length==0 || fastqFiles == null){
            System.out.println("Couldn't find any files that end with \".fq\", \".fq.gz\", \".fastq\", \"_fastq.txt\", \"_fastq.gz\", \"_fastq.txt.gz\", \"_sequence.txt\", or \"_sequence.txt.gz\" in the supplied directory.");
        } else {
            System.out.println("Using the following FASTQ files:");
            countFileNames = new String[fastqFiles.length];
            for (int i=0; i<fastqFiles.length; i++) {
                countFileNames[i] = fastqFiles[i].getName().replaceAll
                    ("(?i)\\.fq$|\\.fq\\.gz$|\\.fastq$|_fastq\\.txt$|_fastq\\.gz$|_fastq\\.txt\\.gz$|_sequence\\.txt$|_sequence\\.txt\\.gz$", ".cnt");
//                        \\. escape . so it doesn't mean 'any char' & escape the backslash
                System.out.println(fastqFiles[i].getAbsolutePath());
            }
        }
        int allReads=0, goodBarcodedReads=0;
        for(int laneNum=0; laneNum<fastqFiles.length; laneNum++) {
            System.out.println("Reading FASTQ file: "+fastqFiles[laneNum]);
            String[] filenameField=fastqFiles[laneNum].getName().split("_");
            ParseBarcodeRead thePBR;  // this reads the key file and store the expected barcodes for this lane
            if(filenameField.length==3) {thePBR=new ParseBarcodeRead(keyFileS, enzyme, filenameField[0], filenameField[1]);}
            else if(filenameField.length==4) {thePBR=new ParseBarcodeRead(keyFileS, enzyme, filenameField[0], filenameField[2]);}  // B08AAABXX_s_1_sequence.txt.gz
            else if(filenameField.length==5) {thePBR=new ParseBarcodeRead(keyFileS, enzyme, filenameField[1], filenameField[3]);}
            else {
                System.out.println("Error in parsing file name:");
                System.out.println("   The filename does not contain either 3, 4, or 5 underscore-delimited values.");
                System.out.println("   Expect: flowcell_lane_fastq.txt.gz OR flowcell_s_lane_fastq.txt.gz OR code_flowcell_s_lane_fastq.txt.gz");
                System.out.println("   Filename: "+fastqFiles[laneNum]);
                continue;
            }
            System.out.println("Total barcodes found in lane:"+thePBR.getBarCodeCount());
            if(thePBR.getBarCodeCount() == 0){
                System.out.println("No barcodes found.  Skipping this flowcell lane."); continue;
            }

            String[] taxaNames=new String[thePBR.getBarCodeCount()];
            for (int i = 0; i < taxaNames.length; i++) {
                taxaNames[i]=thePBR.getTheBarcodes(i).getTaxaName();
            }
			Arrays.sort(taxaNames);
			File tagParseDir = new File(tagParseDirS);
			tagParseDir.mkdir();
			ParseOut[] temOut = new ParseOut[taxaNames.length];
			for (int i = 0; i < taxaNames.length; i++) {
				String tempOutfileS = new File(tagParseDirS, taxaNames[i].replaceAll(":", "_")).getAbsolutePath() + ".par";
				temOut[i] = new ParseOut(tempOutfileS);
			}


            try{
                //Read in qseq file as a gzipped text stream if its name ends in ".gz", otherwise read as text
                if(fastqFiles[laneNum].getName().endsWith(".gz")){
                    br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(fastqFiles[laneNum]))));
                }else{
                    br=new BufferedReader(new FileReader(fastqFiles[laneNum]),65536);
                }
                String sequence="", qualityScore="";
                String temp;
                int currLine=0;
                allReads = 0;
                goodBarcodedReads = 0;
                while (((temp = br.readLine()) != null) && goodBarcodedReads < maxGoodReads) {
                    currLine++;
                    try{
                        //The quality score is every 4th line; the sequence is every 4th line starting from the 2nd.
                        if((currLine+2)%4==0){
                            sequence = temp;
                        }else if(currLine%4==0){
                            qualityScore = temp;
                            allReads++;
                            //After quality score is read, decode barcode using the current sequence & quality  score
                            ReadBarcodeResult rr = thePBR.parseReadIntoTagAndTaxa(sequence, qualityScore, true, 0);
                            if (rr != null){
                                goodBarcodedReads++;
                                int hit = Arrays.binarySearch(taxaNames, rr.getTaxonName());
								temOut[hit].writeTempFile(rr.getRead(), rr.getLength());
                            }
                            if (allReads % 1000000 == 0) {
                                System.out.println("Total Reads:" + allReads + " Reads with barcode and cut site overhang:" + goodBarcodedReads);
                            }
                        }
                    }catch(NullPointerException e){
                        System.out.println("Unable to correctly parse the sequence and "
                        + "quality score from fastq file.  Your fastq file may have been corrupted.");
                        System.exit(0);
                    }
                }
                System.out.println("Total number of reads in lane=" + allReads);
                System.out.println("Total number of good barcoded reads=" + goodBarcodedReads);
                System.out.println("Timing process (sorting, collapsing, and writing TagCount to file).");
                timePoint1 = System.currentTimeMillis();

				for (int i = 0; i < temOut.length; i++) {
					temOut[i].closeOut();
				}
				File[] parsedFiles = tagParseDir.listFiles();
				for (int i = 0; i < parsedFiles.length; i++) {
					ParseIn pi = new ParseIn(parsedFiles[i], tagLengthInLong);
					String outfileS = outputDir + "/" + parsedFiles[i].getName().replace("par", "cnt");
					TagCountMutable tcm = pi.getTagCountMutable();
					tcm.writeTagCountFile(outfileS, FilePacking.Bit, minCount);
				}
				for (int i = 0; i < parsedFiles.length; i++) {
					parsedFiles[i].delete();
				}
				tagParseDir.delete();
                System.out.println("Process took " + (System.currentTimeMillis() - timePoint1) + " milliseconds.");
                br.close();

            } catch(Exception e) {
                System.out.println("Catch testBasicPipeline c="+goodBarcodedReads+" e="+e);
                e.printStackTrace();
            }
            System.out.println("Finished reading "+(laneNum+1)+" of "+fastqFiles.length+" sequence files.");
        }
    }

	private class ParseIn {
		DataInputStream dis;
		int tagNum;
		int tagLengthInLong;
		ParseIn (File infile, int tagLengthInLong) {
			try {
				this.dis = new DataInputStream (new BufferedInputStream (new FileInputStream(infile), 65536));
			}
			catch (Exception e) {
				System.out.println(e.toString());
			}
			tagNum = (int)infile.length()/(8*tagLengthInLong+1);
			this.tagLengthInLong = tagLengthInLong;
		}

		TagCountMutable getTagCountMutable () {
			TagCountMutable tc = new TagCountMutable (tagLengthInLong, tagNum);
			try {
				long[] tag = new long[tagLengthInLong];
				for (int i = 0; i < tagNum; i++) {
					for (int j = 0; j < tagLengthInLong; j++) {
						tag[j] = dis.readLong();
					}
					byte tagLength = dis.readByte();
					tc.addReadCount(tag, tagLength, 1);
				}
				dis.close();
			}
			catch (Exception e) {
				System.out.println(e.toString());
			}
			tc.collapseCounts();
			return tc;
		}
	}

	private class ParseOut {
		DataOutputStream dos;

		ParseOut (String outfileS) {
			try {
				this.dos = new DataOutputStream (new BufferedOutputStream (new FileOutputStream(outfileS), 65536));
			}
			catch (Exception e) {
				System.out.println(e.toString());
			}
		}

		void writeTempFile (long[] seq, byte tagLength) {
			try {
				for (int i = 0; i < seq.length; i++) {
                    dos.writeLong(seq[i]);
                }
				dos.write(tagLength);
			}
			catch (Exception e) {
				System.out.println(e.toString());
			}
		}

		void closeOut () {
			try {
				dos.flush();
				dos.close();
			}
			catch (Exception e) {
				System.out.println(e.toString());
			}
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
