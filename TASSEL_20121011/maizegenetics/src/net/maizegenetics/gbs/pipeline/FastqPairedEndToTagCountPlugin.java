package net.maizegenetics.gbs.pipeline;

import java.awt.Frame;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.text.DecimalFormat;
import java.io.PrintWriter;
import java.io.*;
import java.util.*;

import net.maizegenetics.util.MultiMemberGZIPInputStream;
import javax.swing.ImageIcon;
import net.maizegenetics.gbs.homology.ParseBarcodeRead;
import net.maizegenetics.gbs.homology.ReadBarcodeResult;
import net.maizegenetics.gbs.tagdist.TagCountMutable;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
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
public class FastqPairedEndToTagCountPlugin extends AbstractPlugin {  
    static long timePoint1;
    static long now;
    private ArgsEngine engine = null;
    private Logger logger = Logger.getLogger(FastqPairedEndToTagCountPlugin.class);
    String directoryName=null;
    String keyfile=null;
    String enzyme = null;
    //int maxGoodReads = 200000000;  //standrad default, can't run with 16gb ram
    //int maxGoodReads = 2000000;	// low tester
    int maxGoodReads = 150000000;	// high tester
    int minCount =1;
    String outputDir=null;
    protected static TreeMap <String, Integer> barcodePairs = new TreeMap<String, Integer>();
    protected static ArrayList<String> summaryOutputs = new ArrayList<String>();

    public FastqPairedEndToTagCountPlugin() {
        super(null, false);
    }

    public FastqPairedEndToTagCountPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    private void printUsage(){
        logger.info(
             "\n\nUsage is as follows:\n"
            + " -i  Input directory containing FASTQ files in text or gzipped text.\n"
            + "     NOTE: Directory will be searched recursively and should\n"
            + "     be written WITHOUT a slash after its name.\n\n"
            + " -k  Key file listing barcodes distinguishing the samples\n"
            + " -e  Enzyme used to create the GBS library, if it differs from the one listed in the key file.\n"
            + " -s  Max good reads per lane. (Optional. Default is 200,000,000).\n"
            + " -c  Minimum tag count (default is 1).\n"
            + " -o  Output directory to contain .cnt files (one per FASTQ file, defaults to input directory).\n\n"
        );
    }

    public DataSet performFunction(DataSet input){
        File qseqDirectory = new File(directoryName);
        if (!qseqDirectory.isDirectory()) {
            printUsage();
            throw new IllegalStateException("The input name you supplied is not a directory.");
        }
        countTags(keyfile, enzyme, directoryName, outputDir, maxGoodReads, minCount);  
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
                engine.add("-i", "--input-directory", true);
                engine.add("-k", "--key-file", true);
                engine.add("-e", "--enzyme", true);
                engine.add("-s", "--max-reads", true);
                engine.add("-c", "--min-count", true);
                engine.add("-o", "--output-file", true);
                engine.parse(args);
            }

            if (engine.getBoolean("-i")) { directoryName = engine.getString("-i");}
            else{ printUsage(); throw new IllegalArgumentException("Please specify the location of your FASTQ files."); }

            if(engine.getBoolean("-k")){ keyfile = engine.getString("-k");}
            else{ printUsage(); throw new IllegalArgumentException("Please specify a barcode key file.");}

            if(engine.getBoolean("-e")){ enzyme = engine.getString("-e"); }
            else{ 
                System.out.println("No enzyme specified.  Using enzyme listed in key file.");
            }
        
            if(engine.getBoolean("-s")){ maxGoodReads = Integer.parseInt(engine.getString("-s"));}

            if (engine.getBoolean("-c")) { minCount = Integer.parseInt(engine.getString("-c"));}

            if(engine.getBoolean("-o")){ outputDir = engine.getString("-o");}
            else{outputDir = directoryName;}
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
        
    	BufferedReader br1;
        BufferedReader br2;;

        String[] countFileNames = null;  // counter variable
        TreeMap <String, Integer> pairCount = new TreeMap<String, Integer>(); // stores paired sequences to write to file
        ArrayList <String> hashFileNames = new ArrayList<String>(); // stores names of files resulting from HashMap output
        ArrayList<String> badBarcodeRead2=new ArrayList<String>();
        summaryOutputs.add("Total Reads \t Forward Only \t Reverse Only \t Both");

        
        /* Grab ':' delimited key files */
        String[] tempFileList = keyFileS.split(":");
        String[] keyFileList = new String[2];
        
        if (tempFileList.length == 0){
        	System.out.println("No key file given");
        	keyFileList[0] = "GBS.key"; // = NULL ?
        	keyFileList[1] = "GBS.key"; // = NULL ?
        	
        } else if (tempFileList.length == 1){
        	System.out.println("Only one key file given");
        	keyFileList[0] = tempFileList[0];
        	keyFileList[1] = tempFileList[0];
        	
        } else {
        	keyFileList[0] = tempFileList[0];
        	keyFileList[1] = tempFileList[1];
        }
 
        File inputDirectory = new File(fastqDirectory);
        File[] fastqFiles = DirectoryCrawler.listFiles("(?i).*\\.fq$|.*\\.fq\\.gz$|.*\\.fastq$|.*_fastq\\.txt$|.*_fastq\\.gz$|.*_fastq\\.txt\\.gz$|.*_sequence\\.txt$|.*_sequence\\.txt\\.gz$", inputDirectory.getAbsolutePath());
      
        if(fastqFiles.length !=0 ){
        	Arrays.sort(fastqFiles);
            System.out.println("Using the following FASTQ files:");
            
            
            //COUNTS HOW MANY FILES THERE ARE IN THE INPUT            
            countFileNames = new String[fastqFiles.length];
            for (int i=0; i<fastqFiles.length; i++) {
                countFileNames[i] = fastqFiles[i].getName().replaceAll
                    ("(?i)\\.fq$|\\.fq\\.gz$|\\.fastq$|_fastq\\.txt$|_fastq\\.gz$|_fastq\\.txt\\.gz$|_sequence\\.txt$|_sequence\\.txt\\.gz$", ".cnt");
//                        \\. escape . so it doesn't mean 'any char' & escape the backslash    
                System.out.println(fastqFiles[i].getAbsolutePath());
            }
        }
        
        int allReads=0, goodBarcodedReads=0, goodBarcodedForwardReads=0, goodBarcodedReverseReads=0;
        int numFastqFiles = fastqFiles.length;  //number of files   
        int indexStartOfRead2 = numFastqFiles/2;  // index of where paired file should start, also number of pairs
        
        // check for even number of files
        checkForPairs(indexStartOfRead2);
                         
		/* sets mutildimensional array for 
		 * [x][][] 2 bays for forward and reverse reads
		 * [][x][] number of files expected for each directional read
		 * [][][y] will hold the parsed elements of the file name
		 */
		String [][][] fileReadInfo=new String[2][indexStartOfRead2][5]; 
		
		String[][][] filenameField= new String[2][5][];
		int fileNum=0;
		/* parses file name into array elements that correspond to reads and expected pairing
		 * where the first file would be paired with the file that is half of the total number 
		 * of files.
		 * 
		 * The outter loop controller is set to 2 because there should not be more than two 
		 * directional reads, forward and reverse.
		 */
		for(int read=0; read<2; read++) {
			int loopController, setStart; //control loops and arrays
			int fileController=0; // resets to 0 so files are copied in correct array
			
			//set conditions for the loop
			if(read==0){
				setStart=0;
				loopController=indexStartOfRead2;
			}
			else{
				setStart=indexStartOfRead2;
				loopController=numFastqFiles;
			}
				
				for(fileNum=setStart; fileNum<loopController; fileNum++) {			
					//following block could be set as separate private method
					File outputFile = new File(outputDir+File.separator+countFileNames[fileNum]);
						if(outputFile.isFile()){
				            System.out.println(
				                    "An output file "+countFileNames[fileNum]+"\n"+ 
				                    " already exists in the output directory for file "+fastqFiles[fileNum]+".  Skipping.");
							fileController++;
				            continue;
				        }
					System.out.println("Reading FASTQ file: "+fastqFiles[fileNum]);//print
					filenameField[read][fileNum]=fastqFiles[fileNum].getName().split("_");
					System.arraycopy(filenameField[read][fileNum], 0, fileReadInfo[read][fileController], 0, filenameField[read][fileNum].length);
					fileController++;
					}
		}
		
fileNum=0;


//handle keyfiles and enzymes
//2 arrays for manually inputing multiple enzymes and keys for testing
System.out.println("OLD Key file is:"+ keyFileS);            
System.out.println("OLD enzyme is:"+ enzyme);
String[] hcEnzyme={"PstI-MspI","MspI-PstI"};
String[] hcKeyFiles={"GBS.key","GBS2.key"};
		
			hcKeyFiles[1]=modifyKey2File(hcKeyFiles[1]);


 			/* 
 			 * Reads the key file and store the expected barcodes for a lane.
 			 * Set to a length of 2 to hold up to two key files' worth of information.
 			 * The convention will be that the forward read is [0] and the reverse
 			 * read is[1]
 			 */
            ParseBarcodeRead [] thePBR = new ParseBarcodeRead [2];  
          //  String[][] taxaNames=new String[2][];
            
            for(int b=0;b<indexStartOfRead2;b++){
            	
				if(fileReadInfo[0][b][0]!=null && fileReadInfo[0][b].length==5) {
					thePBR[0]=new ParseBarcodeRead(
							hcKeyFiles[0], hcEnzyme[0], fileReadInfo[0][b][1], fileReadInfo[0][b][3]);
				}
				else {
				 printParsingError();
				continue;
				}
				
				if(fileReadInfo[0][b][0]!=null && fileReadInfo[1][b].length==5) {
					thePBR[1]=new ParseBarcodeRead(
							hcKeyFiles[1], hcEnzyme[1], fileReadInfo[1][b][1], fileReadInfo[1][b][3]);	
				}
				else {
					printParsingError();
				continue;
				}
	
				System.out.println("\nTotal barcodes found in lane:"+thePBR[0].getBarCodeCount());
				if(thePBR[0].getBarCodeCount() == 0){
	                System.out.println("No barcodes found.  Skipping this flowcell lane."); continue;
	            }
				
				System.out.println("\nTotal barcodes found in lane:"+thePBR[1].getBarCodeCount());
				if(thePBR[1].getBarCodeCount() == 0){
	                System.out.println("No barcodes found.  Skipping this flowcell lane."); continue;
	            }
				// as far as I can tell, this bit of code is not used anywhere downstream.
				String[] taxaNamesF=new String[thePBR[0].getBarCodeCount()];
	            for (int i = 0; i < taxaNamesF.length; i++) {
	                taxaNamesF[i]=thePBR[0].getTheBarcodes(i).getTaxaName();
	            }
	            String[] taxaNamesR=new String[thePBR[1].getBarCodeCount()];
	            for (int i = 0; i < taxaNamesR.length; i++) {
	                taxaNamesR[i]=thePBR[1].getTheBarcodes(i).getTaxaName();
	            }
				
				try{
	                //Read in qseq file as a gzipped text stream if its name ends in ".gz", otherwise read as text
	                if(fastqFiles[b].getName().endsWith(".gz")){
	                    br1 = new BufferedReader(new InputStreamReader(
	                    						new MultiMemberGZIPInputStream(
	                    						new FileInputStream(fastqFiles[b]))));
	                    br2 = new BufferedReader(new InputStreamReader(
        										new MultiMemberGZIPInputStream(
        										new FileInputStream(fastqFiles[b+indexStartOfRead2]))));
	                    
	              System.out.println(fastqFiles[b].getName() + "---" + fastqFiles[b+indexStartOfRead2].getName()+"\n");      
	                }else{
	                    br1=new BufferedReader(new FileReader(fastqFiles[b]),65536);
	                    br2=new BufferedReader(new FileReader(fastqFiles[b+indexStartOfRead2]),65536);
	                }
	                String sequenceF="", sequenceR="", qualityScoreF="", qualityScoreR="";
	                String tempF, tempR;
	                                            
	                // clear all counters and controllers for the upcoming section
	                int currLine=0;
	                int bothGood = 0;
	                allReads = 0;
	                goodBarcodedReads = 0;
	                goodBarcodedForwardReads = 0;
	                goodBarcodedReverseReads = 0;
	                ReadBarcodeResult [] rr = new ReadBarcodeResult [2];
	                int pairCount = 0;
	                String tempSeqF=null;
	                String tempSeqR=null;
	                String tempIdF=null;
	                String tempIdR=null;
	                int hashWriteCounter=0;
	                String concatenation=null;
	                String hiseqID=null; // captures reverse ID from raw sequence file
	                int idLine = 1;  // designates lines ID information should be found starting with the first line
	                String qualOverride="BCCFFFFFHHHHHIJJIJJJJJHHHHH#########################################################################";
	                
	                System.out.println("Begin reading raw sequence files");
	                setTime();
	             //original with max controller   
	             //   while ((tempF = br1.readLine()) != null && (tempR = br2.readLine()) != null 
	             //   		&& goodBarcodedReads < maxGoodReads) {
	                // temp read through entire file
	                int rejected=0;
	                int x=0;
	                while ((tempF = br1.readLine()) != null && (tempR = br2.readLine()) != null && x<1000000) {
	                //while ((tempF = br1.readLine()) != null && (tempR = br2.readLine()) != null) {
	                	currLine++;
	                	
	                    try{
	                        //The quality score is every 4th line; the sequence is every 4th line starting from the 2nd.
	                        if(currLine==idLine){
	                        	hiseqID = tempR;
	                        	//System.out.println(currLine);
	                        	idLine+=4;
	                        }else if((currLine+2)%4==0){
	                        	//System.out.println(tempR);
	                            sequenceF = tempF;
	                            sequenceR = tempR;
	                        }else if(currLine%4==0){
	                            qualityScoreF = tempF;
	                            qualityScoreR = tempR;
	                            allReads += 2;
	                            //After quality score is read, decode barcode using the current sequence & quality  score
	                            rr[0] = thePBR[0].parseReadIntoTagAndTaxa(sequenceF, qualityScoreF, true, 0,64);
	                            rr[1] = thePBR[1].forceParseReadIntoTagAndTaxa(sequenceR, qualityScoreR, false, 0,64);
	                            if (rr[0] != null && rr[1] !=null){
	                                goodBarcodedReads+=2;
	                                //goodBarcodedForwardReads++;
	                                //goodBarcodedReverseReads++;
	                                bothGood++;
	                                //System.out.println(qualityScoreR);
	                                tempSeqF=rr[0].toString().substring(0,64);
	                                tempIdF = rr[0].toString().substring(65);
	                                tempSeqR=rr[1].paddedSequence;  //correctly handles Ns present in sequence
	                                tempIdR = rr[1].toString().substring(65);
	                                
	                                concatenation=stitch(tempSeqF, tempSeqR, tempIdF, tempIdR);
	                                
	                                //Check if sequence is part of HashMap
	                                if(pairCount.containsKey(concatenation)){
	                                	// get occurences, increment it, set new value
	                                	pairCount.put(concatenation, pairCount.get(concatenation)+1);
	                                }else{
	                                	// add first occurence
	                                	pairCount.put(concatenation, 1);
	                                }
	                                
	                                pairCount = pairCount.size();
	                                x++;
	                              
	                            }
	                            else if (rr[0] != null){
	                                goodBarcodedForwardReads++;
	                                //x++;
	                     /*           String bb = sequenceR.substring(0,11);
	                                //String bb = sequenceR;
	                                int n = 0;
	                                for(int i=0;i<9;i++){
	                                	if(bb.charAt(i)=='N'){
	                                		n++;
	                                	}
	                                }
	                                if(n>1){
	                                	rejected++;
	                                }else{
	                                	badBarcodeRead2.add(bb);
	                                	//badBarcodeRead2.add(x+"\t"+bb);
	                                }
	                     
	                     */     }
	                            else if (rr[1] != null){
	                                goodBarcodedReverseReads++;
	                                //System.out.println(rr[1].paddedSequence);
	                            }
	                            if (allReads % 10000000 == 0) {
	                            	reportStats(bothGood, goodBarcodedForwardReads, goodBarcodedReverseReads, 
	                            			goodBarcodedReads, allReads, pairCount);
	                            	printListToFile(outputDir+File.separator+"List_bad_read2_barcodes.txt",badBarcodeRead2);
	                     //       	badBarcodeRead2.clear();
	                     //       	System.out.println("REJECTED: >1 N in first 9 bases: "+rejected);
	                            }
	                            
	                            if(allReads % 40000000 ==0){
	                            	processPairsCounted(countFileNames[b],countFileNames[b+indexStartOfRead2],outputDir,
	                            			hashWriteCounter,pairCount, hashFileNames);
	                            	hashWriteCounter++;
	                            }
	                            	
	                        }
	                    }catch(NullPointerException e){
	                        System.out.println("Unable to correctly parse the sequence and "
	                        + "quality score from fastq file.  Your fastq file may have been corrupted.");
	                        System.exit(0);
	                    }
	                    
	                } 
	                
	                //System.out.println("REJECTED: >1 N in first 9 bases: "+rejected);
	           //     printListToFile(outputDir+File.separator+"List_bad_read2_barcodes.txt",badBarcodeRead2);
               // 	badBarcodeRead2.clear();
	                
	                processPairsCounted(countFileNames[b],countFileNames[b+indexStartOfRead2],outputDir,
                			hashWriteCounter,pairCount, hashFileNames);
	                reportTime(now);
	 
                reportStats(bothGood, goodBarcodedForwardReads, goodBarcodedReverseReads, 
            			goodBarcodedReads, allReads, pairCount);
                br1.close();
                br2.close();
              
                fileNum++;
            
		        } catch(Exception e) {
		            System.out.println("Catch testBasicPipeline c="+goodBarcodedReads+" e="+e);
		            e.printStackTrace();
		            System.out.println("Finished reading "+(fileNum+1)+" of "+fastqFiles.length+" sequence files.");
        			}
			}
            combineHashOutputFiles(hashFileNames, outputDir);
        	printIds(outputDir);
        	printListToFile(outputDir+File.separator+"Summary_Stats.txt",summaryOutputs);
        	
    }
    
    
    private static void summarizeCounts(String line){
    	summaryOutputs.add(line);
    }
    
    private static String modifyKey2File(String name){
    	String newName=name+"_mod";
    	String keyContents;
    	ArrayList <String> modified = new ArrayList<String>();
    	int currentLine = 0;
    	int barcodeLength = 0;
    	
    	copyOriginalKeyFile(name, newName);
    	
    	try{
    		BufferedReader inputKey = new BufferedReader(new InputStreamReader(
    				new FileInputStream(name)));
    	
	    	while((keyContents = inputKey.readLine()) != null){
	    		currentLine++;
	    		if(currentLine==1){
	    			System.out.println("Begin modifications to key file");
	    		}else{
	    			String splitKeyLine[] = keyContents.split("\t");
	    			barcodeLength=(splitKeyLine[2].length());
	    			
	    			deletion(splitKeyLine, barcodeLength, modified);
	    			substitution(splitKeyLine, barcodeLength, modified);
	    					
	    				
	    		}

	    	}		
	    	inputKey.close();
	    	}catch(Exception e) {
            System.out.println("modifyKey2File: e="+e);
            e.printStackTrace();
			}
    	
    	printListToFile(newName, modified);
    	
    	return newName;
    }
    
    private static void copyOriginalKeyFile(String name, String newName){
    	String temp;
    	ArrayList <String> hold = new ArrayList<String>();
    	int size=0;
    	
    	try{
    		BufferedReader source = new BufferedReader(new InputStreamReader(
    				new FileInputStream(name)));
    		while((temp=source.readLine())!=null){
    			hold.add(temp);
    		}
    		source.close();
    		}catch(Exception e) {
            System.out.println("copyOriginalKeyFile: e="+e);
            e.printStackTrace();
			}
    	
    	printListToFile(newName, hold);
    	
    }
    
    private static void deletion(String [] lineElement, int length, ArrayList<String>list){
    	String tempBarcode;
    	String barcode = lineElement[2];
    	
    	for(int i=0; i<length; i++){
			if(i==0){
				list.add(lineElement[0]+"\t"+lineElement[1]+"\t"
						+barcode.substring(1,length)+"\t"+lineElement[3]+"_mod\t"
						+lineElement[4]+"\t"+lineElement[5]+"\t"+lineElement[6]+"\t"
						+lineElement[7]);
			}
			else{
				tempBarcode = barcode.substring(0,i)+barcode.substring(i+1,length);
				if(list.contains(tempBarcode)){
					continue;
				}else{
					list.add(lineElement[0]+"\t"+lineElement[1]+"\t"
							+tempBarcode+"\t"+lineElement[3]+"_mod\t"
							+lineElement[4]+"\t"+lineElement[5]+"\t"+lineElement[6]+"\t"
							+lineElement[7]);
				}
			}
    	}
    	
    	for(int i=0; i<length-1; i++){
			if(i==0){
				list.add(lineElement[0]+"\t"+lineElement[1]+"\t"
						+barcode.substring(2,length)+"\t"+lineElement[3]+"_mod\t"
						+lineElement[4]+"\t"+lineElement[5]+"\t"+lineElement[6]+"\t"
						+lineElement[7]);
			}
			else{
				tempBarcode = barcode.substring(0,i)+barcode.substring(i+2,length);
				if(list.contains(tempBarcode)){
					continue;
				}else{
					list.add(lineElement[0]+"\t"+lineElement[1]+"\t"
							+tempBarcode+"\t"+lineElement[3]+"_mod\t"
							+lineElement[4]+"\t"+lineElement[5]+"\t"+lineElement[6]+"\t"
							+lineElement[7]);
				}
			}
    	}
    }
    
    private static void substitution(String [] lineElement, int length, ArrayList<String>list){
    	String originalBarcode = lineElement[2];
    	char [] nucleotides = {'A','C','G','T','N'};
    	char [] barcodeAsArray = originalBarcode.toCharArray();
    	String moddedBarcode;
    	String newLine;
    	
    	for(int i=0;i<length;i++){
    		for(int k=0;k<nucleotides.length;k++){
    			barcodeAsArray = originalBarcode.toCharArray();
    			barcodeAsArray[i]=nucleotides[k];
    			moddedBarcode=new String(barcodeAsArray);
    			newLine = lineElement[0]+"\t"+lineElement[1]+"\t"
						+moddedBarcode+"\t"+lineElement[3]+"_mod\t"
						+lineElement[4]+"\t"+lineElement[5]+"\t"+lineElement[6]+"\t"
						+lineElement[7];
    				
    			if(list.contains(newLine) || moddedBarcode.equals(originalBarcode)){
    				continue;
    			}else{
    				list.add(lineElement[0]+"\t"+lineElement[1]+"\t"
						+moddedBarcode+"\t"+lineElement[3]+"_mod\t"
						+lineElement[4]+"\t"+lineElement[5]+"\t"+lineElement[6]+"\t"
						+lineElement[7]);
    			}
    			
    		}
    	}
    
    }
    
    private static void printListToFile(String name, ArrayList<String> list){
    	int size=0;
    	try{
    		PrintWriter out = new PrintWriter(new BufferedWriter(
    				new FileWriter(name, true)));

    		list.trimToSize();
    		size=list.size();
    		for(int i=0;i<size;i++){
    			out.println(list.get(i).toString());
    		}
    		out.close();
    	}catch (IOException e) {
            System.out.println(e.getMessage());
            }
    }
    
    /**
     * Writes out the contents of a HashMap to a unique file.  The HashMap contains sequences
     * and identification information.  Due to the large amounts of data that is 
     * being process, at the end of the method, the HashMap is cleared
     * to free system resources.
     * @param fileName1 name of the first source data file
     * @param fileName2 name of the second source data file
     * @param dir is the output directory
     * @param order is an integer that is counting the number of files that have been written and is used
     * to help add a uniqueness to the filename
     * @param stored the HashMap that contains the sequences and identification information
     * @param nameLog is the list containing file names processed by this file
     */
    private static void processPairsCounted(String fileName1, String fileName2,String dir,
    		int order, TreeMap<String, Integer> stored, ArrayList <String> nameLog){
    	
    	try {
        	long timeTemp = System.currentTimeMillis();
        	String hashOutName = Integer.toString(order)+"-"+fileName1+"-"+fileName2+".txt";
        	
            PrintWriter out = new PrintWriter(
            		new BufferedWriter(
            				new FileWriter(
            						dir+File.separator+hashOutName, true)));
        
        	nameLog.add(hashOutName);
        	
        	
            for(String h: stored.keySet()){
            	String key = h.toString();
            	String value = stored.get(h).toString();
            	out.println(key+ "\t" + value);
            }
            out.close();		// close PrintWriter
            //reporter
            System.out.println("The number of lines added to " + hashOutName +" is " +stored.size());
            
            stored.clear();  // force memory clear before continuing
            reportTime(timeTemp);
        }catch (IOException e) {
        	System.out.println(e.getMessage());
        	;
        }
    }
    
    /**
     * Checks for an even number of files
     * @param numberOfPairs is the number of pairs
     */
    private static void checkForPairs(int numberOfPairs){
    	if (numberOfPairs % 2 !=0){
    		System.out.println("There are an odd number of files so there won't be correct pairing"); 
    		System.out.println("The number of files detected was "+numberOfPairs); 
        }
    }
    
    /**
     * Reporter method that prints stats to the console
     * @param both - a counter that keeps track of the number of times both sequences register as good reads
     * @param forward - a counter that keeps track of the number of times only the forward sequence registers as a good read
     * @param reverse -  a counter that keeps track of the number of times only the reverse sequence registers as a good read
     * @param allGood - a counter that keeps the current total value of all lines read so far
     * @param totalReads - the total number of lines, good or bad, that the program has read so far
     */
    private static void reportStats(int both, int forward, int reverse, int allGood, int totalReads, int listCount){
    	
    	float percentAll = 100*((float)allGood/totalReads);
    	float percentBoth = 100*((float)both/allGood);
    	float percentForward = 100*((float)forward/allGood);
    	float percentReverse = 100*((float)reverse/allGood);
    	float percentUnique = 100*((float)listCount/allGood);
    	DecimalFormat formatter = new DecimalFormat("00.0");
    	String forSummary = Integer.toString(totalReads)+"\t"+Integer.toString(forward)+"\t"+
    			Integer.toString(reverse)+"\t"+Integer.toString(both);
    	
    	System.out.println("Total Reads:" + totalReads);
    	System.out.println("The number of all good lines encountered so far is "+allGood+" ("+formatter.format(percentAll)+
    			"% of total reads)");
    	System.out.println("The number of only good forward reads in this file is: "+forward+" ("+formatter.format(percentForward)+
    			"% of all good lines read)");
    	System.out.println("The number of only good reverse reads in this file is: "+reverse+" ("+formatter.format(percentReverse)+
    			"% of all good lines read)");
    	System.out.println("The number of good forward and reverse reads in this file is: "+both+" ("+formatter.format(percentBoth)+
    			"% of all good lines read)");
    	System.out.println("The number of unique pairs of good forward and reverse sequence in this file is at most "+listCount+
    			" ("+formatter.format(percentUnique)+"% of all good lines read)\n");
    	
    	summarizeCounts(forSummary);
    }
    
    /**
     * Prints a message to the console indicating there is a problem with the file name structure
     * @param filename - name of file that failed parsing tests
     */
    private static void printParsingError(){
    	
    	System.out.println("Error in parsing file name:");
		System.out.println("   The filename does not contain a 5 underscore-delimited value.");
		System.out.println("   Expected: code_flowcell_s_lane_fastq.txt.gz");
		System.out.println("OR There is already a file in the ouput folder of the same name");
    }
    
    /**
     * Take both sequence and id information from the raw sequence file and concatenate them
     * @param forward - the forward sequence
     * @param reverse - the reverse sequence
     * @param idF - the identification and summary information from the forward sequence
     * @param idR - the identification and summary information from the reverse sequence
     * @return - a string concatenation of the sequences and id information with formatting elements
     */
    private static String stitch(String forward, String reverse, String idF, String idR){
    	//String tempStitch=forward+","+reverse+"\t"+idF+"\t"+idR;
    	String result;
    	String tag = forward+reverse;
    	String id1[]=idF.split(":");
    	String id2[]=idR.split(":");
    	String flowcell = id1[2];
    	String barcodeIDs = id1[1]+":"+id2[1];
    	String lane = id1[3];
    	String tab="\t";
    	
    	result = tag+tab+barcodeIDs+tab+flowcell+tab+lane;
 
    	return result;
    }
    
    /**
     * Combines the HashMap output text files into one summary text file.  The expectation is that
     * incoming files will be tab delimited.
     * @param names - names of the files to process
     * @param directoryInfo - system directory information
     */
    private static void combineHashOutputFiles(ArrayList<String> names, String directoryInfo){
    	
    	String seq;	// sequence
    	String barcodeCombo;	// combination of barcodes that generated a tag
    	String flowcell;	// flowcell
    	String lane;	// lane
    	int numberOfCounts=0; // number of times a sequence has been observed
    	String [] arrayNames = names.toArray(new String[names.size()]); // copy list of files to process
    	
    	TreeMap <String, TreeMap<String, String[]> tagInfo>= new TreeMap<String, TreeMap<String, String[5]>>();
    	
    	// the following HashMap allows for one key (the sequence string) to hold 
    	// many values (id infomation) via the ArrayList
    	//HashMap <String, ArrayList<String>> hma = new HashMap<String, ArrayList<String>>();
    	ArrayList <String> tempArrayList = new ArrayList<String>(); // holds information before addition to HashMap
  		int tempCount = 0;  //reporter variable
    	long tempTime=0;	//reporter variable
    	
    	//Reporter
		tempTime = System.currentTimeMillis();
		System.out.println("\ncombineHashOutputFiles: Gathering information from individual files");
    	
		// process all the files
    	for(int i=0; i<arrayNames.length; i++){
    		
    		BufferedReader hbr=null;
    		String location = directoryInfo+File.separator+arrayNames[i];    		
    		
    		try{
    			String lineRead;
    			hbr=new BufferedReader(new FileReader(location));
    			
    			while ((lineRead = hbr.readLine()) != null){
	    			 
	    			 String splitLine[] = lineRead.split("\t"); // parse infomation from incoming file
	    			 seq = splitLine[0];		
	    			 barcodeCombo = splitLine[1]; 
	    			 logIds(barcodeCombo);
	    			 ids = r1Ids[0]+":"+r1Ids[2]+":"+r1Ids[3]+":"+r1Ids[1]+":"+r2Ids[1]; 
	    			 numberOfCounts=Integer.parseInt(splitline[3]); // copy the sequence count
	    			 
	    				 //Check if sequence is part of HashMap
	    				 if(hma.containsKey(allSeq) ){
	 		            	// get occurences, increment counter, set new value
	 		            	tempArrayList = new ArrayList(hma.get(allSeq));
	 		            	tempArrayList.set(0, Integer.toString(Integer.parseInt(tempArrayList.get(0))+numberOfCounts));
	 		            	tempArrayList.set(1,tempArrayList.get(1)+"\t"+stripRedundant(ids));
	 		            	hma.put(allSeq, tempArrayList);		
	 		            }else{
	 		            	// add new unique line
	 		            	tempArrayList = new ArrayList<String>();
	 		            	tempArrayList.add(Integer.toString(numberOfCounts));
	 		            	tempArrayList.add(ids);
	 		            	hma.put(allSeq, tempArrayList);
	 		            }
	    		 }
    		}catch(IOException io) { 
                System.out.println(io.getMessage());
            }
    		System.out.println(arrayNames[i]+" processed");
    	}
    	// Reporter
    	reportTime(tempTime);    	
    	tempTime = System.currentTimeMillis(); 
    	System.out.println("\ncombineHashOutputFiles: Start writing to output file");
    	System.out.println("The number of lines to be sent to the output file is " + hma.size());
    	
    	try {
    		ArrayList value = new ArrayList();
        	PrintWriter out = new PrintWriter(
            		new BufferedWriter(
            				new FileWriter(
            						directoryInfo+File.separator+"Paired_End_Tags_Info.txt", true)));
        tempCount=0;
        	for(String h: hma.keySet()){
            	String key = h.toString();
            	value= new ArrayList(hma.get(h));
            	out.println(key+ "\t" + value.get(0).toString()+"\t"+value.get(1).toString());
            	tempCount++;
            	// send update to console so user gets status update
            	if(tempCount%1000000 == 0){
            	System.out.println(tempCount+"lines written to file");
            	}
            }
            out.close();		// close PrintWriter
            hma.clear();  //force memory release 
        }catch (IOException e) {
        	 System.out.println(e.getMessage());
        }
    	reportTime(tempTime);
    }
    
    /**
     * Sets a static class variable to current system time.
     */
    private static void setTime(){
    	now = System.currentTimeMillis();
    }
    
    /**
     * A reporter method that takes in a millisecond start time and processes
     * it to a friendlier hour, minute, sec output.  Minutes and seconds are rounded up
     * to the nearest integer.  Method should only be called in code once a process is completed 
     * since it calls on the current system time to act as the end of a process.
     * @param start - a time variable in milliseconds
     */
    private static void reportTime(long start){
    	long end = System.currentTimeMillis();
    	long time = end-start; //milliseconds
    	int sec = 0;
    	int mins = 0;
    	int hours = 0;
    	
    	if(time>=1000){
    		if(time>=60000){
    			if(time>=3600000){
    				hours = (int)time/3600000;
    				mins = Math.round((time % 3600000)*60);
    				System.out.println("Process completed in " +hours + " hours and "+ mins + " minutes\n");
    			}else{
    				mins = (int)time / 60000;
    				System.out.println("Process completed in " +mins + " minutes\n");
    			}
    		}else{
    			sec = Math.round((time / 1000));
				System.out.println("Process completed in " +sec + " seconds\n");
    		}
    	}else{
    		System.out.println("Process completed in " +time + " milliseconds\n");
    	}
    }
    
    private static String stripRedundant(String longID){
    	String tempID[]=longID.split(":");
    	String newID = tempID[2]+":"+tempID[3]+":"+tempID[4];
    	return newID;
    }
    
    /**
     * 
     * @param id
     */
    private static void logIds(String id){
    	
    	if(barcodePairs.containsKey(id)){
    		barcodePairs.put(id, barcodePairs.get(id)+1);
    	}else{
    		barcodePairs.put(id, 1);
    	}
    }
    
    
    private static void printIds(String directoryInfo){
    	long tempTime = System.currentTimeMillis(); 
    	System.out.println("\nprintIds: Start writing to output file");
    	System.out.println("The number of lines to be sent to the output file is " + barcodePairs.size());
    	
    	try {
    		ArrayList value = new ArrayList();
        	PrintWriter out = new PrintWriter(
            		new BufferedWriter(
            				new FileWriter(
            						directoryInfo+File.separator+"Paired_End__ID_Info.txt", true)));
        int tempCount=0;
        	for(String h: barcodePairs.keySet()){
            	String key = h.toString();
            	String val = barcodePairs.get(h).toString();
            	out.println(key+ "\t" + val);
            	tempCount++;
            	// send update to console so user gets status update
            	if(tempCount%1000000 == 0){
            	System.out.println(tempCount+"lines written to file");
            	}
            }
            out.close();		// close PrintWriter
        }catch (IOException e) {
        	 System.out.println(e.getMessage());
        }
    	reportTime(tempTime);
    }
    
    @Override
    public ImageIcon getIcon(){
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
