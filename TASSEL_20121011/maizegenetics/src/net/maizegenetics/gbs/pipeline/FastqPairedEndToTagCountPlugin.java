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
import java.util.HashMap;
import java.util.Set;
import java.util.Map;

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
    private ArgsEngine engine = null;
    private Logger logger = Logger.getLogger(FastqPairedEndToTagCountPlugin.class);
    String directoryName=null;
    String keyfile=null;
    String enzyme = null;
    //int maxGoodReads = 200000000;
    int maxGoodReads = 2000000;
    int minCount =1;
    String outputDir=null;

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
        HashMap <String, Integer> pairCount = new HashMap<String, Integer>(); // stores paired sequences to write to file
        ArrayList <String> hashFileNames = new ArrayList<String>(); // stores names of files resulting from HashMap output
        
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
/*//DEBUG print all array contents
for(int left=0;left<2;left++){
	for(int mid=0;mid<indexStartOfRead2;mid++){
		for(int right=0;right<5;right++){
			System.out.println("fileReadInfo FOR ["+left+"]"+"["+mid+"]"+"["+right+"]"+"IS: "+fileReadInfo[left][mid][right]); //DEBUG 
		}
	}
}
*/
//handle keyfiles and enzymes
//2 arrays for manually inputing multiple enzymes and keys for testing
System.out.println("OLD Key file is:"+ keyFileS);            
System.out.println("OLD enzyme is:"+ enzyme);
String[] hcEnzyme={"PstI-MspI","MspI-PstI"};
String[] hcKeyFiles={"GBS.key","GBS2.key"};
		
 			TagCountMutable [] theTC=new TagCountMutable [2];
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
	                            
	                try{
	                    theTC[0] = new TagCountMutable(2, maxGoodReads);
	                    theTC[1] = new TagCountMutable(2, maxGoodReads);
	                }catch(OutOfMemoryError e){
	                    System.out.println(
	                        "Your system doesn't have enough memory to store the number of sequences"+
	                        " you specified.  Try using a smaller value for the minimum number of reads."
	                    );
	                }
	                
	                // clear all counters and controllers for the upcoming section
	                int currLine=0;
	                int bothGood = 0;
	                allReads = 0;
	                goodBarcodedReads = 0;
	                goodBarcodedForwardReads = 0;
	                goodBarcodedReverseReads = 0;
	                ReadBarcodeResult [] rr = new ReadBarcodeResult [2];
	                int hashCount = 0;
	                String tempSeqF=null;
	                String tempSeqR=null;
	                String tempIdF=null;
	                String tempIdR=null;
	                
	                while ((tempF = br1.readLine()) != null && (tempR = br2.readLine()) != null 
	                		&& goodBarcodedReads < maxGoodReads) {
	                	currLine++;
	                    try{
	                        //The quality score is every 4th line; the sequence is every 4th line starting from the 2nd.
	                        if((currLine+2)%4==0){
	                            sequenceF = tempF;
	                            sequenceR = tempR;
	                        }else if(currLine%4==0){
	                            qualityScoreF = tempF;
	                            qualityScoreR = tempR;
	                            allReads += 2;
	                            //After quality score is read, decode barcode using the current sequence & quality  score
	                            rr[0] = thePBR[0].parseReadIntoTagAndTaxa(sequenceF, qualityScoreF, true, 0,64);
	                            rr[1] = thePBR[1].parseReadIntoTagAndTaxa(sequenceR, qualityScoreR, true, 0,64);
	                            if (rr[0] != null && rr[1] !=null){
	                                goodBarcodedReads+=2;
	                                goodBarcodedForwardReads++;
	                                goodBarcodedReverseReads++;
	                                bothGood++;
	                                
	                                theTC[0].addReadCount(rr[0].getRead(), rr[0].getLength(), 1);
	                                theTC[1].addReadCount(rr[1].getRead(), rr[1].getLength(), 1);
	                                
	                                tempSeqF=rr[0].toString().substring(0,64);
	                                tempIdF = rr[0].toString().substring(65);
	                                tempSeqR=rr[1].toString().substring(0,64);
	                                tempIdR = rr[1].toString().substring(65);
	                                
	                                String concatenation=stitch(tempSeqF, tempSeqR, tempIdF, tempIdR);
	                                
	                                //Check if sequence is part of HashMap
	                                if(pairCount.containsKey(concatenation)){
	                                	// get occurences, increment it, set new value
	                                	pairCount.put(concatenation, pairCount.get(concatenation)+1);
	                                }else{
	                                	// add first occurence
	                                	pairCount.put(concatenation, 1);
	                                }
	                                
	                                hashCount = pairCount.size();
	                              
	                            }
	                            else if (rr[0] != null){
	                                goodBarcodedReads++;
	                                goodBarcodedForwardReads++;
	                               theTC[0].addReadCount(rr[0].getRead(), rr[0].getLength(), 1);
	                            }
	                            else if (rr[1] != null){
	                                goodBarcodedReads++;
	                                goodBarcodedReverseReads++;
	                               theTC[1].addReadCount(rr[1].getRead(), rr[1].getLength(), 1);
	                            }
	                            if (allReads % 10000000 == 0) {
	                            	reportStats(bothGood, goodBarcodedForwardReads, goodBarcodedReverseReads, 
	                            			goodBarcodedReads, allReads, hashCount);
	                            }
	                        }
	                    }catch(NullPointerException e){
	                        System.out.println("Unable to correctly parse the sequence and "
	                        + "quality score from fastq file.  Your fastq file may have been corrupted.");
	                        System.exit(0);
	                    }
	                    
	                } 
	                try {
	                	String hashOutName = countFileNames[b]+"-"+countFileNames[b+indexStartOfRead2]+".txt";
	                	
	                    PrintWriter out = new PrintWriter(
	                    		new BufferedWriter(
	                    				new FileWriter(
	                    						outputDir+File.separator+hashOutName, true)));
	                
	                	hashFileNames.add(hashOutName);
	                	
	                	
		                for(String h: pairCount.keySet()){
		                	String key = h.toString();
		                	String value = pairCount.get(h).toString();
		                	out.println(key+ "\t" + value);
		                }
		                out.close();		// close PrintWriter
		                pairCount.clear();  //force memory release before looping back through
	                }catch (IOException e) {
	                	System.out.println(e.getMessage());
	                	;
	                }
	                
                reportStats(bothGood, goodBarcodedForwardReads, goodBarcodedReverseReads, 
            			goodBarcodedReads, allReads, hashCount);
                System.out.println("Timing process (sorting, collapsing, and writing TagCount to file).");
                timePoint1 = System.currentTimeMillis();
                theTC[0].collapseCounts();
                theTC[1].collapseCounts();
                theTC[0].writeTagCountFile(outputDir+File.separator+countFileNames[b], FilePacking.Bit, minCount);
                theTC[1].writeTagCountFile(outputDir+File.separator+countFileNames[b+indexStartOfRead2], FilePacking.Bit, minCount);
                System.out.println("Process took " + (System.currentTimeMillis() - timePoint1) + " milliseconds.");
                br1.close();
                br2.close();
                //force memory release before looping back through
                theTC[0]=null;
                theTC[1]=null;
                
                fileNum++;
            
		        } catch(Exception e) {
		            System.out.println("Catch testBasicPipeline c="+goodBarcodedReads+" e="+e);
		            e.printStackTrace();
		            System.out.println("Finished reading "+(fileNum+1)+" of "+fastqFiles.length+" sequence files.");
        			}
			}
            combineHashOutputFiles(hashFileNames, outputDir);
        	
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
    	DecimalFormat formatter = new DecimalFormat("00.0");
    	
    	System.out.println("Total Reads:" + totalReads);
    	System.out.println("The number of good lines encountered so far is "+allGood+" (~"+formatter.format(percentAll)+"%)");
    	System.out.println("The number of good forward and reverse reads is: "+both+" (~"+formatter.format(percentBoth)+"%)");
    	System.out.println("The number of total good forward reads is: "+forward+" (~"+formatter.format(percentForward)+"%)");
    	System.out.println("The number of total good reverse reads is: "+reverse+" (~"+formatter.format(percentReverse)+"%)");
    	System.out.println("HashMap size is "+listCount+"\n");
    	System.out.println("Percentages are only an approximation\n");  
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
    
    private static String stitch(String forward, String reverse, String idF, String idR){
    	String tempStitch=forward+","+reverse+"\t"+idF+"\t"+idR;
    	return tempStitch;
    }
    
    private static void combineHashOutputFiles(ArrayList<String> names, String directoryInfo){
    	
    	String allSeq;
    	String ids;
    	int numberOfCounts=0;
    	String [] arrayNames = names.toArray(new String[names.size()]);
    	HashMap <String, ArrayList<String>> hma = new HashMap<String, ArrayList<String>>();
    	ArrayList <String> tempArrayList = new ArrayList<String>();
    	
    	for(int i=0; i<arrayNames.length; i++){
    		
    		BufferedReader hbr=null;
    		String location = directoryInfo+File.separator+arrayNames[i];
    		
    		if(i!=0){
    			//reporter
    			System.out.println("The first file has been processed");
    			System.out.println("The number of unique sequences at this point "+ hma.size());
    		}
    		
    		try{
    			String lineRead;
    			hbr=new BufferedReader(new FileReader(location));
    			while ((lineRead = hbr.readLine()) != null){
	    			 
	    			 String splitline[] = lineRead.split("\t");
	    			 allSeq = splitline[0];
	    			 ids = splitline[1]+"\t"+splitline[2];
	    			 numberOfCounts=Integer.parseInt(splitline[3]);
	    			 
	    			 /*
	    			  * If this is the first file, add everything as is since uniqueness
	    			  * in each separate files is already determined.  This will cut down on
	    			  * needless comparisons
	    			  */
	    			 if(i=0){
	    				// add everything from first file
		            	tempArrayList.add(Integer.toString(numberOfCounts));
		            	tempArrayList.add(ids);
		            	hma.put(allSeq, tempArrayList);
	    			 }
	    			 else{
	    				 //Check if sequence is part of HashMap
	    				 if(hma.containsKey(allSeq) ){
	 		            	// get occurences, increment counter, set new value
	 		            	tempArrayList = hma.get(allSeq);
	 		            	tempArrayList.set(0, tempArrayList.get(0)+numberOfCounts);
	 		            	tempArrayList.set(1,tempArrayList.get(1)+"\t"+ids);
	 		            	hma.put(allSeq, tempArrayList);		            	
	 		            }else{
	 		            	// add new unique line
	 		            	tempArrayList.add(Integer.toString(numberOfCounts));
	 		            	tempArrayList.add(ids);
	 		            	hma.put(allSeq, tempArrayList);
	 		            }
	    			 }
	    		 }
    		}catch(IOException io) { 
                System.out.println(io.getMessage());
            }	 
    	}
    	//Reporter
    	System.out.println("Start writing to output file");
    	System.out.println("The number of lines to be sent to the output file is " + hma.size());
    	try {
        	PrintWriter out = new PrintWriter(
            		new BufferedWriter(
            				new FileWriter(
            						directoryInfo+File.separator+"Paired_End_Info.txt", true)));
        
        	for(String h: hma.keySet()){
            	String key = h.toString();
            	ArrayList value = new ArrayList(hma.get(h));
            	out.println(key+ "\t" + value.get(0)+"\t"+value.get(1));
            }
            out.close();		// close PrintWriter
            hma.clear();  //force memory release 
        }catch (IOException e) {
        	 System.out.println(e.getMessage());
        	;
        }
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
