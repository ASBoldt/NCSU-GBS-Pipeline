/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.tagdist;
import net.maizegenetics.gbs.util.OpenBitSet;
import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.FileReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import net.maizegenetics.genome.BaseEncoder;
import net.maizegenetics.gbs.homology.ParseBarcodeRead;
import net.maizegenetics.gbs.homology.ReadBarcodeResult;

/**
 * This class can read a CBSU TagMapFile into the TagsOnPhysicalMap data
 * structure.
 *
 * @author harriman
 *
 */
public class TagsByTaxaCBSUAdapter {
    boolean cleanCutSites=true;
    static String filePrefix="/Users/edbuckler/SolexaAnal/GBS/test/";
    private int[] origOrder;
    private long[][] sortTags;
    private int tagLengthInLong;

    public static void main(String[] args){

        String infile = "/media/Data/Chih-Wei/rice_chromosomesrenamed.tags";
        String taxonNameFile  = new String("/media/Data/Chih-Wei/RIL_A_B.key");
        String outfile = new String("//media/Data/Chih-Wei/output.tbt");

        TagsByTaxaCBSUAdapter test = new TagsByTaxaCBSUAdapter(infile, taxonNameFile, outfile);
        test.convertCBSUFile(new File(infile), new File(taxonNameFile), new File(outfile));
 
   }

    public TagsByTaxaCBSUAdapter(String infile, String lineNameFile, String outfile) {
        int columns=-1;
        try {
	    BufferedReader br = new BufferedReader(new FileReader(infile), 65536);
            columns=br.readLine().split("\t").length;
            br.close();
        } catch (IOException e) {
	    System.out.println("Catch in reading alignment file: " + e);
	}
        try {
	    BufferedReader br = new BufferedReader(new FileReader(lineNameFile), 65536);
            columns=br.readLine().split("\t").length;
            br.close();
        } catch (IOException e) {
	    System.out.println("Catch in reading line name file: " + e);
	}

//	convertCBSUFile(new File(infile), new File(lineNameFile), new File(outfile));
//	File[] testFiles = DirectoryCrawler.listFiles(".*xa.*", "c:/file_visitor_test");
//	binaryFileTestClass.readBitDistFile(new File("c:/file_visitor_test/outputb.bin"));
//	binaryFileTestClass.writeDistFile(new File("c:/file_visitor_test/outputb.txt"), TagsByTaxa.FilePacking.Text, 0);
//	binaryFileTestClass.writeDistFile(new File("outfile2.txt"), TagsByTaxa.FilePacking.Bit, 0);
    }

//    public void mergeBitDistFiles(File[] inputFiles){
// 	TagsByTaxaBit mergedFile = new TagsByTaxaBit(); //Use void constructor because we don't know anything about the output yet.
//// 	TagsByTaxaBit currFile = new TagsByTaxaBit(); //Use void constructor because we don't know anything about the output yet.
//
////	ArrayList<String> currTaxaNameList = new ArrayList<String>();
////	Hashtable<String, Boolean> taxaHash = new Hashtable<String, Boolean>();
////	Hashtable<Long, Boolean> tagHash = new Hashtable<Long, Boolean>();
////	HashSet currTaxaNameList = new HashSet();
//	HashSet<String> mergedTaxaNameList = new HashSet<String>();
//	int totalTagNumber=0;
//	boolean taxaNamesMatch;
//
//	//Loop through files once to determine total # of tags and taxa from headers.
//	for(File inputFile: inputFiles){
//	    TagsByTaxaBit currFile = new TagsByTaxaBit(inputFile.getAbsolutePath(), FilePacking.Bit);
//	    mergedFile.tagLengthInLong = currFile.tagLengthInLong;
//	    totalTagNumber += currFile.getTagCount();
//	    for(String taxonName: currFile.taxaNames){
//		mergedTaxaNameList.add(taxonName);	//Store name in a master list
//	    }
//	}
//
//	mergedFile.taxaNum = mergedTaxaNameList.size();
//	mergedFile.initMatrices(mergedFile.taxaNum, totalTagNumber); //Now we know how big it should be
//	mergedTaxaNameList.toArray(mergedFile.taxaNames);  //Dump master taxon list into output file's taxon list
//	int totalTagsWritten = 0;
//
//	//Loop through again to merge all records
//        for(File inputFile: inputFiles){
//		TagsByTaxaBit currFile = new TagsByTaxaBit(inputFile.getAbsolutePath(), FilePacking.Bit);
//		OpenBitSet currBitSet = currFile.getTaxaReadBitsForTag(totalTagNumber); //Create bitset for current file
//
//		/*  Begin loop of PAIN */
//		for (int tag = 0; tag < currFile.getTagCount(); tag++) {
//		    for (int j = 0; j < mergedFile.tagLengthInLong; j++) { //Copy sequence from current file to merged file
//			mergedFile.tags[j][tag] = currFile.tags[j][tag];
//		    }
//		    mergedFile.tagLength[tag]=currFile.tagLength[tag];	//Copy read length
//
//		    //If current file contains this taxon, put a reference to it in the master list
//		    //(by setting the corresponding value in the taxon hash to "true").
//		    for (int currTaxon = 0; currTaxon < currFile.taxaNames.length; currTaxon++) {
//			if(currFile.getReadCountForTagTaxon(tag, currTaxon)>0){
//                            mergedFile.setReadCountForTagTaxon(tag, currTaxon, 1);	//increment corresponding read count.
//			}
//		    }
//
////		    for(int currTaxon=0; currTaxon<mergedFile.taxaNum; currTaxon++){
////			if(mergedTaxaNameList.contains(currFile.taxaNames[currTaxon])){		//If master list contains entry for current read,
////                        }
////		    }
//		    //TODO: Write out tag distribution
//		    totalTagsWritten++;
//		}
//		/*End loop */
//            }
// 	System.out.println(totalTagsWritten + " tags written to output file.");
//	System.out.println(mergedFile.taxaNum + " taxa represented in output file.");
// 	mergedFile.writeDistFile(new File("c:/merged.bin"), FilePacking.Bit, 0);
//    }

    public TagsByTaxaCBSUAdapter(String infile, String outfile) {
        TagsByTaxaBitFileMap tbtb2=new TagsByTaxaBitFileMap(infile);
        tagLengthInLong=tbtb2.getTagSizeInLong();
        sortTags=new long[tbtb2.getTagSizeInLong()][tbtb2.getTagCount()];
        origOrder=new int[tbtb2.getTagCount()];
        for (int i = 0; i < tbtb2.getTagCount(); i++) {
            long[] t=tbtb2.getTag(i);
            for (int j = 0; j < tbtb2.getTagSizeInLong(); j++) sortTags[j][i]=t[j];
            origOrder[i]=i;
        }
        Swapper swapperPos = new Swapper() {
           public void swap(int a, int b) {
              int tI;
              tI = origOrder[a]; origOrder[a] = origOrder[b]; origOrder[b] = tI;
              long tL;
              for (int i = 0; i < tagLengthInLong; i++) {
                    tL=sortTags[i][a];sortTags[i][a]=sortTags[i][b]; sortTags[i][b]=tL;
                }
           }
        };
        IntComparator compPos = new IntComparator() {
           public int compare(int a, int b) {
            for (int i = 0; i < tagLengthInLong; i++) {
                if(sortTags[i][a]<sortTags[i][b]) return -1;
                if(sortTags[i][a]>sortTags[i][b]) return 1;
            }
            return 0;
           }
        };
        System.out.println("Position index sort begin.");
        GenericSorting.quickSort(0, origOrder.length, compPos, swapperPos);
        int unique=0;
        for (int i = 1; i < tbtb2.getTagCount(); i++) {
 //           System.out.printf("%d %d %d %n",sortTags[0][i], sortTags[1][i], origOrder[i]);
            if(compPos.compare(i-1,i)!=0) unique++;
        }
        System.out.println("Unique:"+unique);
         System.out.println("Position index sort end.");
         try{
             DataOutputStream fw = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfile), 65536));
            fw.writeInt(unique);
            fw.writeInt(tagLengthInLong);
            fw.writeInt(tbtb2.getTaxaCount());
            for (int t = 0; t < tbtb2.getTaxaCount(); t++) {
                fw.writeUTF(tbtb2.getTaxaName(t));
            }
            int currTag=0;
            OpenBitSet obs;
            for (int i = 0; i < origOrder.length; i++) {
                obs=tbtb2.getTaxaReadBitsForTag(origOrder[i]);
                while(((i+1)<origOrder.length)&&(compPos.compare(i+1,i)==0)) {
                    OpenBitSet obs2=tbtb2.getTaxaReadBitsForTag(origOrder[i+1]);
                    obs.or(obs2);
                    i++;
                }
                for (int j = 0; j < tagLengthInLong; j++){
                    fw.writeLong(sortTags[j][i]); //Write long[] array of sequence to disk
                }
//                byte len=(tbtb2.getTagLength(origOrder[i])<64)?(byte)(tbtb2.getTagLength(origOrder[i])+4):(byte)64;
//                fw.writeByte(len);
                fw.writeByte(tbtb2.getTagLength(origOrder[i]));
                long[] obsInLong=obs.getBits();
                for (int t = 0; t < obsInLong.length; t++) {
                    fw.writeLong(obsInLong[t]);	//Write presence/absence bitmap to file
                }
            }
            fw.close();
    }  catch (Exception e) {
	System.out.println("Catch in writing output file: " + e);
	e.printStackTrace();
    }

    }
 
 
/**@param inFile Name of a tag map file in the format supplied by Harpreet.
 * @param taxonNameFile Name of a file which lists biological line names present in infile, one per line.
 * @param outFile Name of a file in which to write new tag map. */
    public void convertCBSUFile(File inFile, File taxonNameFile, File outFile){
	System.out.println("Reading tag alignment from:" + inFile.toString());
	String[] inputLine={"Not","Started"};
	int numTaxonNames = 0; //Number of distinct taxon names from taxonNameFile
	int currLine = 0;

	try{
	    BufferedReader taxonNameCounter = new BufferedReader(new FileReader(taxonNameFile), 65536);
	    while(taxonNameCounter.readLine() != null){numTaxonNames++;} //1st pass through file counts lines
	    taxonNameCounter.close();
	}catch(Exception e){
	    System.out.println("Catch in counting taxon names: " + e);
	}
	
	String[] taxonNames = new String[numTaxonNames];
	try{
	    BufferedReader taxonNameReader = new BufferedReader(new FileReader(taxonNameFile), 65536);
	    for(int i=0; i<numTaxonNames; i++){taxonNames[i] = taxonNameReader.readLine().trim();} //Trim each line name and add to array
	    taxonNameReader.close();
	}catch(Exception e){
	    System.out.println("Catch in reading taxon names: " + e);
	}
        System.out.println("Taxa names read.  TaxaNumber:"+numTaxonNames);
    try {
	BufferedReader br = new BufferedReader(new FileReader(inFile), 65536);
	BufferedReader tagCounter = new BufferedReader(new FileReader(inFile), 65536);
	DataOutputStream fw = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outFile), 65536));
	int tagNum = 0;
	int tagLengthInLong = 2;
	ReadBarcodeResult tagProcessingResults = null;
	byte tagLength;
	String tagSequence;
	long[] binaryTagSequence = new long[tagLengthInLong];
        System.out.print("Tag counting....");
	while(tagCounter.readLine()!=null){tagNum++;}  tagCounter.close();  //Loop through tag map file, count #tags (i.e. #lines) and then close
        System.out.println("  TagNumber:"+tagNum);
	try{
	    fw.writeInt(tagNum);
	    fw.writeInt(tagLengthInLong);
	    fw.writeInt(numTaxonNames);
	    for (int t = 0; t < numTaxonNames; t++) {
		fw.writeUTF(taxonNames[t]);
	    }
	    for (int row = 0; row < tagNum; row++) { //Fill byte arrays
                if(row%100000==0) System.out.println("Processed Row:"+row);
		inputLine = br.readLine().split("\t");
                String corrSeq=inputLine[5];

                if(
                        !corrSeq.startsWith("CAGC") &&
                        !corrSeq.startsWith("CTGC") &&
                        !corrSeq.endsWith("CGTC") &&
                        !corrSeq.endsWith("CGAC") 
                ) {
                    continue;
                }

		if(cleanCutSites){
                    if(ParseBarcodeRead.getInitialCutSiteRemnant()==null){ParseBarcodeRead.chooseEnzyme("ApeKI");}
                    tagProcessingResults = ParseBarcodeRead.removeSeqAfterSecondCutSite(corrSeq,(byte)(tagLengthInLong*32));}  //Process tag sequence to find cut sites
		if(tagProcessingResults.processedSequence != null){
		    tagLength = tagProcessingResults.getLength(); //If cut site was found, write length of processed tag...
		    tagSequence = tagProcessingResults.paddedSequence;//And sequence of processed tag...
		}else{
		    tagLength = (byte)corrSeq.length(); //...otherwise write length of unprocessed tag.
		    tagSequence = corrSeq; //...and sequence of unprocessed tag.
		}

		binaryTagSequence = BaseEncoder.getLongArrayFromSeq(tagSequence); //Convert sequence from string to long[]

		for (int j = 0; j < tagLengthInLong; j++){
		    fw.writeLong(binaryTagSequence[j]); //Write long[] array of sequence to disk
		}
		fw.writeByte(tagLength);

		OpenBitSet obs = new OpenBitSet(numTaxonNames);	//New binary object to store presence/absence data
                for (int t = 0; t < numTaxonNames; t++) { //Create a marker presence/absence bitmap
//		    if(inputLine[9].substring(t, t+1).matches("1")){
//			obs.set(t);
//		    }
                    if(inputLine[9].charAt(t)=='1'){
			obs.set(t);
		    }
                }
		
                long[] obsInLong=obs.getBits();
                for (int t = 0; t < obsInLong.length; t++) {
                    fw.writeLong(obsInLong[t]);	//Write presence/absence bitmap to file
		}
	    }
            fw.close();
	    br.close();
    }  catch (Exception e) {
	System.out.println("Catch in writing output file: " + e);
	e.printStackTrace();
    }
} catch (Exception e) {
    System.out.println("Catch in reading tag map file: " + e);
}
    }

    boolean isInteger(String s) {
        try{Integer.parseInt(s);}
        catch(NumberFormatException nf) {
            return false;
        }
        return true;
    }




}
