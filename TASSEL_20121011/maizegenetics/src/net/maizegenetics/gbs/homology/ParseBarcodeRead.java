package net.maizegenetics.gbs.homology;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import net.maizegenetics.genome.BaseEncoder;



/**
 * Takes a key file and then sets up the methods to decode a tag
 *
 * Developer: ed
 * 
 */
public class ParseBarcodeRead {
    private static int chunkSize=BaseEncoder.chunkSize;
    //String baseDir="E:/SolexaAnal/";
    private int maximumMismatchInBarcodeAndOverhang=0;
    private static String[] initialCutSiteRemnant=null;
    private static int readEndCutSiteRemnantLength;
    static String nullS="AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    private static String[] likelyReadEnd=null;
    private static String theEnzyme=null;

    //Original design CWGGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG  -- Common Adapter
    //Redesign  CWGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG  -- Common Adapter

    static int maxBarcodeLength=10;
    private Barcode[] theBarcodes;
    private long[] quickBarcodeList;
    private HashMap<Long,Integer> quickMap;

    public ParseBarcodeRead(String keyFile, String enzyme, String flowcell, String lane) {
        if (enzyme != null){ chooseEnzyme(enzyme); }
        else { chooseEnzyme(getKeyFileEnzyme(keyFile)); }

        int totalBarcodes=setupBarcodeFiles(new File(keyFile),  flowcell,  lane);
        System.out.println("Total barcodes found in lane:"+totalBarcodes);
    }
    
    public ParseBarcodeRead(String enzyme, String flowcell, String lane) {
        if (enzyme != null){ chooseEnzyme(enzyme); }
        else { System.out.println("You need to provide an enzyme for your read 2 file. Exiting."); System.exit(0); }

       // int totalBarcodes=setupBarcodeFiles(new File(keyFile),  flowcell,  lane);
       // System.out.println("Total barcodes found in lane:"+totalBarcodes);
    }

    /**Determines which cut sites to look for, and sets them, based on the enzyme used to generate the GBS library.
     * For two-enzyme GBS  both enzymes MUST be specified and separated by a dash "-". e.g. PstI-MspI, SbfI-MspI
     * The enzyme pair "PstI-EcoT22I" uses the Elshire common adapter while PstI-MspI, PstI-TaqI, and SbfI-MspI use a Y adapter (Poland et al. 2012)
     @param enzyme  The name of the enzyme (case insensitive)*/
        //TODO these should all be private static final globals, then just use this set which one is active.
    public static void chooseEnzyme(String enzyme){
        // Check for case-insensitive (?i) match to a known enzyme
        // The common adapter is: [readEndCutSiteRemnant]AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG  
        if(enzyme.matches("(?i)apek[i1]")){
            theEnzyme = "ApeKI";
            initialCutSiteRemnant=new String[]{"CAGC","CTGC"};
            likelyReadEnd= new String[]{"GCAGC","GCTGC","GCAGAGAT","GCTGAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 4;
        }else if(enzyme.matches("(?i)pst[i1]")){  
            theEnzyme = "PstI";
            initialCutSiteRemnant=new String[]{"TGCAG"};
            likelyReadEnd = new String[]{"CTGCAG","CTGCAAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 5;
        }else if(enzyme.matches("(?i)ecot22[i1]")){
            theEnzyme = "EcoT22I";
            initialCutSiteRemnant= new String[]{"TGCAT"};
            likelyReadEnd = new String[]{"ATGCAT","ATGCAAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 5;
        }else if(enzyme.matches("(?i)pas[i1]")){
            theEnzyme = "PasI";
            initialCutSiteRemnant=new String[]{"CAGGG","CTGGG"};
            likelyReadEnd= new String[]{"CCCAGGG","CCCTGGG","CCCTGAGAT","CCCAGAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 5;
        }else if(enzyme.matches("(?i)hpaii|(?i)hpa2")){
            theEnzyme = "HpaII";
            initialCutSiteRemnant=new String[]{"CGG"};
            likelyReadEnd = new String[]{"CCGG","CCGAGATCGG"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 3;
        }else if(enzyme.matches("(?i)msp[i1]")){  
            theEnzyme = "MspI";  // MspI and HpaII are isoschizomers (same recognition seq and overhang)
            initialCutSiteRemnant=new String[]{"CGG"};
            likelyReadEnd = new String[]{"CCGG","CCGAGATCGG"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 3;
        }else if(enzyme.matches("(?i)pst[i1]-ecot22[i1]")){
            theEnzyme = "PstI-EcoT22I";
            initialCutSiteRemnant= new String[]{"TGCAG","TGCAT"};
            likelyReadEnd = new String[]{"ATGCAT","CTGCAG","CTGCAAGAT","ATGCAAGAT"}; // look for EcoT22I site, PstI site, or common adapter for PstI/EcoT22I
            readEndCutSiteRemnantLength = 5;
        }else if(enzyme.matches("(?i)pst[i1]-msp[i1]")){
            theEnzyme = "PstI-MspI";
            initialCutSiteRemnant=new String[]{"TGCAG"};
            // corrected, change from  CCGAGATC to CCGCTCAGG, as Y adapter was used for MspI  -QS
            likelyReadEnd = new String[]{"CCGG","CTGCAG","CCGCTCAGG"}; // look for MspI site, PstI site, or common adapter for MspI
            readEndCutSiteRemnantLength = 3;
        //NCSU added msp-pst, NEED TO VERIFY COMMON ADAPTER FOR PstI    
        }else if(enzyme.matches("(?i)msp[i1]-pst[i1]")){
            theEnzyme = "MspI-PstI";
            initialCutSiteRemnant=new String[]{"CGG"};
            likelyReadEnd = new String[]{"CTGCAG","CCGG","CTGCAAGAT"}; // look for PstI site, MspI site, or common adapter for PstI
            readEndCutSiteRemnantLength = 3;
        }else if(enzyme.matches("(?i)pst[i1]-taq[i1]")){
            theEnzyme = "PstI-TaqI";
            // corrected, change from  TCGAGATC to TCGCTCAGG, as Y adapter was used for TaqI  -QS
            initialCutSiteRemnant=new String[]{"TGCAG"};
            likelyReadEnd = new String[]{"TCGA","CTGCAG","TCGCTCAGG"}; // look for TaqI site, PstI site, or common adapter for TaqI
            readEndCutSiteRemnantLength = 3;
        }else if(enzyme.matches("(?i)sbf[i1]-msp[i1]")){
            theEnzyme = "SbfI-MspI";
            initialCutSiteRemnant=new String[]{"TGCAGG"};
            // corrected, change from  CCGAGATC to CCGCTCAGG, as Y adapter was used for MspI  -QS
            likelyReadEnd = new String[]{"CCGG","CCTGCAGG","CCGCTCAGG"}; // look for MspI site, SbfI site, or common adapter for MspI
            readEndCutSiteRemnantLength = 3;   
        }else if(enzyme.matches("(?i)asis[i1]-msp[i1]")){
            theEnzyme = "AsiSI-MspI";
            initialCutSiteRemnant=new String[]{"ATCGC"};
            // likelyReadEnd for common adapter is CCGCTCAGG, as the Poland et al.(2012) Y adapter was used for MspI
            likelyReadEnd = new String[]{"CCGG","GCGATCGC","CCGCTCAGG"}; // look for MspI site, AsiSI site, or common adapter for MspI
            readEndCutSiteRemnantLength = 3;   
        }else if(enzyme.matches("(?i)bsshii-msp[i1]|(?i)bssh2-msp[i1]")){
            theEnzyme = "BssHII-MspI";
            initialCutSiteRemnant=new String[]{"CGCGC"};
            // likelyReadEnd for common adapter is CCGCTCAGG, as the Poland et al.(2012) Y adapter was used for MspI
            likelyReadEnd = new String[]{"CCGG","GCGCGC","CCGCTCAGG"}; // look for MspI site, BssHII site, or common adapter for MspI
            readEndCutSiteRemnantLength = 3;   
        }else if(enzyme.matches("(?i)fse[i1]-msp[i1]")){
            theEnzyme = "FseI-MspI";
            initialCutSiteRemnant=new String[]{"CCGGCC"};
            // likelyReadEnd for common adapter is CCGCTCAGG, as the Poland et al.(2012) Y adapter was used for MspI
            likelyReadEnd = new String[]{"CCGG","GGCCGGCC","CCGCTCAGG"}; // look for MspI site, FseI site, or common adapter for MspI
            readEndCutSiteRemnantLength = 3;   
        }else if(enzyme.matches("(?i)sal[i1]-msp[i1]")){
            theEnzyme = "SalI-MspI";
            initialCutSiteRemnant=new String[]{"TCGAC"};
            // likelyReadEnd for common adapter is CCGCTCAGG, as the Poland et al.(2012) Y adapter was used for MspI
            likelyReadEnd = new String[]{"CCGG","GTCGAC","CCGCTCAGG"}; // look for MspI site, SalI site, or common adapter for MspI
            readEndCutSiteRemnantLength = 3;   
        } else if(enzyme.matches("(?i)apo[i1]")){
            theEnzyme = "ApoI";
            initialCutSiteRemnant=new String[]{"AATTC","AATTT"};  // corrected from {"AATTG","AATTC"} by Jeff Glaubitz on 2012/09/12
            likelyReadEnd = new String[]{"AAATTC","AAATTT","GAATTC","GAATTT","AAATTAGAT","GAATTAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 5;
        } else if(enzyme.matches("(?i)BamH[i1l]")){
            theEnzyme = "BamHI";
            initialCutSiteRemnant = new String[]{"GATCC"};
            likelyReadEnd = new String[]{"GGATCC", "GGATCAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
//            likelyReadEnd = new String[]{"GGATCC", "AGATCGGAA", "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG"}; // <-- corrected from this by Jeff Glaubitz on 2012/09/12
            readEndCutSiteRemnantLength = 5;
        } else if(enzyme.matches("(?i)mse[i1]")){
            theEnzyme = "MseI";
            initialCutSiteRemnant=new String[]{"TAA"};
            likelyReadEnd = new String[]{"TTAA","TTAAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)RBSTA")) {
            theEnzyme = "RBSTA";
            initialCutSiteRemnant = new String[]{"TA"};
            likelyReadEnd = new String[]{"TTAA", "GTAC", "CTAG", "TTAAGAT", "GTAAGAT", "CTAAGAT"}; // full cut site (from partial digest or chimera) of MseI, CVIQi, XspI or common adapter start
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)RBSCG")) {
            theEnzyme = "RBSCG";
            initialCutSiteRemnant = new String[]{"CG"};
            likelyReadEnd = new String[]{"CCGC", "TCGA", "GCGC", "CCGG", "ACGT", "CCGAGAT", "TCGAGAT", "GCGAGAT", "ACGAGAT"}; // full cut site (from partial digest or chimera) of AciI, TaqaI, HinpI, HpaII, HpyCH4IV or common adapter start
            readEndCutSiteRemnantLength = 3;
        } else{
            System.out.println("The software didn't recognize your cut site.  "
                    +"Currently, only ApeKI, PstI, EcoT22I, PasI, HpaII, MspI, ApoI, BamHI, MseI, RBSTA, and RBSCG are recognized for single enzyme digests, "
                    +"or PstI-EcoT22I, PstI-MspI, PstI-TaqI, SbfI-MspI, AsiSI-MspI, BssHII-MspI, FseI-MspI, or SalI-MspI for two-enzyme digests.");
            System.out.println("For two-enzyme digest, enzyme names should be separated by a dash, e.g. PstI-MspI ");
        }
        System.out.println("Enzyme: " + theEnzyme);
    }

    private String getKeyFileEnzyme(String keyFileName){
        String result=null;
         try{
            BufferedReader br=new BufferedReader(new FileReader(keyFileName), 65536);

            String temp;
            int currLine=0;
            while (((temp = br.readLine()) != null)) {
                String[] s=temp.split("\\t");  //split by whitespace
                if(currLine >0){
                    String enzymeName;
                    if (s.length < 9) {
                        enzymeName="";
                    } else {
                        enzymeName=s[8];
                    }
                    if (!enzymeName.equals("")){
                        result=enzymeName;
                        break;
                    }
                }
                currLine++;
            }
        } catch(Exception e){
            System.out.println("Couldn't open key file to read Enzyme: "+e);
        }
         return result;
    }
    
/* Reads in an Illumina key file, creates a linear array of {@link Barcode} objects
 * representing the barcodes in the key file, then creates a hash map containing
 * indices from the linear array indexed by sequence.  The names of barcode objects
 * follow the pattern samplename:flowcell:lane:well, since sample names alone are not unique.
 *
 * @param keyFile Illumina key file.
 * @param flowcell Only barcodes from this flowcell will be added to the array.
 * @param lane Only barcodes from this lane will be added to the array.
 * @return Number of barcodes in the array.  */
    private int setupBarcodeFiles(File keyFile, String flowcell, String lane) {
        try{
            BufferedReader br=new BufferedReader(new FileReader(keyFile), 65536);
            ArrayList<Barcode> theBarcodesArrayList=new ArrayList<Barcode>();
            String temp;
            while (((temp = br.readLine()) != null)) {
                String[] s=temp.split("\\t");  //split by whitespace
                Barcode theBC=null;
                if(s[0].equals(flowcell)&&s[1].equals(lane)) {
                    String well = (s[6].length()<2)? (s[5]+'0'+s[6]) : s[5]+s[6]; 
                    if (s.length < 8 || s[7] == null || s[7].equals("")) {  // use the plate and well
                        theBC=new     Barcode(s[2],initialCutSiteRemnant,s[3]+":"+s[0]+":"+s[1]+":"+s[4]+":"+well, flowcell, lane);
                    } else {  // use the "libraryPlateWellID" or whatever is in column H of the key file, IF it is an integer
                        try {
                            int libPrepID = Integer.parseInt(s[7]);
                            theBC=new Barcode(s[2],initialCutSiteRemnant,s[3]+":"+s[0]+":"+s[1]+":"+libPrepID, flowcell, lane);
                        } catch (NumberFormatException nfe) {
                            theBC=new Barcode(s[2],initialCutSiteRemnant,s[3]+":"+s[0]+":"+s[1]+":"+s[4]+":"+well, flowcell, lane);
                        }
                    }
                    theBarcodesArrayList.add(theBC);
                    System.out.println(theBC.barcodeS+" "+theBC.taxaName);
                }
            }
            theBarcodes=new Barcode[theBarcodesArrayList.size()];
            theBarcodesArrayList.toArray(theBarcodes);
            Arrays.sort(theBarcodes);
            int nBL = theBarcodes[0].barOverLong.length;
            quickBarcodeList=new long[theBarcodes.length*nBL];
            quickMap = new HashMap();
            for (int i = 0; i < theBarcodes.length; i++) {
            	for (int j=0; j<nBL; j++){
                    quickBarcodeList[i*nBL+j]=theBarcodes[i].barOverLong[j];
                    quickMap.put(theBarcodes[i].barOverLong[j], i);
            	}
            }
            Arrays.sort(quickBarcodeList);
        } catch(Exception e){
            System.out.println("Error with setupBarcodeFiles: "+e);
        }

        return theBarcodes.length;
    }

    private Barcode findBestBarcode(String queryS, int maxDivergence) {
        long query=BaseEncoder.getLongFromSeq(queryS.substring(0, chunkSize));
        //note because the barcodes are polyA after the sequence, they should always
        //sort ahead of the hit, this is the reason for the -(closestHit+2)
        int closestHit=Arrays.binarySearch(quickBarcodeList, query);

/*      THIS IS THE NEW PIPELINE APPROACH THAT DOES NOT WORK
        if(closestHit>-2) return null; //hit or perfect
        if((query&quickBarcodeList[-(closestHit+2)])!=quickBarcodeList[-(closestHit+2)]) return null;
        int index =quickMap.get(quickBarcodeList[-(closestHit+2)]);
  //      System.out.println(theBarcodes[index].barcodeS);
        return theBarcodes[index];
	//note to see if it is a perfect match you can just bit AND
*/

        //  Below is the old pipeline approach, which works (at least for maxDivergence of 0)
        if (closestHit < -1) {  // should always be true, as the barcode+overhang is padded to 32 bases with polyA
            int index = quickMap.get(quickBarcodeList[-(closestHit+2)]);
            if (theBarcodes[index].compareSequence(query, 1)==0) {
                return theBarcodes[index];
            } else if(maxDivergence==0) {
                return null;  // return null if not a perfect match
            }
        }
        else return null;  // should never go to this line

        int maxLength=0, minDiv=maxDivergence+1, countBest=0;
        Barcode bestBC=null;
        for(Barcode bc: theBarcodes) {
            int div=bc.compareSequence(query, maxDivergence+1);
            if(div<=minDiv) {
                if((div<minDiv)||(bc.barOverLength>maxLength)) {
                    minDiv=div;
                    maxLength=bc.barOverLength;
                    bestBC=bc;
                    countBest=1;
                } else {  //it is a tie, so return that not resolvable
                    bestBC=null;
                    countBest++;
                }
            }
        }
        return bestBC;
    }
    
    /**
     * Does a string search for the any exact match barcode.  If found
     * that barcode is returned, otherwise returns null
     * @param headerBC is the index string from the sequence file that should contain a barcode
     * @return	The match barcode found in the header or null if not match is found.
     */
    private Barcode findHeaderBestBarcode(String headerBC) {

    	int counter = theBarcodes.length;
    	
    	for(int i=0; i<counter;i++){
    		//matches a barcode
    		if(theBarcodes[i].getBarcode().contains(headerBC)){
    			//DEBUGGING ELEMENT
    			System.out.println("Barcode detected in header is: " + theBarcodes[i].getBarcode());
    			// END DEBUGGING
    			return theBarcodes[i];
    		}
    	}
    	//DEBUGGING ELEMENT
		System.out.println("NO BARCODE MATCHED IN HEADER: " + headerBC);
		// END DEBUGGING
       return null;
    }
    
    /**
     * Does a string search for the any exact match barcode.  If found
     * that barcode is returned, otherwise returns null
     * @param queryS The sequence string to search
     * @return	The exact match barcode or null if not match is found.
     */
    private Barcode forceFindBestBarcode(String queryS) {
    	//boolean found=false;
    	// assuming chunkSize is 32, look for barcode in first 16 bases
    	String targetArea = queryS.substring(0, (chunkSize/2));	
    	int counter = theBarcodes.length;
    	
    	for(int i=0; i<counter;i++){
    		//matches a barcode completely
    		if(targetArea.indexOf(theBarcodes[i].getBarcode())!=-1){
    			return theBarcodes[i];
    		}
    	}
       return null;
    }


    /**
     * The barcode libraries used for this study can include two types of extraneous sequence
     * at the end of reads.  The first are chimeras created with the free ends.  These will
     * recreate the restriction site.  The second are short regions (less than 64bp), so that will they
     * will contain a portion of site and the universal adapter.
     * This finds the first of site in likelyReadEnd, keeps the restriction site overhang and then sets everything
     * to polyA afterwards
     * @param seq An unprocessed tag sequence.
     * @param maxLength The maximum number of bp in the processed sequence.
     * @return returnValue A ReadBarcodeResult object containing the unprocessed tag, Cut site position, Processed tag, and Poly-A padded tag.
     */
    public static ReadBarcodeResult removeSeqAfterSecondCutSite(String seq, byte maxLength) {
        //this looks for a second restriction site or the common adapter start, and then turns the remaining sequence to AAAA
        int cutSitePosition = 9999;
        ReadBarcodeResult returnValue = new ReadBarcodeResult(seq);

        //Look for cut sites, starting at a point past the length of the initial cut site remnant that all reads begin with
        String match=null;
        for (String potentialCutSite : likelyReadEnd) {
            int p = seq.indexOf(potentialCutSite, 1);
            if( (p>1) && (p<cutSitePosition) )  {
                cutSitePosition=p;
                match = potentialCutSite;
            }
        }
        if(theEnzyme.equalsIgnoreCase("ApeKI") && cutSitePosition==2 
              && ( match.equalsIgnoreCase("GCAGC") || match.equalsIgnoreCase("GCTGC") ) ) {  // overlapping ApeKI cut site: GCWGCWGC
            seq = seq.substring(3,seq.length());  // trim off the initial GCW from GCWGCWGC
            cutSitePosition = 9999;
            returnValue.unprocessedSequence = seq;
            for (String potentialCutSite : likelyReadEnd) {
                int p = seq.indexOf(potentialCutSite, 1);
                if( (p>1) && (p<cutSitePosition) )  cutSitePosition=p;
            }
        }

        if (cutSitePosition < maxLength) {  // Cut site found
            //Trim tag to sequence up to & including the cut site
            returnValue.length = (byte)(cutSitePosition+readEndCutSiteRemnantLength);
            returnValue.processedSequence = seq.substring(0, cutSitePosition+readEndCutSiteRemnantLength);
        } else {
            if (seq.length() <= 0) {
                //If cut site is missing because there is no sequence
                returnValue.processedSequence = "";
                returnValue.length = 0;
            } else {
                //If cut site is missing because it is beyond the end of the sequence (or not present at all)
                returnValue.length = (byte) Math.min(seq.length(), maxLength);
                returnValue.processedSequence = (seq.substring(0, returnValue.length));
            }
        }

        //Pad sequences shorter than max. length with A
        if (returnValue.length < maxLength) {
            returnValue.paddedSequence = returnValue.processedSequence + nullS;
            returnValue.paddedSequence = returnValue.paddedSequence.substring(0, maxLength);
        } else {
        //Truncate sequences longer than max. length
            returnValue.paddedSequence = returnValue.processedSequence.substring(0, maxLength);
            returnValue.length = maxLength;
        }
        return returnValue;
    }

    /**
     * 
     * @param seqS
     * @param qualS
     * @param fastq
     * @param minQual
     * @return If barcode and cut site was found returns the result and processed sequence, if the barcode and cut site were not found return null
     */
    public ReadBarcodeResult parseReadIntoTagAndTaxa(String seqS, String qualS, boolean fastq, int minQual) {
    	 	long[] read=new long[2];
        if((minQual>0)&&(qualS!=null)) {
            int firstBadBase=BaseEncoder.getFirstLowQualityPos(qualS, minQual);
            if(firstBadBase<(maxBarcodeLength+2*chunkSize)) return null;
        }
        int miss = -1;
        if (fastq) { miss=seqS.indexOf('N'); } else { miss=seqS.indexOf('.'); }
        if((miss!=-1)&&(miss<(maxBarcodeLength+2*chunkSize))) return null;  //bad sequence so skip
        Barcode bestBarcode=findBestBarcode(seqS,maximumMismatchInBarcodeAndOverhang);
        if(bestBarcode==null) return null;  //overhang missing so skip
        String genomicSeq=seqS.substring(bestBarcode.barLength, seqS.length());
        ReadBarcodeResult tagProcessingResults = removeSeqAfterSecondCutSite(genomicSeq, (byte)(2*chunkSize));
        String hap=tagProcessingResults.paddedSequence;  //this is slow 20% of total time.   Tag, cut site processed, padded with poly-A

        read=BaseEncoder.getLongArrayFromSeq(hap);
        int pos=tagProcessingResults.length;
        //TODO this instantiation should also include the orginal unprocessedSequence, processedSequence, and paddedSequence - the the object encode it 

        ReadBarcodeResult rbr=new ReadBarcodeResult(read, (byte)pos, bestBarcode.getTaxaName());

        return rbr; 

    }
    
    /**
     * Scan the submitted sequence and quality score to make sure sequence passes
     * criteria of having no missing bases (N) in the target tag region. OVERLOAD: A
     * minimum cut-off value is set here to reject any sequence, regardless of quality,
     * if it is not the desired length 
     * @param seqS	The raw sequence from the data file
     * @param qualS	The quality score from the sequence run
     * @param fastq	A boolean value that controls if a missing value (N) present in the 
     * sequence should be cause to immediately reject the sequence for consideration.
     * @param minQual Logic controller alsong with qualS that is used to check the quality score
     * @param lengthToKeep The exact acceptable sequence length to keep for further processing
     * @return If barcode and cut site was found, and the length was acceptable, returns the result and processed sequence, 
     * 		if the barcode and cut site were not found return null
     */
    public ReadBarcodeResult parseReadIntoTagAndTaxa(String seqS, String qualS, boolean fastq, int minQual, int lengthToKeep) {
    	long[] read=new long[2];
        if((minQual>0)&&(qualS!=null)) {
            int firstBadBase=BaseEncoder.getFirstLowQualityPos(qualS, minQual);
            if(firstBadBase<(maxBarcodeLength+2*chunkSize)) return null;
        }
        int miss = -1;
        if (fastq) { miss=seqS.indexOf('N'); } else { miss=seqS.indexOf('.'); }
        if((miss!=-1)&&(miss<(maxBarcodeLength+2*chunkSize))) return null;  //bad sequence so skip
        Barcode bestBarcode=findBestBarcode(seqS,maximumMismatchInBarcodeAndOverhang);
        if(bestBarcode==null) return null;  //overhang missing so skip
        String genomicSeq=seqS.substring(bestBarcode.barLength, seqS.length());
        ReadBarcodeResult tagProcessingResults = removeSeqAfterSecondCutSite(genomicSeq, (byte)(2*chunkSize));
        if(tagProcessingResults.length != lengthToKeep) return null; // sequence not desired length
        String hap=tagProcessingResults.paddedSequence;  //this is slow 20% of total time.   Tag, cut site processed, padded with poly-A

        read=BaseEncoder.getLongArrayFromSeq(hap);
        int pos=tagProcessingResults.length;
        //TODO this instantiation should also include the orginal unprocessedSequence, processedSequence, and paddedSequence - the the object encode it 

        ReadBarcodeResult rbr=new ReadBarcodeResult(read, (byte)pos, bestBarcode.getTaxaName());

        return rbr; 
    }
    
    /**
     * Modified version of overloaded parseReadIntoTagAndTaxa(String, String, boolean, int, int).  The quality 
     * and bad base detectiong (presense of Ns) have been bypassed.  This calls upon the forceFindBestBarcode looking
     * for an exact barcode match since there are possible barcode variants in the keyfile generated to catch
     * missing or subsituted bases.  Adjustments are made to where the tag starts based on if the enzyme cut side is
     * immediately after the barcode or if there is 1 base away.  Since missing bases (N) may be present and can't be encoded
     * properly by other class methods called upon, the tag sequence assignment to the String that feeds the ReadBarcodeResult
     * methods was changed to handle Strings directly instead of encoding / decoding bits.  The return value is still of type
     * ReadBarcodeResult, but of an overloaded object that has a String sequence elemt.
     * @param seqS	The raw sequence from the data file
     * @param qualS The quality score from the sequence run - BYPASSED
     * @param fastq A boolean value that controls if a missing value (N) present in the 
     * sequence should be cause to immediately reject the sequence for consideration.
     * @param minQual Logic controller alsong with qualS that is used to check the quality score
     * @param lengthToKeep The exact acceptable sequence length to keep for further processing
     * @return If all checks and bypasses are sucessful, returns a ReadBarcodeResult object that contains
     * a direct String variable that holds the tag sequence of interest.
     */
    public ReadBarcodeResult parseReadIntoTagAndTaxaFlexible(String seqS, String qualS, boolean fastq, int minQual, int lengthToKeep) {
    	long[] read=new long[2];
        if((minQual>0)&&(qualS!=null)) {
            int firstBadBase=BaseEncoder.getFirstLowQualityPos(qualS, minQual);
            if(firstBadBase<(maxBarcodeLength+(2*chunkSize))) return null;  
        }
        
        int miss = -1;
        if (fastq) { miss=seqS.indexOf('N'); } else { miss=seqS.indexOf('.'); }
        
       // if(miss!=-1)System.out.println("parseReadIntoTagAndTaxaFlexible "+miss);  //DEBUGGING
        
        if((miss!=-1)&&(miss<(maxBarcodeLength+2*chunkSize))) return null;  //bad sequence
        
        Barcode bestBarcode=forceFindBestBarcode(seqS);
        if(bestBarcode==null) return null;  //barcode missing
        
        String genomicSeq=null;
        
        int cutSiteClose= -1;
        cutSiteClose=checkForCutSite(seqS, bestBarcode);

        if(cutSiteClose == -1){
        	// cut site not found or within designated parameter
        	return null; 
        }else if(cutSiteClose == 0){
        	// take everything after the barcode 
        	genomicSeq=seqS.substring((seqS.indexOf(bestBarcode.getBarcode())+bestBarcode.barLength), seqS.length());
        }else if(cutSiteClose == 1){
        	// take everything after the barcode +1 nucleotide, should be start of cut site remanant
        	genomicSeq=seqS.substring((seqS.indexOf(bestBarcode.getBarcode())+bestBarcode.barLength+1), seqS.length());
        }
        
        //this is slow 20% of total time.   Tag, cut site processed, padded with poly-A
        ReadBarcodeResult tagProcessingResults = removeSeqAfterSecondCutSite(genomicSeq, (byte)(2*chunkSize));

        if(tagProcessingResults.length != lengthToKeep){
    	   return null; // sequence not desired length
       }

        String hap=tagProcessingResults.paddedSequence;  
        
        read=BaseEncoder.getLongArrayFromSeq(hap);
        int pos=tagProcessingResults.length;
        //TODO this instantiation should also include the orginal unprocessedSequence, processedSequence, and paddedSequence - the the object encode it 

        // the key below is the string hap is set directly as an unadulterated string that can be accessed by the
        // calling method directly without conversion to and from bit encoding
        ReadBarcodeResult rbr=new ReadBarcodeResult(read, hap, (byte)pos, bestBarcode.getTaxaName());

        return rbr; 
    }
    
    
    public ReadBarcodeResult parseReadIntoTagAndTaxaNoBarcode(String seqS, String qualS, boolean fastq, int minQual, int lengthToKeep) {
    	long[] read=new long[2];
        if((minQual>0)&&(qualS!=null)) {
            int firstBadBase=BaseEncoder.getFirstLowQualityPos(qualS, minQual);
            if(firstBadBase<(maxBarcodeLength+(2*chunkSize))) return null;  
        }
        
        int miss = -1;
        if (fastq) { miss=seqS.indexOf('N'); } else { miss=seqS.indexOf('.'); }
        
       // if(miss!=-1)System.out.println("parseReadIntoTagAndTaxaNoBarcode "+miss);  // DEBUGGING
        
        if((miss!=-1)&&(miss<(maxBarcodeLength+2*chunkSize))) return null;  //bad sequence
        
        String genomicSeq=null;
        
        int cutSiteClose= -1;
        cutSiteClose=checkForCutSite(seqS);  // look only for cut site near begining of sequence

        if(cutSiteClose == -1){
        	// cut site not found or within designated parameter
        	return null; 
        }else{
        	genomicSeq=seqS.substring(cutSiteClose, seqS.length());
        }
        
        //this is slow 20% of total time.   Tag, cut site processed, padded with poly-A
        ReadBarcodeResult tagProcessingResults = removeSeqAfterSecondCutSite(genomicSeq, (byte)(2*chunkSize));

        if(tagProcessingResults.length != lengthToKeep){
    	   return null; // sequence not desired length
       }

        String hap=tagProcessingResults.paddedSequence;  
        
        read=BaseEncoder.getLongArrayFromSeq(hap);
        int pos=tagProcessingResults.length;
        String taxonID=setCutSiteAsID();

        // the key below is the string hap is set directly as an unadulterated string that can be accessed by the
        // calling method directly without conversion to and from bit encoding, the variable taxonID is used 
        // to placehold an identifier since the data does not have a separate set of barcode identifiers
        ReadBarcodeResult rbr=new ReadBarcodeResult(read, hap, (byte)pos, taxonID);

        return rbr; 
    }
    
    
    public ReadBarcodeResult parseReadIntoTagAndTaxaHeader(String header, String seqS, String qualS, boolean fastq, int minQual, int lengthToKeep) {
    	long[] read=new long[2];
        if((minQual>0)&&(qualS!=null)) {
            int firstBadBase=BaseEncoder.getFirstLowQualityPos(qualS, minQual);
            if(firstBadBase<(maxBarcodeLength+(2*chunkSize))) return null;  
        }
        
        int miss = -1;
        if (fastq) { miss=seqS.indexOf('N'); } else { miss=seqS.indexOf('.'); }
       
        if((miss!=-1)&&(miss<(maxBarcodeLength+2*chunkSize))) return null;  //bad sequence
        
        Barcode bestBarcode=findHeaderBestBarcode(header);
   
        if(bestBarcode==null) return null;  //barcode missing
      
        String genomicSeq=null;
        genomicSeq=seqS.substring(0, seqS.length());		//DEBUGGING
        
        /*  TODO the necessity of cutSiteClose and related methods needs to be
         * reviewed.  It may not be necessary for indexed samples, especially
         * if the stringency is set to high, and if set to low presumably post-processing
         * will draw out miss-called or maladjusted cut site remnants
        int cutSiteClose= -1;
        cutSiteClose=checkForCutSite(seqS);

        if(cutSiteClose == -1){
        	// cut site not found or within designated parameter
        	return null; 
        }else{
        */
        	genomicSeq=seqS.substring(0, seqS.length());
      //  }
        
        //this is slow 20% of total time.   Tag, cut site processed, padded with poly-A
        ReadBarcodeResult tagProcessingResults = removeSeqAfterSecondCutSite(genomicSeq, (byte)(2*chunkSize));

        if(tagProcessingResults.length != lengthToKeep){
    	   return null; // sequence not desired length
       }

        String hap=tagProcessingResults.paddedSequence;  
        
        read=BaseEncoder.getLongArrayFromSeq(hap);
        int pos=tagProcessingResults.length;
        String taxonID=bestBarcode.getBarcode();

        // the key below is the string hap is set directly as an unadulterated string that can be accessed by the
        // calling method directly without conversion to and from bit encoding, the variable taxonID is used 
        // to placehold identification
        ReadBarcodeResult rbr=new ReadBarcodeResult(read, hap, (byte)pos, taxonID);

        return rbr; 
    }
    
    public ReadBarcodeResult parseReadIntoTagAndTaxaHeadPoorSeq(String header, String seqS, String qualS, boolean fastq, int minQual, int lengthToKeep) {
    	long[] read=new long[2];
        if((minQual>0)&&(qualS!=null)) {
            int firstBadBase=BaseEncoder.getFirstLowQualityPos(qualS, minQual);
            if(firstBadBase<(maxBarcodeLength+(2*chunkSize))) return null;  
        }
        
        int miss = -1;
        if (fastq) { miss=seqS.indexOf('N'); } else { miss=seqS.indexOf('.'); }
       
        if((miss!=-1)&&(miss<(maxBarcodeLength+2*chunkSize))) return null;  //bad sequence
        
       // Barcode bestBarcode=findHeaderBestBarcode(seqS, header);   
        //if(bestBarcode==null) return null;  //barcode missing
      
       // String genomicSeq=null;
      //  genomicSeq=seqS.substring(0, seqS.length());		//DEBUGGING
        
        //this is slow 20% of total time.   Tag, cut site processed, padded with poly-A
        ReadBarcodeResult tagProcessingResults = removeSeqAfterSecondCutSite(seqS, (byte)(2*chunkSize));

        if(tagProcessingResults.length != lengthToKeep){return null;} // sequence not desired length

        String hap=tagProcessingResults.paddedSequence;  
        
        read=BaseEncoder.getLongArrayFromSeq(hap);
        int pos=tagProcessingResults.length;
        String taxonID=header;

        // the key below is the string hap is set directly as an unadulterated string that can be accessed by the
        // calling method directly without conversion to and from bit encoding, the variable taxonID is used 
        // to placehold an identifier since the data does not have a separate set of barcode identifiers
        ReadBarcodeResult rbr=new ReadBarcodeResult(read, hap, (byte)pos, taxonID);

        return rbr; 
    }
    
    private int checkForCutSite(String seq, Barcode bar){
    	int response=-1;
    	int areaToScan = initialCutSiteRemnant[0].length()+1;
    	String seqBegin = seq.substring(0,(chunkSize/2)+areaToScan);
    	int indexOfBar = seqBegin.indexOf(bar.getBarcode());
    	String targetArea = seqBegin.substring((indexOfBar+bar.barLength),(indexOfBar+areaToScan+bar.barLength));
    	int present = -1;
    	
    	
    	present = targetArea.indexOf(initialCutSiteRemnant[0]);
    	
    	if (present!=-1 && present==0){
    		response=0;
    	}else if(present!=-1 && present==1){
    		response=1;
    	}

    	return response;
    }
    
    private int checkForCutSite(String seq){
    	int response=-1;
    	int areaToScan = initialCutSiteRemnant[0].length()+1;
    	String seqBegin = seq.substring(0,(chunkSize/2)+areaToScan);
    	int present = -1;
    	
    	present = seqBegin.indexOf(initialCutSiteRemnant[0]);

    	if (present!=-1){
    		response=present;
    	}
    	
    	return response;
    }
    
    
    private String setCutSiteAsID(){
    	String id=null;
    	for(int i=0;i<initialCutSiteRemnant.length;i++){
    		if(id==null){
    			id=initialCutSiteRemnant[i];
    		}else{
    		id=id+initialCutSiteRemnant[i];
    		}
    	}
    	return id;
    }
    
    public int getBarCodeCount() {
        return theBarcodes.length;
    }

    public Barcode getTheBarcodes(int index) {
        return theBarcodes[index];
    }

    public String[] getTaxaNames(){
        String[] result = new String[getBarCodeCount()];
            for (int i = 0; i < result.length; i++) {
                result[i]=getTheBarcodes(i).getTaxaName();
            }
            return result;
    }
    
    static public String[] getInitialCutSiteRemnant() {
        return initialCutSiteRemnant;
    }

}



