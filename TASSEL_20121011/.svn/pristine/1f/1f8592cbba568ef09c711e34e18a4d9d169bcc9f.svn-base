/*
 * ImportUtils
 */
package net.maizegenetics.pal.alignment;

import net.maizegenetics.pal.datatype.DataType;
import net.maizegenetics.pal.datatype.IUPACNucleotides;
import net.maizegenetics.pal.datatype.TextDataType;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.util.ExceptionUtils;
import net.maizegenetics.util.Utils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.TreeSet;
import java.util.regex.Pattern;


/**
 * The class imports PAL alignment datatypes from
 * various file formats.
 *
 * @author terry
 */
public class ImportUtils {
    
    public static final int NUM_HAPMAP_NON_TAXA_HEADERS = 11;

    private ImportUtils() {
        // Utility Class - do not instantiate.
    }

    /**
     * Creates an alignment from a Hapmap file
     * @param filename
     * @return alignment
     */
    public static Alignment readFromHapmap(String filename, String chrom) {
        try {
            ArrayList<byte[]>[] theBLOBList = new ArrayList[GdpdmBLOBUtils.totalBLOBtypes];
            for (int i = 0; i < theBLOBList.length; i++) {
                theBLOBList[i] = new ArrayList<byte[]>();
            }
            byte[][] alleleBLOB = null;
            //array for BLOBs holding coverage plus header [taxa][SNPs]
            byte[] positionBLOB = null;
            byte[] SNPidBLOB = null;
//            String[][] chromsAvailable = GdpdmBLOBUtils.getFileInfo(Utils.addSuffixIfNeeded(filename, ".txt"), 2, 0);
            String[][] chromsAvailable = GdpdmBLOBUtils.getFileInfo(filename, 2, 0);
            int[] chromInfo = GdpdmBLOBUtils.getChromInFileInfo(chromsAvailable, chrom);
//            theBLOBList = AllelePositionBLOBUtils.createAllelePositionBLOBsFromHapmap(Utils.addSuffixIfNeeded(filename, ".txt"), theBLOBList, chrom, chromInfo);
            theBLOBList = AllelePositionBLOBUtils.createAllelePositionBLOBsFromHapmap(filename, theBLOBList, chrom, chromInfo);
            if (theBLOBList[0].size() > 1) {
                alleleBLOB = new byte[theBLOBList[0].size()][];
                //array for BLOBs holding SNPs plus header [taxa][SNPs]
                for (int i = 0; i < alleleBLOB.length; i++) {
                    alleleBLOB[i] = theBLOBList[0].get(i);
                }
                positionBLOB = (byte[])theBLOBList[1].get(0);
                SNPidBLOB = theBLOBList[4].get(0);
            }
            return new Pack1Alignment(alleleBLOB, positionBLOB, SNPidBLOB);
        }
        catch (Exception e) {
            throw new IllegalArgumentException("Problem creating Alignment: " + ExceptionUtils.getExceptionCauses(e));
        }
    }

    /**
     * Creates an alignment from all chromosomes in Hapmap file
     * @param filename
     * @return Alignment one combine alignment
     */
    public static Alignment readFromHapmap(String filename) {
        try {
            long currentTime = System.currentTimeMillis();
//            String[][] chromsAvailable = GdpdmBLOBUtils.getFileInfo(Utils.addSuffixIfNeeded(filename, ".txt"), 2, 0);
            String[][] chromsAvailable = GdpdmBLOBUtils.getFileInfo(filename, 2, 0);
            //System.out.println(chromsAvailable[1][0]);
            String[] chroms = new String[chromsAvailable.length - 1];
//            int[] chromInfo;
            Alignment[] align = new Alignment[chromsAvailable.length - 1];
            for (int i = 1; i < chromsAvailable.length; i++) {
                chroms[i - 1] = chromsAvailable[i][0];
            }
            long prevTime = currentTime;
        currentTime = System.currentTimeMillis();
        System.out.println("Time to count lines: " + ((currentTime - prevTime) / 1000));
            ArrayList<byte[]>[][] theBLOBList = new ArrayList[chroms.length][GdpdmBLOBUtils.totalBLOBtypes];
            for (int i = 0; i < theBLOBList.length; i++) {
                for (int j = 0; j < theBLOBList[i].length; j++) {
                    theBLOBList[i][j] = new ArrayList<byte[]>();
                }
            }
//            chromInfo = GdpdmBLOBUtils.getChromInFileInfo(chromsAvailable, chrom);
//            theBLOBList = AllelePositionBLOBUtils.createAllelePositionBLOBsFromHapmap(Utils.addSuffixIfNeeded(filename, ".txt"), theBLOBList, chroms, chromsAvailable);
            theBLOBList = AllelePositionBLOBUtils.createAllelePositionBLOBsFromHapmap(filename, theBLOBList, chroms, chromsAvailable);
            prevTime = currentTime;
            currentTime = System.currentTimeMillis();
            System.out.println("Time to read file: " + ((currentTime - prevTime) / 1000));
            for (int i = 0; i < theBLOBList.length; i++) {
                byte[][] alleleBLOB = null;
                //array for BLOBs holding coverage plus header [taxa][SNPs]
                byte[] positionBLOB = null;
                byte[] SNPidBLOB = null;
                if (theBLOBList[i][0].size() > 1) {
                    alleleBLOB = new byte[theBLOBList[i][0].size()][];
                    //array for BLOBs holding SNPs plus header [taxa][SNPs]
                    for (int j = 0; j < alleleBLOB.length; j++) {
                        alleleBLOB[j] = theBLOBList[i][0].get(j);
                    }
                    positionBLOB = (byte[])theBLOBList[i][1].get(0);
                    SNPidBLOB = theBLOBList[i][4].get(0);
                }
                align[i] = new Pack1Alignment(alleleBLOB, positionBLOB, SNPidBLOB);
                prevTime = currentTime;
            currentTime = System.currentTimeMillis();
            System.out.println("Time to create Alignment: " + ((currentTime - prevTime) / 1000));
            }
            //return align;
            if(align.length==1) return align[0];
            Alignment test= CombineAlignment.getInstance(align);
            return test;
        }
        catch (Exception e) {
            e.printStackTrace();
            throw new IllegalArgumentException("Problem creating Alignment: " + ExceptionUtils.getExceptionCauses(e));
        }
    }

    public static Pack1Alignment readFromLZMA(String filenameRoot) {
        throw new UnsupportedOperationException("Unsupported Operation");
    }

    /**
     * Reads data BLOBs from a given zip file
     * @param filenameRoot
     * @return Alignment
     */
    public static Alignment readFromZip(String filenameRoot) {
        try {
            ArrayList<byte[]>[] theBLOBList = new ArrayList[GdpdmBLOBUtils.totalBLOBtypes];
            for (int i = 0; i < theBLOBList.length; i++) {
                theBLOBList[i] = new ArrayList<byte[]>();
            }
            byte[][] alleleBLOB = null;
            //array for BLOBs holding coverage plus header [taxa][SNPs]
            byte[] positionBLOB = null;
            byte[] SNPidBLOB = null;
            theBLOBList = GdpdmBLOBUtils.readBLOBFromZip(Utils.addSuffixIfNeeded(filenameRoot, ".zip"), theBLOBList);
            if (theBLOBList[0].size() > 0) {
                alleleBLOB = new byte[theBLOBList[0].size()][];
                //array for BLOBs holding SNPs plus header [taxa][SNPs]
                for (int i = 0; i < alleleBLOB.length; i++) {
                    alleleBLOB[i] = theBLOBList[0].get(i);
                }
                positionBLOB = (byte[])theBLOBList[1].get(0);
                if (theBLOBList[4].size() == 0) {
                    SNPidBLOB = new byte[1030];
                    SNPidBLOB[GdpdmBLOBUtils.blobTypeField[0]] = 0x35;
                }
                else {
                    SNPidBLOB = theBLOBList[4].get(0);
                }
            }
//            System.out.println(alleleBLOB.length);
//            System.out.println(SNPidBLOB.length);
            return new Pack1Alignment(alleleBLOB, positionBLOB, SNPidBLOB);
        }
        catch (Exception e) {
            throw new IllegalArgumentException("Problem creating Alignment: " + ExceptionUtils.getExceptionCauses(e));
        }
    }

    /**
     * Creates an alignment from BLOBs in a gzip compressed file
     * @param fileName name of gzip file
     * @return Alignment
     */

public static Alignment readFromGZIP(String fileName) {
        try {
            ArrayList<byte[]>[] theBLOBList = new ArrayList[GdpdmBLOBUtils.totalBLOBtypes];
            for (int i = 0; i < theBLOBList.length; i++) {
                theBLOBList[i] = new ArrayList<byte[]>();
            }
            byte[][] alleleBLOB = null;
            byte[] positionBLOB = null;
            byte[] SNPidBLOB = null;
            theBLOBList = GdpdmBLOBUtils.readBLOBfromGZIP(Utils.addSuffixIfNeeded(fileName, ".gz"), theBLOBList);
            if (theBLOBList[0].size() > 0) {
                alleleBLOB = new byte[theBLOBList[GdpdmBLOBUtils.alleleBLOBtype - 49].size()][];
                for (int i = 0; i < alleleBLOB.length; i++) {
                    alleleBLOB[i] = theBLOBList[GdpdmBLOBUtils.alleleBLOBtype - 49].get(i);
                }
                positionBLOB = (byte[])theBLOBList[GdpdmBLOBUtils.allelePositionBLOBtype - 49].get(0);
                if (theBLOBList[GdpdmBLOBUtils.SNPIdBLOBtype - 49].size() == 0) {
                    SNPidBLOB = new byte[1030];
                    SNPidBLOB[GdpdmBLOBUtils.blobTypeField[0]] = 0x35;
                }
                else {
                    SNPidBLOB = theBLOBList[GdpdmBLOBUtils.SNPIdBLOBtype - 49].get(0);
                }
            }
            return new Pack1Alignment(alleleBLOB, positionBLOB, SNPidBLOB);
        }
        catch (Exception e) {
            throw new IllegalArgumentException("Problem creating Alignment: " + ExceptionUtils.getExceptionCauses(e));
        }
    }

    /**
     * Creates an alignment from a set of Plink files and desired chromosome
     * @param PEDfileName
     * @param MAPfileName
     * @param chrom
     * @return an alignment for desired chromosome
     */
    public static Alignment readFromPLINK(String PEDfileName, String MAPfileName, String chrom) {
        try {
            ArrayList<byte[]>[] theBLOBList = new ArrayList[GdpdmBLOBUtils.totalBLOBtypes];
            for (int i = 0; i < theBLOBList.length; i++) {
                theBLOBList[i] = new ArrayList<byte[]>();
            }
            byte[][] alleleBLOB = null;
            //array for BLOBs holding coverage plus header [taxa][SNPs]
            byte[] positionBLOB = null;
            byte[] SNPidBLOB = null;
            String[][] chromsAvailable = GdpdmBLOBUtils.getFileInfo(Utils.addSuffixIfNeeded(MAPfileName, ".map"), 0, 1);
            int[] chromInfo = GdpdmBLOBUtils.getChromInFileInfo(chromsAvailable, chrom);
            int numTaxa = GdpdmBLOBUtils.countLinesInFile(Utils.addSuffixIfNeeded(PEDfileName, ".ped"));
            theBLOBList = AllelePositionBLOBUtils.createAllelePositionBLOBsFromPlink(Utils.addSuffixIfNeeded(PEDfileName, ".ped"), Utils.addSuffixIfNeeded(MAPfileName, ".map"), theBLOBList, chrom, chromInfo, numTaxa);
            if (theBLOBList[0].size() > 1) {
                alleleBLOB = new byte[theBLOBList[0].size()][];
                //array for BLOBs holding SNPs plus header [taxa][SNPs]
                for (int i = 0; i < alleleBLOB.length; i++) {
                    alleleBLOB[i] = theBLOBList[0].get(i);
                }
                positionBLOB = (byte[])theBLOBList[1].get(0);
                SNPidBLOB = theBLOBList[4].get(0);
            }
            return new Pack1Alignment(alleleBLOB, positionBLOB, SNPidBLOB);
        }
        catch (Exception e) {
            throw new IllegalArgumentException("Problem creating Alignment: " + ExceptionUtils.getExceptionCauses(e));
        }
    }

    /**
     * Creates an alignment for each chromosome in a .map/.ped file pair
     * @param PEDfileName
     * @param MAPfileName
     * @return an alignment for each chromosome
     */
    public static Alignment readFromPLINK(String PEDfileName, String MAPfileName) {
        try {
            String[][] chromsAvailable = GdpdmBLOBUtils.getFileInfo(Utils.addSuffixIfNeeded(MAPfileName, ".map"), 0, 1);
//            int[] chromInfo;
            int numTaxa = GdpdmBLOBUtils.countLinesInFile(Utils.addSuffixIfNeeded(PEDfileName, ".ped"));
            String[] chroms = new String[chromsAvailable.length];
            Alignment[] align = new Alignment[chromsAvailable.length];
            for (int i = 0; i < chromsAvailable.length; i++) {
                chroms[i] = chromsAvailable[i][0];
            }
//            chrom = chromsAvailable[i][0];
//            chromInfo = GdpdmBLOBUtils.getChromInFileInfo(chromsAvailable, chrom);
            ArrayList<byte[]>[][] theBLOBList = new ArrayList[chroms.length][GdpdmBLOBUtils.totalBLOBtypes];
            for (int i = 0; i < theBLOBList.length; i++) {
                for (int j = 0; j < theBLOBList[i].length; j++) {
                    theBLOBList[i][j] = new ArrayList<byte[]>();
                }
            }
            byte[][] alleleBLOB = null;
            //array for BLOBs holding coverage plus header [taxa][SNPs]
            byte[] positionBLOB = null;
            byte[] SNPidBLOB = null;
            theBLOBList = AllelePositionBLOBUtils.createAllelePositionBLOBsFromPlink(Utils.addSuffixIfNeeded(PEDfileName, ".ped"), Utils.addSuffixIfNeeded(MAPfileName, ".map"), theBLOBList, chroms, chromsAvailable, numTaxa);
            for (int i = 0; i < theBLOBList.length; i++) {
                if (theBLOBList[i][0].size() > 1) {
                    alleleBLOB = new byte[theBLOBList[i][0].size()][];
                    //array for BLOBs holding SNPs plus header [taxa][SNPs]
                    for (int j = 0; j < alleleBLOB.length; j++) {
                        alleleBLOB[j] = theBLOBList[i][0].get(j);
                    }
                    positionBLOB = (byte[])theBLOBList[i][1].get(0);
                    SNPidBLOB = theBLOBList[i][4].get(0);
                }
                align[i] = new Pack1Alignment(alleleBLOB, positionBLOB, SNPidBLOB);
            }
            //return align;
            return CombineAlignment.getInstance(align);
        }
        catch (Exception e) {
            e.printStackTrace();
            throw new IllegalArgumentException("Problem creating Alignment: " + ExceptionUtils.getExceptionCauses(e));
        }
    }

    /**
     * Creates an alignment from a set of Flapjack files and desired chromosome
     * @param genotypeFileName
     * @param MAPfileName
     * @param chrom
     * @return an alignment for desired chromosome
     */
    public static Alignment readFromFlapjack(String genotypeFileName, String MAPfileName, String chrom) {
        try {
            ArrayList<byte[]>[] theBLOBList = new ArrayList[GdpdmBLOBUtils.totalBLOBtypes];
            for (int i = 0; i < theBLOBList.length; i++) {
                theBLOBList[i] = new ArrayList<byte[]>();
            }
            byte[][] alleleBLOB = null;
            //array for BLOBs holding coverage plus header [taxa][SNPs]
            byte[] positionBLOB = null;
            byte[] SNPidBLOB = null;
            String[][] chromsAvailable = GdpdmBLOBUtils.getFileInfo(Utils.addSuffixIfNeeded(MAPfileName, ".map"), 1, 0);
            int[] chromInfo = GdpdmBLOBUtils.getChromInFileInfo(chromsAvailable, chrom);
            int numTaxa = GdpdmBLOBUtils.countLinesInFile(Utils.addSuffixIfNeeded(genotypeFileName, ".geno")) - 1;
            theBLOBList = AllelePositionBLOBUtils.createAllelePositionBLOBsFromFlapjack(Utils.addSuffixIfNeeded(genotypeFileName, ".geno"), Utils.addSuffixIfNeeded(MAPfileName, ".map"), theBLOBList, chrom, chromInfo, numTaxa);
            if (theBLOBList == null) {
                theBLOBList = AllelePositionBLOBUtils.createAllelePositionBLOBsFromFlapjack(Utils.addSuffixIfNeeded(genotypeFileName, ".geno"), Utils.addSuffixIfNeeded(MAPfileName, ".map"), theBLOBList, chrom, chromInfo, numTaxa);
            }
            if (theBLOBList[0].size() > 1) {
                alleleBLOB = new byte[theBLOBList[0].size()][];
                //array for BLOBs holding SNPs plus header [taxa][SNPs]
                for (int i = 0; i < alleleBLOB.length; i++) {
                    alleleBLOB[i] = theBLOBList[0].get(i);
                }
                positionBLOB = (byte[])theBLOBList[1].get(0);
                SNPidBLOB = theBLOBList[4].get(0);
            }
            return new Pack1Alignment(alleleBLOB, positionBLOB, SNPidBLOB);
        }
        catch (Exception e) {
            throw new IllegalArgumentException("Problem creating Alignment: " + ExceptionUtils.getExceptionCauses(e));
        }
    }

    /**
     * Creates a CombineAlignment from Flapjack files
     * Does NOT assume the file is sorted by position and/or chromosome
     * @param genotypeFileName
     * @param MAPfileName
     * @return CombineAlignment
     */
    public static Alignment readFromFlapjack(String genotypeFileName, String MAPfileName) {
        try {
            Object[][] fileInfo = GdpdmBLOBUtils.getFlapjackInfo(MAPfileName);

            // initializing theBLOBList
            ArrayList<byte[]>[][] theBLOBList = new ArrayList[fileInfo.length][GdpdmBLOBUtils.totalBLOBtypes];
            for (int i =0; i < theBLOBList.length; i++) {
                for (int j = 0; j < theBLOBList[i].length; j++) {
                    theBLOBList[i][j] = new ArrayList<byte[]>();
                }
            }

            String[] chromNames = new String[fileInfo.length];
            Hashtable[] tableArray = new Hashtable[fileInfo.length];
            int[] idLengths = new int[fileInfo.length];
            int[] numSites = new int[fileInfo.length];

            int numTaxa = GdpdmBLOBUtils.countLinesInFile(Utils.addSuffixIfNeeded(genotypeFileName, ".geno")) - 1;

            // defining fields describing flapjack file
            for (int i =0; i < fileInfo.length; i++) {
                chromNames[i] = (String)fileInfo[i][0];
                tableArray[i] = (Hashtable)fileInfo[i][1];
                idLengths[i] = ((Integer)fileInfo[i][2]).intValue();
                numSites[i] = tableArray[i].size();
            }

            // populate theBLOBList
            theBLOBList = AllelePositionBLOBUtils.createAllelePositionBLOBsFromFlapjack(Utils.addSuffixIfNeeded(genotypeFileName, ".geno"), Utils.addSuffixIfNeeded(MAPfileName, ".map"), theBLOBList, chromNames, tableArray, idLengths, numSites, numTaxa);

            Alignment[] align = new Alignment[chromNames.length];

            // create alignments
            for (int i = 0; i < theBLOBList.length; i++) {
                byte[][] alleleBLOB = null;
                //array for BLOBs holding coverage plus header [taxa][SNPs]
                byte[] positionBLOB = null;
                byte[] SNPidBLOB = null;
                if (theBLOBList[i][0].size() > 1) {
                    alleleBLOB = new byte[theBLOBList[i][0].size()][];
                    //array for BLOBs holding SNPs plus header [taxa][SNPs]
                    for (int j = 0; j < alleleBLOB.length; j++) {
                        alleleBLOB[j] = theBLOBList[i][0].get(j);
                    }
                    positionBLOB = (byte[])theBLOBList[i][1].get(0);
                    SNPidBLOB = theBLOBList[i][4].get(0);
                }
                align[i] = new Pack1Alignment(alleleBLOB, positionBLOB, SNPidBLOB);
            }
            return CombineAlignment.getInstance(align);
        }
        catch (Exception e) {
            throw new IllegalArgumentException("Problem creating Alignment: " + ExceptionUtils.getExceptionCauses(e));
        }
    }

    /**
     * Reads files in Flapjack format for maps with physical position.
     * @param genotypeFileName	the name and path of the Flapjack genotype file
     * @param MAPfileName	the name and path of the Flapjack map file
     * @param hasHetSeparator	boolean indicating that heterozygotes are separated by a character
     * @param hasNucleotides	boolean indicating that the marker values are nucleotides
     * @param missing	the missing value character
     * @param hetsep	the character separating the marker genotypes for heterozygotes 
     * @return an Alignment
     * @throws IOException
     */
    public static Alignment readFromFlapjackPhysical(String genotypeFileName, String MAPfileName, 
    		boolean hasHetSeparator, boolean hasNucleotides, String missing, String hetsep) throws IOException {

    	String input;
    	String[] info;
    	Pattern tab = Pattern.compile("\t");
    	Pattern hetpat;
    	if (hasHetSeparator) {
    		hetpat = Pattern.compile(hetsep);
    	} else {
    		hetpat = null;
    	}

    	class Snp {
    		final String id;
    		final String chr;
    		final int pos;
    		int col;

    		Snp(String id, String chr, int pos) {
    			this.id = id;
    			this.chr = chr;
    			this.pos = pos;
    		}
    	}

    	//read the map information (do not include snps without a valid position)
    	HashMap<String, Snp> SnpMap = new HashMap<String, Snp>();

    	BufferedReader br = new BufferedReader(new FileReader(MAPfileName));
    	while ((input = br.readLine()) != null) {
    		info = tab.split(input);
    		if (info.length > 2) {
    			try { 
    				int pos = Integer.parseInt(info[2]);
    				SnpMap.put(info[0], new Snp(info[0], info[1], pos));
    			} catch (NumberFormatException e) {}
    		}
    	}
    	br.close();

    	//read the genotype file

    	//create an array of snps and read in the taxa
    	TreeSet<String> locusSet = new TreeSet<String>();
    	br = new BufferedReader(new FileReader(genotypeFileName));
    	input = br.readLine(); //gets the snps
    	info = tab.split(input.trim());
    	int numberOfSnps = info.length;
    	Snp[] snps = new Snp[numberOfSnps];
    	for (int s = 0; s < numberOfSnps; s++) {
    		snps[s] = SnpMap.get(info[s]);
    		if (snps[s] != null) {
        		snps[s].col = s + 1;
    			locusSet.add(snps[s].chr);
    		}
    	}
    	int numberOfLoci = locusSet.size();

    	//index the snps by locus, position
    	//index is int[][] where locus is the first dimension and order by position is the second
    	ArrayList<Snp> snpList = new ArrayList<Snp>();
    	for (Snp s:snps) snpList.add(s);
    	Collections.sort(snpList, new Comparator<Snp>(){

    		@Override
    		public int compare(Snp s1, Snp s2) {
    			int locuscomp = s1.chr.compareTo(s2.chr);
    			if (locuscomp != 0) return locuscomp;
    			if (s1.pos < s2.pos) return -1;
    			if (s1.pos > s2.pos) return 1;
    			return 0;
    		}
    	});

    	//count the number of snps per locus
    	int[] snpcountPerLocus = new int[numberOfLoci];
    	for (int i = 0; i < numberOfLoci; i++) snpcountPerLocus[i] = 0;
    	int locusNumber = -1;
    	String prevchr = "";
    	for (Snp snp:snpList) {
    		if (!snp.chr.equals(prevchr)) {
    			locusNumber++;
    			prevchr = snp.chr;
    		}
    		snpcountPerLocus[locusNumber]++;
    	}

    	//set the snp index
    	int[][] snpIndex = new int[numberOfLoci][];
    	for (int locus = 0; locus < numberOfLoci; locus++) {
    		snpIndex[locus] = new int[snpcountPerLocus[locus]];
    	}
    	locusNumber = -1;
    	prevchr = "";
    	int snpNumber = 0;
    	for (Snp snp:snpList) {
    		if (!snp.chr.equals(prevchr)) {
    			locusNumber++;
    			prevchr = snp.chr;
    			snpNumber = 0;
    		}
    		snpIndex[locusNumber][snpNumber++] = snp.col;
    	}

    	//read in taxa
    	ArrayList<String> taxaList = new ArrayList<String>();
    	while((input = br.readLine()) != null) {
    		String taxon = input.substring(0, input.indexOf('\t'));
    		taxaList.add(taxon);
    	}
    	int numberOfTaxa = taxaList.size();
    	String[] taxanames = new String[numberOfTaxa];
    	taxaList.toArray(taxanames);
    	br.close();

    	DataType dt;
    	if (hasNucleotides) {
    		dt = new IUPACNucleotides();
    	} else {
    		dt = new TextDataType();
    	}

    	//set up the arrays to hold the alignment data
    	String[][] sequences = new String[numberOfLoci][numberOfTaxa];
    	int[][] positions = new int[numberOfLoci][];
    	String[][] siteNames = new String[numberOfLoci][];
    	for (int locus = 0; locus < numberOfLoci; locus++) {
    		String[] locusSiteNames = new String[snpcountPerLocus[locus]];
    		int[] locusPositions = new int[snpcountPerLocus[locus]];
    		for (int s = 0; s < snpcountPerLocus[locus]; s++) {
    			locusSiteNames[s] = snps[snpIndex[locus][s] - 1].id; //minus 1 because snpIndex is on snp column in the input file
    			locusPositions[s] = snps[snpIndex[locus][s] - 1].pos;
    		}
    		positions[locus] = locusPositions;
    		siteNames[locus] = locusSiteNames;
    	}

    	//read in the genotypes
    	br = new BufferedReader(new FileReader(genotypeFileName));
    	br.readLine();
    	int linecount = 0;
    	if (hasNucleotides) {
    		while ((input = br.readLine()) != null) {
    			info = tab.split(input);
    			for (int locus = 0; locus < numberOfLoci; locus++) {
    				StringBuilder sb = new StringBuilder();
    				for (int s = 0; s < snpcountPerLocus[locus]; s++) {
    					int col = snpIndex[locus][s];
    					String val = info[col];
    					if (val.equals(missing)) {
    						sb.append(IUPACNucleotides.UNKNOWN_CHARACTER);
    					} else if (val.length() == 1) {
    						sb.append(val);
    					} else if (hasHetSeparator) {
    						sb.append((char)IUPACNucleotides.getDegerateSNPByteFromTwoSNPs((byte) val.charAt(0), (byte) val.charAt(2)));
    					} else {
    						sb.append((char)IUPACNucleotides.getDegerateSNPByteFromTwoSNPs((byte) val.charAt(0), (byte) val.charAt(1)));
    					}
    				}
    				sequences[locus][linecount] = sb.toString();
    			}
    			linecount++;
    		}

    	} else {
    		TextDataType tdt = (TextDataType) dt;
    		while ((input = br.readLine()) != null) {

    			info = tab.split(input);
    			for (int locus = 0; locus < numberOfLoci; locus++) {
    				StringBuilder sb = new StringBuilder();
    				for (int s = 0; s < snpcountPerLocus[locus]; s++) {
    					int col = snpIndex[locus][s];
    					String val = info[col];
    					if (val.equals(missing)) {
    						sb.append(TextDataType.UNKNOWN_STRING);
    					} else if (val.length() == 1) {
    						sb.append(val);
    					} else if (hasHetSeparator) {
    						val = val.replace(hetsep, ":"); 
    						sb.append(tdt.getCharFromTextRepresentation(val));
    					} else {
    						val = new String(new char[]{val.charAt(0), ':', val.charAt(1)});
    						sb.append(tdt.getCharFromTextRepresentation(val));
    					}
    				}
    				sequences[locus][linecount] = sb.toString();
    			}
    			linecount++;
    		}
    	}
    	br.close();


    	//create the Alignments
    	Alignment[] theAlignments = new Alignment[numberOfLoci];
    	locusNumber = 0;
    	for (String locusName:locusSet) {
    		int lastsite = snpcountPerLocus[locusNumber] - 1;
    		Locus thisLocus = new Locus(locusName, locusName, positions[locusNumber][0], positions[locusNumber][lastsite],null, null);
    		theAlignments[locusNumber] = new SimpleAlignment(new SimpleIdGroup(taxanames), sequences[locusNumber], dt, null, positions[locusNumber], null, thisLocus,null, siteNames[locusNumber], false);
    		locusNumber++;
    	}

    	return CombineAlignment.getInstance(theAlignments);

    }
    
    /**
     * Reads files in Flapjack format for maps with genetic positions.
     * @param genotypeFileName	the name and path of the Flapjack genotype file
     * @param MAPfileName	the name and path of the Flapjack map file
     * @param hasHetSeparator	boolean indicating that heterozygotes are separated by a character
     * @param hasNucleotides	boolean indicating that the marker values are nucleotides
     * @param missing	the missing value character
     * @param hetsep	the character separating the marker genotypes for heterozygotes 
     * @return an Alignment
     * @throws IOException
     */
    public static Alignment readFromFlapjackGenetic(String genotypeFileName, String MAPfileName, 
    		boolean hasHetSeparator, boolean hasNucleotides, String missing, String hetsep) throws IOException {

    	String input;
    	String[] info;
    	Pattern tab = Pattern.compile("\t");
    	Pattern hetpat;
    	if (hasHetSeparator) {
    		hetpat = Pattern.compile(hetsep);
    	} else {
    		hetpat = null;
    	}

    	class Snp {
    		final String id;
    		final String chr;
    		final double pos;
    		int col;

    		Snp(String id, String chr, double pos) {
    			this.id = id;
    			this.chr = chr;
    			this.pos = pos;
    		}
    	}

    	//read the map information (do not include snps without a valid position)
    	HashMap<String, Snp> SnpMap = new HashMap<String, Snp>();

    	BufferedReader br = new BufferedReader(new FileReader(MAPfileName));
    	while ((input = br.readLine()) != null) {
    		info = tab.split(input);
    		if (info.length > 2) {
    			try { 
    				double pos = Double.parseDouble(info[2]);
    				SnpMap.put(info[0], new Snp(info[0], info[1], pos));
    			} catch (NumberFormatException e) {}
    		}
    	}
    	br.close();

    	//read the genotype file

    	//create an array of snps and read in the taxa
    	TreeSet<String> locusSet = new TreeSet<String>();
    	br = new BufferedReader(new FileReader(genotypeFileName));
    	input = br.readLine(); //gets the snps
    	info = tab.split(input.trim());
    	int numberOfSnps = info.length;
    	Snp[] snps = new Snp[numberOfSnps];
    	for (int s = 0; s < numberOfSnps; s++) {
    		snps[s] = SnpMap.get(info[s]);
    		if (snps[s] != null) {
        		snps[s].col = s + 1;
    			locusSet.add(snps[s].chr);
    		}
    	}
    	int numberOfLoci = locusSet.size();

    	//index the snps by locus, position
    	//index is int[][] where locus is the first dimension and order by position is the second
    	ArrayList<Snp> snpList = new ArrayList<Snp>();
    	for (Snp s:snps) snpList.add(s);
    	Collections.sort(snpList, new Comparator<Snp>(){

    		@Override
    		public int compare(Snp s1, Snp s2) {
    			int locuscomp = s1.chr.compareTo(s2.chr);
    			if (locuscomp != 0) return locuscomp;
    			if (s1.pos < s2.pos) return -1;
    			if (s1.pos > s2.pos) return 1;
    			return 0;
    		}
    	});

    	//count the number of snps per locus
    	int[] snpcountPerLocus = new int[numberOfLoci];
    	for (int i = 0; i < numberOfLoci; i++) snpcountPerLocus[i] = 0;
    	int locusNumber = -1;
    	String prevchr = "";
    	for (Snp snp:snpList) {
    		if (!snp.chr.equals(prevchr)) {
    			locusNumber++;
    			prevchr = snp.chr;
    		}
    		snpcountPerLocus[locusNumber]++;
    	}

    	//set the snp index
    	int[][] snpIndex = new int[numberOfLoci][];
    	for (int locus = 0; locus < numberOfLoci; locus++) {
    		snpIndex[locus] = new int[snpcountPerLocus[locus]];
    	}
    	locusNumber = -1;
    	prevchr = "";
    	int snpNumber = 0;
    	for (Snp snp:snpList) {
    		if (!snp.chr.equals(prevchr)) {
    			locusNumber++;
    			prevchr = snp.chr;
    			snpNumber = 0;
    		}
    		snpIndex[locusNumber][snpNumber++] = snp.col;
    	}

    	//read in taxa
    	ArrayList<String> taxaList = new ArrayList<String>();
    	while((input = br.readLine()) != null) {
    		String taxon = input.substring(0, input.indexOf('\t'));
    		taxaList.add(taxon);
    	}
    	int numberOfTaxa = taxaList.size();
    	String[] taxanames = new String[numberOfTaxa];
    	taxaList.toArray(taxanames);
    	br.close();

    	DataType dt;
    	if (hasNucleotides) {
    		dt = new IUPACNucleotides();
    	} else {
    		dt = new TextDataType();
    	}

    	//set up the arrays to hold the alignment data
    	String[][] sequences = new String[numberOfLoci][numberOfTaxa];
    	int[][] positions = new int[numberOfLoci][];
    	String[][] siteNames = new String[numberOfLoci][];
    	for (int locus = 0; locus < numberOfLoci; locus++) {
    		String[] locusSiteNames = new String[snpcountPerLocus[locus]];
    		int[] locusPositions = new int[snpcountPerLocus[locus]];
    		for (int s = 0; s < snpcountPerLocus[locus]; s++) {
    			locusSiteNames[s] = snps[snpIndex[locus][s] - 1].id; //minus 1 because snpIndex is on snp column in the input file
    			locusPositions[s] = (int)(1000 * snps[snpIndex[locus][s] - 1].pos);
    		}
    		positions[locus] = locusPositions;
    		siteNames[locus] = locusSiteNames;
    	}

    	//read in the genotypes
    	br = new BufferedReader(new FileReader(genotypeFileName));
    	br.readLine();
    	int linecount = 0;
    	if (hasNucleotides) {
    		while ((input = br.readLine()) != null) {
    			info = tab.split(input);
    			for (int locus = 0; locus < numberOfLoci; locus++) {
    				StringBuilder sb = new StringBuilder();
    				for (int s = 0; s < snpcountPerLocus[locus]; s++) {
    					int col = snpIndex[locus][s];
    					String val = info[col];
    					if (val.equals(missing)) {
    						sb.append(IUPACNucleotides.UNKNOWN_CHARACTER);
    					} else if (val.length() == 1) {
    						sb.append(val);
    					} else if (hasHetSeparator) {
    						sb.append((char)IUPACNucleotides.getDegerateSNPByteFromTwoSNPs((byte) val.charAt(0), (byte) val.charAt(2)));
    					} else {
    						sb.append((char)IUPACNucleotides.getDegerateSNPByteFromTwoSNPs((byte) val.charAt(0), (byte) val.charAt(1)));
    					}
    				}
    				sequences[locus][linecount] = sb.toString();
    			}
    			linecount++;
    		}

    	} else {
    		TextDataType tdt = (TextDataType) dt;
    		while ((input = br.readLine()) != null) {

    			info = tab.split(input);
    			for (int locus = 0; locus < numberOfLoci; locus++) {
    				StringBuilder sb = new StringBuilder();
    				for (int s = 0; s < snpcountPerLocus[locus]; s++) {
    					int col = snpIndex[locus][s];
    					String val = info[col];
    					if (val.equals(missing)) {
    						sb.append(TextDataType.UNKNOWN_STRING);
    					} else if (val.length() == 1) {
    						sb.append(val);
    					} else if (hasHetSeparator) {
    						val = val.replace(hetsep, ":"); 
    						sb.append(tdt.getCharFromTextRepresentation(val));
    					} else {
    						val = new String(new char[]{val.charAt(0), ':', val.charAt(1)});
    						sb.append(tdt.getCharFromTextRepresentation(val));
    					}
    				}
    				sequences[locus][linecount] = sb.toString();
    			}
    			linecount++;
    		}
    	}
    	br.close();


    	//create the Alignments
    	Alignment[] theAlignments = new Alignment[numberOfLoci];
    	locusNumber = 0;
    	for (String locusName:locusSet) {
    		int lastsite = snpcountPerLocus[locusNumber] - 1;
    		Locus thisLocus = new Locus(locusName, locusName, positions[locusNumber][0], positions[locusNumber][lastsite],null, null);
    		theAlignments[locusNumber] = new SimpleAlignment(new SimpleIdGroup(taxanames), sequences[locusNumber], dt, null, positions[locusNumber], null, thisLocus,null, siteNames[locusNumber], false);
    		locusNumber++;
    	}

    	return CombineAlignment.getInstance(theAlignments);

    }
    
    /**
     * Creates an array of alignments from a set of Flapjack files
     * @param genotypeFileName
     * @param MAPfileName
     * @return align an array of alignments from selected Flapjack files
     */
    public static Alignment readFromFlapjackOld(String genotypeFileName, String MAPfileName) {
        try {
            String[][] chromsAvailable = GdpdmBLOBUtils.getFileInfo(Utils.addSuffixIfNeeded(MAPfileName, ".map"), 1, 0);
//            int[] chromInfo;
            int numTaxa = GdpdmBLOBUtils.countLinesInFile(Utils.addSuffixIfNeeded(genotypeFileName, ".geno")) - 1;
            String[] chroms = new String[chromsAvailable.length];
            Alignment[] align = new Alignment[chromsAvailable.length];
//            ArrayList<Alignment> align = new ArrayList<Alignment>();
            //ArrayList<byte[]>[] theBLOBList = new ArrayList[GdpdmBLOBUtils.totalBLOBtypes];
            for (int i = 0; i < chromsAvailable.length; i++) {
                chroms[i] = chromsAvailable[i][0];
            }
//            chromInfo = GdpdmBLOBUtils.getChromInFileInfo(chromsAvailable, chrom);
            ArrayList<byte[]>[][] theBLOBList = new ArrayList[chroms.length][GdpdmBLOBUtils.totalBLOBtypes];
            for (int i = 0; i < theBLOBList.length; i++) {
                for (int j = 0; j < theBLOBList[i].length; j++) {
                    theBLOBList[i][j] = new ArrayList<byte[]>();
                }
            }
            theBLOBList = AllelePositionBLOBUtils.createAllelePositionBLOBsFromFlapjack(Utils.addSuffixIfNeeded(genotypeFileName, ".geno"), Utils.addSuffixIfNeeded(MAPfileName, ".map"), theBLOBList, chroms, chromsAvailable, numTaxa);
//            if (theBLOBList == null) {
//                theBLOBList = AllelePositionBLOBUtils.createAllelePositionBLOBsFromFlapjack(Utils.addSuffixIfNeeded(genotypeFileName, ".geno"), Utils.addSuffixIfNeeded(MAPfileName, ".map"), theBLOBList, chroms, chromsAvailable, numTaxa);
//            }
            for (int i = 0; i < theBLOBList.length; i++) {
                byte[][] alleleBLOB = null;
                //array for BLOBs holding coverage plus header [taxa][SNPs]
                byte[] positionBLOB = null;
                byte[] SNPidBLOB = null;
                if (theBLOBList[i][0].size() > 1) {
                    alleleBLOB = new byte[theBLOBList[i][0].size()][];
                    //array for BLOBs holding SNPs plus header [taxa][SNPs]
                    for (int j = 0; j < alleleBLOB.length; j++) {
                        alleleBLOB[j] = theBLOBList[i][0].get(j);
                    }
                    positionBLOB = (byte[])theBLOBList[i][1].get(0);
                    SNPidBLOB = theBLOBList[i][4].get(0);
                }
                align[i] = new Pack1Alignment(alleleBLOB, positionBLOB, SNPidBLOB);
//                Alignment a = new Pack1Alignment(alleleBLOB, positionBLOB, SNPidBLOB);
//                if (a != null) {
//                    align.add(a);
//                }
            }
//            Alignment[] aligns = new Alignment[align.size()];
//            for (int i = 0; i < align.size(); i++) {
//                aligns[i] = align.get(i);
//            }
            //return align;
            return CombineAlignment.getInstance(align);
        }
        catch (Exception e) {
            throw new IllegalArgumentException("Problem creating Alignment: " + ExceptionUtils.getExceptionCauses(e));
        }
    }

    public static Alignment createPack1AlignmentFromFile(ArrayList<String> fileList) {
        ArrayList<byte[]>[] theBLOBList = new ArrayList[GdpdmBLOBUtils.totalBLOBtypes];
        for (int i = 0; i < theBLOBList.length; i++) {
            theBLOBList[i] = new ArrayList<byte[]>();
        }
        byte[][] coverageBLOB = null, alleleBLOB = null;
        //array for BLOBs holding coverage plus header [taxa][SNPs]
        byte[] refSeqBLOB = null, positionBLOB = null, SNPidBLOB = null, qualityBLOB = null;
        try {
            for (String theFileName : fileList) {
                if (theFileName.endsWith(".gz")) {
                    theBLOBList = GdpdmBLOBUtils.readBLOBfromGZIP(theFileName, theBLOBList);
                } else if (theFileName.endsWith(".7z")) {
                    theBLOBList = GdpdmBLOBUtils.readBLOBfromLZMA(theFileName, theBLOBList);
                } else if (theFileName.endsWith(".txt")) {
                    String[][] chromsAvailable = GdpdmBLOBUtils.getFileInfo(theFileName, 2, 0);
                    String chrom = chromsAvailable[1][0];
                    int[] chromInfo = GdpdmBLOBUtils.getChromInFileInfo(chromsAvailable, chrom);
                    theBLOBList = AllelePositionBLOBUtils.createAllelePositionBLOBsFromHapmap(theFileName, theBLOBList, chrom, chromInfo);
                } else {
                    throw new IllegalArgumentException("Unknown file type.  Must be .map, .ped, .txt, .zip, or .7z.");
                }
            }
        } catch (Exception e) {
            throw new IllegalArgumentException("Problem reading Alignment: " + ExceptionUtils.getExceptionCauses(e));
        }
        if ((theBLOBList == null)) {
            throw new IllegalArgumentException("File empty.");
        }
        else if ((theBLOBList.length == 0)) {
            throw new IllegalArgumentException("No data in this file.");
        }
        //todo put reading coverage here.
        if (theBLOBList[0].size() > 1) {
            alleleBLOB = new byte[theBLOBList[0].size()][];
            //array for BLOBs holding SNPs plus header [taxa][SNPs]
            for (int i = 0; i < alleleBLOB.length; i++) {
                alleleBLOB[i] = theBLOBList[0].get(i);
            }
            positionBLOB = (byte[]) theBLOBList[1].get(0);
            if (theBLOBList[4].size() == 0) {
                SNPidBLOB = new byte[1030];
                SNPidBLOB[GdpdmBLOBUtils.blobTypeField[0]] = 0x35;
            }
            else {
                SNPidBLOB = theBLOBList[4].get(0);
            }
          //  SNPidBLOB = theBLOBList[4].get(0);
        }
        if (theBLOBList[2].size() > 0) {
            coverageBLOB = new byte[theBLOBList[2].size()][];
            for (int i = 0; i < alleleBLOB.length; i++) {
                coverageBLOB[i] = theBLOBList[2].get(i);
            }
        }
        if(!theBLOBList[6].isEmpty()) {
            qualityBLOB=(byte[]) theBLOBList[6].get(0);
        }
        try {
            if (coverageBLOB == null) {
                return new Pack1Alignment(alleleBLOB, positionBLOB, SNPidBLOB, qualityBLOB);
            } else {
                return new Pack1AlignmentWithCoverage(alleleBLOB, positionBLOB, SNPidBLOB, coverageBLOB, refSeqBLOB);
            }
        } catch (Exception e) {
            throw new IllegalArgumentException("Problem creating Alignment: " + ExceptionUtils.getExceptionCauses(e));
        }

    }

    public static Alignment createPack1AlignmentFromFile(String alleleFileName, String coverageFileName, String refseqFileName) {
        ArrayList<String> fileList = new ArrayList<String>();
        if (alleleFileName.length() != 0) {
            fileList.add(alleleFileName);
        }
        if ((coverageFileName != null) && (coverageFileName.length() != 0)) {
            fileList.add(coverageFileName);
        }
        if ((refseqFileName != null) && (refseqFileName.length() != 0)) {
            fileList.add(refseqFileName);
        }
        return createPack1AlignmentFromFile(fileList);
    }
}
