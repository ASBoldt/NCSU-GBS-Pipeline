/*
 * AllelePositionBLOBUtils
 */
package net.maizegenetics.pal.alignment;

import java.util.Hashtable;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.Random;
import java.util.regex.Pattern;
import net.maizegenetics.util.ExceptionUtils;
import net.maizegenetics.util.Utils;

/**
 * These utils are used for the creation create a set of BLOB for genotypic data including methods for:
 * 1. converting a HapMap text file to zipped alleleBLOBs
 * 2. creating a pack1alignment from alleleBLOBs
 * 3. Methods for extracting individuals bases and ranges of bases
 * 4. Encoding information for the IUPAC code into half byte
 *
 * @author ed
 */
public class AllelePositionBLOBUtils {

    //HapMap format constants
    public static final int totalNonTaxaHeaders = 11;
    //Plink .ped format constants
    public static final int totalNonSNPHeaders = 6;
    //IUPAC code equivalents
    public static final byte[] AIUPAC = {(byte) 'A'};
    public static final byte[] CIUPAC = {(byte) 'C'};
    public static final byte[] GIUPAC = {(byte) 'G'};
    public static final byte[] TIUPAC = {(byte) 'T'};
    public static final byte[] RIUPAC = {(byte) 'A', (byte) 'G'};
    public static final byte[] YIUPAC = {(byte) 'C', (byte) 'T'};
    public static final byte[] SIUPAC = {(byte) 'G', (byte) 'C'};
    public static final byte[] WIUPAC = {(byte) 'A', (byte) 'T'};
    public static final byte[] KIUPAC = {(byte) 'G', (byte) 'T'};
    public static final byte[] MIUPAC = {(byte) 'A', (byte) 'C'};
    public static final byte[] BIUPAC = {(byte) '+'};
    public static final byte[] DIUPAC = {(byte) '+', (byte) '-'};
    public static final byte[] HIUPAC = {};
    public static final byte[] VIUPAC = {};
    public static final byte[] NIUPAC = {(byte) 'N'};
    public static final byte[] gapIUPAC = {(byte) '-'};
    // Compiled Whitespace Pattern
    private static Pattern WHITESPACE_PATTERN = Pattern.compile("\\s");

    private AllelePositionBLOBUtils() {
    }

    public static IdGroup getTaxa(byte[][] alleleBLOB) {
        Identifier[] ids = new Identifier[alleleBLOB.length];
        for (int i = 0; i < alleleBLOB.length; i++) {
            String tn = new String(alleleBLOB[i], GdpdmBLOBUtils.taxaField[0], GdpdmBLOBUtils.taxaField[1]);
            ids[i] = new Identifier(tn.trim());
        }
        return new SimpleIdGroup(ids);
    }

    /**
     * Unpacks a base from an allele BLOB which uses half byte coding
     * @param alleleBLOB
     * @param site
     * @return an ASCII byte representation of IUPAC character
     */
    public static byte getBaseFromAlleleBLOB(byte[] alleleBLOB, int site) {
        int index = site / 2;
        // int currBase=alleleBLOB[index];
        if (site % 2 == 0) {
            return GdpdmBLOBUtils.bases[(alleleBLOB[index + GdpdmBLOBUtils.totalHeaderPadding] >> 4) & 0xf];
        } else {
            return GdpdmBLOBUtils.bases[alleleBLOB[index + GdpdmBLOBUtils.totalHeaderPadding] & 0xf];
        }
    }

    /**
     * Unpacks a base from an allele BLOB which uses half byte coding
     * @param alleleBLOB
     * @param startSite
     * @param endSite
     * @return an ASCII byte representation of IUPAC character
     */
    public static byte[] getBaseRangeFromAlleleBLOB(byte[] alleleBLOB, int startSite, int endSite) {
        //this is roughly 3x times faster at get bases then individually.
        int totalSites = endSite - startSite + 1;
        int currPos = 0;
        byte[] allele = new byte[totalSites];
        if (startSite % 2 == 1) {
            allele[0] = getBaseFromAlleleBLOB(alleleBLOB, startSite);
            startSite++;
            currPos++;
        }
        if (endSite % 2 == 0) {
            allele[totalSites - 1] = getBaseFromAlleleBLOB(alleleBLOB, endSite);
            endSite--;
            totalSites--;
        }
        ByteBuffer bb = ByteBuffer.wrap(alleleBLOB);
        bb.position((startSite / 2) + GdpdmBLOBUtils.totalHeaderPadding);
        while (currPos < totalSites) {
            //todo need to check whether i am missing the last bases
            byte b = bb.get();
            allele[currPos] = GdpdmBLOBUtils.bases[(b >> 4) & 0xF];
            currPos++;
            allele[currPos] = GdpdmBLOBUtils.bases[b & 0xF];
            currPos++;
        }

        return allele;
    }

    public static Pack1Alignment permuteAlignment(Pack1Alignment p1a) {
        Random generator = new Random();
        byte[][] permAlleleBLOBs = new byte[p1a.alleleBLOB.length][];
        for (int i = 0; i < p1a.getSequenceCount(); i++) {
            permAlleleBLOBs[i] = (byte[]) p1a.alleleBLOB[i].clone();
        }
        int randTaxa = 0;
        for (int b = 0; b < p1a.getSiteCount(); b++) {
            for (int i = 0; i < p1a.getSequenceCount(); i++) {
                randTaxa = generator.nextInt(p1a.getSequenceCount());
                byte b1 = getBaseFromAlleleBLOB(permAlleleBLOBs[i], b);
                byte b2 = getBaseFromAlleleBLOB(permAlleleBLOBs[randTaxa], b);
                setHalfByteInAlleleBLOB(permAlleleBLOBs[i], b, (char) b2);
                setHalfByteInAlleleBLOB(permAlleleBLOBs[randTaxa], b, (char) b1);
            }
        }
        Pack1Alignment newP1A = new Pack1Alignment(permAlleleBLOBs, p1a.getVariableSitesBLOB(), p1a.getSNPidBLOB());
        return newP1A;

    }

    /**
     * This will make a hard copy of an alignment.
     * @param p1a Input aniglment
     * @param newIdGroup
     * @return
     */
    public static Pack1Alignment subsetCopyAlignment(Pack1Alignment p1a, IdGroup newIdGroup) {
        byte[][] permAlleleBLOBs = new byte[newIdGroup.getIdCount()][];
        for (int i = 0; i < newIdGroup.getIdCount(); i++) {
            int p1aSeqNumber = p1a.getIdGroup().whichIdNumber(newIdGroup.getIdentifier(i).getFullName());
            permAlleleBLOBs[i] = (byte[]) p1a.alleleBLOB[p1aSeqNumber].clone();
        }
        Pack1Alignment newP1A = new Pack1Alignment(permAlleleBLOBs, p1a.getVariableSitesBLOB(), p1a.getSNPidBLOB());
        return newP1A;

    }

    /**
     * This will make a hard copy of an alignment with all taxa.
     * @param p1a Input aniglment
     * @return
     */
    public static Pack1Alignment copyAlignment(Pack1Alignment p1a) {
        return subsetCopyAlignment(p1a, p1a.getIdGroup());
    }

    /**
     * Create a set of packed allele BLOBs and one position BLOB from a HapMap text file
     * @param infileName -  HapMap text file
     * @param theBLOBList - the BLOBs being added to
     * @param chrom - the desired chromosome
     * @param chromInfo - [0] = number of sites, [1] = depth in file, [2] = id length
     * @return createZip - the arrays of blobs needed to build an alignment
     */
    public static ArrayList<byte[]>[] createAllelePositionBLOBsFromHapmap(String infileName, ArrayList<byte[]>[] theBLOBList, String chrom, int[] chromInfo) {
        int numSites = 0;  //number of size in the alignment
        int numTaxa; //number of taxa or individuals
        byte[][] alleleB; //array for BLOBs holding SNPs plus header [taxa][SNPs]
        byte[] locusPositionB;  //array for BLOB holding the positions
        byte[] SNPidB; //array for BLOB holding SNP ids
        byte[] qualityB; //array for BLOB holding SNP ids
        int minPosition = Integer.MAX_VALUE, maxPosition = Integer.MIN_VALUE;  //smallest and largest QTL position
        String currLocus = "ERROR", currGenomeBuild = "ERROR"; //name of the locus, name of the genome build

        try {
            //int[] sites = GdpdmBLOBUtils.numSitesInChrom(GdpdmBLOBUtils.getChromsAvailableCounts(infileName, 2), chrom);
            numSites = chromInfo[0];
            if (numSites == 0) {
                throw new Exception("Desired chromosome not available in " + infileName + "."); //chrom not available
            }
            int depthInFile = chromInfo[1]; //how many lines into .txt, data for desired chromosome begins
            int SNPidLength = chromInfo[2];
            //numSites = GdpdmBLOBUtils.countLinesInFile(infileName) - 1;  //assume there is a one line header
            //int SNPidLength = GdpdmBLOBUtils.getSNPIDLength(infileName, depthInFile, numSites, 0, true); //find max SNP id length
            //BufferedReader fileIn = new BufferedReader(new FileReader(infileName), 1000000);
            BufferedReader fileIn = GdpdmBLOBUtils.getBufferedReader(infileName, 1000000);
            //BufferedReader fileIn = new BufferedReader(new InputStreamReader((new URL(infileName)).openStream()), 1000000);
            String[] header = WHITESPACE_PATTERN.split(fileIn.readLine());
            //String[] header = fileIn.readLine().split("\t");
            numTaxa = header.length - totalNonTaxaHeaders;
            int packSiteSize = (numSites % 2 == 1) ? ((numSites + 1) / 2) : (numSites / 2); //number of bytes to pack bases, add an extra site if odd
            alleleB = new byte[numTaxa][GdpdmBLOBUtils.totalHeaderPadding + packSiteSize];
            locusPositionB = new byte[GdpdmBLOBUtils.totalHeaderPadding + ((Integer.SIZE / Byte.SIZE) * numSites)];
            SNPidB = new byte[GdpdmBLOBUtils.totalHeaderPadding + (SNPidLength * numSites)];
            qualityB = new byte[GdpdmBLOBUtils.totalHeaderPadding + numSites];
            ByteBuffer lpB = ByteBuffer.wrap(locusPositionB);
            ByteBuffer sidB = ByteBuffer.wrap(SNPidB);
            ByteBuffer qBB = ByteBuffer.wrap(qualityB);
            lpB.position(GdpdmBLOBUtils.totalHeaderPadding);
            sidB.position(GdpdmBLOBUtils.totalHeaderPadding);
            qBB.position(GdpdmBLOBUtils.totalHeaderPadding);
            int prevPosition = -1;
            for (int depth = 0; depth < depthInFile - 1; depth++) {
                fileIn.readLine(); //shift .txt reader to desired position in .txt file
            }
            for (int site = 0; site < numSites; site++) {
                String[] s = WHITESPACE_PATTERN.split(fileIn.readLine());
                int position = Integer.parseInt(s[3]);
                if (position < prevPosition) {
                    throw new Exception("Sites are not properly sorted at " + position + " and " + prevPosition);
                }
                if (site == 0) {
                    currLocus = s[2];
                    currGenomeBuild = s[5];
                    minPosition = position;
                }
                for (int i = 0; i < numTaxa; i++) {
                    if (s[totalNonTaxaHeaders + i].length() == 2) {
                        setHalfByteInAlleleBLOB(alleleB[i], site, (char) getBaseFromSNPValue(s[totalNonTaxaHeaders + i].getBytes()));
                    } else {
                        setHalfByteInAlleleBLOB(alleleB[i], site, s[totalNonTaxaHeaders + i].charAt(0));
                    }
                }
                sidB.put(s[0].getBytes()); //get bytes representing SNP id, insert into SNP id BLOB
                lpB.putInt(position);
                try {
                    qBB.put(Byte.parseByte(s[10]));
                } catch (NumberFormatException nfe) {
                    qBB.put((byte) -1);
                }
                sidB.position(GdpdmBLOBUtils.totalHeaderPadding + (SNPidLength * (site + 1))); //increment snp id BLOB buffer correctly
                prevPosition = position;
            }
            maxPosition = prevPosition;
            for (int i = 0; i < numTaxa; i++) {
                GdpdmBLOBUtils.setHeaderOfBLOB(alleleB[i], GdpdmBLOBUtils.compressionAlgorithm, GdpdmBLOBUtils.alleleBLOBtype,
                        numSites, currGenomeBuild,
                        currLocus, minPosition, maxPosition, header[i + totalNonTaxaHeaders], 4);
            }
            GdpdmBLOBUtils.setHeaderOfBLOB(locusPositionB, GdpdmBLOBUtils.compressionAlgorithm, GdpdmBLOBUtils.allelePositionBLOBtype,
                    numSites, currGenomeBuild,
                    currLocus, minPosition, maxPosition, null, 32);
            GdpdmBLOBUtils.setHeaderOfBLOB(SNPidB, GdpdmBLOBUtils.compressionAlgorithm, GdpdmBLOBUtils.SNPIdBLOBtype,
                    numSites, currGenomeBuild,
                    currLocus, minPosition, maxPosition, null, SNPidLength * 8);
//                    currLocus, minPosition, maxPosition, null, SNPidLength);
            GdpdmBLOBUtils.setHeaderOfBLOB(qualityB, GdpdmBLOBUtils.compressionAlgorithm, GdpdmBLOBUtils.siteQualityBLOBtype,
                    numSites, currGenomeBuild, currLocus, minPosition, maxPosition, null, 1);
            for (byte[] ab : alleleB) {
                int blobType = Integer.parseInt("" + (char) ab[GdpdmBLOBUtils.blobTypeField[0]]);
                theBLOBList[blobType - 1].add(ab);
            }
            int lpBlobType = Integer.parseInt("" + (char) locusPositionB[GdpdmBLOBUtils.blobTypeField[0]]);
            theBLOBList[lpBlobType - 1].add(locusPositionB);
            int snpidBLOBType = Integer.parseInt("" + (char) SNPidB[GdpdmBLOBUtils.blobTypeField[0]]);
            theBLOBList[snpidBLOBType - 1].add(SNPidB);
            int qualtyBLOBType = Integer.parseInt("" + (char) qualityB[GdpdmBLOBUtils.blobTypeField[0]]);
            theBLOBList[qualtyBLOBType - 1].add(qualityB);
            return theBLOBList;
        } catch (Exception e) {
            System.err.println("File IO in createAlleleBLOBs: " + e);
        }
        return null;
    }

    /**
     * Create a set of packed allele BLOBs and one position BLOB from a HapMap text file
     * @param infileName -  HapMap text file
     * @param theBLOBList - the BLOBs being added to
     * @param chrom - the desired chromosome
     * @param chromInfo - [0] = number of sites, [1] = depth in file, [2] = id length
     * @return createZip - the arrays of blobs needed to build an alignment
     */
    public static ArrayList<byte[]>[][] createAllelePositionBLOBsFromHapmap(String infileName, ArrayList<byte[]>[][] theBLOBList, String[] chroms, String[][] chromsAvailable) throws IOException{
        int numSites = 0;  //number of size in the alignment
        int numTaxa; //number of taxa or individuals
        byte[][] alleleB; //array for BLOBs holding SNPs plus header [taxa][SNPs]
        byte[] locusPositionB;  //array for BLOB holding the positions
        byte[] SNPidB; //array for BLOB holding SNP ids
        byte[] qualityB; //array for BLOB holding SNP ids
        int minPosition = Integer.MAX_VALUE, maxPosition = Integer.MIN_VALUE;  //smallest and largest QTL position
        String currLocus = "ERROR", currGenomeBuild = "ERROR"; //name of the locus, name of the genome build

        BufferedReader fileIn = null;
        try {
            //BufferedReader fileIn = new BufferedReader(new FileReader(infileName), 1000000);
            //fileIn = GdpdmBLOBUtils.getBufferedReader(infileName, 1000000);
            fileIn = Utils.getBufferedReader(infileName, 1000000);
            int lineInFile = 0;
            //BufferedReader fileIn = new BufferedReader(new InputStreamReader((new URL(infileName)).openStream()), 1000000);

            String[] header = WHITESPACE_PATTERN.split(fileIn.readLine());
            lineInFile++;
            //int[] sites = GdpdmBLOBUtils.numSitesInChrom(GdpdmBLOBUtils.getChromsAvailableCounts(infileName, 2), chrom);
            int depthInFile = 1;
            for (int x = 0; x < chroms.length; x++) {
                int[] chromInfo = GdpdmBLOBUtils.getChromInFileInfo(chromsAvailable, chroms[x]);
                numSites = chromInfo[0];
                if (numSites == 0) {
                    throw new IllegalArgumentException("Desired chromosome not available in " + infileName + "."); //chrom not available
                }
//                depthInFile += chromInfo[1]; //how many lines into .txt, data for desired chromosome begins
                int SNPidLength = chromInfo[2];
                //numSites = GdpdmBLOBUtils.countLinesInFile(infileName) - 1;  //assume there is a one line header
                //int SNPidLength = GdpdmBLOBUtils.getSNPIDLength(infileName, depthInFile, numSites, 0, true); //find max SNP id length
                //String[] header = fileIn.readLine().split("\t");
                numTaxa = header.length - totalNonTaxaHeaders;
                int packSiteSize = (numSites % 2 == 1) ? ((numSites + 1) / 2) : (numSites / 2); //number of bytes to pack bases, add an extra site if odd
                alleleB = new byte[numTaxa][GdpdmBLOBUtils.totalHeaderPadding + packSiteSize];
                locusPositionB = new byte[GdpdmBLOBUtils.totalHeaderPadding + ((Integer.SIZE / Byte.SIZE) * numSites)];
                SNPidB = new byte[GdpdmBLOBUtils.totalHeaderPadding + (SNPidLength * numSites)];
                qualityB = new byte[GdpdmBLOBUtils.totalHeaderPadding + numSites];
                ByteBuffer lpB = ByteBuffer.wrap(locusPositionB);
                ByteBuffer sidB = ByteBuffer.wrap(SNPidB);
                ByteBuffer qBB = ByteBuffer.wrap(qualityB);
                lpB.position(GdpdmBLOBUtils.totalHeaderPadding);
                sidB.position(GdpdmBLOBUtils.totalHeaderPadding);
                qBB.position(GdpdmBLOBUtils.totalHeaderPadding);
                int prevPosition = -1;
                while (depthInFile < chromInfo[1]) {
                    fileIn.readLine();
                    depthInFile++;
                    lineInFile++;
                }
//                for (int depth = 0; depth < depthInFile - 1; depth++) {
//                    fileIn.readLine(); //shift .txt reader to desired position in .txt file
//                }
                for (int site = 0; site < numSites; site++) {
                    String[] s = WHITESPACE_PATTERN.split(fileIn.readLine());
                    depthInFile++;
                    lineInFile++;
                    int position = Integer.parseInt(s[3]);
                    if (position < prevPosition) {
                        throw new IllegalArgumentException("Sites are not properly sorted at " + position + " and " + prevPosition);
                    }
                    if (site == 0) {
                        currLocus = s[2];
                        currGenomeBuild = s[5];
                        minPosition = position;
                    }
                    if (numTaxa + totalNonTaxaHeaders != s.length) {
                        throw new IllegalStateException("Number of Taxa: " + numTaxa + " does not match number of values: " + (s.length - totalNonTaxaHeaders) + " at line in file: " + lineInFile + " chromosome: " + chroms[x] + " site: " + site);
                    }
                    for (int i = 0; i < numTaxa; i++) {
                        if (s[totalNonTaxaHeaders + i].length() == 2) {
                            setHalfByteInAlleleBLOB(alleleB[i], site, (char) getBaseFromSNPValue(s[totalNonTaxaHeaders + i].getBytes()));
                        } else {
                            setHalfByteInAlleleBLOB(alleleB[i], site, s[totalNonTaxaHeaders + i].charAt(0));
                        }
                    }
                    sidB.put(s[0].getBytes()); //get bytes representing SNP id, insert into SNP id BLOB
                    lpB.putInt(position);
                    sidB.position(GdpdmBLOBUtils.totalHeaderPadding + (SNPidLength * (site + 1))); //increment snp id BLOB buffer correctly
                    prevPosition = position;
                }
                maxPosition = prevPosition;
                for (int i = 0; i < numTaxa; i++) {
                    GdpdmBLOBUtils.setHeaderOfBLOB(alleleB[i], GdpdmBLOBUtils.compressionAlgorithm, GdpdmBLOBUtils.alleleBLOBtype,
                            numSites, currGenomeBuild,
                            currLocus, minPosition, maxPosition, header[i + totalNonTaxaHeaders], 4);
                }
                GdpdmBLOBUtils.setHeaderOfBLOB(locusPositionB, GdpdmBLOBUtils.compressionAlgorithm, GdpdmBLOBUtils.allelePositionBLOBtype,
                        numSites, currGenomeBuild,
                        currLocus, minPosition, maxPosition, null, 32);
                GdpdmBLOBUtils.setHeaderOfBLOB(SNPidB, GdpdmBLOBUtils.compressionAlgorithm, GdpdmBLOBUtils.SNPIdBLOBtype,
                        numSites, currGenomeBuild,
                        currLocus, minPosition, maxPosition, null, SNPidLength * 8);
//                        currLocus, minPosition, maxPosition, null, SNPidLength);
                for (byte[] ab : alleleB) {
                    int blobType = Integer.parseInt("" + (char) ab[GdpdmBLOBUtils.blobTypeField[0]]);
                    theBLOBList[x][blobType - 1].add(ab);
                }
                int lpBlobType = Integer.parseInt("" + (char) locusPositionB[GdpdmBLOBUtils.blobTypeField[0]]);
                theBLOBList[x][lpBlobType - 1].add(locusPositionB);
                int snpidBLOBType = Integer.parseInt("" + (char) SNPidB[GdpdmBLOBUtils.blobTypeField[0]]);
                theBLOBList[x][snpidBLOBType - 1].add(SNPidB);
            }
            return theBLOBList;
        } finally {
            try {
                fileIn.close();
            } catch (Exception e) {
                // do nothing
            }
        }
    }

    /**
     * Create a set of packed allele BLOBs and one position BLOB from Plink files, .PED and .MAP.
     * Currently user has to input desired chromosome from file, otherwise error is thrown.
     * @param PEDfileName - .PED file
     * @param MAPfileName - .MAP file
     * @param theBLOBList - the BLOBs being added to
     * @param chrom - name of desired chromosome
     * @param chromInfo - [0] = number of sites, [1] = depth in file, [2] = id length
     * @param numTaxa - number of individuals
     * @return theBLOBList - the arrays of blobs needed to build an alignment
     */
    public static ArrayList<byte[]>[] createAllelePositionBLOBsFromPlink(String PEDfileName, String MAPfileName, ArrayList<byte[]>[] theBLOBList, String chrom, int[] chromInfo, int numTaxa) {
        int numSites;  //number of size in the alignment
        byte[][] alleleB; //array for BLOBs holding SNPs plus header [taxa][SNPs]
        byte[] locusPositionB;  //array for BLOB holding the positions
        byte[] SNPidB; //array for BLOB holding SNP ids
        int minPosition = Integer.MAX_VALUE, maxPosition = Integer.MIN_VALUE;  //smallest and largest QTL position
        String currLocus = "ERROR", currGenomeBuild = "ERROR"; //name of the locus, name of the genome build

        try {
            //BufferedReader PEDfileIn = new BufferedReader(new FileReader(PEDfileName), 1000000);
            //BufferedReader MAPfileIn = new BufferedReader(new FileReader(MAPfileName), 1000000);
            BufferedReader PEDfileIn = GdpdmBLOBUtils.getBufferedReader(PEDfileName, 1000000);
            BufferedReader MAPfileIn = GdpdmBLOBUtils.getBufferedReader(MAPfileName, 1000000);
            //BufferedReader PEDfileIn = new BufferedReader(new InputStreamReader((new URL(PEDfileName)).openStream()), 1000000);
            //BufferedReader MAPfileIn = new BufferedReader(new InputStreamReader((new URL(MAPfileName)).openStream()), 1000000);
            numSites = chromInfo[0];
            int depthInFile = chromInfo[1]; //how many lines into .MAP, data for desired chromosome begins
            int SNPidLength = chromInfo[2];
            if (numSites == 0) {
                throw new Exception("Desired chromosome not available in " + MAPfileName + "."); //chrom not available
            }
            int packSiteSize = (numSites % 2 == 1) ? ((numSites + 1) / 2) : (numSites / 2); //number of bytes to pack bases, add an extra site if odd
            alleleB = new byte[numTaxa][GdpdmBLOBUtils.totalHeaderPadding + packSiteSize];
            locusPositionB = new byte[GdpdmBLOBUtils.totalHeaderPadding + ((Integer.SIZE / Byte.SIZE) * numSites)];
            SNPidB = new byte[GdpdmBLOBUtils.totalHeaderPadding + (SNPidLength * numSites)];
            ByteBuffer lpB = ByteBuffer.wrap(locusPositionB);
            ByteBuffer sidB = ByteBuffer.wrap(SNPidB);
            lpB.position(GdpdmBLOBUtils.totalHeaderPadding);
            sidB.position(GdpdmBLOBUtils.totalHeaderPadding);
            for (int depth = 0; depth < depthInFile; depth++) {
                MAPfileIn.readLine(); //shift .map reader to desired position in .map file
            }
            int prevPosition = -1;
            for (int site = 0; site < numSites; site++) {
                String[] mapLine = WHITESPACE_PATTERN.split(MAPfileIn.readLine());
                int position = Integer.parseInt(mapLine[3]);
                if (position < prevPosition) {
                    throw new Exception("Sites are not properly sorted at " + position + " and " + prevPosition);
                }
                if (site == 0) {
                    currLocus = mapLine[0]; //chrom name is first column in .map file
                    currGenomeBuild = "NA"; //genome build not available in Plink format
                    minPosition = position;
                }
                sidB.put(mapLine[1].getBytes()); //get bytes representing SNP id, insert into SNP id BLOB
                lpB.putInt(position);
                sidB.position(GdpdmBLOBUtils.totalHeaderPadding + (SNPidLength * (site + 1))); //increment snp id BLOB buffer correctly
                prevPosition = position;
            }
            maxPosition = prevPosition;
            GdpdmBLOBUtils.setHeaderOfBLOB(locusPositionB, GdpdmBLOBUtils.compressionAlgorithm,
                    GdpdmBLOBUtils.allelePositionBLOBtype, numSites, currGenomeBuild, currLocus,
                    minPosition, maxPosition, null, 32);
            GdpdmBLOBUtils.setHeaderOfBLOB(SNPidB, GdpdmBLOBUtils.compressionAlgorithm,
                    GdpdmBLOBUtils.SNPIdBLOBtype, numSites, currGenomeBuild, currLocus,
                    minPosition, maxPosition, null, SNPidLength * 8);
//                    minPosition, maxPosition, null, SNPidLength);
            byte[] SNPValue = new byte[2];
            for (int taxa = 0; taxa < numTaxa; taxa++) {
                String[] pedLine = WHITESPACE_PATTERN.split(PEDfileIn.readLine());
                StringBuilder taxonName = new StringBuilder(pedLine[1]); //creating taxonName, namelvl0 = IndividID
                taxonName.append(":");
                taxonName.append(pedLine[0]); //creating taxonName, namelvl1 = FamilyID
                GdpdmBLOBUtils.setHeaderOfBLOB(alleleB[taxa], GdpdmBLOBUtils.compressionAlgorithm,
                        GdpdmBLOBUtils.alleleBLOBtype, numSites, currGenomeBuild, currLocus,
                        minPosition, maxPosition, taxonName.toString(), 4);
                for (int site = 0; site < numSites; site++) {
                    //SNP data in .ped file is in pairs (data always diploid) however pairs are still
                    //whitespace delimited, pairs need to be partnered up after splitting
                    SNPValue[0] = getBaseFromNumeric(pedLine[totalNonSNPHeaders + depthInFile * 2 + site * 2].charAt(0));
                    SNPValue[1] = getBaseFromNumeric(pedLine[totalNonSNPHeaders + depthInFile * 2 + site * 2 + 1].charAt(0));
                    setHalfByteInAlleleBLOB(alleleB[taxa], site, (char) getBaseFromSNPValue(SNPValue));
                }
            }
            for (byte[] ab : alleleB) {
                int blobType = Integer.parseInt("" + (char) ab[GdpdmBLOBUtils.blobTypeField[0]]);
                theBLOBList[blobType - 1].add(ab);
            }
            int lpBlobType = Integer.parseInt("" + (char) locusPositionB[GdpdmBLOBUtils.blobTypeField[0]]);
            theBLOBList[lpBlobType - 1].add(locusPositionB);
            int snpidBLOBType = Integer.parseInt("" + (char) SNPidB[GdpdmBLOBUtils.blobTypeField[0]]);
            theBLOBList[snpidBLOBType - 1].add(SNPidB);
            PEDfileIn.close();
            MAPfileIn.close();
            return theBLOBList;
        } catch (Exception e) {
            System.err.println("File IO in createAlleleBLOBsFromPlink: " + e);
        }
        return null;
    }

    /**
     * Create a set of packed allele BLOBs and one position BLOB from Plink files, .PED and .MAP.
     * Currently user has to input desired chromosome from file, otherwise error is thrown.
     * @param PEDfileName - .PED file
     * @param MAPfileName - .MAP file
     * @param theBLOBList - the BLOBs being added to
     * @param chrom - name of desired chromosome
     * @param chromInfo - [0] = number of sites, [1] = depth in file, [2] = id length
     * @param numTaxa - number of individuals
     * @return theBLOBList - the arrays of blobs needed to build an alignment
     */
    public static ArrayList<byte[]>[][] createAllelePositionBLOBsFromPlink(String PEDfileName, String MAPfileName, ArrayList<byte[]>[][] theBLOBList, String[] chroms, String[][] chromsAvailable, int numTaxa) throws IOException{
        int numSites;  //number of size in the alignment
        byte[][] alleleB; //array for BLOBs holding SNPs plus header [taxa][SNPs]
        byte[] locusPositionB;  //array for BLOB holding the positions
        byte[] SNPidB; //array for BLOB holding SNP ids
        int minPosition = Integer.MAX_VALUE, maxPosition = Integer.MIN_VALUE;  //smallest and largest QTL position
        String currLocus = "ERROR", currGenomeBuild = "ERROR"; //name of the locus, name of the genome build

        //BufferedReader PEDfileIn = new BufferedReader(new FileReader(PEDfileName), 1000000);
        //BufferedReader MAPfileIn = new BufferedReader(new FileReader(MAPfileName), 1000000);
//            BufferedReader PEDfileIn = GdpdmBLOBUtils.getBufferedReader(PEDfileName, 1000000);
        BufferedReader MAPfileIn = GdpdmBLOBUtils.getBufferedReader(MAPfileName, 1000000);
        //BufferedReader PEDfileIn = new BufferedReader(new InputStreamReader((new URL(PEDfileName)).openStream()), 1000000);
        //BufferedReader MAPfileIn = new BufferedReader(new InputStreamReader((new URL(MAPfileName)).openStream()), 1000000);
        int depthInFile = 0;
        for (int x = 0; x < chroms.length; x++) {
            BufferedReader PEDfileIn = GdpdmBLOBUtils.getBufferedReader(PEDfileName, 1000000);
            int[] chromInfo = GdpdmBLOBUtils.getChromInFileInfo(chromsAvailable, chroms[x]);
            numSites = chromInfo[0];
//                int depthInFile = chromInfo[1]; //how many lines into .MAP, data for desired chromosome begins
            int SNPidLength = chromInfo[2];
            if (numSites == 0) {
                throw new IllegalStateException("Desired chromosome not available in " + MAPfileName + "."); //chrom not available
            }
            int packSiteSize = (numSites % 2 == 1) ? ((numSites + 1) / 2) : (numSites / 2); //number of bytes to pack bases, add an extra site if odd
            alleleB = new byte[numTaxa][GdpdmBLOBUtils.totalHeaderPadding + packSiteSize];
            locusPositionB = new byte[GdpdmBLOBUtils.totalHeaderPadding + ((Integer.SIZE / Byte.SIZE) * numSites)];
            SNPidB = new byte[GdpdmBLOBUtils.totalHeaderPadding + (SNPidLength * numSites)];
            ByteBuffer lpB = ByteBuffer.wrap(locusPositionB);
            ByteBuffer sidB = ByteBuffer.wrap(SNPidB);
            lpB.position(GdpdmBLOBUtils.totalHeaderPadding);
            sidB.position(GdpdmBLOBUtils.totalHeaderPadding);
            while (depthInFile < chromInfo[1]) {
                MAPfileIn.readLine();
                depthInFile++;
            }
//                for (int depth = 0; depth < depthInFile; depth++) {
//                    MAPfileIn.readLine(); //shift .map reader to desired position in .map file
//                }
            int prevPosition = -1;
            for (int site = 0; site < numSites; site++) {
                String[] mapLine = WHITESPACE_PATTERN.split(MAPfileIn.readLine());
                depthInFile++;
                int position = Integer.parseInt(mapLine[3]);
                if (position < prevPosition) {
                    throw new IllegalStateException("Sites are not properly sorted at " + position + " and " + prevPosition);
                }
                if (site == 0) {
                    currLocus = mapLine[0]; //chrom name is first column in .map file
                    currGenomeBuild = "NA"; //genome build not available in Plink format
                    minPosition = position;
                }
                sidB.put(mapLine[1].getBytes()); //get bytes representing SNP id, insert into SNP id BLOB
                lpB.putInt(position);
                sidB.position(GdpdmBLOBUtils.totalHeaderPadding + (SNPidLength * (site + 1))); //increment snp id BLOB buffer correctly
                prevPosition = position;
            }
            maxPosition = prevPosition;
            GdpdmBLOBUtils.setHeaderOfBLOB(locusPositionB, GdpdmBLOBUtils.compressionAlgorithm,
                    GdpdmBLOBUtils.allelePositionBLOBtype, numSites, currGenomeBuild, currLocus,
                    minPosition, maxPosition, null, 32);
            GdpdmBLOBUtils.setHeaderOfBLOB(SNPidB, GdpdmBLOBUtils.compressionAlgorithm,
                    GdpdmBLOBUtils.SNPIdBLOBtype, numSites, currGenomeBuild, currLocus,
                    minPosition, maxPosition, null, SNPidLength * 8);
//                        minPosition, maxPosition, null, SNPidLength);
            byte[] SNPValue = new byte[2];
//                System.out.println(numTaxa);
            for (int taxa = 0; taxa < numTaxa; taxa++) {
//                    System.out.println("testing1");
                String[] pedLine = WHITESPACE_PATTERN.split(PEDfileIn.readLine());
//                    System.out.println(pedLine[1]);
                StringBuilder taxonName = new StringBuilder(pedLine[1]); //creating taxonName, namelvl0 = IndividID
                taxonName.append(":");
                taxonName.append(pedLine[0]); //creating taxonName, namelvl1 = FamilyID
                GdpdmBLOBUtils.setHeaderOfBLOB(alleleB[taxa], GdpdmBLOBUtils.compressionAlgorithm,
                        GdpdmBLOBUtils.alleleBLOBtype, numSites, currGenomeBuild, currLocus,
                        minPosition, maxPosition, taxonName.toString(), 4);
//                    System.out.println(numSites);
//                    System.out.println(depthInFile);
                for (int site = 0; site < numSites; site++) {
                    //SNP data in .ped file is in pairs (data always diploid) however pairs are still
                    //whitespace delimited, pairs need to be partnered up after splitting
                    SNPValue[0] = getBaseFromNumeric(pedLine[totalNonSNPHeaders + (depthInFile - numSites) * 2 + site * 2].charAt(0));
                    SNPValue[1] = getBaseFromNumeric(pedLine[totalNonSNPHeaders + (depthInFile - numSites) * 2 + site * 2 + 1].charAt(0));
                    setHalfByteInAlleleBLOB(alleleB[taxa], site, (char) getBaseFromSNPValue(SNPValue));
                }
            }
            for (byte[] ab : alleleB) {
                int blobType = Integer.parseInt("" + (char) ab[GdpdmBLOBUtils.blobTypeField[0]]);
                theBLOBList[x][blobType - 1].add(ab);
            }
            int lpBlobType = Integer.parseInt("" + (char) locusPositionB[GdpdmBLOBUtils.blobTypeField[0]]);
            theBLOBList[x][lpBlobType - 1].add(locusPositionB);
            int snpidBLOBType = Integer.parseInt("" + (char) SNPidB[GdpdmBLOBUtils.blobTypeField[0]]);
            theBLOBList[x][snpidBLOBType - 1].add(SNPidB);
            PEDfileIn.close();
        }
//            PEDfileIn.close();
        MAPfileIn.close();
        return theBLOBList;

    }

    /**
     * Create a set of packed allele BLOBs and one position BLOB from a set of Flapjack files.
     * Currently user has to input desired chromosome from file, otherwise error is thrown.
     * @param genotypeFileName - genotype file
     * @param MAPfileName - .MAP file
     * @param theBLOBList - the BLOBs being added to
     * @param chrom - name of desired chromosome
     * @param chromInfo - [0] = number of sites, [1] = depth in file, [2] = id length
     * @param numTaxa - number of individuals
     * @return theBLOBList - the arrays of blobs needed to build an alignment
     */
    public static ArrayList<byte[]>[] createAllelePositionBLOBsFromFlapjack(String genotypeFileName, String MAPfileName, ArrayList<byte[]>[] theBLOBList, String chrom, int[] chromInfo, int numTaxa) {
        int numSites;  //number of size in the alignment
        byte[][] alleleB; //array for BLOBs holding SNPs plus header [taxa][SNPs]
        byte[] locusPositionB;  //array for BLOB holding the positions
        byte[] SNPidB; //array for BLOB holding SNP ids
        int minPosition = Integer.MAX_VALUE, maxPosition = Integer.MIN_VALUE;  //smallest and largest QTL position
        String currLocus = "ERROR", currGenomeBuild = "ERROR"; //name of the locus, name of the genome build

        try {
            //BufferedReader genotypeFileIn = new BufferedReader(new FileReader(genotypeFileName), 1000000);
            //BufferedReader MAPfileIn = new BufferedReader(new FileReader(MAPfileName), 1000000);
            BufferedReader genotypeFileIn = GdpdmBLOBUtils.getBufferedReader(genotypeFileName, 1000000);
            BufferedReader MAPfileIn = GdpdmBLOBUtils.getBufferedReader(MAPfileName, 1000000);
            //BufferedReader genotypeFileIn = new BufferedReader(new InputStreamReader((new URL(genotypeFileName)).openStream()), 1000000);
            //BufferedReader MAPfileIn = new BufferedReader(new InputStreamReader((new URL(MAPfileName)).openStream()), 1000000);
            //numTaxa = GdpdmBLOBUtils.countLinesInFile(genotypeFileName) - 1;
            //int[] chroms = GdpdmBLOBUtils.numSitesInChrom(GdpdmBLOBUtils.getChromsAvailableCounts(MAPfileName, 1), chrom);
            numSites = chromInfo[0];
            if (numSites == 0) {
                throw new Exception("Desired chromosome not available in " + MAPfileName + "."); //chrom not available in .map file
            }
            int depthInFile = chromInfo[1]; //how many lines into .MAP, data for desired chromosome begins
            int SNPidLength = chromInfo[2]; //GdpdmBLOBUtils.getSNPIDLength(MAPfileName, depthInFile, numSites, 0, false); //find max SNP id length
            for (int depth = 0; depth < depthInFile; depth++) {
                MAPfileIn.readLine(); //shift .map reader to desired position in .map file
            }
            int packSiteSize = (numSites % 2 == 1) ? ((numSites + 1) / 2) : (numSites / 2); //number of bytes to pack bases, add an extra site if odd
            alleleB = new byte[numTaxa][GdpdmBLOBUtils.totalHeaderPadding + packSiteSize];
            locusPositionB = new byte[GdpdmBLOBUtils.totalHeaderPadding + ((Integer.SIZE / Byte.SIZE) * numSites)];
            SNPidB = new byte[GdpdmBLOBUtils.totalHeaderPadding + (SNPidLength * numSites)];
            ByteBuffer lpB = ByteBuffer.wrap(locusPositionB);
            ByteBuffer sidB = ByteBuffer.wrap(SNPidB);
            lpB.position(GdpdmBLOBUtils.totalHeaderPadding);
            sidB.position(GdpdmBLOBUtils.totalHeaderPadding);
            String[] sites = WHITESPACE_PATTERN.split(genotypeFileIn.readLine());
            String[] mapLine = WHITESPACE_PATTERN.split(MAPfileIn.readLine());
            int startLoc = 1;
            int prevPosition = -1;
            for (int i = 0; i < sites.length - 1; i++) { // find correct column in genotype file desired data starts
                if (sites[i + 1].equals(mapLine[0])) {
                    for (int site = 0; site < numSites; site++) {
                        int position = (int) Double.parseDouble(mapLine[2]); // possible third column is genetic distance not position
                        if (position < prevPosition) {
                            throw new Exception("Sites are not properly sorted at " + position + " and " + prevPosition);
                        }
                        if (site == 0) {
                            currLocus = mapLine[1]; //chrom name is second column in .map file
                            currGenomeBuild = "NA"; //genome build not available in Flapjack format
                            minPosition = position;
                        }
                        sidB.put(mapLine[0].getBytes()); //get bytes representing SNP id, insert into SNP id BLOB
                        lpB.putInt(position);
                        sidB.position(GdpdmBLOBUtils.totalHeaderPadding + (SNPidLength * (site + 1))); //increment snp id BLOB buffer correctly
                        prevPosition = position;
                        if (site != numSites - 1) {
                            mapLine = WHITESPACE_PATTERN.split(MAPfileIn.readLine());
                        }
                    }
                    break;
                }
                if (i == sites.length - 2) {
//                    if (getSingleChrom) {
                    throw new Exception("Desired chromosome not available in " + genotypeFileName + "."); //chrom not available in genotype file
//                    }
//                    else {
//                        return null;
//                    }
                }
                startLoc++;
            }
            maxPosition = prevPosition;
            GdpdmBLOBUtils.setHeaderOfBLOB(locusPositionB, GdpdmBLOBUtils.compressionAlgorithm,
                    GdpdmBLOBUtils.allelePositionBLOBtype, numSites, currGenomeBuild, currLocus,
                    minPosition, maxPosition, null, 32);
            GdpdmBLOBUtils.setHeaderOfBLOB(SNPidB, GdpdmBLOBUtils.compressionAlgorithm,
                    GdpdmBLOBUtils.SNPIdBLOBtype, numSites, currGenomeBuild, currLocus,
                    minPosition, maxPosition, null, SNPidLength * 8);
//                    minPosition, maxPosition, null, SNPidLength);
            byte[] SNPValue = new byte[2];
            for (int taxa = 0; taxa < numTaxa; taxa++) {
                String[] genotypeLine = WHITESPACE_PATTERN.split(genotypeFileIn.readLine());
                GdpdmBLOBUtils.setHeaderOfBLOB(alleleB[taxa], GdpdmBLOBUtils.compressionAlgorithm,
                        GdpdmBLOBUtils.alleleBLOBtype, numSites, currGenomeBuild, currLocus,
                        minPosition, maxPosition, genotypeLine[0], 4);
                for (int site = 0; site < numSites; site++) {
                    //SNP data in flapjack either, single letter, X/X, or empty string for no data
                    String snp = genotypeLine[startLoc + site];
                    byte[] SNPs = snp.getBytes();
                    if (SNPs.length == 1) {
                        SNPValue[0] = SNPs[0];
                        SNPValue[1] = SNPs[0];
                    } else if (SNPs.length == 3) {
                        SNPValue[0] = SNPs[0];
                        SNPValue[1] = SNPs[2];
                    } else if (SNPs.length == 0) {
                        SNPValue[0] = (byte) '-';
                        SNPValue[1] = (byte) '-';
                    } else {
                        throw new Exception(snp + " is not in a supported SNP format.");
                    }
                    setHalfByteInAlleleBLOB(alleleB[taxa], site, (char) getBaseFromSNPValue(SNPValue));
                }
            }
            for (byte[] ab : alleleB) {
                int blobType = Integer.parseInt("" + (char) ab[GdpdmBLOBUtils.blobTypeField[0]]);
                theBLOBList[blobType - 1].add(ab);
            }
            int lpBlobType = Integer.parseInt("" + (char) locusPositionB[GdpdmBLOBUtils.blobTypeField[0]]);
            theBLOBList[lpBlobType - 1].add(locusPositionB);
            int snpidBLOBType = Integer.parseInt("" + (char) SNPidB[GdpdmBLOBUtils.blobTypeField[0]]);
            theBLOBList[snpidBLOBType - 1].add(SNPidB);
            genotypeFileIn.close();
            MAPfileIn.close();
            return theBLOBList;
        } catch (Exception e) {
            System.err.println("File IO in createAlleleBLOBsFromFlapjack: " + e);
        }
        return null;
    }

    /**
     * Create allele, position, and id BLOBs from Flapjack file, one set of BLOBs per chromosome
     * Flapjack file does not need to be sorted by position or chromosome
     * @param genotypeFileName
     * @param MAPfileName
     * @param theBLOBList
     * @param chromNames
     * @param tableArray
     * @param idLengths
     * @param numSites
     * @param numTaxa
     * @return theBLOBList
     */
    public static ArrayList<byte[]>[][] createAllelePositionBLOBsFromFlapjack(String genotypeFileName, String MAPfileName, ArrayList<byte[]>[][] theBLOBList, String[] chromNames, Hashtable[] tableArray, int[] idLengths, int[] numSites, int numTaxa) {
        int numChroms = chromNames.length;
        int[] minPositions = new int[numChroms];
        int[] maxPositions = new int[numChroms];
        String[] taxaNames = new String[numTaxa];
        byte[][][] alleleB = new byte[numChroms][numTaxa][]; //array for BLOBs holding SNPs plus header [chromosome][taxa][SNPs]
        byte[][] locusPositionB = new byte[numChroms][];  //array for BLOB holding the positions [chromosome][position]
        byte[][] SNPidB = new byte[numChroms][]; //array for BLOB holding SNP ids [chromosome][IDs]

        try {
            int[] packSiteSize = new int[numChroms];
            ByteBuffer[] lpB = new ByteBuffer[numChroms];
            ByteBuffer[] sidB = new ByteBuffer[numChroms];

            // initialzing byte arrays
            for (int i = 0; i < numChroms; i++) {
                packSiteSize[i] = (numSites[i] % 2 == 1) ? ((numSites[i] + 1) / 2) : (numSites[i] / 2); //number of bytes to pack bases, add an extra site if odd
                for (int j = 0; j < numTaxa; j++) {
                    alleleB[i][j] = new byte[GdpdmBLOBUtils.totalHeaderPadding + packSiteSize[i]];
                }
                locusPositionB[i] = new byte[GdpdmBLOBUtils.totalHeaderPadding + ((Integer.SIZE / Byte.SIZE) * numSites[i])];
                SNPidB[i] = new byte[GdpdmBLOBUtils.totalHeaderPadding + (idLengths[i] * numSites[i])];
                lpB[i] = ByteBuffer.wrap(locusPositionB[i]);
                sidB[i] = ByteBuffer.wrap(SNPidB[i]);
                lpB[i].position(GdpdmBLOBUtils.totalHeaderPadding);
                sidB[i].position(GdpdmBLOBUtils.totalHeaderPadding);
            }

            BufferedReader MAPfileIn = GdpdmBLOBUtils.getBufferedReader(MAPfileName, 1000000);
            String mapFileLine = "";

            while ((mapFileLine = MAPfileIn.readLine()) != null) {
                String[] fileLineArray = WHITESPACE_PATTERN.split(mapFileLine);
                int chromIndex = -1;
                for (int i = 0; i < numChroms; i++) {
                    if (fileLineArray[1].equals(chromNames[i])) {
                        chromIndex = i;
                        i = numChroms;
                    }
                }

                // finding order of SNPs by position
                int position = Integer.parseInt(fileLineArray[2]);
                int positionIndex = (Integer) tableArray[chromIndex].get(position);

                // setting min and max positions for BLOB headers
                if (positionIndex == 0) {
                    minPositions[chromIndex] = position;
                } else if (positionIndex == (numSites[chromIndex] - 1)) {
                    maxPositions[chromIndex] = position;
                }

                BufferedReader genotypeFileIn = GdpdmBLOBUtils.getBufferedReader(genotypeFileName, 1000000);
                String[] sites = WHITESPACE_PATTERN.split(genotypeFileIn.readLine());

                // find SNP location in file
                int snpIndex = -1;
                for (int i = 1; i < sites.length; i++) {
                    if (sites[i].equals(fileLineArray[0])) {
                        snpIndex = i;
                        i = sites.length;
                    }
                }

                // populate position and id BLOB body
                lpB[chromIndex].position(GdpdmBLOBUtils.totalHeaderPadding + (4 * (positionIndex)));
                lpB[chromIndex].putInt(position);
                sidB[chromIndex].position(GdpdmBLOBUtils.totalHeaderPadding + (idLengths[chromIndex] * (positionIndex)));
                sidB[chromIndex].put(fileLineArray[0].getBytes());

                // populate SNP value BLOB body
                for (int i = 0; i < numTaxa; i++) {
                    String[] genoFileLine = WHITESPACE_PATTERN.split(genotypeFileIn.readLine());

                    // saving taxa names
                    taxaNames[i] = genoFileLine[0];

                    byte[] SNPValue = new byte[2];
                    byte[] SNPs = genoFileLine[snpIndex].getBytes();
                    if (SNPs.length == 1) {
                        SNPValue[0] = SNPs[0];
                        SNPValue[1] = SNPs[0];
                    } else if (SNPs.length == 3) {
                        SNPValue[0] = SNPs[0];
                        SNPValue[1] = SNPs[2];
                    } else if (SNPs.length == 0) {
                        SNPValue[0] = (byte) '-';
                        SNPValue[1] = (byte) '-';
                    } else {
                        throw new Exception(genoFileLine[snpIndex] + " is not in a supported SNP format.");
                    }
                    setHalfByteInAlleleBLOB(alleleB[chromIndex][i], positionIndex, (char) getBaseFromSNPValue(SNPValue));
                }

                genotypeFileIn.close();
            }

            MAPfileIn.close();

            // populate BLOB headers
            for (int i = 0; i < numChroms; i++) {
                GdpdmBLOBUtils.setHeaderOfBLOB(locusPositionB[i], GdpdmBLOBUtils.compressionAlgorithm,
                        GdpdmBLOBUtils.allelePositionBLOBtype, numSites[i], "NA", chromNames[i],
                        minPositions[i], maxPositions[i], null, 32);
                GdpdmBLOBUtils.setHeaderOfBLOB(SNPidB[i], GdpdmBLOBUtils.compressionAlgorithm,
                        GdpdmBLOBUtils.SNPIdBLOBtype, numSites[i], "NA", chromNames[i],
                        minPositions[i], maxPositions[i], null, idLengths[i] * 8);
                for (int j = 0; j < numTaxa; j++) {
                    GdpdmBLOBUtils.setHeaderOfBLOB(alleleB[i][j], GdpdmBLOBUtils.compressionAlgorithm,
                            GdpdmBLOBUtils.alleleBLOBtype, numSites[i], "NA", chromNames[i],
                            minPositions[i], maxPositions[i], taxaNames[j], 4);
                }
            }

            // populating theBLOBList
            for (int i = 0; i < numChroms; i++) {
                int lpBlobType = Integer.parseInt("" + (char) locusPositionB[i][GdpdmBLOBUtils.blobTypeField[0]]);
                theBLOBList[i][lpBlobType - 1].add(locusPositionB[i]);
                int snpidBLOBType = Integer.parseInt("" + (char) SNPidB[i][GdpdmBLOBUtils.blobTypeField[0]]);
                theBLOBList[i][snpidBLOBType - 1].add(SNPidB[i]);
                for (int j = 0; j < numTaxa; j++) {
                    int blobType = Integer.parseInt("" + (char) alleleB[i][j][GdpdmBLOBUtils.blobTypeField[0]]);
                    theBLOBList[i][blobType - 1].add(alleleB[i][j]);
                }
            }

            return theBLOBList;
        } catch (Exception e) {
            throw new IllegalArgumentException("Problem creating Alignment from Flapjack files: " + ExceptionUtils.getExceptionCauses(e));
        }
    }

    /**
     * Create a set of packed allele BLOBs and one position BLOB from a set of Flapjack files.
     * Currently user has to input desired chromosome from file, otherwise error is thrown.
     * @param genotypeFileName - genotype file
     * @param MAPfileName - .MAP file
     * @param theBLOBList - the BLOBs being added to
     * @param chrom - name of desired chromosome
     * @param chromInfo - [0] = number of sites, [1] = depth in file, [2] = id length
     * @param numTaxa - number of individuals
     * @return theBLOBList - the arrays of blobs needed to build an alignment
     */
    public static ArrayList<byte[]>[][] createAllelePositionBLOBsFromFlapjack(String genotypeFileName, String MAPfileName, ArrayList<byte[]>[][] theBLOBList, String[] chroms, String[][] chromsAvailable, int numTaxa) {
        int numSites;  //number of size in the alignment
        byte[][] alleleB; //array for BLOBs holding SNPs plus header [taxa][SNPs]
        byte[] locusPositionB;  //array for BLOB holding the positions
        byte[] SNPidB; //array for BLOB holding SNP ids
        int minPosition = Integer.MAX_VALUE, maxPosition = Integer.MIN_VALUE;  //smallest and largest QTL position
        String currLocus = "ERROR", currGenomeBuild = "ERROR"; //name of the locus, name of the genome build

        try {
            //BufferedReader genotypeFileIn = new BufferedReader(new FileReader(genotypeFileName), 1000000);
            //BufferedReader MAPfileIn = new BufferedReader(new FileReader(MAPfileName), 1000000);
//            BufferedReader genotypeFileIn = GdpdmBLOBUtils.getBufferedReader(genotypeFileName, 1000000);
            BufferedReader MAPfileIn = GdpdmBLOBUtils.getBufferedReader(MAPfileName, 1000000);
            //BufferedReader genotypeFileIn = new BufferedReader(new InputStreamReader((new URL(genotypeFileName)).openStream()), 1000000);
            //BufferedReader MAPfileIn = new BufferedReader(new InputStreamReader((new URL(MAPfileName)).openStream()), 1000000);
            //numTaxa = GdpdmBLOBUtils.countLinesInFile(genotypeFileName) - 1;
            //int[] chroms = GdpdmBLOBUtils.numSitesInChrom(GdpdmBLOBUtils.getChromsAvailableCounts(MAPfileName, 1), chrom);
            int depthInFile = 0;
            for (int x = 0; x < chroms.length; x++) {
                BufferedReader genotypeFileIn = GdpdmBLOBUtils.getBufferedReader(genotypeFileName, 1000000);
                int[] chromInfo = GdpdmBLOBUtils.getChromInFileInfo(chromsAvailable, chroms[x]);
                numSites = chromInfo[0];
                if (numSites == 0) {
                    throw new Exception("Desired chromosome not available in " + MAPfileName + "."); //chrom not available in .map file
                }
//                int depthInFile = chromInfo[1]; //how many lines into .MAP, data for desired chromosome begins
                int SNPidLength = chromInfo[2]; //GdpdmBLOBUtils.getSNPIDLength(MAPfileName, depthInFile, numSites, 0, false); //find max SNP id length
//                for (int depth = 0; depth < depthInFile; depth++) {
//                    MAPfileIn.readLine(); //shift .map reader to desired position in .map file
//                }
                while (depthInFile < chromInfo[1]) {
                    MAPfileIn.readLine();
                    depthInFile++;
                }
                int packSiteSize = (numSites % 2 == 1) ? ((numSites + 1) / 2) : (numSites / 2); //number of bytes to pack bases, add an extra site if odd
                alleleB = new byte[numTaxa][GdpdmBLOBUtils.totalHeaderPadding + packSiteSize];
                locusPositionB = new byte[GdpdmBLOBUtils.totalHeaderPadding + ((Integer.SIZE / Byte.SIZE) * numSites)];
                SNPidB = new byte[GdpdmBLOBUtils.totalHeaderPadding + (SNPidLength * numSites)];
                ByteBuffer lpB = ByteBuffer.wrap(locusPositionB);
                ByteBuffer sidB = ByteBuffer.wrap(SNPidB);
                lpB.position(GdpdmBLOBUtils.totalHeaderPadding);
                sidB.position(GdpdmBLOBUtils.totalHeaderPadding);
                String[] sites = WHITESPACE_PATTERN.split(genotypeFileIn.readLine());
                String[] mapLine = WHITESPACE_PATTERN.split(MAPfileIn.readLine());
                depthInFile++;
                int startLoc = 1;
                int prevPosition = -1;
                for (int i = 0; i < sites.length - 1; i++) { // find correct column in genotype file desired data starts
                    if (sites[i + 1].equals(mapLine[0])) {
                        for (int site = 0; site < numSites; site++) {
                            int position = (int) Double.parseDouble(mapLine[2]); // possible third column is genetic distance not position
                            if (position < prevPosition) {
                                throw new Exception("Sites are not properly sorted at " + position + " and " + prevPosition);
                            }
                            if (site == 0) {
                                currLocus = mapLine[1]; //chrom name is second column in .map file
                                currGenomeBuild = "NA"; //genome build not available in Flapjack format
                                minPosition = position;
                            }
                            sidB.put(mapLine[0].getBytes()); //get bytes representing SNP id, insert into SNP id BLOB
                            lpB.putInt(position);
                            sidB.position(GdpdmBLOBUtils.totalHeaderPadding + (SNPidLength * (site + 1))); //increment snp id BLOB buffer correctly
                            prevPosition = position;
                            if (site != numSites - 1) {
                                mapLine = WHITESPACE_PATTERN.split(MAPfileIn.readLine());
                                depthInFile++;
                            }
                        }
                        break;
                    }
                    if (i == sites.length - 2) {
                        //                    if (getSingleChrom) {
                        throw new Exception("Desired chromosome not available in " + genotypeFileName + "."); //chrom not available in genotype file
                        //                    }
                        //                    else {
                        //                        return null;
                        //                    }
                    }
                    startLoc++;
                }
                maxPosition = prevPosition;
                GdpdmBLOBUtils.setHeaderOfBLOB(locusPositionB, GdpdmBLOBUtils.compressionAlgorithm,
                        GdpdmBLOBUtils.allelePositionBLOBtype, numSites, currGenomeBuild, currLocus,
                        minPosition, maxPosition, null, 32);
                GdpdmBLOBUtils.setHeaderOfBLOB(SNPidB, GdpdmBLOBUtils.compressionAlgorithm,
                        GdpdmBLOBUtils.SNPIdBLOBtype, numSites, currGenomeBuild, currLocus,
                        minPosition, maxPosition, null, SNPidLength * 8);
//                        minPosition, maxPosition, null, SNPidLength);
                byte[] SNPValue = new byte[2];
                for (int taxa = 0; taxa < numTaxa; taxa++) {
                    String[] genotypeLine = WHITESPACE_PATTERN.split(genotypeFileIn.readLine());
                    GdpdmBLOBUtils.setHeaderOfBLOB(alleleB[taxa], GdpdmBLOBUtils.compressionAlgorithm,
                            GdpdmBLOBUtils.alleleBLOBtype, numSites, currGenomeBuild, currLocus,
                            minPosition, maxPosition, genotypeLine[0], 4);
                    for (int site = 0; site < numSites; site++) {
                        //SNP data in flapjack either, single letter, X/X, or empty string for no data
                        String snp = genotypeLine[startLoc + site];
                        byte[] SNPs = snp.getBytes();
                        if (SNPs.length == 1) {
                            SNPValue[0] = SNPs[0];
                            SNPValue[1] = SNPs[0];
                        } else if (SNPs.length == 3) {
                            SNPValue[0] = SNPs[0];
                            SNPValue[1] = SNPs[2];
                        } else if (SNPs.length == 0) {
                            SNPValue[0] = (byte) '-';
                            SNPValue[1] = (byte) '-';
                        } else {
                            throw new Exception(snp + " is not in a supported SNP format.");
                        }
                        setHalfByteInAlleleBLOB(alleleB[taxa], site, (char) getBaseFromSNPValue(SNPValue));
                    }
                }
                for (byte[] ab : alleleB) {
                    int blobType = Integer.parseInt("" + (char) ab[GdpdmBLOBUtils.blobTypeField[0]]);
                    theBLOBList[x][blobType - 1].add(ab);
                }
                int lpBlobType = Integer.parseInt("" + (char) locusPositionB[GdpdmBLOBUtils.blobTypeField[0]]);
                theBLOBList[x][lpBlobType - 1].add(locusPositionB);
                int snpidBLOBType = Integer.parseInt("" + (char) SNPidB[GdpdmBLOBUtils.blobTypeField[0]]);
                theBLOBList[x][snpidBLOBType - 1].add(SNPidB);
                genotypeFileIn.close();
            }
//            genotypeFileIn.close();
            MAPfileIn.close();
            return theBLOBList;
        } catch (Exception e) {
            System.err.println("File IO in createAlleleBLOBsFromFlapjack: " + e);
        }
        return null;
    }

    /**
     * Converstion of possible use of 1/2/3/4 into A/C/G/T in .PED file
     * @param num
     * @return correct A/C/G/T
     */
    public static byte getBaseFromNumeric(char num) {
        switch (num) {
            case '0':
                return (byte) 'N';
            case '1':
                return (byte) 'A';
            case '2':
                return (byte) 'C';
            case '3':
                return (byte) 'G';
            case '4':
                return (byte) 'T';
            case 'D':
                return (byte) '0';
            default:
                return (byte) num; //in case .PED file uses A/C/G/T
        }
    }

    public static byte[] getSNPValueForPlink(byte[] base) {
        byte[] SNPValues = new byte[base.length];
        for (int i = 0; i < base.length; i++) {
            if ((char) base[i] == 'N') {
                SNPValues[i] = (byte) '0';
            } else if ((char) base[i] == '0') {
                SNPValues[i] = (byte) 'D';
            } else {
                SNPValues[i] = base[i];
            }
        }
        return SNPValues;
    }

    /**
     * Conversion of sequence character into SNP Values.
     * This uese the standard IUPAC, except gap is 15 Hex F
     * @param base
     * @return character array of length 1 to 3 based on sequence character
     */
    public static byte[] getSNPValueFromBase(char base) {
        return getSNPValueFromHalfByte(getHalfByteFromBase(base));
    }

    /**
     * Conversion of SNP Values into a sequence character.
     * This uses the standard IUPAC, except gap is 15 or Hex F
     * @param snpValue
     * @return sequence character
     */
    public static byte getBaseFromSNPValue(byte[] snpValue) {
        return getBaseFromHalfByte(getHalfByteFromSNPValue(snpValue));
    }

    /**
     * Overloaded method calling getBaseFromSNPValue twice to process two data states
     * @param snpValue1
     * @param snpValue2
     * @return byte array, first, second element is conversion of first, second
     * SNP Value respectively
     */
    public static byte[] getBaseFromSNPValue(byte[] snpValue1, byte[] snpValue2) {
        byte[] twoBases = {getBaseFromSNPValue(snpValue1), getBaseFromSNPValue(snpValue2)};
        return twoBases;
    }

    /**
     * Overloaded method calling getBaseFromHalfByte twice to process two data states
     * @param halfByte1
     * @param halfByte2
     * @return byte array, first, second element is conversion of first, second
     * HalfBytes respectively
     */
    public static byte[] getBaseFromHalfByte(byte halfByte1, byte halfByte2) {
        byte[] twoBases = {getBaseFromHalfByte(halfByte1), getBaseFromHalfByte(halfByte2)};
        return twoBases;
    }

    /**
     * Conversion of half byte (4bit) encoding into SNP Values.
     * This uses the standard IUPAC, except gap is 15 or Hex F
     * @param halfByte
     * @return character array of length 1 to 3 based on halfByte
     */
    public static byte[] getSNPValueFromHalfByte(byte halfByte) {
        switch (halfByte) {
            case 0xE:
                return NIUPAC;
            case 0x0:
                return AIUPAC;
            case 0x1:
                return CIUPAC;
            case 0x2:
                return GIUPAC;
            case 0x3:
                return TIUPAC;
            case 0xF:
                return gapIUPAC;
            case 0x4:
                return RIUPAC;
            case 0x5:
                return YIUPAC;
            case 0x6:
                return SIUPAC;
            case 0x7:
                return WIUPAC;
            case 0x8:
                return KIUPAC;
            case 0x9:
                return MIUPAC;
            case 0xA:
                return BIUPAC;
            case 0xB:
                return DIUPAC;
            case 0xC:
                return HIUPAC;
            case 0xD:
                return VIUPAC;
            default:
                return NIUPAC;
        }
    }

    /**
     * Conversion of SNP Values to a half byte (4bit) encoding.
     * This uses the standard IUPAC, except gap is 15 or Hex F
     * possible TODO: currently assumes snpValues characters are capital letters
     * mistaking one character for two is not possible as standard ASCII values
     * go up to 126 and lowest two character summation in use is 132
     * @param snpValue
     * @return two bit encoding of base in a byte
     */
    public static byte getHalfByteFromSNPValue(byte[] snpValue) {
        int valueSum = 0;
        for (int i = 0; i < snpValue.length; i++) {
            valueSum += snpValue[i];
        }
        switch (valueSum) {
            case (int) GdpdmBLOBUtils.UNKNOWN_CHARACTER:
                return 0xE;
            case 65: //'A'
                return 0x0;
            case 67: //'C'
                return 0x1;
            case 71: //'G'
                return 0x2;
            case 84: //'T'
                return 0x3;
            case (int) GdpdmBLOBUtils.GAP_CHARACTER:
                return 0xF;
            case 90: //'-' + '-'
                return 0xF;
            case 136: //'A+G'
                return 0x4;
            case 151: //'C+T'
                return 0x5;
            case 138: //'G+C'
                return 0x6;
            case 149: //'A+T'
                return 0x7;
            case 155: //'G+T'
                return 0x8;
            case 132: //'A+C'
                return 0x9;
            case 43: //'+'
                return 0xA;
            case 86: //'+' + '+'
                return 0xA;
            case 48: //'0'
                return 0xB;
            case 96: //'0' + '0'
                return 0xB;
            case 88: //'+' + '-'
                return 0xB;
            case 130: //'A+A'
                return 0x0;
            case 134: //'C+C'
                return 0x1;
            case 142: //'G+G'
                return 0x2;
            case 168: //'T+T'
                return 0x3;
            default:
                return 0xE;
        }
    }

    /**
     * Conversion of half byte (4bit) encoding to a sequence character (2 bytes).
     * This uses the standard IUPAC, except gap is 15 or Hex F
     * @param halfByte
     * @return 1 byte sequence character
     */
    public static byte getBaseFromHalfByte(byte halfByte) {
        switch (halfByte) {
            case 0xE:
                return GdpdmBLOBUtils.UNKNOWN_CHARACTER;
            case 0x0:
                return 'A';
            case 0x1:
                return 'C';
            case 0x2:
                return 'G';
            case 0x3:
                return 'T';
            case 0xF:
                return GdpdmBLOBUtils.GAP_CHARACTER;
            case 0x4:
                return 'R';
            case 0x5:
                return 'Y';
            case 0x6:
                return 'S';
            case 0x7:
                return 'W';
            case 0x8:
                return 'K';
            case 0x9:
                return 'M';
            case 0xA:
                return '+';
            case 0xB:
                return '0';
            case 0xC:
                return GdpdmBLOBUtils.UNKNOWN_CHARACTER;  //unused
            case 0xD:
                return GdpdmBLOBUtils.UNKNOWN_CHARACTER;  //unused
            default:
                return GdpdmBLOBUtils.UNKNOWN_CHARACTER;
        }
    }

    /**
     * Conversion of sequence character (2 bytes) to a half byte (4bit) encoding.
     * This uses the standard IUPAC, except gap is 15 or Hex F
     * @param base
     * @return two bit encoding of base in a byte
     */
    public static byte getHalfByteFromBase(char base) {
        switch (base) {
            case GdpdmBLOBUtils.UNKNOWN_CHARACTER:
                return 0xE;
            case 'A':
                return 0x0;
            case 'C':
                return 0x1;
            case 'G':
                return 0x2;
            case 'T':
                return 0x3;
            case GdpdmBLOBUtils.GAP_CHARACTER:
                return 0xF;  //these are out of order because they are more common then below
            case 'R':
                return 0x4;
            case 'Y':
                return 0x5;
            case 'S':
                return 0x6;
            case 'W':
                return 0x7;
            case 'K':
                return 0x8;
            case 'M':
                return 0x9;
            case '+':
                return 0xA;
            case '0':
                return 0xB;
            case 'H':
                return 0xC;
            case 'V':
                return 0xD;
            default:
                return 0xE; //N
        }
    }

    public static final boolean isBaseHomozygous(byte base) {
        switch (base) {
            case 'A':
                return true;
            case 'C':
                return true;
            case 'G':
                return true;
            case 'T':
                return true;
            case GdpdmBLOBUtils.GAP_CHARACTER:
                return true;  //these are out of order because they are more common then below
            case '+':
                return true;
            default:
                return false; //N
        }
    }

    public static void setHalfByteInAlleleBLOB(byte[] seqBA, int site, char base) {
        int bbase = getHalfByteFromBase(base);
        int index = (site / 2) + GdpdmBLOBUtils.totalHeaderPadding;
        if (site % 2 == 0) {
            seqBA[index] = (byte) ((bbase << 4) | (seqBA[index] & 0x0f));
        } else {
            seqBA[index] = (byte) (bbase | (seqBA[index] & 0xf0));
        }
    }

    /**
     *
     * @param args
     * @throws java.lang.Exception
     */
    public static void main(String[] args) throws Exception {
        String pedName = "flapjack multiple test.flpjk.geno";
        String mapName = "flapjack multiple test.flpjk.map";
//        String name = "hapmap multiple test.txt";

//        net.maizegenetics.pal.alignment.Alignment[] a = net.maizegenetics.pal.alignment.ImportUtils.readFromFlapjack(pedName, mapName);
//
//        net.maizegenetics.pal.alignment.ExportUtils.writeToFlapjack(a, "outtest", '\t');
//
//        a.getDataType();
//        a.getIdGroup();

//        net.maizegenetics.pal.alignment.ExportUtils.writeToHapmap(a, true, name, '\t');

//        String[][] chroms = GdpdmBLOBUtils.getFileInfo(name, 2, 0);
//
//        System.out.println(chroms[0][0]);
//        System.out.println(GdpdmBLOBUtils.getChromInFileInfo(chroms, "chrom")[1]);
//        System.out.println(GdpdmBLOBUtils.countLinesInFile(chroms));
//        System.out.println(GdpdmBLOBUtils.countLinesInFile(name));

//        net.maizegenetics.pal.alignment.Alignment a = net.maizegenetics.pal.alignment.ImportUtils.readFromZip(name);
//        int chr = 8;
////           String hapmapInputFile="C:/EdStuff/Solexa/test/"+chr+".snps.haploid.1.13.log2ml.HP1.txt", testDirectory="C:/EdStuff/Solexa/test/";
//        String hapmapInputFile = "E:/SolexaAnal/AGPLog2ML/combined1_13HP/" + chr + ".snps.combined.1.13.log2ml.HP1.txt", testDirectory = "E:/SolexaAnal/test/";
////           ArrayList<byte[]> theBLOBs=createAllelePositionBLOBs(hapmapInputFile);
////           GdpdmBLOBUtils.writeBLOBtoFiles(theBLOBs, testDirectory+"sepf/");
////           GdpdmBLOBUtils.writeBLOBtoZip(theBLOBs, testDirectory + "all_" +chr+ "pos.zip");
////           GdpdmBLOBUtils.writeBLOBtoLZMA(theBLOBs, testDirectory + "all_" +chr+ "pos.7z");
//        //       Pack1Alignment p1a=(Pack1Alignment)GdpdmBLOBUtils.createPack1AlignmentFromFile("E:/SolexaAnal/test/all_"+chr+"pos.7z", "", "");
//        Pack1Alignment p1a = (Pack1Alignment) ImportUtils.createPack1AlignmentFromFile("C:/EdStuff/Solexa/test/all_8pos.7z", "", "");
//        //        Pack1Alignment p1a=(Pack1Alignment)GdpdmBLOBUtils.createPack1AlignmentFromFile("C:/EdStuff/Solexa/test/all_8pos.7z", "C:/EdStuff/Solexa/test/all_8cov.7z", "");
//        //   Pack1Alignment p1a=(Pack1Alignment)GdpdmBLOBUtils.createPack1AlignmentFromFile("C:/EdStuff/Solexa/test/all_8pos.7z", "C:/EdStuff/Solexa/test/all_8cov.7z", "");
//
//        int closestSNP = Math.abs(p1a.getSiteOfPhysicalPosition(179738793, null));
//        System.out.println("Site closest to position 179738793 is: " + closestSNP);
//        for (int j = closestSNP - 100; j < closestSNP + 100; j++) {
//            if (p1a.getSiteSummary(j).getNumberMissing() > 8) {
//                continue;
//            }
//            System.out.print(j + "\t" + p1a.getPositionInLocus(j) + "\t");
//            for (int i = 0; i < p1a.getSequenceCount(); i++) {
//                System.out.print(p1a.getBaseChar(i, j) + "\t");
//            }
//            System.out.println();
//        }
//        for (int j = 0; j < 1; j++) {
//            for (int i = 0; i < p1a.getSequenceCount(); i++) {
//                // p1a.getAlignedSequenceString(i);
//                // System.out.println(p1a.getIdentifier(i).getName()+"  "+ p1a.getAlignedSequenceString(i).substring(0, 10));
//                // System.out.println(p1a.getIdentifier(i).getName()+"  "+ p1a.getAlignedSequenceString(i).substring(p1a.getSiteCount()-10));
//                byte[] allele = new byte[p1a.getSiteCount()];
//                // for(int b=0; b<p1a.getSiteCount(); b++) {allele[b]=p1a.getData( i, b);}
//
//                //
//                allele = p1a.getBaseRange(i, 0, p1a.getSiteCount() - 1);
//
//                System.out.print(p1a.getIdGroup().getIdentifier(i).getName() + "\t");
//                System.out.print(new String(allele, closestSNP - 100, 100));
//                System.out.print("  ");
//                System.out.println(new String(allele, closestSNP, 100));
//            }
//        }
//        System.out.println(p1a.getPositionInLocus(100));
//        System.out.println(new String(p1a.getBaseRange(0, 100, 150)));
//        System.out.println(p1a.getPositionInLocus(150));
    }
}
