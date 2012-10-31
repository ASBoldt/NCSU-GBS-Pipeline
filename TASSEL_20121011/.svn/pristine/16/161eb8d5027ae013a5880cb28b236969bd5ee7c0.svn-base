/*
 * GdpdmBLOBUtils
 */
package net.maizegenetics.pal.alignment;

import SevenZip.Compression.LZMA.Decoder;
import SevenZip.Compression.LZMA.Encoder;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;

import java.net.MalformedURLException;
import java.net.URL;

import java.nio.ByteBuffer;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;

import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;
import java.util.zip.ZipOutputStream;

import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import net.maizegenetics.pal.datatype.DataType;

import java.util.regex.Pattern;
import net.maizegenetics.util.Utils;

/**
 *
 * @author ed
 */
public class GdpdmBLOBUtils {

    public static final char UNKNOWN_CHARACTER = DataType.UNKNOWN_CHARACTER; //'N';
    public static final char GAP_CHARACTER = DataType.GAP_CHARACTER; //'-';
    public static byte[] bases = {'A', 'C', 'G', 'T', 'R', 'Y', 'S', 'W', 'K', 'M', '+', '0', 'H', 'V', UNKNOWN_CHARACTER, GAP_CHARACTER};  //the question
    /**
     * Constants that describe the GDPDM blobs used for storing diversity data
     * Field variables define the start and end position in the blob for the variables
     */
    public static final int totalHeaderPadding = 1024;
    public static final byte[] compressionAlgorithm = {'0', '0', '1'};
    public static final int[] compressionAlgorithmField = {0, 3}; //start position in header, maxLength
    public static final int[] blobTypeField = {3, 1}; //start position in header, maxLength
    public static final int[] numSitesField = {4, 4}; //start position in header, maxLength
    public static final int[] genomeBuildField = {8, 10};  //start position in header, maxLength
    public static final int[] locusField = {18, 25}; //start position in header, maxLength
    public static final int[] startField = {43, 4};  //start position in header, maxLength
    public static final int[] endField = {47, 4};  //start position in header, maxLength
    public static final int[] taxaField = {51, 150}; //start position in header, maxLength
    public static final int idLengthField = 201; //start position in header, maxLength determined dynamically
    public static final String[] blobTypeSuffixForZip = {".bc01", ".bc02", ".bc03", ".bc04", ".bc05"};
    public static final String[] blobTypeNameForZip = {"AlleleValueBLOB", "PositionBLOB", "CoverageBLOB", "RefSeqBLOBs", "IdBLOB", "IndelSizeBLOB", "SiteQualityBLOB"};
    public static final byte alleleBLOBtype = '1'; //div_allele.binary_value
    public static final byte allelePositionBLOBtype = '2'; //div_allele_assay.binary_position
    public static final byte coverageBLOBtype = '3';
    public static final byte refSeqBLOBtype = '4'; //same structure as alleleBLOBtype except contiguous bases
    public static final byte SNPIdBLOBtype = '5'; //div_allele_assay.binary_id
    public static final byte indelSizeBLOBtype = '6';  //div_allele_assay.binary_annotation
    public static final byte siteQualityBLOBtype = '7';  //div_allele_assay.binary_annotation
    public static final int totalBLOBtypes = 7;

    // Compiled Whitespace Pattern
    private static Pattern WHITESPACE_PATTERN = Pattern.compile("\\s");
    private static int bufferSize = 67108864;

    /**
     * Method for setting the headers for the blob.  Taxon is only need for certain types.
     * @param blob
     * @param compressionAlgorithm
     * @param blobType
     * @param numSites
     * @param genomeBuild
     * @param locus
     * @param minPosition
     * @param maxPosition
     * @param taxonName
     */
    public static void setHeaderOfBLOB(byte[] blob, byte[] compressionAlgorithm, byte blobType, int numSites, String genomeBuild, String locus, int minPosition, int maxPosition, String taxonName) {
        ByteBuffer bbb = ByteBuffer.wrap(blob);
        bbb.position(compressionAlgorithmField[0]);
        bbb.put(compressionAlgorithm);
        bbb.position(blobTypeField[0]);
        bbb.put(blobType);
        bbb.position(numSitesField[0]);
        bbb.putInt(numSites);
        bbb.position(genomeBuildField[0]);
        // Byte 9-18: Bit 65-144: Genome Version (10 byte) (ASCII)
        bbb.put(genomeBuild.getBytes(), 0, Math.min(genomeBuild.length(), genomeBuildField[1]));
        bbb.position(locusField[0]);
        //Byte 19-43: Bit 145-344: Chromosome Name (25 byte) (ASCII)
        bbb.put(locus.getBytes(), 0, Math.min(locus.length(), locusField[1]));
        bbb.position(startField[0]);
        //Byte 44-47: Bit 345-376: Start Position (4 bytes). Smallest Absolute Physical Position.
        bbb.putInt(minPosition);
        bbb.position(endField[0]);
        //Byte 48-51: Bit 377-408: End Position (4 bytes). Largest Absolute Physical Position
        bbb.putInt(maxPosition);
        if (taxonName != null) {
            bbb.position(taxaField[0]);
            //Byte 52-201: Accession Name (25 bytes) (ASCII)
            bbb.put(taxonName.getBytes(), 0, Math.min(taxonName.length(), taxaField[1]));
        }
    }

    /**
     * Method for setting the headers for the blob.  Taxon is only need for certain types.
     * @param blob
     * @param compressionAlgorithm
     * @param blobType
     * @param numSites
     * @param genomeBuild
     * @param locus
     * @param minPosition
     * @param maxPosition
     * @param taxonName
     * @param idLength
     */
    public static void setHeaderOfBLOB(byte[] blob, byte[] compressionAlgorithm, byte blobType, int numSites, String genomeBuild, String locus, int minPosition, int maxPosition, String taxonName, int idLength) {
        ByteBuffer bbb = ByteBuffer.wrap(blob);
        bbb.position(compressionAlgorithmField[0]);
        bbb.put(compressionAlgorithm);
        bbb.position(blobTypeField[0]);
        bbb.put(blobType);
        bbb.position(numSitesField[0]);
        bbb.putInt(numSites);
        bbb.position(genomeBuildField[0]);
        // Byte 9-18: Bit 65-144: Genome Version (10 byte) (ASCII)
        bbb.put(genomeBuild.getBytes(), 0, Math.min(genomeBuild.length(), genomeBuildField[1]));
        bbb.position(locusField[0]);
        //Byte 19-43: Bit 145-344: Chromosome Name (25 byte) (ASCII)
        bbb.put(locus.getBytes(), 0, Math.min(locus.length(), locusField[1]));
        bbb.position(startField[0]);
        //Byte 44-47: Bit 345-376: Start Position (4 bytes). Smallest Absolute Physical Position.
        bbb.putInt(minPosition);
        bbb.position(endField[0]);
        //Byte 48-51: Bit 377-408: End Position (4 bytes). Largest Absolute Physical Position
        bbb.putInt(maxPosition);
        if (taxonName != null) {
            bbb.position(taxaField[0]);
            //Byte 52-201: Accession Name (25 bytes) (ASCII)
            bbb.put(taxonName.getBytes(), 0, Math.min(taxonName.length(), taxaField[1]));
        }
        bbb.position(idLengthField);
        //Byte 202-205: SNP id (4 bytes)
        bbb.putInt(idLength);
    }

    /**
     * Count the number of lines in a file
     * @param infileName
     * @return lines
     */
    public static int countLinesInFile(String infileName) {
        int lines = 0;
        try {
            //BufferedReader fileIn = new BufferedReader(new FileReader(infileName), 1000000);
            BufferedReader fileIn = getBufferedReader(infileName, 1000000);
            while (fileIn.readLine() != null) {
                lines++;
            }
            fileIn.close();
        } catch (Exception e) {
            System.err.println("File IO in countLinesInFile: " + e);
        }
        return lines;
    }

    public static int countLinesInFile(String[][] chromCounts) {
        int lines = 0;
        for (int i = 0; i < chromCounts.length; i++) {
            lines += Integer.parseInt(chromCounts[i][1]);
        }
        return lines;
    }

    /**
     * Creates appropriate BufferedReader depending on data loaction
     * @param inSourceName
     * @return BufferedReader
     */
    public static BufferedReader getBufferedReader(String inSourceName, int bufSize) {
        try {
            if (inSourceName.startsWith("http")) {
                return new BufferedReader(new InputStreamReader((new URL(inSourceName)).openStream()), bufSize);
            } else {
                return new BufferedReader(new FileReader(inSourceName), bufSize);
            }
        } catch (Exception e) {
            System.err.println("File IO in getBufferedReader: " + e);
        }
        return null;
    }

//    /**
//     * Finds the length of the longest SNP ID in the given file
//     * @param infileName
//     * @param depthInFile
//     * @param sites
//     * @param type
//     * Hapmap type = 0
//     * Plink type = 1
//     * Flapjack type = 0
//     * @param hasHeader
//     * @return SNPIDLength
//     */
//    public static int getSNPIDLength(String infileName, int depthInFile, int numSites, int type, boolean hasHeader) {
//        int SNPIDLength = 0;
//        try {
//            //BufferedReader fileIn = new BufferedReader(new FileReader(infileName), 1000000);
//            BufferedReader fileIn = getBufferedReader(infileName, 1000000);
//            //BufferedReader fileIn = new BufferedReader(new InputStreamReader((new URL(infileName)).openStream()), 1000000);
//            if(hasHeader) {
//                fileIn.readLine();
//            }
//            for (int depth = 0; depth < depthInFile - 1; depth++) {
//                fileIn.readLine(); //shift file to desired position
//            }
//            for (int i = 0; i < numSites; i++) {
//                String SNPID = WHITESPACE_PATTERN.split(fileIn.readLine())[type];
//                if (SNPID.length() > SNPIDLength) {
//                    SNPIDLength = SNPID.length();
//                }
//            }
//            fileIn.close();
//        } catch (Exception e) {
//            System.err.println("File IO in countLinesInFile: " + e);
//        }
//        return SNPIDLength;
//    }
    /**
     * Get chromosome IDs for all represented chromosomes in a given .map file
     * along with number of sites on each chromosome
     * currently unique to Plink files
     * @param infileName
     * @param chromFileType, column number where chromosome ID are stored, varies based on format
     * Hapmap type = 2
     * Plink type = 0
     * FlapJack type = 1
     * @param idFileType, column number where SNP ID are stored, varies based on format
     * Hapmap type = 0
     * Plink type = 1
     * Flapjack type = 0
     * @return stringArray
     */
    public static String[][] getFileInfo(String infileName, int chromFileType, int idFileType) {
        try {
            ArrayList<String> stringList = new ArrayList<String>();
            ArrayList<Integer> countList = new ArrayList<Integer>();
            ArrayList<Integer> idLengthList = new ArrayList<Integer>();
            int count = 0;
            String previousChrom = "";
            int SNPidLength = 0;
            //BufferedReader fileIn = new BufferedReader(new FileReader(infileName), 1000000);
            //BufferedReader fileIn = getBufferedReader(infileName, 1000000);
            BufferedReader fileIn = Utils.getBufferedReader(infileName, 1000000);
            //BufferedReader fileIn = new BufferedReader(new InputStreamReader((new URL(infileName)).openStream()), 1000000);
            String fileInLine = fileIn.readLine();
            while (fileInLine != null) {
                String[] fileLine = WHITESPACE_PATTERN.split(fileInLine);
                String currentChrom = fileLine[chromFileType].trim();
                String SNPID = fileLine[idFileType].trim();
                if (!currentChrom.equals(previousChrom)) {
                    if (stringList.size() != 0) {
                        countList.add(new Integer(count));
                        idLengthList.add(new Integer(SNPidLength));
                        count = 0;
                        SNPidLength = 0;
                    }
                    stringList.add(currentChrom);
                }
                count++;
                if (SNPID.length() > SNPidLength) {
                    SNPidLength = SNPID.length();
                }
                previousChrom = currentChrom;
                fileInLine = fileIn.readLine();
            }
            countList.add(new Integer(count));
            idLengthList.add(new Integer(SNPidLength));
            String[][] stringArray = new String[stringList.size()][3];
            for (int i = 0; i < stringList.size(); i++) {
                stringArray[i][0] = stringList.get(i);
                //System.out.println(stringList.get(i));
                stringArray[i][1] = countList.get(i).toString();
                //System.out.println(countList.get(i));
                stringArray[i][2] = idLengthList.get(i).toString();
            //System.out.println(stringArray[i][2]);
            }
            return stringArray;
        } catch (Exception e) {
            System.err.println("File IO in getFileInfo: " + e);
        }
        return null;
    }

    /**
     * Returns chromosome names, sorted positions, id lengths for each chromosome
     * for Flapjack files
     * @param infileName
     * @return
     */
    public static Object[][] getFlapjackInfo(String infileName) {
        try {
            ArrayList<String> chromNameList = new ArrayList<String>();
            ArrayList<ArrayList<Integer>> chromPositionList = new ArrayList<ArrayList<Integer>>();
            ArrayList<Integer> idLengthList = new ArrayList<Integer>();

            String previousChrom = "";
            BufferedReader fileIn = getBufferedReader(infileName, 1000000);

            String fileInLine;
            while ((fileInLine = fileIn.readLine()) != null) {
                String[] fileLine = WHITESPACE_PATTERN.split(fileInLine);
                //check chromosome
                int i = chromNameList.lastIndexOf(fileLine[1]);
                if (i == -1) {
                    chromNameList.add(fileLine[1]);
                    ArrayList<Integer> tempPosition = new ArrayList<Integer>();
                    tempPosition.add(new Integer(fileLine[2]));
                    chromPositionList.add(tempPosition);
                    idLengthList.add(fileLine[0].length());
                }
                else {
                    chromPositionList.get(i).add(new Integer(fileLine[2]));
                    if (fileLine[0].length() > idLengthList.get(i)) {
                        idLengthList.set(i, fileLine[0].length());
                    }
                }
            }
            return sortChroms(chromNameList, chromPositionList, idLengthList);
        }
        catch (Exception e) {
            System.err.println("File IO in getFlapjackinfo: " + e);
            return null;
        }
    }

    /**
     * Sorts by chromosome if chromosomes are numeric and sorts within chromosome by position
     * returns sorted arrays
     * @param nameList
     * @param positionList
     * @param idLengthList
     * @return sortedArray
     */
    public static Object[][] sortChroms(ArrayList<String> nameList, ArrayList<ArrayList<Integer>> positionList, ArrayList<Integer> idLengthList) {
        try {
            int[] chroms = new int[nameList.size()];
            int[] sortedChroms = new int[nameList.size()];
            for (int i = 0; i < nameList.size(); i++) {
                chroms[i] = Integer.parseInt(nameList.get(i));
                sortedChroms[i] = Integer.parseInt(nameList.get(i));
            }
            // sort chromosomes
            Arrays.sort(sortedChroms);
            // sort position and idLength lists accordingly
            for (int i = 0; i < chroms.length - 1; i++) {
                for (int j = 0; j < sortedChroms.length; j++) {
                    if (chroms[i] == sortedChroms[j]) {
                        ArrayList<Integer> temp = positionList.get(j);
                        positionList.set(j, positionList.get(i));
                        positionList.set(i, temp);
                        Integer tempInt = idLengthList.get(j);
                        idLengthList.set(j, idLengthList.get(i));
                        idLengthList.set(i, tempInt);
                        j = sortedChroms.length;
                    }
                }
            }
            // populate return array
            Object[][] sortedArray = new Object[sortedChroms.length][3];
            Hashtable[] tableArray = sortPositions(positionList);
            for (int i = 0; i < sortedChroms.length; i++) {
                sortedArray[i][0] = "" + sortedChroms[i];
                sortedArray[i][1] = tableArray[i];
                sortedArray[i][2] = idLengthList.get(i);
            }
            return sortedArray;
        }
        // if chromosomes are not numeric, only sort by position
        catch (NumberFormatException nfe) {
            Object[][] sortedArray = new Object[nameList.size()][3];
            Hashtable[] tableArray = sortPositions(positionList);
            for (int i = 0; i < nameList.size(); i++) {
                sortedArray[i][0] = nameList.get(i);
                sortedArray[i][1] = tableArray[i];
                sortedArray[i][2] = idLengthList.get(i);
            }
            return sortedArray;
        }
    }

    /**
     * Sorts by position and stores data in hashtables, with position being the key and index being the value
     * @param positionList
     * @return tableArray
     */
    public static Hashtable[] sortPositions(ArrayList<ArrayList<Integer>> positionList) {
        Hashtable[] tableArray = new Hashtable[positionList.size()];
        for (int i = 0; i < positionList.size(); i++) {
            tableArray[i] = new Hashtable(positionList.get(i).size());
            int[] positions = new int[positionList.get(i).size()];
            for (int j = 0; j < positions.length; j++) {
                positions[j] = positionList.get(i).get(j);
            }
            Arrays.sort(positions);
            for (int j = 0; j < positions.length; j++) {
                tableArray[i].put(positions[j], j);
            }
        }
        return tableArray;
    }

    /**
     * Returns the number of sites on a given chromosome as well as how far into a .map file that the
     * sites for a given chromosome begins, works in conjunction with getChromsAvailableCounts
     * @param chromCounts
     * @param chrom
     * @return chromInfo[0] = numSites, chromInfo[1] = depth in .map file the sites begin
     * chromInfo[2] = length of longest SNP id within the chromosome
     */
    public static int[] getChromInFileInfo(String[][] chromCounts, String chrom) {
        int[] chromInfo = {0, 0, 0};
        for (int i = 0; i < chromCounts.length; i++) {
            if (chromCounts[i][0].equals(chrom)) {
                chromInfo[0] = Integer.parseInt(chromCounts[i][1]);
                chromInfo[2] = Integer.parseInt(chromCounts[i][2]);
                break;
            }
            chromInfo[1] += Integer.parseInt(chromCounts[i][1]);
        }
        return chromInfo;
    }

    /**
     * Reads alleleBLOB and siteBLOB from LZMA compressed file
     * @param infileName the name of the compressed file
     * @param theBLOBList blobs in the order allele[][], site[], reference[], and coverage[][]
     * @return ArrayList<byte[]>[] the same blob list as above but with the added blobs
     */
    public static ArrayList<byte[]>[] readBLOBfromLZMA(String infileName, ArrayList<byte[]>[] theBLOBList) {
        try {
            FileInputStream fis = new FileInputStream(infileName);
            int propertiesSize = 5;
            byte[] properties = new byte[propertiesSize];
            if (fis.read(properties, 0, propertiesSize) != propertiesSize) {
                throw new IOException("input .lzma file is too short");
            }
            Decoder decoder = new Decoder();
            if (!decoder.SetDecoderProperties(properties)) {
                throw new IOException("Incorrect stream properties");
            }
            long outSize = 0;
            for (int i = 0; i < 8; i++) {
                int v = fis.read();
                if (v < 0) {
                    throw new IOException("Can\'t read stream size");
                }
                outSize |= ((long) v) << (8 * i);
            }
            ByteArrayOutputStream baos = new ByteArrayOutputStream((int) outSize + 10);
            if (!decoder.Code(fis, baos, outSize)) {
                throw new IOException("Error in data stream");
            }
            fis.close();
            ByteBuffer bb = ByteBuffer.wrap(baos.toByteArray());
            baos.close();
            int numBLOBs = bb.getInt();
            for (int i = 0; i < numBLOBs; i++) {
                int entrySize = bb.getInt();
                byte[] theBLOB = new byte[entrySize];
                bb.get(theBLOB);
                int blobType = Integer.parseInt("" + (char) theBLOB[GdpdmBLOBUtils.blobTypeField[0]]);
                theBLOBList[blobType - 1].add(theBLOB);
            }
            System.out.println("Compressed LZMA read");
        } catch (IOException ee) {
            System.err.println("Data could not be saved: " + ee);
        }
        return theBLOBList;
    }

    public static void writeBLOBtoFiles(ArrayList<byte[]> theBLOBList, String outFileDirectory) {
        try {
            for (byte[] b : theBLOBList) {
                String s = new String(b, taxaField[0], taxaField[1]).trim();
                if (s.equals("")) {
                    s = "all";
                }
                //this is a likely a position file
                String[] taxaParts = s.split(":");
                String locusName = (new String(b, locusField[0], locusField[1])).trim();
                String suffix = ".unk";
                int blobType = Integer.parseInt("" + (char) b[GdpdmBLOBUtils.blobTypeField[0]]);
                suffix = blobTypeSuffixForZip[blobType - 1];
                FileOutputStream fos = new FileOutputStream(outFileDirectory + taxaParts[0] + "_" + locusName + suffix);
                fos.write(b);
                fos.close();
            }
        } catch (IOException ee) {
            System.err.println("Data could not be saved: " + ee);
        }
    }

    /**
     * Writes alleleBLOB and siteBLOB to a LZMA compressed file
     * The headers for the file include file of the file, parameters of compression,
     * Format is:
     * Encoding properties
     * uncompressed fileSize
     * number of BLOBs
     * blobSize1 blob1
     * blobSize2 blob2
     * etc.
     * @param theBLOBList the list of BLOBs to be written
     * @param outfileName the name of the compressed file
     */
    public static void writeBLOBtoLZMA(ArrayList<byte[]> theBLOBList, String outfileName) {
        try {
            FileOutputStream fos = new FileOutputStream(outfileName);
            //  ZipOutputStream zos = new ZipOutputStream(fos);
            int totalSizeHeaders = theBLOBList.size() + 1;
            //the 1 is for the header telling how many blobs there are
            int size = 4 * totalSizeHeaders;
            for (byte[] b : theBLOBList) {
                size += b.length;
            }
            ByteBuffer bb = ByteBuffer.allocate(size);
            bb.putInt(theBLOBList.size());
            // number of taxon
            for (byte[] b : theBLOBList) {
                bb.putInt(b.length);
                bb.put(b);
            }
            ByteArrayInputStream bais = new ByteArrayInputStream(bb.array());
            Encoder encoder = new Encoder();
            encoder.SetAlgorithm(1);
            int dictionary = 1 << 21;
            encoder.SetDictionarySize(dictionary);
            encoder.WriteCoderProperties(fos);
            long fileSize = bb.position();
            for (int i = 0; i < 8; i++) {
                fos.write((int) (fileSize >>> (8 * i)) & 255);
            }
            encoder.Code(bais, fos, -1, -1, null);
            fos.flush();
            fos.close();
            bais.close();
            System.out.println("Compressed LZMA written");
        } catch (IOException ee) {
            System.err.println("Data could not be saved: " + ee);
        }
    }

    /**
     * Writes alleleBLOB, siteBLOB, and idBLOB to a gzip compressed file
     * @param theBLOBList the list of BLOBs to be written
     * @param outfileName the name of the compressed file
     */
    public static void writeBLOBtoGZIP(ArrayList<byte[]>[] theBLOBList, String outfileName) {
        try {
            FileOutputStream fos = new FileOutputStream(outfileName);
            GZIPOutputStream gos = new GZIPOutputStream(fos, bufferSize);
            for (int i = 0; i < theBLOBList.length; i++) {
                for (byte[] b : theBLOBList[i]) {
                    gos.write(b);
                }
            }
            gos.close();
        } catch (IOException ee) {
            System.err.println("Data could not be save: " + ee);
        }
    }

    /**
     * Reads alleleBLOB, siteBLOB, and idBLOB from a gzip compressed file
     * @param infileName name of the compressed file
     * @param theBLOBList blobs
     * @return ArrayList<byte[]>[] theBLOBList populated
     */
    public static ArrayList<byte[]>[] readBLOBfromGZIP(String infileName, ArrayList<byte[]>[] theBLOBList) throws MalformedURLException, IOException, FileNotFoundException {

        System.out.println("Reading:" + infileName);

        GZIPInputStream gis = null;

        if (infileName.startsWith("http")) {
            gis = new GZIPInputStream((new URL(infileName)).openStream(), bufferSize);
        } else {
            gis = new GZIPInputStream(new FileInputStream(infileName), bufferSize);
        }

        //ZipEntry ze;
        //fis = new FileInputStream(infileName);
        //zis = new ZipInputStream(fis);
        ByteBuffer bb = ByteBuffer.allocate(50000000);
        //byte[] buf = new byte[40960];
        byte[] bufHeader = new byte[1024];
        byte[] bufBody;
        byte[] bufTemp;
        int blobType = 0;
        int c = 0;
        while ((c = gis.read(bufHeader)) != -1) {
            bb.clear();
            bb.put(bufHeader, 0, c);
            if (c != 1024) {
                int bytesLeft = 1024 - c;
                while (bytesLeft > 0) {
                    bufTemp = new byte[bytesLeft];
                    if ((c = gis.read(bufTemp)) != -1) {
                        bb.put(bufTemp, 0, c);
                        bytesLeft -= c;
                    } else {
                        bytesLeft = 0;
                    }
                }
            }
            blobType = bb.get(GdpdmBLOBUtils.blobTypeField[0]) - 48;
            int numSites = bb.getInt(GdpdmBLOBUtils.numSitesField[0]);
            int bodyLength;
            if (blobType == 1) {
                int tempFieldLength = bb.getInt(GdpdmBLOBUtils.idLengthField);
                if (tempFieldLength != 4) {
                    bb.position(idLengthField);
                    bb.putInt(4);
                    bb.position(1024);
                }
                bodyLength = (numSites / 2) + (numSites % 2);
                bufBody = new byte[bodyLength];
                if ((c = gis.read(bufBody)) != -1) {
                    bb.put(bufBody, 0, c);
                    if (c != bodyLength) {
                        int bytesLeft = bodyLength - c;
                        while (bytesLeft > 0) {
                            bufTemp = new byte[bytesLeft];
                            if ((c = gis.read(bufTemp)) != -1) {
                                bb.put(bufTemp, 0, c);
                                bytesLeft -= c;
                            } else {
                                bytesLeft = 0;
                            }
                        }
                    }
                }
            } else if (blobType == 2) {
                int tempFieldLength = bb.getInt(GdpdmBLOBUtils.idLengthField);
                if (tempFieldLength != 32) {
                    bb.position(idLengthField);
                    bb.putInt(32);
                    bb.position(1024);
                }
                bodyLength = numSites * 4;
                bufBody = new byte[bodyLength];
                if ((c = gis.read(bufBody)) != -1) {
                    bb.put(bufBody, 0, c);
                    if (c != bodyLength) {
                        int bytesLeft = bodyLength - c;
                        while (bytesLeft > 0) {
                            bufTemp = new byte[bytesLeft];
                            if ((c = gis.read(bufTemp)) != -1) {
                                bb.put(bufTemp, 0, c);
                                bytesLeft -= c;
                            } else {
                                bytesLeft = 0;
                            }
                        }
                    }
                }
            } else if (blobType == 5) {
                int tempIdLength = bb.getInt(GdpdmBLOBUtils.idLengthField);
                int idLength = 0;
                if (tempIdLength < 24) {
                    idLength = tempIdLength;
                    bb.position(idLengthField);
                    bb.putInt(idLength * 8);
                    bb.position(1024);
                } else {
                    idLength = tempIdLength / 8;
                }
                bodyLength = numSites * idLength;
                bufBody = new byte[bodyLength];
                if ((c = gis.read(bufBody)) != -1) {
                    bb.put(bufBody, 0, c);
                    if (c != bodyLength) {
                        int bytesLeft = bodyLength - c;
                        while (bytesLeft > 0) {
                            bufTemp = new byte[bytesLeft];
                            if ((c = gis.read(bufTemp)) != -1) {
                                bb.put(bufTemp, 0, c);
                                bytesLeft -= c;
                            } else {
                                bytesLeft = 0;
                            }
                        }
                    }
                }
            }
            int entrySize = bb.position();
            byte[] theBLOB = new byte[entrySize];
            bb.position(0);
            bb.get(theBLOB);
            //int blobType = Integer.parseInt("" + (char) theBLOB[GdpdmBLOBUtils.blobTypeField[0]]);
            theBLOBList[blobType - 1].add(theBLOB);
            blobType = 0;
        }
        gis.close();
        System.out.println(infileName + " read");

        return theBLOBList;
    }

    public static void writeBLOBtoZip(ArrayList<byte[]>[] theBLOBList, String outfileName) {
        try {
            FileOutputStream fos = new FileOutputStream(outfileName);
            ZipOutputStream zos = new ZipOutputStream(fos);
            zos.setLevel(9);
            for (int i = 0; i < theBLOBList.length; i++) {
                if (theBLOBList[i].size() != 0) {
                    ZipEntry thisEntry = new ZipEntry(blobTypeNameForZip[i] + blobTypeSuffixForZip[i]);
                    zos.putNextEntry(thisEntry);
                }
                for (byte[] b : theBLOBList[i]) {
                    zos.write(b);
                }
                zos.closeEntry();
            }
            zos.close();
            System.out.println("Compressed zip written");
        } catch (IOException ee) {
            System.err.println("Data could not be saved: " + ee);
        }
    }

    /**
     * Create Pack1Alignment from the given infile.  The infile should be a zipped file for allele
     * blobs and one position blob
     * @param infileName
     * @param theBLOBList blobs in the order allele[][], site[], reference[], and coverage[][]
     * @return ArrayList<byte[]>[] the same blob list as above but with the added blobs
     */
    public static ArrayList<byte[]>[] readBLOBFromZip(String infileName, ArrayList<byte[]>[] theBLOBList) {
        try {
            System.out.println("Reading:" + infileName);

            ZipInputStream zis = null;

            if (infileName.startsWith("http")) {
                zis = new ZipInputStream((new URL(infileName)).openStream());
            } else {
                zis = new ZipInputStream(new FileInputStream(infileName));
            }

            //ZipEntry ze;
            //fis = new FileInputStream(infileName);
            //zis = new ZipInputStream(fis);
            ByteBuffer bb = ByteBuffer.allocate(50000000);
            //byte[] buf = new byte[40960];
            byte[] bufHeader = new byte[1024];
            byte[] bufBody;
            byte[] bufTemp;
            int blobType = 0;
            int c = 0;
            while (zis.getNextEntry() != null) {
                while ((c = zis.read(bufHeader)) != -1) {
                    bb.clear();
                    bb.put(bufHeader, 0, c);
                    if (c != 1024) {
                        int bytesLeft = 1024 - c;
                        while (bytesLeft > 0) {
                            bufTemp = new byte[bytesLeft];
                            if ((c = zis.read(bufTemp)) != -1) {
                                bb.put(bufTemp, 0, c);
                                bytesLeft -= c;
                            } else {
                                bytesLeft = 0;
                            }
                        }
                    }
                    blobType = bb.get(GdpdmBLOBUtils.blobTypeField[0]) - 48;
//                    System.out.println(blobType);
                    int numSites = bb.getInt(GdpdmBLOBUtils.numSitesField[0]);
                    int bodyLength;
                    if (blobType == 1) {
                        bb.position(idLengthField);
                        bb.putInt(4);
                        bb.position(1024);
                        bodyLength = (numSites / 2) + (numSites % 2);
                        bufBody = new byte[bodyLength];
                        if ((c = zis.read(bufBody)) != -1) {
                            bb.put(bufBody, 0, c);
                            if (c != bodyLength) {
                                int bytesLeft = bodyLength - c;
                                while (bytesLeft > 0) {
                                    bufTemp = new byte[bytesLeft];
                                    if ((c = zis.read(bufTemp)) != -1) {
                                        bb.put(bufTemp, 0, c);
                                        bytesLeft -= c;
                                    } else {
                                        bytesLeft = 0;
                                    }
                                }
                            }
                        }
                    } else if (blobType == 2) {
                        bb.position(idLengthField);
                        bb.putInt(32);
                        bb.position(1024);
                        bodyLength = numSites * 4;
                        bufBody = new byte[bodyLength];
                        if ((c = zis.read(bufBody)) != -1) {
                            bb.put(bufBody, 0, c);
                            if (c != bodyLength) {
                                int bytesLeft = bodyLength - c;
                                while (bytesLeft > 0) {
                                    bufTemp = new byte[bytesLeft];
                                    if ((c = zis.read(bufTemp)) != -1) {
                                        bb.put(bufTemp, 0, c);
                                        bytesLeft -= c;
                                    } else {
                                        bytesLeft = 0;
                                    }
                                }
                            }
                        }
                    } else if (blobType == 5) {
                        int tempIdLength = bb.getInt(GdpdmBLOBUtils.idLengthField);
                        int idLength = 0;
                        if (tempIdLength < 24) {
                            idLength = tempIdLength;
                            bb.position(idLengthField);
                            bb.putInt(idLength * 8);
                            bb.position(1024);
                        } else {
                            idLength = tempIdLength / 8;
                        }
                        bodyLength = numSites * idLength;
                        bufBody = new byte[bodyLength];
                        if ((c = zis.read(bufBody)) != -1) {
                            bb.put(bufBody, 0, c);
                            if (c != bodyLength) {
                                int bytesLeft = bodyLength - c;
                                while (bytesLeft > 0) {
                                    bufTemp = new byte[bytesLeft];
                                    if ((c = zis.read(bufTemp)) != -1) {
                                        bb.put(bufTemp, 0, c);
                                        bytesLeft -= c;
                                    } else {
                                        bytesLeft = 0;
                                    }
                                }
                            }
                        }
                    }
                    int entrySize = bb.position();
                    byte[] theBLOB = new byte[entrySize];
                    bb.position(0);
                    bb.get(theBLOB);
                    //int blobType = Integer.parseInt("" + (char) theBLOB[GdpdmBLOBUtils.blobTypeField[0]]);
                    theBLOBList[blobType - 1].add(theBLOB);
                }
                zis.closeEntry();
                blobType = 0;
            }
            zis.close();
            System.out.println(infileName + " read");
        } catch (Exception ee) {
            System.err.println("Data could not be read: " + ee);
        } catch (Error ee) {
            System.err.println("Data could not be read: " + ee);
        }
        return theBLOBList;
    }

//    public static void writePack1AlignmentToZip(Pack1Alignment p1a, String outfileName) {
//        ArrayList<byte[]> theBLOBList = new ArrayList<byte[]>();
//        for (int i = 0; i < p1a.getSequenceCount(); i++) {
//            theBLOBList.add(p1a.getAlleleBLOBs(i));
//        }
//        theBLOBList.add(p1a.getVariableSitesBLOB());
//        writeBLOBtoZip(theBLOBList, outfileName);
//    }
}
