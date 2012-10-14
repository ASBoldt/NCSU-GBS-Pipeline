/*
 * ExportUtils
 */
package net.maizegenetics.pal.alignment;

import java.io.FileWriter;
import java.io.BufferedWriter;

import java.util.ArrayList;

import net.maizegenetics.util.Utils;

import java.util.regex.Pattern;

import org.apache.log4j.Logger;




/**
 * The class exports PAL alignment datatypes to
 * various file formats.
 *
 * @author Jon
 */
public class ExportUtils {

    private static final Logger myLogger = Logger.getLogger(ExportUtils.class);

    private ExportUtils() {
        // Utility Class - do not instantiate.
    }

    /**
     * Writes multiple alignments to single Hapmap file. Currently no error checking
     * @param alignment array of alignemnts
     * @param diploid
     * @param filename
     * @param delimChar
     */
    public static void writeToHapmap(Alignment alignment, boolean diploid, String filename, char delimChar) {
        if (delimChar != ' ' && delimChar != '\t') {
            throw new IllegalArgumentException("Delimiter charater must be either a blank space or a tab.");
        }
//        for (int i = 0; i < alignment.length; i++) {
//            for (int j = i+1; j < alignment.length; j++) {
//                if (alignment[j].getSequenceCount() != alignment[i].getSequenceCount()) {
//                    throw new IllegalArgumentException("Selected chromosomes must have data for the same number of plant lines.");
//                }
//                for (int taxa = 0; taxa < alignment[i].getSequenceCount(); taxa++) {
//                    if (!alignment[i].getIdGroup().getIdentifier(taxa).getFullName().trim().equals(alignment[j].getIdGroup().getIdentifier(taxa).getFullName().trim())) {
//                        throw new IllegalArgumentException("Selected chromosomes must have data for the same plant lines.");
//                    }
//                }
//            }
//        }
        try {
            String fullFileName = Utils.addSuffixIfNeeded(filename, ".hmp.txt", new String[]{".hmp.txt", ".hmp.txt.gz"});
            BufferedWriter bw = Utils.getBufferedWriter(fullFileName);
            bw.write("rs#");
            bw.write(delimChar);
            bw.write("alleles");
            bw.write(delimChar);
            bw.write("chrom");
            bw.write(delimChar);
            bw.write("pos");
            bw.write(delimChar);
            bw.write("strand");
            bw.write(delimChar);
            bw.write("assembly#");
            bw.write(delimChar);
            bw.write("center");
            bw.write(delimChar);
            bw.write("protLSID");
            bw.write(delimChar);
            bw.write("assayLSID");
            bw.write(delimChar);
            bw.write("panelLSID");
            bw.write(delimChar);
            bw.write("QCcode");
            bw.write(delimChar);
//            int numTaxa = alignment[0].getSequenceCount();
            int numTaxa = alignment.getSequenceCount();
            for (int taxa = 0; taxa < numTaxa; taxa++) {//finish filling out first row
                //not completely sure this does what I want, I need to access the
                //accession name from every alleleBLOB in bytes [52-201] but there
                //doesn't seem to be a method to access that in Alignment
//                String sequenceID = alignment[0].getIdGroup().getIdentifier(taxa).getFullName().trim();
                String sequenceID = alignment.getIdGroup().getIdentifier(taxa).getFullName().trim();
                sequenceID = sequenceID.replaceAll("\\s", "_");
                bw.write(sequenceID);
                if (taxa != numTaxa - 1) {
                    bw.write(delimChar);
                }
            }
            bw.write("\n");
//            for (int n = 0; n < alignment.length; n++) {
//            int numSites = alignment[n].getSiteCount();
            int numSites = alignment.getSiteCount();
            for (int site = 0; site < numSites; site++) {
//                StringBuilder sb = new StringBuilder(alignment[n].getSNPID(site)); //rs#
                StringBuilder sb = new StringBuilder(alignment.getSNPID(site));
//                    if (alignment[n].getSNPID(site).length() == 0) {
//                        sb.append("chr");
//                        sb.append(alignment[n].getLocusName(site));
//                        sb.append("_");
//                        sb.append(alignment[n].getPositionInLocus(site));
//                    }
                sb.append(delimChar);
//                SiteSummary ss = alignment[n].getSiteSummary(site);
                SiteSummary ss = alignment.getSiteSummary(site);
                int numAlleles = ss.getAlleles().length;
                //currently does not correctly display if numAlleles > 2
                if (numAlleles == 0) {
                    sb.append("NA"); //if data does not exist
                }
                for (int i = 0; i < Math.min(2, numAlleles) ; i++) {
                    if (i > 0) {
                        sb.append('/');
                    }
                    sb.append(ss.getAlleles()[i]); //alleles
                }
                sb.append(delimChar);
//                sb.append(alignment[n].getLocusName(site)); //chrom
                sb.append(alignment.getLocusName(site));
                sb.append(delimChar);
//                sb.append(alignment[n].getPositionInLocus(site)); // pos
                sb.append(alignment.getPositionInLocus(site));
                sb.append(delimChar);
                sb.append("+"); //strand
                sb.append(delimChar);
                sb.append("NA"); //assembly# not supported
                sb.append(delimChar);
                sb.append("NA"); //center unavailable
                sb.append(delimChar);
                sb.append("NA"); //protLSID unavailable
                sb.append(delimChar);
                sb.append("NA"); //assayLSID unavailable
                sb.append(delimChar);
                sb.append("NA"); //panelLSID unavailable
                sb.append(delimChar);
                sb.append("NA"); //QCcode unavailable
                sb.append(delimChar);
                for (int taxa = 0; taxa < numTaxa; taxa++) {
                    if (diploid == false) {
//                        byte baseIUPAC = alignment[n].getBase(taxa, site);
                        byte baseIUPAC = alignment.getBase(taxa, site);
                        sb.append((char)baseIUPAC);
                    }else {
//                            byte[] b = AllelePositionBLOBUtils.getSNPValueFromBase((char)alignment[n].getBase(taxa, site));
                            byte[] b = AllelePositionBLOBUtils.getSNPValueFromBase((char)alignment.getBase(taxa, site));
                            if(b.length==1) {
                                sb.append((char)b[0]);
                                sb.append((char)b[0]);
                            }
                            else {
                                sb.append((char)b[0]);
                                sb.append((char)b[1]);
                            }
                    }
                    if (taxa != (numTaxa - 1)) {
                        sb.append(delimChar);
                    }
                }
                bw.write(sb.toString());
                bw.write("\n");
            }
//            }
            bw.close();
        }
        catch (Exception e) {
            System.out.println("Error writing writeToHapmap: "+e);
            e.printStackTrace();
        }
    }

    /**
     * Writes multiple alignments to single Hapmap file. Currently no error checking
     * @param alignment array of alignemnts
     * @param diploid
     * @param filename
     * @param delimChar
     */
    public static void writeToHapmap(Alignment alignment, AlignmentMask mask, boolean diploid, String filename, char delimChar) {
        if (delimChar != ' ' && delimChar != '\t') {
            throw new IllegalArgumentException("Delimiter charater must be either a blank space or a tab.");
        }
//        for (int i = 0; i < alignment.length; i++) {
//            for (int j = i+1; j < alignment.length; j++) {
//                if (alignment[j].getSequenceCount() != alignment[i].getSequenceCount()) {
//                    throw new IllegalArgumentException("Selected chromosomes must have data for the same number of plant lines.");
//                }
//                for (int taxa = 0; taxa < alignment[i].getSequenceCount(); taxa++) {
//                    if (!alignment[i].getIdGroup().getIdentifier(taxa).getFullName().trim().equals(alignment[j].getIdGroup().getIdentifier(taxa).getFullName().trim())) {
//                        throw new IllegalArgumentException("Selected chromosomes must have data for the same plant lines.");
//                    }
//                }
//            }
//        }
        try {
            String fullFileName = Utils.addSuffixIfNeeded(filename, ".hmp.txt", new String[]{".hmp.txt", ".hmp.txt.gz"});
            BufferedWriter bw = Utils.getBufferedWriter(fullFileName);
            bw.write("rs#");
            bw.write(delimChar);
            bw.write("alleles");
            bw.write(delimChar);
            bw.write("chrom");
            bw.write(delimChar);
            bw.write("pos");
            bw.write(delimChar);
            bw.write("strand");
            bw.write(delimChar);
            bw.write("assembly#");
            bw.write(delimChar);
            bw.write("center");
            bw.write(delimChar);
            bw.write("protLSID");
            bw.write(delimChar);
            bw.write("assayLSID");
            bw.write(delimChar);
            bw.write("panelLSID");
            bw.write(delimChar);
            bw.write("QCcode");
            bw.write(delimChar);
//            int numTaxa = alignment[0].getSequenceCount();
            int numTaxa = alignment.getSequenceCount();
            for (int taxa = 0; taxa < numTaxa; taxa++) {//finish filling out first row
                //not completely sure this does what I want, I need to access the
                //accession name from every alleleBLOB in bytes [52-201] but there
                //doesn't seem to be a method to access that in Alignment
//                String sequenceID = alignment[0].getIdGroup().getIdentifier(taxa).getFullName().trim();
                String sequenceID = alignment.getIdGroup().getIdentifier(taxa).getFullName().trim();
                sequenceID = sequenceID.replaceAll("\\s", "_");
                bw.write(sequenceID);
                if (taxa != numTaxa - 1) {
                    bw.write(delimChar);
                }
            }
            bw.write("\n");
//            for (int n = 0; n < alignment.length; n++) {
//            int numSites = alignment[n].getSiteCount();
            int numSites = alignment.getSiteCount();
            for (int site = 0; site < numSites; site++) {
//                StringBuilder sb = new StringBuilder(alignment[n].getSNPID(site)); //rs#
                StringBuilder sb = new StringBuilder(alignment.getSNPID(site));
//                    if (alignment[n].getSNPID(site).length() == 0) {
//                        sb.append("chr");
//                        sb.append(alignment[n].getLocusName(site));
//                        sb.append("_");
//                        sb.append(alignment[n].getPositionInLocus(site));
//                    }
                sb.append(delimChar);
//                SiteSummary ss = alignment[n].getSiteSummary(site);
                SiteSummary ss = alignment.getSiteSummary(site);
                int numAlleles = ss.getAlleles().length;
                //currently does not correctly display if numAlleles > 2
                if (numAlleles == 0) {
                    sb.append("NA"); //if data does not exist
                }
                for (int i = 0; i < Math.min(2, numAlleles) ; i++) {
                    if (i > 0) {
                        sb.append('/');
                    }
                    sb.append(ss.getAlleles()[i]); //alleles
                }
                sb.append(delimChar);
//                sb.append(alignment[n].getLocusName(site)); //chrom
                sb.append(alignment.getLocusName(site));
                sb.append(delimChar);
//                sb.append(alignment[n].getPositionInLocus(site)); // pos
                sb.append(alignment.getPositionInLocus(site));
                sb.append(delimChar);
                sb.append("+"); //strand
                sb.append(delimChar);
                sb.append("NA"); //assembly# not supported
                sb.append(delimChar);
                sb.append("NA"); //center unavailable
                sb.append(delimChar);
                sb.append("NA"); //protLSID unavailable
                sb.append(delimChar);
                sb.append("NA"); //assayLSID unavailable
                sb.append(delimChar);
                sb.append("NA"); //panelLSID unavailable
                sb.append(delimChar);
                sb.append("NA"); //QCcode unavailable
                sb.append(delimChar);
                for (int taxa = 0; taxa < numTaxa; taxa++) {
                    if (diploid == false) {
//                        byte baseIUPAC = alignment[n].getBase(taxa, site);
                        byte baseIUPAC = alignment.getBase(taxa, site);
                        if (mask.getMask(taxa, site) == 0x0) {
                            sb.append((char)baseIUPAC);
                        }
                        else if (mask.getMask(taxa, site) == 0x1) {
                            sb.append((char)(baseIUPAC + 32));
                        }
                    }
                    else {
//                            byte[] b = AllelePositionBLOBUtils.getSNPValueFromBase((char)alignment[n].getBase(taxa, site));
                            byte[] b = AllelePositionBLOBUtils.getSNPValueFromBase((char)alignment.getBase(taxa, site));
                            if(b.length==1) {
                                if (mask.getMask(taxa, site) == 0x0) {
                                    sb.append((char)b[0]);
                                    sb.append((char)b[0]);
                                }
                                else if (mask.getMask(taxa, site) == 0x1) {
                                    sb.append((char)(b[0] + 32));
                                    sb.append((char)(b[0] + 32));
                                }
                            }
                            else {
                                if (mask.getMask(taxa, site) == 0x0) {
                                    sb.append((char)b[0]);
                                    sb.append((char)b[1]);
                                }
                                else if (mask.getMask(taxa, site) == 0x1) {
                                    sb.append((char)(b[0] + 32));
                                    sb.append((char)(b[1] + 32));
                                }
                            }
                    }
                    if (taxa != (numTaxa - 1)) {
                        sb.append(delimChar);
                    }
                }
                bw.write(sb.toString());
                bw.write("\n");
            }
//            }
            bw.close();
        }
        catch (Exception e) {
            System.out.println("Error writing writeToHapmap: "+e);
            e.printStackTrace();
        }
    }

//    public static void writeToHapmap(Alignment alignment, boolean diploid, String filename, char delimChar) {
//        writeToHapmap(alignment.getAlignments(), diploid, filename, delimChar);
//    }

    public static void writeToLZMA(Pack1Alignment alignment, String filenameRoot) {
    }

    /**
     * Writes single set of value and pos BLOBs to a set of zip files
     * @param alignment
     * @param filenameRoot
     */
    public static void writeToZip(Alignment alignment, String filenameRoot) {
        Pack1Alignment align = (Pack1Alignment)alignment;
        ArrayList<byte[]>[] theBLOBList = new ArrayList[GdpdmBLOBUtils.totalBLOBtypes];
        for (int i = 0; i < theBLOBList.length; i++) {
            theBLOBList[i] = new ArrayList<byte[]>();
        }
        byte[][] alleleBLOBs = align.getAlleleBLOBs();
        for (byte[] b : alleleBLOBs) {
            theBLOBList[GdpdmBLOBUtils.alleleBLOBtype - 49].add(b);
        }
        theBLOBList[GdpdmBLOBUtils.allelePositionBLOBtype - 49].add(align.getVariableSitesBLOB());
        theBLOBList[GdpdmBLOBUtils.SNPIdBLOBtype - 49].add(align.getSNPidBLOB());
        GdpdmBLOBUtils.writeBLOBtoZip(theBLOBList, Utils.addSuffixIfNeeded(filenameRoot, ".BLOB.zip"));
    }

    /**
     * Writes an alignment to a gzip compressed file
     * @param alignment alignment to write to file
     * @param fileName name of gzip file
     */

    public static void writeToGZIP(Alignment alignment, String fileName) {
        Pack1Alignment align = (Pack1Alignment)alignment;
        ArrayList<byte[]>[] theBLOBList = new ArrayList[GdpdmBLOBUtils.totalBLOBtypes];
        for (int i = 0; i < theBLOBList.length; i++) {
            theBLOBList[i] = new ArrayList<byte[]>();
        }
        byte[][] alleleBLOBs = align.getAlleleBLOBs();
        for (byte[] b : alleleBLOBs) {
            theBLOBList[GdpdmBLOBUtils.alleleBLOBtype - 49].add(b);
        }
        theBLOBList[GdpdmBLOBUtils.allelePositionBLOBtype - 49].add(align.getVariableSitesBLOB());
        theBLOBList[GdpdmBLOBUtils.SNPIdBLOBtype - 49].add(align.getSNPidBLOB());
        GdpdmBLOBUtils.writeBLOBtoGZIP(theBLOBList, Utils.addSuffixIfNeeded(fileName, ".BLOB.gz"));
    }

    /**
     * Writes given set of alignments to a set of Plink files
     * @param alignment
     * @param filename
     * @param delimChar
     */
    public static void writeToPlink(Alignment alignment, String filename, char delimChar) {
        if (delimChar != ' ' && delimChar != '\t') {
            throw new IllegalArgumentException("Delimiter charater must be either a blank space or a tab.");
        }
//        for (int i = 0; i < alignment.length; i++) {
//            for (int j = i+1; j < alignment.length; j++) {
//                if (alignment[j].getSequenceCount() != alignment[i].getSequenceCount()) {
//                    throw new IllegalArgumentException("Selected chromosomes must have data for the same number of plant lines.");
//                }
//                for (int taxa = 0; taxa < alignment[i].getSequenceCount(); taxa++) {
//                    if (!alignment[i].getIdGroup().getIdentifier(taxa).getFullName().trim().equals(alignment[j].getIdGroup().getIdentifier(taxa).getFullName().trim())) {
//                        throw new IllegalArgumentException("Selected chromosomes must have data for the same plant lines.");
//                    }
//                }
//            }
//        }
        try {
            BufferedWriter MAPbw = new BufferedWriter(new FileWriter(Utils.addSuffixIfNeeded(filename, ".plk.map")), 1000000);
            int numSites = alignment.getSiteCount();
            for (int site = 0; site < numSites; site++) {
                MAPbw.write(alignment.getLocusName(site)); // chromosome name
                MAPbw.write(delimChar);
//                    if (alignment[i].getSNPID(site).length() == 0) {
//                        MAPbw.write("chr");
//                        MAPbw.write(alignment[i].getLocusName(site));
//                        MAPbw.write("_");
//                        MAPbw.write(alignment[i].getPositionInLocus(site));
//                    }
//                    else {
                    MAPbw.write(alignment.getSNPID(site)); // rs#
//                    }
                MAPbw.write(delimChar);
                MAPbw.write("-9"); // genetic distance unavailable
                MAPbw.write(delimChar);
                MAPbw.write(Integer.toString(alignment.getPositionInLocus(site))); // position
                MAPbw.write("\n");
            }
            MAPbw.close();
            BufferedWriter PEDbw = new BufferedWriter(new FileWriter(Utils.addSuffixIfNeeded(filename, ".plk.ped")), 1000000);
            // Compiled : Pattern
            Pattern splitter = Pattern.compile(":");
            int numTaxa = alignment.getSequenceCount();
            for (int taxa = 0; taxa < numTaxa; taxa++) {
                String[] name = splitter.split(alignment.getIdGroup().getIdentifier(taxa).getFullName().trim());
                if (name.length != 1) {
                    PEDbw.write(name[1]); // namelvl 1 if is available
                }
                else {
                    PEDbw.write("-9");
                }
                PEDbw.write(delimChar);
                //PEDbw.write(name[0]); // namelvl 0
                PEDbw.write(alignment.getIdGroup().getIdentifier(taxa).getFullName().trim()); // namelvl 0
                PEDbw.write(delimChar);
                PEDbw.write("-9"); // paternal ID unavailable
                PEDbw.write(delimChar);
                PEDbw.write("-9"); // maternal ID unavailable
                PEDbw.write(delimChar);
                PEDbw.write("-9"); // gender is both
                PEDbw.write(delimChar);
                PEDbw.write("-9"); // phenotype unavailable, changed to "-9" from "NA" due to error with Haploview compatibility
                PEDbw.write(delimChar);
                for (int site = 0; site < numSites; site++) {
                    byte[] b = AllelePositionBLOBUtils.getSNPValueForPlink(AllelePositionBLOBUtils.getSNPValueFromBase((char)alignment.getBase(taxa, site)));
                    PEDbw.write((char)b[0]);
                    PEDbw.write(delimChar);
                    PEDbw.write((char)b[b.length - 1]);
                    if (site != numSites - 1) {
                        PEDbw.write(delimChar);
                    }
                }
                PEDbw.write("\n");
            }
            PEDbw.close();
        }
        catch (Exception e) {
            myLogger.error("Error writing writeToPlink: "+e);
            throw new IllegalStateException("writeToPlink: " + e.getMessage());
        }
    }

//    public static void writeToPlink(Alignment alignment, String filename, char delimChar) {
//        writeToPlink(alignment.getAlignments(), filename, delimChar);
//    }

    /**
     * Writes given set of alignments to a set of Flapjack files
     * @param alignment
     * @param filename
     * @param delimChar
     */
    public static void writeToFlapjack(Alignment alignment, String filename, char delimChar) {
        if (delimChar != ' ' && delimChar != '\t') {
            throw new IllegalArgumentException("Delimiter charater must be either a blank space or a tab.");
        }
//        for (int i = 0; i < alignment.length; i++) {
//            for (int j = i+1; j < alignment.length; j++) {
//                if (alignment[j].getSequenceCount() != alignment[i].getSequenceCount()) {
//                    throw new IllegalArgumentException("Selected chromosomes must have data for the same number of plant lines.");
//                }
//                for (int taxa = 0; taxa < alignment[i].getSequenceCount(); taxa++) {
//                    if (!alignment[i].getIdGroup().getIdentifier(taxa).getFullName().trim().equals(alignment[j].getIdGroup().getIdentifier(taxa).getFullName().trim())) {
//                        throw new IllegalArgumentException("Selected chromosomes must have data for the same plant lines.");
//                    }
//                }
//            }
//        }
        try {
            BufferedWriter MAPbw = new BufferedWriter(new FileWriter(Utils.addSuffixIfNeeded(filename, ".flpjk.map")), 1000000);
            BufferedWriter DATbw = new BufferedWriter(new FileWriter(Utils.addSuffixIfNeeded(filename, ".flpjk.geno")), 1000000);
            int numSites = alignment.getSiteCount();
            for (int site = 0; site < numSites; site++) {
//                    if (alignment[i].getSNPID(site).length() == 0) {
//                        MAPbw.write("chr");
//                        MAPbw.write(alignment[i].getLocusName(site));
//                        MAPbw.write("_");
//                        MAPbw.write(alignment[i].getPositionInLocus(site));
//                    }
//                    else {
                    MAPbw.write(alignment.getSNPID(site)); // rs#
//                    }
                MAPbw.write(delimChar);
                MAPbw.write(alignment.getLocusName(site)); // chromosome name
                MAPbw.write(delimChar);
                MAPbw.write(Integer.toString(alignment.getPositionInLocus(site))); // position
                MAPbw.write("\n");
                DATbw.write(delimChar);
                DATbw.write(alignment.getSNPID(site));
            }
            MAPbw.close();
            DATbw.write("\n");
            int numTaxa = alignment.getSequenceCount();
            for (int taxa = 0; taxa < numTaxa; taxa++) {
                DATbw.write(alignment.getIdGroup().getIdentifier(taxa).getFullName().trim());
                DATbw.write(delimChar);
                for (int site = 0; site < numSites; site++) {
                    byte[] b = AllelePositionBLOBUtils.getSNPValueFromBase((char)alignment.getBase(taxa, site));
                    if (b.length == 1) {
                        if ((char)b[0] == 'N') {
                            DATbw.write('-');
                        }
                        else {
                            DATbw.write((char)b[0]);
                        }
                    }
                    else if (b.length == 2) {
                        DATbw.write((char)b[0]);
                        DATbw.write('/');
                        DATbw.write((char)b[1]);
                    }
                    else if (b.length == 3) {
                        DATbw.write((char)b[0]);
                        DATbw.write('/');
                        DATbw.write((char)b[1]);
                        DATbw.write('/');
                        DATbw.write((char)b[2]);
                    }
                    if (site != numSites - 1) {
                        DATbw.write(delimChar);
                    }
                }
                DATbw.write("\n");
            }
            DATbw.close();
        }
        catch (Exception e) {
            System.out.println("Error writing writeToFlapjack: "+e);
        }
    }

//    public static void writeToFlapjack(Alignment alignment, String filename, char delimChar) {
//        writeToFlapjack(alignment.getAlignments(), filename, delimChar);
//    }
}
