/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.HashSet;
import java.util.Set;
import net.maizegenetics.pal.alignment.AbstractAlignment;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AllelePositionBLOBUtils;
import net.maizegenetics.pal.alignment.GdpdmBLOBUtils;
import net.maizegenetics.pal.alignment.Locus;
import net.maizegenetics.pal.datatype.IUPACNucleotides;

/**
 *
 * @author edbuckler
 */
public class MutableAlignmentForGBS extends AbstractAlignment {
    public static final byte refAllele=GdpdmBLOBUtils.bases[0]; // A
    public static final byte altAllele=GdpdmBLOBUtils.bases[1]; // C
    public static final byte hetAllele=GdpdmBLOBUtils.bases[9]; // M = AC
    public static final byte missingAllele=GdpdmBLOBUtils.UNKNOWN_CHARACTER; // N
    private byte[][] seq;
    byte[][] imputedSeq;
    private int nDCOs;

    private Locus myLocus;
    private int[] variableSites;
    private byte[] strand;
    private int siteNumber;

    public enum FlankType {EMPTY, PURE_SAME_ALLELE, PURE_DIFF_ALLELE, MIX_OR_PURE_HET, PURE_REF_ALLELE, PURE_ALT_ALLELE}

    public MutableAlignmentForGBS(Alignment a) {
       super(a.getIdGroup(),new IUPACNucleotides());
       nDCOs = 0;
       myLocus=a.getLocus(0);
       siteNumber=a.getSiteCount();
       seq=new byte[getSequenceCount()][siteNumber];
       variableSites=new int[siteNumber];
       strand = new byte[siteNumber];
       for(int t=0; t<getSequenceCount(); t++) seq[t]=a.getBaseRange(t, 0, siteNumber-1);
       for(int s = 0; s <siteNumber; s++) {
           variableSites[s]=a.getPositionInLocus(s);
//           strand[s] = a.isPositiveStrand(s) ? (byte) '+' : (byte) '-';
       }
       System.out.println("New mutable alignment has " + siteNumber + " sites, and " + getSequenceCount() + " taxa");
    }

    public int imputeMissingData() {
        int nImputed = 0;
        for(int t=0; t<getSequenceCount(); t++) {
            byte startBase=-1;
            int startSite=0;
            for(int s = 0; s <getSiteCount(); s++) {
                if((getBase(t,s)==refAllele)||(getBase(t,s)==altAllele)) {
                    if(startBase==-1) {
                        startBase=getBase(t,s); //fill in the ends
                        if (s>0) ++nImputed;  // count the intial site, if imputed
                    }
                    if(startBase==getBase(t,s)) nImputed += fillRange(t, startSite,s,startBase);
                    startSite=s;
                    startBase=getBase(t,s);
                }
            }
            if (startBase==-1) {
                fillRange(t, startSite,getSiteCount()-1,missingAllele);  // nothing but missing data
            } else {
                nImputed += fillRange(t, startSite,getSiteCount()-1,startBase); //fill the ends
                if (startSite < (getSiteCount()-1)) ++nImputed;  // count the final site, if imputed
            }
        }
        return nImputed;
    }

    public int imputeMissingDataIncludingHets() {
        int nImputed = 0;
        for(int t=0; t<getSequenceCount(); t++) {
            byte startBase=-1;
            int startSite=0;
            for(int s = 0; s <getSiteCount(); s++) {
                if((getBase(t,s)==refAllele)||(getBase(t,s)==altAllele)||(getBase(t,s)==hetAllele)) {
                    if(startBase==-1) {
                        startBase=getBase(t,s); //fill in the ends
                        if (s>0) ++nImputed;  // count the intial site, if imputed
                    }
                    if(startBase==getBase(t,s)) nImputed += fillRange(t, startSite,s,startBase);
                    startSite=s;
                    startBase=getBase(t,s);
                }
            }
            if (startBase==-1) { 
                fillRange(t, startSite,getSiteCount()-1,missingAllele);  // nothing but missing data
            } else {
                nImputed += fillRange(t, startSite,getSiteCount()-1,startBase); //fill the ends
                if (startSite < (getSiteCount()-1)) ++nImputed;  // count the final site, if imputed
            }
        }
        return nImputed;
    }

    public int removeDCOs(int window) {
        for(int t=0; t<getSequenceCount(); t++) {
            for(int s = 0; s <getSiteCount(); s++) {
                setBase(t, s, callDCOsBasedOnFlanks(t,s,window));
            }
        }
        return nDCOs;
    }

    private byte callDCOsBasedOnFlanks(int taxon, int site, int window) {
        byte focusCall = getBase(taxon,site);
        if (focusCall == missingAllele) {
            return missingAllele;
        }

        int prevSite = site-1;
        int nextSite = site+1;
        byte[] upstreamFlank = new byte[window];
        for (int base=0; base<window; base++) {upstreamFlank[base]=missingAllele;}
        byte[] downstreamFlank = new byte[window];
        for (int base=0; base<window; base++) {downstreamFlank[base]=missingAllele;}

        // get the flankType of the upstream window of bases
        boolean collecting = true;
        int nFlankBases = 0;
        while (collecting) {
            if (prevSite>=0 && getBase(taxon,prevSite)!=missingAllele) {
                upstreamFlank[nFlankBases] = getBase(taxon,prevSite);
                ++nFlankBases;
            }
            --prevSite;
            if (nFlankBases==window || prevSite<0) {
                collecting = false;
            }
        }
        FlankType upstreamFlankType = getFlankType(focusCall, upstreamFlank);

        // get the flankType of the downstream window of bases
        collecting = true;
        nFlankBases = 0;
        while (collecting) {
            if (nextSite<seq[taxon].length && getBase(taxon,nextSite)!=missingAllele) {
                    downstreamFlank[nFlankBases] = getBase(taxon,nextSite);
                    ++nFlankBases;
            }
            ++nextSite;
            if (nFlankBases==window || nextSite>=seq[taxon].length) {
                collecting = false;
            }
        }
        FlankType downstreamFlankType = getFlankType(focusCall, downstreamFlank);

        // find DCOs, count them, and convert them to missing
        if (focusCall==hetAllele) {
            if (    (upstreamFlankType==FlankType.PURE_REF_ALLELE && downstreamFlankType==FlankType.PURE_REF_ALLELE)
                 || (upstreamFlankType==FlankType.PURE_ALT_ALLELE && downstreamFlankType==FlankType.PURE_ALT_ALLELE) ) {
                ++nDCOs;
                return missingAllele;
            }
            else {
                return hetAllele;
            }
        }
        else {
            if (upstreamFlankType==FlankType.PURE_DIFF_ALLELE && downstreamFlankType==FlankType.PURE_DIFF_ALLELE) {
                ++nDCOs;
                return missingAllele;
            }
            return focusCall;
        }
    }

    public void callHetSegments(int window) {
        imputedSeq = new byte[getSequenceCount()][siteNumber];
        for(int t=0; t<getSequenceCount(); t++) {
            for(int s = 0; s <getSiteCount(); s++) {
                imputedSeq[t][s] = callHetsBasedOnFlanks(t,s,window);
            }
        }
    }

    private byte callHetsBasedOnFlanks(int taxon, int site, int window) {
        // assumes that removeDCOs run beforehand (so we don't set 7 bases on either side of a solitary het to N)
        byte focusCall = getBase(taxon,site);
        if (focusCall == missingAllele) {
            return missingAllele;
        }

        int prevSite = site-1;
        int nextSite = site+1;
        byte[] upstreamFlank = new byte[window];
        for (int base=0; base<window; base++) {upstreamFlank[base]=missingAllele;}
        byte[] downstreamFlank = new byte[window];
        for (int base=0; base<window; base++) {downstreamFlank[base]=missingAllele;}

        // get the flankType of the upstream window of bases
        boolean collecting = true;
        int nFlankBases = 0;
        while (collecting) {
            if (prevSite>=0 && getBase(taxon,prevSite)!=missingAllele) {
                upstreamFlank[nFlankBases] = getBase(taxon,prevSite);
                ++nFlankBases;
            }
            --prevSite;
            if (nFlankBases==window || prevSite<0) {
                collecting = false;
            }
        }
        FlankType upstreamFlankType = getFlankType(focusCall, upstreamFlank);

        // get the flankType of the downstream window of bases
        collecting = true;
        nFlankBases = 0;
        while (collecting) {
            if (nextSite<seq[taxon].length && getBase(taxon,nextSite)!=missingAllele) {
                    downstreamFlank[nFlankBases] = getBase(taxon,nextSite);
                    ++nFlankBases;
            }
            ++nextSite;
            if (nFlankBases==window || nextSite>=seq[taxon].length) {
                collecting = false;
            }
        }
        FlankType downstreamFlankType = getFlankType(focusCall, downstreamFlank);

        // make the call
        if (focusCall==hetAllele) {
            if (    (upstreamFlankType==FlankType.PURE_REF_ALLELE && downstreamFlankType==FlankType.PURE_REF_ALLELE)
                 || (upstreamFlankType==FlankType.PURE_ALT_ALLELE && downstreamFlankType==FlankType.PURE_ALT_ALLELE) ) {
                ++nDCOs;
                return missingAllele;
            }
            else {
                return hetAllele;
            }
        }
        else {
            if (upstreamFlankType==FlankType.PURE_SAME_ALLELE && downstreamFlankType==FlankType.PURE_SAME_ALLELE) {
                return focusCall;
            }
            else if (upstreamFlankType==FlankType.PURE_SAME_ALLELE && downstreamFlankType==FlankType.PURE_DIFF_ALLELE) {
                return focusCall;
            }
            else if (upstreamFlankType==FlankType.PURE_DIFF_ALLELE && downstreamFlankType==FlankType.PURE_SAME_ALLELE) {
                return focusCall;
            }
            else if (upstreamFlankType==FlankType.PURE_DIFF_ALLELE && downstreamFlankType==FlankType.PURE_DIFF_ALLELE) {
                ++nDCOs;
                return missingAllele;
            }
            else if (upstreamFlankType==FlankType.PURE_SAME_ALLELE && downstreamFlankType==FlankType.MIX_OR_PURE_HET) {
                return missingAllele;
            }
            else if (upstreamFlankType==FlankType.MIX_OR_PURE_HET && downstreamFlankType==FlankType.PURE_SAME_ALLELE) {
                return missingAllele;
            }
            else if (upstreamFlankType==FlankType.EMPTY && downstreamFlankType==FlankType.MIX_OR_PURE_HET) {
                return missingAllele;
            }
            else if (upstreamFlankType==FlankType.MIX_OR_PURE_HET && downstreamFlankType==FlankType.EMPTY) {
                return missingAllele;
            }
            else if (upstreamFlankType==FlankType.PURE_DIFF_ALLELE && downstreamFlankType==FlankType.MIX_OR_PURE_HET) {
                return hetAllele;
            }
            else if (upstreamFlankType==FlankType.MIX_OR_PURE_HET && downstreamFlankType==FlankType.PURE_DIFF_ALLELE) {
                return hetAllele;
            }
            else if (upstreamFlankType==FlankType.MIX_OR_PURE_HET && downstreamFlankType==FlankType.MIX_OR_PURE_HET) {
                return hetAllele;
            }
            return focusCall;
        }
    }

    private FlankType getFlankType(byte focus, byte[] flankAlleles) {
        Set Alleles = new HashSet();
        for (int flankSite = 0; flankSite < flankAlleles.length; ++flankSite) {
            if (flankAlleles[flankSite] != missingAllele) {
                Alleles.add(flankAlleles[flankSite]);
            }
        }
        if (Alleles.size()==1) {
            if (Alleles.contains(hetAllele)) {
                return FlankType.MIX_OR_PURE_HET;
            }
            else if (Alleles.contains(focus)) {
                return FlankType.PURE_SAME_ALLELE;
            }
            else if (focus==hetAllele) {
                if (Alleles.contains(refAllele)) {
                    return FlankType.PURE_REF_ALLELE;
                }
                else {
                    return FlankType.PURE_ALT_ALLELE;
                }
            }
            else {
                return FlankType.PURE_DIFF_ALLELE;
            }
        }
        else if (Alleles.size()>1) {
            return FlankType.MIX_OR_PURE_HET;
        }
        return FlankType.EMPTY;
    }

    public int getnDCOs() {
        return nDCOs;
    }

    /**
     * @param taxa inclusive starting base
     * @param start inclusive starting base
     * @param end inclusive starting base
     * @param fillBase
     */
    private int fillRange(int taxa, int start, int end, byte fillBase) {
        int nImputed = 0;
        if(end<=start) return nImputed;
        for(int i=start; i<=end; i++) {
            setBase(taxa,i,fillBase);
            ++nImputed;
        }
        return nImputed-2;
    }

    public void writeToHapmapBC2S3(boolean diploid, boolean hetsCalled, String filename, char delimChar, int chr) {
        if (delimChar != ' ' && delimChar != '\t') {
            throw new IllegalArgumentException("Delimiter charater must be either a blank space or a tab.");
        }
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(filename + ".hmp.txt"), 1000000);
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
            int numTaxa = getSequenceCount();
            for (int taxa = 0; taxa < numTaxa; taxa++) {//finish filling out first row
                String sequenceID = getIdGroup().getIdentifier(taxa).getFullName().trim();
                bw.write(sequenceID);
                if (taxa != numTaxa - 1) {
                    bw.write(delimChar);
                }
            }
            bw.write("\n");
            int numSites = getSiteCount();
            for (int site = 0; site < numSites; site++) {
                int numAlleles = getAlleles(site).length;
                if (numAlleles>0 && getAlleles(site)[0] == 'C') {
                    StringBuilder sb = new StringBuilder(getSNPID(site)); //rs#
                    if (sb.substring(0, sb.length()).startsWith("null")) {
                        sb.replace(0, 4, ""+chr);
                    }
                    sb.append(delimChar);
                    //currently does not correctly display if numAlleles > 2
                    if (numAlleles == 0) {
                        sb.append("NA"); //if data does not exist
                    }
                    for (int i = 0; i < numAlleles ; i++) {
                        if (i > 0) {
                            sb.append('/');
                        }
                        sb.append((char) getAlleles(site)[i]); //alleles
                    }
                    sb.append(delimChar);
                    sb.append(getLocusName(site)); //chrom
                    sb.append(delimChar);
                    sb.append(getPositionInLocus(site)); // pos
                    sb.append(delimChar);
                    sb.append( (char) this.strand[site] );
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
                            byte baseIUPAC;
                            if (hetsCalled == false) {
                                baseIUPAC = getBase(taxa, site);
                                sb.append((char)baseIUPAC);
                            }
                            else {
                                baseIUPAC = imputedSeq[taxa][site];
                                sb.append((char)baseIUPAC);
                            }
                        }
                        else {
                            byte[] b;
                            if (hetsCalled == false) {
                                b = AllelePositionBLOBUtils.getSNPValueFromBase((char) getBase(taxa, site));
                            }
                            else {
                                b = AllelePositionBLOBUtils.getSNPValueFromBase((char) imputedSeq[taxa][site]);
                            }
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
            }
            bw.close();
        }
        catch (Exception e) {
            System.out.println("Error writing writeToHapmap: "+e);
        }
    }

    public int writeToHapmap(boolean diploid, boolean hetsCalled, String filename, char delimChar, int chr) {
        int nGenos = 0;
        if (delimChar != ' ' && delimChar != '\t') {
            throw new IllegalArgumentException("Delimiter charater must be either a blank space or a tab.");
        }
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(filename + ".hmp.txt"), 1000000);
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
            int numTaxa = getSequenceCount();
            for (int taxa = 0; taxa < numTaxa; taxa++) {//finish filling out first row
                String sequenceID = getIdGroup().getIdentifier(taxa).getFullName().trim();
                bw.write(sequenceID);
                if (taxa != numTaxa - 1) {
                    bw.write(delimChar);
                }
            }
            bw.write("\n");
            int numSites = getSiteCount();
            for (int site = 0; site < numSites; site++) {
                int numAlleles = getAlleles(site).length;
                if (numAlleles==2 || numAlleles==3) {  // getAlleles doesn't support hets
                    StringBuilder sb = new StringBuilder(getSNPID(site)); //rs#
                    if (sb.substring(0, sb.length()).startsWith("null")) {
                        sb.replace(0, 4, ""+chr);
                    }
                    sb.append(delimChar);
                    //currently does not correctly display if numAlleles > 2
                    if (numAlleles == 0) {
                        sb.append("NA"); //if data does not exist
                    }
                    for (int i = 0; i < numAlleles ; i++) {
                        if (i > 0) {
                            sb.append('/');
                        }
                        sb.append((char) getAlleles(site)[i]); //alleles
                    }
                    sb.append(delimChar);
                    sb.append(getLocusName(site)); //chrom
                    sb.append(delimChar);
                    sb.append(getPositionInLocus(site)); // pos
                    sb.append(delimChar);
                    sb.append( (char) this.strand[site] );
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
                            byte baseIUPAC;
                            if (hetsCalled == false) {
                                baseIUPAC = getBase(taxa, site);
                                sb.append((char)baseIUPAC);
                            }
                            else {
                                baseIUPAC = imputedSeq[taxa][site];
                                sb.append((char)baseIUPAC);
                            }
                            ++nGenos;
                        }
                        else {
                            byte[] b;
                            if (hetsCalled == false) {
                                b = AllelePositionBLOBUtils.getSNPValueFromBase((char) getBase(taxa, site));
                            }
                            else {
                                b = AllelePositionBLOBUtils.getSNPValueFromBase((char) imputedSeq[taxa][site]);
                            }
                            if(b.length==1) {
                                sb.append((char)b[0]);
                                sb.append((char)b[0]);
                            }
                            else {
                                sb.append((char)b[0]);
                                sb.append((char)b[1]);
                            }
                            ++nGenos;
                        }
                        if (taxa != (numTaxa - 1)) {
                            sb.append(delimChar);
                        }
                    }
                    bw.write(sb.toString());
                    bw.write("\n");
                }
            }
            bw.close();
        }
        catch (Exception e) {
            System.out.println("Error writing writeToHapmap: "+e);
        }
        return nGenos;
    }

    public void writeToHapmapBC2S3CA(boolean diploid, boolean hetsCalled, String filename, char delimChar, int chr) {
        if (delimChar != ' ' && delimChar != '\t') {
            throw new IllegalArgumentException("Delimiter charater must be either a blank space or a tab.");
        }
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(filename + ".hmp.txt"), 1000000);
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
            int numTaxa = getSequenceCount();
            for (int taxa = 0; taxa < numTaxa; taxa++) {//finish filling out first row
                String sequenceID = getIdGroup().getIdentifier(taxa).getFullName().trim();
                bw.write(sequenceID);
                if (taxa != numTaxa - 1) {
                    bw.write(delimChar);
                }
            }
            bw.write("\n");
            int numSites = getSiteCount();
            for (int site = 0; site < numSites; site++) {
                int numAlleles = getAlleles(site).length;
                boolean noVAllele = true;
                if (numAlleles == 3 && getAlleles(site)[2] == 'V') {noVAllele=false;}
                if (numAlleles == 4 && (getAlleles(site)[2] == 'V' || getAlleles(site)[3] == 'V')) {noVAllele=false;}
                if (numAlleles>1 && getAlleles(site)[0] == 'C' && getAlleles(site)[1] == 'A' && noVAllele) {
                    StringBuilder sb = new StringBuilder(getSNPID(site)); //rs#
                    if (sb.substring(0, sb.length()).startsWith("null")) {
                        sb.replace(0, 4, ""+chr);
                    }
                    sb.append(delimChar);
                    //currently does not correctly display if numAlleles > 2
                    if (numAlleles == 0) {
                        sb.append("NA"); //if data does not exist
                    }
                    for (int i = 0; i < numAlleles ; i++) {
                        if (i > 0) {
                            sb.append('/');
                        }
                        sb.append((char) getAlleles(site)[i]); //alleles
                    }
                    sb.append(delimChar);
                    sb.append(getLocusName(site)); //chrom
                    sb.append(delimChar);
                    sb.append(getPositionInLocus(site)); // pos
                    sb.append(delimChar);
                    sb.append( (char) this.strand[site] );
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
                            byte baseIUPAC;
                            if (hetsCalled == false) {
                                baseIUPAC = getBase(taxa, site);
                                sb.append((char)baseIUPAC);
                            }
                            else {
                                baseIUPAC = imputedSeq[taxa][site];
                                sb.append((char)baseIUPAC);
                            }
                        }
                        else {
                            byte[] b;
                            if (hetsCalled == false) {
                                b = AllelePositionBLOBUtils.getSNPValueFromBase((char) getBase(taxa, site));
                            }
                            else {
                                b = AllelePositionBLOBUtils.getSNPValueFromBase((char) imputedSeq[taxa][site]);
                            }
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
            }
            bw.close();
        }
        catch (Exception e) {
            System.out.println("Error writing writeToHapmap: "+e);
        }
    }



    @Override
    public char getBaseChar(int taxon, int site) {
        return (char)getBase(taxon,site);
    }

    @Override
    public byte getBase(int taxon, int site) {
        return seq[taxon][site];
    }

    public void setBase(int taxon, int site, byte base) {
        seq[taxon][site]=base;
    }

    @Override
    public byte getBase(int taxon, int site, int allele) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte getBase(int taxon, Locus locus, int physicalPosition, int allele) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public String[] getSNPIDs() {
        int sites = getSiteCount();
        String[] SNPids = new String[sites];
        for (int i = 0; i < sites; i++) {
            SNPids[i] = getSNPID(i);
        }
        return SNPids;
    }

    @Override
    public String getSNPID(int site) {
        return getLocus(site).getChromosomeName()+"_"+variableSites[site];
    }

    @Override
    public int getSiteCount() {
        return siteNumber;
    }

    @Override
    public int getLocusSiteCount(Locus locus) {
        return getSiteCount();
    }

    @Override
    public int getPositionInLocus(int site) {
        return variableSites[site];
    }

    @Override
    public int getSiteOfPhysicalPosition(int physicalPosition, Locus locus) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte getPositionType(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte[] getPositionTypes() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public Locus getLocus(int site) {
        return myLocus;
    }

    @Override
    public Locus[] getLoci() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getNumLoci() {
        return 1;
    }

}
