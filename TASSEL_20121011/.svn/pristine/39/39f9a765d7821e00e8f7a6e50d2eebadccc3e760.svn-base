/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.gbs.util;

import java.util.Collections;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import net.maizegenetics.pal.alignment.AbstractAlignment;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.Locus;
import net.maizegenetics.pal.datatype.DataType;
import net.maizegenetics.pal.datatype.IUPACNucleotides;

/**
 *  A proposal for storing alignment.  This data alignment is optimized for operations
 * involving lots of SNPs from a taxon - imputation, genetic distance, kinship, diversity, etc.
 * It is not optimized for LD or association mapping - use SBitAlignmentTest.
 *
 * Additionally, this structure should only used used when there are very few alleles (e.g. <=3).
 * SBitAlignmentTest is much more efficient for alignments with many alleles.
 *
 * @author edbuckler
 */
public class SBitAlignment extends AbstractAlignment {
    OpenBitSet[][] sData;
    protected int[] variableSites;  //positions of the variable sites
    protected byte[] SNPid; //SNP IDs
    protected byte[][] allele;
    protected Locus[] locus;  //this will be changed to locus object and eventually an list of loci
    private boolean isOnlyPolymorphic = true;
    protected int numSites, numTaxa;
    protected int maxAlleles=2;  //2 or 3 are reasonable numbers  i.e. only 2 or 3 most common alleles are used.  others are set to missing.


    public SBitAlignment(Alignment a) {
        super(a);
        numTaxa=getIdGroup().getIdCount();
        initVariableSites(a);
        initAlleles(a);
        loadAlleles(a);
        locus=a.getLoci();
    }

    private void initVariableSites(Alignment a) {
        this.numSites = a.getSiteCount();
        variableSites = new int[getSiteCount()];
        for (int i = 0; i < getSiteCount(); i++) {
            variableSites[i] = a.getPositionInLocus(i);
        }
    }

    private void initAlleles(Alignment a) {
        int homoHetAlleles=3;
        allele = new byte[homoHetAlleles][getSiteCount()];
        for (int i = 0; i < getSiteCount(); i++) {
            byte[] alleles = a.getAlleles(i);  //need to make sure this is getting haploid alleles
            for (int j = 0; j < maxAlleles; j++) {
                allele[j][i]=(j<alleles.length) ? alleles[j] : DataType.UNKNOWN_BYTE;
            }
            //this is hardwired for two alleles and third state is the het
            allele[2][i]=IUPACNucleotides.getDegerateSNPByteFromTwoSNPs(allele[0][i],allele[1][i]);
        }
    }

    private void loadAlleles(Alignment a) {
        sData=new OpenBitSet[maxAlleles][numSites];
        for (int al = 0; al < maxAlleles; al++) {
            for (int s = 0; s < numSites; s++) {
                sData[al][s]=new OpenBitSet(numTaxa);
            }
        }
        for (int s = 0; s < numSites; s++) {
            for (int t = 0; t < numTaxa; t++) {
                byte cb=a.getBase(t, s);
                //this is hard wired for a two allele model, shift when getBase returns the real alleles
                if((cb==allele[0][s])||(cb==allele[2][s])) sData[0][s].fastSet(t);
                if((cb==allele[1][s])||(cb==allele[2][s])) sData[1][s].fastSet(t);
            }
        }

    }

    public OpenBitSet getSiteBitsWithClone(int site, int allele) {
        return (OpenBitSet)sData[allele][site].clone();
    }

    public OpenBitSet getSiteBitsNoClone(int site, int allele) {
        return (OpenBitSet)sData[allele][site];
    }

    @Override
    public char getBaseChar(int taxon, int site) {
        return (char)getBase(taxon, site);
    }

    @Override
    public byte getBase(int taxon, int site) {
        int a=sData[0][site].fastGet(taxon)?1:0;
        a+=sData[1][site].fastGet(taxon)?2:0;
        if(a==0) return DataType.UNKNOWN_BYTE;
        return allele[--a][site];
    }

    @Override
    public byte getBase(int taxon, int site, int allele) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte getBase(int taxon, Locus locus, int physicalPosition, int allele) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

            /**
     * Return sorted list of alleles from highest frequency to lowest
     * at given site in alignment. Assumes IUPAC and diploid.
      * Resulting double dimension array
     * holds alleles (actually chars) in result[0].  And the counts
     * are in result[1]. All alleles are counted once.  No support
     * for higher ploidys.
     *
     * @param a alignment
     * @param site site
     * @param includeGAP whether to include GAP
     * @return sorted list of alleles and counts
     */
     @Override
    public int[][] getAllelesSortedByFrequency(int site, boolean includeGAP) {
        int[] stateCnt=new int[Byte.MAX_VALUE];
        for (int i = 0; i < getSequenceCount(); i++) {
            byte b=this.getBase(i, site);
            if(b==DataType.UNKNOWN_BYTE) continue;
            byte[] dipB=IUPACNucleotides.getDiploidValueFromIUPACCode(b);
            stateCnt[dipB[0]]++;
            stateCnt[dipB[1]]++;
        }
        TreeMap<Integer,Byte> alleleMap=new TreeMap<Integer,Byte>(Collections.reverseOrder());
        for (byte i = 0; i < stateCnt.length; i++) {
            if(stateCnt[i]>0) alleleMap.put(stateCnt[i],i);
        }
        int[][] result = new int[2][alleleMap.size()];
        Set<Map.Entry<Integer,Byte>> am=alleleMap.entrySet();
        int e=0;
        for (Map.Entry<Integer,Byte> eam: am) {
            result[0][e]=eam.getValue();
            result[1][e]=eam.getKey();
            e++;
        }
        return result;
    }

    @Override
    public String[] getSNPIDs() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getSNPID(int site) {
        return "S"+getLocus(site).getChromosomeName()+"_"+variableSites[site];
    }

    @Override
    public int getSiteCount() {
        return numSites;
    }

    @Override
    public int getLocusSiteCount(Locus locus) {
        throw new UnsupportedOperationException("Not supported yet.");
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
        return locus[0];
    }

    @Override
    public Locus[] getLoci() {
        return locus;
    }

    @Override
    public int getNumLoci() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

}
