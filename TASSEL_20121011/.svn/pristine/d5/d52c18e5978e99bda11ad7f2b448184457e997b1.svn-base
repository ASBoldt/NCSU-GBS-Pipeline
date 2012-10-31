/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.solexa;

import cern.colt.list.LongArrayList;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AnnotatedAlignmentUtils;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.FilterAlignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.Pack1Alignment;
import net.maizegenetics.pal.datatype.DataType;
import net.maizegenetics.pal.datatype.IUPACNucleotides;
import net.maizegenetics.pal.distance.IBSDistanceMatrix;
import net.maizegenetics.pal.distance.WriteDistanceMatrix;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.IdGroupUtils;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.popgen.BasicImputation;
import net.maizegenetics.stats.MLM.Kinship;

/**
 *
 * @author edbuckler
 */
public class AlignmentFusion {
    long[] nameLong;
    int[] pos;
    int[] a1site;
    int[] a2site;

    public AlignmentFusion(String infile1, String infile2, String outfile) {
        System.out.println("Alignment Loading:"+infile1);
        Pack1Alignment p1a = (Pack1Alignment) ImportUtils.createPack1AlignmentFromFile(infile1, "", "");
        System.out.println("Alignment Loaded:"+infile1);
        System.out.printf("Sites %d Taxa %d %n", p1a.getSiteCount(),p1a.getSequenceCount());
        
        System.out.println("Alignment Loading:"+infile2);
        Pack1Alignment p2a = (Pack1Alignment) ImportUtils.createPack1AlignmentFromFile(infile2, "", "");
        System.out.println("Alignment Loaded:"+infile2);
        System.out.printf("Sites %d Taxa %d %n", p2a.getSiteCount(),p2a.getSequenceCount());

        int siteInP1AandP2A=0;
        int alleleAgree=0, alleleDisagree=0;
        IUPACNucleotides nuc=new IUPACNucleotides();

        int[] p12aRedirect = new int[p1a.getSequenceCount()];
        for (int t = 0; t < p1a.getSequenceCount(); t++) {
            p12aRedirect[t] = p2a.getIdGroup().whichIdNumber(p1a.getTaxaName(t));
        }
        for (int i = 0; i < p1a.getSiteCount(); i++) {
           int currPosition=p1a.getPositionInLocus(i);
           int p2pos=p2a.getSiteOfPhysicalPosition(currPosition, null);
           if(p2pos<0) {continue;}
           if((p1a.getMajorAllele(i)!=p2a.getMajorAllele(p2pos))||(p1a.getMinorAllele(i)!=p2a.getMinorAllele(p2pos))) {continue;}
           siteInP1AandP2A++;
           for (int t = 0; t < p1a.getSequenceCount(); t++) {
                if(p12aRedirect[t]<0) continue;
                byte a1=p1a.getBase(t, i);
                byte a2=p2a.getBase(p12aRedirect[t], p2pos);
                if(nuc.isHeterozygote(a1)||nuc.isHeterozygote(a2)||(a1==DataType.UNKNOWN_BYTE)||(a2==DataType.UNKNOWN_BYTE)) continue;
                if((a1==(byte)'+')||(a2==(byte)'+')||(a1==(byte)'-')||(a2==(byte)'-')) continue;
                if(a1==a2) {alleleAgree++;}
                else {alleleDisagree++;
                    System.out.print(p1a.getTaxaName(t)+" "+p2a.getTaxaName(p12aRedirect[t])+":");
                    System.out.println(i+" "+(char)a1+" v "+(char)a2);
                }
           }
        }
        double disagreeRate=(double)(alleleDisagree)/(double)(alleleAgree+alleleDisagree);
        System.out.println("Shared Sites:"+siteInP1AandP2A);
        System.out.printf("Alleles Agree: %d Disagree: %d Rate: %g %n", alleleAgree, alleleDisagree, disagreeRate);
        
    }

     public AlignmentFusion(String infile1, int minQ1, String infile2, int minQ2, boolean union, String outfile) {
        System.out.println("Alignment Loading:"+infile1);
        Pack1Alignment p1a = (Pack1Alignment) ImportUtils.createPack1AlignmentFromFile(infile1, "", "");
        System.out.println("Alignment Loaded:"+infile1);
        System.out.printf("Sites %d Taxa %d %n", p1a.getSiteCount(),p1a.getSequenceCount());
        int redudantP1APositions=0;
        for (int i = 1; i < p1a.getSiteCount(); i++) {
            if(p1a.getPositionInLocus(i)==p1a.getPositionInLocus(i-1)) {redudantP1APositions++;
               // System.out.println(p1a.getPositionInLocus(i));
            }

        }
        System.out.println("SitesWithMultiplePoly:"+redudantP1APositions);

        System.out.println("Alignment Loading:"+infile2);
        Pack1Alignment p2a = (Pack1Alignment) ImportUtils.createPack1AlignmentFromFile(infile2, "", "");
        System.out.println("Alignment Loaded:"+infile2);
        System.out.printf("Sites %d Taxa %d %n", p2a.getSiteCount(),p2a.getSequenceCount());
        redudantP1APositions=0;
        for (int i = 1; i < p2a.getSiteCount(); i++) {
            if(p2a.getPositionInLocus(i)==p2a.getPositionInLocus(i-1)) {redudantP1APositions++;
                //System.out.println(p2a.getPositionInLocus(i));
            }
        }
        System.out.println("SitesWithMultiplePoly:"+redudantP1APositions);

        createSitePolyList(p1a, minQ1, p2a, minQ2);

        int siteInP1AandP2A=0;
        int alleleAgree=0, alleleDisagree=0;
        IUPACNucleotides nuc=new IUPACNucleotides();


        int[] p12aRedirect = new int[p1a.getSequenceCount()];
        for (int t = 0; t < p1a.getSequenceCount(); t++) {
            p12aRedirect[t] = p2a.getIdGroup().whichIdNumber(p1a.getTaxaName(t));
        }
        try{
        BufferedWriter fileOut = new BufferedWriter(new FileWriter(outfile, false), 100000);
        fileOut.write("rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode\t");
        for (int t = 0; t < p1a.getSequenceCount(); t++) {
            fileOut.write(p1a.getIdGroup().getIdentifier(t).getFullName() + "\t");
        }
        fileOut.write("\n");
        for (int i = 0; i < pos.length; i++) {
           if((a1site[i]==-1)&&(a2site[i]==-1)) continue;
           if(((a1site[i]==-1)||(a2site[i]==-1))&&(union==false)) continue;
           StringBuilder sb = new StringBuilder();
           Alignment ca=(a1site[i]>-1)?p1a:p2a;
           int cSite=(a1site[i]>-1)?a1site[i]:a2site[i];
           //System.out.printf("%d %d %d %d %n",i, pos[i],a1site[i],a2site[i]);
           sb.append(ca.getSNPID(cSite)+"\t");
           sb.append((char)ca.getMajorAllele(cSite) + "/" + (char)ca.getMinorAllele(cSite) + "\t");
           sb.append(ca.getLocus(cSite) + "\t");
           sb.append((int) ca.getPositionInLocus(cSite) + "\t");
           sb.append("+\tAGPv1\tMaizeDiv\tSBS\tMHPv2");
           if(a1site[i]>-1) sb.append("C");
           if(a2site[i]>-1) sb.append("B");
           sb.append("\tNAMfnd\t");
           if(ca.hasSiteScores()) {sb.append((int)ca.getSiteScore(0, cSite)+"\t");}
                else {sb.append("QC+\t");}

           for (int t = 0; t < p1a.getSequenceCount(); t++) {
                char a1score=(a1site[i]>-1)?(char)p1a.getBaseChar(t, a1site[i]):DataType.UNKNOWN_CHARACTER;
                char a2score=((a2site[i]>-1)&&(p12aRedirect[t]>-1))?(char)p2a.getBaseChar(p12aRedirect[t], a2site[i]):DataType.UNKNOWN_CHARACTER;
                if(a1score==a2score) {sb.append(a1score+"\t");  alleleAgree++;}
                else if((a1score==DataType.UNKNOWN_CHARACTER)&&union) {sb.append(a2score+"\t");}
                else if((a2score==DataType.UNKNOWN_CHARACTER)&&union) {sb.append(a1score+"\t");}
                else  {sb.append(DataType.UNKNOWN_CHARACTER + "\t"); alleleDisagree++;}
            }
           //System.out.println(i+" "+sb.toString());
           fileOut.write(sb.toString());
           fileOut.write("\n");
        }
        fileOut.flush();
        fileOut.close();
         }
        catch(IOException ioe) {
            System.err.println("IO Error in AlignmentFusion(String infile1, int minQ1, String infile2, int minQ2, boolean union, String outfile)");
            ioe.printStackTrace();
        }
        double disagreeRate=(double)(alleleDisagree)/(double)(alleleAgree+alleleDisagree);
        System.out.println("Shared Sites:"+siteInP1AandP2A);
        System.out.printf("Alleles Agree: %d Disagree: %d Rate: %g %n", alleleAgree, alleleDisagree, disagreeRate);

    }

    private void createSitePolyList(Alignment a1, int minQ1, Alignment a2, int minQ2) {
        LongArrayList sp=new LongArrayList(a1.getSiteCount());
        System.out.println("Starting to create the combined - alignment 1");
        int c=a1.getSiteCount()-1;
        for (int i = 0; i < a1.getSiteCount(); i++) {
            if(a1.getSiteScore(0, i)<minQ1) continue;
            long x=((long)a1.getPositionInLocus(i)<<16)|((long)a1.getMajorAllele(i)<<8)|((long)a1.getMinorAllele(i));
            sp.add(x);
        }
        c=sp.size();
        sp.sort();
        System.out.println("Starting to create the combined - alignment 2");
        for (int i = 0; i < a2.getSiteCount(); i++) {
            if(a2.getSiteScore(0, i)<minQ2) continue;
            long x=((long)a2.getPositionInLocus(i)<<16)|((long)a2.getMajorAllele(i)<<8)|((long)a2.getMinorAllele(i));
            long y=((long)a2.getPositionInLocus(i)<<16)|((long)a2.getMinorAllele(i)<<8)|((long)a2.getMajorAllele(i));
            if((sp.binarySearchFromTo(x,0,c)<0)&&(sp.binarySearchFromTo(y,0,c)<0)) sp.add(x);
        }
        sp.trimToSize();
        sp.sort();
        System.out.println("CombinedListSize:"+sp.size());
        nameLong=sp.elements();
        pos=new int[sp.size()];
        a1site=new int[sp.size()]; Arrays.fill(a1site, -1);
        a2site=new int[sp.size()]; Arrays.fill(a2site, -1);
        for (int i = 0; i < a1.getSiteCount(); i++) {
            long x=((long)a1.getPositionInLocus(i)<<16)|((long)a1.getMajorAllele(i)<<8)|((long)a1.getMinorAllele(i));
            int index=Arrays.binarySearch(nameLong, x);
            if(index<0) continue;
            a1site[index]=i;
            pos[index]=a1.getPositionInLocus(i);
        }
        for (int i = 0; i < a2.getSiteCount(); i++) {
            long x=((long)a2.getPositionInLocus(i)<<16)|((long)a2.getMajorAllele(i)<<8)|((long)a2.getMinorAllele(i));
            int index=Arrays.binarySearch(nameLong, x);
            if(index<0) {
                long y=((long)a2.getPositionInLocus(i)<<16)|((long)a2.getMinorAllele(i)<<8)|((long)a2.getMajorAllele(i));
                index=Arrays.binarySearch(nameLong, y);
            }
            if(index<0) continue;
            a2site[index]=i;
            pos[index]=a2.getPositionInLocus(i);
        }
//        for (int i = 0; i < 100; i++) {
//            System.out.printf("%d %d %d %n",pos[i],a1site[i],a2site[i]);
//
//        }
    }

    /**
     * Returns the number of lines with homozygous states for the ref allele,
     * alternate allele, and combined for both.
     * @return
     */
    public int[] getHomozygousCounts(Alignment a1, int site) {
        int[] cnt=new int[3];
        for (int i=0; i<a1.getSequenceCount(); i++) {
            byte b=a1.getBase(i, site);
             if(b==DataType.UNKNOWN_BYTE) {}
             else if(b == a1.getMajorAllele(site)) {cnt[0]++;}
             else if(b == (byte) a1.getMinorAllele(site)) {cnt[1]++;}
             else {cnt[2]++;}//heterozygous
            
        }
        return cnt;
    }


    public static void filterImputeAlignments() {
         String outFile1="/Users/edbuckler/SolexaAnal/HapMapV2/jermpipe/SNPS201010/imputed/chr#.CSHLALLBGI.h90_f50.Q87Q87Union.hmp.txt";
         String outFile2="/Users/edbuckler/SolexaAnal/HapMapV2/jermpipe/SNPS201010/imputed/chr#.CSHLALLBGI.h90_f50.Q87Q87Union.imp.hmp.txt";
         String inDir1="/Users/edbuckler/SolexaAnal/HapMapV2/jermpipe/SNPS201010/fusion/chr#.CSHLALLBGI.h90_f12.Q87Q87Union.hmp.txt";
          for (int i = 3; i >= 1; i--) {
            String in1=inDir1.replaceFirst("#", ""+i);
            String out1=outFile1.replaceFirst("#", ""+i);
            String out2=outFile2.replaceFirst("#", ""+i);
              System.out.println("Reading:"+in1);
            Pack1Alignment p1a = (Pack1Alignment) ImportUtils.createPack1AlignmentFromFile(in1, "", "");
            System.out.println("Filtering:"+in1);
            Alignment a=AnnotatedAlignmentUtils.removeSitesBasedOnFreqIgnoreGapsMissing(p1a, 0.0001, 50);
            System.out.println("writing:"+out1);
            ExportUtils.writeToHapmap(a, false, out1, '\t');
            System.out.println("written:"+out1);
            a=p1a=null;
            System.out.println("Reading:"+out1);
            p1a = (Pack1Alignment) ImportUtils.createPack1AlignmentFromFile(out1, "", "");
            System.out.println("Imputing:"+out1);
            Pack1Alignment p2a=BasicImputation.imputeBySite(p1a, 31, 0);
            System.out.println("Writing:"+out2);
            ExportUtils.writeToHapmap(p2a, false, out2, '\t');
            System.out.println("Written:"+out2);

        }

    }

     public static void filterAlignmentsEstK() {
         String outFile2="/Users/edbuckler/SolexaAnal/HapMapV2/jermpipe/SNPS201010/K/chr#.CSHLALLBGI.h90_f50.Q87Q87Union.DistAll.hmp.txt";
         String outFile1="/Users/edbuckler/SolexaAnal/HapMapV2/jermpipe/SNPS201010/K/chr#.CSHLALLBGI.h90_f12.Q87Q87Union.KinAll.hmp.txt";
         String inDir1="/Users/edbuckler/SolexaAnal/HapMapV2/jermpipe/SNPS201010/fusion/chr#.CSHLALLBGI.h90_f12.Q87Q87Union.hmp.txt";
          for (int i = 3; i >= 1; i--) {
            String in1=inDir1.replaceFirst("#", ""+i);
            String out1=outFile1.replaceFirst("#", ""+i);
            String out2=outFile2.replaceFirst("#", ""+i);
              System.out.println("Reading:"+in1);
            Pack1Alignment p1a = (Pack1Alignment) ImportUtils.createPack1AlignmentFromFile(in1, "", "");
            for(int ep=0; ep<p1a.getPositionInLocus(p1a.getSiteCount()-1); ep+=10000000) {
                int start=p1a.getSiteOfPhysicalPosition(ep, null);
                if(start<0) start=-(start+1);
                int end=p1a.getSiteOfPhysicalPosition(ep+10000000, null);
                if(end<0) end=-(end+1);
                if(end>=p1a.getSiteCount()) end=p1a.getSiteCount()-1;
                Alignment a=AnnotatedAlignmentUtils.removeSitesOutsideRange(p1a, start, end);
                int twoDigitep=ep/10000000;
                out1=outFile1.replaceFirst("#", i+"_"+twoDigitep);
                out2=outFile2.replaceFirst("#", i+"_"+twoDigitep);
                System.out.println(out1+" s:"+start+" e:"+end);
                Kinship kin2 = new Kinship(a);
                WriteDistanceMatrix.saveDelimitedDistanceMatrix(kin2.getDm(), out1);
                IBSDistanceMatrix theIBSDM2= new IBSDistanceMatrix(a);
                WriteDistanceMatrix.saveDelimitedDistanceMatrix(theIBSDM2, out2);
              }

        }

    }

     public static void filterAlignmentsForMaize() {
         String outFile1="/Users/edbuckler/SolexaAnal/HapMapV2/jermpipe/SNPS201010/test/chr#.CSHLALLBGI.h90_f50.Q87Q87Union.maizeteo.hmp.txt";
         String inDir1="/Users/edbuckler/SolexaAnal/HapMapV2/jermpipe/SNPS201010/fusion/chr#.CSHLALLBGI.h90_f12.Q87Q87Union.hmp.txt";
         for (int i = 10; i >= 1; i--) {
            String in1=inDir1.replaceFirst("#", ""+i);
            String out1=outFile1.replaceFirst("#", ""+i);
            System.out.println("Reading:"+in1);
            Pack1Alignment p1a = (Pack1Alignment) ImportUtils.createPack1AlignmentFromFile(in1, "", "");
            boolean[] include=new boolean[p1a.getSequenceCount()];
            for (int t = 0; t < p1a.getSequenceCount(); t++) {
                  Identifier id=p1a.getIdGroup().getIdentifier(t);
                  if(id.getFullName().contains("MZ")) {include[t]=true;}
                  else if(id.getFullName().contains("TEO")) {include[t]=true;}
                    else {include[t]=false;}
              }

            IdGroup newGroup=IdGroupUtils.idGroupSubset(p1a.getIdGroup(),include);
            Alignment a=FilterAlignment.getInstance(p1a, newGroup);
            System.out.println("writing:"+out1);
            ExportUtils.writeToHapmap(a, false, out1, '\t');
            System.out.println("written:"+out1);


        }

    }

    public static void main(String[] args) throws Exception {
//        filterImputeAlignments();
 //       filterAlignmentsEstK();
        filterAlignmentsForMaize();
        System.exit(0);
//        String infile1="/Users/edbuckler/SolexaAnal/HapMapV2/test/comps/chr10.08232010.genotypes.log0ct_h0_f12.LR95H9.hmp.txt";
//        String infile1="/Users/edbuckler/SolexaAnal/HapMapV2/test/comps/chr10.08232010.genotypes.log0ct_h0_f12.LR8H9.hmp.txt";
//        String infile2="/Users/edbuckler/SolexaAnal/HapMapV2/test/comps/chr10.genotypes.bgi.log0ct_h0_f12.LR95H9.hmp.txt";
//        String infile2="/Users/edbuckler/SolexaAnal/HapMapV2/test/comps/chr10.genotypes.bgi.log0ct_h0_f12.LR8H9.hmp.txt";

//        String infile1="/Users/edbuckler/SolexaAnal/HapMapV2/jermpipe/20101010/cshlpe/test.cshl.hmp.txt";
//        String infile2="/Users/edbuckler/SolexaAnal/HapMapV2/jermpipe/20101010/cshlpe/test.bgi.hmp.txt";
/*
    String infile2="/Users/edbuckler/SolexaAnal/HapMapV2/jermpipe/20101010/cshlpe/chr10.BGI-Phase2.genotypes.log0ct_h0_f12.HET90.qscore.hmp.txt";
    String infile1="/Users/edbuckler/SolexaAnal/HapMapV2/jermpipe/20101010/cshlpe/chr10.CSHL-Phase2.101010.genotypes.log0ct_h0_f12.HET90.qscore.hmp.txt";


        String outfile="/Users/edbuckler/SolexaAnal/HapMapV2/jermpipe/20101010/fusion/chr10fusionQ88Q95int.hmp.txt";
        String parentDirectoryS="/Users/edbuckler/SolexaAnal/HapMapV2/test/comps/";


//        File parentDirectory = new File(parentDirectoryS);
//        File[] fileList = parentDirectory.listFiles();
//        for(File fn: fileList) {
//            if(!fn.getName().endsWith(".txt")) continue;
        //     AlignmentFusion hld3=new AlignmentFusion(infile1, infile2, outfile);
             AlignmentFusion hld3=new AlignmentFusion(infile1, 88, infile2, 95, false, outfile);
//        }
  */
          String inDir1="/Users/edbuckler/SolexaAnal/HapMapV2/jermpipe/SNPS201010/CSHL-ALL/HP1/chr#.CSHL-ALL.101010.genotypes.log1ct_h90_f12.hmp.txt";
          String inDir2="/Users/edbuckler/SolexaAnal/HapMapV2/jermpipe/SNPS201010/BGI/HP1/chr#.BGI-Phase2.genotypes.log1ct_h90_f12.hmp.txt";
          String outDir1="/Users/edbuckler/SolexaAnal/HapMapV2/jermpipe/SNPS201010/fusion/chr#.CSHLALLBGI.h90_f12.Q87Q87Union.hmp.txt";
          for (int i = 1; i <= 10; i++) {
            String in1=inDir1.replaceFirst("#", ""+i);
            String in2=inDir2.replaceFirst("#", ""+i);
            String out1=outDir1.replaceFirst("#", ""+i);
            AlignmentFusion hld3=new AlignmentFusion(in1, 87, in2, 87, true, out1);

        }

    }


}
