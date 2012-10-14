/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS;

import java.io.File;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AnnotatedAlignmentUtils;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.FilterAlignment;
import net.maizegenetics.pal.ids.IdGroup;

/**
 *
 * @author edbuckler
 */
public class GBSPipelinesKen {
    static String baseDir="/usr/local/maizediv/";
    static String qseqDir=baseDir+"illumina/433V4AAXX/Ken/fastq";
    static String qseqKey=baseDir+"illumina/433V4AAXX/Ken/433V4AAXX_key.txt";
    static String parsedTaxaDir=baseDir+"illumina/433V4AAXX/Ken/taxa";
    static String collapsedTaxaDir=baseDir+"illumina/433V4AAXX/Ken/taxacollapse/yeast";
    static String combinedFlowcellsReads=baseDir+"illumina/433V4AAXX/Ken/yeast/yeastReadCountsMin2_20100519.bin";
    static String myReadsByTaxa=baseDir+"illumina/433V4AAXX/Ken/yeast/yeastReadByTaxaMin2_20100519.txt";

    public GBSPipelinesKen() {
    }

    public static void main(String[] args) {
        createRefMarkerMap();

//         ReadsByTaxa theRBT=new ReadsByTaxa(theReadsByTaxaMinCount,true);
//         ReadsWPhysicalMap rwpmTest4=new ReadsWPhysicalMap(obsCommonReadsWithGenetic2,true);
//        // rwpmTest4.printRows(1000);
//         rwpmTest4.sortTable(false);
////
//
//         IdGroup goodTaxa=findRobustTaxa(theRBT,rwpmTest4);
//         for(int i=1; i<=10; i++) {
// //           createRobustHapMapfile(goodTaxa,theRBT,rwpmTest4,hapMap,(byte)i);
//         }


///**
//  CompareReadDistribution crd=new CompareReadDistribution(CompareReadDistribution.Analysis.EvalVirtualDigest,
//          refGenomeDigested, readsWorthMapping);
// */
//    CompareReadDistribution crd=new CompareReadDistribution(CompareReadDistribution.Analysis.EvalVirtualDigest, refGenomeDigested, combinedFlowcellsReads);
//

    }

    public static void createRobustHapMapfile(IdGroup goodTaxa, ReadsByTaxa theRBT, ReadsWPhysicalMap theRWMP,
            String baseFile, byte chromosome) {
        System.out.println("Creating Chromosome:"+chromosome);
        Alignment a = new ConvertGBSToAlignment(theRBT, theRWMP, 0.01, 1, chromosome);
        System.out.println("Base Alignment  Taxa:" + a.getSequenceCount() + " Sites:" + a.getSiteCount());

        AlignmentFilterByGBSUtils.hetsByLine(a, true);
        AlignmentFilterByGBSUtils.countCrossoversByLine(a);
        AlignmentFilterByGBSUtils.countDCO(a,false);
        System.out.println();

        if(goodTaxa==null) goodTaxa=a.getIdGroup();
        a = FilterAlignment.getInstance(a, goodTaxa);
        System.out.println("Taxa Filtered");
        System.out.println("TaxaFilter Alignment  Taxa:" + a.getSequenceCount() + " Sites:" + a.getSiteCount());
        AlignmentFilterByGBSUtils.hetsByLine(a, true);
        AlignmentFilterByGBSUtils.countCrossoversByLine(a);
        AlignmentFilterByGBSUtils.countDCO(a,false);
        System.out.println();

        System.out.println("Filtering Sites");
        a = AnnotatedAlignmentUtils.removeConstantSitesIgnoreGapsMissing(a);
        System.out.println("Codominant Alignment  Taxa:" + a.getSequenceCount() + " Sites:" + a.getSiteCount());
        int[] goodLowHetSites = AlignmentFilterByGBSUtils.getLowHetSNPs(a, 0.05, 20);

        a = FilterAlignment.getInstance(a, goodLowHetSites);
        System.out.println("SiteFilter Alignment  Taxa:" + a.getSequenceCount() + " Sites:" + a.getSiteCount());
        AlignmentFilterByGBSUtils.hetsByLine(a, true);
  
        int[] goodHighLDSites = AlignmentFilterByGBSUtils.getGoodSitesByLD(a, 0.9, 20, 10, true);
        a = FilterAlignment.getInstance(a, goodHighLDSites);
        System.out.println("LDFilter Alignment  Taxa:" + a.getSequenceCount() + " Sites:" + a.getSiteCount());

        int[] goodLowDCOSites = AlignmentFilterByGBSUtils.getLowDCOSNPs(a, 0.03, 0);
        a = FilterAlignment.getInstance(a, goodLowDCOSites);
        System.out.println("DCOFilter Alignment  Taxa:" + a.getSequenceCount() + " Sites:" + a.getSiteCount());
        AlignmentFilterByGBSUtils.hetsByLine(a, true);
        AlignmentFilterByGBSUtils.countCrossoversByLine(a);
        AlignmentFilterByGBSUtils.countDCO(a,false);
        System.out.println();

        ExportUtils.writeToHapmap(a, false, baseFile + ".c" + chromosome, '\t');
        
        MutableAlignmentForGBS mutAlign=new MutableAlignmentForGBS(a);
        mutAlign.imputeMissingData();
        System.out.println("Imputed Alignment  Taxa:" + mutAlign.getSequenceCount() + " Sites:" + mutAlign.getSiteCount());
        AlignmentFilterByGBSUtils.countDCO(mutAlign,false);
        AlignmentFilterByGBSUtils.countCrossoversByLine(mutAlign);

        ExportUtils.writeToHapmap(mutAlign, false, baseFile + ".imc" + chromosome, '\t');

    }

    public static IdGroup findRobustTaxa(ReadsByTaxa theRBT, ReadsWPhysicalMap theRWMP) {
        Alignment a = new ConvertGBSToAlignment(theRBT, theRWMP, 0.01, 1, (byte)99);
        System.out.println("Base Alignment  Taxa:" + a.getSequenceCount() + " Sites:" + a.getSiteCount());
        AlignmentFilterByGBSUtils.hetsByLine(a, true);
        AlignmentFilterByGBSUtils.countDCO(a,true);
        System.out.println("");

        a = AnnotatedAlignmentUtils.removeConstantSitesIgnoreGapsMissing(a);
        System.out.println("Codominant Alignment  Taxa:" + a.getSequenceCount() + " Sites:" + a.getSiteCount());

        System.out.println("Het Filtering");
        IdGroup goodLineSubset = AlignmentFilterByGBSUtils.getLowHetIdGroup(a, 0.03, a.getSiteCount()/10);
        a = FilterAlignment.getInstance(a, goodLineSubset);
        System.out.println("TaxaFilter Alignment  Taxa:" + a.getSequenceCount() + " Sites:" + a.getSiteCount());
        AlignmentFilterByGBSUtils.hetsByLine(a, true);

        int[] goodLowHetSites = AlignmentFilterByGBSUtils.getLowHetSNPs(a, 0.05, 20);
        a = FilterAlignment.getInstance(a, goodLowHetSites);
        System.out.println("SiteFilter Alignment  Taxa:" + a.getSequenceCount() + " Sites:" + a.getSiteCount());
        AlignmentFilterByGBSUtils.hetsByLine(a, true);
        System.out.println("");

        System.out.println("DCO Filtering");
        int[] goodLowDCOSites = AlignmentFilterByGBSUtils.getLowDCOSNPs(a, 0.03, 0);
        a = FilterAlignment.getInstance(a, goodLowDCOSites);
        System.out.println("DCOSiteFilter Alignment  Taxa:" + a.getSequenceCount() + " Sites:" + a.getSiteCount());
        AlignmentFilterByGBSUtils.countDCO(a,true);

        IdGroup goodLineSubsetByDCO = AlignmentFilterByGBSUtils.getLowDCOIdGroup(a, 0.02, a.getSiteCount()/10);
        a = FilterAlignment.getInstance(a, goodLineSubsetByDCO);
        System.out.println("DCOTaxaFilter Alignment  Taxa:" + a.getSequenceCount() + " Sites:" + a.getSiteCount());

        AlignmentFilterByGBSUtils.countDCO(a,false);

        return a.getIdGroup();
    }

    public static void createRefMarkerMap() {

       //Read and parse qseq or fastq file
//       ParseBarcodeFiles pbf=new ParseBarcodeFiles(qseqDir,qseqKey,parsedTaxaDir,0,true);
//       System.gc();

       //Take a GBSReads folder and collapse the folder by sorting and collapse identical reads
//       ReadCounts rc=new ReadCounts(parsedTaxaDir, collapsedTaxaDir, 1, true, false);
//       System.gc();

       //Fuse multiple taxa files into one, minimum frequency is 2, binary=true
       CombineReadCounts shff=new CombineReadCounts(collapsedTaxaDir,combinedFlowcellsReads, 2, true);
       System.gc();

       CreateReadsByTaxa crbt=new CreateReadsByTaxa(combinedFlowcellsReads, collapsedTaxaDir, myReadsByTaxa, false);


    }

}
