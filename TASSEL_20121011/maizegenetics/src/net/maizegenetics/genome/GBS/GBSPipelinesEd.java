/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AnnotatedAlignmentUtils;
import net.maizegenetics.pal.alignment.CombineAlignment;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.FilterAlignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.ids.IdGroup;

/**
 *
 * @author edbuckler
 */
public class GBSPipelinesEd {
    static String baseDir="/Users/edbuckler/SolexaAnal/";
  //  static String baseDir="E:/SolexaAnal/";
//    static String qseqDir=baseDir+"NGG_IBM/434LFAAXX/qseq";
//    static String qseqKey=baseDir+"NGG_IBM/434LFAAXX/434LFAAXX.key.txt";
//    static String qseqDir=baseDir+"NGG_IBM/433V4AAXX/qseq";
//    static String qseqKey=baseDir+"NGG_IBM/433V4AAXX/433V4AAXX_key.txt";
    static String qseqDir=baseDir+"NGG_IBM/434GFAAXX/qseq";
    static String qseqKey=baseDir+"NGG_IBM/434GFAAXX/434GFAAXX_key.txt";
    static String taxaNameFile=baseDir+"NGG_IBM/IBMNames.txt";
    static String parsedTaxaDir=baseDir+"NGG_IBM/434LFAAXX/taxaQual";
    static String collapsedTaxaDir=baseDir+"NGG_IBM/434LFAAXX/taxacollapse";
    static String qualcollapsedTaxaDir=baseDir+"NGG_IBM/434LFAAXX/taxacollapsequal";
   // static String combinedHighQuality=baseDir+"NGG_IBM/434LFAAXX/counts/ReadCountsMin2_100430qual.bin";
    static String combinedFlowcellsReads=baseDir+"NGG_IBM/434LFAAXX/counts/RCBaseAGPV2_QualMin2_100503.bin";
    static String readsWorthMapping=baseDir+"NGG_IBM/434LFAAXX/counts/RCBaseAGPV2_QualMin5_100511.bin";
//    static String readsWorthMapping=baseDir+"NGG_IBM/434LFAAXX/counts/ReadCountsMin5_100430qual.bin";
    static String theReadsByTaxaAll=baseDir+"NGG_IBM/434LFAAXX/counts/ReadByTaxaMin2QualwNonQual_100511.bin";
    static String theReadsByTaxaLowQual=baseDir+"NGG_IBM/434LFAAXX/counts/ReadByTaxaLowQual_100518.bin";
    static String theReadsByTaxaLowQualMultiCells=baseDir+"NGG_IBM/434LFAAXX/counts/ReadByTaxaLowQualMulti_100518.bin";
    static String theReadsByTaxaMinLowQual=baseDir+"NGG_IBM/434LFAAXX/counts/ReadByTaxaMin20LowQual_100518.bin";

    static String theReadsByTaxaMinCount=baseDir+"NGG_IBM/434LFAAXX/counts/ReadByTaxaMin5QualwNonQual_100511.bin";
    static String refGenomeFasta=baseDir+"AGPv2/maize_agp_v2.fasta";
    static String refGenomeDigested=baseDir+"AGPv2/B73refVirtualV2.rwpm.bin";
    static String obsCommonReadsOnMap=baseDir+"NGG_IBM/434LFAAXX/counts/AGP_QualReadMap_100511.rwpm.bin";
     static String baseGeneticMap=baseDir+"AGPv2/IBM_644_NAM_SNP_genos_283_RILs_AGPv2_hapmap.txt";

    static String obsCommonReadsOnMapTest=baseDir+"NGG_IBM/434LFAAXX/counts/commonReadMap_100424test.rwpm.bin";
    static String obsCommonReadsWithGenetic1=baseDir+"NGG_IBM/434LFAAXX/counts/AGP_QualReadMapwGenetic_100503.rwpm.bin";
    static String obsCommonReadsWithGenetic2=baseDir+"NGG_IBM/434LFAAXX/counts/AGP_QualReadMapwGenetic_100511.rwpm.bin";
    static String obsCommonReadsWithGenetic3=baseDir+"NGG_IBM/434LFAAXX/counts/AGP_QualReadMapwLRFGenetic_100513.rwpm.bin";
    static String hapMap=baseDir+"NGG_IBM/434LFAAXX/counts/commonReadMap_100508t";
    static String hapMap2=baseDir+"SNP55K/SNP55hapmap05172010.txt";
    static String fullGeneticOuput=baseDir+"NGG_IBM/434LFAAXX/counts/geneticOutputWSkim_100528c3_1_55k.txt";

    public GBSPipelinesEd() {
    }

    public static void main(String[] args) {
 //        createRefMarkerMap();

//        parseQseqAgainstRefRBT();
//        appendParseQseqAgainstRefRBT();

//         ReadsByTaxa theRBT=new ReadsByTaxa(theReadsByTaxaMinCount,true);
//         ReadsWPhysicalMap rwpm=new ReadsWPhysicalMap(obsCommonReadsWithGenetic2,true);
//         rwpm.sortTable(false);
//         IdGroup goodTaxa=findRobustTaxa(theRBT,rwpm);
//         for(int i=1; i<=10; i++) {
//            createRobustHapMapfile(goodTaxa,theRBT,rwpm,hapMap,(byte)i);
//         }

         mapMarkersWithGenetics();


///**
//  CompareReadDistribution crd=new CompareReadDistribution(CompareReadDistribution.Analysis.EvalVirtualDigest,
//          refGenomeDigested, readsWorthMapping);
// */
//    CompareReadDistribution crd=new CompareReadDistribution(CompareReadDistribution.Analysis.EvalVirtualDigest, refGenomeDigested, combinedFlowcellsReads);
//

    }

    public static void parseQseqAgainstRefRBT() {
       String[] taxaNames=readStringArrayFromFile(new File(taxaNameFile));
       System.out.println(Arrays.toString(taxaNames));
       ReadCounts rc3=new ReadCounts(combinedFlowcellsReads, true);
       System.out.println("Building the Reads by Taxa table");
       ReadsByTaxa theRBT=new ReadsByTaxa(taxaNames, rc3);
       QseqToReadByTaxa.processDirectoryForCounting(new File(qseqDir), new File(qseqKey), theRBT, -1, false);
       theRBT.writeDistFile(new File(theReadsByTaxaLowQual), true, -1);
       // readsWorthMapping, collapsedTaxaDir, theReadsByTaxaAll
    }

    public static void appendParseQseqAgainstRefRBT() {
       System.out.println("Building the Reads by Taxa table");
//       ReadsByTaxa theRBT=new ReadsByTaxa(theReadsByTaxaLowQual, true);
       ReadsByTaxa theRBT=new ReadsByTaxa(theReadsByTaxaLowQualMultiCells, true);
       QseqToReadByTaxa.processDirectoryForCounting(new File(qseqDir), new File(qseqKey), theRBT, -1, false);
       theRBT.writeDistFile(new File(theReadsByTaxaLowQualMultiCells), true, -1);
       // readsWorthMapping, collapsedTaxaDir, theReadsByTaxaAll
    }

    public static void mapMarkersWithGenetics() {
  //      ReadsByTaxa theRBT=new ReadsByTaxa(theReadsByTaxaMinCount,true);
//          ReadsByTaxa theRBT=new ReadsByTaxa(theReadsByTaxaLowQual,true);
//          ReadsByTaxa theRBT=new ReadsByTaxa(theReadsByTaxaMinLowQual,true);
          ReadsByTaxa theRBT=new ReadsByTaxa(theReadsByTaxaLowQualMultiCells,true);


       //this BLAST and basic map only needs to be run when obsCommonReadsWithGenetic3 not present
//        ReadsWPhysicalMap rwpmTest2=MapHaplotypes.blastUsingPhysicalMap(new ReadsWPhysicalMap(refGenomeDigested,true),
//                  new ReadsWPhysicalMap(obsCommonReadsOnMap,true),34);
//         System.gc();
//         //this genetic check is fast an only evaluates flanking markers
//         rwpmTest2=MapHaplotypes.checkWorkingMapWithGenetic(baseGeneticMap,
//                  rwpmTest2, theRBT,false, 0.01);
//         rwpmTest2.writeCountFile(new File(obsCommonReadsWithGenetic3), Integer.MAX_VALUE, false,
//           false, Float.MAX_VALUE, true);
//
//        System.exit(0);

        ReadsWPhysicalMap rwpm=new ReadsWPhysicalMap(obsCommonReadsWithGenetic3,true);
//        Alignment[] gMap=new Alignment[10];
//        for(int chromosome=1; chromosome<=10; chromosome++) {
//            gMap[chromosome-1] = ImportUtils.readFromHapmap(hapMap+".imc" + chromosome+".hmp.txt");
//            System.out.println("Chromosome "+chromosome+" GeneticMap  Taxa:" + gMap[chromosome-1].getSequenceCount()
//                    + " Sites:" + gMap[chromosome-1].getSiteCount());
//         }
        Alignment[] gMap= ((CombineAlignment)ImportUtils.readFromHapmap(hapMap2)).getAlignments();
//        MapHaplotypes.mapWorkingMapWithGenetic(gMap, rwpm, theRBT, 0.05, 20, 10, null);
        MapHaplotypes.mapWorkingMapWithGenetic(gMap, rwpm, theRBT, 0.05, 20, 10, fullGeneticOuput);
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
        a=FilterAlignment.getInstance(a, goodTaxa);
     //   a = new FilterAlignment(a, goodTaxa);
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

         a=FilterAlignment.getInstance(a, goodLowHetSites);
     //   a = new FilterAlignment(a, goodLowHetSites);
        System.out.println("SiteFilter Alignment  Taxa:" + a.getSequenceCount() + " Sites:" + a.getSiteCount());
        AlignmentFilterByGBSUtils.hetsByLine(a, true);
  
        int[] goodHighLDSites = AlignmentFilterByGBSUtils.getGoodSitesByLD(a, 0.9, 20, 10, true);
        a=FilterAlignment.getInstance(a, goodHighLDSites);
     //   a = new FilterAlignment(a, goodHighLDSites);
        System.out.println("LDFilter Alignment  Taxa:" + a.getSequenceCount() + " Sites:" + a.getSiteCount());

        int[] goodLowDCOSites = AlignmentFilterByGBSUtils.getLowDCOSNPs(a, 0.03, 0);
        a=FilterAlignment.getInstance(a, goodLowDCOSites);
      //  a = new FilterAlignment(a, goodLowDCOSites);
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
        /**
     * Virtual digest of a reference
     */
/*        VirtualDigester vd=new VirtualDigester(new File(refGenomeFasta),new File(refGenomeDigested));
        System.gc();

      //Sort the physical map for quick searching and save
        ReadsWPhysicalMap rwpmVD=new ReadsWPhysicalMap(refGenomeDigested,true);
        rwpmVD.sortTable(true);
        rwpmVD.writeCountFile(new File(refGenomeDigested));
        System.gc();

    //   //Read and parse qseq file
       ParseBarcodeFiles pbf=new ParseBarcodeFiles(qseqDir,qseqKey,parsedTaxaDir,10);
       System.gc();

    //Take a GBSReads folder and collapse the folder by sorting and collapse identical reads
       ReadCounts rc=new ReadCounts(parsedTaxaDir, collapsedTaxaDir, 1, true, false);
       System.gc();

    //Fuse multiple taxa files into one, minimum frequency is 2, binary=true, simpleFilter=false
    //Use the physical genome as a base, and then add high quality reads to the top
       ReadsWPhysicalMap theRef=new ReadsWPhysicalMap(refGenomeDigested,true);
       theRef.sortTable(true);
       CombineReadCounts shff=new CombineReadCounts(theRef,
               qualcollapsedTaxaDir,combinedFlowcellsReads, 2, true, false);
       System.gc();
 
       //Found that using the reference reads as all greater than 2 in frequency was the best
       //so readsWorthMapping=combinedFlowcellsReads
       //now use all the reads, but only populate with those with perfect matches
       CreateReadsByTaxa crbt=new CreateReadsByTaxa(readsWorthMapping, collapsedTaxaDir, theReadsByTaxaAll, true);
       System.gc();

       ReadsByTaxa theRBT=new ReadsByTaxa(theReadsByTaxaAll,true);
       theRBT.writeDistFile(new File(theReadsByTaxaMinCount), true, 5);
       theRBT.writeReadCountFile(new File(readsWorthMapping), true, 5);
       System.gc();
*/
      //Create an empty mapping list from the readList
 
       ReadCounts rc3=new ReadCounts(readsWorthMapping, true);
       ReadsWPhysicalMap rwpm1=new ReadsWPhysicalMap(rc3);
       rwpm1.writeCountFile(new File(obsCommonReadsOnMap));
       System.gc();

         ReadsWPhysicalMap rwpmTest2=MapHaplotypes.blastUsingPhysicalMap(new ReadsWPhysicalMap(refGenomeDigested,true),
                  new ReadsWPhysicalMap(obsCommonReadsOnMap,true),3);
        rwpmTest2.writeCountFile(new File(obsCommonReadsOnMap));


         ReadsByTaxa theRBT2=new ReadsByTaxa(theReadsByTaxaMinCount,true);
         System.gc();
         //this genetic check is fast an only evaluates flanking markers
         rwpmTest2=MapHaplotypes.checkWorkingMapWithGenetic(baseGeneticMap,
                  rwpmTest2, theRBT2,false, 0.00001);
         rwpmTest2.writeCountFile(new File(obsCommonReadsWithGenetic2), Integer.MAX_VALUE, true,
           true, 0.01f, true);


    }

    public static String[] readStringArrayFromFile(File infile) {
        ArrayList<String> theStrings=new ArrayList<String>();
        try{
            BufferedReader br=new BufferedReader(new FileReader(infile));
            while(br.ready()) {
                String s=br.readLine().trim();
                if(s.isEmpty()==false) {theStrings.add(s);}
            }
            String[] results=new String[theStrings.size()];
            theStrings.toArray(results);
            return results;

        }catch(IOException e) {
            System.out.println("Error reading plain text file\n"+e.getMessage());
            return null;
        }
    }

}
