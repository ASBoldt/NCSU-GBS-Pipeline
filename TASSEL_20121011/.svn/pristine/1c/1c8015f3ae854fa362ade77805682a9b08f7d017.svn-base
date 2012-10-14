/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.gbs.pipeline;

import cern.jet.random.Binomial;
import edu.cornell.lassp.houle.RngPack.RandomJava;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import net.maizegenetics.gbs.homology.TagMatchFinder;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.gbs.util.OpenBitSet;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.tagdist.TagsByTaxaBitFileMap;
import net.maizegenetics.gbs.util.SBitAlignment;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.statistics.FisherExact;
import org.biojava3.core.util.ConcurrencyTools;

/**
 * Find sets of nearly homologuous tags, and then tests whether they are likely alleles
 * or systematic sequencing errors.  This really only works with inbred materials.
 * This also requires a TagsOnPhysicalMap file to report information.
 * @author edbuckler
 */
public class TagHomologyPhaseMT {
    static int minTaxaCnt=0;
    static int maxSize=600000;
    SBitAlignment[] theAnchorChrInBits;
    int[] tbt2anchRedirect;

    public TagHomologyPhaseMT(TagsOnPhysicalMap theTOPM, TagsByTaxa theTBT, String anchorFileName, String outHapMap) {
        final int cores=Runtime.getRuntime().availableProcessors();
        System.out.println("TagHomologyPhaseMT processors available:"+cores);
        TagMatchFinder theTMF=new TagMatchFinder(theTBT);
        loadAnchorMaps(anchorFileName);
        IdGroup anchorIdGroup=theAnchorChrInBits[0].getIdGroup();
        tbt2anchRedirect = new int[theTBT.getTaxaNames().length];
        for (int t = 0; t < theTBT.getTaxaCount(); t++) {
            tbt2anchRedirect[t] = anchorIdGroup.whichIdNumber(theTBT.getTaxaName(t));
        }
        System.out.println("Alignment Loaded:"+anchorFileName);
        ExecutorService pool = Executors.newCachedThreadPool();
        ThreadPoolExecutor tpe=(ThreadPoolExecutor)pool;
        FisherExact theFE=new FisherExact(theTBT.getTaxaCount()+10);

       int countSNPs=0, countLoci=0;
       long time=System.currentTimeMillis();
        System.out.println("Max Pool Size "+tpe.getMaximumPoolSize());
        System.out.println("Core Pool Size "+tpe.getCorePoolSize());
       int pauses=0;
       for (int i = 0; (i < theTBT.getTagCount())&&(countSNPs<maxSize); i++) {
           long[] cTag=theTBT.getTag(i);
           TreeMap<Integer,Integer> al=theTMF.findMatchesWithIntLengthWords(cTag, 1, false);
           int blastChr=-1, blastPos=-1;
           int indexInTOPM=0;
            if(theTOPM!=null) {
                indexInTOPM=theTOPM.getTagIndex(cTag);
                if(indexInTOPM>-1) {
                    blastChr=theTOPM.getChromosome(indexInTOPM);
                    blastPos=theTOPM.getStartPosition(indexInTOPM);
                }
            }
           OpenBitSet cBitDist=theTBT.getTaxaReadBitsForTag(i);
           int cCnt=(int)cBitDist.cardinality();
           for(Entry<Integer,Integer> ht: al.entrySet()) {
                int hitIndex=ht.getKey();
                if(hitIndex==i) continue;
                int div=ht.getValue();
                int hitIndexInTOPM=theTOPM.getTagIndex(theTBT.getTag(hitIndex));
                OpenBitSet hBitDist=theTBT.getTaxaReadBitsForTag(hitIndex);
                int hCnt=(int)hBitDist.cardinality();
                int cnt11=(int)OpenBitSet.intersectionCount(cBitDist, hBitDist);
                int cnt01=cCnt-cnt11;
                int cnt10=hCnt-cnt11;
                int cnt00=theTBT.getTaxaCount()-cnt11-cnt01-cnt10;
                double exp11=(double)hCnt*(double)cCnt/(double)theTBT.getTaxaCount();
                double relExp11=(double)cnt11/exp11;
                double p=theFE.getCumlativeP(cnt11, cnt01, cnt10, cnt00);
                if(p<0.99) {
//                    System.out.println(BaseEncoder.getSequenceFromLong(cTag));
//                    System.out.print(BaseEncoder.getSequenceFromLong(theTOPM.getTag(hitIndex)));
                    System.out.printf("con %d %d %d %d %d %d %g %g %g ",i, hitIndex, div, cCnt, hCnt, cnt11, exp11, relExp11, p);
                    OpenBitSet reCBitDist=getTagsInBits(cBitDist, tbt2anchRedirect, anchorIdGroup.getIdCount());
                    double[] bestR=getBestHit(reCBitDist);
                    String s=String.format("c %d %d %d %d %d %d %g %d %d ",indexInTOPM,blastChr, blastPos,
                        (int)bestR[0],(int)bestR[1],(int)bestR[2],bestR[3], (int)bestR[4], reCBitDist.cardinality());
                    if(bestR[3]<1) {System.out.print(s);}
                    reCBitDist=getTagsInBits(hBitDist, tbt2anchRedirect, anchorIdGroup.getIdCount());
                    bestR=getBestHit(reCBitDist);
                    s=String.format("h %d %d %d %d %d %d %g %d %d ",hitIndexInTOPM,theTOPM.getChromosome(hitIndexInTOPM), theTOPM.getStartPosition(hitIndexInTOPM),
                        (int)bestR[0],(int)bestR[1],(int)bestR[2],bestR[3], (int)bestR[4], reCBitDist.cardinality());
                    if(bestR[3]<1) {System.out.print(s);}
                    OpenBitSet orBD=(OpenBitSet)cBitDist.clone();
                    orBD.or(hBitDist);
                    reCBitDist=getTagsInBits(orBD, tbt2anchRedirect, anchorIdGroup.getIdCount());
                    bestR=getBestHit(reCBitDist);
                    s=String.format("o %d %d %d %d %d %d %g %d %d ",i,blastChr, blastPos,
                        (int)bestR[0],(int)bestR[1],(int)bestR[2],bestR[3], (int)bestR[4], reCBitDist.cardinality());
                    if(bestR[3]<1) {System.out.print(s);}
                    System.out.println();
                }
           }
       }
       System.out.println("Total Pauses or Yields:"+pauses);
       System.out.println("Main ThreadsCnt:"+Thread.activeCount()+" AlignmentThreadsCnt:"+ConcurrencyTools.getThreadPool().getActiveCount());
        try {
                // Wait a while for existing tasks to terminate
             if (!pool.awaitTermination(5, TimeUnit.SECONDS)) {
               pool.shutdownNow(); // Cancel currently executing tasks
               // Wait a while for tasks to respond to being cancelled
               if (!pool.awaitTermination(5, TimeUnit.SECONDS)) {System.err.println("Pool did not terminate");}
               else {System.out.println("Pool did terminate");}
             }
           } catch (InterruptedException ie) {
               System.err.println("Pool did not terminate");
             // (Re-)Cancel if current thread also interrupted
             pool.shutdownNow();
             // Preserve interrupt status
             Thread.currentThread().interrupt();
        }
       System.out.println("TC:"+Thread.activeCount()+" BJC"+ConcurrencyTools.getThreadPool().getActiveCount());
       ConcurrencyTools.shutdown();

       System.out.println("Wrote to hapmap:"+outHapMap);
        System.out.println("TC:"+Thread.activeCount());
    }

    private double[] getBestHit(OpenBitSet testBitSet) {
        double[] bestR={-1,-1,-1, 1, -1};
        double[][] theResults=new double[theAnchorChrInBits.length][];
        ScanChromosomeMTR[] threads = new ScanChromosomeMTR[theAnchorChrInBits.length];
        for (int cn=0; cn<theAnchorChrInBits.length; cn++) {
            threads[cn] = new ScanChromosomeMTR(theAnchorChrInBits[cn],testBitSet, cn, theResults);
            threads[cn].start();
        }
        for (int t=0; t<theAnchorChrInBits.length; t++) {
            try {threads[t].join();}
            catch (Exception e) {e.printStackTrace();}
        }
        for (int cn = 0; cn < theAnchorChrInBits.length; cn++) {
            double[] r=theResults[cn];
            if(r[3]<bestR[3]) bestR=r.clone();
        }
        return bestR;
    }

    public static OpenBitSet getTagsInBits(OpenBitSet origBitSet, int[] reDirect, int anchorTaxa) {
        OpenBitSet resultBitSet=new OpenBitSet(anchorTaxa);
        for (int j = 0; j < reDirect.length; j++) {
            if(reDirect[j]<0) continue;
            if (origBitSet.fastGet(j)) {resultBitSet.fastSet(reDirect[j]);}
        }
        return resultBitSet;
    }

    private void loadAnchorMaps(String anchorFileName) {
        if(anchorFileName.contains("*")) {
            int chrN=6;
            theAnchorChrInBits=new SBitAlignment[chrN];
            for (int i = 0; i < chrN; i++) {
                String file=anchorFileName.replace("*",""+(i+1));
                System.out.println("Reading:"+file);
                Alignment a = ImportUtils.readFromHapmap(file,""+(i+1));
                theAnchorChrInBits[i] = new SBitAlignment(a);
                System.out.printf("Chr %s Sites %d Taxa %d %n", theAnchorChrInBits[i].getLocus(0),
                        theAnchorChrInBits[i].getSiteCount(),theAnchorChrInBits[i].getSequenceCount());
            }
        } else {
            Alignment[] a = ImportUtils.readFromHapmap(anchorFileName).getAlignments();
            theAnchorChrInBits=new SBitAlignment[a.length];
            for (int i = 0; i < a.length; i++) {
                theAnchorChrInBits[i] = new SBitAlignment(a[i]);
    //            System.out.printf("Chr %d Sites %d Taxa %d %n", theAnchorChrInBits[i].getChr(),
    //                    theAnchorChrInBits[i].getNumSites(),theAnchorChrInBits[i].getNumTaxa());
                System.out.printf("Chr %s Sites %d Taxa %d %n", theAnchorChrInBits[i].getLocus(0),
                        theAnchorChrInBits[i].getSiteCount(),theAnchorChrInBits[i].getSequenceCount());
            }
        }
    }


    public static List<Object> getKeysFromValue(Map<?, ?> hm, Object value){
        List <Object>list = new ArrayList<Object>();
        for(Object o:hm.keySet()){
            if(hm.get(o).equals(value)) {
                list.add(o);
            }
        }
        return list;
    }

    private class ScanChromosomeMTR extends Thread {
        SBitAlignment refAlignment;
        OpenBitSet obsymj;
        Binomial binomFunc=new Binomial(5, 0.5, new RandomJava());
        int resultIndex;
        double[][] resultReport;

        public ScanChromosomeMTR(SBitAlignment refAlignment, OpenBitSet ymj, int resultIndex, double[][] resultReport) {
            this.refAlignment=refAlignment;
            this.obsymj=ymj;
            //obsymj=new OpenBitSet(ymj,ymj.length);
            this.resultIndex=resultIndex;
            this.resultReport=resultReport;
        }

        public void run() {
            //in the future this may call other threads to split the effort up more
            long tests=0;
            int bestSite=-1, countSig=0;
            double bestP=2;
            for (int i = 0; i < refAlignment.getSiteCount(); i++) {
             //   if(i%10000==0) System.out.println("scanChromosome chr:"+refAlignment.getChr()+"+site:"+i+" with Chr:"+refAlignment.getChr()+" test:"+tests);
                OpenBitSet mj=refAlignment.getSiteBitsNoClone(i, 0);
                OpenBitSet mn=refAlignment.getSiteBitsNoClone(i, 1);
                if(mn.lastCalculatedCardinality()>4) {
                    double p=TagCallerAgainstAnchorMT.testSites(obsymj, mj, mn, binomFunc);
                    if(p<bestP) {bestP=p; bestSite=i;}
                    if(p<0.0001) countSig++;
                }
                tests++;
            }
            int chr=Integer.parseInt(refAlignment.getLocus(bestSite).getChromosomeName());
            double[] result={chr, bestSite, refAlignment.getPositionInLocus(bestSite),bestP, countSig};
          //  if(bestP<0.000001) System.out.println(Arrays.toString(result));
            resultReport[resultIndex]=result;
        }

    }
    

     public static void main(String[] args) {
        System.out.println("Starting TagsToSNPByAlignmentMT");
        String tagMapFileS="/Users/edbuckler/SolexaAnal/GBS/reftags/14FCGBS.tg.ndup.bin";
 //       String tagsByTaxaS="/Users/edbuckler/SolexaAnal/GBS/test/bitIBM110210.dist.txt";
  //      String tagsByTaxaS="/Users/edbuckler/SolexaAnal/GBS/test/HiSeq282_110214.dist.txt";
        String tagsByTaxaS="/Users/edbuckler/SolexaAnal/GBS/dist/NAMsort_110303.bibin";
     //   String tagsByTaxaS="/Users/edbuckler/SolexaAnal/GBS/test/IBM110210.dist.txt";
    //    String outHapMap="/Users/edbuckler/SolexaAnal/GBS/test/HiSeq282wHet_110215.hmp.txt";
        String outHapMap="/Users/edbuckler/SolexaAnal/GBS/test/NAMwHet_110315mt.hmp.txt";
        String anchorMapFile="/Users/edbuckler/SolexaAnal/GBS/hmp/ch*_NAMwLowHetF10LD_110303imp.hmp.txt";

       TagsOnPhysicalMap theTOPM=new TagsOnPhysicalMap(tagMapFileS,true);  //Reads tag map file into an object
       //TagsByTaxa theTBT=new TagsByTaxaBit(tagsByTaxaS,FilePacking.Bit);
 //      TagsByTaxa theTBT=new TagsByTaxaBit(tagsByTaxaS,FilePacking.Text);
       TagsByTaxa theTBT=new TagsByTaxaBitFileMap(tagsByTaxaS);
   //    TagsByTaxa theTBT=null;
       theTOPM.sortTable(true);
 //      theTOPM.printRows(5, true, true);
       TagHomologyPhaseMT theTSBAMT=new TagHomologyPhaseMT(theTOPM, theTBT, anchorMapFile, outHapMap);
       
       
    }


}
