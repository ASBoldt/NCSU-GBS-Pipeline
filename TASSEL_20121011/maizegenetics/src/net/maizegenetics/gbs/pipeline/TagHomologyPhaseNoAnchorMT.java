/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.gbs.pipeline;

import cern.jet.random.Binomial;
import edu.cornell.lassp.houle.RngPack.RandomJava;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import net.maizegenetics.gbs.homology.TagMatchFinder;
import net.maizegenetics.gbs.util.OpenBitSet;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.tagdist.TagsByTaxaBitFileMap;
import net.maizegenetics.gbs.util.SBitAlignment;
import net.maizegenetics.pal.datatype.DataType;
import net.maizegenetics.pal.statistics.FisherExact;
import org.biojava3.core.util.ConcurrencyTools;

/**
 *Find sets of nearly homologous tags, and then tests whether they are likely alleles
 * or systematic sequencing errors.  This really only works with inbred materials.
 *
 * @author edbuckler
 */
public class TagHomologyPhaseNoAnchorMT {
    static int minTaxaCnt=0;
    static int maxSize=600000;


    public TagHomologyPhaseNoAnchorMT(TagsByTaxa theTBT, String outHapMap) {
        final int cores=Runtime.getRuntime().availableProcessors();
        System.out.println("TagHomologyPhaseMT processors available:"+cores);
        TagMatchFinder theTMF=new TagMatchFinder(theTBT);
        ExecutorService pool = Executors.newCachedThreadPool();
        ThreadPoolExecutor tpe=(ThreadPoolExecutor)pool;
        FisherExact theFE=new FisherExact(theTBT.getTaxaCount()+10);

       int countSNPs=0, countLoci=0;
       long time=System.currentTimeMillis();
        System.out.println("Max Pool Size "+tpe.getMaximumPoolSize());
        System.out.println("Core Pool Size "+tpe.getCorePoolSize());
        System.out.println("Comp	Tag1	Tag2	Divergence	Tag1Cnt1	Tag2Cnt2	Cnt11	exp11	Ratio11OE	P");
       int pauses=0;
	   ArrayList<filteredSnpPair> filSnpList = new ArrayList();
       for (int i = 0; (i < theTBT.getTagCount())&&(countSNPs<maxSize); i++) {
           long[] cTag=theTBT.getTag(i);
           TreeMap<Integer,Integer> al=theTMF.findMatchesWithIntLengthWords(cTag, 1, false);
           OpenBitSet cBitDist=theTBT.getTaxaReadBitsForTag(i);
           int cCnt=(int)cBitDist.cardinality();
           for(Entry<Integer,Integer> ht: al.entrySet()) {
                int hitIndex=ht.getKey();
                if(hitIndex==i) continue;  //prevent self match
                int div=ht.getValue();
                OpenBitSet hBitDist=theTBT.getTaxaReadBitsForTag(hitIndex);
                int hCnt=(int)hBitDist.cardinality();
                int cnt11=(int)OpenBitSet.intersectionCount(cBitDist, hBitDist);
                int cnt01=cCnt-cnt11;
                int cnt10=hCnt-cnt11;
                int cnt00=theTBT.getTaxaCount()-cnt11-cnt01-cnt10;
                double exp11=(double)hCnt*(double)cCnt/(double)theTBT.getTaxaCount();
                double relExp11=(double)cnt11/exp11;
                double p=theFE.getCumlativeP(cnt11, cnt01, cnt10, cnt00);
                double propPresent=(double)(cCnt+hCnt)/(double)theTBT.getTaxaCount();
                if((p<0.01)&&(relExp11<0.2)&&(propPresent>0.2)) {
					filteredSnpPair fsp = new filteredSnpPair(i, hitIndex, div, cCnt, hCnt, cnt11, exp11, relExp11, p, cBitDist, hBitDist);
					fsp.swap();
					filSnpList.add(fsp);
//                    System.out.println(BaseEncoder.getSequenceFromLong(cTag));
//                    System.out.print(BaseEncoder.getSequenceFromLong(theTOPM.getTag(hitIndex)));
  //                  System.out.printf("con %d %d %d %d %d %d %g %g %g ",i, hitIndex, div, cCnt, hCnt, cnt11, exp11, relExp11, p);
//                    System.out.print(" ");
//                    System.out.println(bitsToPseudoSeq(theTBT.getTaxaCount(),cBitDist,hBitDist));
                }
           }
       }
	   filteredSnpPair[] filSnps = filSnpList.toArray(new filteredSnpPair[filSnpList.size()]);
	   System.out.println(filSnpList.size());
	   System.out.println(filSnps.length);
	   Arrays.sort(filSnps);
	   System.out.println(filSnps.length);
	   try {
			BufferedWriter bw = new BufferedWriter (new FileWriter(outHapMap), 65536);
			bw.write("rs	alleles	chrom	pos	strand	assembly	center	protLSID	assayLSID	panelLSID	QCcode	");
			for (int i = 0; i < theTBT.getTaxaCount()-1; i++) {
				bw.write(theTBT.getTaxaNames()[i]+ "\t");
			}
			bw.write(theTBT.getTaxaNames()[theTBT.getTaxaCount()-1]);
			bw.newLine();
			int count = 0;
			for (int i = 0; i < filSnps.length; i += 2) {
				count++;
				bw.write(filSnps[i].mkStr(count, 1, count, theTBT.getTaxaCount()));
				bw.newLine();
			}
			bw.flush();
			bw.close();
	   }
	   catch (Exception e) {
			System.out.println(e.toString());
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
	class filteredSnpPair implements Comparable <filteredSnpPair> {
		int queryIndex;
		int hitIndex;
		int div;
		int cCnt;
		int hCnt;
		int cnt11;
		double exp11;
		double relExp11;
		double p;
		OpenBitSet cBitDist;
		OpenBitSet hBitDist;
		filteredSnpPair(int queryIndex, int hitIndex, int div, int cCnt, int hCnt, int cnt11, double exp11, double relExp11, double p, OpenBitSet cBitDist,  OpenBitSet hBitDist) {
			this.queryIndex = queryIndex;
			this.hitIndex = hitIndex;
			this.div = div;
			this.cCnt = cCnt;
			this.hCnt = hCnt;
			this.cnt11 = cnt11;
			this.exp11 = exp11;
			this.relExp11 = relExp11;
			this.p = p;
			this.cBitDist = cBitDist;
			this.hBitDist = hBitDist;
		}
		public String mkStr (int ID, int chr, int posi, int taxaCount) {
			StringBuilder sb = new StringBuilder();
			sb.append(ID).append("\tA/C\t").append(chr).append("\t").append(posi).append("\t+\tNA\tSWGDiv\tGBS\tSWGV1\tSWGPop\tQC+\t");
			sb.append(bitsToPseudoSeq(taxaCount,cBitDist,hBitDist));
			sb.deleteCharAt(sb.length()-1);
			String str = sb.toString();
			return str;
		}
		public void swap () {
			if (queryIndex > hitIndex) {
				int medium;
				OpenBitSet obs;
				medium = queryIndex;
				queryIndex = hitIndex;
				hitIndex = medium;
				medium = cCnt;
				hCnt = cCnt;
				cCnt = medium;
				obs = cBitDist;
				cBitDist = hBitDist;
				hBitDist = obs;
			}
		}
		public int compareTo (filteredSnpPair o) {
			return queryIndex - o.queryIndex;
		}
	}
    public static String bitsToPseudoSeq(int numTaxa, OpenBitSet t1, OpenBitSet t2) {
        StringBuilder sb=new StringBuilder();
        String[] allele={"A","C","R"};
        for (int i = 0; i < numTaxa; i++) {
            int a=t1.fastGet(i)?1:0;
            a+=t2.fastGet(i)?2:0;
            if(a==0) {sb.append(DataType.UNKNOWN_CHARACTER);}
            else {sb.append(allele[--a]);}
            sb.append("\t");
        }
        return sb.toString();
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

        //String tagsByTaxaS="/Users/edbuckler/SolexaAnal/GBS/dist/NAMsort_110303.bibin";
		String tagsByTaxaS = "M:/switch/ReadByTaxaBit/TagByTaxaBitAll.txt";
     //   String tagsByTaxaS="/Users/edbuckler/SolexaAnal/GBS/test/IBM110210.dist.txt";
    //    String outHapMap="/Users/edbuckler/SolexaAnal/GBS/test/HiSeq282wHet_110215.hmp.txt";
        //String outHapMap="/Users/edbuckler/SolexaAnal/GBS/test/NAMwHet_110315mt.hmp.txt";
		String outHapMap = "M:/switch/Hapmap/outHapMAP.txt";


       //TagsByTaxa theTBT=new TagsByTaxaBit(tagsByTaxaS,FilePacking.Bit);
 //      TagsByTaxa theTBT=new TagsByTaxaBit(tagsByTaxaS,FilePacking.Text);
       TagsByTaxa theTBT=new TagsByTaxaBitFileMap(tagsByTaxaS);

       TagHomologyPhaseNoAnchorMT theTSBAMT=new TagHomologyPhaseNoAnchorMT(theTBT, outHapMap);
       
       
    }


}
