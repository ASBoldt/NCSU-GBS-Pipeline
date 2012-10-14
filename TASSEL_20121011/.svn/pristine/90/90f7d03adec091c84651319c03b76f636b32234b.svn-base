package net.maizegenetics.gbs.pipeline;

import cern.jet.random.Binomial;
import edu.cornell.lassp.houle.RngPack.RandomJava;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.lang.reflect.Array;
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
//import net.maizegenetics.gbs.maps.TagsOnContigs;
import net.maizegenetics.gbs.util.OpenBitSet;
import net.maizegenetics.gbs.tagdist.Tags;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.tagdist.TagsByTaxaBitFileMap;
import net.maizegenetics.gbs.util.SBitAlignment;
import net.maizegenetics.genome.BaseEncoder;
import net.maizegenetics.genome.GBS.deNovoContigs.deNovoContigTags;
import net.maizegenetics.pal.datatype.DataType;
import net.maizegenetics.pal.statistics.FisherExact;
import org.biojava3.core.util.ConcurrencyTools;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.BinomialDistributionImpl;
import org.apache.commons.math.stat.inference.ChiSquareTestImpl;

/**
 *
 * @author jessepoland
 */
public class TagsToSNPsNoAnchor {
    static int minTaxaCnt=0;
    static int maxSize=100000;
    static boolean numericOutput = false;
    


    /**
     * @param theTBT
     * @param outHapMap
     * @param minDiff
     * @param maxMissingData
     * @param minAlleleFreq
     * @param isDHpopulation, asking if population is double haploid.  If so, will reject any SNPs that are heterozygous, 
     */
    public TagsToSNPsNoAnchor(TagsByTaxa theTBT, String outHapMap, int divergence, double maxMissingData, double minAlleleFreq, double maxHet, boolean isDHpopulation, boolean isBiparental, boolean callHets, double pVal) {
        
    	System.gc();
    	//do some checks
    	if(maxMissingData>1 || maxMissingData<0){
    		System.out.println("Max Missing Data (maxMissingData) must be between 0 and 1.  Entered Value: " + maxMissingData);
    		System.exit(0);
    	}
    	if(minAlleleFreq>0.5 || minAlleleFreq<0){
    		System.out.println("Minimum Allele Frequency (minAlleleFreq) must be between 0 and 0.5.  Entered Value: " + minAlleleFreq);
    		System.exit(0);
    	}
    	
    	
    	final int cores=Runtime.getRuntime().availableProcessors();
        System.out.println("TagHomologyPhaseMT processors available:"+cores);
        TagMatchFinder theTMF=new TagMatchFinder(theTBT);
        ExecutorService pool = Executors.newCachedThreadPool();
        ThreadPoolExecutor tpe=(ThreadPoolExecutor)pool;
        FisherExact theFE=new FisherExact(theTBT.getTaxaCount()+10);
        ChiSquareTestImpl theChiSq = new ChiSquareTestImpl();

        int nSNPs=0;
       int countSNPs=0, countLoci=0;
       long time=System.currentTimeMillis();
       System.out.println("Max Pool Size "+tpe.getMaximumPoolSize());
       System.out.println("Core Pool Size "+tpe.getCorePoolSize());
       System.out.println("Comp	Tag1	Tag2	Divergence	Tag1Cnt1	Tag2Cnt2	Cnt11	exp11	Ratio11OE	P");
       
       int pauses=0;
       ArrayList<filteredSnpPair> filSnpList = new ArrayList<filteredSnpPair>(0);
       
       //deNovoContigTags tagsOnCtg = new deNovoContigTags("/Users/jpoland/Genome/Barley/morex_TagsFromContigs.bin", true);
       //tagsOnCtg.getTags()[0].getRead();
       //TagsOnContigs theCtgTags = new TagsOnContigs("/Users/jpoland/Genome/Barley/tagsOnContig.bin", true, false);
       //TagMatchFinder theCtgTMF=new TagMatchFinder(theCtgTags);
       
       boolean[] tested = new boolean[theTBT.getTagCount()]; 
       
//	   for (int i = 0; (i < 100000)&&(countSNPs<maxSize); i++) {
	 
			   

        	   try {
        			BufferedWriter bw = new BufferedWriter (new FileWriter(outHapMap), 65536);
        			bw.write("rs	alleles	chrom	pos	strand	assembly	center	protLSID	assayLSID	panelLSID	QCcode	");
        			for (int j = 0; j < theTBT.getTaxaCount()-1; j++) {
        				bw.write(theTBT.getTaxaNames()[j]+ "\t");
        			}
        			bw.write(theTBT.getTaxaNames()[theTBT.getTaxaCount()-1]);
        			bw.newLine();
        			int count = 0;

        			int nTaxa = theTBT.getTaxaCount();
        			OpenBitSet reference = new OpenBitSet(theTBT.getTagCount());
        			
        			for (int i = 0; (i < theTBT.getTagCount())&&(nSNPs<maxSize); i++) {   
        			//	for (int i = 0; i < 100000; i++) { 
        				   long[] tagA=theTBT.getTag(i);
        				   OpenBitSet bitDistA=theTBT.getTaxaReadBitsForTag(i); // the taxa with tagA
        				   
        				   //if(bitDistA.cardinality() < nTaxa*0.3*(1-maxMissingData)) continue;  // skip the rare tags, use tags that show up at reasonable frequency (major allele)
        				   //if(bitDistA.cardinality() < 8) continue;  // skip the rare tags, use tags that show up at reasonable frequency (major allele)
        				   if(reference.get(i)) continue;  // check if this tag has been used as a reference, skip if so
        				   reference.set(i); // set this tags a being a reference
 
        				   TreeMap<Integer,Integer> al=theTMF.findMatchesWithIntLengthWords(tagA, divergence, false);  //TreeMap is a set of Tag pairs(query Tag, hit Tag). There is only 1 mismatch between each pair
        				   if(al.size()<2) continue; // skip stuff that only has one match (self)

         				   int refCnt=0;
         				   for(Entry<Integer,Integer> ht: al.entrySet()) {
        					   if(reference.get(ht.getKey())) refCnt++;
        				   }
        				   if(refCnt>1) continue;  // skip as there was already a reference tag for this TreeMap
        				   
        				   
        				   // populate an array of the tags with hits and their index
        				   ArrayList<byte[]> tags = new ArrayList<byte[]>(0);
        				   int[] hitIdx = new int[al.size()]; 
        				   //int refTagIdx = Integer.MAX_VALUE;
        				   int hi = 0;
        				   for(Entry<Integer,Integer> ht: al.entrySet()) {
        					   hitIdx[hi] = ht.getKey();
        					   hi++;
        					   long[] t = theTBT.getTag(ht.getKey());
        					   tags.add(BaseEncoder.getByteSeqFromLong(t));  
        					   //if(ht.getKey()<refTagIdx){refTagIdx = ht.getKey();} //set the first tag as the reference
        				   }
        				   //if(refTagIdx<i) continue; //skip tags that were previously tested
        				   
        				   

        				   
        				   // test each base in the set of tags for differences
        				   for(int idx=0; idx<tags.get(0).length; idx++){
        					   
        					   byte a1 = tags.get(0)[idx]; 
        					   byte a2 = tags.get(0)[idx]; // this will fail on tri-allelic sites
        					   
        					   for(int t=0; t<tags.size(); t++){
        						   if(tags.get(t)[idx]!=a1){
        							   a2 = tags.get(t)[idx];
        						   }  
        					   } 
        					   
        					   if(a1==a2) continue;  // skip over non-variable sites
        					   
        					   OpenBitSet bitDistB=theTBT.getTaxaReadBitsForTag(0);
        					   for(int z=0; z<bitDistB.size(); z++) { bitDistB.fastClear(z); }
        					   
        					   for(int t=0; t<tags.size(); t++){
        						   if(tags.get(t)[idx]==a1){
        							   OpenBitSet bitDistTmp = theTBT.getTaxaReadBitsForTag(hitIdx[t]);
        							   for(int b=0; b<bitDistA.size(); b++){
        								   if(bitDistTmp.fastGet(b)) bitDistA.fastSet(b);
        							   }
        						   }  
        						   if(tags.get(t)[idx]==a2){
        							   OpenBitSet bitDistTmp = theTBT.getTaxaReadBitsForTag(hitIdx[t]);
        							   for(int b=0; b<bitDistB.size(); b++){
        								   if(bitDistTmp.fastGet(b)) bitDistB.fastSet(b);
        							   }
        						   }
        		   
        					   } 

        		               int cntTaxa = theTBT.getTaxaCount();
        		               int cnt10 = (int)bitDistA.cardinality();
        		               int cnt01 = (int)bitDistB.cardinality();
        		               int cnt11 = (int)OpenBitSet.intersectionCount(bitDistA, bitDistB); //number of lines with both tags
        		               int cntWTag = (int)OpenBitSet.unionCount(bitDistA, bitDistB); // number of lines with one of the two tags
        		               int cnt00 = cntTaxa - cntWTag;
        					   
        					   
        		               double propPresent= (double)cntWTag/(double)cntTaxa;
        		               if(propPresent<(1-maxMissingData)) continue;  // skip if missing too much data
        		               
        		               
        		               
        		               // filter for heterozygousity
        		               if(isDHpopulation && cnt11>0) continue;  // skip if there are any heterozygous calls
        		               //if(cnt11>10) continue; // filter that the two tags don't show up in the same individuals
        		               
        		               
        		               double percentHet = (double)cnt11/(double)cntWTag;  // TODO figure out if divisor should be total taxa or #taxa with one of the tags (currently)
        		               if(percentHet>maxHet) continue; // skip if heterozygousity is > maxHet
        		            
        		               // check if the snp pair has segregation distortion, this basically does the same thing as the minor allele freq
  
        		               //if(isBiparental && (Math.abs(percent10-0.5)>0.3 || Math.abs(percent01-0.5)>0.3)) continue; // skip if one allele is >80% or <20%
        		 
        		               double freqA = (double)cnt10/(double)cntWTag; 
        		               double freqB = (double)cnt01/(double)cntWTag;
        		               if(Math.min(freqA, freqB )<minAlleleFreq) continue;
        		               
        					   

        		               

        		               
          		               double a = (double)cnt10/(double)cntTaxa;
        		               double b = (double)cnt01/(double)cntTaxa; 
        		               //double percent11 = (double)cnt11/(double)cntTaxa;
        		               //System.out.println(cnt10 + "\t" + cnt01 + "\t" + cnt11);
        		               if(cnt11 > a*b*(double)cntTaxa*0.6) continue; // skip stuff that has too much hets
        		               
        		               
        		               // filter for two alleles are not associated (inbred)
        		               //double exp11=(double)cnt01*(double)cnt10/(double)theTBT.getTaxaCount();
        		               //double relExp11=(double)cnt11/exp11;
        		               double p=theFE.getCumlativeP(cnt11, cnt01, cnt10, cnt00);
        		               if(p>pVal) continue;
        		               
        		               
        		               /*
        		               int trials = theTBT.getTaxaCount();
        		               double tProb = a*b;
        		               double p = 1;
        		               
        		               BinomialDistributionImpl x = new BinomialDistributionImpl(trials, tProb);
        		               try {
        		                   p = x.cumulativeProbability(cnt11);
        		               } catch (Exception e) {
        		                   System.err.println("Error in the BinomialDistributionImpl");
        		               }
        		               
        		               if(p>pChiSq) continue;
        		               */
        		               
        		               /*
        		               double pvalChiSq = 1;
        		               double[] expected = {a*(double)cntTaxa, b*(double)cntTaxa, a*b*(double)cntTaxa};
        		               long[] observed = {cnt10-cnt11, cnt01-cnt11, cnt11};
        		               try {
        							pvalChiSq = theChiSq.chiSquareTest(expected, observed);
        						} catch (IllegalArgumentException e) {
        							// TODO Auto-generated catch block
        							e.printStackTrace();
        						} catch (MathException e) {
        							// TODO Auto-generated catch block
        							e.printStackTrace();
        						}
        						if(pvalChiSq>pChiSq) continue;
        		           		*/
        		               
        						nSNPs++;
        						//String snps[] = makeSNPCalls(BaseEncoder.getSequenceFromLong(tagA), BaseEncoder.getSequenceFromLong(tagB));
        						String snps[] = new String[2];
        						snps[1] = BaseEncoder.getSequenceFromLong(tagA);             
        						snps[0] = BaseEncoder.getCharBase(a1) + "/" + BaseEncoder.getCharBase(a2);          
        						            
        						int hitIndex = hitIdx[0];
        						int div = 0;
        						
//        						filteredSnpPair fsp = new filteredSnpPair(snps[1], snps[0], i, hitIndex, div, cnt10, cnt01, cnt11, exp11, relExp11, p, bitDistA, bitDistB, callHets);
//        						fsp.swap();
//        						filSnpList.add(fsp);
        						
        						//"rs	alleles	chrom	pos	strand	assembly	center	protLSID	assayLSID	panelLSID	QCcode	"
        						//TODO need to check if this shoudl be  getTagLength(i)-1. That was outputting 63 bp strings
        						bw.write( BaseEncoder.getSequenceFromLong(tagA).substring(0,theTBT.getTagLength(i)) + "\t" + snps[0] + "\t" + "0" + "\t" + nSNPs + "\t" +"-1"+"\t" + idx + "\tNA\t" + cnt10 +"\t"+ cnt01 +"\t"+ cnt11 +"\t"+ propPresent +"\t" );
        						bw.write(bitsToPseudoSeq(theTBT.getTaxaCount(),bitDistA,bitDistB, snps[0], callHets, numericOutput)+ "\n");   
        		                //System.out.println(BaseEncoder.getSequenceFromLong(tagA));
        		                
        						
        						if(nSNPs%100!=0) continue;
//        		                System.out.print(BaseEncoder.getSequenceFromLong(theTOPM.getTag(hitIndex)));
        		                System.out.println("# SNPs: " + nSNPs);
        						//System.out.printf("con %d %d %d %d %d %d %g %g %g ",i, hitIndex, div, cnt10, cnt01, cnt11, exp11, relExp11, p);
        		                //System.out.print("\n");
        		                System.out.printf("con %d %d %d %d %d %d %g ",i, hitIndex, div, cnt10, cnt01, cnt11, p);
        		                System.out.println();
        		                
        						   
        		                for(int t=0; t<tags.size(); t++){
        							   System.out.print(BaseEncoder.getSequenceFromLong(theTBT.getTag(hitIdx[t])) + "\t");
        							   OpenBitSet bitDist=theTBT.getTaxaReadBitsForTag(hitIdx[t]);
        							   for(int bt=0; bt<bitDist.size(); bt++){ System.out.print(bitDist.getBit(bt) + "\t");}
        							   System.out.println();
        						   } 
        		                System.out.println("\t\t\t\t\t\t\t\t\t" + bitsToPseudoSeq(theTBT.getTaxaCount(),bitDistA,bitDistB, snps[0], callHets, numericOutput));
        		                System.out.println(); 
        		                
        			
        				   }
        			
        			}
        			
        			
        			
        			
        			bw.flush();
        			bw.close();
        			System.out.println("Finished writing to: " + outHapMap);
        			System.out.println("\n# SNPs: " + nSNPs);
        			System.out.println("Divergence: " + divergence);
        			System.out.println("Het: " + maxHet);
        			System.out.println("Missing Data: " + maxMissingData);
        			System.out.println("Min Allele Freq: " + minAlleleFreq);
        			System.out.println("p-value: " + pVal);
        			System.out.println("\n" + "Finished running TagsToSNPsNoAnchor");
        			
        			
        			
        			
        			
        			
        	   }
        	   catch (Exception e) {
        			System.out.println(e.toString());
        	   }
        	   
        	   
        	   
        	   
        	   
        	   
        	   
		   }

		   
    
    /**
     * @param theTBT
     * @param outHapMap
     * @param minDiff
     * @param maxMissingData
     * @param minAlleleFreq
     * @param isDHpopulation, asking if population is double haploid.  If so, will reject any SNPs that are heterozygous, 
     */
    public TagsToSNPsNoAnchor(TagsByTaxa theTBT, String outHapMap, int divergence, double maxMissingData, double minAlleleFreq, double maxHet, boolean callHets, double pVal) {
        
    	System.gc();
    	//do some checks
    	if(maxMissingData>1 || maxMissingData<0){
    		System.out.println("Max Missing Data (maxMissingData) must be between 0 and 1.  Entered Value: " + maxMissingData);
    		System.exit(0);
    	}
    	if(minAlleleFreq>0.5 || minAlleleFreq<0){
    		System.out.println("Minimum Allele Frequency (minAlleleFreq) must be between 0 and 0.5.  Entered Value: " + minAlleleFreq);
    		System.exit(0);
    	}
    	
    	
    	final int cores=Runtime.getRuntime().availableProcessors();
        System.out.println("TagHomologyPhaseMT processors available:"+cores);
        TagMatchFinder theTMF=new TagMatchFinder(theTBT);
        ExecutorService pool = Executors.newCachedThreadPool();
        ThreadPoolExecutor tpe=(ThreadPoolExecutor)pool;
        FisherExact theFE=new FisherExact(theTBT.getTaxaCount()+10);
        ChiSquareTestImpl theChiSq = new ChiSquareTestImpl();

        int nSNPs=0;
       int countSNPs=0, countLoci=0;
       long time=System.currentTimeMillis();
       System.out.println("Max Pool Size "+tpe.getMaximumPoolSize());
       System.out.println("Core Pool Size "+tpe.getCorePoolSize());
       System.out.println("Comp	Tag1	Tag2	Divergence	Tag1Cnt1	Tag2Cnt2	Cnt11	exp11	Ratio11OE	P");
       
       int pauses=0;
       ArrayList<filteredSnpPair> filSnpList = new ArrayList<filteredSnpPair>(0);
       
       //deNovoContigTags tagsOnCtg = new deNovoContigTags("/Users/jpoland/Genome/Barley/morex_TagsFromContigs.bin", true);
       //tagsOnCtg.getTags()[0].getRead();
       //TagsOnContigs theCtgTags = new TagsOnContigs("/Users/jpoland/Genome/Barley/tagsOnContig.bin", true, false);
       //TagMatchFinder theCtgTMF=new TagMatchFinder(theCtgTags);
       
       boolean[] tested = new boolean[theTBT.getTagCount()]; 
       
//	   for (int i = 0; (i < 100000)&&(countSNPs<maxSize); i++) {
	 
			   

        	   try {
        			BufferedWriter bw = new BufferedWriter (new FileWriter(outHapMap), 65536);
        			bw.write("rs	alleles	chrom	pos	strand	assembly	center	protLSID	assayLSID	panelLSID	QCcode	");
        			for (int j = 0; j < theTBT.getTaxaCount()-1; j++) {
        				bw.write(theTBT.getTaxaNames()[j]+ "\t");
        			}
        			bw.write(theTBT.getTaxaNames()[theTBT.getTaxaCount()-1]);
        			bw.newLine();
        			int count = 0;

        			int nTaxa = theTBT.getTaxaCount();
        			OpenBitSet reference = new OpenBitSet(theTBT.getTagCount());
        			
        			for (int i = 0; (i < theTBT.getTagCount())&&(nSNPs<maxSize); i++) {   
        			//	for (int i = 0; i < 100000; i++) { 
        				   long[] tagA=theTBT.getTag(i);
        				   OpenBitSet bitDistA=theTBT.getTaxaReadBitsForTag(i); // the taxa with tagA
        				   
        				   //if(bitDistA.cardinality() < nTaxa*0.3*(1-maxMissingData)) continue;  // skip the rare tags, use tags that show up at reasonable frequency (major allele)
        				   //if(bitDistA.cardinality() < 8) continue;  // skip the rare tags, use tags that show up at reasonable frequency (major allele)
        				   if(reference.get(i)) continue;  // check if this tag has been used as a reference, skip if so
        				   reference.set(i); // set this tags a being a reference
 
        				   TreeMap<Integer,Integer> al=theTMF.findMatchesWithIntLengthWords(tagA, divergence, false);  //TreeMap is a set of Tag pairs(query Tag, hit Tag). There is only 1 mismatch between each pair
        				   if(al.size()<2) continue; // skip stuff that only has one match (self)

         				   int refCnt=0;
         				   for(Entry<Integer,Integer> ht: al.entrySet()) {
        					   if(reference.get(ht.getKey())) refCnt++;
        				   }
        				   if(refCnt>1) continue;  // skip as there was already a reference tag for this TreeMap
        				   
        				   
        				   // populate an array of the tags with hits and their index
        				   ArrayList<byte[]> tags = new ArrayList<byte[]>(0);
        				   int[] hitIdx = new int[al.size()]; 
        				   //int refTagIdx = Integer.MAX_VALUE;
        				   int hi = 0;
        				   for(Entry<Integer,Integer> ht: al.entrySet()) {
        					   hitIdx[hi] = ht.getKey();
        					   hi++;
        					   long[] t = theTBT.getTag(ht.getKey());
        					   tags.add(BaseEncoder.getByteSeqFromLong(t));  
        					   //if(ht.getKey()<refTagIdx){refTagIdx = ht.getKey();} //set the first tag as the reference
        				   }
        				   //if(refTagIdx<i) continue; //skip tags that were previously tested
        				   
        				   

        				   
        				   // test each base in the set of tags for differences
        				   for(int idx=0; idx<tags.get(0).length; idx++){
        					   
        					   byte a1 = tags.get(0)[idx]; 
        					   byte a2 = tags.get(0)[idx]; // this will fail on tri-allelic sites
        					   
        					   for(int t=0; t<tags.size(); t++){
        						   if(tags.get(t)[idx]!=a1){
        							   a2 = tags.get(t)[idx];
        						   }  
        					   } 
        					   
        					   if(a1==a2) continue;  // skip over non-variable sites
        					   
        					   OpenBitSet bitDistB=theTBT.getTaxaReadBitsForTag(0);
        					   for(int z=0; z<bitDistB.size(); z++) { bitDistB.fastClear(z); }
        					   
        					   for(int t=0; t<tags.size(); t++){
        						   if(tags.get(t)[idx]==a1){
        							   OpenBitSet bitDistTmp = theTBT.getTaxaReadBitsForTag(hitIdx[t]);
        							   for(int b=0; b<bitDistA.size(); b++){
        								   if(bitDistTmp.fastGet(b)) bitDistA.fastSet(b);
        							   }
        						   }  
        						   if(tags.get(t)[idx]==a2){
        							   OpenBitSet bitDistTmp = theTBT.getTaxaReadBitsForTag(hitIdx[t]);
        							   for(int b=0; b<bitDistB.size(); b++){
        								   if(bitDistTmp.fastGet(b)) bitDistB.fastSet(b);
        							   }
        						   }
        		   
        					   } 

        		               int cntTaxa = theTBT.getTaxaCount();
        		               int cnt10 = (int)bitDistA.cardinality();
        		               int cnt01 = (int)bitDistB.cardinality();
        		               int cnt11 = (int)OpenBitSet.intersectionCount(bitDistA, bitDistB); //number of lines with both tags
        		               int cntWTag = (int)OpenBitSet.unionCount(bitDistA, bitDistB); // number of lines with one of the two tags
        		               int cnt00 = cntTaxa - cntWTag;
        					   
        					   
        		               double propPresent= (double)cntWTag/(double)cntTaxa;
        		               if(propPresent<(1-maxMissingData)) continue;  // skip if missing too much data
        		               
        		               
        		               
        		               // filter for heterozygousity
        		             
        		               //if(cnt11>10) continue; // filter that the two tags don't show up in the same individuals
        		               
        		               
        		               double percentHet = (double)cnt11/(double)cntWTag;  // TODO figure out if divisor should be total taxa or #taxa with one of the tags (currently)
        		               if(percentHet>maxHet) continue; // skip if heterozygousity is > maxHet
        		            
        		               // check if the snp pair has segregation distortion, this basically does the same thing as the minor allele freq
  
        		               //if(isBiparental && (Math.abs(percent10-0.5)>0.3 || Math.abs(percent01-0.5)>0.3)) continue; // skip if one allele is >80% or <20%
        		 
        		               double freqA = (double)cnt10/(double)cntWTag; 
        		               double freqB = (double)cnt01/(double)cntWTag;
        		               if(Math.min(freqA, freqB )<minAlleleFreq) continue;
        		               
        					   

        		               

        		               
          		               double a = (double)cnt10/(double)cntTaxa;
        		               double b = (double)cnt01/(double)cntTaxa; 
        		               //double percent11 = (double)cnt11/(double)cntTaxa;
        		               //System.out.println(cnt10 + "\t" + cnt01 + "\t" + cnt11);
        		               if(cnt11 > a*b*(double)cntTaxa*0.6) continue; // skip stuff that has too much hets
        		               
        		               
        		               // filter for two alleles are not associated (inbred)
        		               //double exp11=(double)cnt01*(double)cnt10/(double)theTBT.getTaxaCount();
        		               //double relExp11=(double)cnt11/exp11;
        		               double p=theFE.getCumlativeP(cnt11, cnt01, cnt10, cnt00);
        		               if(p>pVal) continue;
        		               
        		               
        		               /*
        		               int trials = theTBT.getTaxaCount();
        		               double tProb = a*b;
        		               double p = 1;
        		               
        		               BinomialDistributionImpl x = new BinomialDistributionImpl(trials, tProb);
        		               try {
        		                   p = x.cumulativeProbability(cnt11);
        		               } catch (Exception e) {
        		                   System.err.println("Error in the BinomialDistributionImpl");
        		               }
        		               
        		               if(p>pChiSq) continue;
        		               */
        		               
        		               /*
        		               double pvalChiSq = 1;
        		               double[] expected = {a*(double)cntTaxa, b*(double)cntTaxa, a*b*(double)cntTaxa};
        		               long[] observed = {cnt10-cnt11, cnt01-cnt11, cnt11};
        		               try {
        							pvalChiSq = theChiSq.chiSquareTest(expected, observed);
        						} catch (IllegalArgumentException e) {
        							// TODO Auto-generated catch block
        							e.printStackTrace();
        						} catch (MathException e) {
        							// TODO Auto-generated catch block
        							e.printStackTrace();
        						}
        						if(pvalChiSq>pChiSq) continue;
        		           		*/
        		               
        						nSNPs++;
        						//String snps[] = makeSNPCalls(BaseEncoder.getSequenceFromLong(tagA), BaseEncoder.getSequenceFromLong(tagB));
        						String snps[] = new String[2];
        						snps[1] = BaseEncoder.getSequenceFromLong(tagA);             
        						snps[0] = BaseEncoder.getCharBase(a1) + "/" + BaseEncoder.getCharBase(a2);          
        						            
        						int hitIndex = hitIdx[0];
        						int div = 0;
        						
//        						filteredSnpPair fsp = new filteredSnpPair(snps[1], snps[0], i, hitIndex, div, cnt10, cnt01, cnt11, exp11, relExp11, p, bitDistA, bitDistB, callHets);
//        						fsp.swap();
//        						filSnpList.add(fsp);
        						
        						//"rs	alleles	chrom	pos	strand	assembly	center	protLSID	assayLSID	panelLSID	QCcode	"
        						//TODO need to check if this shoudl be  getTagLength(i)-1. That was outputting 63 bp strings
        						bw.write( BaseEncoder.getSequenceFromLong(tagA).substring(0,theTBT.getTagLength(i)) + "\t" + snps[0] + "\t" + "0" + "\t" + nSNPs + "\t" +"-1"+"\t" + idx + "\tNA\t" + cnt10 +"\t"+ cnt01 +"\t"+ cnt11 +"\t"+ propPresent +"\t" );
        						bw.write(bitsToPseudoSeq(theTBT.getTaxaCount(),bitDistA,bitDistB, snps[0], callHets, numericOutput)+ "\n");   
        		                //System.out.println(BaseEncoder.getSequenceFromLong(tagA));
        		                
        						
        						if(nSNPs%100!=0) continue;
//        		                System.out.print(BaseEncoder.getSequenceFromLong(theTOPM.getTag(hitIndex)));
        		                System.out.println("# SNPs: " + nSNPs);
        						//System.out.printf("con %d %d %d %d %d %d %g %g %g ",i, hitIndex, div, cnt10, cnt01, cnt11, exp11, relExp11, p);
        		                //System.out.print("\n");
        		                System.out.printf("con %d %d %d %d %d %d %g ",i, hitIndex, div, cnt10, cnt01, cnt11, p);
        		                System.out.println();
        		                
        						   
        		                for(int t=0; t<tags.size(); t++){
        							   System.out.print(BaseEncoder.getSequenceFromLong(theTBT.getTag(hitIdx[t])) + "\t");
        							   OpenBitSet bitDist=theTBT.getTaxaReadBitsForTag(hitIdx[t]);
        							   for(int bt=0; bt<bitDist.size(); bt++){ System.out.print(bitDist.getBit(bt) + "\t");}
        							   System.out.println();
        						   } 
        		                System.out.println("\t\t\t\t\t\t\t\t\t" + bitsToPseudoSeq(theTBT.getTaxaCount(),bitDistA,bitDistB, snps[0], callHets, numericOutput));
        		                System.out.println(); 
        		                
        			
        				   }
        			
        			}
        			
        			
        			
        			
        			bw.flush();
        			bw.close();
        			System.out.println("Finished writing to: " + outHapMap);
        			System.out.println("\n# SNPs: " + nSNPs);
        			System.out.println("Divergence: " + divergence);
        			System.out.println("Het: " + maxHet);
        			System.out.println("Missing Data: " + maxMissingData);
        			System.out.println("Min Allele Freq: " + minAlleleFreq);
        			System.out.println("p-value: " + pVal);
        			System.out.println("\n" + "Finished running TagsToSNPsNoAnchor");
        			
        			
        			
        			
        			
        			
        	   }
        	   catch (Exception e) {
        			System.out.println(e.toString());
        	   }
        	   
        	   
        	   
        	   
        	   
        	   
        	   
		   }

		   
		   /*
		   
           for(Entry<Integer,Integer> ht: al.entrySet()) {
                int hitIndex=ht.getKey();
                //if(hitIndex<i) continue;  //prevent self match 
                if(hitIndex<=i) continue;  //prevent self match and searches only the tags that were not previously searched (>i)
                int div=ht.getValue();
                 
                long[] tagB = theTBT.getTag(hitIndex);
                //System.out.println(BaseEncoder.getSequenceFromLong(tagA));
                //System.out.println(BaseEncoder.getSequenceFromLong(tagB));

                
                //hBitDist is the count of hTag in each Taxa.
                OpenBitSet bitDistB=theTBT.getTaxaReadBitsForTag(hitIndex);
                //hCnt is the number the taxa with presence(1) of hTag
                int cntTaxa = theTBT.getTaxaCount();
                int cnt01 = (int)bitDistB.cardinality();
                int cnt11 = (int)OpenBitSet.intersectionCount(bitDistA, bitDistB); //number of lines with both tags
                int cntWTag = (int)OpenBitSet.unionCount(bitDistA, bitDistB); // number of lines with one of the two tags
                int cnt00 = cntTaxa - cntWTag;
                
                
                double propPresent= (double)cntWTag/(double)cntTaxa;
                if(propPresent<(1-maxMissingData)) continue;  // skip if missing too much data
                
                
                
                // filter for heterozygousity
                if(isDHpopulation && cnt11>0) continue;  // skip if there are any heterozygous calls
                //if(cnt11>10) continue; // filter that the two tags don't show up in the same individuals
                
                
                double percentHet = (double)cnt11/(double)cntWTag;  // TODO figure out if divisor should be total taxa or #taxa with one of the tags (currently)
                if(percentHet>maxHet) continue; // skip if heterozygousity is > maxHet
             
                // check if the snp pair has segregation distortion, this basically does the same thing as the minor allele freq
                double percent10 = (double)cnt10/(double)cntWTag;
                double percent01 = (double)cnt01/(double)cntWTag; 
                double percent11 = (double)cnt11/(double)cntWTag;
                if(isBiparental && (Math.abs(percent10-0.5)>0.3 || Math.abs(percent01-0.5)>0.3)) continue; // skip if one allele is >80% or <20%
  
                
                if(Math.min(percent10, percent01)<minAlleleFreq) continue;
                
                
                //int cnt11=(int)OpenBitSet.intersectionCount(cBitDist, hBitDist);
                
                //caculate the contingency table
                // determine if tag pairs are not independent, i.e on the same chromosome, this will only work for inbreds
                //int cnt01=cntA-cnt11;
                //int cnt10=cntB-cnt11;
                //int cnt00=theTBT.getTaxaCount()-cnt11-cnt01-cnt10;
                
                //if(cnt10!=0 && cnt01!=0) continue; // skip if the two tags are 
                
                // filter for two alleles are not associated (inbred)
                double exp11=(double)cnt01*(double)cnt10/(double)theTBT.getTaxaCount();
                double relExp11=(double)cnt11/exp11;
                double p=theFE.getCumlativeP(cnt11, cnt01, cnt10, cnt00);
                //if(p>0.01) continue;
                
                double pvalChiSq = 1;
                double[] expected = {percent10*(double)cntTaxa, percent01*(double)cntTaxa, percent10*percent01*(double)cntTaxa};
                long[] observed = {cnt10, cnt01, cnt11};
                try {
					pvalChiSq = theChiSq.chiSquareTest(expected, observed);
				} catch (IllegalArgumentException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (MathException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				if(pvalChiSq>pChiSq) continue;
            

                //set some filters
                
                //if(propPresent>1-(percentMissingData/100) && p<0.05) {
              
                //if((p<0.01)&&(propPresent>0.1)) {
                //if((p<0.01)&&(relExp11<0.1)&&(propPresent>0.1)) {
                //if((cnt11==0 && cnt00==0 && cnt01+cnt10>10) || (cnt11+cnt00>10 && cnt01==0 && cnt10==0) ) {
					String snps[] = makeSNPCalls(BaseEncoder.getSequenceFromLong(tagA), BaseEncoder.getSequenceFromLong(tagB));
                	filteredSnpPair fsp = new filteredSnpPair(snps[1], snps[0], i, hitIndex, div, cnt10, cnt01, cnt11, exp11, relExp11, p, bitDistA, bitDistB, callHets);
					fsp.swap();
					filSnpList.add(fsp);
                    System.out.println(BaseEncoder.getSequenceFromLong(tagA));
                    System.out.println(BaseEncoder.getSequenceFromLong(tagB));
//                    System.out.print(BaseEncoder.getSequenceFromLong(theTOPM.getTag(hitIndex)));
                    System.out.printf("con %d %d %d %d %d %d %g %g %g ",i, hitIndex, div, cnt10, cnt01, cnt11, exp11, relExp11, p);
                    System.out.print(" ");
                    System.out.println(bitsToPseudoSeq(theTBT.getTaxaCount(),bitDistA,bitDistB, snps[0], callHets));
                
           }
           
           */
       
       
	   
	         
	           
	           
	    /*       
	           
	   filteredSnpPair[] filSnps = filSnpList.toArray(new filteredSnpPair[filSnpList.size()]);
	   
	   System.out.println(filSnpList.size());
	   System.out.println("SNPs: " + filSnps.length);
	   Arrays.sort(filSnps);
	   System.out.println(filSnps.length);
	   try {
			BufferedWriter bw = new BufferedWriter (new FileWriter(outHapMap), 65536);
			bw.write("rs	alleles	chrom	pos	strand	assembly	center	protLSID	assayLSID	panelLSID	QCcode	");
			for (int j = 0; j < theTBT.getTaxaCount()-1; j++) {
				bw.write(theTBT.getTaxaNames()[j]+ "\t");
			}
			bw.write(theTBT.getTaxaNames()[theTBT.getTaxaCount()-1]);
			bw.newLine();
			int count = 0;
			for (int j = 0; j < filSnps.length; j++) {
				count++;
				bw.write(filSnps[j].mkStr(filSnps[j].sequence, 0, count, theTBT.getTaxaCount(), callHets));
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
    
    */
    
    


    /**
     * @param theTBT
     * @param outHapMap
     * @param minDiff
     * @param maxMissingData
     * @param minAlleleFreq
     * @param isDHpopulation, asking if population is double haploid.  If so, will reject any SNPs that are heterozygous, 
     */
    public TagsToSNPsNoAnchor(TagsByTaxa theTBT, TagsByTaxa theTBT2, String outHapMap, int divergence, double maxMissingData, double minAlleleFreq, double maxHet, boolean callHets, double pVal) {
        
    	boolean numericOut = false;
    	System.gc();
    	//do some checks
    	if(maxMissingData>1 || maxMissingData<0){
    		System.out.println("Max Missing Data (maxMissingData) must be between 0 and 1.  Entered Value: " + maxMissingData);
    		System.exit(0);
    	}
    	if(minAlleleFreq>0.5 || minAlleleFreq<0){
    		System.out.println("Minimum Allele Frequency (minAlleleFreq) must be between 0 and 0.5.  Entered Value: " + minAlleleFreq);
    		System.exit(0);
    	}
    	
    	
    	final int cores=Runtime.getRuntime().availableProcessors();
        System.out.println("TagHomologyPhaseMT processors available:"+cores);
        TagMatchFinder theTMF=new TagMatchFinder(theTBT);
        ExecutorService pool = Executors.newCachedThreadPool();
        ThreadPoolExecutor tpe=(ThreadPoolExecutor)pool;
        FisherExact theFE=new FisherExact(theTBT.getTaxaCount()+10);
        ChiSquareTestImpl theChiSq = new ChiSquareTestImpl();

        int nSNPs=0;
       int countSNPs=0, countLoci=0;
       long time=System.currentTimeMillis();
       System.out.println("Max Pool Size "+tpe.getMaximumPoolSize());
       System.out.println("Core Pool Size "+tpe.getCorePoolSize());
       System.out.println("Comp	Tag1	Tag2	Divergence	Tag1Cnt1	Tag2Cnt2	Cnt11	exp11	Ratio11OE	P");
       
       int pauses=0;
       ArrayList<filteredSnpPair> filSnpList = new ArrayList<filteredSnpPair>(0);
       
       //deNovoContigTags tagsOnCtg = new deNovoContigTags("/Users/jpoland/Genome/Barley/morex_TagsFromContigs.bin", true);
       //tagsOnCtg.getTags()[0].getRead();
       //TagsOnContigs theCtgTags = new TagsOnContigs("/Users/jpoland/Genome/Barley/tagsOnContig.bin", true, false);
       //TagMatchFinder theCtgTMF=new TagMatchFinder(theCtgTags);
       
       boolean[] tested = new boolean[theTBT.getTagCount()]; 
       
//	   for (int i = 0; (i < 100000)&&(countSNPs<maxSize); i++) {
	 
			   

        	   try {
        			BufferedWriter bw = new BufferedWriter (new FileWriter(outHapMap), 65536);
        			bw.write("rs	alleles	chrom	pos	strand	assembly	center	protLSID	assayLSID	panelLSID	QCcode");
        			//for (int j = 0; j < theTBT.getTaxaCount(); j++) {
        			//	bw.write("\t" + theTBT.getTaxaNames()[j]);
        			//}
        			for (int j = 0; j < theTBT2.getTaxaCount(); j++) {
        				bw.write("\t" + theTBT2.getTaxaNames()[j]);
        			}
        			
        			//bw.write(theTBT.getTaxaNames()[theTBT.getTaxaCount()-1]);
        			//bw.write(theTBT2.getTaxaNames()[theTBT2.getTaxaCount()-1]); // write the second tbt set
        			bw.newLine();
        			int count = 0;
        			
        			OpenBitSet reference = new OpenBitSet(theTBT.getTagCount());

        			for (int i = 0; (i < theTBT.getTagCount())&&(nSNPs<maxSize); i++) {   
        				/*   
        				//cTag is the a query Tag
        				   long[] tagA=theTBT.getTag(i);
        				   //TreeMap is a set of Tag pairs(query Tag, hit Tag). There is only 1 mismatch between each pair
        				   TreeMap<Integer,Integer> al=theTMF.findMatchesWithIntLengthWords(tagA, divergence, false);
        				   if(al.size()<2) continue;
        				   //cBitDist is the count of the cTag in each Taxa. Here, we use the Bit Verision. It means presence(1) and absence(0). Just like 0	1	0  1 1 1 0 0
        				   OpenBitSet bitDistA=theTBT.getTaxaReadBitsForTag(i);
        				   //cCnt is the number the taxa with presence(1) of cTag
        				   
        		           
        				   //Here, we try to Fisher test each Tag pairs with mismatch = 1
        				   
        				   //Byte[][] tags = new Byte[al.size()][64];
        				  
        				   ArrayList<byte[]> tags = new ArrayList<byte[]>(0);
        				   int[] hitIdx = new int[al.size()]; 
        				   int refTagIdx = Integer.MAX_VALUE;
        				   int hi = 0;
        				   for(Entry<Integer,Integer> ht: al.entrySet()) {
        					   hitIdx[hi] = ht.getKey();
        					   hi++;
        					   long[] t = theTBT.getTag(ht.getKey());
        					   tags.add(BaseEncoder.getByteSeqFromLong(t));  
        					   if(ht.getKey()<refTagIdx){refTagIdx = ht.getKey();} //set the first tag as the reference
        				   }
        				   if(refTagIdx<i) continue; //skip tags that were previously tested
        				   */
        				   
        				   long[] tagA=theTBT.getTag(i);
        				   OpenBitSet bitDistA=theTBT.getTaxaReadBitsForTag(i); // the taxa with tagA
        				   
        				   //if(bitDistA.cardinality() < nTaxa*0.3*(1-maxMissingData)) continue;  // skip the rare tags, use tags that show up at reasonable frequency (major allele)
        				   //if(bitDistA.cardinality() < 8) continue;  // skip the rare tags, use tags that show up at reasonable frequency (major allele)
        				   if(reference.get(i)) continue;  // check if this tag has been used as a reference, skip if so
        				   reference.set(i); // set this tags a being a reference
 
        				   TreeMap<Integer,Integer> al=theTMF.findMatchesWithIntLengthWords(tagA, divergence, false);  //TreeMap is a set of Tag pairs(query Tag, hit Tag). There is only 1 mismatch between each pair
        				   if(al.size()<2) continue; // skip stuff that only has one match (self)

         				   int refCnt=0;
         				   for(Entry<Integer,Integer> ht: al.entrySet()) {
        					   if(reference.get(ht.getKey())) refCnt++;
        				   }
        				   if(refCnt>1) continue;  // skip as there was already a reference tag for this TreeMap
        				   
        				   
        				   // populate an array of the tags with hits and their index
        				   ArrayList<byte[]> tags = new ArrayList<byte[]>(0);
        				   int[] hitIdx = new int[al.size()]; 
        				   //int refTagIdx = Integer.MAX_VALUE;
        				   int hi = 0;
        				   for(Entry<Integer,Integer> ht: al.entrySet()) {
        					   hitIdx[hi] = ht.getKey();
        					   hi++;
        					   long[] t = theTBT.getTag(ht.getKey());
        					   tags.add(BaseEncoder.getByteSeqFromLong(t));  
        					   //if(ht.getKey()<refTagIdx){refTagIdx = ht.getKey();} //set the first tag as the reference
        				   }
        				   //if(refTagIdx<i) continue; //skip tags that were previously tested
        				   
        				   
        				   // test each base in the set of tags for differences
        				   for(int idx=0; idx<tags.get(0).length; idx++){
        					   
        					   byte a1 = tags.get(0)[idx]; 
        					   byte a2 = tags.get(0)[idx]; // this will fail on tri-allelic sites
        					   
        					   for(int t=0; t<tags.size(); t++){
        						   if(tags.get(t)[idx]!=a1){
        							   a2 = tags.get(t)[idx];
        						   }  
        					   } 
        					   
        					   if(a1==a2) continue;  // skip over non-variable sites
        					   
        					   OpenBitSet bitDistB=theTBT.getTaxaReadBitsForTag(0);
        					   for(int z=0; z<bitDistB.size(); z++) { bitDistB.fastClear(z); }
        					   
        					   // make bit sets for the set of selection candidates
        					   OpenBitSet bitDistA2=theTBT2.getTaxaReadBitsForTag(0);
        					   for(int z=0; z<bitDistA2.size(); z++) { bitDistA2.fastClear(z); }
        					   
        					   OpenBitSet bitDistB2=theTBT2.getTaxaReadBitsForTag(0);
        					   for(int z=0; z<bitDistB2.size(); z++) { bitDistB2.fastClear(z); }
        					   
        					   if(theTBT2.getTagIndex(tagA)>-1){
        						   bitDistA2 = theTBT2.getTaxaReadBitsForTag(theTBT2.getTagIndex(tagA));
        					   }
        					   
        					   for(int t=0; t<tags.size(); t++){
        						   if(tags.get(t)[idx]==a1){
        							   OpenBitSet bitDistTmp = theTBT.getTaxaReadBitsForTag(hitIdx[t]);
        							   for(int b=0; b<bitDistA.size(); b++){
        								   if(bitDistTmp.fastGet(b)) bitDistA.fastSet(b);
        							   }
        							   // look throught the second tbt (selection candidates) and make calls
        							   int selIdx = theTBT2.getTagIndex(theTBT.getTag(hitIdx[t]));
                					   if(selIdx>-1){
                						   OpenBitSet bitDistTmp2 = theTBT2.getTaxaReadBitsForTag(selIdx);
            							   for(int b=0; b<bitDistA2.size(); b++){
            								   if(bitDistTmp2.fastGet(b)) bitDistA2.fastSet(b);
            							   }
                					   }
        							   
        						   }  
        						   if(tags.get(t)[idx]==a2){
        							   OpenBitSet bitDistTmp = theTBT.getTaxaReadBitsForTag(hitIdx[t]);
        							   for(int b=0; b<bitDistB.size(); b++){
        								   if(bitDistTmp.fastGet(b)) bitDistB.fastSet(b);
        							   }
        							   // look throught the second tbt (selection candidates) and make calls
        							   int selIdx = theTBT2.getTagIndex(theTBT.getTag(hitIdx[t]));
                					   if(selIdx>-1){
                						   OpenBitSet bitDistTmp2 = theTBT2.getTaxaReadBitsForTag(selIdx);
            							   for(int b=0; b<bitDistB2.size(); b++){
            								   if(bitDistTmp2.fastGet(b)) bitDistB2.fastSet(b);
            							   }
                					   }
        						   }
        		   
        					   } 

        		               int cntTaxa = theTBT.getTaxaCount();
        		               int cnt10=(int)bitDistA.cardinality();
        		               int cnt01 = (int)bitDistB.cardinality();
        		               int cnt11 = (int)OpenBitSet.intersectionCount(bitDistA, bitDistB); //number of lines with both tags
        		               int cntWTag = (int)OpenBitSet.unionCount(bitDistA, bitDistB); // number of lines with one of the two tags
        		               int cnt00 = cntTaxa - cntWTag;
        					   
        					   
        		               double propPresent= (double)cntWTag/(double)cntTaxa;
        		               if(propPresent<(1-maxMissingData)) continue;  // skip if missing too much data
        		               
        		               
        		               
        		               // filter for heterozygousity
        		               //if(isDHpopulation && cnt11>0) continue;  // skip if there are any heterozygous calls
        		               //if(cnt11>10) continue; // filter that the two tags don't show up in the same individuals
        		               
        		               
        		               double percentHet = (double)cnt11/(double)cntWTag;  // TODO figure out if divisor should be total taxa or #taxa with one of the tags (currently)
        		               if(percentHet>maxHet) continue; // skip if heterozygousity is > maxHet
        		            
        		               // check if the snp pair has segregation distortion, this basically does the same thing as the minor allele freq
  
        		               //if(isBiparental && (Math.abs(percent10-0.5)>0.3 || Math.abs(percent01-0.5)>0.3)) continue; // skip if one allele is >80% or <20%
        		 
        		               double freqA = (double)cnt10/(double)cntWTag; 
        		               double freqB = (double)cnt01/(double)cntWTag;
        		               if(Math.min(freqA, freqB )<minAlleleFreq) continue;
        		               
        					   

        		               

        		               
          		               double a = (double)cnt10/(double)cntTaxa;
        		               double b = (double)cnt01/(double)cntTaxa; 
        		               //double percent11 = (double)cnt11/(double)cntTaxa;
        		               //System.out.println(cnt10 + "\t" + cnt01 + "\t" + cnt11);
        		               if(cnt11 > a*b*(double)cntTaxa*0.6) continue; // skip stuff that has too much hets
        		               
        		               
        		               // filter for two alleles are not associated (inbred)
        		               //double exp11=(double)cnt01*(double)cnt10/(double)theTBT.getTaxaCount();
        		               //double relExp11=(double)cnt11/exp11;
        		               double p=theFE.getCumlativeP(cnt11, cnt01, cnt10, cnt00);
        		               if(p>pVal) continue;
        		               
        		               
        		               /*
        		               int trials = theTBT.getTaxaCount();
        		               double tProb = a*b;
        		               double p = 1;
        		               
        		               BinomialDistributionImpl x = new BinomialDistributionImpl(trials, tProb);
        		               try {
        		                   p = x.cumulativeProbability(cnt11);
        		               } catch (Exception e) {
        		                   System.err.println("Error in the BinomialDistributionImpl");
        		               }
        		               
        		               if(p>pChiSq) continue;
        		               */
        		               
        		               /*
        		               double pvalChiSq = 1;
        		               double[] expected = {a*(double)cntTaxa, b*(double)cntTaxa, a*b*(double)cntTaxa};
        		               long[] observed = {cnt10-cnt11, cnt01-cnt11, cnt11};
        		               try {
        							pvalChiSq = theChiSq.chiSquareTest(expected, observed);
        						} catch (IllegalArgumentException e) {
        							// TODO Auto-generated catch block
        							e.printStackTrace();
        						} catch (MathException e) {
        							// TODO Auto-generated catch block
        							e.printStackTrace();
        						}
        						if(pvalChiSq>pChiSq) continue;
        		           		*/
        		               
        						nSNPs++;
        						//String snps[] = makeSNPCalls(BaseEncoder.getSequenceFromLong(tagA), BaseEncoder.getSequenceFromLong(tagB));
        						String snps[] = new String[2];
        						snps[1] = BaseEncoder.getSequenceFromLong(tagA);             
        						snps[0] = BaseEncoder.getCharBase(a1) + "/" + BaseEncoder.getCharBase(a2);          
        						            
        						int hitIndex = hitIdx[0];
        						int div = 0;
        						
//        						filteredSnpPair fsp = new filteredSnpPair(snps[1], snps[0], i, hitIndex, div, cnt10, cnt01, cnt11, exp11, relExp11, p, bitDistA, bitDistB, callHets);
//        						fsp.swap();
//        						filSnpList.add(fsp);
        						
        						//"rs	alleles	chrom	pos	strand	assembly	center	protLSID	assayLSID	panelLSID	QCcode	"
        						bw.write( BaseEncoder.getSequenceFromLong(tagA).substring(0,theTBT.getTagLength(i)-1) + "\t" + snps[0] + "\t" + "0" + "\t" + nSNPs + "\t" +"-1"+"\t" + idx + "\tNA\t" + cnt10 +"\t"+ cnt01 +"\t"+ cnt11 +"\t"+ propPresent +"\t" );
        						//bw.write(bitsToPseudoSeq(theTBT.getTaxaCount(),bitDistA,bitDistB, snps[0], callHets, numericOut));   
        		                //System.out.println(BaseEncoder.getSequenceFromLong(tagA));
        		                // write the second set
        						bw.write(bitsToPseudoSeq(theTBT2.getTaxaCount(),bitDistA2,bitDistB2, snps[0], callHets, numericOut)+ "\n");   
        						//bw.newLine();
        						
        						
        						
        						if(nSNPs%1000!=0) continue;
//        		                System.out.print(BaseEncoder.getSequenceFromLong(theTOPM.getTag(hitIndex)));
        		                System.out.println("# SNPs: " + nSNPs);
        						//System.out.printf("con %d %d %d %d %d %d %g %g %g ",i, hitIndex, div, cnt10, cnt01, cnt11, exp11, relExp11, p);
        		                //System.out.print("\n");
        		                System.out.printf("con %d %d %d %d %d %d %g ",i, hitIndex, div, cnt10, cnt01, cnt11, p);
        		                //System.out.println("\t\t\t\t\t\t" + bitsToPseudoSeq(theTBT.getTaxaCount(),bitDistA,bitDistB, snps[0], callHets, numericOut));
        		                System.out.println("\t\t\t\t\t\t" + bitsToPseudoSeq(theTBT2.getTaxaCount(),bitDistA2,bitDistB2, snps[0], callHets, numericOut));  
        		                
        		                for(int t=0; t<tags.size(); t++){
        							   System.out.print(BaseEncoder.getSequenceFromLong(theTBT.getTag(hitIdx[t])) + "\t");
        							   OpenBitSet bitDist=theTBT.getTaxaReadBitsForTag(hitIdx[t]);
        							   for(int bt=0; bt<bitDist.size(); bt++){ System.out.print(bitDist.getBit(bt) + "\t");}
        							   System.out.println();
        						   } 
        		                System.out.println(); 

        			
        				   }
        			
        			}
        			
        			
        			
        			
        			bw.flush();
        			bw.close();
        			System.out.println("Finished writing to: " + outHapMap);
        			System.out.println("\n# SNPs: " + nSNPs);
        			System.out.println("Divergence: " + divergence);
        			System.out.println("Het: " + maxHet);
        			System.out.println("Missing Data: " + maxMissingData);
        			System.out.println("Min Allele Freq: " + minAlleleFreq);
        			System.out.println("p-value: " + pVal);
        			System.out.println("\n" + "Finished running TagsToSNPsNoAnchor");
        			
        			
        			
        			
        			
        			
        	   }
        	   catch (Exception e) {
        			System.out.println(e.toString());
        	   }
        	   
        	   
        	   
        	   
        	   
        	   
        	   
		   }

	
    
    
    
    
    
	public String[] makeSNPCalls(String allele1, String allele2){
		
		String[] snpString = new String[2];
		snpString[0]="";
		snpString[1]="";
		
		//for(int i=0; i<snps.length; i++){
			for(int a=0; a<allele1.length(); a++){
				//char c1 = allele1.charAt(a);
				//char c2 = allele2.charAt(a);
				//System.out.println(c1 == c2);
				if(allele1.charAt(a) == allele2.charAt(a)){
					snpString[1] = snpString[1] + allele1.charAt(a);
				}
				if(!(allele1.charAt(a) == allele2.charAt(a))){
					/*
					if(allele1.charAt(a) < allele2.charAt(a)){
						snpString[0] = allele1.charAt(a) + "/" + allele2.charAt(a);
						snpString[1] = snpString[1] + "[" + allele1.charAt(a) + "/" + allele2.charAt(a) + "]";
					}
					else{
						snpString[0] = allele2.charAt(a) + "/" + allele1.charAt(a);
						snpString[1] = snpString[1] + "[" + allele2.charAt(a) + "/" + allele1.charAt(a) + "]";
					}
					*/
					
					snpString[0] = allele1.charAt(a) + "/" + allele2.charAt(a);
					snpString[1] = snpString[1] + "[" + allele1.charAt(a) + "/" + allele2.charAt(a) + "]";
				}
				
				//if(!(allele1.charAt(a) == allele2.charAt(a))){
				//	snpString[0] = allele1.charAt(a) + "/" + allele2.charAt(a);
				//	snpString[1] = allele1.substring(0,a) + "[" + allele1.charAt(a) + "/" + allele2.charAt(a) + "]" + allele1.substring(a+1);
				//	break;
				//}
				
			}
		//}
		return snpString;
	}
    
    
    
    
    
	class filteredSnpPair implements Comparable <filteredSnpPair> {
		String sequence; //store sequence of SNP marker
		String snp;
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
		
		filteredSnpPair(String sequence, String snp, int queryIndex, int hitIndex, int div, int cCnt, int hCnt, int cnt11, double exp11, double relExp11, double p, OpenBitSet cBitDist,  OpenBitSet hBitDist, boolean callHets) {
			this.sequence = sequence;
			this.snp = snp;
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
		
		public String mkStr (String name, int chr, int posi, int taxaCount, boolean callHets, boolean numericOut) {
			StringBuilder sb = new StringBuilder();
			sb.append(name).append("\t").append(snp).append("\t").append(chr).append("\t").append(posi).append("\t+\tNA\tSWGDiv\tGBS\tSWGV1\tSWGPop\tQC+\t");
			sb.append(bitsToPseudoSeq(taxaCount,cBitDist,hBitDist, snp, callHets, numericOut));
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
	
	
	
	
	public static ArrayList<filteredSnpPair> collapseSeqError(TagsByTaxa theTBT, ArrayList<filteredSnpPair> filSnpList){
		
		TagMatchFinder theTMF=new TagMatchFinder(theTBT);
		
		int i = 0;
		while(i<filSnpList.size()){
			int idx = filSnpList.get(i).queryIndex;
			long[] tagA=theTBT.getTag(idx);
			TreeMap<Integer,Integer> al=theTMF.findMatchesWithIntLengthWords(tagA, 2, false);
			
			
			
			OpenBitSet bitDistA=theTBT.getTaxaReadBitsForTag(i);
			
			
		}
		
		
		
		
		
		return filSnpList;
	}
	
	
    public static String bitsToPseudoSeq(int numTaxa, OpenBitSet t1, OpenBitSet t2, String snp, boolean callHets, boolean numericOut) {
        StringBuilder sb=new StringBuilder();
        //String[] allele={"A","C","R"};
        String[] tmp = snp.split("/");
        String[] allele = {tmp[0], tmp[1], "N"};
        if(callHets) allele[2] = "H";
        if(numericOut) allele = new String[]{"-1", "1", "0"};
        for (int i = 0; i < numTaxa; i++) {
            int a=t1.fastGet(i)?1:0;
            a+=t2.fastGet(i)?2:0;
            if(a==0) {
            	if(!numericOut) sb.append(DataType.UNKNOWN_CHARACTER);
            	else sb.append("-9");
            		}
            //if(a==0) {sb.append("-9");}
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
        System.out.println("Starting TagsToSNPsNoAnchor");


		TagsByTaxa theTBT = new TagsByTaxaBitFileMap("/Users/jpoland/gbs/barley/mb/tbt_merge_20111024.bin");
		

		String outHapMap = "/Users/jpoland/gbs/barley/mb/MB_HapMap_Dif3_Na60_20111121.hap";
		
        System.out.println("Starting TagsToSNPByAlignmentMT");


	    int nDiff = 3;
	    double maxMissingData = 0.6;
	    double minorAlleleFreq = 0.2;
	    double maxHet = 0.03;
	    boolean isDHpopulation = false;
	    boolean isBiparental = false;
	    boolean callHets = true;
	    boolean numericOut = false;
	    double chiSqPVal = 0.001;
	         
		new TagsToSNPsNoAnchor(theTBT, outHapMap, nDiff, maxMissingData, minorAlleleFreq, maxHet, isDHpopulation, isBiparental, callHets, chiSqPVal);
		System.gc();
		
		
		
		
		System.out.println("\n" + "DONE!");
	
       
    }


}
