package net.maizegenetics.genome.GBS;

import net.maizegenetics.genome.GBS.ReadCounts;
import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map.Entry;
import java.util.TreeMap;
import net.maizegenetics.genome.BaseEncoder;
import net.maizegenetics.genome.PairwiseAlignment;


/**
 * Created by IntelliJ IDEA.
 * User: ed
 * Date: Jan 29, 2008
 * Time: 1:05:23 PM
 * To change this template use File | Settings | File Templates.
 */
public class TagMatchFinder  {
    private static final int chunkSize=BaseEncoder.chunkSize;
   // static long[][] genomeTags;
    private long[][] refTags;
 //   private long query;
 //   private int genomePosition;

    private static int[][] lookups;
    public static final int wordLength=16;        //best results 8 for sensitivity
    public static final int maxDivergence=15;

//    int[] refTagNumbers=null;  //used if snp synonym table is used and need to return
//    boolean useSNPSynonyms=false;
//    int matches=0, synmatches=0;
    private long sTime, cTime, aTime;

    private PairwiseAlignment pa;
    private Thread liveThread=null; //use to keep this alive while another thread is alive


    public int[][] getTagLookTable(){
        return lookups;
    }


    public TagMatchFinder(long[][] refTags) {
        this.refTags=refTags;
       // genomeTags=refTags.haplotype;
        init();
         pa=new PairwiseAlignment(64,32);
         long[] query={refTags[0][10],refTags[1][10]};
         TreeMap<Integer,Integer> al=findMatchesWithIntLengthWords(query, 1000,false);
         System.out.println(al.toString());

    }

    public static void main(String[] args) {
       ReadCounts refTags=new ReadCounts("C:/EdStuff/Solexa/NGG_IBM/counts/countMerge091223_20.txt",true);
       TagMatchFinder tmf=new TagMatchFinder(refTags.haplotype);

    }

    private void init() {
        lookups=new int[2][2*refTags.length*refTags[0].length];  //4 x 16base words
        int count=0;
        for(int i=0; i<refTags[0].length; i++) {
            for(int j=0; j<refTags.length; j++) {
                int word[]=BaseEncoder.getIntFromLong(refTags[j][i]);
                lookups[0][count]=word[0];
                lookups[1][count]=i;
                count++;
                lookups[0][count]=word[1];
                lookups[1][count]=i;
                count++;
            //    System.out.println(BaseEncoder.getSequenceFromLong(refTags[0][i]));
//                if(word[0]==-2147483648) {System.out.println(BaseEncoder.getSequenceFromLong(refTags[0][i])+BaseEncoder.getSequenceFromLong(refTags[1][i]));
//                System.out.println(BaseEncoder.getSequenceFromInt(word[0])+" "+BaseEncoder.getSequenceFromInt(word[1]));
//                }
//                if(word[0]==0) {System.out.println(BaseEncoder.getSequenceFromLong(refTags[0][i])+BaseEncoder.getSequenceFromLong(refTags[1][i]));
//                System.out.println(BaseEncoder.getSequenceFromInt(word[0])+" "+BaseEncoder.getSequenceFromInt(word[1]));
//                }
            }
        }
       GenericSorting.quickSort(0, lookups[0].length, comp, swapper);
       reduceDuplicates();
       pa=new PairwiseAlignment(32,32);
//          System.out.println("shiftLookup size="+shiftLookup.length+"  shiftQuery size="+shiftQuery.length);
    }

  private void reduceDuplicates() {
      int start=0, end=-1, currHap=lookups[0][0], duplicated=0;
      System.out.println(BaseEncoder.getSequenceFromInt(currHap)+" "+currHap);
      for (int i = 0; i < lookups[0].length; i++) {
          if(lookups[0][i]==currHap) {end=i;}
          else {
            if(((end-start)>1000)) {
                //System.out.println(BaseEncoder.getSequenceFromInt(currHap)+" "+(end-start));
                for(int j=start; j<=end; j++) {
                    lookups[0][j]=Integer.MAX_VALUE;
                    duplicated++;
                }
            }
            currHap=lookups[0][i];
            start=end=i;
          }
      }
      GenericSorting.quickSort(0, lookups[0].length, comp, swapper);

      int[][] newlookups=new int[2][lookups[0].length-duplicated];
      System.arraycopy(lookups[0], 0, newlookups[0], 0, lookups[0].length-duplicated);
      System.arraycopy(lookups[1], 0, newlookups[1], 0, lookups[1].length-duplicated);
      System.out.println("Old Lookup Size:"+lookups[0].length+"  new size:"+newlookups[0].length);
      lookups=newlookups;
  }

  /**
   *
   * @param query
   * @param maxDiv
   * @param keepOnlyBest
   * @return
   */
  public TreeMap<Integer,Integer> findMatchesWithIntLengthWords(long[] query, int maxDiv, boolean keepOnlyBest) {
      TreeMap<Integer,Integer> hitsAndDiv=new TreeMap<Integer,Integer>();
      int bestDiv=maxDiv;
      for (int i = 0; i < query.length; i++) {
          int words[]=BaseEncoder.getIntFromLong(query[i]);
          for(int word: words) {
      //      System.out.println("word"+word+" string:"+BaseEncoder.getSequenceFromInt(word));
            int hit=Arrays.binarySearch(lookups[0], word);
            if(hit<0) continue;  //no hit for this word
            while((hit>0)&&(word==lookups[0][hit-1])) {hit--;}  //backup to the first hit
            while((hit<lookups[0].length)&&(lookups[0][hit]==word)) {
                int indexPossHapHit=lookups[1][hit];
                if(!hitsAndDiv.containsKey(indexPossHapHit)) {
                    int div=BaseEncoder.seqDifferences(refTags[0][indexPossHapHit], query[0], maxDiv)+
                        BaseEncoder.seqDifferences(refTags[1][indexPossHapHit], query[1], maxDiv);
                    if(div<=maxDiv) {hitsAndDiv.put(indexPossHapHit, div); //System.out.println("Hit"+indexPossHapHit+" div:"+div);
                    }
                    if(div<bestDiv) bestDiv=div;
                }
                hit++;
            }
          }
      }
      if(keepOnlyBest&&hitsAndDiv.size()>1) {
          TreeMap<Integer,Integer> bestHitsAndDiv=new TreeMap<Integer,Integer>();
          for(Entry<Integer,Integer> e: hitsAndDiv.entrySet()) {
              if(e.getValue()==bestDiv) bestHitsAndDiv.put(e.getKey(), e.getValue());
          }
          return bestHitsAndDiv;
      } else return hitsAndDiv;
      //1bp indels can be dealt with by shift
  }

     Swapper swapper = new Swapper() {
   public void swap(int a, int b) {
      int t1, t2;
      t1 = lookups[0][a]; lookups[0][a] = lookups[0][b];	lookups[0][b] = t1;
      t2 = lookups[1][a]; lookups[1][a] = lookups[1][b]; lookups[1][b] = t2;
   }
};

IntComparator comp = new IntComparator() {
   public int compare(int a, int b) {
      if (lookups[0][a]==lookups[0][b]) return lookups[1][a]==lookups[1][b] ? 0 : (lookups[1][a]<lookups[1][b] ? -1 : 1);
      return lookups[0][a]<lookups[0][b] ? -1 : 1;
   }
};
}
