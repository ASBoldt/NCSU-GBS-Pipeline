package net.maizegenetics.genome.GBS;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Arrays;
import net.maizegenetics.genome.BaseEncoder;

/**
 * Digest the genome and the creates output file of potential haplotypes of 64bp in
 * length.  Includes location, direction, and size of fragments.
 * User: ed
 */
public class VirtualDigester {
    private enum Ends {Left, Both, Right};
    private String[] cuts={"GCAGC","GCTGC"}; // TseI or ApeKI

    private int reOffset1=1, reOffset2=4;  // normally 1 for CCGG
    private int maxTags=9000000+1;
    private static int chunkSize=BaseEncoder.chunkSize;
    private String polyA="AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    private String polyT="TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";


    public VirtualDigester(File inputSeqFile, File outputSeqFile) {
        System.out.println("-1="+BaseEncoder.getSequenceFromLong(-1));
        System.out.println("Max Value="+BaseEncoder.getSequenceFromLong(Long.MAX_VALUE));
        System.out.println("MIN Value="+BaseEncoder.getSequenceFromLong(Long.MIN_VALUE));
        ReadsWPhysicalMap theSHM=new ReadsWPhysicalMap(maxTags);
        readGenome(inputSeqFile, theSHM);
        theSHM.printRows(50);  
        theSHM.writeCountFile(outputSeqFile);
    }


    long[] getEncoded64Sequence(String s, boolean reverseComp) {
        if(reverseComp) s=BaseEncoder.getReverseComplement(s);
        if(s.length()<2*chunkSize) {
            s=s+polyA.substring(s.length());
        }
        long[] parts=new long[2];
        parts[0]=BaseEncoder.getLongFromSeq(s.substring(0, chunkSize));
        parts[1]=BaseEncoder.getLongFromSeq(s.substring(chunkSize, 2*chunkSize));
        if((parts[0]==-1)||(parts[1]==-1)) return null;
        return parts;
    }

    public void readGenome(File inFile, ReadsWPhysicalMap theSHM) {
      // only read in a fragment that is flanked by the restriction sites
      int nTags=0;
      int cutSiteTotal=0;
      int currBase=0, csbStart=0;
      String currBAC="";
      byte chromosome=-1;
      try{
          FileReader fr=new FileReader(inFile);
          BufferedReader br=new BufferedReader(fr);
          String temp="";
          StringBuilder csb, nsb; // current string builder, next string builder
          csb=new StringBuilder();
          nsb=new StringBuilder();
          Ends currEnd=Ends.Right;
          while((br.ready())&&(theSHM.getReadTotal()<(maxTags-2))) {
              temp=br.readLine().trim();
              if(theSHM.getReadTotal()%10000<1) {
                  System.out.println(currBase+":" + theSHM.getReadTotal());
              }
              if(temp.startsWith(">")) {
                  //Flush the system
                  csb=new StringBuilder();
                  nsb=new StringBuilder();
                  csbStart=currBase=0;
                  //start a new chromosome
                  currBAC="CHR:"+temp;
                  System.out.println(currBAC);
                  System.out.println("Tags at start f:"+theSHM.getReadTotal());
                  chromosome=Byte.parseByte(temp.split("chr")[1]);
                  currEnd=Ends.Right;
               } else {
                  csb.append(temp);
                  int nextSite=this.nextCutSite(csb);
                  if(nextSite>-1) {
                      ++cutSiteTotal;
                      nsb.append(csb.substring(nextSite+reOffset1, csb.length()));
                      csb.setLength(nextSite+reOffset2);
                      if (currEnd == Ends.Both) {  // telomeres will not be included in the library
                        addHaplotypes(csb.toString(), theSHM, currEnd, csbStart,chromosome);
                        nTags+=2;  // one tag for each strand
                      }
                      csb=nsb;
                      csbStart=currBase+temp.length()-csb.length();
                      nsb=new StringBuilder();
                      currEnd=Ends.Both;
                  }
                  currBase+=temp.length();
               }
          }
          System.out.println("Total tags for "+Arrays.deepToString(cuts) +" tags = " + nTags+"  RE cut sites:" + cutSiteTotal);
          fr.close();
      }
      catch(Exception e) {
          System.out.println("Catch c="+nTags+" e="+e);
          e.printStackTrace();
          System.exit(1);

      }

    }

     private void addHaplotypes(String seq, ReadsWPhysicalMap theSHM, Ends theEnds, int start, byte chromosome) {
         int distance=0;
         int length=(seq.length()<(2*chunkSize))?seq.length():2*chunkSize;
         if(theEnds==Ends.Both) distance=seq.length();
         if((theEnds==Ends.Both)||(theEnds==Ends.Left)) {
             long[] encodhaps=getEncoded64Sequence(seq.substring(0,length), false);
              if(encodhaps!=null) {  // Nb: virtual reads containing N will encode as null (18,059 reads in AGPv1; 17,943 in AGPv2)
                  theSHM.addHaplotype(encodhaps, chromosome, (byte)'+', start+1, start+(2*chunkSize),
                          (short)distance,(byte)0);
              }
         }
         if((theEnds==Ends.Both)||(theEnds==Ends.Right)) {
             long[] encodhaps=getEncoded64Sequence(seq.substring(seq.length()-length,seq.length()), true);
              if(encodhaps!=null) {  // Nb: virtual reads containing N will encode as null (18,059 reads in AGPv1; 17,943 in AGPv2)
                  theSHM.addHaplotype(encodhaps, chromosome, (byte)'-', start+seq.length()-(2*chunkSize)+1, start+seq.length(),
                          (short)distance, (byte)0);
                  if(encodhaps[0]==0) {System.out.println(seq);}
              }
         }
     }




    public void writeDistribution(File outFile, ReadCounts shf) {
    //  int f=0;
      try{
         FileWriter fw=new FileWriter(outFile,true);
          int[] freq=new int[1000];
          int[] distFreq=new int[5000];
          int c=1;
          for(int i=0; i<shf.haplotype[0].length-1; i++) {
              if(shf.hapcount[i]<distFreq.length) {
                  distFreq[shf.hapcount[i]]++;
              } else {
                  distFreq[distFreq.length-1]++;
              }
              if((shf.haplotype[0][i]==shf.haplotype[0][i+1])&&(shf.haplotype[1][i]==shf.haplotype[1][i+1])) {
               c++;
              }
              else {
              if(c>freq.length-1) {c=freq.length-1;}
              freq[c]=freq[c]+1;
              c=1;
              }
            }
          fw.write(Arrays.deepToString(cuts)+"\n");
          for(int i=0; i<freq.length; i++) {
              fw.write(i+"\t"+freq[i]+"\n");
          }
          fw.write("Distance\n");
          for(int i=0; i<distFreq.length; i++) {
              fw.write(i+"\t"+distFreq[i]+"\n");
          }
          System.out.println(Arrays.deepToString(cuts)+" ProportionSingle= "+(double)freq[1]/(double)shf.haplotype[0].length+
                  " totaltags="+shf.haplotype[0].length);
          fw.close();
      }
      catch(Exception e) {
          System.out.println("Catch c="+e);
          e.printStackTrace();
      }

    }


   int nextCutSite(int currSite, String s) {
       int minCut=999999999, t;
       for(int i=0; i<cuts.length; i++) {
           t=s.indexOf(cuts[i],currSite);
           if((t>-1)&&(t<minCut)) minCut=t;
        }
       if(minCut<999999999) return minCut;
       else return -1;
   }

      int nextCutSite(StringBuilder s) {
       int minCut=999999999, t;
       for(String cs: cuts) {
           t=s.indexOf(cs);
           if((t>-1)&&(t<minCut)) minCut=t;
        }
       if(minCut<999999999) return minCut;
       else return -1;
   }


  public static void main(String[] args) {
 //      VirtualDigester be=new VirtualDigester(null, null);
//       VirtualDigester be=new VirtualDigester(new File("E:/SolexaAnal/maize_pseudo/maize_pseudo.seq"),
//               new File("E:/SolexaAnal/test/maize_pseudo.cut.bin2"));
//       VirtualDigester be=new VirtualDigester(new File("/Users/edbuckler/SolexaAnal/maize_pseudo/maize_pseudo.seq"),
//               new File("/Users/edbuckler/SolexaAnal/maize_pseudo/maize_pseudo.cut.bin"));
       VirtualDigester be=new VirtualDigester(new File(args[0]),
               new File(args[1]));
  }

//Haplotypes positions written to:E:\SolexaAnal\test\maize_pseudo.cut.bin
//Number of Haplotypes in file:8423380

  //    String[] cuts={"CCGG","AGCT", "CTAG", "CGCG","CATG","GTAC", "GATC",
//                   "GGCC","ACGT", "TGCA", "TTAA", "TCGA", "AATT",
//                   "CCGC", "GCGC"};
//    String[] cuts={"CCGG","ACGT","CCGC"};
//    String[] cuts={"CCTC","CCCG"};
//    String[] cuts={"GGATCC","GAATTC", "AAGCTT","CTGCAG","GATATC","CATATG","TCTAGA","CTCGAG","GGTACC","GAGCTC","CCGCGG",
//"ACTAGT", "GCATGC","CCCGGG"};
  //  String[] cuts={"TCTAGA", "GCATGC","ACTAGT"};   //good 6bp cutters for enrichment
/*   String[] cuts= {"AACGTT",	"AAGCTT",	"AAGGAG",	"AATATT",	"ACATGT",	"ACCGGT",	"ACCTGC",	"ACGCGT",	"ACTAGT",
            "ACTGGG",	"AGATCT",	"AGCGCT",	"AGGCCT",	"AGTACT",	"ATCGAT",	"ATGCAT",	"ATTAAT",
            "CAACAC",	"CAATTG",	"CACGAG",	"CACGTC",	"CACGTG",	"CAGCTC",	"CAGCTG",	"CATATG",
            "CCATGG",	"CCCACA",	"CCCAGC",	"CCCGGG",	"CCGCGG",	"CCGCTC",	"CCTAGG",	"CGAACG",
            "CGATCG",	"CGGCCG",	"CGTACG",	"CGTCTC",	"CTCGAG",	"CTCTTC",	"CTGAAG",	"CTGCAG",
            "CTGGAC",	"CTGGAG",	"CTTAAG",	"CTTGAG",	"GAACCA",	"GAAGAC",	"GAATGC",	"GAATTC",
            "GACCGA",	"GACGTC",	"GAGCTC",	"GAGGAG",	"GATATC",	"GCAATG",	"GCAGTG",	"GCATGC",
            "GCCGAG",	"GCCGGC",	"GCGATG",	"GCGCGC",	"GCTAGC",	"GGATCC",	"GGCGCC",	"GGCGGA",
            "GGGCCC",	"GGTACC",	"GGTCTC",	"GTATAC",	"GTATCC",	"GTCGAC",	"GTGCAC",	"GTGCAG",
            "GTTAAC",	"TACGTA",	"TCATGA",	"TCCGGA",	"TCGCGA",	"TCGTAG",	"TCTAGA",	"TGATCA",
            "TGCGCA",	"TGGCCA",	"TGTACA",	"TTATAA",	"TTCGAA",	"TTTAAA"};
*/
//    String[] cuts= {"AACGTT", "TACGTA", "GTTAAC", "CGTACG", "CGATCG", "CACGAG"};    //best 6bp cutters

 /*   String[] cuts={"GCCGC",	"GCAGC",	"CCATC",	"ACGGC",	"GGATC",	"CCCGT",	"CTCAG",
                   "GTCTC",	"CCAGA",	"ACTGG",	"CATCG",	"CCCGC",	"GGGAC",	"GGATG",
                   "GACGC",	"CCTTC",	"GGTGA",	"GAAGA",	"GAGTC",	"GCATC",	"GGGTC",
                   "ATGAA",	"ACGGA",	"GCGAC",	"TCGTA"};
*/
//    String[] cuts={"TACGA","TCGTA"};    //top 5bp cutter and revcomp
//    String[] cuts={"GCAGC","GCTGC"};      //most cutting 5bp cutter, requires 15bp cutback
//    String[] cuts={"ATGAA","TTCAT"};      //most cutting 5bp cutter, requires 15bp cutback
 //   String[] cuts={"CCAGG","CCTGG"};     //EcoRII
//   String[] cuts={"GAATC","GATTC"};     //TfiI
//    String[] cuts={"TTGAA","TTCAA"};     //AgsI
//     String[] cuts={"GTGAC","GTCAC"};     //Tsp45I
//    String[] cuts={"GGACC","GGTCC"};     //AvaII

}
