/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.gbs.pipeline;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Map;
import java.util.List;
import java.io.FileReader;
import java.io.BufferedReader;
import java.text.DecimalFormat;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.TreeMap;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.genome.BaseEncoder;
import net.maizegenetics.pal.alignment.Alignment;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava3.alignment.template.Profile;
import org.biojava3.alignment.Alignments;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import net.maizegenetics.pal.alignment.SimpleAlignment;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.pal.datatype.IUPACNucleotides;
import net.maizegenetics.pal.alignment.AnnotatedAlignmentUtils;
import net.maizegenetics.pal.datatype.DataType;

/**
 *
 * @author Qi
 */
public class vcf {
    public String[] TAXA;
    public ArrayList <vcf_record> sites;
    int alleleskept=3;
    private final static DecimalFormat df = new DecimalFormat("#.##");
    public vcf(String importfile, int ak) 
    {
        alleleskept = ak;
        System.out.println ("Process " + importfile);
        sites = new ArrayList <vcf_record>();
           try 
           {
                BufferedReader br = new BufferedReader(new FileReader(importfile), 65536);
                String strLine;
                String[] inputLine;
                 while ((strLine = br.readLine()) != null) 
                 {
                     if (strLine.startsWith("##")) continue;
                     if (strLine.startsWith ("#CHROM"))
                     {
                         inputLine = strLine.trim().split("\\s");
                         TAXA = Arrays.copyOfRange(inputLine, 9, inputLine.length);
                     //System.out.println(TAXA.length + " " + inputLine.length + " " + inputLine[9]);
                     //System.exit(0);
                     
                         break;
                     }
                 }
                 //System.out.println(TAXA.length);
                 //System.exit(0);
                 
                 String patternStr = "NS=([0-9]+)";
                 Pattern NSpattern = Pattern.compile(patternStr);
                 patternStr = "DP=([0-9]+)";
                 Pattern DPpattern = Pattern.compile(patternStr);
                 patternStr = "AF=([0-9.,]+)";
                 Pattern AFpattern = Pattern.compile(patternStr);
                 
                 Matcher matcher;
                 
                 while ((strLine = br.readLine()) != null) 
                 {
                     inputLine = strLine.trim().split("\\s");
                     //public  vcf_record (String chrom, int pos, String id, char ref, String alt, int qual, String filter, int ns, int dp, double[] af, ArrayList <genotype> GTs)
                     if (inputLine.length<10)
                     {
                         continue;
                     }

                     
                     String CHROM = inputLine[0];
                     int POS = Integer.parseInt(inputLine[1]);
                     //String ID = inputLine[2];	
                     char REF = inputLine[3].charAt(0);
                     String ALT	= inputLine[4];
                     int QUAL = Integer.parseInt(inputLine[5]); 
                     String FILTER = inputLine[6];
                     String INFO = inputLine[7];
                     
                     matcher = NSpattern.matcher(INFO);
                     matcher.find();
                     int NS = Integer.parseInt(matcher.group(1));
                     matcher = DPpattern.matcher(INFO);
                     matcher.find();
                     int DP = Integer.parseInt(matcher.group(1));
                     matcher = AFpattern.matcher(INFO);
                     matcher.find();
                     String[] AFstr = matcher.group(1).split(",");
                     double[] AFarray = new double[AFstr.length];
                     //double ref_frq = 1.0;
                     for (int i=0; i<AFstr.length; i++)
                     {
                         AFarray[i] = Double.parseDouble(AFstr[i]);
                         //ref_frq -= AFarray[i+1];
                     }
                     //AFarray[0] = ref_frq;
                     
                     byte[][] gts = new byte[alleleskept+1][TAXA.length];
                     String[] genoytpeStrs = Arrays.copyOfRange(inputLine, 9, inputLine.length);

                     for (int i=0; i<genoytpeStrs.length; i++)
                     {
                         String gt =genoytpeStrs[i];
                         genotype4vcf gt_record = new genotype4vcf(gt, REF, ALT);
                         byte[] vcftyptes = gt_record.make_vcf_bytes(alleleskept);
                         //System.out.println(gt_record.allele2reads.values().toString());
                         for (int k=0; k<=alleleskept; k++)
                         {
                             gts[k][i] = vcftyptes[k];
                         }
                         
                     }
                     
                     //System.out.println("pos " + POS);
                     boolean filter = false;
                     if (FILTER.startsWith("P")) {filter=true;}
                     vcf_record v = new vcf_record(CHROM, POS, REF, ALT, QUAL, filter, NS, DP, AFarray, gts);
                     sites.add(v);
                 }

            br.close();
           } 
           catch (Exception e) {
            System.out.println("Catch in reading TagCount file e=" + e);
            e.printStackTrace();


        }
    }
    
    public vcf(TagsByTaxa theTBT, TagsOnPhysicalMap theTOPM, int chr, int minTaxaWithLocus, double minMAF, double mxMAF,  int ak, int blocksize)
    {
        TAXA = theTBT.getTaxaNames();
        sites = new ArrayList <vcf_record>();
        alleleskept = ak;
        int tempfilecount = 0;
        theTOPM.sortTable(true);
        
        AlignedTags currTAL = new AlignedTags(Integer.MIN_VALUE, Byte.MIN_VALUE, Integer.MIN_VALUE);
        int[] currPos = null;
        int countLoci=0;
        
        for (int i = 0; i < theTOPM.getSize(); i++) 
        {
            int ri=theTOPM.getReadIndexForPositionIndex(i);  // process tags in order of physical position

        
            if(theTOPM.getChromosome(ri)!=chr) continue;    //Skip tags from other chromosomes
            if(Arrays.equals(currPos, theTOPM.getPositionArray(ri))) 
            {  // add tag to current TAL
                currTAL.addTag(ri, theTOPM, theTBT); 
            } 
            else 
            {  
                if(currTAL.getSize()>1) 
                {  // finish the current TAL
                    ArrayList <vcf_record> vcfrecords = currTAL.getSNPCallsQuant4vcf( minTaxaWithLocus, minMAF, alleleskept);
                    if (vcfrecords!=null)
                    {
                        sites.addAll(vcfrecords);
                        countLoci++;
                        //if (countLoci>100) break;
                        if(countLoci%100==0) {
                            System.out.println("AlignedTags: " + countLoci + " Sites:" + sites.size() + " Position:" + currPos[2]);
                        }
                        if (blocksize> 0 && countLoci>blocksize)
                        {
                            tempfilecount ++;
                            this.write_vcf("t" + tempfilecount + ".vcf");
                            sites = new ArrayList <vcf_record>();
                            countLoci = 0;
                        }
                    }
                }
                // start a new TAL with the current tag
                currPos=theTOPM.getPositionArray(ri);  // currPos[0]=chr (int)  // currPos[1]=strand (byte)  // currPos[2]=minPosition (int)
                if((currPos[1]!=TagsOnPhysicalMap.byteMissing)&&(currPos[2]!=TagsOnPhysicalMap.intMissing)) {  // we already know that currPos[0]==targetChromo
                    currTAL = new AlignedTags(currPos[0], (byte) currPos[1], currPos[2]);
                    currTAL.addTag(ri, theTOPM, theTBT);
                } else currPos=null;  // invalid position
            }
        }
        this.sort_by_position();
        this.filter_vcf (minTaxaWithLocus, minMAF, mxMAF, false);
    }
    
    public void write_vcf (String outputfile)
    {
        try {
                BufferedWriter bw = new BufferedWriter(new FileWriter(outputfile), 65536);
                bw.write("##fileformat=VCFv4.0"); bw.newLine();
                bw.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" ); bw.newLine();
                bw.write("##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">" ); bw.newLine();
                bw.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth (only filtered reads used for calling)\">" ); bw.newLine();
                bw.write("##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality\">" ); bw.newLine();
                bw.write("####FORMAT=<ID=PL,Number=3,Type=Float,Description=\"Normalized, Phred-scaled likelihoods for AA,AB,BB genotypes where A=ref and B=alt; not applicable if site is not biallelic\">" ); bw.newLine();
                bw.write("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">" ); bw.newLine();
                bw.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">"); bw.newLine();
                bw.write("##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">"); bw.newLine();
                bw.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
                //System.out.println("total site:" + sites.size());
                for (int i = 0; i < TAXA.length; i++) {
                        bw.write("\t" + TAXA[i]);
                }
                bw.newLine();
                for (int i = 0; i < sites.size(); i++) {
                    //System.out.println("III" + i);
                        vcf_record site = sites.get(i);
                        if (site.FILTER)
                        {
                            
                            String AFstr = "";
                            for( double a :site.AF)
                            {
                               AFstr += df.format(a) + ",";
                            }
                            if (AFstr.length()>0)
                            {
                                AFstr = AFstr.substring(0, AFstr.length()-1);
                            }
                            bw.write(site.CHROM + "\t" + site.POS + "\t" + "S"+site.CHROM + "_" + site.POS + "\t" +  site.REF + "\t" +  site.ALT + "\t" +  site.QUAL + "\tPASS\tNS=" + site.NS + ";DP=" + site.DP + ";AF=" + AFstr + "\tGT:AD:DP:GQ:PL") ;
                            byte[][] genotypes_bytes = site.genotypes_bytes;
                            ArrayList <Byte> allelelist = new ArrayList();
                            allelelist.add((byte) site.REF);
                            String[] ALTstr = site.ALT.split(",");
                            for (String a: ALTstr)
                            {
                                allelelist.add((byte) a.charAt(0));
                            }
                            

                            for (int k=0; k<TAXA.length; k++)
                            {
                                byte[] vcf_byte = new byte[alleleskept+1];
                                for (int j=0; j<=alleleskept; j++)
                                {
                                    vcf_byte[j]= genotypes_bytes[j][k];
                                }
                                //System.out.println("III" + i +" " + k  + " " + vcf_byte.toString());
                                genotype4vcf gv = new genotype4vcf(allelelist, vcf_byte);
                                
                                bw.write("\t" + gv.vcf_genotype_str);
                                
                            }
                            
  
                            bw.newLine();  
                        }
                              
                }
                bw.flush();
                bw.close();

        } catch (Exception e) {
                System.out.println(e.toString());
        }
    }
    
    public String chr2site(int siteindex)
    {
        return sites.get(siteindex).CHROM;
    }
    public void sort_by_position ()
    {
        System.out.println("initPhysicalSort");
        Collections.sort(sites, new PositionComparator());
        
    }
    public void filter_vcf (int minTaxaWithLocus, double minMAF, double mxMAF, boolean biallelic)
    {
        System.out.println(minTaxaWithLocus);
        //Collections.sort(sites, new PositionComparator());
        System.out.println ("minMAF: " + minMAF);
        System.out.println ("mxMAF: " + mxMAF);
        int left = 0;
        for (int i=0; i<sites.size(); i++)
        {
            if (sites.get(i).NS<minTaxaWithLocus)
            {
                sites.get(i).FILTER = false;
                continue;
            }
            double[] AF = sites.get(i).AF;
            if ((AF != null) && (AF.length>0))
            {
                if (AF[0]<minMAF)
                {
                    //System.out.println ("minMAF: filtered" + AF[0]);
                    sites.get(i).FILTER = false;
                    continue;
                }
            }
            if ((AF != null) && (AF.length>0))
            {
                if (AF[0]>mxMAF)
                {
                    sites.get(i).FILTER = false;
                    continue;
                }
            }
            if (biallelic)
            {
                if (AF.length>1)
                {
                    sites.get(i).FILTER = false;
                    continue;
                }
            }
            left ++;
        }
        System.out.println ("left " + left);
    }
        
    public void remove_duplicate()
    {
        System.out.println("remove duplicates");
        for (int i=1; i<sites.size(); i++)
        {
            
            if ((sites.get(i).POS == sites.get(i-1).POS)  && (sites.get(i).CHROM.equals(sites.get(i-1).CHROM)))
            {
                //merge the two records
                
                vcf_record vcf1 = sites.get(i);
                vcf_record vcf2 = sites.get(i-1);
                ArrayList <Byte> allelelist1 = new ArrayList();
                ArrayList <Byte> allelelist2 = new ArrayList();
                HashMap <Byte, Integer> all_AllelesHash = new HashMap<Byte, Integer>();
                all_AllelesHash.put((byte) vcf1.REF, 0);
                allelelist1.add((byte) vcf1.REF);
                all_AllelesHash.put((byte) vcf2.REF, 0);
                allelelist2.add((byte) vcf2.REF);
                String[] aa = vcf1.ALT.split(",");
                for (String a : aa)
                {
                    all_AllelesHash.put((byte)a.charAt(0), 0);
                    allelelist1.add((byte)a.charAt(0));
                }
                aa = vcf2.ALT.split(",");
                for (String a : aa)
                {
                    all_AllelesHash.put((byte)a.charAt(0), 0);
                    allelelist2.add((byte)a.charAt(0));
                }
                
                byte[][] genotypes1 = vcf1.genotypes_bytes;
                ArrayList<genotype4vcf> genotypelist1 = new ArrayList<genotype4vcf>();
                for (int k=0; k<TAXA.length; k ++)
                {
                    byte[] vbyte = new byte[alleleskept+1];
                    for (int a=0; a<=alleleskept; a++)
                    {
                        vbyte[a] = genotypes1[a][k];
                    }
                    genotype4vcf g = new genotype4vcf(allelelist1, vbyte);
                    genotypelist1.add(g);
                }
                
                byte[][] genotypes2 = vcf2.genotypes_bytes;
                ArrayList<genotype4vcf> genotypelist2 = new ArrayList<genotype4vcf>();
                for (int k=0; k<TAXA.length; k ++)
                {
                    byte[] vbyte = new byte[alleleskept+1];
                    for (int a=0; a<=alleleskept; a++)
                    {
                        vbyte[a] = genotypes2[a][k];
                    }
                    genotype4vcf g = new genotype4vcf(allelelist2, vbyte);
                    genotypelist2.add(g);
                }
                
                ArrayList<HashMap<Byte, Integer>> taxa2AlleleCounts = new ArrayList<HashMap<Byte, Integer>>();               
                for (int tx = 0; tx < TAXA.length; tx++) 
                {
                    HashMap<Byte,Integer> combined_alleleCounts=new HashMap <Byte,Integer>();
                    for (byte allele: all_AllelesHash.keySet())
                    {
                        combined_alleleCounts.put(allele, 0);
                    }
        
                    //System.out.println(genotypelist1.get(tx).DP);
                    
                    
                    
                    if (genotypelist1.get(tx).DP>0)
                    {
                        HashMap<Byte,Integer> alleleCounts1 = genotypelist1.get(tx).allele2reads;
                        for (byte allele: alleleCounts1.keySet())
                        {
                        if(!combined_alleleCounts.containsKey(allele)){combined_alleleCounts.put(allele, 0);}
                        int sum = combined_alleleCounts.get(allele) + alleleCounts1.get(allele);
                        if (sum>127) sum = 127;
                       combined_alleleCounts.put(allele, sum);
                        }
                    }

                    if (genotypelist2.get(tx).DP>0)
                    {
                        HashMap<Byte,Integer> alleleCounts2 = genotypelist2.get(tx).allele2reads;
                        for (byte allele: alleleCounts2.keySet())
                        {
                            //System.out.println(sites.get(i).POS);
                            if(!combined_alleleCounts.containsKey(allele)){combined_alleleCounts.put(allele, 0);}
                            int sum = combined_alleleCounts.get(allele) + alleleCounts2.get(allele);
                            if (sum>127) sum = 127;
                           combined_alleleCounts.put(allele, sum);
                        }
                    }
                    taxa2AlleleCounts.add(combined_alleleCounts);
                    //System.out.println(sites.get(i).POS + combined_alleleCounts.keySet().toString() + combined_alleleCounts.values().toString());
                }
                
                vcf_record myrecord = AlignedTags.getVCFRecord (taxa2AlleleCounts, Integer.parseInt(sites.get(i).CHROM), sites.get(i).POS, alleleskept);
                if (myrecord == null)
                {
                    sites.get(i).FILTER = false;
                    //
                    sites.get(i-1).FILTER = false;
                }
                else
                {
                    sites.set(i, myrecord);
                    //
                    sites.get(i-1).FILTER = false;
                }
                
            }
        }
    }
    
    public void merge_duplicate_taxa()
    {
        System.out.println("remove duplicate taxa");
        TreeMap<String,List<Integer>> sortedIds2=new TreeMap<String,List<Integer>>();
        int uniqueTaxa =0;
        ArrayList <String> newtaxalist= new ArrayList();
        for (int i=0; i<TAXA.length-1; i++)
        {
            String name = TAXA[i].split(":")[0].toUpperCase();
                List<Integer> l = sortedIds2.get(name);
                if (l == null) {
                    sortedIds2.put(name, l=new ArrayList<Integer>());
                    newtaxalist.add(name);
                    uniqueTaxa++;
                }
                l.add(i);
        }
        
        String []newTAXA = new String[uniqueTaxa];
        ArrayList <vcf_record> newsites = new ArrayList();
        newtaxalist.toArray(newTAXA);
        Arrays.sort(newTAXA);
        
        for (int i=1; i<sites.size(); i++)
        {
            vcf_record vcf1 = sites.get(i);
            byte[][] genotypes1 = vcf1.genotypes_bytes;
            
            byte[][] genotypes2 = new byte[alleleskept][uniqueTaxa];
            int samplewithdata =0;
            for (int k=0; k<uniqueTaxa; k ++)
            {
                List<Integer> listoftaxaid = sortedIds2.get(newTAXA[k]);
                int totalreads_per_site = 0;
                for (int taxaid: listoftaxaid)
                {
                    for (int alleleid=0; alleleid<alleleskept; alleleid++)
                    {
                        totalreads_per_site += genotypes2[alleleid][taxaid];
                        if (genotypes2[alleleid][k]>0)
                        {
                           genotypes2[alleleid][k] += genotypes2[alleleid][taxaid];
                        }
                        else
                        {
                            genotypes2[alleleid][k] = genotypes2[alleleid][taxaid];
                        }
                    }

                }
                if (totalreads_per_site>0)
                {
                    samplewithdata ++;
                }
            }
            vcf1.NS = samplewithdata;
            vcf1.genotypes_bytes = genotypes2;


        }
            TAXA = newTAXA;
    }

    
    
}
class vcf_record
{
    public String CHROM;
    public int POS;
    //public String ID;
    public char REF;
    public String ALT;
    public int QUAL;
    boolean FILTER;
    public int NS;
    public int DP;
    public double[] AF;
    //public String AFstr;
    public byte[][] genotypes_bytes;

    public  vcf_record (String chrom, int pos, char ref, String alt, int qual, boolean filter, int ns, int dp, double[] af, byte[][] GTs)
    {
         CHROM = chrom;
         POS = pos;
         //ID = id;
         REF = ref;
         ALT= alt;
        QUAL = qual;
        FILTER = filter;
        NS=ns;
        DP = dp;
        AF= af;
        genotypes_bytes = GTs;
//        AFstr = "";
//        for( double a :AF)
//        {
//           AFstr += df.format(a) + ",";
//        }
//        if (AFstr.length()>0)
//        {
//            AFstr = AFstr.substring(0, AFstr.length()-1);
//        }
    }

}
class genotype4vcf
{
    //GT:AD:DP:GQ:PL
    public final static genotype_score genoscore  = new genotype_score();
    HashMap<Byte, Integer> allele2reads;
    public ArrayList <Byte> alleleslist;
    public int[] PL= new int[3];; //[l(aa ab bb)]
    public byte[] GT; //biallelic genotypes
    public int[] GTbyIndex = new int[2]; //biallelic genotypes
    public int DP=0;
    public int GQ;
    public String vcf_genotype_str = "./.";
    
    public genotype4vcf (HashMap<Byte, Integer> genohash)
    {
        GT = new byte[2];
        GT[0] = DataType.UNKNOWN_BYTE;
        GT[1] = DataType.UNKNOWN_BYTE;
        allele2reads = genohash;
        
        for (byte a: allele2reads.keySet())
        {
            if (allele2reads.get(a)>127) allele2reads.put(a,127);
        }
        alleleslist = new ArrayList <Byte>();
        
        if (genohash.isEmpty())
        {
            return;
        }
                
        for (int readcount : genohash.values())
        {
            DP += readcount;
        }
        if (DP==0) return;
        alleleslist.addAll(genohash.keySet());
        Collections.sort(alleleslist, new HashValueComparator(genohash));
        
        if (alleleslist.size()==1)
        {
            PL = genoscore.getscore(allele2reads.get(alleleslist.get(0)), 0);
            GT[0] = alleleslist.get(0);
            GT[1] = alleleslist.get(0);
        }
        else
        {
            int[] scores = genoscore.getscore(allele2reads.get(alleleslist.get(0)), allele2reads.get(alleleslist.get(1)));
            PL[0] = scores[0];  PL[1] = scores[1]; PL[2] = scores[2];
            GQ = scores[3];
            if ((PL[1] <= PL[0]) && (PL[1] <= PL[2]))
            {
                GT[0] = alleleslist.get(0);
                GT[1] = alleleslist.get(1);
            }
            else if ((PL[0] <= PL[1]) && (PL[0] <= PL[2]))
            {
                GT[0] = alleleslist.get(0);
                GT[1] = alleleslist.get(0);
            }
            else
            {
                GT[0] = alleleslist.get(1);
                GT[1] = alleleslist.get(1);
            }
        }
//     
//        if (GT[0] == 78) 
//        {
//             System.out.println(alleleslist.get(0));
//             System.exit(0);
//             
//        }
        
    }
    public genotype4vcf (ArrayList <Byte> inputallelelist, byte[] vcf_bytes)
    {
         alleleslist = inputallelelist;
         GT = new byte[2];
         GT[0] = DataType.UNKNOWN_BYTE;
         GT[1] = DataType.UNKNOWN_BYTE;
         DP=0;
         vcf_genotype_str ="./.";
         if (vcf_bytes[0] == 0)
         {
            return;     
         }
         String gtstr = genotype_score.getgenotypefromcode(vcf_bytes[0]);

         GT[0] = (byte) gtstr.charAt(0);
         GT[1] = (byte) gtstr.charAt(1);
         GTbyIndex[0] = alleleslist.indexOf(GT[0]);
         GTbyIndex[1] = alleleslist.indexOf(GT[1]);
         if ((GTbyIndex[0]<0) || (GTbyIndex[1]<0)) 
         {
             return;
         }
         
         DP=0;
        allele2reads = new HashMap<Byte, Integer> ();
        //System.out.println("III" +  gtstr + allelelist.toString() + " " + vcf_bytes[0] + " " + vcf_bytes[1] + " " + vcf_bytes[2]);

        for (int i=0; i<alleleslist.size(); i++)
        {
           allele2reads.put(alleleslist.get(i), 0);
        }
        //System.out.println("III" +  gtstr + allelelist.toString() + " " + vcf_bytes[0] + " " + vcf_bytes[1] + " " + vcf_bytes[2]);

        for (int i=1; i<vcf_bytes.length  && i<=alleleslist.size(); i++)
        {
            allele2reads.put(alleleslist.get(i-1), (int)vcf_bytes[i]);
        }
        //allele2reads.put(GT[0], (int)vcf_bytes[1]);
        //allele2reads.put(GT[1], (int)vcf_bytes[2]);
        //System.out.println("III" +  gtstr + allelelist.toString() + " " + vcf_bytes[0] + " " + vcf_bytes[1] + " " + vcf_bytes[2]);

        String ADstr = "";
        for (int i=0; i<alleleslist.size(); i++)
        {
            int c = allele2reads.get(alleleslist.get(i));
           ADstr += c + ",";
           DP += c;
        }
         //System.out.println("III" +  gtstr + allelelist.toString() + " " + vcf_bytes[0] + " " + vcf_bytes[1] + " " + vcf_bytes[2]);

        if (ADstr.length()>0)
        {
            ADstr = ADstr.substring(0, ADstr.length()-1);
        }
       int gt1count;
       int gt2count;
       if (GT[0] == GT[1])
       {
           if (GTbyIndex[0]>0)
           {
                gt1count = allele2reads.get(alleleslist.get(0));
                gt2count = allele2reads.get(GT[0]);           
           }
           else
           {
                gt1count = allele2reads.get(GT[0]);
                gt2count =0;
                if (alleleslist.size()>1)
                {
                    gt2count = allele2reads.get(alleleslist.get(1));
                }
           }
       }
       else
       {
           gt1count = allele2reads.get(GT[0]);
           gt2count = allele2reads.get(GT[1]);
       }
       int[] scores = genoscore.getscore(gt1count, gt2count);
       PL[0] = scores[0];  PL[1] = scores[1]; PL[2] = scores[2]; GQ = scores[3];
        vcf_genotype_str = Integer.toString(GTbyIndex[0]) + "/" + Integer.toString(GTbyIndex[1]) + ":" + ADstr + ":" + DP + ":" + GQ+ ":" + PL[0] + "," + PL[1] + "," + PL[2];
        
    }
    public genotype4vcf (String vcf_str, char ref, String alt)
    {
        vcf_genotype_str = vcf_str;
        
        
        alleleslist = new ArrayList <Byte>();
        alleleslist.add((byte) ref);
        String[] alts = alt.split(",");
        for (String a : alts)
        {
            alleleslist.add((byte) a.charAt(0));
        }
        
        GT = new byte[2];
        GT[0] = DataType.UNKNOWN_BYTE;
        GT[1] = DataType.UNKNOWN_BYTE;
        DP=0;
        allele2reads = new HashMap<Byte, Integer> ();
        for (byte allele: alleleslist)
        {
           allele2reads.put(allele, 0);
        }
        
        if (vcf_str.equals("./."))
        {
           return;
        }
        
        String[] genotypeinfo = vcf_genotype_str.split(":");
        String[] genotype_str = genotypeinfo[0].split("/");
        //System.out.println ("pos" + vcf_str + " " +  alleleslist.size());
        GT[0] = alleleslist.get(Integer.parseInt(genotype_str[0]));
        GT[1] = alleleslist.get(Integer.parseInt(genotype_str[1]));
        GTbyIndex[0] = Integer.parseInt(genotype_str[0]);
        GTbyIndex[1] = Integer.parseInt(genotype_str[1]);
        DP = Integer.parseInt(genotypeinfo[2]);
        GQ = Integer.parseInt(genotypeinfo[3]);
        
        String[] PLstrs = genotypeinfo[4].split(",");
        PL[0] = Integer.parseInt(PLstrs[0]);
        PL[1] = Integer.parseInt(PLstrs[1]);
        PL[2] = Integer.parseInt(PLstrs[2]);
        
        //allele2reads = new HashMap<Byte, Integer>();
        String[] counts = genotypeinfo[1].split(",");
        for (int i=0; i<counts.length; i++)
        {
            int count = Integer.parseInt(counts[i]);
            if (count>127) count = 127;
            allele2reads.put(alleleslist.get(i), count );
        }
        //System.out.println(allele2reads.values().toString());
        
    }
    public void assignVCFIndex (ArrayList <Byte> inputallelelist)
    {
        if (allele2reads.isEmpty())
        {
            return;
        }
        alleleslist = inputallelelist;
        
        GTbyIndex[0] = alleleslist.indexOf(GT[0]);
        GTbyIndex[1] = alleleslist.indexOf(GT[1]);
//        
//        if (GTbyIndex[0]<0) 
//        {
//            System.out.println ((char) allelelist.get(0).byteValue());
//            System.out.println ((char) allelelist.get(1).byteValue());
//             System.out.println((char)GT[0]);
//             System.exit(0);
//             
//        }
        if (GTbyIndex[0] > GTbyIndex[1])
        {
            int t = PL[0];
            PL[0] = PL[2];
            PL[2] = t;
            
            t = GTbyIndex[0];
            GTbyIndex[0] = GTbyIndex[1];
            GTbyIndex[1] = t;
            
            byte tb = GT[0];
            GT[0] = GT[1];
            GT[1] = tb;
        }
        else if ((GTbyIndex[0] == GTbyIndex[1]) && (GTbyIndex[0]>0)) 
        {
            int t = PL[0];
            PL[0] = PL[2];
            PL[2] = t;
        }
        String ADstr = "";
        int total_used_reads = 0;
        for (int a=0; a<alleleslist.size(); a++)
        {
            byte allele = alleleslist.get(a);
            if (allele2reads.containsKey(allele))
            {
                int reads = allele2reads.get(allele);
                total_used_reads += reads;
                if (reads>127) reads=127;
                ADstr += Integer.toString(reads) + ",";
            }
        }
        if (total_used_reads==0)
        {
            return;
        }
                
        if (ADstr.length()>0) ADstr=ADstr.substring(0, ADstr.length()-1);
        vcf_genotype_str = Integer.toString(GTbyIndex[0]) + "/" + Integer.toString(GTbyIndex[1]) + ":" + ADstr + ":" + DP + ":" + GQ+ ":" + PL[0] + "," + PL[1] + "," + PL[2];
    }
    
    public byte[] make_vcf_bytes (int alleleskept)  //first byte: genotyping code; second to last byte: count of reads for each allele, up to the number of kept alleles
    {
        byte[] vcfbytes = new byte[alleleskept+1];
        if (GT[0] == DataType.UNKNOWN_BYTE)
        {
            vcfbytes[0]= 0;
            for (int i=1; i<=alleleskept; i++)
            {
                vcfbytes[i]= 0;
            }
        }
        else
        {
            //System.out.println("GG " + (char)GT[0] + " " + (char)GT[0]);
            vcfbytes[0] =genoscore.getcodefromgenotype((char)GT[0], (char)GT[1]);
            //byte b = (byte) (a-128);
            for (int i=1; i<=alleleskept; i++)
            {
                vcfbytes[i]= 0;
            }
            for (int i=0; i<alleleslist.size()  && i<alleleskept; i++)
            {
                vcfbytes[i+1]= (byte)(int)(allele2reads.get(alleleslist.get(i)));
            }
        }
        //System.out.println("III" + vcf_genotype_str + " " + vcfbytes[0] + " " + vcfbytes[1] + " " + vcfbytes[2]);
        //System.exit(0);
        return vcfbytes;
    }
}
class AlignedTags {
    public ArrayList<SingleTBT> theTags = new ArrayList<SingleTBT>(); 
    int startPosition;
    int chromosome;
    byte strand;
    int indexOfRef;
    int[] positionsOfVariableSites;
    //byte[][] callsAtVariableSitesByTag;
    private final static int maxSNPsPerLocus = 64;
    private final static int maxAlignmentSize = 10000;



    public AlignedTags(int chromosome, byte strand, int startPosition) {
        this.chromosome = chromosome;
        this.strand = strand;
        this.startPosition = startPosition;
        positionsOfVariableSites = null;
        //callsAtVariableSitesByTag = null;
    }
    
    public void addTag(int tagTOPMIndex, TagsOnPhysicalMap theTOPM, TagsByTaxa theTBT) {
        SingleTBT singleTBT = new SingleTBT(tagTOPMIndex, theTOPM, theTBT);
        if (singleTBT.taxaWithTag > 0) {
            theTags.add(singleTBT);
        }
    }
    
    public int getSize() {
        return theTags.size();
    }
    
    public int getChromosome() {
        return chromosome;
    }
    
    public byte getStrand() {
        return strand;
    }
    
    public int getStartPosition() {
        return startPosition;
    }
    
    public int getTOPMIndexOfTag(int tagIndex) {
        return theTags.get(tagIndex).tagTOPMIndex;
    }
    
    //public byte getCallAtVariableSiteForTag(int site, int tagIndex) {
    //    return callsAtVariableSitesByTag[site][tagIndex];
    //}
    
    public int getNumberTaxaCovered() {
        if (theTags.size() < 2) return 0;
        int nTaxaCovered = 0;
        boolean[] covered = new boolean[theTags.get(0).tagDist.length];  // initializes to false
        for (SingleTBT sTBT : theTags) {
            for (int tx=0; tx < covered.length; ++tx) {
                if (sTBT.tagDist[tx] > 0) covered[tx] = true;
            }
        }
        for (int tx=0; tx < covered.length; ++tx) {
            if (covered[tx]) ++nTaxaCovered;
        }
        return nTaxaCovered;
    }
    
    
    public ArrayList<vcf_record> getSNPCallsQuant4vcf(int minTaxaWithLocus, double minMAF, int alleleskept) {
        
        if (theTags.size() < 1) return null;
        //System.out.println("tags size: " + theTags.size());
//        if ((top2tags==true) && (theTags.size() > 2))
//        {
//           Collections.sort(theTags, new Comparator(){
//            public int compare (Object t1, Object t2)
//            {
//                SingleTBT tt1 = (SingleTBT) t1;
//                SingleTBT tt2 = (SingleTBT) t2;
//               if (tt1.readCount < tt2.readCount) {return 1;}
//               else if (tt1.readCount == tt2.readCount) {return 0;}
//                       else {return -1;}
//            }
//           });
//           int secondCount =  theTags.get(1).readCount;
//           int thirdCount =  theTags.get(2).readCount;
//           //System.out.println("tags size: " + theTags.get(1).readCount + " " + secondCount + " " + thirdCount);
//           if (((float)thirdCount/(float)secondCount)>0.5) {return null;}
//           for (int i=2; i<theTags.size(); i++) theTags.remove(i);
//        }
        Alignment tagAlignment = this.getVariableSites();
        
        if(tagAlignment==null) return null;
        int nTaxa = theTags.get(0).tagDist.length;
        int nSites = tagAlignment.getSiteCount();
        if (nSites<1 || nTaxa<1) {
            return null;
        }
        
        positionsOfVariableSites = new int[nSites];
        ArrayList <vcf_record> vcf_list = new  ArrayList<vcf_record>();
        //callsAtVariableSitesByTag = new byte[tagAlignment.getSiteCount()][theTags.size()];
        for (int s = 0; s < nSites; s++) {
            int positionsInLocus = tagAlignment.getPositionInLocus(s);
            int positionInChr = (strand == -1) ? startPosition-positionsInLocus : startPosition+positionsInLocus;
            positionsOfVariableSites[s] = positionInChr;
           ArrayList<HashMap<Byte, Integer>> taxa2AlleleCounts = new ArrayList<HashMap<Byte, Integer>>();
            for (int tx = 0; tx < nTaxa; tx++) {
                HashMap<Byte,Integer> alleleCounts = new HashMap<Byte,Integer>();
                //int readtotal =0;
                for (int tg = 0; tg < tagAlignment.getSequenceCount(); tg++) {
                    int tagIndex=Integer.parseInt(tagAlignment.getTaxaName(tg).split("_")[0]);  // taxaName in tagAlignment is set to indexInTheTags_"refSeq"|"no"
                    byte baseToAdd=tagAlignment.getBase(tg, s);
                    if (baseToAdd == 78)
                    {
                        continue;
                    }
                                        
                    if (strand==-1) baseToAdd = complement(baseToAdd);

                    //callsAtVariableSitesByTag[s][tagIndex] = baseToAdd;
                    if (alleleCounts.containsKey(baseToAdd)) {
                        alleleCounts.put(baseToAdd, alleleCounts.get(baseToAdd)+ theTags.get(tagIndex).tagDist[tx]);
                    } else {
                        alleleCounts.put(baseToAdd, (int) theTags.get(tagIndex).tagDist[tx]);
                    }
                    //readtotal += (int) theTags.get(tagIndex).tagDist[tx];
                }
                taxa2AlleleCounts.add(alleleCounts);
                //if (positionInChr==524471){System.out.println(tx + " " + readtotal);}
            }
            
            vcf_record myrecord = getVCFRecord (taxa2AlleleCounts, this.chromosome, positionsOfVariableSites[s], alleleskept);
            if ((myrecord !=null) && (myrecord.NS>minTaxaWithLocus) && (myrecord.AF[0]>minMAF))
            {
                vcf_list.add(myrecord);
            }
        }
        //vcf_record[] vcfarrays = vcf_list.toArray(new String[list.size()]);
        
        return vcf_list;
    }
    
    public int[] getPositionsOfVariableSites() {
        return positionsOfVariableSites;
    }
   private void assignRefTag() {
        indexOfRef = Integer.MIN_VALUE;
        int lengthOfRef = Integer.MIN_VALUE;
        int tagIndex = 0;
        for (SingleTBT sTBT : theTags) {
            if (sTBT.divergence==0  && sTBT.tagLength>lengthOfRef) {
                indexOfRef = tagIndex;
                lengthOfRef = sTBT.tagLength;
            }
            ++tagIndex;
        }
    }
    
     private Alignment getVariableSites() {
        this.assignRefTag();
        boolean checkReplicateProfiles = false;
        boolean printOutRefWithGaps = false;
        List<DNASequence> lst = new ArrayList<DNASequence>();
        if (theTags.size() < 2) return null;
        if(theTags.size()>maxAlignmentSize) return null;   // should we use a maxAlignmentSize (upper limit of tags) here?
        
        int tagIndex = 0;
        for (SingleTBT sTBT : theTags) {
            DNASequence ds=new DNASequence(sTBT.tagTrimmed);
            String refMark = (tagIndex==indexOfRef) ? "refSeq" : "no";
            ds.setOriginalHeader(tagIndex+"_"+refMark);    // OriginalHeader set to indexInTheTags_'refSeq'|'no'
            ds.setCompoundSet(AmbiguityDNACompoundSet.getDNACompoundSet());
            lst.add(ds);
            ++tagIndex;
        }
        Profile<DNASequence, NucleotideCompound> profile = Alignments.getMultipleSequenceAlignment(lst);
        if (checkReplicateProfiles) {
            System.out.printf("Clustal1:%d%n%s%n", startPosition,profile);
            Profile<DNASequence, NucleotideCompound> profile2 = Alignments.getMultipleSequenceAlignment(lst);
            System.out.printf("Clustal2:%d%n%s%n", startPosition, profile2);
        }
        String[] aseqs=new String[theTags.size()];
        String[] names=new String[theTags.size()];
        boolean refTagWithGaps = false;
        int[] positions = null;
        for (int i = 0; i < aseqs.length; i++) {
            aseqs[i]=profile.getAlignedSequence(i+1).getSequenceAsString();
            names[i]=profile.getAlignedSequence(i+1).getOriginalSequence().getOriginalHeader();
            if (names[i].split("_")[1].equals("refSeq")) {  // names were set to indexInTheTags_"refSeq"|"no"
                if (aseqs[i].contains("-")) {
                    refTagWithGaps = true;
                    positions = new int[aseqs[i].length()];
                    positions[0] = 0;
                    for (int site=1; site<aseqs[i].length(); site++) {
                        positions[site] = (aseqs[i].charAt(site)=='-')? (positions[site-1]) : (positions[site-1]+1);
                    }
                }
            }
        }
        profile=null;
        Alignment aa = null;
        
        if (refTagWithGaps) {
            aa = SimpleAlignment.getInstance(new SimpleIdGroup(names),aseqs, new IUPACNucleotides(), positions);
        } else {
            aa = SimpleAlignment.getInstance(new SimpleIdGroup(names),aseqs, new IUPACNucleotides());
        }
        Alignment faa = AnnotatedAlignmentUtils.removeConstantSitesIgnoreMissing(aa);
        if (printOutRefWithGaps && refTagWithGaps) { 
            System.out.println("chr"+chromosome+"  pos:"+startPosition+"  strand:"+strand+"  AA\n"+aa.toString().trim());
            System.out.println("chr"+chromosome+"  pos:"+startPosition+"  strand:"+strand+"  FAA\n"+faa.toString()); 
        }
        if(faa.getSiteCount()>maxSNPsPerLocus || faa.getSiteCount()<1 || faa.getSequenceCount()<2) return null;
        SimpleAlignment sfaa = SimpleAlignment.getInstance(faa);
        return sfaa;
    }
     
     public static vcf_record getVCFRecord (ArrayList<HashMap<Byte, Integer>> taxa2AlleleCounts, int chrom, int position, int alleleskept) 
    {        
        //DecimalFormat df = new DecimalFormat("#.##");
        
        int readdepth = 0;
        int taxaWithData = 0;
        ArrayList<genotype4vcf> genotypes= new ArrayList<genotype4vcf>();
        HashMap <Byte, Integer> allalleles2taxacount = new HashMap <Byte, Integer>();
        for (int i =0; i<taxa2AlleleCounts.size(); i++)
        {
           HashMap <Byte, Integer> allele2reads = taxa2AlleleCounts.get(i);
           // System.out.println(allele2reads.values().toString());
           genotype4vcf theGT = new genotype4vcf(allele2reads);
           //System.out.println("ddda" + allele2reads.keySet().toString() + allele2reads.values().toString() + theGT.GT[0]);
           genotypes.add(theGT);
           if (theGT.GT[0] != DataType.UNKNOWN_BYTE)
           {
               taxaWithData ++;
               readdepth += theGT.DP;
               if (allalleles2taxacount.containsKey(theGT.GT[0])){allalleles2taxacount.put(theGT.GT[0], allalleles2taxacount.get(theGT.GT[0]) + 1);}
               else {allalleles2taxacount.put(theGT.GT[0], 1);}
               if (allalleles2taxacount.containsKey(theGT.GT[1])){allalleles2taxacount.put(theGT.GT[1], allalleles2taxacount.get(theGT.GT[1]) + 1);}
               else {allalleles2taxacount.put(theGT.GT[1], 1);}
           }
        }
        //why allele count =1
        //System.out.println("DDDDDDD" + position + " " + allalleles2taxacount.size());
        if (allalleles2taxacount.size()<2){return null;}
        if (allalleles2taxacount.isEmpty()){return null;}
        
            
        ArrayList <Byte> allalleleslist = new ArrayList();
        allalleleslist.addAll(allalleles2taxacount.keySet());
        Collections.sort(allalleleslist, new HashValueComparator(allalleles2taxacount));
        
        if (allalleleslist.size()>alleleskept)
        {
            
            allalleleslist.subList(alleleskept, allalleleslist.size()).clear();
        }
        // assign teh reference allele: the major alleles
        char refallele = (char)allalleleslist.get(0).byteValue();
        String altalleles = "";
        
        //prepare the minor allele string and allele frequency
        double[] allfrq = new double[allalleleslist.size()-1];
        //System.out.println(allalleleslist.size());
        if (allalleleslist.size()>1)
        {
            //int[] MAF = new int[allalleleslist.size()];
            //Arrays.fill(MAF, 0);
            for (int i=1; i<allalleleslist.size(); i++)
            {
               altalleles = altalleles + (char)allalleleslist.get(i).byteValue() + ",";
               allfrq[i-1] = (double)allalleles2taxacount.get(allalleleslist.get(i)) / ((double)taxaWithData * 2.0);  
            }
            altalleles = altalleles.substring(0, altalleles.length() -1);
        }
        
        for (genotype4vcf theGT: genotypes)
        {
            theGT.assignVCFIndex(allalleleslist);
        }
        //public void vcf_record (String chrom, int pos, String id, String ref, String[] alt, int qual, String filter, int ns, int dp, double[] af)  
        if(allfrq.length > 0)
        {
            byte[][] gts = new byte[alleleskept+1][genotypes.size()];
            for (int i=0; i<genotypes.size(); i++)
            {
                byte[] vcfbytes = genotypes.get(i).make_vcf_bytes(alleleskept);
                for (int t =0; t<=alleleskept; t++)
                {
                    gts[t][i] = vcfbytes[t];
                }
            }
            vcf_record thisvcf_record = new vcf_record(Integer.toString(chrom), position, refallele, altalleles, 20, true, taxaWithData, readdepth, allfrq, gts);  
            return thisvcf_record;
        }
        else
        {
            return null;
        }
        
    }
   
     public static byte complement(byte geno) {
        byte comp = Byte.MIN_VALUE;
        switch (geno) {  
            case 'A':  comp = 'T';  break;
            case 'C':  comp = 'G';  break;
            case 'G':  comp = 'C';  break;
            case 'T':  comp = 'A';  break;
            case 'K':  comp = 'M';  break;
            case 'M':  comp = 'K';  break;
            case 'R':  comp = 'Y';  break;
            case 'S':  comp = 'S';  break;
            case 'W':  comp = 'W';  break;
            case 'Y':  comp = 'R';  break;
            case '-':  comp = '-';  break;  // both strands have the deletion
            case '+':  comp = '+';  break;  // both strands have the insertion
            case '0':  comp = '0';  break;
            case 'N':  comp = 'N';  break;
            default:   comp = 'N';  break;
        }
        return comp;
    }
    
    
}
class SingleTBT {
    int tagTOPMIndex;
    int readCount;
    long[] tag;
    byte tagLength;
    int divergence;
    String tagTrimmed;
    int tagTBTIndex; //index in the TBT
    int taxaWithTag;
    byte[] tagDist;  // observed count of the tag for each taxon
    private static String polyN="NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";

    SingleTBT(int tagTOPMIndex, TagsOnPhysicalMap theTOPM, TagsByTaxa theTBT) {

        this.tagTOPMIndex = tagTOPMIndex;
        tag = theTOPM.getTag(tagTOPMIndex);
        tagLength = (byte) theTOPM.getTagLength(tagTOPMIndex);
        divergence = theTOPM.getDivergence(tagTOPMIndex);
        tagTrimmed = BaseEncoder.getSequenceFromLong(tag).substring(0, tagLength);
        if(tagLength<64) tagTrimmed=tagTrimmed+polyN.substring(0, 64-tagLength);
        tagTBTIndex = theTBT.getTagIndex(tag);
        taxaWithTag = (tagTBTIndex>-1) ? theTBT.getNumberOfTaxaWithTag(tagTBTIndex) : 0;  // tags with 0 taxaWithTag will not be added to TagsAtLocus
        tagDist = (tagTBTIndex>-1) ? theTBT.getTaxaReadCountsForTag(tagTBTIndex) : null;
        readCount =0;
        if (tagDist != null)
        {
            for (int i=0; i<tagDist.length; i++)
            {
                readCount += tagDist[i];
            }
        }
    }
}

class HashValueComparator implements Comparator {

  Map base;
  public HashValueComparator(Map base) {
      this.base = base;
  }

  public int compare(Object a, Object b) {

    if((Integer)base.get(a) < (Integer)base.get(b)) {
      return 1;
    } else if((Integer)base.get(a) == (Integer)base.get(b)) {
      return 0;
    } else {
      return -1;
    }
  }
}

class PositionComparator implements Comparator<vcf_record> {
    public int compare(vcf_record v1, vcf_record v2) 
    {
        if (Integer.parseInt(v1.CHROM) < Integer.parseInt(v2.CHROM))
        {
            return -1;
        }
        else if (Integer.parseInt(v1.CHROM) == Integer.parseInt(v2.CHROM))
        {
            if (v1.POS< v2.POS)
            {
                return -1;
            }
            else if (v1.POS==v2.POS)
            {   
                return 0;
            }
            else
            {
                return 1;
            }
        }
        else
        {
            return 1; 
        }
    }
}