package net.maizegenetics.genome.solexa;

import net.maizegenetics.pal.statistics.ChiSquareDistribution;
import net.maizegenetics.pal.statistics.ContigencyTable;
import net.maizegenetics.pal.statistics.FisherExact;

import java.text.DecimalFormat;
import java.util.Arrays;
import net.maizegenetics.pal.alignment.AllelePositionBLOBUtils;
import net.maizegenetics.pal.datatype.DataType;
/**
 * 
 * User: ed
 * Date: Feb 17, 2009
 * Time: 1:43:38 PM
 * 
 */
public class SNPDistV2 {
    ContigencyTable ct;
    FisherExact fe;
    String snpText, varType;
    String[] alleleDef=new String[3];
    double startPos, endPos,snpLogP, snpX2P, maxThreshold, minorAlleleProp[],pi,pivar,fst,fstLogP,
           tropHets,tempHets,tropN,tempN,popNum;
    double maxMLScore=-1;  //maximum machine learn score
    int chromosome, taxaNumWithReads, numberOfTaxa,totalFreq,k, permutations, df, afreq[],
        totalReads4Taxa[], taxaFreq[], dist[][];
    double[] avgQual; //store the average quality scores for different bases
    double[] mlQual; //store the machine learning quality score

    int[] tropical={2,3,4,5,6,7,9,10,13,14,17,19,21,22,26,27}, //CML52 left out, CML52R is the real thing
          temperate={0,1,11,12,15,16,18,20,23,24,25},
          all;
    static String[] taxaOrder;
    int altAlleleNumber;  //the number of the alternate allele - either 1 or 2.

    //attributes good for inbreds


    /**=Parse the Allele read counts from split line of file representing one snp into int [][]*/
    public SNPDistV2(int maxTaxa, String snpText, ContigencyTable ct, FisherExact fe, int altAllele) {
       this.ct=ct;
       this.fe=fe;
       this.numberOfTaxa=maxTaxa;
       this.snpText=snpText;
       this.altAlleleNumber=altAllele;
       chromosome=-1;
       startPos=endPos=-1;
       afreq=new int[3];
       dist=new int[3][maxTaxa];
       avgQual=new double[3];
       minorAlleleProp=new  double[maxTaxa];
       totalReads4Taxa=new int[maxTaxa];
       all=new int[maxTaxa];
       permutations=1000;
       snpLogP=snpX2P=maxThreshold=pi=pivar=fst=fstLogP=-1;
       popNum=2;
       k=1;//segregating sites to be analyzed in calcpi
       df=1;//df for Chi2 test of Fst values
       taxaNumWithReads=totalFreq=0;
       tropHets=tempHets=tropN=tempN=0;
 //      String s[]=snpText.split("\t");
       String s[]=snpText.split("\\s");
       if(s.length==(maxTaxa+8)) parseAlleleFreqNew(s);
       else parseAlleleProp(s);
    }

   public SNPDistV2(int maxTaxa, String snpText, ContigencyTable ct, FisherExact fe) {
       this(maxTaxa, snpText, ct, fe, 1);
   }

   private void parseAlleleFreqNew(String s[]) {
        try{
           //begin reading row from file and assigning appropriate variables
           int colPos=0;
           String chrString=s[colPos++];
           if(chrString.contains("r")) {chrString=chrString.split("r")[1];}
           chromosome=Integer.parseInt(chrString);
           startPos=Double.parseDouble(s[colPos++]);
           endPos=Double.parseDouble(s[colPos++]);
           varType=s[colPos++];
           String[] a=s[colPos++].split("/");
           for(int i=0; i<a.length; i++) alleleDef[i]=a[i];
           a=s[colPos++].split("/");
           if(!a[0].equals("NA")) {
               avgQual[0]=Double.parseDouble(a[0]);
               avgQual[1]=Double.parseDouble(a[1]);
               if(a.length==3) avgQual[2]=Double.parseDouble(a[2]);
                }
             else {
                avgQual[0]=avgQual[1]=avgQual[2]=99;
             }

           //total counts
           a=s[colPos++].split("/");
           afreq[0]=Integer.parseInt(a[0]);
           afreq[1]=Integer.parseInt(a[1]);
           afreq[2]=Integer.parseInt(a[2]);
           totalFreq=afreq[0]+afreq[1]+afreq[2];

           colPos++; //ignore total taxa counts
           for (int i = 0; i < numberOfTaxa; i++) {
                a=s[colPos++].split("/");
                for(int j=0; j<3; j++) {
                    totalReads4Taxa[i]+=dist[j][i] = Integer.parseInt(a[j]);
                }
           }
           //end reading from row in file

           //calculate proportion of reads coming from alt 1
           for (int i=0; i<numberOfTaxa; i++) {
                 all[i]=i;
                 minorAlleleProp[i]=((dist[0][i]+dist[altAlleleNumber][i])>0)?
                   ((double)dist[altAlleleNumber][i]/((double)dist[0][i]+(double)dist[altAlleleNumber][i])):Double.NaN;
                 if(Double.isNaN(minorAlleleProp[i])==false) taxaNumWithReads++;
              //   System.out.println(startPos+" "+i+" "+dist[0][i]+" "+dist[1][i]);
           }
        } catch (Exception e)
            {
            System.out.println("parseAlleleFreq ERROR: " + e +"\t" +snpText);
            e.printStackTrace();
        }
   }

   private void parseAlleleProp(String s[]) {
        try{

           //begin reading row from file and assigning appropriate variables
           int colPos=0;
           chromosome=Integer.parseInt(s[colPos++]);
           startPos=Double.parseDouble(s[colPos++]);
           endPos=Double.parseDouble(s[colPos++]);
           varType=s[colPos++];
           String[] a=s[colPos++].split("/");
           for(int i=0; i<a.length; i++) alleleDef[i]=a[i];
           avgQual[0]=Double.parseDouble(s[colPos++]);
           avgQual[1]=Double.parseDouble(s[colPos++]);
           totalFreq=Integer.parseInt(s[colPos++]);
           afreq[0]=Integer.parseInt(s[colPos++]);
           afreq[1]=Integer.parseInt(s[colPos++]);
           for (int i = 0; i < numberOfTaxa; i++) { minorAlleleProp[i]= Double.parseDouble(s[colPos++]);}
           taxaNumWithReads=Integer.parseInt(s[colPos++]);
           tropN=Double.parseDouble(s[colPos++]);
           tempN=Double.parseDouble(s[colPos++]);
           snpLogP=Double.parseDouble(s[colPos++]);
           snpX2P=Double.parseDouble(s[colPos++]);
           maxThreshold=Double.parseDouble(s[colPos++]);
           pi=Double.parseDouble(s[colPos++]);
           pivar=Double.parseDouble(s[colPos++]);
           fst=Double.parseDouble(s[colPos++]);
           fstLogP=Double.parseDouble(s[colPos++]);
           maxMLScore=Double.parseDouble(s[colPos++]);
        }
        catch (Exception e)
            {System.err.println("parseAlleleProp ERROR: " + e +"\t" +snpText);
            e.printStackTrace();
        }
   }


   public String toStringNuc() {
        StringBuilder sb=new StringBuilder();
        sb.append(chromosome+"\t"+startPos+"\t"+endPos+"\t"+varType+"\t"+alleleDef+"\t");
        for (int i=0; i<numberOfTaxa; i++)
            sb.append(minorAlleleProp[i]+"\t"+dist[0][i]+"\t"+dist[1][i]+"\t");
        sb.append(taxaNumWithReads+"\t"+tropN+"\t"+tempN+"\t"+snpLogP+"\t"+maxMLScore+"\t"+maxThreshold+
                "\t"+pi+"\t"+pivar+"\t"+fst+"\t"+fstLogP);
       return sb.toString();
    }

   public static String toStringMachineHeader(int line) {
      StringBuilder sb=new StringBuilder();
      sb.append("CHROMOSOME\tSTARTPOS\tENDPOS\tPOLYMORPHISM_TYPE\tALLELE_DEFINITION\tLINE\t");
        sb.append("QualRef\tQualAlt1\tQualAlt2\t");
        sb.append("FreqAll\tFreqAll_Ref\tFreqAll_Alt1\tFreqAll_Alt2\t");
        sb.append("AllMinorAlleleProp\t");
        sb.append("FreqLine\tFreqLine_Ref\tFreqLine_Alt1\tFreqLine_Alt2\t");
        sb.append("LineMinorAlleleProp\t");
        sb.append("taxaNum\tsnpLogP\tsnpLogPX2\tmaxThreshold\t");
      //  sb.append("maxdiscrimScore\tlinediscrimScore\t");
      return sb.toString();

   }

      public String toStringMachine(int line, boolean returnRefLike) {
        if((dist[1][line]==0)&&(returnRefLike==false)) return null;
        StringBuilder sb=new StringBuilder();
        sb.append(chromosome+"\t"+startPos+"\t"+endPos+"\t"+varType+"\t"+alleleDef+"\t"+taxaOrder[line]+"\t");
        sb.append(avgQual[0]+"\t"+avgQual[1]+"\t"+avgQual[2]+"\t");
        sb.append(totalFreq+"\t"+afreq[0]+"\t"+afreq[1]+"\t"+afreq[2]+"\t");
        sb.append((double)afreq[1]/(double)totalFreq+"\t");
        sb.append(totalReads4Taxa[line]+"\t"+dist[0][line]+"\t"+dist[1][line]+"\t"+dist[2][line]+"\t");
        sb.append(minorAlleleProp[line]+"\t");
        sb.append(taxaNumWithReads+"\t"+snpLogP+"\t"+snpX2P+"\t"+maxThreshold+"\t");
       // sb.append(maxMLScore+"\t"+getMachineLearnScore(line)+"\t");
       return sb.toString();
    }


      //BAC position B73orCML52 TotalRead RefReadTotal Alt1ReadTotal Alt2ReadTotal alt1Prop Alt2prop CML52RefCnt CML52AltCnt CML52Ratio FisherLODTest

public static String toStringPropHeader(String[] taxaOrder) {
        StringBuilder sb=new StringBuilder();
        sb.append("CHROMOSOME\tSTARTPOS\tENDPOS\tPOLYMORPHISM_TYPE\tALLELE_DEFINITION\t");
         sb.append("RefQual\tAlt1Qual\tTotalCnt\tRefCnt\tAlt1Cnt\t");
        for(int i=0; i<taxaOrder.length;i++){sb.append(taxaOrder[i]+"\t");}
        sb.append("TAXA_W/_READS\tTROPICAL_W/_READS\tTEMPERATE_W/_READS\tSNP_LOG_P\tSNP_LOG_PX2\t" +
                     "MAX_THRESHOLD\tPI\tPI_VAR\tFST\tFST_LOG_P\tmaxMLScore\t");
        return sb.toString();
          }

   public String toStringProp() {
        DecimalFormat f = new DecimalFormat("#0.##");
        DecimalFormat f2 = new DecimalFormat("#0.0");
        StringBuilder sb=new StringBuilder();
        sb.append(chromosome+"\t"+f2.format(startPos)+"\t"+f2.format(endPos)+"\t"+varType+"\t"+alleleDef[0]+"/"+alleleDef[altAlleleNumber]+"\t");
        sb.append(f2.format(avgQual[0])+"\t"+f2.format(avgQual[altAlleleNumber])+"\t"+totalFreq+"\t"+afreq[0]+"\t"+afreq[altAlleleNumber]+"\t");
        for (int i=0; i<numberOfTaxa; i++)
            {if(Double.isNaN(minorAlleleProp[i])) sb.append("NaN\t");
            else sb.append(f.format(minorAlleleProp[i])+"\t");}
        sb.append(taxaNumWithReads+"\t"+tropN+"\t"+tempN+"\t"+snpLogP+"\t"+snpX2P+"\t"+maxThreshold+
                "\t"+pi+"\t"+pivar+"\t"+fst+"\t"+fstLogP+"\t"+maxMLScore+"\t");
       return sb.toString();
    }

   public static String toStringAttributesHeader() {
        StringBuilder sb=new StringBuilder();
        sb.append("CHROMOSOME\tSTARTPOS\tENDPOS\tPOLYMORPHISM_TYPE\tALLELE_DEFINITION\t");
         sb.append("RefQual\tAlt1Qual\tTotalCnt\tRefCnt\tAlt1Cnt\t");
        sb.append("TAXA_W/_READS\tSNP_LOG_P\tSNP_LOG_PX2\t");
        return sb.toString();
          }

   public String toStringAttributes() {
        DecimalFormat f = new DecimalFormat("#0.##");
        DecimalFormat f2 = new DecimalFormat("#0.0");
        StringBuilder sb=new StringBuilder();
        sb.append(chromosome+"\t"+f2.format(startPos)+"\t"+f2.format(endPos)+"\t"+varType+"\t"+alleleDef[0]+"/"+alleleDef[altAlleleNumber]+"\t");
        sb.append(f2.format(avgQual[0])+"\t"+f2.format(avgQual[altAlleleNumber])+"\t"+totalFreq+"\t"+afreq[0]+"\t"+afreq[altAlleleNumber]+"\t");
        sb.append(taxaNumWithReads+"\t"+snpLogP+"\t"+snpX2P+"\t");
       return sb.toString();
    }


    public static String toStringPropQualHeader() {
        StringBuilder sb=new StringBuilder();
        sb.append("BAC_ACCESSION\tPOSITION_IN_BAC\tPOLYMORPHISM_TYPE\tALLELE_DEFINITION\t");
        for(int i=0; i<taxaOrder.length;i++){sb.append(taxaOrder[i]+"\t");}
         for(int i=0; i<taxaOrder.length;i++){sb.append(taxaOrder[i]+"_Q\t");}
             sb.append("TAXA_W/_READS\tTROPICAL_W/_READS\tTEMPERATE_W/_READS\tSNP_LOG_P\t" +
                     "MAX_THRESHOLD\tPI\tPI_VAR\tFST\tFST_LOG_P\tmaxMLScore\t");
        return sb.toString();
          }


    public String toStringPropQual() {
        DecimalFormat f = new DecimalFormat("#0.00");
        f.getDecimalFormatSymbols().setNaN("X");
        StringBuilder sb=new StringBuilder();
        sb.append(chromosome+"\t"+startPos+"\t"+endPos+"\t"+varType+"\t"+alleleDef+"\t");
        for (int i=0; i<numberOfTaxa; i++)
        {if(Double.isNaN(minorAlleleProp[i])) sb.append("NaN\t");
            else sb.append(f.format(minorAlleleProp[i])+"\t");}
        for (int i=0; i<numberOfTaxa; i++)
        {if(Double.isNaN(mlQual[i])) sb.append("NaN\t");
            else sb.append(f.format(mlQual[i])+"\t");}
        sb.append(taxaNumWithReads+"\t"+tropN+"\t"+tempN+"\t"+snpLogP+"\t"+maxThreshold+
                "\t"+pi+"\t"+pivar+"\t"+fst+"\t"+fstLogP+"\t");
       sb.append(maxMLScore+"\t");
       return sb.toString();
    }

   public void scoreSNP() {
       if(afreq[altAlleleNumber]<3) return;
       int[][] nzcounts=eliminateZeroColumns(dist);
       try{
           ct.setMatrix(nzcounts);
           snpLogP=-Math.log10(ct.calcRapidMonteCarloExactTest(permutations)+(1.0/(double)permutations));
       }catch(Exception e) {
           System.err.println(e);
           snpLogP=-1.0;
           System.out.println(Arrays.deepToString(nzcounts));
           e.printStackTrace();
       }
   }

   public void scoreSNPX2ThenContigency() {
      scoreSNPFastX2();
      if(snpX2P>2) scoreSNP();

   }

    public void scoreSNPFastX2() {
      // snpLogP=-1.0;
       if((afreq[0]<2)||(afreq[altAlleleNumber]<2)) return;
       int[][] nzcounts=eliminateZeroColumns(dist);
       double sum=0, sumr[], sumc[];
       sumr=new double[nzcounts.length];
       sumc=new double[nzcounts[0].length];
       for(int r=0; r<nzcounts.length; r++) {
           for(int c=0; c<nzcounts[0].length; c++) {
               sum+=nzcounts[r][c];
               sumr[r]+=nzcounts[r][c];
               sumc[c]+=nzcounts[r][c];
           }
       }
       double g=0, exp, dev;
       for(int r=0; r<nzcounts.length; r++) {
           for(int c=0; c<nzcounts[0].length; c++) {
               exp=sumr[r]*sumc[c]/sum;
               dev=(nzcounts[r][c]-exp);
               g+=dev*dev/exp;
           }
       }
      double p=1-ChiSquareDistribution.cdf(g, nzcounts[0].length-1);
      snpX2P=-Math.log10(p+0.000001);
   }

   public int[][] eliminateZeroColumns(int[][] c) {
        int gc=0;
        for(int i=0; i<c[0].length; i++) {
            if((c[0][i]+c[altAlleleNumber][i])>0) gc++;}
        int[][] nc=new int[2][gc];
        gc=0;
        for(int i=0; i<c[0].length; i++) {
            if((c[0][i]+c[altAlleleNumber][i])>0) {
                nc[0][gc]=c[0][i];
                nc[1][gc]=c[altAlleleNumber][i];
                gc++; }}
           return nc; }

   public void calcMaxThresholdPiFst() {
        double minP=2, currP, bestThreshold=0.5;
        minP=findFisherP(bestThreshold);
        for(double currThres=.01; currThres<1; currThres+=0.05) {
           currP=findFisherP(currThres);
           if(currP<minP) {minP=currP;
               bestThreshold=currThres;}}
        maxThreshold=bestThreshold;
        calcPiMeanVar();
        calcFstFstLogP();
    }

   private double findFisherP(double givenThreshold) {
        int[][] contig=new int[2][2];
        for (int i=0; i<numberOfTaxa; i++) {
                if(minorAlleleProp[i]>givenThreshold) {
                    contig[0][0]+=dist[0][i];
                    contig[1][0]+=dist[1][i];
                } else {
                   contig[0][1]+=dist[0][i];
                   contig[1][1]+=dist[1][i];}}
         return fe.getCumlativeP(contig[0][0], contig[1][0], contig[0][1], contig[1][1]);
   }

    private void calcPiMeanVar() {
        int taxaAboveThreshold=0, taxaBelowThreshold=0;
        for (int i=0; i<all.length; i++)
            if(Double.isNaN(minorAlleleProp[i])==false)
                if(minorAlleleProp[all[i]]>maxThreshold) taxaAboveThreshold++;
                    else taxaBelowThreshold++;
        double n=taxaAboveThreshold+taxaBelowThreshold;
        pi=(double)taxaAboveThreshold*(double)taxaBelowThreshold/(double)(n*(n-1)/2);
        double b1=(n+1)/(3*(n-1));
        double b2=(2*((n*n)+n+3))/(9*n*(n-1));
        pivar=((b1/(double)k)*pi)+(b2*pi*pi);
    }

   private void calcFstFstLogP() {
   double tropSNP = 0, tempSNP=0;
   double tropfreq=0, tempfreq=0;
     for (int i=0; i<all.length; i++){
            if(Double.isNaN(minorAlleleProp[i])==false){
                for(int j=0;j<tropical.length;j++) {
                    if(all[i]==tropical[j]) {tropN++;
                        if(minorAlleleProp[i]>maxThreshold) tropSNP++;}}
                for(int j=0;j<temperate.length;j++) {
                    if(all[i]==temperate[j]) {tempN++;
                        if(minorAlleleProp[i]>maxThreshold) tempSNP++;}}
            }}
    tropfreq=tropSNP/tropN;
    tempfreq=tempSNP/tempN;
    double n_mean = tropN/popNum + tempN/popNum;
    double n_var = (popNum*n_mean - (((tropN*tropN)/
                   (popNum*n_mean))+((tempN*tempN)/(popNum*n_mean))))/(popNum-1);
    double p_mean = ((tropN*tropfreq)/(popNum*n_mean))+
                    ((tempN*tempfreq)/(popNum*n_mean));
    double s_var = ((tropN*((tropfreq-p_mean)*(tropfreq-p_mean)))/
                   ((popNum-1)*n_mean))+((tempN*((tempfreq-p_mean)*
                   (tempfreq-p_mean)))/((popNum-1)*n_mean));
    double h_mean = ((tropN*tropHets)/(popNum*n_mean))+((tempN*tempHets)/
                    (popNum*n_mean));
    double a = n_mean/n_var*(s_var - ((1/(n_mean-1))*((p_mean*(1-p_mean))
               - ((s_var*(popNum-1))/popNum) - (h_mean/4))));
    double b = (n_mean/(n_mean-1))*((p_mean*(1-p_mean))-(((popNum-1)*s_var)/
               popNum)-((((2*n_mean)-1)*h_mean)/(4*n_mean)));
    double c = h_mean/2;
    fst = a/(a+b+c);
    double chi2=(popNum-1)*n_mean*fst;
    fstLogP=-Math.log10(ChiSquareDistribution.pdf(chi2, df));
   }

  public double scoreMaxMachineLearnScore() {
       double maxScore=0;
       mlQual=new double[numberOfTaxa];
       for (int i=0; i<numberOfTaxa; i++) {
           mlQual[i]=getMachineLearnScore(i);
           if(mlQual[i]>maxScore) maxScore=mlQual[i];
       }
       this.maxMLScore=maxScore;
       return maxScore;
   }

   private double getMachineLearnScoreOld(int line) {
       if(minorAlleleProp[line]==0) return 0.0;
       double score=-0.70017;
       score+=(2.58224*minorAlleleProp[line]);
       score+=(0.29189*snpLogP);
       score+=(-1.47105*maxThreshold);
       score+=(0.15049*((double)afreq[1]/(double)totalFreq));
       score+=(0.05771*dist[1][line]);
       score+=(-0.00166*afreq[1]);
       return score;
   }

    private double getMachineLearnScore(int line) {
           //this returns p-values
       double score=-Math.log10(0.5);
       if((minorAlleleProp[line]>0.9)&&(afreq[1]>1)&&(avgQual[1]>15)) {
           if(snpLogP>2) {score=-Math.log10(0.02);}
           else {score=-Math.log10(0.12);}
       }
       return score;
   }

    /**
     * Returns the number of lines with homozgous states for the ref allele,
     * alternate allele, and combined for both.
     * @return
     */
    public int[] getHomozygousCounts() {
        int[] cnt=new int[3];
        for (int i=0; i<numberOfTaxa; i++) {
             if(Double.isNaN(minorAlleleProp[i])==false) {
                 cnt[2]++;
                 if(minorAlleleProp[i]==0) cnt[0]++;
                 if(minorAlleleProp[i]==1) cnt[1]++;
             }
        }
        return cnt;
    }

    public double getMedianRelativeDepth(double[] chrAvgDepths) {
        if(chrAvgDepths.length!=numberOfTaxa) return Double.NaN;
        double[] relativeDepth=new double[numberOfTaxa];
        for (int t = 0; t < numberOfTaxa; t++) {
            relativeDepth[t]=(double)totalReads4Taxa[t]/chrAvgDepths[t];
        }
        Arrays.sort(relativeDepth);
        return relativeDepth[numberOfTaxa/2];
    }

    public double getMedianRelativeDepthOfPresent(double[] chrAvgDepths) {
        if(chrAvgDepths.length!=numberOfTaxa) return Double.NaN;
        double[] relativeDepth=new double[numberOfTaxa];
        int missing=0;
        for (int t = 0; t < numberOfTaxa; t++) {
            relativeDepth[t]=(double)totalReads4Taxa[t]/chrAvgDepths[t];
            if(totalReads4Taxa[t]==0) missing++;
        }
  //      System.out.print(this.startPos+" "+totalReads4Taxa[0]+" "+relativeDepth[0]+" ");
        Arrays.sort(relativeDepth);
        int index=((numberOfTaxa-missing)/2)+missing;
   //     System.out.println(Arrays.toString(dist));
     //   if(relativeDepth[index]>1.5) System.out.println(Arrays.toString(relativeDepth));
        return relativeDepth[index];
    }

    public long[][] getAllelesInBits(boolean ignoreHeterozygous) {
        int lgPerSite = (numberOfTaxa / 64) + 1;
        long[][] seq = new long[2][lgPerSite];
        for (int j = 0; j < numberOfTaxa; j++) {
            if(ignoreHeterozygous&&(minorAlleleProp[j]<1)&&(minorAlleleProp[j]>0)) continue;
            int index=j/64;
            int offset=j%64;
            if (minorAlleleProp[j]<1) {  //reference alleles
                seq[0][index]=seq[0][index]|(1L<<offset);
            } 
            if (minorAlleleProp[j]>0) {  //alt alleles
                seq[1][index]=seq[1][index]|(1L<<offset);
            }
        }
        return seq;
    }

    public long[][] getAllelesInBits(boolean ignoreHeterozygous, int[] reDirect, int redirectTaxaCnt) {
        int lgPerSite = (redirectTaxaCnt / 64) + 1;
        long[][] seq = new long[2][lgPerSite];
        for (int j = 0; j < numberOfTaxa; j++) {
            if(reDirect[j]<0) continue;
            if(ignoreHeterozygous&&(minorAlleleProp[j]<1)&&(minorAlleleProp[j]>0)) continue;
            int index=reDirect[j]/64;
            int offset=reDirect[j]%64;
            if (minorAlleleProp[j]<1) {  //reference alleles
                seq[0][index]=seq[0][index]|(1L<<offset);
            }
            if (minorAlleleProp[j]>0) {  //alt alleles
                seq[1][index]=seq[1][index]|(1L<<offset);
            }
        }
        return seq;
    }

    public byte[] getAllelesInBytes(boolean includeHets, double refToHetThreshold, double hetToAltThreshold) {
        byte altBase = 'S', refBase = 'S', hetBase='S';
        refBase = (byte)this.alleleDef[0].charAt(0);
        altBase = (byte)this.alleleDef[altAlleleNumber].charAt(0);
        byte[] stateByte={refBase, altBase};
        hetBase=AllelePositionBLOBUtils.getBaseFromHalfByte(AllelePositionBLOBUtils.getHalfByteFromSNPValue(stateByte));
        if(this.varType.equals("IDP")) {
            altBase=(byte)this.alleleDef[altAlleleNumber].charAt(0);
            refBase=(altBase==(byte)'-')?(byte)'+':(byte)'-';
            hetBase=(byte)'0';
        }
        if(includeHets==false) {hetBase=DataType.UNKNOWN_CHARACTER;}
        byte[] seq = new byte[numberOfTaxa];
        for (int j = 0; j < numberOfTaxa; j++) {
            if (Double.isNaN(minorAlleleProp[j])) {seq[j]=DataType.UNKNOWN_CHARACTER;}
            else if(minorAlleleProp[j]<refToHetThreshold) {seq[j] = refBase;}
            else if(minorAlleleProp[j]>hetToAltThreshold) {seq[j]=altBase;}
            else {seq[j]=hetBase;}
        }
        return seq;
    }

    public byte[] getAllelesInBytes(boolean includeHets, double refToHetThreshold, double hetToAltThreshold, int[] reDirect, int redirectTaxaCnt) {
        byte altBase = 'S', refBase = 'S', hetBase='S';
        refBase = (byte)this.alleleDef[0].charAt(0);
        altBase = (byte)this.alleleDef[altAlleleNumber].charAt(0);
        byte[] stateByte={refBase, altBase};
        hetBase=AllelePositionBLOBUtils.getBaseFromHalfByte(AllelePositionBLOBUtils.getHalfByteFromSNPValue(stateByte));
        if(this.varType.equals("IDP")) {
            altBase=(byte)this.alleleDef[altAlleleNumber].charAt(0);
            refBase=(altBase==(byte)'-')?(byte)'+':(byte)'-';
            hetBase=(byte)'0';
        }
        if(includeHets==false) {hetBase=DataType.UNKNOWN_CHARACTER;}
        byte[] seq = new byte[redirectTaxaCnt];
        for (int j = 0; j < numberOfTaxa; j++) {
            if(reDirect[j]<0) continue;
            if (Double.isNaN(minorAlleleProp[j])) {seq[reDirect[j]]=DataType.UNKNOWN_CHARACTER;}
            else if(minorAlleleProp[j]<refToHetThreshold) {seq[reDirect[j]] = refBase;}
            else if(minorAlleleProp[j]>hetToAltThreshold) {seq[reDirect[j]]=altBase;}
            else {seq[reDirect[j]]=hetBase;}
        }
        return seq;
    }


    /*
    Variable                          0             1

                           Constant                   -0.70017      -0.19431
                           minorAlleleProp_line_       2.58224       0.80922
                           snpLogP                     0.29189       0.22904
                           maxThreshold               -1.47105      -1.45200
                           AllminorAlleleProp          0.15049       0.54025
                           dist_1__line_               0.05771       0.05157
                           afreq_1_                   -0.00166      -0.00140
      
    */
}
