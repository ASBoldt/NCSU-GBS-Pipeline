package net.maizegenetics.genome.solexa;

import net.maizegenetics.pal.statistics.ChiSquareDistribution;
import net.maizegenetics.pal.statistics.ContigencyTable;
import net.maizegenetics.pal.statistics.FisherExact;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.text.DecimalFormatSymbols;
/**
 * Created by IntelliJ IDEA.
 * User: ed
 * Date: Feb 17, 2009
 * Time: 1:43:38 PM
 * To change this template use File | Settings | File Templates.
 */
public class SNPDist {
    ContigencyTable ct;
    FisherExact fe;
    String bacName, snpText, varType, alleleDef;
    double bacPos,snpLogP,maxThreshold, minorAlleleProp[],pi,pivar,fst,fstLogP,
           tropHets,tempHets,tropN,tempN,popNum;
    double maxMLScore=-1;  //maximum machine learn score
    int taxaNum, maxTaxa,totalFreq,k, permutations, df, afreq[],
        totalReads4Taxa[], taxaFreq[], dist[][];
    double[] qual; //store the machine learning quality score

    int[] tropical={2,3,4,5,6,7,8,9,12,13,16,18,20,21,25,26},
          temperate={0,1,10,11,14,15,17,19,22,23,24},
          all;
    static String[] taxaOrder={"B73","B97","CML103","CML228","CML247","CML277","CML322",
                                         "CML333","CML52","CML69","HP301","IL14H","KI11","KI3","KY21","M162W",
                                         "M37W","MO17","MO18W","MS71","NC350","NC358","OH43","OH7B","P39","TX303","TZI8"};

    /**=Parse the Allele read counts from split line of file representing one snp into int [][]*/
    public SNPDist(int maxTaxa, String snpText, ContigencyTable ct, FisherExact fe) {
       this.ct=ct;
       this.fe=fe;
       this.maxTaxa=maxTaxa;
       this.snpText=snpText;
       afreq=new int[3];
       dist=new int[3][maxTaxa];
       minorAlleleProp=new  double[maxTaxa];
       totalReads4Taxa=new int[maxTaxa];
       all=new int[maxTaxa];
       permutations=10000;
       snpLogP=maxThreshold=pi=pivar=fst=fstLogP=-1;
       popNum=2;
       k=1;//segregating sites to be analyzed in calcpi
       df=1;//df for Chi2 test of Fst values
       taxaNum=totalFreq=0;
       tropHets=tempHets=tropN=tempN=0;
       String s[]=snpText.split("\t");
       if(s.length==((maxTaxa*3)+10)) parseAlleleFreq(s);
       else if(s.length==(maxTaxa+13)) parseAlleleProp(s);
    }


   private void parseAlleleFreq(String s[]) {
        try{
           //begin reading row from file and assigning appropriate variables
           int colPos=0;
           bacName=s[colPos++];
           bacPos=Double.parseDouble(s[colPos++]);
           varType=s[colPos++];
           alleleDef=s[colPos++];
           afreq[0]=Integer.parseInt(s[colPos++]);
           afreq[1]=Integer.parseInt(s[colPos++]);
           afreq[2]=Integer.parseInt(s[colPos++]);
           totalFreq=afreq[0]+afreq[1]+afreq[2];
           colPos+=3;  //ignore total taxa counts
           for (int i = 0; i < maxTaxa * 3; i++)
                totalReads4Taxa[i/3]+=dist[i%3][i/3] = Integer.parseInt(s[colPos++]);
           //end reading from row in file
          //this corrects an error in the base file.
           if(afreq[1]<afreq[2]) {
              int temp;
              for (int i=0; i<maxTaxa; i++) {temp=dist[1][i]; dist[1][i]=dist[2][i]; dist[2][i]=temp;}
              temp=afreq[1];  afreq[1]=afreq[2];  afreq[2]=temp;
           }
           //calculate proportion of reads coming from alt 1
           for (int i=0; i<maxTaxa; i++) {
                 all[i]=i;
                 minorAlleleProp[i]=((dist[0][i]+dist[1][i])>0)?
                   ((double)dist[1][i]/((double)dist[0][i]+(double)dist[1][i])):Double.NaN;
                 if(Double.isNaN(minorAlleleProp[i])==false) taxaNum++;
           }
        }catch (Exception e)
            {System.err.println("parseAlleleFreq ERROR: " + e +"\t" +snpText);}
   }

   private void parseAlleleProp(String s[]) {
        try{
           //begin reading row from file and assigning appropriate variables
           int colPos=0;
           bacName=s[colPos++];
           bacPos=Double.parseDouble(s[colPos++]);
           varType=s[colPos++];
           alleleDef=s[colPos++];
       //these should be output in the future.
//           afreq[0]=Integer.parseInt(s[colPos++]);
//           afreq[1]=Integer.parseInt(s[colPos++]);
//           afreq[2]=Integer.parseInt(s[colPos++]);
//           totalFreq=afreq[0]+afreq[1]+afreq[2];
           for (int i = 0; i < maxTaxa; i++) { minorAlleleProp[i]= Double.parseDouble(s[colPos++]);}
           taxaNum=Integer.parseInt(s[colPos++]);
           tropN=Double.parseDouble(s[colPos++]);
           tempN=Double.parseDouble(s[colPos++]);
           snpLogP=Double.parseDouble(s[colPos++]);
           maxThreshold=Double.parseDouble(s[colPos++]);
           pi=Double.parseDouble(s[colPos++]);
           pivar=Double.parseDouble(s[colPos++]);
           fst=Double.parseDouble(s[colPos++]);
           fstLogP=Double.parseDouble(s[colPos++]);
        }
        catch (Exception e)
            {System.err.println("ERROR: " + e +"\t" +snpText);}
   }


   public String toStringNuc() {
        StringBuilder sb=new StringBuilder();
        sb.append(bacName+"\t"+bacPos+"\t"+varType+"\t"+alleleDef+"\t");
        for (int i=0; i<maxTaxa; i++)
            sb.append(minorAlleleProp[i]+"\t"+dist[0][i]+"\t"+dist[1][i]+"\t");
        sb.append(taxaNum+"\t"+tropN+"\t"+tempN+"\t"+snpLogP+"\t"+maxMLScore+"\t"+maxThreshold+
                "\t"+pi+"\t"+pivar+"\t"+fst+"\t"+fstLogP);
       return sb.toString();
    }

   public static String toStringNucHeader(int line) {
      StringBuilder sb=new StringBuilder();
      sb.append("bacName\tbacPos\tvarType\talleleDef\t");
        sb.append("totalFreq\tafreq[0]\tafreq[1]\tafreq[2]\t");
        sb.append("AllminorAlleleProp\t");
        sb.append("totalReads4Taxa[line]\tdist[0][line]\tdist[1][line]\tdist[2][line]\t");
        sb.append("minorAlleleProp[line]\t");
        sb.append("taxaNum\tsnpLogP\tmaxThreshold\t");
        sb.append("maxdiscrimScore\tlinediscrimScore\t");
      return sb.toString();

   }

      public String toStringNuc(int line) {
        StringBuilder sb=new StringBuilder();
        sb.append(bacName+"\t"+bacPos+"\t"+varType+"\t"+alleleDef+"\t");
        sb.append(totalFreq+"\t"+afreq[0]+"\t"+afreq[1]+"\t"+afreq[2]+"\t");
        sb.append((double)afreq[1]/(double)totalFreq+"\t");
        sb.append(totalReads4Taxa[line]+"\t"+dist[0][line]+"\t"+dist[1][line]+"\t"+dist[2][line]+"\t");
        sb.append(minorAlleleProp[line]+"\t");
        sb.append(taxaNum+"\t"+snpLogP+"\t"+maxThreshold+"\t");
        sb.append(maxMLScore+"\t"+getMachineLearnScore(line)+"\t");
       return sb.toString();
    }

      //BAC position B73orCML52 TotalRead RefReadTotal Alt1ReadTotal Alt2ReadTotal alt1Prop Alt2prop CML52RefCnt CML52AltCnt CML52Ratio FisherLODTest

   public String toStringProp() {
        StringBuilder sb=new StringBuilder();
        sb.append(bacName+"\t"+bacPos+"\t"+varType+"\t"+alleleDef+"\t");
        for (int i=0; i<maxTaxa; i++)
            sb.append(minorAlleleProp[i]+"\t");
        sb.append(taxaNum+"\t"+tropN+"\t"+tempN+"\t"+snpLogP+"\t"+maxThreshold+
                "\t"+pi+"\t"+pivar+"\t"+fst+"\t"+fstLogP);
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
        sb.append(bacName+"\t"+bacPos+"\t"+varType+"\t"+alleleDef+"\t");
        for (int i=0; i<maxTaxa; i++)
        {if(Double.isNaN(minorAlleleProp[i])) sb.append("NaN\t");
            else sb.append(f.format(minorAlleleProp[i])+"\t");}
        for (int i=0; i<maxTaxa; i++)
        {if(Double.isNaN(qual[i])) sb.append("NaN\t");
            else sb.append(f.format(qual[i])+"\t");}
        sb.append(taxaNum+"\t"+tropN+"\t"+tempN+"\t"+snpLogP+"\t"+maxThreshold+
                "\t"+pi+"\t"+pivar+"\t"+fst+"\t"+fstLogP+"\t");
       sb.append(maxMLScore+"\t");
       return sb.toString();
    }

   public void scoreSNP() {
       if(afreq[1]<3) return;
       int[][] nzcounts=eliminateZeroColumns(dist);
       try{
           ct.setMatrix(nzcounts);
           snpLogP=-Math.log10(ct.calcRapidMonteCarloExactTest(permutations)+.0001);
       }catch(Exception e) {System.err.println(e);
           snpLogP=-1.0;}
   }

   public void scoreSNPFast() {
      // snpLogP=-1.0;
       if(afreq[1]<3) return;
       int[][] nzcounts=eliminateZeroColumns(dist);
       double[][] nzcountd=new double[2][nzcounts[0].length];
       double sum=0, sumr[], sumc[];
       sumr=new double[nzcounts.length];
       sumc=new double[nzcounts[0].length];
       //continuity correction
       for(int c=0; c<nzcounts[0].length; c++) {
          if(nzcounts[0][c]==0) {nzcountd[0][c]=nzcounts[0][c]+0.5;  nzcountd[1][c]=nzcounts[1][c]-0.5;}
            else if(nzcounts[1][c]==0) {nzcountd[0][c]=nzcounts[0][c]-0.5;  nzcountd[1][c]=nzcounts[1][c]+0.5;}
            else {nzcountd[0][c]=nzcounts[0][c];  nzcountd[1][c]=nzcounts[1][c];}
           }
       for(int r=0; r<nzcountd.length; r++) {
           for(int c=0; c<nzcountd[0].length; c++) {
               sum+=nzcountd[r][c];
               sumr[r]+=nzcountd[r][c];
               sumc[c]+=nzcountd[r][c];
           }
       }
       double g=0, exp;
       for(int r=0; r<nzcounts.length; r++) {
           for(int c=0; c<nzcounts[0].length; c++) {
               exp=sumr[r]*sumc[c]/sum;
               g+=nzcountd[r][c]*Math.log((nzcountd[r][c]/exp));
           }
       }
      g=g*2;
      double p=1-ChiSquareDistribution.cdf(g, nzcounts[0].length-1);
      snpLogP=-Math.log10(p+0.000001);
   }

    public void scoreSNPFastX2() {
      // snpLogP=-1.0;
       if(afreq[1]<3) return;
       int[][] nzcounts=eliminateZeroColumns(dist);
       double sum=0, sumr[], sumc[];
       sumr=new double[nzcounts.length];
       sumc=new double[nzcounts[0].length];
       //continuity correction
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
      snpLogP=-Math.log10(p+0.000001);
   }

   public static int[][] eliminateZeroColumns(int[][] c) {
        int gc=0;
        for(int i=0; i<c[0].length; i++) {
            if((c[0][i]+c[1][i])>0) gc++;}
        int[][] nc=new int[2][gc];
        gc=0;
        for(int i=0; i<c[0].length; i++) {
            if((c[0][i]+c[1][i])>0) {
                nc[0][gc]=c[0][i];
                nc[1][gc]=c[1][i];
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
        for (int i=0; i<maxTaxa; i++) {
                if(minorAlleleProp[i]>givenThreshold) {
                    contig[0][0]+=dist[0][i];
                    contig[1][0]+=dist[1][i];
                } else {
                   contig[0][1]+=dist[0][i];
                   contig[1][1]+=dist[1][i];}}
         return fe.getCumlativeP(contig[0][0], contig[1][0], contig[0][1], contig[1][1]); }

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
       qual=new double[maxTaxa];
       for (int i=0; i<maxTaxa; i++) {
           qual[i]=getMachineLearnScore(i);
           if(qual[i]>maxScore) maxScore=qual[i];
       }
       this.maxMLScore=maxScore;
       return maxScore;
   }

   private double getMachineLearnScore(int line) {
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
