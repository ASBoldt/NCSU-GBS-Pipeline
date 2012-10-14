/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome;

//import cern.colt.Arrays;
import java.util.Arrays;
import java.io.FileWriter;
import java.util.Random;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AllelePositionBLOBUtils;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.Pack1Alignment;
import net.maizegenetics.pal.datatype.DataType;

/**
 *
 * @author ed
 */
public class HaplotypeLengthDistribution {
    static byte UNKNOWN=(byte)DataType.UNKNOWN_CHARACTER;

    
    //TODO Find SNPs that break long haplotypes
    //TODO Impute missing data based on longest shared hit
    //TODO Make haplotype groups 

    public static int[] maxHaplotypeLength(Alignment align, int initialSite, int taxa1, int taxa2, boolean ignoreInitialSite) {
        int[] hapDescription=new  int[3];  //coded length, left start site, right start site
        //hapDescription[0]=-1;
        int b;
        byte s1b, s2b;
        boolean ibd;
        s1b=align.getBase(taxa1,initialSite);
        s2b=align.getBase(taxa2,initialSite);
        if(((s1b!=UNKNOWN)&&(s2b!=UNKNOWN))&&(ignoreInitialSite==false)) {
            if(s1b==s2b) {
                hapDescription[0]++;
                hapDescription[1]=hapDescription[2]=initialSite;
            } else return hapDescription;
        }
        ibd=true;
        for(b=initialSite-1; (b>=0)&&ibd; b--){
            s1b=align.getBase(taxa1,b);
            s2b=align.getBase(taxa2,b);
            if((s1b!=UNKNOWN)&&(s2b!=UNKNOWN)) {
                if(s1b==s2b) {
                    hapDescription[0]++;
                    hapDescription[1]=b;
                } else {ibd=false;}
            }
        }
        ibd=true;
        for(b=initialSite+1; (b<align.getSiteCount())&&ibd; b++){
            s1b=align.getBase(taxa1,b);
            s2b=align.getBase(taxa2,b);
            if((s1b!=UNKNOWN)&&(s2b!=UNKNOWN)) {
                if(s1b==s2b) {
                    hapDescription[0]++;
                    hapDescription[2]=b;
                } else {ibd=false;}
            }
        }
//        if (hapDescription[0]==-1) hapDescription[0]=0;

        return hapDescription;
    }


    public static int[] maxHaplotypeLength(Alignment align, int initialSite, int taxa1, int taxa2, int maxMismatch) {
        int[] hapDescription=new  int[3];  //coded length, left start site, right start site
        //hapDescription[0]=-1;
        int b, rightMissCnt=0, leftMissCnt=0;
        byte s1b, s2b;
        boolean ibd;
        int[][] leftMis=new int[2][maxMismatch+1];
        int[][] rightMis=new int[2][maxMismatch+1];
        ibd=true;
        for(b=initialSite; (b>=0)&&ibd; b--){
            s1b=align.getBase(taxa1,b);
            s2b=align.getBase(taxa2,b);
//            if(AllelePositionBLOBUtils.isBaseHomozygous(s1b)&&AllelePositionBLOBUtils.isBaseHomozygous(s2b)) {
            if((s1b!=UNKNOWN)&&(s2b!=UNKNOWN)) {
                if(s1b==s2b) {
                    leftMis[0][leftMissCnt]++;
                    //hapDescription[1]=b;
                } else {
                    leftMis[1][leftMissCnt]=b;
                    leftMissCnt++;
                    if(leftMissCnt>maxMismatch) {ibd=false;} else {
                        leftMis[0][leftMissCnt]=leftMis[0][leftMissCnt-1];
                    }
                }
            }
        }
        ibd=true;
        for(b=initialSite; (b<align.getSiteCount())&&ibd; b++){
            s1b=align.getBase(taxa1,b);
            s2b=align.getBase(taxa2,b);
//            if(AllelePositionBLOBUtils.isBaseHomozygous(s1b)&&AllelePositionBLOBUtils.isBaseHomozygous(s2b)) {
            if((s1b!=UNKNOWN)&&(s2b!=UNKNOWN)) {
                if(s1b==s2b) {
                    if(b!=initialSite) rightMis[0][rightMissCnt]++;
                    //hapDescription[2]=b;
                } else {
                    rightMis[1][rightMissCnt]=b;
                    rightMissCnt++;
                    if(rightMissCnt>maxMismatch) {ibd=false;}
                    else {rightMis[0][rightMissCnt]=rightMis[0][rightMissCnt-1];
                    }
                }
            }
        }
//        System.out.println("Left:"+Arrays.toString(leftMis[0])+Arrays.toString(leftMis[1]));
//        System.out.println("Right:"+Arrays.toString(rightMis[0])+Arrays.toString(rightMis[1]));
        for(int i=0; i<=maxMismatch; i++) {
            int length=leftMis[0][maxMismatch-i]+rightMis[0][i];
         //   if(leftMis[1][0]==initialSite) length=0;  //mismatch at the initial site
            if(hapDescription[0]<length) {
                hapDescription[0]=length;
                hapDescription[1]=leftMis[1][maxMismatch-i];
                hapDescription[2]=rightMis[1][maxMismatch-i];
            }
        }
        return hapDescription;
    }


    public static int[][] maxHaplotypeLengthMatrixPrint(Alignment align, int initialSite) {
        int[][] m=new int[align.getSequenceCount()][align.getSequenceCount()];
        for(int i=0; i<align.getSequenceCount(); i++) {
            System.out.print(align.getIdGroup().getIdentifier(i).getName()+"\t");
            for(int j=0; j<align.getSequenceCount(); j++) {
               m[i][j]=maxHaplotypeLength(align, initialSite, i, j, 1)[0];
               System.out.print(m[i][j]+"\t");
           }
           System.out.print("\n");
        }
        return m;
    }
    /**
     * Note in the return matrix I extend the length of the haplotype by 1 for all, so there is a difference between missing for zero.
     * @param align
     * @param initialSite
     * @param maxMismatch
     * @return
     */
     public static int[][] maxHaplotypeLengthMatrix(Alignment align, int initialSite, int maxMismatch) {
        int[][] m=new int[align.getSequenceCount()][align.getSequenceCount()];
        for(int i=0; i<align.getSequenceCount(); i++) {
            for(int j=0; j<i; j++) {
               m[j][i]=m[i][j]=maxHaplotypeLength(align, initialSite, i, j, maxMismatch)[0]+1;
               if(align.getBase(j, initialSite)==UNKNOWN) m[i][j]*=-1;
               if(align.getBase(i, initialSite)==UNKNOWN) m[j][i]*=-1;
           }
        }
        return m;
    }

       /**
     * Note in the return vector I extend the length of the haplotype by 1 for all, so there is a difference between missing for zero.
     * @param align
     * @param initialSite
     * @param maxMismatch
     * @return
     */
     public static int[] maxHaplotypeLengthMatrix(Alignment align, int initialSite, int seq, int maxMismatch) {
        int[] m=new int[align.getSequenceCount()];
        for(int i=0; i<align.getSequenceCount(); i++) {
                if(i==seq) {m[seq]=0;}
               else {m[i]=maxHaplotypeLength(align, initialSite, i, seq, maxMismatch)[0]+1;
                    if(align.getBase(i, initialSite)==UNKNOWN) m[i]*=-1;
                }
        }
        return m;
    }

    public static void lengthDistPerSite(Alignment align) {
        int skip=100;
         int propLength=align.getSiteCount();
         int[][] lengthDist=new int[align.getSequenceCount()+1][propLength/skip+1];
         int[][] freqDist=new int[align.getSequenceCount()+1][propLength/skip+1];
         for(int i=0; i<align.getSequenceCount(); i++) { //=align.getSequenceCount(); i++) {
            for(int j=0; j<align.getSequenceCount(); j++) {
                if(i==j) continue;
                for (int b = 0; b < propLength; b+=skip) {
             //       int[] hapDist=maxHaplotypeLength(align,b,i,j,false);
                    int[] hapDist=maxHaplotypeLength(align,b,i,j,1);
             //       System.out.println("initbase="+b+" oldstyle="+hapDist[0]+" newstyle="+hapDist2[0]);
                    if(hapDist[0]>lengthDist[i][b/skip])
                        {lengthDist[i][b/skip]=hapDist[0];
                     //   System.out.println("initbase="+b+" oldstyle="+hapDist[0]+" newstyle="+hapDist2[0]);
                        }
                    if(hapDist[0]>30) freqDist[i][b/skip]++;
            //        int[] hapDistIgnore=maxHaplotypeLength(align,b,i,j,true);
            //        if((hapDistIgnore[0]-hapDist[0])>lengthDist[i][b]) {lengthDist[i][b]=hapDistIgnore[0]-hapDist[0];}

                }
            }
         }
        StringBuilder sb=new StringBuilder();
        sb.append("Site:Line\t");
        for(int i=0; i<align.getSequenceCount(); i++) {sb.append(align.getIdGroup().getIdentifier(i).getName()+"\t");}
        sb.append("RandomTaxa\n");
        for(int j=0; j<lengthDist[0].length; j++) {
            sb.append(j+"\t");
            for(int i=0; i<=align.getSequenceCount(); i++) {sb.append(lengthDist[i][j]+"\t");}
            sb.append("\n");
        }
//        System.out.println(sb.toString());
//        sb=new StringBuilder();
//        for(int i=0; i<align.getSequenceCount(); i++) {sb.append(align.getIdGroup().getIdentifier(i).getName()+"\t");}
//        sb.append("RandomTaxa\n");
//        for(int j=0; j<freqDist[0].length; j++) {
//            sb.append(j+"\t");
//            for(int i=0; i<=align.getSequenceCount(); i++) {sb.append(freqDist[i][j]+"\t");}
//            sb.append("\n");
//        }
       // System.out.println(sb.toString());
        try{
 //           FileWriter fw=new FileWriter("E:/SolexaAnal/test/chr8Hap30bpLengthDistImpPerm.txt");
            FileWriter fw=new FileWriter("/Users/edbuckler/SolexaAnal/HapMapV2/chr10Hap30bpLengthDistMM1perm.txt");
            fw.write(sb.toString());
            fw.close();

        } catch (Exception e) {
            System.out.println("Error writing HapFreqDist.txt"+e);
        }

    }

    public static int[] evaluateImputation(Pack1Alignment align, int everyXsiteImpute, int maxMismatch) {
        int[] impAccTest=new int[2];
        //copy of initial data so it can set to missing
        byte[][] impAlleleBLOBs=new byte[align.getSequenceCount()][];
        for(int i=0; i<align.getSequenceCount(); i++) {
            impAlleleBLOBs[i]=(byte[])align.getAlleleBLOBs(i).clone();
        }
        //record the answers and then set to missing.
        int maxSite=align.getSiteCount();//10000;
      //  int maxSite=10000;
        int totalTestSites=maxSite/everyXsiteImpute+1;
        int count=0;
        int[] testSites=new int[totalTestSites], testLines=new int[totalTestSites];
        byte[] answers=new byte[totalTestSites];
        for (int b = 0; b <maxSite; b+=everyXsiteImpute) {
            testSites[count]=b;
           // int i=count%align.getSequenceCount();
            testLines[count]=count%align.getSequenceCount();
        //    testLines[count]=count%32;  //just test the maize
         //   for(int i=0; i<align.getSequenceCount(); i++) {
                answers[count]=AllelePositionBLOBUtils.getBaseFromAlleleBLOB(impAlleleBLOBs[testLines[count]],b);
                AllelePositionBLOBUtils.setHalfByteInAlleleBLOB(impAlleleBLOBs[testLines[count]],b,(char)UNKNOWN);
           // }
            count++;
        }

        Pack1Alignment testP1A=new Pack1Alignment(impAlleleBLOBs,align.getVariableSitesBLOB(),align.getSNPidBLOB());
//        align=null;  //allow to be garbage collected
        //impute the sites using the standard method
        Pack1Alignment impP1A=imputeBySite(testP1A,maxMismatch);
        //evaluate imputation
        int testedBases=0, correctBases=0;
        //count=0;
        for(int t=0; t<totalTestSites; t++) {
           // for(int i=0; i<impP1A.getSequenceCount(); i++) {
                if((answers[t]!=UNKNOWN)&&AllelePositionBLOBUtils.isBaseHomozygous(answers[t])&&
                        AllelePositionBLOBUtils.isBaseHomozygous(impP1A.getBase(testLines[t], testSites[t]))) {
                   testedBases++;
                   if(answers[t]==impP1A.getBase(testLines[t], testSites[t])) correctBases++;
                   System.out.print("Site:"+testSites[t]+" Line:"+testLines[t]+" answer:"+(char)answers[t]+" imp:"+(char)impP1A.getBase(testLines[t], testSites[t]));
                   System.out.println(" CommonBase:"+impP1A.getSiteSummary(testSites[t]).getAlleles()[0]+
                           " MAF:"+impP1A.getSiteSummary(testSites[t]).getAlleleFrequency(0));
                }
           // }
        }
        System.out.println("Imputed testedBases:"+testedBases+" correctBases:"+correctBases);
        impAccTest[0]=testedBases; impAccTest[1]=correctBases;
        return impAccTest;
    }

    /**
     * Mask a site to unknown and then attempts to impute it.  It plays the bad form game
     * of playing with the BLOBs in align, but this is the only fast way to do things millions of
     * times.
     * @param align
     * @param maxMismatch
     * @param seq
     * @param site
     * @return
     */
    public static byte imputeOneSite(Pack1Alignment align, int maxMismatch, int seq, int site) {
        byte answer=align.getBase(seq, site);
        AllelePositionBLOBUtils.setHalfByteInAlleleBLOB(align.getAlleleBLOBs(seq),site,(char)UNKNOWN);
        int[] dm=maxHaplotypeLengthMatrix(align, site, seq, maxMismatch);  //fast to just calc the vector
    //    int[] dm=new int[align.getSequenceCount()];
        int maxMatch=0, bestLine=-1;
        for(int j=0; j<align.getSequenceCount(); j++) {
            if(dm[j]>maxMatch) {maxMatch=dm[j]; bestLine=j;}
        }
        byte impAnswer=UNKNOWN;
        if(bestLine>-1) impAnswer=align.getBase(bestLine, site);
        AllelePositionBLOBUtils.setHalfByteInAlleleBLOB(align.getAlleleBLOBs(seq),site,(char)answer);
        return impAnswer;
    }

    public static void summarizeImputationBySite(Pack1Alignment align, int maxMismatch, String outfile) {
        int totalKnown=0, totalCorrect=0;
        try{
        FileWriter fw=new FileWriter(outfile);
        fw.write("site\tMinorAlleleFrequency\tknown\tcorrectImp\n");
        int propLength=align.getSiteCount();//1000;
        for (int site = 0; site < propLength; site++) {
            int known=0, correctImp=0;
            for(int i=0; i<align.getSequenceCount(); i++) {
                byte answer=align.getBase(i, site);
                if(AllelePositionBLOBUtils.isBaseHomozygous(answer)) {
                    known++;
                    byte impAnswer=imputeOneSite(align, maxMismatch, i, site);
                    if(answer==impAnswer) correctImp++;
                }
            }
//            System.out.println(site+"\t"+align.getMinorAlleleFrequency(site)+"\t"+known+"\t"+correctImp);
            fw.write(site+"\t"+align.getPositionInLocus(site)+"\t"+align.getMinorAlleleFrequency(site)+"\t"+known+"\t"+correctImp+"\n");
        }
        fw.close();
            System.out.println("Known %d");
        } catch (Exception e) {
            System.out.println("Error writing summarizeImputationBySite"+e);
        }
    }

    public static Pack1Alignment imputeBySite(Pack1Alignment align, int maxMismatch) {
        int propLength=align.getSiteCount();
        byte[][] impAlleleBLOBs=new byte[align.getSequenceCount()][];
        for(int i=0; i<align.getSequenceCount(); i++) {
            impAlleleBLOBs[i]=(byte[])align.getAlleleBLOBs(i).clone();
        }
        int knownSNPs=0, unknownSNPs=0;
         for (int b = 0; b < propLength; b++) {
            if(b%1000==0)
                {System.out.println("Imputed base:"+b+" known:"+knownSNPs+" unknownSNPs:"+unknownSNPs);}
            int[][] dm=maxHaplotypeLengthMatrix(align, b, maxMismatch);
            for(int i=0; i<align.getSequenceCount(); i++) {
                 byte focusBase=align.getBase(i, b);
                 if(focusBase==UNKNOWN) {
                    int maxMatch=0, bestLine=-1;
                    for(int j=0; j<align.getSequenceCount(); j++) {
                        if(dm[i][j]>maxMatch) {maxMatch=dm[i][j]; bestLine=j;}
                    }
                    if(bestLine==-1) {
                        System.out.println("Error Imputed base:"+b+" focusBase"+(char)focusBase+" known:"+knownSNPs+" unknownSNPs:"+unknownSNPs);
                    }
                    unknownSNPs++;
                    if(maxMatch>0) {AllelePositionBLOBUtils.setHalfByteInAlleleBLOB(impAlleleBLOBs[i],b,
                            align.getBaseChar(bestLine, b));}
                    else {
                        AllelePositionBLOBUtils.setHalfByteInAlleleBLOB(impAlleleBLOBs[i],b,
                            align.getSiteSummary(b).getAlleles()[0]);
                    }
                }  else {knownSNPs++;}
            }
         }
        Pack1Alignment newAlign=new Pack1Alignment(impAlleleBLOBs,align.getVariableSitesBLOB(),align.getSNPidBLOB());
        return newAlign;
    }

    public static void createHaplotypeFile(Alignment align, int minHapLength, int maxMismatch, String outfile) {
         int propLength=align.getSiteCount();
         int[] hapCntSum=new int[align.getSequenceCount()];
         byte[][] hapName=new byte[align.getSequenceCount()][propLength];
         int[] cntHaps=new int[propLength];
         for (int b = 0; b < propLength; b++) {
             for(int i=0; i<align.getSequenceCount(); i++) {hapName[i][b]=(byte)i;}
             int[][] dm=maxHaplotypeLengthMatrix(align, b, maxMismatch);
             for(int i=0; i<align.getSequenceCount(); i++) {
                for(int j=i+1; j<align.getSequenceCount(); j++) {
                    if(dm[i][j]>minHapLength) hapName[j][b]=hapName[i][b];
                }
            }
            int[] hapCnt=new int[align.getSequenceCount()];
            for(int i=0; i<align.getSequenceCount(); i++) {hapCnt[hapName[i][b]]++;}
            Arrays.sort(hapCnt);
            for(int i=0; i<align.getSequenceCount(); i++) {
                if (hapCnt[i]>0) cntHaps[b]++;
                hapCntSum[i]+=hapCnt[i];
            }
           //System.out.println(b+"\t"+Arrays.toString(hapCntSum));
         }
         System.out.println(align.getSiteCount()+"\t"+Arrays.toString(hapCntSum));
        StringBuilder sb=new StringBuilder();
        sb=new StringBuilder();
        sb.append("site\tchr\tposition\t");
        for(int i=0; i<align.getSequenceCount(); i++) {sb.append(align.getIdGroup().getIdentifier(i).getName()+"\t");}
        sb.append("hapcount\t");
        sb.append("\n");
        
        try{
            FileWriter fw=new FileWriter(outfile);
            fw.write(sb.toString());
            sb=new StringBuilder();
            for(int j=0; j<hapName[0].length; j++) {
                sb.append(j+"\t"+align.getLocusName(j)+"\t"+align.getPositionInLocus(j)+"\t");
                for(int i=0; i<align.getSequenceCount(); i++) {sb.append(hapName[i][j]+"\t");}
                sb.append(cntHaps[j]+"\t");
                sb.append("\n");
                fw.write(sb.toString());
                sb=new StringBuilder();
            }
            fw.close();
        } catch (Exception e) {
            System.out.println("Error writing HapFreqDist.txt"+e);
        }

    }


    public static void lengthDistribution(Alignment align) {
        int[][] lengthDist=new int[align.getSequenceCount()+1][1000];
        Random generator = new Random();
      //  int i=0;
        int count=0;
        for(int i=0; i<=align.getSequenceCount(); i++) {
            for(int j=0; j<align.getSequenceCount(); j++) {
                if(i==j) continue;
                for (int b = 0; b < align.getSiteCount(); b++) {
                    byte s1b;//=p1a.getBase(i,b);
                    byte s2b=align.getBase(j,b);
                    if(i==align.getSequenceCount()) {
                        int randTaxa;
                        do{randTaxa=generator.nextInt(align.getSequenceCount());}
                        while(randTaxa==j);
                        s1b=align.getBase(randTaxa,b);
                    } else {
                        s1b=align.getBase(i,b);
                    }
                    if((s1b!=(byte)DataType.UNKNOWN_CHARACTER)&&(s2b!=(byte)DataType.UNKNOWN_CHARACTER)) {
                        if(s1b==s2b) {count++;}
                        else {
                            if(count>=lengthDist[i].length) {lengthDist[i][lengthDist.length-1]++;} else {lengthDist[i][count]++;}
                            count=0;
                        }
                    }
                }
//                System.out.println(i+" "+j+" "+count);
            }
        }
        System.out.println(Arrays.toString(lengthDist));
        StringBuilder sb=new StringBuilder();
        for(int i=0; i<align.getSequenceCount(); i++) {sb.append(align.getIdGroup().getIdentifier(i).getName()+"\t");}
        sb.append("RandomTaxa\n");
        for(int j=0; j<lengthDist[0].length; j++) {
            sb.append(j+"\t");
            for(int i=0; i<=align.getSequenceCount(); i++) {sb.append(lengthDist[i][j]+"\t");}
            sb.append("\n");
        }
        System.out.println(sb.toString());
    }

    /**
     *
     * @param args
     * @throws java.lang.Exception
     */
    public static void main(String[] args) throws Exception {
        String basePath=System.getProperty("user.dir")+"/maizegenetics/test/datafiles/"; 
//        String infile=basePath+"all_8pos.zip";
    //    String infile=basePath+"all_8posNoCML52.zip";
        String infile="/Users/edbuckler/SolexaAnal/HapMapV2/anchor/h100kchr10.genotypes.log2ct_h90.good.hmp.txt";
     //   String infile="/Users/edbuckler/SolexaAnal/HapMapV2/anchor/chr10.R1_BWT_NV.small.log2_h90.good.hmp.txt";

        /**
         * Change the hets to missing as all sort of this class are having problems with this.
         *
         * Also need to examine the LD of the missing data
         */

        Pack1Alignment p1a = (Pack1Alignment) ImportUtils.createPack1AlignmentFromFile(infile, "", "");
        System.out.println("Alignment Loaded:"+infile);
        System.out.printf("Sites %d Taxa %d %n:", p1a.getSiteCount(),p1a.getSequenceCount());
//        IdGroup idWithoutCML54=new SimpleIdGroup(p1a.getIdGroup(),8);
//        Pack1Alignment p1aV2 = AllelePositionBLOBUtils.subsetCopyAlignment(p1a, idWithoutCML54);
//        GdpdmBLOBUtils.writePack1AlignmentToZip(p1aV2, infile2);
      //  lengthDistribution(p1a);
        p1a=AllelePositionBLOBUtils.permuteAlignment(p1a);
        lengthDistPerSite(p1a);
//        maxHaplotypeLengthMatrix(p1a,0,1);
//        Pack1Alignment newP1A=imputeBySite(p1a, 1);
//        GdpdmBLOBUtils.writePack1AlignmentToZip(newP1A, "C:/EdStuff/Solexa/test/chr8Imp.zip");
 //       GdpdmBLOBUtils.writePack1AlignmentToZip(newP1A, "E:/SolexaAnal/test/chr8Imp.zip");
    //    
  //      evaluateImputation(p1a,30,1);
 //       summarizeImputationBySite(p1a,4,"/Users/edbuckler/SolexaAnal/HapMapV2/test/chr10ImpByErrorRateMis0.txt");
//        summarizeImputationBySite(p1a,0,"E:/SolexaAnal/test/chr8ImpByErrorRateMis0.txt");
        //Create haplotype file from imputed data
        String impInfile="C:/EdStuff/Solexa/test/chr8Imp.zip";
//        String impInfile="E:/SolexaAnal/test/chr8Imp.zip";
//        Pack1Alignment p1a = (Pack1Alignment) ImportUtils.createPack1AlignmentFromFile(impInfile, "", "");
    //    maxHaplotypeLengthMatrix(p1a,66696);
      String outfile="/Users/edbuckler/SolexaAnal/0_10.snpsindels.combined.1.13.log2ml.HP1/1hap.txt";

 //      createHaplotypeFile(newP1A,30,1,outfile);
  //      createHaplotypeFile(p1a,30,1,"E:/SolexaAnal/test/chr8Hap30bpHapNameNum.txt");
    }
}
