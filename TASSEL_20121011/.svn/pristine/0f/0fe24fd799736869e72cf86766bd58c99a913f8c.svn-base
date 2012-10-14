/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.pal.alignment;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;
import cern.colt.list.IntArrayList;
import cern.colt.list.ShortArrayList;
import java.nio.IntBuffer;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * These utils are used for the creation create a set of coverage BLOBs for genotypic data including methods for:
 * 1. converting a Jer-Ming coverage text file to zipped coverageBLOBs
 * 2. subsetting a Jer-Ming coverage text file based on lines covered
 * 3. methods for getting coverage for a single bases or for a range of bases
 * 
 *
 * @author ed
 */
public class CoverageBLOBUtils {
    private static final char UNKNOWN_CHARACTER='N';  //Todo change to DataType.UNKNOWN_CHARACTER when fixed
    /**
     * Constants that decribe the GDPDM blobs used for storing diversity data
     * Field variables define the start and end position in the blob for the variables
     */
  ;


    //Jer-ming Chia coverage format constants
    public static final int totalNonTaxaHeaders=4;


    private CoverageBLOBUtils() {}


//        private static void writeBLOBtoZip(byte[][] coverageBLOB, String outfileName) {
//         try {
//            FileOutputStream fos = new FileOutputStream(outfileName);
//            java.util.zip.ZipOutputStream zos = new ZipOutputStream(fos);
//            zos.setLevel(5);
//            for(int i=0; i<coverageBLOB.length; i++) {
//                String s=new String(coverageBLOB[i],GdpdmBLOBUtils.taxaField[0],GdpdmBLOBUtils.taxaField[1]);
//                String[] taxaParts=s.trim().split(":");
////       System.out.println(s+" size="+alleleBLOB[i].length);
//                String locusName=(new String(coverageBLOB[i],GdpdmBLOBUtils.locusField[0],GdpdmBLOBUtils.locusField[1])).trim();
//                ZipEntry thisEntry = new ZipEntry(taxaParts[0]+"_"+locusName+".cvb");
//                zos.putNextEntry(thisEntry);
//                zos.write(coverageBLOB[i]);
//                zos.closeEntry();
//            }
//            fos.close();
//
//        } catch (IOException ee) {
//            System.err.println("Data could not be saved: " + ee);
//        }
//     }

    /**
     * Create Pack1Alignment from the given infile.  The infile should be a zipped file for allele
     * blobs and one position blob
     * @param infileName
     */
/*    public static byte[][] getCoverageFromZipBLOB(String infileName) {
        byte[][] coverageB=null; //array for BLOBs holding SNPs plus header [taxa][SNPs]
        Pack1Alignment p1a=null;
        try {
          System.out.println("Reading:"+infileName);
          FileInputStream fis = new FileInputStream(infileName);
          ZipInputStream zis = new ZipInputStream(fis);
          ZipEntry ze;
          int numTaxa=0;
          while((ze=zis.getNextEntry())!=null){
            if(ze.getName().contains(".cvb")) {
                numTaxa++;
            }
          }
          fis.close();
          coverageB=new byte[numTaxa][];
          fis = new FileInputStream(infileName);
          zis = new ZipInputStream(fis);
          numTaxa=0;
          ByteBuffer bb=ByteBuffer.allocate(10000000);
          byte[] buf=new byte[40960];
          while((ze=zis.getNextEntry())!=null){
            bb.clear();
            //System.out.println(ze.getName());
            int c;
                 while ((c=zis.read(buf)) != -1) {
                    bb.put(buf, 0, c);
                }
            int entrySize=bb.position();
           // System.out.println(entrySize);
            if(ze.getName().contains(".cvb")) {
                coverageB[numTaxa]=new byte[entrySize];
                bb.position(0);
                bb.get(coverageB[numTaxa]);
//      System.out.println(new String(alleleB[numTaxa]));
                numTaxa++;
            }
            zis.closeEntry();
          }
          zis.close();
          System.out.println("Coverage data zip read created");
        }catch (IOException ee) {
            System.err.println("Coverage Data in zip could not be read: " + ee);
        }
        return coverageB;
    }

*/
        /**
     * Count the number of lines in a file
     * @param infileName
     * @return lines
     */
    public static int subsetRowsInCoverageFile(String infileName, String outfileName, int minTaxa, boolean addHeader) {
        int lines=0;
        try{
            BufferedReader fileIn   = new BufferedReader(new FileReader(infileName),1000000);
            BufferedWriter fileOut   = new BufferedWriter(new FileWriter(outfileName),1000000);
            String s;
            while ((s=fileIn.readLine()) != null){
                String[] sp=s.split("\t");
                int count=0;
                for(int i=4; i<sp.length; i++) {
                    if(Integer.parseInt(sp[i])>0) count++;
                }
                if(count>minTaxa) fileOut.write(s+"\n");
                lines++;
            }
            fileIn.close();
            fileOut.close();
        } catch (Exception e) {
            System.err.println("File IO in countLinesInFile: " + e );
        }
        return lines;
    }

    //"B73:temperate:mays:mays:Zea	B97:temperate:mays:mays:Zea	CML103:tropical:mays:mays:Zea	CML228:tropical:mays:mays:Zea	CML247:tropical:mays:mays:Zea	CML277:tropical:mays:mays:Zea	CML322:tropical:mays:mays:Zea	CML333:tropical:mays:mays:Zea	CML52:BC1:mays:mays:Zea	CML52R:tropical:mays:mays:Zea	CML69:tropical:mays:mays:Zea	HP301:popcorn:mays:mays:Zea	IL14H:sweet:mays:mays:Zea	KI11:tropical:mays:mays:Zea	KI3:tropical:mays:mays:Zea	KY21:temperate:mays:mays:Zea	M162W:temperate:mays:mays:Zea	M37W:tropical:mays:mays:Zea	MO17:temperate:mays:mays:Zea	MO18W:tropical:mays:mays:Zea	MS71:temperate:mays:mays:Zea	NC350:tropical:mays:mays:Zea	NC358:tropical:mays:mays:Zea	OH43:temperate:mays:mays:Zea	OH7B:temperate:mays:mays:Zea	P39:sweet:mays:mays:Zea	TX303:tropical:mays:mays:Zea	TZI8:tropical:mays:mays:Zea"

    /**
     * Create a coverage file from text coverage file.  The first show should have .
     * @param infileName name of text infile with coveage information
     * @param outfileName name of the compressed file for writing to.
     */
    public static ArrayList<byte[]> createCoverageBLOBs(String infileName) {
        int lines=0;
        IntArrayList[] startPosList;
        ShortArrayList[] lengthList;
        ArrayList<byte[]> theBLOBs=new ArrayList<byte[]>();
        String[] taxaNames;
        try{
            BufferedReader fileIn   = new BufferedReader(new FileReader(infileName),1000000);
            String[] spp=fileIn.readLine().split("\t");
            lines=spp.length-4;
            taxaNames=new String[lines];
            startPosList=new IntArrayList[lines];
            lengthList=new ShortArrayList[lines];
            for(int i=0; i<lines; i++) {
                startPosList[i]=new IntArrayList();
                lengthList[i]=new ShortArrayList();
                taxaNames[i]=spp[i+4];
            }
//            fileIn.close();
//            fileIn   = new BufferedReader(new FileReader(infileName),1000000);
         //   BufferedWriter fileOut   = new BufferedWriter(new FileWriter(outfileName),1000000);
            String s;
            while ((s=fileIn.readLine()) != null){
               // System.out.println(s);
                String[] sp=s.split("\t");
                if(sp.length<lines) continue;
               // int count=0;
                int start=Integer.parseInt(sp[1]);
        if(start%100000==0) System.out.println(s);
                int end=Integer.parseInt(sp[2]);
                short size=(short)(end-start+1);
                for(int i=0; i<lines; i++) {
                    if(Integer.parseInt(sp[i+4])>0) {
                        int lastElement=startPosList[i].size()-1;
                   // System.out.println("lastElement="+lastElement);
                        if((lastElement>-1)&&(start==startPosList[i].getQuick(lastElement)+lengthList[i].getQuick(lastElement))) {
                            //adjacent
                            short newLength=(short)(lengthList[i].getQuick(lastElement)+size);
                            if(newLength<255)
                            {lengthList[i].setQuick(lastElement, newLength);}
                            else {
                               startPosList[i].add(start);
                               lengthList[i].add(size);
                            }
                        } else {
                            startPosList[i].add(start);
                            lengthList[i].add(size);
                        }
                    }
                }
            }
            fileIn.close();
            //coverageBLOB=new byte[lines][];
            for(int i=0; i<lines; i++) {
                System.out.println("Line:"+taxaNames[i]+" size="+startPosList[i].size());
                byte[] theBLOB=new byte[GdpdmBLOBUtils.totalHeaderPadding+(startPosList[i].size()*8)];
             //   String taxonName="undefined"+i;
                GdpdmBLOBUtils.setHeaderOfBLOB(theBLOB, GdpdmBLOBUtils.compressionAlgorithm, GdpdmBLOBUtils.coverageBLOBtype,
                        startPosList[i].size(), "undefined", "locus", startPosList[i].getQuick(0),
                        startPosList[i].getQuick(startPosList[i].size()-1), taxaNames[i]);
                ByteBuffer cBB=ByteBuffer.wrap(theBLOB);
                cBB.position(GdpdmBLOBUtils.totalHeaderPadding);
                for (int j = 0; j < startPosList[i].size(); j++) {
                    cBB.putInt(startPosList[i].getQuick(j)); //start position
                    cBB.putInt(startPosList[i].getQuick(j)+lengthList[i].getQuick(j)-1);  //end position
                }
                theBLOBs.add(theBLOB);
            }
            return theBLOBs;
            //fileOut.close();
        } catch (Exception e) {
            System.err.println("File IO in createCoverageBLOBs: " + e );
            return null;
        }
    }

    /**
     * Determines whether the given site is covered in the sequencing
     * @param startEndArray - must be ordered by StartOfRegion1 EndOfRegion1 StartOfRegion2 EndOfRegion2
     * @param position
     * @return true is the site was in the range of covered sites, false otherwise
     */
    public static boolean isCovered(int[] startEndArray, int position) {
        if(position<startEndArray[0]||position>startEndArray[startEndArray.length-1]) return false;
        int r=Arrays.binarySearch(startEndArray, position);
        if(r>0) return true;  //direct hit on start or end so covered
        r=-r;
        if(r%2==0) return true;  //if the insert site is even the it must fall between start and end
        return false;
    }

    
      /**
     *
     * @param args
     * @throws java.lang.Exception
     */
    public static void main(String[] args) throws Exception {
         int chr=8;
           String coverageInputFile="C:/EdStuff/Solexa/test/"+chr+"_12_coverageh.txt", testDirectory="C:/EdStuff/Solexa/test/";
//           String hapmapInputFile="E:/SolexaAnal/AGPLog2ML/combined1_13HP/"+chr+".snps.haploid.1.13.log2ml.HP1.txt", testDirectory="E:/SolexaAnal/test/";
     //   subsetRowsInCoverageFile("C:/EdStuff/Solexa/AGPCoverage/maxmap1/8.coverage.1.combined.txt","C:/EdStuff/Solexa/test/8_12_coverage.txt", 12, true);
//        createCoverageBLOBs("C:/EdStuff/Solexa/test/8_12_coverage.txt","C:/EdStuff/Solexa/test/junk1.zip");
       // byte[][] b=new byte[10][10];
                //getCoverageFromZipBLOB("C:/EdStuff/Solexa/test/junk1.zip");
//        ArrayList<byte[]> theBLOBs=createCoverageBLOBs(coverageInputFile);
//        GdpdmBLOBUtils.writeBLOBtoFiles(theBLOBs, testDirectory+"sepf/");
//        GdpdmBLOBUtils.writeBLOBtoZip(theBLOBs, testDirectory + "all_" +chr+ "cov.zip");
//        GdpdmBLOBUtils.writeBLOBtoLZMA(theBLOBs, testDirectory + "all_" +chr+ "cov.7z");
/*        ArrayList<byte[]>[] theBLOBList = new ArrayList[GdpdmBLOBUtils.totalBLOBtypes];
        for (int i = 0; i < theBLOBList.length; i++) {
            theBLOBList[i] = new ArrayList<byte[]>();
        }
        theBLOBList=GdpdmBLOBUtils.readBLOBfromLZMA(testDirectory + "all_" +chr+ "cov.7z", theBLOBList);
        byte[] b=theBLOBList[2].get(0);
        ByteBuffer bb=ByteBuffer.wrap(b);
        bb.position(GdpdmBLOBUtils.numSitesField[0]);
        int segments=bb.getInt();
        IntBuffer ib=ByteBuffer.wrap(b,GdpdmBLOBUtils.totalHeaderPadding,
                    b.length-GdpdmBLOBUtils.totalHeaderPadding).asIntBuffer();
        int[] ii=new int[segments*2];
        for(int i=0; i<ii.length; i++) {
             ii[i]=ib.get();
        }
        int r=0;
        long time=System.currentTimeMillis();
        for(int i=0; i<1000000; i++) {
             boolean bc=isCovered(ii,i);
          //  if(bc) System.out.println(i);
        }
        long time2=System.currentTimeMillis();
        System.out.println(r+" "+(time2-time));

*/
        int sitesToPull=2600;
        Pack1Alignment p1a=(Pack1Alignment)ImportUtils.createPack1AlignmentFromFile("C:/EdStuff/Solexa/test/all_8pos.7z", "C:/EdStuff/Solexa/test/all_8cov.7z", "");
        for (int j = 0; j < 1; j++) {
            for (int i = 0; i < p1a.getSequenceCount(); i++) {
                // p1a.getAlignedSequenceString(i);
                // System.out.println(p1a.getIdentifier(i).getName()+"  "+ p1a.getAlignedSequenceString(i).substring(0, 10));
                // System.out.println(p1a.getIdentifier(i).getName()+"  "+ p1a.getAlignedSequenceString(i).substring(p1a.getSiteCount()-10));
                byte[] allele = new byte[sitesToPull];
                 for(int bp=0; bp<sitesToPull; bp++) {allele[bp]=p1a.getBase( i, bp);}

                //
               // allele = p1a.getDataRange(i, 0, p1a.getSiteCount() - 1);

                 System.out.print(p1a.getIdGroup().getIdentifier(i).getName());
                System.out.println(new String(allele, 0, 1000));
              //  System.out.print("  ");
               // System.out.println(new String(allele, allele.length - 10, 10));
            }
        }
    }
}
