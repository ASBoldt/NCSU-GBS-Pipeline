/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.gbs.pipeline;

import java.awt.Frame;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import javax.swing.ImageIcon;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.gbs.tagdist.AbstractTags;
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.gbs.tagdist.Tags;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.util.ArgsEngine;
import net.maizegenetics.genome.BaseEncoder;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.util.DirectoryCrawler;
import org.apache.log4j.Logger;

/* Implements an external mergesort to combine multiple tag-count files.
 * @author edbuckler
 */
public class CompareMultipleTagCountPlugin extends AbstractPlugin {

    private static long[][] ctags;  //array of current tags, [indexOfFile][indexOfTagLong], emptyLong=0, fileFinish=Long.MAX
    private static int[] ctagCnt; //array of tag counts [indexOfFile]
    private static byte[] ctagLength;  //array of tag length [indexOfFile]
    private static int[] chunkTagSizes;  // array of number of tags in each file  [indexOfFile]
    private static boolean outputSeqStr = false;
    private static ArrayList <TagWithCounts> tagCountsMatrix = new ArrayList<TagWithCounts>();
    private static int[] readCountPerFile;
    private static int[][] tagDistPerFile;
    private static int maxDepth = 1000;
    
//    private int[] ctagRemaining; //How many tags still remain to be read from each chunk file
    private static int tagLengthInLong=2;
    private static int tagsRead=0, outCnt=0;
    private static int rwOpen=0;
    private static DataInputStream[] rw;
    private static ArgsEngine myArgsEngine = null;
        static     File inputDirectory=null;
       static  String[] inputFileNames=null;
       static  String outputFileName = null;
       static  int minCount = 0;
        private static final Logger myLogger=Logger.getLogger(CompareMultipleTagCountPlugin.class);

    public CompareMultipleTagCountPlugin() {
        super(null, false);
    }

    public CompareMultipleTagCountPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    @Override
    public void setParameters(String[] args){
        if(args.length==0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        
        if (myArgsEngine == null) {
            myArgsEngine = new ArgsEngine();
            myArgsEngine.add("-i", "--input-directory", true);
            myArgsEngine.add("-o", "--output-directory", true);
            myArgsEngine.add("-c", "--min-count", true);
            myArgsEngine.add("-t", "--fastq");
        }

        myArgsEngine.parse(args);
        if (myArgsEngine.getBoolean("-o")) {
            outputFileName = myArgsEngine.getString("-o");
        }else{
            printUsage();
            throw new IllegalArgumentException("Please specify an output file.");
        }

        if (myArgsEngine.getBoolean("-i")) {
            inputDirectory = new File(myArgsEngine.getString("-i"));
            if (!inputDirectory.isDirectory()) {
                printUsage();
                throw new IllegalArgumentException("The input name you supplied is not a directory.");
            }
            inputFileNames = DirectoryCrawler.listFileNames(".*\\.cnt", inputDirectory.getAbsolutePath());
            if(inputFileNames == null || inputFileNames.length == 0){
                printUsage();
                throw new IllegalArgumentException("Couldn't find any files ending in \".cnt\" in the directory you specified.");
            }
            myLogger.info("Merging the following .cnt files...");
            for (String filename : inputFileNames) {
                myLogger.info(filename);
            }
            myLogger.info("...to \""+outputFileName+"\"."   );
        }else{
            printUsage();
            throw new IllegalArgumentException("Please specify an input file.\n");

        }

        if(myArgsEngine.getBoolean("-c")){
            minCount = Integer.parseInt(myArgsEngine.getString("-c"));
        }else{
            minCount =1;
        }

    }


private void printUsage(){
            myLogger.info(
                    "\nUsage is as follows:\n"
                    + "-i  Input directory containing .cnt files\n"
                    + "-o  Output directory name\n"
                    + "-c Minimum summed depth of reads to be output (default 1)\n"
            );
}
 
    @Override
    public DataSet performFunction(DataSet input) {
        mergeChunks(inputFileNames, inputDirectory, outputFileName, minCount);
        return null;
    }
    

    public  void mergeChunks(String[] chunkFileNames, File inputDirectory, String outputFileName, int minCount) {
        rw=new DataInputStream[chunkFileNames.length];
        readCountPerFile = new int[chunkFileNames.length];
        tagDistPerFile=new int[chunkFileNames.length][maxDepth+1];
        chunkTagSizes=new int[rw.length];
        ctagCnt=new int[rw.length];
        ctagLength=new byte[rw.length];
        rwOpen=rw.length;
        
        Arrays.fill(readCountPerFile, 0);
        try {
            for (int f = 0; f < rw.length; f++) {
                String infile=chunkFileNames[f];
                rw[f]=new DataInputStream(new BufferedInputStream(new FileInputStream(infile), 4000000));
                chunkTagSizes[f]=rw[f].readInt();
                tagLengthInLong=rw[f].readInt();
                myLogger.info("Opened :"+infile+" tags="+chunkTagSizes[f]);
            }
//            outStream = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFileName+".fq"), 655360));
            ctags=new long[rw.length][tagLengthInLong];
            outStack("B:");
            int t=0;
            while(rwOpen>0) {
                long[] minTag=updateCurrentTags();
           //     myLogger.info(BaseEncoder.getSequenceFromLong(minTag));
           //     outStack("U:"+t);
                writeTags(minTag, minCount);
           //     outStack("A:"+t);
                if(t%1000000==0) {
                    System.out.printf("t=%d tagsRead=%d outCnt=%d rwOpen=%d %n",t,tagsRead,outCnt,rwOpen);
                    outStack("AA:"+t);
                    myLogger.info(BaseEncoder.getSequenceFromLong(minTag));
                    }
                t++;
            }
//            outStream.flush();
//            outStream.close();
        } catch (Exception e) {
            myLogger.info("Catch in reading TagCount file e=" + e);
            e.printStackTrace();
        }
        
        //start sorting
        System.out.println("start sorting ...");
        Collections.sort(tagCountsMatrix, new CompareIntArray());

        //String filenamepattern = ".+[\\/]";
        try
        {
            //write matrix
            DataOutputStream outStream = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFileName + "\\matrix.txt"), 655360));
//            if (outputSeqStr)
//            {
//                outStream.writeBytes("TAG\tLength\tSum");
//            }
//            else
//            {
                outStream.writeBytes("Tag\tLength\tSum");
//            }
            for (int j=0; j<chunkFileNames.length; j++)
            {
                outStream.writeBytes("\t" +  removePath(chunkFileNames[j]));
            }
            outStream.writeBytes("\n");
            for (int i=0; i<tagCountsMatrix.size(); i++)
            {
                String tagSequence = BaseEncoder.getSequenceFromLong(tagCountsMatrix.get(i).tagseq);
                int tagLen = tagCountsMatrix.get(i).taglen;
                int[] mycounts = tagCountsMatrix.get(i).counts;
                int counttotal = mycounts[0];
                

                outStream.writeBytes(tagSequence + "\t" + tagLen + "\t" + counttotal);


                for (int j=1; j<=chunkFileNames.length; j++)
                {
                    outStream.writeBytes("\t" + mycounts[j]);
                }
                outStream.writeBytes("\n");
            }
            outStream.close();
           
           //write read distribution
           outStream = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFileName + "\\dep_dist.txt"), 2000));
           outStream.writeBytes("Depth");
            for (int j=0; j<chunkFileNames.length; j++)
            {
                outStream.writeBytes("\t" + removePath(chunkFileNames[j]));
            }
            outStream.writeBytes("\n");  
            
            for (int i=1; i<maxDepth; i++)
            {
                outStream.writeBytes("" +  i);
                for (int j=0; j<chunkFileNames.length; j++)
                {
                    outStream.writeBytes("\t" + tagDistPerFile[j][i]);
                }
                outStream.writeBytes("\n");  
            }
            
            outStream.close();
            
            
            //write shared tags
            outStream = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFileName + "\\shared_tags.txt"), 2000));
           outStream.writeBytes("Range\tShared");
            for (int j=0; j<chunkFileNames.length; j++)
            {
                outStream.writeBytes("\t" + removePath(chunkFileNames[j]));
            }
            outStream.writeBytes("\n"); 
            int binsize = 100000;
            int[] bincount = new int[chunkFileNames.length];
            int shared =0;
            Arrays.fill(bincount, 0);
             for (int i=0; i<tagCountsMatrix.size(); i++)
            {
                int[] mycounts = tagCountsMatrix.get(i).counts;
                int t=0;
                for (int j=1; j<=chunkFileNames.length; j++)
                {
                    if (mycounts[j]>0)
                    {
                        bincount[j-1] ++;
                        t++;
                    }
                }
                if (t == chunkFileNames.length)
                {
                    shared ++;
                }
                if   (((i%binsize)==0) && (i>0))
                {
                    outStream.writeBytes("" + i);
                    outStream.writeBytes("\t" + shared);
                    for (int j=0; j<chunkFileNames.length; j++)
                    {
                        outStream.writeBytes("\t" + bincount[j]);
                    }
                    outStream.writeBytes("\n");
                    Arrays.fill(bincount, 0);
                    shared =0;
                }

            }
            outStream.close();          
                        
        }
        catch (Exception e) {
            myLogger.info("Catch in write out  file e=" + e);
            e.printStackTrace();
        }

    }

    private  void outStack(String prefix) {
        System.out.print(prefix+":");
         for (int f = 0; f < rw.length; f++) {
             System.out.print(ctags[f][0]+":");
         }
         myLogger.info("");
    }

    private  void writeTags(long[] writeTag, int minCount) {
        int count=0;
        int[] counts = new int[rw.length + 1];
        Arrays.fill(counts, 0);
        int tagLength=-1;

        //Loop through input buffers, compare current records, write smallest to output buffer
        for (int f = 0; f < rw.length; f++) {
            if(AbstractTags.compareTags(ctags[f],writeTag)==0) {
                count+=ctagCnt[f];
                counts[f+1] = ctagCnt[f];
                tagLength=ctagLength[f];
                for (int j = 0; j < tagLengthInLong; j++) ctags[f][j]=0;
            }
        }
        if(count >= minCount) 
        {
            outCnt++;
            counts[0] = count;
            //counts[rw.length+1] = tagLength;
            TagWithCounts twc = new TagWithCounts(writeTag, tagLength, counts);
            tagCountsMatrix.add(twc);
//            if (outputSeqStr)
//            {
//                tagList.add(writeTag);
//            }
        }
    }


  private  long[] updateCurrentTags() {
        long[] minTag=new long[tagLengthInLong];
        minTag[0]=Long.MAX_VALUE;
        for (int f = 0; f < rw.length; f++) {
            if(ctags[f][0]==0) readNextTag(f);
            if(AbstractTags.compareTags(ctags[f],minTag)<0) minTag=ctags[f].clone();
        }
        return minTag;
    }

    private  void readNextTag(int f) {
        if(rw[f]==null) return;
        try{
            for (int j = 0; j < tagLengthInLong; j++) {
                ctags[f][j] = rw[f].readLong();
            }
            ctagLength[f]=rw[f].readByte();
            int cnt = rw[f].readInt();
            ctagCnt[f] = cnt;
            readCountPerFile[f] += cnt;
            
            int cntbin = cnt;
            
            if (cntbin>maxDepth)
            {
                cntbin = maxDepth;
            }
            
            tagDistPerFile[f][cntbin] ++;
            tagsRead++;
            

        } catch(IOException eof) {
            try{
                myLogger.info("Finished reading file "+f+".");
                rw[f].close();
                rw[f]=null;
                for (int i= 0; i < tagLengthInLong; i++) {

                    ctags[f][i]=Long.MAX_VALUE;
                }
                rwOpen--;
            }
            catch(IOException eof2) {
                myLogger.info("Catch closing" + eof2);
                rw[f]=null;
            }
        }
    }
    
    private String removePath (String filePathName)
    {
        if ( filePathName == null )
        return null;

        int dotPos = filePathName.lastIndexOf( '.' );
        int slashPos = filePathName.lastIndexOf( '\\' );
        if ( slashPos == -1 )
        slashPos = filePathName.lastIndexOf( '/' );

        if ( dotPos > slashPos )
        {
        return filePathName.substring( slashPos > 0 ? slashPos + 1 : 0,
            dotPos );
        }

        return filePathName.substring( slashPos > 0 ? slashPos + 1 : 0 );
    }
            
    @Override
    public String getToolTipText() {
        return null;
    }
    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return null;
    }

}

class TagWithCounts
{
    public int[] counts;
    public int taglen;
    public long[] tagseq;
    public TagWithCounts(long[] seq, int thelen, int[] thecounts)
    {
        tagseq =seq;
        taglen =thelen;
        counts = thecounts;
    }

}
class CompareIntArray implements Comparator <TagWithCounts>{
    public int compare(TagWithCounts v1, TagWithCounts v2) 
    {
        if (v1.counts[0] < v2.counts[0])
        {
            return 1;
        }
        else if (v1.counts[0] == v2.counts[0])
        {
            return 0;
        }
        else
        {
            return -1; 
        }
    }
}