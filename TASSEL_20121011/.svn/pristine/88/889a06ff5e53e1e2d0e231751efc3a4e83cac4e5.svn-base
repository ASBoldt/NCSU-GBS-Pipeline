package net.maizegenetics.genome.GBS;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Map.Entry;
import java.util.TreeMap;

/**
 * Holds haplotypes data compressed in long and a counts.  Has basic filtering methods.
 *
 * User: ed
 * Date: Jan 26, 2008
 * Time: 8:12:12 AM
 * To change this template use File | Settings | File Templates.
 */
public class ReadsWPhysicalMap implements Reads {
    private ReadWithPositionInfo[] shp;
    private int currentSize=0, maxSize=0;


    public ReadsWPhysicalMap() {

    }

    public ReadsWPhysicalMap(String inFile, boolean binary) {
        if(binary) {
            File fn=new File(inFile);
            int rows=(int)(fn.length()/ReadWithPositionInfo.totalByteSize);
            initMatrices(rows);
            readBinarySolexaHaplotypeMap(new File(inFile),rows);
       //     readBinaryReadsWPositionInfo(new File(inFile),rows);
            currentSize=rows;
        } else throw new UnsupportedOperationException("Not supported yet.");
    }

    public ReadsWPhysicalMap(String inFile, boolean binary, int maxCapacity) {
        if(binary) {
            File fn=new File(inFile);
            int rows=(int)(fn.length()/ReadWithPositionInfo.totalByteSize);
            initMatrices(maxCapacity);
            readBinarySolexaHaplotypeMap(new File(inFile),rows);
       //     readBinaryReadsWPositionInfo(new File(inFile),rows);
            currentSize=rows;
        } else throw new UnsupportedOperationException("Not supported yet.");
    }
    
    public ReadsWPhysicalMap(int rows) {
        initMatrices(rows);
    }

    public ReadsWPhysicalMap(Reads readList) {
        initMatrices(readList.getReadTotal());
        for (int i = 0; i < readList.getReadTotal(); i++) {
            addHaplotype(readList.getRead(i));
        }
    }

    private void initMatrices(int rows) {
        shp=new ReadWithPositionInfo[rows];
        maxSize=rows;
    }

    protected void printRows(int numRows) {
        for(int i=0; i<numRows; i++) {
            System.out.println(shp[i].toString());
        }
    }

   protected void printRows(int numRows, boolean requirePhysPosition) {
       int outCount=0;
       for(int i=0; outCount<numRows; i++) {
            if((requirePhysPosition==true)&&(shp[i].isChromosomeKnown()!=true)) continue;
            System.out.println(shp[i].toString());
            outCount++;
        }
    }

    public long sortTable(boolean byHaplotype) {
        System.out.print("Starting Read Table Sort ...");
        long time=System.currentTimeMillis();
        if(byHaplotype) {
            Arrays.sort(shp);
        } else {
            Arrays.sort(shp, new PositionComparator());
        }
        long totalTime=System.currentTimeMillis()-time;
        System.out.println("Done in "+totalTime+"ms");
        return totalTime;
    }

    public int addHaplotype(long[] haplotype, byte chromosome, byte strand, int positionMin,
            int positionMax, short nextCutDistance, byte divergence) {
        if(currentSize>=maxSize) return -1;  //failed to add
        ReadWithPositionInfo theSHP=new ReadWithPositionInfo(haplotype, chromosome, strand, positionMin,
             positionMax, nextCutDistance, divergence);
        shp[currentSize]=theSHP;
        currentSize++;
        return currentSize;
    }

    public int addHaplotype(long[] haplotype, byte chromosome, byte strand, int positionMin,
            int positionMax, short nextCutDistance, byte divergence, float dcoP, float mapP,
            byte multimaps) {
        if(currentSize>=maxSize) return -1;  //failed to add
        ReadWithPositionInfo theSHP=new ReadWithPositionInfo(haplotype, chromosome, strand, positionMin,
             positionMax, nextCutDistance, divergence, dcoP, mapP, multimaps);
        shp[currentSize]=theSHP;
        currentSize++;
        return currentSize;
    }

    public int addHaplotype(long[] haplotype) {
        if(currentSize>=maxSize) return -1;  //failed to add
        ReadWithPositionInfo theSHP=new ReadWithPositionInfo(haplotype);
        shp[currentSize]=theSHP;
        currentSize++;
        return currentSize;
    }

    public int addHaplotype(ByteBuffer bb) {
        if(currentSize>=maxSize) return -1;  //failed to add
        ReadWithPositionInfo theSHP=new ReadWithPositionInfo(bb);
        shp[currentSize]=theSHP;
        currentSize++;
        return currentSize;
    }

    public int addHaplotype(ReadWithPositionInfo rwpi) {
        if(currentSize>=maxSize) return -1;  //failed to add
        shp[currentSize]=rwpi;
        currentSize++;
        return currentSize;
    }

    protected void readBinarySolexaHaplotypeMap(File currentFile, int hapsToRead) {
    //    String chromosome;
        long[] haplotype=new long[2];
        byte chrB;
        try{
         DataInputStream dis = new DataInputStream(new BufferedInputStream (new FileInputStream(currentFile),65536));
         System.out.println("File = " + currentFile);
            String temp, sl;
            //lonhaplotype=new long[2];
            for(int i=0; i<hapsToRead; i++) {
                haplotype[0]=dis.readLong();
                haplotype[1]=dis.readLong();
                chrB=dis.readByte();
              //  chrB[1]=dis.readByte();
                byte strand=dis.readByte();
                int positionMin=dis.readInt();
                int positionMax=dis.readInt();
                short nextCutDistance=dis.readShort();
                byte divergence=dis.readByte();
                float dcoP=dis.readFloat();
                float mapP=dis.readFloat();
                byte multimaps=dis.readByte();
                addHaplotype(haplotype, chrB, strand, positionMin,
                    positionMax, nextCutDistance, divergence, dcoP, mapP, multimaps);
            }
        dis.close();
        }
        catch(Exception e) {
            System.out.println("Error c="+currentSize+" e="+e);
        }
       System.out.println("Count of Tags="+currentSize);
    }

    protected void readBinaryReadsWPositionInfo(File currentFile, int hapsToRead) {
    //    String chromosome;
        byte[] ba=new byte[ReadWithPositionInfo.totalByteSize];
        try{
         FileInputStream fis = new FileInputStream(currentFile);
         FileChannel fc=fis.getChannel();
     //    ByteBuffer bb=ByteBuffer.allocate(ReadWithPositionInfo.totalByteSize);
         ByteBuffer bb=ByteBuffer.wrap(ba);
         System.out.println("File = " + currentFile);
       //  byte[] b;
            for(int i=0; i<hapsToRead; i++) {
                int rl=fc.read(bb);
                System.out.println("bb:"+bb.toString()+" rl:"+rl);
                if(rl!=ReadWithPositionInfo.totalByteSize) System.out.println("Error reading at row:"+i);
                addHaplotype(bb);
            }
        fis.close();
        }
        catch(Exception e) {
            System.out.println("Error c="+currentSize+" e="+e);
        }
       System.out.println("Count of Tags="+currentSize);
    }

        @Override
    public int getReadIndex(long[] read) {
        ReadWithPositionInfo querySHP=new ReadWithPositionInfo(read);
        querySHP.chromosome=Byte.MIN_VALUE;  //this sets search so that it will
          //be before the first matching read
        int hit1=Arrays.binarySearch(shp, querySHP);
        int firstHit=(-hit1)-1;
        if((firstHit<shp.length)&&Arrays.equals(shp[firstHit].getRead(),read)) return firstHit;
//        int misMatch=BaseEncoder.seqDifferences(read[0], shp[firstHit].readSet[0], 3)+
//           BaseEncoder.seqDifferences(read[1], shp[firstHit].readSet[1], 3);
//        if(misMatch==0) return firstHit;
        return hit1;
    }

    public int[] getReadIndexSet(long[] query) {
        int hit=getReadIndex(query);
        if(hit<0) return null;
        ArrayList<Integer> hits=new ArrayList<Integer>();
        for(int i=hit; (i<shp.length)&&(Arrays.equals(shp[i].getRead(),query)); i++) {
           hits.add(i);
        }
        int[] rhits=new int[hits.size()];
        for (int i = 0; i < rhits.length; i++) rhits[i]=hits.get(i);
        return rhits;
    }

    public int getReadIndex(byte chr, byte strand, int posMin) {
        ReadWithPositionInfo querySHP=new ReadWithPositionInfo(chr, strand, posMin);
        PositionComparator posCompare = new PositionComparator();
        int hit1=Arrays.binarySearch(shp, querySHP, posCompare);
        return hit1;
    }

   protected void writeCountFile(File outFile, int minResolution, boolean requirePhysPosition,
           boolean requireDCOMap, float minDCOP, boolean binary) {
       int hapsOutput=0;
       try{
       DataOutputStream fw=new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outFile),4000000));       
       for(int i=0; i<currentSize; i++) {
         if(shp[i]==null) continue;
         if((requirePhysPosition==true)&&(shp[i].isChromosomeKnown()!=true)) continue;
         if((shp[i].positionMax-shp[i].positionMin)>minResolution) continue;
         if(requireDCOMap) {
             if(Float.isNaN(shp[i].getDcoP())||(shp[i].getDcoP()>minDCOP)||(shp[i].getDcoP()<0)) continue;

         }
         if(!binary) {fw.writeBytes(shp[i].toString()+"\n");}
         else {fw.write(shp[i].toByte());}
         hapsOutput++;
       }
       fw.flush();
       fw.close();
       System.out.println("Haplotypes positions written to:"+outFile.toString());
       System.out.println("Number of Haplotypes in file:"+hapsOutput);
       }
       catch(Exception e) {
             System.out.println("Catch in writing output file e="+e);
       }
   }

   protected void writeCountFile(File outFile) {
       writeCountFile(outFile, Integer.MAX_VALUE, false, false, Float.NaN, true);
   }

   protected void writeCountFile(File outFile, boolean binary) {
       writeCountFile(outFile, Integer.MAX_VALUE, false, false, Float.NaN, binary);
   }

    public int getSize() {
        return shp.length;
    }

   public void getChromosomeDistribution(int binSize) {
       TreeMap<String,Integer> dist=new TreeMap<String,Integer>();
       for(ReadWithPositionInfo theSHP: shp) {
           String bin;
           if(binSize<1) {bin="C\t"+theSHP.getChromosome();}
           else {bin="C\t"+theSHP.getChromosome()+"\t"+Math.floor(theSHP.positionMin/binSize);}
           if((bin.equals("C\t0_\t149.0"))) System.out.println("Error:"+theSHP.toString());
           Integer count = dist.get(bin);
           if(count == null){
                    count = 0;
            }
           dist.put(bin, count + 1);
       }
       for(Entry<String,Integer> e: dist.entrySet()) {
           System.out.println(e.getKey()+"\t"+e.getValue());
       }
       System.out.println("TreeMap"+dist.toString());
   }

    public void getLengthDistribution(int maxSize) {
       int[] sizes=new int[maxSize+1];
       for(ReadWithPositionInfo theSHP: shp) {
           if(theSHP.nextCutDistance<0) {sizes[maxSize-1]++;}
           else if(theSHP.nextCutDistance>maxSize) {sizes[maxSize]++;}
           else {sizes[theSHP.nextCutDistance]++;}
       }
       System.out.println("Size\tCount");
       for(int i=0; i<sizes.length; i++) {
           System.out.println(i+"\t"+sizes[i]);
       }
   }


   public static void main(String[] args) {
       if((args.length!=1)||(args[0].equals("-h"))) {
           System.out.println("Require args: infile ");
           System.exit(1);
       }
        else if (args.length==1) {
           ReadsWPhysicalMap be=new ReadsWPhysicalMap(args[0],true);
           be.printRows(100);
           System.out.print("Sort beginning: ");
           long t=be.sortTable(true);
           System.out.println("Done in "+t+" milliseconds");
           be.printRows(100);
//           be.writeCountFile(new File(args[0]), 10000000, true);
//           System.out.print("Sort beginning: ");
//           t=be.sortTable(false);
//           System.out.println("Done in "+t+" milliseconds");
//           be.printRows(100);
        //   be.getChromosomeDistribution(1000000);
        //   be.getLengthDistribution(5000);
       }
  }

    @Override
    public long[] getRead(int index) {
        if((index<0)||(index>=shp.length)) return null;
        return shp[index].getRead();
    }

    @Override
    public int getReadCount(int index) {
        return Integer.MIN_VALUE;
    }

    @Override
    public boolean areReadsUnique() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getReadTotal() {
        return this.currentSize;
    }

    public ReadWithPositionInfo getReadWithPosition(int index) {
        if((index<0)||(index>=shp.length)) return null;
        return shp[index];
    }
}

class PositionComparator implements Comparator{
    public int compare(Object o1, Object o2){
        ReadWithPositionInfo shp1 = (ReadWithPositionInfo)o1;
        ReadWithPositionInfo shp2 = (ReadWithPositionInfo)o2;
        if(shp1.equals(shp2)) return 0;
        if (shp1.chromosome < shp2.chromosome) {
                return -1;
        } else if (shp1.chromosome > shp2.chromosome) {
            return 1;
        }
        if (shp1.positionMin < shp2.positionMin) {
            return -1;
        } else if  (shp1.positionMin > shp2.positionMin) {
            return 1;
        }
        if (shp1.strand < shp2.strand) {
            return -1;
        } else if  (shp1.strand > shp2.strand) {
            return 1;
        }
        return 0;
	}

}

class HaplotypeComparator implements Comparator{
    public int compare(Object o1, Object o2){
        ReadWithPositionInfo shp1 = (ReadWithPositionInfo)o1;
        ReadWithPositionInfo shp2 = (ReadWithPositionInfo)o2;
        if(shp1.equals(shp2)) return 0;
        if (shp1.read0 < shp2.read0) {
                return -1;
        } else if (shp1.read0 > shp2.read0) {
            return 1;
        }
        if (shp1.read1 < shp2.read1) {
            return -1;
        } else if (shp1.read1 > shp2.read1) {
            return 1;
        }
        return 0;
	}

}

class DCOPComparator implements Comparator{
    public int compare(Object o1, Object o2){
        ReadWithPositionInfo shp1 = (ReadWithPositionInfo)o1;
        ReadWithPositionInfo shp2 = (ReadWithPositionInfo)o2;
        if (Float.isNaN(shp1.dcoP) && Float.isNaN(shp2.dcoP)) {
            return 0;
        }
        else if (Float.isNaN(shp2.dcoP)) {
            return -1;
        }
        else if (Float.isNaN(shp1.dcoP)) {
            return 1;
        }

        if (shp1.dcoP < shp2.dcoP) {
                return -1;
        } else if (shp1.dcoP > shp2.dcoP) {
            return 1;
        }
        return 0;
    }

}

class DivergenceComparator implements Comparator{
    public int compare(Object o1, Object o2){
        ReadWithPositionInfo shp1 = (ReadWithPositionInfo)o1;
        ReadWithPositionInfo shp2 = (ReadWithPositionInfo)o2;
        if (       (shp1.divergence == Byte.MIN_VALUE || shp1.divergence < 0)
                && (shp2.divergence == Byte.MIN_VALUE || shp2.divergence < 0) ) {
            return 0;
        }
        else if (shp2.divergence == Byte.MIN_VALUE || shp2.divergence < 0) {
            return -1;
        }
        else if (shp1.divergence == Byte.MIN_VALUE || shp1.divergence < 0) {
            return 1;
        }

        if (shp1.divergence < shp2.divergence) {
                return -1;
        } else if (shp1.divergence > shp2.divergence) {
            return 1;
        }
        return 0;
    }
}
