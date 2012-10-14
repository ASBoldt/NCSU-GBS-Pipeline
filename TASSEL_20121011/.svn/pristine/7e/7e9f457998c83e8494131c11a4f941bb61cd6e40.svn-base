package net.maizegenetics.genome.GBS.deNovoContigs;

import net.maizegenetics.genome.GBS.*;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.TreeMap;
import net.maizegenetics.genome.BaseEncoder;

/**
 * User: glaubitz
 * Date: Sept 28, 2010
 * To change this template use File | Settings | File Templates.
 */
public class deNovoContigTags {
    private TagFromDeNovoContig[] theTags;
    private int currentSize=0, maxSize=0;

    public deNovoContigTags(String inFile, boolean binary) {
        if(binary) {
            File fn=new File(inFile);
//            int rows=(int)(fn.length()/54);   // temp fix
            int rows=(int)(fn.length()/TagFromDeNovoContig.totalByteSize);
            initMatrices(rows);
            readBinaryNovelContigTags(fn, rows);
            currentSize=rows;
        } else throw new UnsupportedOperationException("Not supported yet.");
    }
    
    public deNovoContigTags(int rows) {
        initMatrices(rows);
    }

    private void initMatrices(int rows) {
        theTags=new TagFromDeNovoContig[rows];
        maxSize=rows;
    }

    protected void printRows(int numRows) {
        System.out.println(TagFromDeNovoContig.getHeaderString());
        for(int i=0; i<numRows; i++) {
            System.out.println(theTags[i].toString());
        }
    }

   protected void printRows(int numRows, short refCntCutoff, float bestRCutoff) {
       int outCount=0;
       System.out.println(TagFromDeNovoContig.getHeaderString());
       for(int i=0; outCount<numRows; i++) {
            if(theTags[i].getRefCnt()<refCntCutoff || theTags[i].getBestR()>bestRCutoff) continue;
            System.out.println(theTags[i].toString());
            outCount++;
        }
    }

    public enum deNovoCtgTagsSortType { byCONTIG, byGENCHR, byTAG }

    protected long sortTable(deNovoCtgTagsSortType sortType) {
        System.out.print("Sorting tags from de novo contigs...");
        long time=System.currentTimeMillis();
        switch (sortType) {
            case byCONTIG : Arrays.sort(theTags);
                            break;
            case byGENCHR : Arrays.sort(theTags, new PositionComparator());
                            break;
            case byTAG    : Arrays.sort(theTags, new SeqComparator());
                            break;
        }
        long totalTime=System.currentTimeMillis()-time;
        System.out.println("Done in "+totalTime+"ms");
        return totalTime;
    }

    public int addTag(long[] tag, int novelContig, byte tagNum, byte nTagsInCtg, byte inNovelCtg, byte genChr,
            int startPos, int endPos, short refCnt, short altCnt, float bestR, byte chromosome,
            byte strand, int positionMin, int positionMax, short nextCutDistance, byte divergence,
            byte multimaps, byte phyConsenChr, int phyConsenPosition) {
        if(currentSize>=maxSize) return -1;  //failed to add
        TagFromDeNovoContig newTag = new TagFromDeNovoContig(tag, novelContig, tagNum, nTagsInCtg, inNovelCtg, genChr,
            startPos, endPos, refCnt, altCnt, bestR, chromosome, strand, positionMin, positionMax,
            nextCutDistance, divergence, multimaps, phyConsenChr, phyConsenPosition);  
        theTags[currentSize] = newTag;
        currentSize++;
        return currentSize;
    }

    public int addTag(TagFromDeNovoContig tfnc) {
        if(currentSize>=maxSize) return -1;  //failed to add
        theTags[currentSize]=tfnc;
        currentSize++;
        return currentSize;
    }

    protected void readBinaryNovelContigTags(File currentFile, int nTagsToRead) {
        long[] tag=new long[2];
        try{
         DataInputStream dis = new DataInputStream(new BufferedInputStream (new FileInputStream(currentFile),65536));
         System.out.println("File = " + currentFile);
            String temp, sl;
            for(int i=0; i<nTagsToRead; i++) {
                tag[0]=dis.readLong();
                tag[1]=dis.readLong();
                int novelContig = dis.readInt();
                byte tagNum = dis.readByte();
                byte nTagsInCtg = dis.readByte();
                byte inNovelCtg = dis.readByte();    
                byte genChr = dis.readByte();
                int startPos = dis.readInt();
                int endPos = dis.readInt();
                short refCnt = dis.readShort();
                short altCnt = dis.readShort();
                float bestR = dis.readFloat();
                byte chrByte=dis.readByte();
                byte strand=dis.readByte();
                int positionMin=dis.readInt();
                int positionMax=dis.readInt();
                short nextCutDistance=dis.readShort();
                byte divergence=dis.readByte();
                byte multimaps=dis.readByte();
                byte phyConsenChr=dis.readByte();     // Uncomment | comment  (temp fix)
                int phyConsenPosition=dis.readInt();  // Uncomment | comment  (temp fix)
                addTag(tag, novelContig, tagNum, nTagsInCtg, inNovelCtg, genChr, startPos, endPos, refCnt, altCnt, bestR,
                        chrByte, strand, positionMin, positionMax, nextCutDistance, divergence, multimaps, phyConsenChr, phyConsenPosition);
                                                                                      // change temp fix  "(byte) -1, Integer.MIN_VALUE" to "phyConsenChr, phyConsenPosition"
            }
        dis.close();
        }
        catch(Exception e) {
            System.out.println("Error c="+currentSize+" e="+e);
        }
       System.out.println("Count of Tags="+currentSize);
    }

   protected void writeToFile(File outFile, float bestRCutoff, boolean novelContigsOnly, boolean wPhyOrGenPositionOnly, boolean wGenPositionOnly, boolean ref, boolean binary) {
        System.out.println();
        System.out.println("writeToFile");
        System.out.println("-----------");

       int hapsOutput=0;
       try{
       DataOutputStream fw=new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outFile),4000000));       
       if(!binary) {fw.writeBytes(TagFromDeNovoContig.getHeaderString()+"\n");}
       for(int i=0; i<currentSize; i++) {
         if(theTags[i]==null) continue;
         if(novelContigsOnly && theTags[i].inNovelCtg < 1 ) continue;
         if(wPhyOrGenPositionOnly && theTags[i].phyChr < 0 && theTags[i].genChr < 1) continue;
         if(wGenPositionOnly && theTags[i].genChr < 1) continue;
         if(wGenPositionOnly && ref && theTags[i].altCnt > theTags[i].refCnt) continue;
         if(wGenPositionOnly && !ref && theTags[i].refCnt > theTags[i].altCnt) continue;
         if(theTags[i].getBestR()>bestRCutoff) continue;
         if(!binary) {fw.writeBytes(theTags[i].toString()+"\n");}
         else {fw.write(theTags[i].toByte());}
         hapsOutput++;
       }
       fw.flush();
       fw.close();
       System.out.println("Tag positions written to:"+outFile.toString());
       System.out.println("Number of tags in file:"+hapsOutput);
       }
       catch(Exception e) {
             System.out.println("Catch in writing output file e="+e);
       }
   }

   protected void writeToFile(File outFile) {
       writeToFile(outFile, Float.MAX_VALUE, false, false, false, true, true);
   }

   protected void writeToFile(File outFile, boolean binary) {
       writeToFile(outFile, Float.MAX_VALUE, false, false, false, true, binary);
   }

   protected void writeNovelContigsOnlyToFile(File outFile, boolean binary) {
       writeToFile(outFile, Float.MAX_VALUE, true, false, false, true, binary);
   }

   protected void writeTagsWGeneticPositionToFile(File outFile, boolean ref, boolean binary) {
       writeToFile(outFile, Float.MAX_VALUE, false, false, true, ref, binary);
   }

   protected void writeTagsWPhyOrGenPositionToFile(File outFile, boolean binary) {
       writeToFile(outFile, Float.MAX_VALUE, false, true, false, true, binary);
   }

    public int getSize() {
        return theTags.length;
    }

    public void findNovelContigs(int maxDiv, double matchProportion, double consenseThresh) {
        System.out.println();
        System.out.println("findNovelContigs");
        System.out.println("----------------");
        System.out.println("/n/tmaxDiv = " + maxDiv);
        System.out.println("/tmatchProportion = " + matchProportion);
        System.out.println("/tconsenseThresh = " + consenseThresh + "/n");

        // Nb: theTags must be sorted by contig for this to work
        this.sortTable(deNovoCtgTagsSortType.byCONTIG);

        int nNovelCtgs = 0;
        ContigPhyMappingInfo contigInfo = new ContigPhyMappingInfo(theTags[0].contig);
        contigInfo.contigStartIndex = 0;
        for(int i=0; i<currentSize; i++) {
            if (theTags[i].contig == contigInfo.contigNum) {
                contigInfo.recordPhysicalPosition(theTags[i], maxDiv);
            } else {
                contigInfo.contigEndIndex = i-1;
                contigInfo.determineIfContigIsNovel(matchProportion, consenseThresh);
                if (contigInfo.isNovel) {
                    for (int t = contigInfo.contigStartIndex; t <= contigInfo.contigEndIndex; ++t) theTags[t].inNovelCtg = 1;
                    ++nNovelCtgs;
                } else {
                    for (int t = contigInfo.contigStartIndex; t <= contigInfo.contigEndIndex; ++t) theTags[t].inNovelCtg = 0;
                }
                contigInfo = new ContigPhyMappingInfo(theTags[i].contig);
                contigInfo.contigStartIndex = i;
                contigInfo.recordPhysicalPosition(theTags[i], maxDiv);
            }
            if (i%100000 == 0) System.out.println("...finished tag " + i);
        }
        contigInfo.contigEndIndex = theTags.length-1;
        contigInfo.determineIfContigIsNovel(matchProportion, consenseThresh);  // deal with the final contig
        if (contigInfo.isNovel) {
            for (int t = contigInfo.contigStartIndex; t <= contigInfo.contigEndIndex; ++t) theTags[t].inNovelCtg = 1;
            ++nNovelCtgs;
        } else {
            for (int t = contigInfo.contigStartIndex; t <= contigInfo.contigEndIndex; ++t) theTags[t].inNovelCtg = 0;
        }
        System.out.println(nNovelCtgs + " novel contigs found");
    }

    public void determineConsensusGenPositions(boolean ref) {
        System.out.println();
        System.out.println("determineConsensusGenPositions");
        System.out.println("------------------------------");

        // Nb: theTags must be sorted by contig for this to work
        // If the contig does not have consensus, genChr will be set to -1 for all its tags
        // If the contig has consensus, genChr will be set to -1 for all tags within the contig that do not share the consensus position
        // use writeContigsWGeneticPositionsToFile(File outFile, boolean ref, boolean binary) to filter out the non-consensus tags
        this.sortTable(deNovoCtgTagsSortType.byCONTIG);

        int nContigs = 0;
        int nCtgsWGenConsensus = 0;
        ContigGenMappingInfo contigGenInfo = new ContigGenMappingInfo(theTags[0].contig);
        contigGenInfo.contigStartIndex = 0;
        for(int i=0; i<currentSize; i++) {
            if (theTags[i].contig == contigGenInfo.contigNum) {
                contigGenInfo.recordGeneticPosition(theTags[i], ref);
            } else {
                contigGenInfo.contigEndIndex = i-1;
                if (contigGenInfo.checkGenConsensus()) {
                    ++nCtgsWGenConsensus;
                    for (int t = contigGenInfo.contigStartIndex; t <= contigGenInfo.contigEndIndex; ++t) {
                        if (!contigGenInfo.tagMatchesConsensus(theTags[t])) theTags[t].genChr = -1;
                    }
                } else {
                    for (int t = contigGenInfo.contigStartIndex; t <= contigGenInfo.contigEndIndex; ++t) theTags[t].genChr = -1;
                }
                ++nContigs;
                contigGenInfo = new ContigGenMappingInfo(theTags[i].contig);
                contigGenInfo.contigStartIndex = i;
                contigGenInfo.recordGeneticPosition(theTags[i], ref);
            }
            if (i%100000 == 0) System.out.println("...finished tag " + i);
        }
        // deal with the final contig
        contigGenInfo.contigEndIndex = theTags.length-1;
        if (contigGenInfo.checkGenConsensus()) {
            ++nCtgsWGenConsensus;
            for (int t = contigGenInfo.contigStartIndex; t <= contigGenInfo.contigEndIndex; ++t) {
                if (!contigGenInfo.tagMatchesConsensus(theTags[t])) theTags[t].genChr = -1;
            }
        } else {
            for (int t = contigGenInfo.contigStartIndex; t <= contigGenInfo.contigEndIndex; ++t) theTags[t].genChr = -1;
        }
        ++nContigs;
        System.out.println(nCtgsWGenConsensus + " contigs with genetic consensus found out of " + nContigs + " contigs total");
    }

    public void determineConsensusPhyPositions(int maxDiv) {
        // Nb: theTags must be sorted by contig for this to work (tags from each contig must be grouped together)
        // Initially thought that this should be run AFTER determineConsensusGenPositions and writeContigsWGeneticPositionsToFile
        //    - wanted to get the consensus physical positions for the tags in the contig that agree genetically
        //          = physical position of the contig haplotype (merged tags)
        // BUT, I found that many non-novel contigs ended up without a consensus physical position, so I decided
        //       to run this immediately after physicallyMapAllTags() in the main pipeline class

        System.out.println("\ndetermineConsensusPhyPositions for all de novo contigs");
        System.out.println("------------------------------------------------------");

        this.sortTable(deNovoCtgTagsSortType.byCONTIG);

        int nCtgs = 0;
        ContigPhyMappingInfo contigInfo = new ContigPhyMappingInfo(theTags[0].contig);
        contigInfo.contigStartIndex = 0;
        for(int i=0; i<currentSize; i++) {
            if (theTags[i].contig == contigInfo.contigNum) {
                contigInfo.recordPhysicalPosition(theTags[i], maxDiv);
            } else {
                contigInfo.contigEndIndex = i-1;
                contigInfo.determineConsensusPhyPosition();
                if (contigInfo.bestChr > -1) {
                    if (contigInfo.avgPositionOnBestChr > 0) {
                        for (int t = contigInfo.contigStartIndex; t <= contigInfo.contigEndIndex; ++t) {
                            theTags[t].setPhyConsenChr(contigInfo.bestChr);
                            theTags[t].setPhyConsenPosition(contigInfo.avgPositionOnBestChr);
                        }
                    }
                    else {
                        System.out.println("WARNING: consensus physical position of <1 found for contig " + contigInfo.contigNum);  // this helped detect an overflow bug when int used to calc the sum of positions en route to avg
                    }
                    ++nCtgs;
                }
                contigInfo = new ContigPhyMappingInfo(theTags[i].contig);
                contigInfo.contigStartIndex = i;
                contigInfo.recordPhysicalPosition(theTags[i], maxDiv);
            }
            if (i%100000 == 0) System.out.println("...finished tag " + i);
        }
        contigInfo.contigEndIndex = theTags.length-1;
        contigInfo.determineConsensusPhyPosition();  // deal with the final contig
        if (contigInfo.bestChr > -1) {
            for (int t = contigInfo.contigStartIndex; t <= contigInfo.contigEndIndex; ++t) {
                theTags[t].setPhyConsenChr(contigInfo.bestChr);
                theTags[t].setPhyConsenPosition(contigInfo.avgPositionOnBestChr);
            }
            ++nCtgs;
        }
        System.out.println("Consensus physical positions obtained for " + nCtgs + " contigs");
    }

    public void lookupGeneticMappingResult(File GenMapResultFile, int nGenMappedTags, boolean ref)  {
        // the boolean ref indicates if we expect the tag to segregate with the reference genome (B73) or not (e.g., false = expect to seg with Mo17)
        System.out.println();
        System.out.println("lookupGeneticMappingResult");
        System.out.println("--------------------------");

        String[] GMRTags = new String[nGenMappedTags];
        int[] genChr = new int[nGenMappedTags];
        int[] genStart = new int[nGenMappedTags];
        int[] genEnd = new int[nGenMappedTags];
        int[] refCnt = new int[nGenMappedTags];
        int[] altCnt = new int[nGenMappedTags];
        float[] bestR = new float[nGenMappedTags];

        // tailor these to the GenMapResultFile -- THIS FILE MUST BE SORTED BY TAG
        int tagField = 0,
            genChrField = 1,
            genStartField = 2,
            genEndField = 3,
            refCntField = 5,
            altCntField = 6,
            bestRField = 15;

        try {
            BufferedReader br = new BufferedReader(new FileReader(GenMapResultFile), 65536);
            String inputLine = br.readLine();  // skip the header line
            int i = 0;
            while ( (inputLine = br.readLine()) != null) {
                String[] fields = inputLine.split("\t");
                GMRTags[i] = fields[tagField];
                genChr[i] = Integer.parseInt(fields[genChrField]);
                genStart[i] = Integer.parseInt(fields[genStartField]);
                genEnd[i] = Integer.parseInt(fields[genEndField]);
                refCnt[i] = Integer.parseInt(fields[refCntField]);
                altCnt[i] = Integer.parseInt(fields[altCntField]);
                bestR[i] = Float.parseFloat(fields[bestRField]);
                ++i;
            }
            if (i != nGenMappedTags) {
                System.out.println("WARNING: found " + i + " tags in GenMapResultFile but expected " + nGenMappedTags);
            } else {
                System.out.println("Read " + i + " tags from " + GenMapResultFile);
            }
        } catch (Exception e) {
            System.out.println("Error reading GenMapResultFile: " + e);
        }
        int hasGMR = 0;
        for(int i=0; i<currentSize; i++) {
            long[] dNCTag = new long[2];
            dNCTag[0] = theTags[i].tag0;
            dNCTag[1] = theTags[i].tag1;
            int hit=Arrays.binarySearch(GMRTags, BaseEncoder.getSequenceFromLong(dNCTag));
            if (hit > -1 && 
                        ( (ref && (refCnt[hit] > altCnt[hit])  )  ||         // B73 sequence
                          (!ref && (altCnt[hit] > refCnt[hit]) )       )     // Mo17 sequence
                      ) {
                theTags[i].setGenChr(genChr[hit]);
                theTags[i].setStartGenPos(genStart[hit]);
                theTags[i].setEndGenPos(genEnd[hit]);
                theTags[i].setRefCnt(refCnt[hit]);
                theTags[i].setAltCnt(altCnt[hit]);
                theTags[i].setBestR(bestR[hit]);
                ++hasGMR;
            }
            if (i % 100000 == 0) {
                System.out.println("...finished adding genetic mapping results to " + hasGMR + " out of " + i + " de novo contig tags");
            }
        }
        System.out.println("\nFinished! " + hasGMR + " out of " + currentSize + " de novo contig tags were genetically mapped");
    }

    public void mergeIntoHaplotypes(ReadsByTaxa fullRBT, String outFile) {
        System.out.println();
        System.out.println("mergeIntoHaplotypes");
        System.out.println("-------------------");

        // Nb: theTags must be sorted by contig for this to work
        this.sortTable(deNovoCtgTagsSortType.byCONTIG);

        byte[][] countByContigAndTaxon = new byte[fullRBT.haplotypeNum][fullRBT.taxaNum];  // this has the capacity for many more haplotypes than is needed
        for (int tagIndex=0; tagIndex<fullRBT.haplotypeNum; ++tagIndex) {
            for (int taxon=0; taxon<fullRBT.taxaNum; ++taxon) {
                countByContigAndTaxon[tagIndex][taxon] = 0;
            }
        }
        TreeMap <Integer, Integer> contigNumToIndexTreeMap = new TreeMap<Integer, Integer>();  // TreeMap b/c I initially used the first tag (String) to identify the ctg (problematic b/c of duplicate tags)

        // collect tags from the same contig to be merged into a haplotype
        int currentContig = theTags[0].contig;
        int contigCount = 0;
        ArrayList<String> TagsToMerge = new ArrayList<String>();
        for(int i=0; i<currentSize; i++) {
            if (theTags[i].contig == currentContig) {
                TagsToMerge.add(BaseEncoder.getSequenceFromLong(theTags[i].tag0) + BaseEncoder.getSequenceFromLong(theTags[i].tag1));
            }
            else {
                // process the previous contig
                contigNumToIndexTreeMap.put(currentContig, contigCount);
                for (int tagIndex=0; tagIndex < TagsToMerge.size(); ++tagIndex) {
                    for (int taxon=0; taxon<fullRBT.taxaNum; ++taxon) {
                        int newCount = countByContigAndTaxon[contigCount][taxon]
                                + fullRBT.getReadCountForTaxa(fullRBT.getReadIndex(BaseEncoder.getLongArrayFromSeq(TagsToMerge.get(tagIndex))), taxon);
                        if (newCount > Byte.MAX_VALUE) {
                            countByContigAndTaxon[contigCount][taxon] = Byte.MAX_VALUE;
                        } else {
                            countByContigAndTaxon[contigCount][taxon] = (byte) newCount;
                        }
                    }
                }
                // start a new contig
                currentContig = theTags[i].contig;
                ++contigCount;
                TagsToMerge = new ArrayList<String>();
                TagsToMerge.add(BaseEncoder.getSequenceFromLong(theTags[i].tag0) + BaseEncoder.getSequenceFromLong(theTags[i].tag1));
            }
        }
        // process the final contig
        contigNumToIndexTreeMap.put(currentContig, contigCount);
        for (int tagIndex=0; tagIndex < TagsToMerge.size(); ++tagIndex) {
            for (int taxon=0; taxon<fullRBT.taxaNum; ++taxon) {
                int newCount = countByContigAndTaxon[contigCount][taxon]
                        + fullRBT.getReadCountForTaxa(fullRBT.getReadIndex(BaseEncoder.getLongArrayFromSeq(TagsToMerge.get(tagIndex))), taxon);
                if (newCount > Byte.MAX_VALUE) {
                    countByContigAndTaxon[contigCount][taxon] = Byte.MAX_VALUE;
                } else {
                    countByContigAndTaxon[contigCount][taxon] = (byte) newCount;
                }
            }
        }

        // write the new "CHbT" (contig haplotypes by taxa) to a text file
        int ctgsOutput = 0;
        try {
            DataOutputStream fw = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outFile), 65536));
            fw.writeBytes(fullRBT.taxaNum + "\t" + contigNumToIndexTreeMap.size() + "\n");
            for (int t = 0; t < fullRBT.taxaNum; t++) {
                fw.writeBytes("\t" + fullRBT.getTaxaName(t));
            }
            fw.writeBytes("\n");

            for(Integer ctg : contigNumToIndexTreeMap.keySet()) {
                fw.writeBytes("" + ctg);
                for (int t = 0; t < fullRBT.taxaNum; t++) {
                    fw.writeBytes("\t" + countByContigAndTaxon[contigNumToIndexTreeMap.get(ctg)][t]);
                }
                fw.writeBytes("\n");
                ctgsOutput++;
            }
            fw.flush();
            fw.close();
            System.out.println("Counts per contig by taxa written to:" + outFile.toString());
            System.out.println("Number of contigs in file:" + ctgsOutput);
        } catch (Exception e) {
            System.out.println("Catch in writing output file e=" + e);
        }
    }

    private static class ContigPhyMappingInfo {
        int contigNum;
        int contigStartIndex, contigEndIndex;
        boolean isNovel;
        int nTags;
        int nTagsUniqueMatch;
        int[] nUniqueMatchesByChr;  // index 0 to 10; value = count of tags with unique match
        ArrayList<Integer>[] UniqueMatchPositions;
        int bestChr;
        int medianPositionOnBestChr;
        int avgPositionOnBestChr;
        int nTagsWBestPosition;

        ContigPhyMappingInfo(int contigNum) {
            this.contigNum = contigNum;
            contigStartIndex = -1;
            contigEndIndex = -1;
            isNovel = false;
            nTags = 0;
            nTagsUniqueMatch = 0;
            nUniqueMatchesByChr = new int[11];
            for (int chr=0; chr<11; ++chr) nUniqueMatchesByChr[chr] = 0;   // index 0 to 10; value = count of tags with unique match (divergence <= maxDiv)
            UniqueMatchPositions = (ArrayList<Integer>[])new ArrayList[11]; // index = chr 0 to 10; value = ArrayList of positions on that chr
            for (int chr=0; chr<11; ++chr) UniqueMatchPositions[chr] = new ArrayList<Integer>();
            bestChr = -1;
            medianPositionOnBestChr = Integer.MIN_VALUE;
            avgPositionOnBestChr = Integer.MIN_VALUE;
            nTagsWBestPosition = 0;
        }

        void recordPhysicalPosition(TagFromDeNovoContig tag, int maxDiv) {
            ++nTags;
            if (tag.divergence <= maxDiv && tag.multimaps == 1) {
                ++nTagsUniqueMatch;
                ++nUniqueMatchesByChr[tag.phyChr];
                UniqueMatchPositions[tag.phyChr].add(tag.positionMin);
            }
        }
        
        void determineIfContigIsNovel(double matchProportion, double consenseThresh) {
//            if ( nTags < 4 ) { isNovel = false; return; }  # use this if you want to filter out small contigs
            if ( (double) nTagsUniqueMatch/nTags < matchProportion ) { isNovel = true; return; }
            int maxChrN = 0;
            for (int chr=0; chr<11;  ++chr) {
                if (nUniqueMatchesByChr[chr] > maxChrN) {
                    maxChrN = nUniqueMatchesByChr[chr];
                    bestChr = chr;   // NB: in the (very rare?) event of a tie, the lowest numbered chrom wins
                }
            }
            if (bestChr < 0) { isNovel = true; return; }
            Integer[] PositionsOnBestChr = new Integer[UniqueMatchPositions[bestChr].size()];
            UniqueMatchPositions[bestChr].toArray(PositionsOnBestChr);
            Arrays.sort(PositionsOnBestChr);
            medianPositionOnBestChr = PositionsOnBestChr[PositionsOnBestChr.length/2];
            for (Integer pos: PositionsOnBestChr) {
                if (Math.abs(pos-medianPositionOnBestChr) <= 1000000) { ++nTagsWBestPosition; }
            }
            if ( (double) nTagsWBestPosition/nTagsUniqueMatch < consenseThresh ) {
                isNovel = true; return;
            }
            isNovel = false;
            return;
        }

        void determineConsensusPhyPosition() {
            int maxPerfectSinglesByChr = 0;
            for (int chr=0; chr<11; ++chr) {
                if (nUniqueMatchesByChr[chr] > maxPerfectSinglesByChr) {
                    maxPerfectSinglesByChr = nUniqueMatchesByChr[chr];
                    bestChr = chr;   // NB: in the (rare?) event of a tie, the lowest numbered chrom wins
                }
            }
            if (bestChr < 0) { return; }  // chr0 is a legitimate physical position
            Integer phyPositsOnBestChr[] = new Integer[UniqueMatchPositions[bestChr].size()];
            UniqueMatchPositions[bestChr].toArray(phyPositsOnBestChr);
            long sumPositionsOnBestChr = 0;  // previously used an int which caused an overflow bug for ~20 contigs
            for (Integer phyPosit : phyPositsOnBestChr) {
                sumPositionsOnBestChr += phyPosit;
            }
            sumPositionsOnBestChr = sumPositionsOnBestChr/phyPositsOnBestChr.length;
            avgPositionOnBestChr = (int) sumPositionsOnBestChr;
        }
    }

    private static class ContigGenMappingInfo {
        int contigNum;
        int contigStartIndex, contigEndIndex;
        int nTagsGenMapped;
        int[] totRefCntsByChr;  // index 0 to 10; value = sum of refCnt's (or altCnt's if ref=false) per chr
        ArrayList<Integer>[] GenPositions;  // index = chr 0 to 10; int average of genStart and genEnd
        ArrayList<Integer>[] refCnts;  // index = chr 0 to 10; refCnts (or altCnt's if ref=false) for each position to use for weighting GenPositions
        int bestChr;
        double weightedAvgMbPositionOnBestChr;
        int totRefCntsOnBestChr;
        int totRefCntsForTagsMatchingConsensus;

        ContigGenMappingInfo(int contigNum) {
            this.contigNum = contigNum;
            contigStartIndex = -1;
            contigEndIndex = -1;
            nTagsGenMapped = 0;
            totRefCntsByChr = new int[11];
            for (int chr=0; chr<11; ++chr) totRefCntsByChr[chr] = 0;
            GenPositions = (ArrayList<Integer>[])new ArrayList[11]; 
            for (int chr=0; chr<11; ++chr) GenPositions[chr] = new ArrayList<Integer>();
            refCnts = (ArrayList<Integer>[])new ArrayList[11];
            for (int chr=0; chr<11; ++chr) refCnts[chr] = new ArrayList<Integer>();
            bestChr = -1;
            weightedAvgMbPositionOnBestChr = 0.0;
            totRefCntsOnBestChr = 0;
            totRefCntsForTagsMatchingConsensus = 0;
        }

        void recordGeneticPosition(TagFromDeNovoContig tag, boolean ref) {
            if (tag.genChr < 1) return;
            if (ref && tag.altCnt > tag.refCnt) return;
            if (!ref && tag.refCnt > tag.altCnt) return;
            ++nTagsGenMapped;
            if (ref) {
                totRefCntsByChr[tag.genChr] += tag.refCnt;
                refCnts[tag.genChr].add( (int) tag.refCnt);
            } else {
                totRefCntsByChr[tag.genChr] += tag.altCnt;
                refCnts[tag.genChr].add( (int) tag.altCnt);
            }
            GenPositions[tag.genChr].add( (tag.endGenPos+tag.startGenPos)/2 );
        }

        boolean checkGenConsensus() {
            int maxRefCntByChr = 0;
            for (int chr=1; chr<11; ++chr) {
                if (totRefCntsByChr[chr] > maxRefCntByChr) {
                    maxRefCntByChr = totRefCntsByChr[chr];
                    bestChr = chr;   // NB: in the (very rare?) event of a tie, the lowest numbered chrom wins
                }
            }
            if (bestChr < 1) { return false; }
            Integer genPositsOnBestChr[] = new Integer[GenPositions[bestChr].size()];
            GenPositions[bestChr].toArray(genPositsOnBestChr);
            Integer refCntsOnBestChr[] = new Integer[refCnts[bestChr].size()];
            refCnts[bestChr].toArray(refCntsOnBestChr);
            for (int p = 0; p < genPositsOnBestChr.length; ++p) {
                weightedAvgMbPositionOnBestChr += refCntsOnBestChr[p] * ( (double) genPositsOnBestChr[p] / 1000000 );
                totRefCntsOnBestChr += refCntsOnBestChr[p];
            }
            weightedAvgMbPositionOnBestChr /= totRefCntsOnBestChr;

            int nTagsMatchingConsensus = 0;
            for (int p = 0; p < genPositsOnBestChr.length; ++p) {
                if ( Math.abs( (double) genPositsOnBestChr[p]/1000000 - weightedAvgMbPositionOnBestChr) < 50) {
                    ++nTagsMatchingConsensus;
                    totRefCntsForTagsMatchingConsensus += refCntsOnBestChr[p];
                }
            }
            if ( totRefCntsForTagsMatchingConsensus > 29 && (double) nTagsMatchingConsensus/nTagsGenMapped > 0.5) return true;
            return false;
        }
        
        boolean tagMatchesConsensus(TagFromDeNovoContig tag) {
            if (tag.genChr==bestChr && Math.abs((double)((tag.endGenPos+tag.startGenPos)/2)/1000000-weightedAvgMbPositionOnBestChr) < 50) {
                return true;
            }
            return false;
        }

    }
}

class PositionComparator implements Comparator{
    public int compare(Object o1, Object o2){
        TagFromDeNovoContig tag1 = (TagFromDeNovoContig)o1;
        TagFromDeNovoContig tag2 = (TagFromDeNovoContig)o2;
        if(tag1.equals(tag2)) return 0;
        if (tag1.genChr < tag2.genChr) {
                return -1;
        } else if (tag1.genChr > tag2.genChr) {
            return 1;
        }
        if (tag1.startGenPos < tag2.startGenPos) {
            return -1;
        } else if  (tag1.startGenPos > tag2.startGenPos) {
            return 1;
        }
        if (tag1.endGenPos < tag2.endGenPos) {
            return -1;
        } else if  (tag1.endGenPos > tag2.endGenPos) {
            return 1;
        }
        return 0;
    }
}

class SeqComparator implements Comparator {
    public int compare(Object o1, Object o2){
        TagFromDeNovoContig novTag1 = (TagFromDeNovoContig)o1;
        TagFromDeNovoContig novTag2 = (TagFromDeNovoContig)o2;
        if(novTag1.equals(novTag2)) return 0;
        if (novTag1.tag0 < novTag2.tag0) {
            return -1;
        } else if (novTag1.tag0 > novTag2.tag0) {
            return 1;
        }
        if (novTag1.tag1 < novTag2.tag1) {
            return -1;
        } else if (novTag1.tag1 > novTag2.tag1) {
            return 1;
        }
        if (novTag1.contig < novTag2.contig) {
            return -1;
        } else if (novTag1.contig > novTag2.contig) {
            return 1;
        }
        if (novTag1.tagNum < novTag2.tagNum) {
            return -1;
        } else if (novTag1.tagNum > novTag2.tagNum) {
            return 1;
        }
        return 0;
    }
}
