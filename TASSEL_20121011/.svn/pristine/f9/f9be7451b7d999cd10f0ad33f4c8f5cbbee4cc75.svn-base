/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.genome.GBS;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.util.Map.Entry;
import java.util.Random;
import java.util.TreeMap;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AllelePositionBLOBUtils;
import net.maizegenetics.pal.alignment.CombineAlignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.Locus;
import net.maizegenetics.pal.statistics.FisherExact;
import org.apache.commons.math.distribution.BinomialDistributionImpl;

/**
 *
 * @author ed
 */
public class MapHaplotypes {

   
    private MapHaplotypes() {
    }

    static public ReadsWPhysicalMap blastUsingPhysicalMap(ReadsWPhysicalMap physicalMap,
            ReadsWPhysicalMap workingMap, int maxMismatch) {
        physicalMap.sortTable(true);
        ReadBLASTer theReadBLASTer = new ReadBLASTer(physicalMap);
        int perfectSingleHit = 0, perfectMultiHit = 0, newSingleMap = 0, newMultiMap = 0, noMap = 0;
        int count = 0;
        for (int i = 0; i < workingMap.getReadTotal(); i++) {
            ReadWithPositionInfo workReadPos = workingMap.getReadWithPosition(i);
            long[] currRead = workReadPos.getRead();
            count++;
            int[] phyHits = physicalMap.getReadIndexSet(currRead);
            TreeMap<Integer, Integer> al;//=theReadBLASTer.findMatchesWithIntLengthWords(currRead, 3,true);
//            System.out.println(i+"\t"+Arrays.toString(phyHits)+"  "+al.toString());
//            if(count%100000==0) System.out.println("currRead:"+BaseEncoder.getSequenceFromLong(currRead));
            if (phyHits != null) {
                if (phyHits.length == 1) {
                    workReadPos.setPosition(physicalMap.getReadWithPosition(phyHits[0]));
                    workReadPos.setMultimaps(1);
////                    System.out.println("currRead:"+BaseEncoder.getSequenceFromLong(currRead)+" "+"physRead:"+BaseEncoder.getSequenceFromLong(physicalMap.shp[phyHits[0]].readSet));
                    workReadPos.setDivergence(0);
                    workReadPos.nextCutDistance = physicalMap.getReadWithPosition(phyHits[0]).nextCutDistance;
                    perfectSingleHit++;
                } else {
                    workReadPos.setMultimaps(phyHits.length);
                    workReadPos.setDivergence(0);
                    //                   System.out.println("MM  currRead:"+BaseEncoder.getSequenceFromLong(currRead)+" "+"physRead:"+BaseEncoder.getSequenceFromLong(physicalMap.shp[phyHits[1]].readSet));
                    perfectMultiHit++;
                }
            } else {
                al = theReadBLASTer.findMatchesWithIntLengthWords(currRead, maxMismatch, true);
                if (al.size() == 1) {
                    workReadPos.setPosition(physicalMap.getReadWithPosition(al.firstEntry().getKey()));
                    workReadPos.setMultimaps(1);
                    workReadPos.setDivergence(al.firstEntry().getValue());
                    workReadPos.nextCutDistance = physicalMap.getReadWithPosition(al.firstEntry().getKey()).nextCutDistance;
                    newSingleMap++;
                } else if (al.size() > 1) {
                    workReadPos.setMultimaps(al.size());
                    workReadPos.setDivergence(al.firstEntry().getValue());
                    newMultiMap++;
                } else {
                    noMap++;
                }
            }
//            if(count%100000==0){
//               System.out.print("Total:"+workingMap.getReadTotal()+"PerfectSingle:" + perfectSingleHit + "  perfMulti:" + perfectMultiHit);
//                System.out.println("newSingleMap:"+newSingleMap+"\tnewMultiMap:"+newMultiMap+"\tnoMap"+noMap);
//
//            }
        }
        System.out.print("MaxMismatch:" + maxMismatch + "\tTotal:" + workingMap.getReadTotal() + "\tPerfectSingle:" + perfectSingleHit + "\tperfMulti:" + perfectMultiHit);
        System.out.println("\tnewSingleMap:" + newSingleMap + "\tnewMultiMap:" + newMultiMap + "\tnoMap:" + noMap);
        theReadBLASTer = null;
        physicalMap = null;
        return workingMap;
    }

    static public ReadsWPhysicalMap checkWorkingMapWithGenetic(String geneticMapFile, ReadsWPhysicalMap workingMap,
            ReadsByTaxa theRBT, boolean checkRestOfChromosomes, double sigThres) {
        FisherExact theFisherExact = new FisherExact(1000);
        //load genetic map into Treemap with chromosome names and alignment
        TreeMap<String, Alignment> alignmentTree = loadGeneticAlignments(geneticMapFile);
        int[] taxaMapAlignmentToRBT = new int[theRBT.getTaxaCount()];
        Alignment a = alignmentTree.firstEntry().getValue();
        for (int t = 0; t < taxaMapAlignmentToRBT.length; t++) {
            taxaMapAlignmentToRBT[t] = a.getIdGroup().whichIdNumber(theRBT.getTaxaName(t));
        }
//        taxaMapAlignmentToRBT=permuteTaxaRedirectDistToHap(taxaMapAlignmentToRBT);
        int mapsSomewhere = 0, perfectSingleHit = 0, nonFlankBetter = 0, nonFlank100XBetter = 0;
        int count = 0;
        for (int i = 0; i < workingMap.getReadTotal(); i++) {
            ReadWithPositionInfo workReadPos = workingMap.getReadWithPosition(i);
            int rbtIndex = theRBT.getReadIndex(workReadPos.getRead());
            if (rbtIndex < 0) {
                continue; //for some reason the read is no in the read by taxa list
            }
            Alignment test;
            //best chromosome P - record the best of the two flanking markers
            if (workReadPos.isChromosomeKnown() == false) {
                continue;
            }
//  if(workReadPos.isChromosomeKnown()==true) continue;
//  if(workReadPos.getDivergence()!=0) continue;
//  test=alignmentTree.get("1");
            test = alignmentTree.get("" + workReadPos.getChromosome());
            if (test == null) {
                //System.out.println("Chromosome "+workReadPos.getChromosomeString()+" missing");
                continue;
            }
            int searchSite = test.getSiteOfPhysicalPosition(workReadPos.positionMin, null);
            int leftSite, rightSite;
            if (-(searchSite + 1) >= test.getSiteCount() - 1) {
                searchSite = test.getSiteCount() - 1;
            }
            if (searchSite < 0) {
                leftSite = -(searchSite + 1);
                rightSite = leftSite + 1;
            } else {
                leftSite = rightSite = searchSite;
            }
//            double leftP=scoreByFisherExact(test, leftSite, theRBT, rbtIndex, taxaMapAlignmentToRBT, theFisherExact);
//            double rightP=scoreByFisherExact(test, rightSite, theRBT, rbtIndex, taxaMapAlignmentToRBT, theFisherExact);
            double flankP = scoreByDCO(test, leftSite, theRBT, rbtIndex, taxaMapAlignmentToRBT);
//            double flankP=Math.min(leftP, rightP);
            //compare to a scan against of the rest of the genome
            double bestNonChrP = 1;
            String bestNonChr = "  ";
            int bestNonChrPos = -1;
            if (checkRestOfChromosomes) {
                for (Entry<String, Alignment> nonAEntry : alignmentTree.entrySet()) {
                    if (nonAEntry.getKey().toString().equals("" + workReadPos.getChromosome())) {
                        continue;
                    }
                    test = nonAEntry.getValue();
                    for (int s = 0; s < test.getSiteCount(); s++) {
                        // double p=scoreByFisherExact(test, s, theRBT, rbtIndex, taxaMapAlignmentToRBT, theFisherExact);
                        double p = scoreByDCO(test, s, theRBT, rbtIndex, taxaMapAlignmentToRBT);
                        if (p < bestNonChrP) {
                            bestNonChrP = p;
                            bestNonChr = test.getLocusName(s);
                            bestNonChrPos = test.getPositionInLocus(s);
                        }
                    }
                }
            }
            //if the flank is the best set P-value to flank, otherwise set to negative of rest of genome
            if (Double.MAX_VALUE == flankP) {
                continue;  //a problem with the mapping
            }
            if ((flankP < sigThres) || (bestNonChrP < sigThres)) {
                mapsSomewhere++;
            }
            if ((flankP < sigThres)) {
                perfectSingleHit++;
            }
            if (bestNonChrP < flankP) {
//                System.out.println(flankP+"\t"+bestNonChrP+"\t"+bestNonChr+"\t"+bestNonChrPos+"\t"+workReadPos.toString());
                nonFlankBetter++;
                if ((flankP / 100) > bestNonChrP) {
//                    System.out.println(flankP+"\t"+bestNonChrP+"\t"+bestNonChr+"\t"+bestNonChrPos+"\t"+workReadPos.toString());
                    nonFlank100XBetter++;
                }
                flankP = -1 * flankP;
            }


            count++;
            workReadPos.setDcoP((float) flankP);
            //  System.out.println(flankP+"\t"+bestNonChrP+"\t"+bestNonChr+"\t"+bestNonChrPos+"\t"+workReadPos.toString());
            if (count % 1000 == 0) {
                System.out.println("Tests of flanking markers:" + count + " mapping:" + mapsSomewhere
                        + " flankBest:" + perfectSingleHit + " nonFlankBest:" + nonFlankBetter + " nonFlank100XBetter:" + nonFlank100XBetter);
            }
        }
        System.out.println("Total:" + workingMap.getReadTotal() + "PerfectSingle:" + perfectSingleHit);
        System.out.println("Tests of flanking markers:" + count + " mapping:" + mapsSomewhere
                + " flankBest:" + perfectSingleHit + " nonFlankBest:" + nonFlankBetter + " nonFlank100XBetter:" + nonFlank100XBetter);
        return workingMap;
    }

    /**
     * Tests each read (tag) in a ReadsWPhysicalMap for linkage to physically closest framework marker
     *
     * Framework marker alleles must be coded as A and C in HapMap format. Recomb with A recorded in dcoP.
     * Recomb with C recored in mapP.
     *
     * @param frameworkGeneticMapFile Path to external framework genetic map file (HapMap format)
     * @param workingMap              ReadsWPhysicalMap containing the reads that you want to map.  These should have first been physically mapped.
     * @param theRBT                  ReadsByTaxa indicating the distribution each read (tag) across taxa
     * @param sigThres                Minimum recomb rate to bother recording (dcoP or mapP set to 1.0 otherwise)
     * @param minCountForTesting      Minimum sample size to bother calculating recomb (dcoP or mapP set to 1.0 otherwise)
     * @return                        ReadsWPhysicalMap that now contains the recomb rates (dcoP and mapP)
     */
    static public ReadsWPhysicalMap mapGeneticallyVsFramework(String frameworkGeneticMapFile, ReadsWPhysicalMap workingMap,
            ReadsByTaxa theRBT, double sigThres, int minCountForTesting, boolean multipleChrs) {
        //load genetic map into Treemap with chromosome names and alignment
        TreeMap<String, Alignment> alignmentTree = loadGeneticAlignments(frameworkGeneticMapFile, multipleChrs);
        int[] taxaMapAlignmentToRBT = new int[theRBT.getTaxaCount()];
        Alignment a = alignmentTree.firstEntry().getValue();
        for (int t = 0; t < taxaMapAlignmentToRBT.length; t++) {
            taxaMapAlignmentToRBT[t] = a.getIdGroup().whichIdNumber(theRBT.getTaxaName(t));
        }
        int mapsToFlank = 0, mapsToA = 0, mapsToC = 0, doesNotMapToFlank = 0;
        int notInRBT = 0, noPhyPosit = 0, noChr = 0;
        int count = 0;
        for (int i = 0; i < workingMap.getReadTotal(); i++) {
            ReadWithPositionInfo workReadPos = workingMap.getReadWithPosition(i);
            int rbtIndex = theRBT.getReadIndex(workReadPos.getRead());
            if (rbtIndex < 0) { ++notInRBT; continue;} //for some reason the read is not in the read by taxa list
            if (workReadPos.isChromosomeKnown() == false) { ++noPhyPosit; continue; }
            Alignment test = alignmentTree.get("" + workReadPos.getChromosome());
            if (test == null) { ++noChr; continue; }

            int searchSite = test.getSiteOfPhysicalPosition(workReadPos.positionMin, null);

            //  searchSite is the rightSite: "The insertion point is defined as the point at which the key
            //     would be inserted into the array: the index of the first element greater than the key"
            int testSite, leftSite, rightSite;
            if (-(searchSite + 1) >= test.getSiteCount()) {
                searchSite = test.getSiteCount();
            }
            if (searchSite < 0) {
                leftSite = -(searchSite + 1)-1;
                if (leftSite < 0) {leftSite = 0;}
                rightSite = -(searchSite + 1);
                if (rightSite < 0) {rightSite = 0;}
            } else {
                leftSite = searchSite-1;
                rightSite = searchSite;
            }
            if (Math.abs(workReadPos.positionMin-test.getPositionInLocus(rightSite))
                    < Math.abs(workReadPos.positionMin)-test.getPositionInLocus(leftSite)) {
                testSite = rightSite;
            } else {
                testSite = leftSite;
            }
            double flankRecombAlleleA = testFlankingAlleleA(testSite, test, theRBT, taxaMapAlignmentToRBT, rbtIndex, minCountForTesting);
            double flankRecombAlleleC = testFlankingAlleleC(testSite, test, theRBT, taxaMapAlignmentToRBT, rbtIndex, minCountForTesting);
            if (flankRecombAlleleA < sigThres || flankRecombAlleleC < sigThres) {
                mapsToFlank++;
                if (flankRecombAlleleA < sigThres) {
                    workReadPos.setDcoP((float) flankRecombAlleleA);  // use DcoP to record recomb rate with the A allele
                    workReadPos.setMapP((float) 1.0);
                    ++mapsToA;
                }
                if (flankRecombAlleleC < sigThres) {
                    workReadPos.setMapP((float) flankRecombAlleleC);  // use MapP to record recomb rate with the C allele
                    workReadPos.setDcoP((float) 1.0);
                    ++mapsToC;
                }
            } else {
                ++doesNotMapToFlank;
                workReadPos.setDcoP((float) 1.0);
                workReadPos.setMapP((float) 1.0);
            }
            count++;
            if (count % 1000 == 0) {
                System.out.println("Tests of flanking markers:" + count + "  mapsToFlank:" + mapsToFlank + "  mapsToA:" + mapsToA + "  mapsToC:" + mapsToC);
            }
        }
        System.out.println("Total:" + workingMap.getReadTotal() + "  MapsToFlank:" + mapsToFlank + "  doesNotMapToFlank:" + doesNotMapToFlank
                + "  notInRBT:" + notInRBT + "  noPhyPosit:" + noPhyPosit + "  noChr:" + noChr);
        System.out.println("Tests of flanking markers:" + count + "  mapsToFlank:" + mapsToFlank + "  mapsToA:" + mapsToA + "  mapsToC:" + mapsToC);
        return workingMap;
    }

    static public ReadsWPhysicalMap checkWorkingMapWithGeneticBC2(String geneticMapFile, ReadsWPhysicalMap workingMap,
            ReadsByTaxa theRBT, double sigThres) {
        //load genetic map into Treemap with chromosome names and alignment
        TreeMap<String, Alignment> alignmentTree = loadGeneticAlignments(geneticMapFile);
        int[] taxaMapAlignmentToRBT = new int[theRBT.getTaxaCount()];
        Alignment a = alignmentTree.firstEntry().getValue();
        for (int t = 0; t < taxaMapAlignmentToRBT.length; t++) {
            taxaMapAlignmentToRBT[t] = a.getIdGroup().whichIdNumber(theRBT.getTaxaName(t));
        }
        int mapsToFlank = 0, mapsToA = 0, mapsToC = 0, doesNotMapToFlank = 0;
        int notInRBT = 0, noPhyPosit = 0, noChr = 0;
        int count = 0;
        for (int i = 0; i < workingMap.getReadTotal(); i++) {
            ReadWithPositionInfo workReadPos = workingMap.getReadWithPosition(i);
            int rbtIndex = theRBT.getReadIndex(workReadPos.getRead());
            if (rbtIndex < 0) { ++notInRBT; continue;} //for some reason the read is not in the read by taxa list
            if (workReadPos.isChromosomeKnown() == false) { ++noPhyPosit; continue; }
            Alignment test = alignmentTree.get("" + workReadPos.getChromosome());
            if (test == null) { ++noChr; continue; }

            int searchSite = test.getSiteOfPhysicalPosition(workReadPos.positionMin, null);

            //  searchSite is the rightSite: "The insertion point is defined as the point at which the key
            //     would be inserted into the array: the index of the first element greater than the key"
            int testSite, leftSite, rightSite;
            if (-(searchSite + 1) >= test.getSiteCount()) {
                searchSite = test.getSiteCount();
            }
            if (searchSite < 0) {
                leftSite = -(searchSite + 1)-1;
                if (leftSite < 0) {leftSite = 0;}
                rightSite = -(searchSite + 1);
                if (rightSite < 0) {rightSite = 0;}
            } else {
                leftSite = searchSite-1;
                rightSite = searchSite;
            }
            if (Math.abs(workReadPos.positionMin-test.getPositionInLocus(rightSite))
                    < Math.abs(workReadPos.positionMin)-test.getPositionInLocus(leftSite)) {
                testSite = rightSite;
            } else {
                testSite = leftSite;
            }
            double flankRecombMinorAllele = testFlankingAlleleA(testSite, test, theRBT, taxaMapAlignmentToRBT, rbtIndex, 5);  // assumes the minor (teosinte) allele is A = refAllele
            double flankRecombMajorAllele = testFlankingAlleleC(testSite, test, theRBT, taxaMapAlignmentToRBT, rbtIndex, 50);  // assumes the major (W22) allele is C (for "corn")
            if (flankRecombMinorAllele < sigThres || flankRecombMajorAllele < sigThres) {
                mapsToFlank++;
                if (flankRecombMinorAllele < sigThres) {
                    workReadPos.setDcoP((float) flankRecombMinorAllele);
                    workReadPos.setMapP((float) 1.0);
                    ++mapsToA;
                }
                if (flankRecombMajorAllele < sigThres) {
                    workReadPos.setMapP((float) flankRecombMajorAllele);  // use MapP to record recomb rate with W22 (recurrent parent) allele
                    workReadPos.setDcoP((float) 1.0);
                    ++mapsToC;
                }
            }
            else {
                ++doesNotMapToFlank;
                workReadPos.setDcoP((float) 1.0);
                workReadPos.setMapP((float) 1.0);
            }
            count++;
            if (count % 1000 == 0) {
                System.out.println("Tests of flanking markers:" + count + "  mapsToFlank:" + mapsToFlank + "  mapsToTeo:" + mapsToA + "  mapsToW22:" + mapsToC);
            }
        }
        System.out.println("Total:" + workingMap.getReadTotal() + "  MapsToFlank:" + mapsToFlank + "  doesNotMapToFlank:" + doesNotMapToFlank
                + "  notInRBT:" + notInRBT + "  noPhyPosit:" + noPhyPosit + "  noChr:" + noChr);
        System.out.println("Tests of flanking markers:" + count + "  mapsToFlank:" + mapsToFlank + "  mapsToTeo:" + mapsToA + "  mapsToW22:" + mapsToC);
        return workingMap;
    }

    static public ReadsWPhysicalMap mapWorkingMapWithGenetic(Alignment[] gMap, ReadsWPhysicalMap workingMap,
            ReadsByTaxa theRBT, double sigThres, int minCountForTesting, int initSkimRate, String outFile) {
        int[] taxaMapAlignmentToRBT = new int[theRBT.getTaxaCount()];
        for (int t = 0; t < taxaMapAlignmentToRBT.length; t++) {
            taxaMapAlignmentToRBT[t] = gMap[0].getIdGroup().whichIdNumber(theRBT.getTaxaName(t));
        }
 //       taxaMapAlignmentToRBT=permuteTaxaRedirectDistToHap(taxaMapAlignmentToRBT);
        int mapsSomewhere = 0, perfectSingleHit = 0, lessThanMin = 0, greaterThanMin = 0;
        int count = 0;
        try{
        DataOutputStream fw=null;
        if(outFile!=null) fw=new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outFile),4000000));
        for (int i = 0; i < workingMap.getReadTotal(); i++) {
            ReadWithPositionInfo workReadPos = workingMap.getReadWithPosition(i);
            int rbtIndex = theRBT.getReadIndex(workReadPos.getRead());   // this is a binary search, so requires that theRBT is sorted by read (tag)
            if (rbtIndex < 0) {
                continue; //for some reason the read is no in the read by taxa list
            }
            if (theRBT.getTaxaCountForRead(rbtIndex) < minCountForTesting) {  // changed by Jeff: the number of Taxa the read was seen in
                lessThanMin++;
                continue;
            }
            greaterThanMin++;
            int bestSkimAlign=-1;
            double bestSkimP=1;
            for (int currAlign = 0; currAlign < gMap.length; currAlign++) {
                GeneticMapResult gmr=findBestHitOnChromosome(gMap[currAlign], theRBT,
                        taxaMapAlignmentToRBT, rbtIndex, minCountForTesting, initSkimRate);
                if(gmr.bestP<bestSkimP) {bestSkimP=gmr.bestP; bestSkimAlign=currAlign;}
            }  //end of chromosomes
            if(bestSkimAlign<0) continue; //nothing was found, probably sample size
            GeneticMapResult gmr=findBestHitOnChromosome(gMap[bestSkimAlign], theRBT,
                        taxaMapAlignmentToRBT, rbtIndex, minCountForTesting, 1);
            if ((gmr.bestP < sigThres)) {
                mapsSomewhere++;
                workReadPos.setMapP((float) gmr.bestP);
                String blastGeneticAgree = workReadPos.isChromosomeKnown() ? "NO" : "NA";
                if ((gMap[bestSkimAlign].getLocus(gmr.startPos).toString().equals("" + workReadPos.getChromosome())) &&
                        (workReadPos.positionMin >= gMap[bestSkimAlign].getPositionInLocus(gmr.startPos) - 5e6)
                        && (workReadPos.positionMin <= gMap[bestSkimAlign].getPositionInLocus(gmr.endPos) + 5e6)) {
                    perfectSingleHit++;
                    blastGeneticAgree = "YES";
                }
                StringBuilder sb=new StringBuilder();
                sb.append(gMap[bestSkimAlign].getLocus(gmr.startPos) + " "
                        + gMap[bestSkimAlign].getPositionInLocus(gmr.startPos) + " " + gMap[bestSkimAlign].getPositionInLocus(gmr.endPos) + " ");
                String source = (gmr.refCnt > gmr.altCnt) ? "REF" : "ALT";
                sb.append(blastGeneticAgree + " " + gmr.refCnt + " " + gmr.altCnt + " " + source + " ");
                sb.append(workReadPos.toString());
                if(fw!=null) {fw.writeBytes(sb.toString()); fw.writeBytes("\n");}
                else {System.out.println(sb.toString());}
            }

            count++;
            //  System.out.println(flankP+"\t"+bestNonChrP+"\t"+bestNonChr+"\t"+bestNonChrPos+"\t"+workReadPos.toString());
            if (count % 1000 == 0) {
                System.out.println("Total reads:" + i + "   Tested reads:" + count + "   mapping:" + mapsSomewhere
                        + "   BlastWithinGenetic:" + perfectSingleHit);
                if(fw!=null) fw.flush();
            }
        } //end of reads
            if(fw!=null) fw.close();
        }  //end of try
        catch(Exception e) {
            System.out.println("Error writing to file:"+e.getMessage());
            e.printStackTrace();
        }
        System.out.println("Total:" + workingMap.getReadTotal() + "PerfectSingle:" + perfectSingleHit);
        System.out.println("Tests of flanking markers:" + count + " mapping:" + mapsSomewhere
                + " flankBest:" + perfectSingleHit);
        return workingMap;
    }


    private static GeneticMapResult findBestHitOnChromosome(Alignment gMap, ReadsByTaxa theRBT,
            int[] taxaMapAlignmentToRBT, int rbtIndex, int minCountForTesting, int stepSize) {

        int startBest = -1, endBest = -1, altCntBest = -1, refCntBest = -1;
        double bestP = 1;
        for (int geneticSite = 0; geneticSite < gMap.getSiteCount(); geneticSite+=stepSize) {
            int refCnt = 0, altCnt = 0;
            for (int t = 0; t < taxaMapAlignmentToRBT.length; t++) {
                if (taxaMapAlignmentToRBT[t] < 0) {
                    continue;
                }
                if (theRBT.getReadCountForTaxa(rbtIndex, t) < 1) {
                    continue;
                }
                byte currBase = gMap.getBase(taxaMapAlignmentToRBT[t], geneticSite);
                if (currBase == AlignmentFilterByGBSUtils.refAllele) {
                    refCnt++;
                } else if (currBase == AlignmentFilterByGBSUtils.altAllele) {
                    altCnt++;
                }
            }
            int cnt = refCnt + altCnt;
            if (cnt < minCountForTesting) {
                continue;
            }
            double currP = (double) Math.min(refCnt, altCnt) / (double) cnt;
            if (currP < bestP) {
                bestP = currP;
                startBest = endBest = geneticSite;
                altCntBest = altCnt;
                refCntBest = refCnt;
            } else if (currP == bestP) {
                endBest = geneticSite;
            }
        }  //end of one chromosome
        GeneticMapResult gmr=new GeneticMapResult(startBest, endBest, refCntBest,
                altCntBest, bestP, gMap.getLocus(startBest));
        return gmr;
    }

    private static double testFlankingAlleleA(int FlankSite, Alignment gMap, ReadsByTaxa theRBT,
            int[] taxaMapAlignmentToRBT, int rbtIndex, int minCountForTesting) {

        double recomb = 1;
        int refCnt = 0, altCnt = 0;
        for (int t = 0; t < taxaMapAlignmentToRBT.length; t++) {
            if (taxaMapAlignmentToRBT[t] < 0) {
                continue;
            }
            if (theRBT.getReadCountForTaxa(rbtIndex, t) < 1) {
                continue;
            }
            byte currBase = gMap.getBase(taxaMapAlignmentToRBT[t], FlankSite);
            if (currBase == AlignmentFilterByGBSUtils.refAllele) {
                refCnt++;
            } else if (currBase == AlignmentFilterByGBSUtils.altAllele) {
                altCnt++;
            }
        }
        int cnt = refCnt + altCnt;
        if (cnt < minCountForTesting) {
            recomb = 1;
        } else {
            recomb = altCnt / (double) cnt;
        }
        return recomb;
    }

    private static float testFlankingAlleleC(int FlankSite, Alignment gMap, ReadsByTaxa theRBT,
            int[] taxaMapAlignmentToRBT, int rbtIndex, int minCountForTesting) {

        float recomb = 1;
        int refCnt = 0, altCnt = 0;
        for (int t = 0; t < taxaMapAlignmentToRBT.length; t++) {
            if (taxaMapAlignmentToRBT[t] < 0) {
                continue;
            }
            if (theRBT.getReadCountForTaxa(rbtIndex, t) < 1) {
                continue;
            }
            byte currBase = gMap.getBase(taxaMapAlignmentToRBT[t], FlankSite);
            if (currBase == AlignmentFilterByGBSUtils.refAllele) {
                refCnt++;
            } else if (currBase == AlignmentFilterByGBSUtils.altAllele) {
                altCnt++;
            }
        }
        int cnt = refCnt + altCnt;
        if (cnt < minCountForTesting) {
            recomb = 1;
        } else {
            recomb = (float) refCnt / (float) cnt;
        }
        return recomb;
    }


    private static TreeMap<String, Alignment> loadGeneticAlignments(String geneticMapFile) {
        TreeMap<String, Alignment> alignmentMap = new TreeMap<String, Alignment>();
        Alignment[] refMarkers = ((CombineAlignment) ImportUtils.readFromHapmap(geneticMapFile)).getAlignments(); // temp fix, need to be changed for actual functionality
        System.out.println("Genetic Map File Read:" + geneticMapFile);
        for (Alignment a : refMarkers) {
            System.out.println("Locus:" + a.getLocus(0) + " Sites:" + a.getSiteCount() + " Taxa:" + a.getSequenceCount());
            alignmentMap.put(a.getLocus(0).getName(), a);
        }
        return alignmentMap;
    }

    private static TreeMap<String, Alignment> loadGeneticAlignments(String geneticMapFile, boolean multipleChrs) {
        TreeMap<String, Alignment> alignmentMap = new TreeMap<String, Alignment>();
        Alignment[] refMarkers;
        if (multipleChrs) {
            refMarkers = ((CombineAlignment) ImportUtils.readFromHapmap(geneticMapFile)).getAlignments(); // temp fix, need to be changed for actual functionality
        } else {
            refMarkers = new Alignment[1];
            refMarkers[0] = ImportUtils.readFromHapmap(geneticMapFile);
        }
        System.out.println("Genetic Map File Read:" + geneticMapFile);
        for (Alignment a : refMarkers) {
            System.out.println("Locus:" + a.getLocus(0) + " Sites:" + a.getSiteCount() + " Taxa:" + a.getSequenceCount());
            alignmentMap.put(a.getLocus(0).getName(), a);
        }
        return alignmentMap;
    }

    /**
     * Tests whether the presense versus absense of a haplotypes is associated with
     * the score of another marker.
     * @param a Reference alignment
     * @param testSiteInAlignment site within the alignment
     * @param testHaplotype the index of the test haplotype
     * @return p-value of the Fisher exact test
     */
    private static double scoreByFisherExact(Alignment a, int testSiteInAlignment, ReadsByTaxa theRBT,
            int testHaplotype, int[] taxaMapAlignmentToRBT, FisherExact theFisherExact) {
        //This test whether the presenese
        byte refBase = a.getMajorAllele(testSiteInAlignment);
        int[][] states = new int[2][2];
        for (int t = 0; t < taxaMapAlignmentToRBT.length; t++) {
            if (taxaMapAlignmentToRBT[t] < 0) {
                continue;  //missing taxa
            }
            byte b = a.getBase(taxaMapAlignmentToRBT[t], testSiteInAlignment);
            if (isHomozygousGood((char) b) == false) {
                continue;
            }
            if (refBase == b) {
                if (theRBT.getReadCountForTaxa(testHaplotype, t) > 0) {
                    states[0][1]++;
                } else {
                    states[0][0]++;
                }
            } else if (theRBT.getReadCountForTaxa(testHaplotype, t) > 0) {
                states[1][1]++;
            } else {
                states[1][0]++;
            }
        }
        return theFisherExact.getTwoTailedP(states[0][0], states[1][0], states[0][1], states[1][1]);
    }

    /**
     * Tests whether the presense of a haplotypes is associated with
     * the phase of a marker interval.
     * @param a Reference alignment
     * @param testSiteInAlignment site within the alignment
     * @param testHaplotype the index of the test haplotype
     * @return p-value of the Fisher exact test
     */
    private static double scoreByDCO(Alignment a, int testSiteInAlignment, ReadsByTaxa theRBT,
            int testHaplotype, int[] taxaMapAlignmentToRBT) {
        //score the interval with test site on the left
        if (testSiteInAlignment < 0 || testSiteInAlignment >= (a.getSiteCount() - 1)) {
            return Double.MAX_VALUE;  //end of the chromosome so no interval
        }
        //Deals with segregation distortion, no sure a good idea
//        if ((a.getMinorAlleleFrequency(testSiteInAlignment) < .3) || (a.getMinorAlleleFrequency(testSiteInAlignment + 1) < .3)) {
//            return Double.MAX_VALUE;
//        }
        byte refBase1 = a.getMajorAllele(testSiteInAlignment);
        byte refBase2 = a.getMajorAllele(testSiteInAlignment + 1);
        int[][] states = new int[2][2];
        int[][] statesAll = new int[2][2];
        double p = 1.0;
        double dcoProp = 1.0;
        //Create 2X2 matrix for all the data and just for those with haplotypes
        for (int t = 0; t < taxaMapAlignmentToRBT.length; t++) {
            if (taxaMapAlignmentToRBT[t] < 1) {
                continue;  //missing taxa
            }     
            byte b1 = a.getBase(taxaMapAlignmentToRBT[t], testSiteInAlignment);
            byte b2 = a.getBase(taxaMapAlignmentToRBT[t], testSiteInAlignment + 1);
            if ((isHomozygousGood((char) b1) == false) || (isHomozygousGood((char) b2) == false)) {
                continue;  //unresolved for hets
            }
            if (refBase1 == b1) {
                if (refBase2 == b2) {
                    statesAll[0][0]++;
                } else {
                    statesAll[0][1]++;
                }
            } else if (refBase2 == b2) {
                statesAll[1][0]++;
            } else {
                statesAll[1][1]++;
            }
            if (theRBT.getReadCountForTaxa(testHaplotype, t) < 1) {
                continue; //only focus on present haplotypes (a dominant test)
            }
            if (refBase1 == b1) {
                if (refBase2 == b2) {
                    states[0][0]++;
                } else {
                    states[0][1]++;
                }
            } else if (refBase2 == b2) {
                states[1][0]++;
            } else {
                states[1][1]++;
            }
        }
        //Determine phase direction for a biparental cross
        double sumd1 = statesAll[0][0] + statesAll[1][1];
        double sumd2 = statesAll[1][0] + statesAll[0][1];
        int majorAll, majorDom, minorAll, minorDom;
        if (sumd1 > sumd2) {
            if (states[0][0] > states[1][1]) {
                majorAll = statesAll[0][0];
                minorAll = statesAll[1][1];
                majorDom = states[0][0];
                minorDom = states[1][1];
            } else {
                majorAll = statesAll[1][1];
                minorAll = statesAll[0][0];
                majorDom = states[1][1];
                minorDom = states[0][0];
            }
        } else {
            if (states[1][0] > states[0][1]) {
                majorAll = statesAll[1][0];
                minorAll = statesAll[0][1];
                majorDom = states[1][0];
                minorDom = states[0][1];
            } else {
                majorAll = statesAll[0][1];
                minorAll = statesAll[1][0];
                majorDom = states[0][1];
                minorDom = states[1][0];
            }
        }
        int sumDom = minorDom + majorDom;
        dcoProp = (double) minorDom / (double) sumDom;
        double minorProb = (double) minorAll / (double) (minorAll + majorAll);
        BinomialDistributionImpl x = new BinomialDistributionImpl(sumDom, minorProb);
        try {
            p = x.cumulativeProbability(minorDom);
        } catch (Exception e) {
            System.err.println("Error in the BinomialDistributionImpl");
        }

//       if(p<0.000001) System.out.println("scoreByDCO="+Arrays.deepToString(states)+
//                Arrays.deepToString(statesAll)+" dcoProp="+dcoProp+
//                " minorProb="+minorProb+" p="+p+" sites="+testSiteInAlignment
//                +"-"+(testSiteInAlignment+1));

        return p;
    }

    private static boolean isHomozygousGood(char base) {
        byte hf = AllelePositionBLOBUtils.getHalfByteFromBase(base);
        if (hf < 0x4) {
            return true;
        } else {
            return false;
        }
    }

    private static int[] permuteTaxaRedirectDistToHap(int[] taxaRedirectDistToHap) {
        Random x = new Random();
        for (int t = 0; t < taxaRedirectDistToHap.length; t++) {
            if (taxaRedirectDistToHap[t] > -1) {
                int pull = x.nextInt(taxaRedirectDistToHap.length);
                if (taxaRedirectDistToHap[pull] > -1) {
                    int temp = taxaRedirectDistToHap[t];
                    taxaRedirectDistToHap[t] = taxaRedirectDistToHap[pull];
                    taxaRedirectDistToHap[pull] = temp;
                }
            }
        }
        return taxaRedirectDistToHap;
    }

    private static class GeneticMapResult {
        int startPos, endPos, altCnt, refCnt;
        double bestP;
        Locus theChromosome;

        GeneticMapResult(int startPos, int endPos, int refCnt, int altCnt, double bestP, Locus theChromosome) {
            this.startPos=startPos;
            this.endPos=endPos;
            this.refCnt=refCnt;
            this.altCnt=altCnt;
            this.bestP=bestP;
            this.theChromosome=theChromosome;
        }
    }
}

