/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS.deNovoContigs;
import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import net.maizegenetics.genome.GBS.*;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.Locus;


/**
 *
 * @author jcg233
 */
public class mapContigHaplotypes {

    private mapContigHaplotypes() {
    }

    static public void geneticallyMapContigHaplos(Alignment[] gMap, ContigHaplotypesByTaxa theCHBT, double recombThres, int minCountForTesting, int initSkimRate, String outFile) {
        System.out.println();
        System.out.println("geneticallyMapContigHaplos");
        System.out.println("--------------------------");
        System.out.println("Parameters:");
        System.out.println("\trecombThresh: " + recombThres);
        System.out.println("\tminCountForTesting: " + minCountForTesting);
        System.out.println("\tinitSkimRate: " + initSkimRate);
        System.out.println();

        int[] taxaMapAlignmentToCHBT = new int[theCHBT.getTaxaCount()];
        for (int t = 0; t < taxaMapAlignmentToCHBT.length; t++) {
            taxaMapAlignmentToCHBT[t] = gMap[0].getIdGroup().whichIdNumber(theCHBT.getTaxaName(t));
        }
        int mapsSomewhere = 0;
        int count = 0;
        int tooFewTaxa = 0;
        int failedSkim = 0;
        try{
        DataOutputStream fw=null;
        if(outFile!=null) fw=new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outFile),4000000));
        String header = "ContigNum\tGeneticChr\tGeneticStart\tGeneticEnd\tGeneticMean\tGeneticSpan\tB73count\tMo17count\tSegnType\tRecomb";
        if(fw!=null) {fw.writeBytes(header); fw.writeBytes("\n");}
        else {System.out.println(header);}
        for (int i = 0; i < theCHBT.getTotalNumCtgs(); i++) {
            if (theCHBT.getTaxaCountForHaplotype(i) < minCountForTesting) {  // the number of Taxa the contigHaplotype was seen in
                count++;  tooFewTaxa++;  continue;
            }
            int bestSkimAlign=-1;
            double bestSkimP=1;
            for (int currAlign = 0; currAlign < gMap.length; currAlign++) {
                GeneticMapResult gmr=findBestHitOnChromosome(gMap[currAlign], theCHBT, taxaMapAlignmentToCHBT, i, minCountForTesting, initSkimRate);
                if(gmr.bestP<bestSkimP) {bestSkimP=gmr.bestP; bestSkimAlign=currAlign;}
            }  //end of chromosomes
            if(bestSkimAlign<0) {count++;  ++failedSkim;  continue;}  //nothing was found, probably sample size
            GeneticMapResult gmr=findBestHitOnChromosome(gMap[bestSkimAlign], theCHBT, taxaMapAlignmentToCHBT, i, minCountForTesting, 1);
            if (gmr.bestP < recombThres) {
                mapsSomewhere++;
                StringBuilder sb=new StringBuilder();
                sb.append(theCHBT.getContigID(i) + "\t" + gMap[bestSkimAlign].getLocus(gmr.startPos)
                        + "\t" + gMap[bestSkimAlign].getPositionInLocus(gmr.startPos)
                        + "\t" + gMap[bestSkimAlign].getPositionInLocus(gmr.endPos)
                        + "\t" + (gMap[bestSkimAlign].getPositionInLocus(gmr.endPos)+gMap[bestSkimAlign].getPositionInLocus(gmr.startPos))/2
                        + "\t" + (gMap[bestSkimAlign].getPositionInLocus(gmr.endPos)-gMap[bestSkimAlign].getPositionInLocus(gmr.startPos))
                        + "\t");
                String source = (gmr.refCnt > gmr.altCnt) ? "REF" : "ALT";
                sb.append(gmr.refCnt + "\t" + gmr.altCnt + "\t" + source + "\t" + gmr.bestP);
                if(fw!=null) {fw.writeBytes(sb.toString()); fw.writeBytes("\n");}
                else {System.out.println(sb.toString());}
            }

            count++;
            //  System.out.println(flankP+"\t"+bestNonChrP+"\t"+bestNonChr+"\t"+bestNonChrPos+"\t"+workReadPos.toString());
            if (count % 1000 == 0) {
                System.out.println("  Tested contigs: " + count + "    mapping: " + mapsSomewhere);
                if(fw!=null) fw.flush();
            }
        } //end of reads
            if(fw!=null) fw.close();
        }  //end of try
        catch(Exception e) {
            System.out.println("Error writing to file:"+e.getMessage());
            e.printStackTrace();
        }
        System.out.println("Total: " + count + "   mapping: " + mapsSomewhere + "   failedSkim: " + failedSkim + "   tooFewTaxa: " + tooFewTaxa);
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

    private static GeneticMapResult findBestHitOnChromosome(Alignment gMap, ContigHaplotypesByTaxa theCHBT,
            int[] taxaMapAlignmentToRBT, int chbtIndex, int minCountForTesting, int stepSize) {

        int startBest = -1, endBest = -1, altCntBest = -1, refCntBest = -1;
        double bestP = 1;
        for (int geneticSite = 0; geneticSite < gMap.getSiteCount(); geneticSite+=stepSize) {
            int refCnt = 0, altCnt = 0;
            for (int t = 0; t < taxaMapAlignmentToRBT.length; t++) {
                if (taxaMapAlignmentToRBT[t] < 0) {
                    continue;
                }
                if (theCHBT.getHaplotypeCountForTaxa(chbtIndex, t) < 1) {
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
        GeneticMapResult gmr=new GeneticMapResult(startBest, endBest, refCntBest, altCntBest, bestP, gMap.getLocus(startBest));
        return gmr;
    }

}
