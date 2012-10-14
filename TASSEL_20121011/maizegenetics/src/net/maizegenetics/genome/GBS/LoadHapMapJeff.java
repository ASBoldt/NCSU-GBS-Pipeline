/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.CombineAlignment;
import net.maizegenetics.pal.alignment.ImportUtils;

/**
 *
 * @author glaubitz
 */
public class LoadHapMapJeff {
    static String baseDir="/usr/local/maizediv/";
    static String inFrameworkMap=baseDir+"illumina/434LFAAXX/GBS_framework/commonReadMap_26058_SNPs_256_RILs_hmp.txt";
    static String in55KMap=baseDir+"illumina/434LFAAXX/GBS_framework/IBM_55K_hapmap.txt";
    static String in55KMap282="/home/glaubitz/data/55K/SNP55K_maize282_hapmap_20100513.txt";
    static String in3KSNPs282="Q:/data/SNPs/281_vs_3122_SNPs_AGPv1/SNP3093_maize281_hapmap_20100614.txt";
    static String in55KHapMapV2Samples="Q:/data/SNPs/55K_chip/AGPv2/SNP55K_hapmapV2Samples_AGPv1_20100823.txt";
    static String in55kIBMframe="C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/ShiranVDs/IBM_55K_08192010.hmp.txt";

    public static void main(String[] args) {
        String poly, genoStr;
        Alignment[] framework;
        Alignment[] illumina55K;
//        framework = ((CombineAlignment)ImportUtils.readFromHapmap(inFrameworkMap)).getAlignments();
//        for (Alignment a : framework) {
//            poly = a.isAllPolymorphic() ? "true" : "false";
//            System.out.println("Locus: " + a.getLocus(0) + " Sites: " + a.getSiteCount() + " Taxa: " + a.getSequenceCount() + " allPoly: " + poly);
//            for (int i = a.getSiteCount()-1; i < a.getSiteCount(); i++) {  // ;
////                       System.out.println(i+" : "+(char)a.getMajorAllele(i)+" "+a.getMinorAllele(i)+" "+a.getMinorAlleleFrequency(i));
//                genoStr = a.getLocusName(i) + " " + a.getPositionInLocus(i) + " ";
//                for (int j = 0; j< a.getSequenceCount(); j++)  {
//                    genoStr += a.getBaseChar(j, i);
//                }
//                System.out.println(genoStr);
//            }
//        }
        illumina55K = ((CombineAlignment)ImportUtils.readFromHapmap(in55kIBMframe)).getAlignments();
        for (Alignment a : illumina55K) {
            poly = a.isAllPolymorphic() ? "true" : "false";
            System.out.println("Locus: " + a.getLocus(0) + " Sites: " + a.getSiteCount() + " Taxa: " + a.getSequenceCount() + " allPoly: " + poly);
            for (int i = a.getSiteCount()-1; i < a.getSiteCount(); i++) {  // ;
                       System.out.println(i+" : "+(char)a.getMajorAllele(i)+" "+a.getMinorAllele(i)+" "+a.getMinorAlleleFrequency(i));
                genoStr = a.getLocusName(i) + " " + a.getPositionInLocus(i) + " ";
                for (int j = 0; j< a.getSequenceCount(); j++)  {
                    genoStr += a.getBaseChar(j, i);
                }
                System.out.println(genoStr);
            }
        }
    }
}
