/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS.switchgrass;

import java.util.Arrays;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AnnotatedAlignmentUtils;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.distance.DistanceMatrix;
import net.maizegenetics.pal.distance.WriteDistanceMatrix;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.pal.popgen.LinkageDisequilibrium;
import net.maizegenetics.pal.tree.Tree;



/**
 *
 * @author fl262
 */
public class LDMatrixOutput {
	static double minPresense=0.5;
        static double minMAF=0.05;
        static double maxHeterozygosity=0.7;

	public static void main (String[] args) {
		 Alignment a = ImportUtils.readFromHapmap("M:/Snps_SWGPop1&Pop2_removed.txt");
		 System.out.println("File Loaded");
		 System.out.println("SiteFiltered Alignment  Taxa:" + a.getSequenceCount() + " Sites:" + a.getSiteCount());
		 int minCount=(int)Math.round(a.getSequenceCount()*minPresense);
		 a=AnnotatedAlignmentUtils.removeSitesBasedOnFreqIgnoreGapsMissing(a,minMAF,minCount);
//		 int[] goodLowHetSites = AlignmentFilterByGBSUtils.getLowHetSNPs(a, false, maxHeterozygosity, minCount);
//        String[] taxaFilter={"BLANK"};
//       // IdGroup keepTaxa=AlignmentFilterByGBSUtils.getFilteredIdGroupByName(a.getIdGroup(), taxaFilter, false);
//         a = FilterAlignment.getInstance(a, goodLowHetSites);
		 System.out.println("SiteFiltered Alignment  Taxa:" + a.getSequenceCount() + " Sites:" + a.getSiteCount());
		  System.out.println("Filter Sites");
		  ExportUtils.writeToHapmap(a, false, "M:/temp.hmp.txt", '\t');
		  System.out.println("HapMap written");
		  a = ImportUtils.readFromHapmap("M:/temp.hmp.txt");
		   System.out.println("LD Starting");
		 LinkageDisequilibrium theLD=new LinkageDisequilibrium( a,
            20, 50, LinkageDisequilibrium.testDesign.All);
		 theLD.run();
		  System.out.println("LD Complete");
		 double[][] r2mat=new double[a.getSiteCount()][a.getSiteCount()];
		 String[] siteNames=a.getSNPIDs();
		 double[] minP=new double[a.getSiteCount()];
		 Arrays.fill(minP, 1);
		 for (int i = 0; i < a.getSiteCount(); i++) {
			 for (int j = 0; j < i; j++) {
				 double val=1-theLD.getRSqr(i, j);
				 if(Double.isNaN(val)) val=1.0;
				 r2mat[i][j]=r2mat[j][i]=val;
				 double p=theLD.getPVal(i, j);
				 if(Double.isNaN(p)) continue;
				 if(p<minP[i]) minP[i]=p;
				  if(p<minP[j]) minP[j]=p;
			 }
		}
		for (int i = 0; i < a.getSiteCount(); i++) {
			System.out.printf("%d %g %n",i,minP[i]);
		}
		DistanceMatrix dm=new DistanceMatrix(r2mat, new SimpleIdGroup(siteNames));
		Tree theTree = new net.maizegenetics.pal.tree.NeighborJoiningTree(dm);
		System.out.println(theTree.toString());
		theTree = new net.maizegenetics.pal.tree.UPGMATree(dm);
//		System.out.println(theTree.toString());
		WriteDistanceMatrix.saveDelimitedDistanceMatrix(dm, "M:/sitedist.txt");
	}
}
