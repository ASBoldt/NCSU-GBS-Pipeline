/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.gbs.pipeline;

import java.awt.Frame;
import java.io.File;
import java.io.FilenameFilter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.TreeSet;
import javax.swing.ImageIcon;
import net.maizegenetics.gbs.util.ArgsEngine;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.gbs.tagdist.TagCountMutable;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import org.apache.log4j.Logger;

/**
 *
 * @author Fei Lu
 */
public class UMergeTaxaTagCountPlugin extends AbstractPlugin {
	private ArgsEngine engine = null;
	private String parentDir = null;
	private boolean ifMerge = true;
	private int minCount = 5;
	private int maxTagOfTaxa = 10000000;
	private int maxTagOfAll = 60000000;

	private Logger logger = Logger.getLogger(UMergeTaxaTagCountPlugin.class);

	public UMergeTaxaTagCountPlugin() {
        super(null, false);
    }

    public UMergeTaxaTagCountPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

	private void printUsage (){
        logger.info(
             "\n\nUsage is as follows:\n"
			+ " -w  Working directory to contain subdirectories\n"
			+ " -t  Option to merge taxa (y/n). Default: y\n"
			+ " -m  Maximum tag number in the merged TagCount file. Default: 60000000"
			+ " -c  Minimum count of a tag must be present to be output. Default: 5"
        );
    }

	public DataSet performFunction (DataSet input) {
		File pd = new File (parentDir);
		String tagCountOfTaxaDirS = new File (pd, UCreatWorkingDirPlugin.childDir[2]).getAbsolutePath();
		if (ifMerge) mergeTaxa(tagCountOfTaxaDirS);
		String mergedTagCountOfAllS = new File (pd, UCreatWorkingDirPlugin.childDir[3]).getAbsolutePath() + "/mergedAll.cnt";
		this.mergeAllTagCount(tagCountOfTaxaDirS, mergedTagCountOfAllS);
		return null;
	}

	private class sortFilebyTaxaName implements Comparator <File> {
		public int compare (File o1, File o2) {
			return o1.getName().split("_")[0].compareTo(o2.getName().split("_")[0]);
		}
	}
	
	public void mergeTaxa (String tagCountOfTaxaDirS) {
		File[] allfile = new File(tagCountOfTaxaDirS).listFiles(new ListFilter("cnt", 2));
		Arrays.sort(allfile, new sortFilebyTaxaName());
		TreeSet<String> uniqueTaxaNameSet = new TreeSet();
		String[] taxaNameArray = new String[allfile.length];
		for (int i = 0; i < allfile.length; i++) {
			String[] temp = allfile[i].getName().split("_");
			taxaNameArray[i] = temp[0];
			System.out.println(temp[0]);
			uniqueTaxaNameSet.add(temp[0]);
		}
		String[] uniqueTaxaNameArray = uniqueTaxaNameSet.toArray(new String[uniqueTaxaNameSet.size()]);
		for (int i = 0; i < uniqueTaxaNameArray.length; i++) {
			ArrayList<Integer> idenList = new ArrayList();
			int hit = Arrays.binarySearch(taxaNameArray, uniqueTaxaNameArray[i]);
			while (hit > 0 && uniqueTaxaNameArray[i].equals(taxaNameArray[hit-1])) {
				hit--;
			}

			while (hit < taxaNameArray.length && uniqueTaxaNameArray[i].equals(taxaNameArray[hit])) {
				idenList.add(hit);
				hit++;
			}
			if (idenList.size() > 1) {
				Integer[] iden = idenList.toArray(new Integer[idenList.size()]);
				TagCounts primeTC = new TagCounts (allfile[iden[0]].getAbsolutePath(), FilePacking.Bit);
				TagCountMutable primeTCM = new TagCountMutable (primeTC, maxTagOfTaxa);

				for (int j = 1; j < iden.length; j++) {
					TagCounts tc = new TagCounts(allfile[iden[j]].getAbsolutePath(), FilePacking.Bit);
					primeTCM.addReadCounts(tc);
					if (j % 5 == 0) {
						primeTCM.collapseCounts();
						System.out.println(uniqueTaxaNameArray[i] + " currently contains " + primeTCM.getCurrentSize() + " tags.");
					}
				}
				primeTCM.collapseCounts();
				System.out.println(uniqueTaxaNameArray[i] + " finally contains " + primeTCM.getCurrentSize() + " tags.");
				String mergedFileName = new File (allfile[iden[0]].getParent(), uniqueTaxaNameArray[i] + "_merged.cnt").toString();
				primeTCM.writeTagCountFile(mergedFileName, FilePacking.Bit, 1);
				for (int j = 0; j < iden.length; j++) allfile[iden[j]].delete();
				System.out.println(idenList.size() + " tagCount files of " + uniqueTaxaNameArray[i] + " are merged in " + mergedFileName);
			}
		}
	}

	public void mergeAllTagCount (String tagCountOfTaxaDirS, String mergedTagCountOfAllS) {
		File[] all = new File(tagCountOfTaxaDirS).listFiles(new ListFilter("cnt", 2));
		TagCounts primeTC = new TagCounts(all[0].getAbsolutePath(), FilePacking.Bit);
		TagCountMutable primeTCM = new TagCountMutable (primeTC, maxTagOfAll);
		for (int i = 1; i < all.length; i++) {
			TagCounts tc = new TagCounts(all[i].getAbsolutePath(), FilePacking.Bit);
			if (primeTCM.getCurrentSize()+tc.getSize() > maxTagOfAll) primeTCM.removeRareTag(2);
			primeTCM.addReadCounts(tc);
			if (i % 10 == 0) {
				primeTCM.collapseCounts();
				System.out.println("Merged tagCount file of all taxa currently contains " + primeTCM.getCurrentSize() + " tags.");
			}
		}
		primeTCM.collapseCounts();
		System.out.println("Merged tagCount file of all taxa finally contains " + primeTCM.getCurrentSize() + " tags.");
		primeTCM.writeTagCountFile(mergedTagCountOfAllS, FilePacking.Bit, minCount);
	}

	@Override
    public void setParameters(String [] args) {
		if(args.length==0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
		if(engine == null){
			engine = new ArgsEngine();
            engine.add("-w", "--working-directory", true);
			engine.add("-t", "--merging-taxa", true);
			engine.add("-c", "--min-count", true);
			engine.add("-m", "--max-tagNumber", true);
            engine.parse(args);
        }
		if (engine.getBoolean("-w")) { parentDir = engine.getString("-w");}
        else{ printUsage(); throw new IllegalArgumentException("Please specify the working directory."); }
		if (engine.getBoolean("-t")) {
			String m = engine.getString("-t");
			if (m.equalsIgnoreCase("n")) ifMerge = false;
			else ifMerge = true;
		}
		if (engine.getBoolean("-c")) minCount = Integer.parseInt(engine.getString("-c"));
		if (engine.getBoolean("-m")) maxTagOfAll = Integer.parseInt(engine.getString("-m"));
    }

	private class ListFilter implements FilenameFilter {
		String filterS;
		int mode;
		ListFilter (String filterS, int mode) {
			this.filterS = filterS;
			this.mode = mode;
		}
		public boolean accept(File dir, String name) {
			if (mode == 0) {
				return name.startsWith(filterS);
			}
			else if (mode == 1) {
				return name.contains(filterS);
			}
			else if (mode == 2) {
				return name.endsWith(filterS);
			}
			return false;
		}
	}

	public ImageIcon getIcon() {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	public String getButtonName() {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	public String getToolTipText() {
		throw new UnsupportedOperationException("Not supported yet.");
	}

}
