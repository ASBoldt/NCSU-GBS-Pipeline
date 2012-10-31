/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.gbs.pipeline;

import java.awt.Frame;
import java.io.File;
import java.io.FilenameFilter;
import java.util.Arrays;
import java.util.Comparator;
import javax.swing.ImageIcon;
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.gbs.tagdist.Tags;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByte;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.gbs.util.ArgsEngine;
import org.apache.log4j.Logger;

/**
 *
 * @author Fei Lu
 */
public class UTagPairToTBTPlugin extends AbstractPlugin {
	private ArgsEngine engine = null;
	private Logger logger = Logger.getLogger(UTagPairToTBTPlugin.class);
	private String parentDir=null;

	public UTagPairToTBTPlugin () {
		super(null, false);
	}

	public UTagPairToTBTPlugin (Frame parentFrame) {
		super(parentFrame, false);
    }

	private void printUsage (){
        logger.info(
             "\n\nUsage is as follows:\n"
			+ " -w  Working directory to contain subdirectories\n"
        );
    }

	public DataSet performFunction(DataSet input) {
		File pd = new File (parentDir);
		String tagCountOfTaxaDirS = new File (pd, UCreatWorkingDirPlugin.childDir[2]).getAbsolutePath();
		String tagPairFileS = new File (pd, UCreatWorkingDirPlugin.childDir[4]).getAbsolutePath() + "/tagPair.tps";
		String TBTFileS = new File (pd, UCreatWorkingDirPlugin.childDir[5]).getAbsolutePath() + "/tbt.bin";
		this.creatTBT(tagCountOfTaxaDirS, tagPairFileS, TBTFileS);
		return null;
	}

	private class sortFilebyTaxaName implements Comparator <File> {
		public int compare (File o1, File o2) {
			return o1.getName().split("_")[0].compareTo(o2.getName().split("_")[0]);
		}
	}

	private void creatTBT (String tagCountOfTaxaDirS, String tagPairFileS, String TBTFileS) {
		TagCounts masterTag = new TagCounts (tagPairFileS,  FilePacking.Bit);
		masterTag.sort();
		File[] tcFiles = new File(tagCountOfTaxaDirS).listFiles(new ListFilter("cnt", 2));
		Arrays.sort(tcFiles, new sortFilebyTaxaName());
		String[] taxaNames = new String[tcFiles.length];
		for (int i = 0; i < tcFiles.length; i++) {
			taxaNames[i] = tcFiles[i].getName().split("\\.")[0];
		}
		TagsByTaxa theTBT = new TagsByTaxaByte(taxaNames, masterTag);
		for (int i = 0; i < tcFiles.length; i++) {
			TagCounts tc = new TagCounts(tcFiles[i].toString(), FilePacking.Bit);
			for (int j = 0; j < tc.getSize(); j++) {
				long[] query = tc.getTag(j);
				int tagIndex = theTBT.getTagIndex(query);
				if (tagIndex > -1) {
					int cnt = tc.getReadCount(j);
					theTBT.addReadsToTagTaxon(tagIndex, i, cnt);
				}
			}
		}
		theTBT.writeDistFile(new File(TBTFileS),FilePacking.Byte, 1);
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
            engine.parse(args);
        }
		if (engine.getBoolean("-w")) { parentDir = engine.getString("-w");}
        else{ printUsage(); throw new IllegalArgumentException("Please specify the working directory."); }
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
