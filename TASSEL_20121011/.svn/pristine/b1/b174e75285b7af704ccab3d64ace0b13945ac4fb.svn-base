/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.gbs.pipeline;

import net.maizegenetics.gbs.tagdist.UTagPairs;
import java.awt.Frame;
import java.io.File;
import javax.swing.ImageIcon;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.gbs.util.ArgsEngine;
import org.apache.log4j.Logger;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByte;

/**
 *
 * @author Fei Lu
 */
public class UTBTToMapInfoPlugin extends AbstractPlugin {
	private ArgsEngine engine = null;
	private Logger logger = Logger.getLogger(UTBTToMapInfoPlugin.class);
	private String parentDir=null;

	public UTBTToMapInfoPlugin () {
		super(null, false);
	}

	public UTBTToMapInfoPlugin (Frame parentFrame) {
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
		String tagPairFileS = new File (pd, UCreatWorkingDirPlugin.childDir[4]).getAbsolutePath() + "/tagPair.tps";
		String TBTFileS = new File (pd, UCreatWorkingDirPlugin.childDir[5]).getAbsolutePath() + "/tbt.bin";
		String mapInfoFileS = new File (pd, UCreatWorkingDirPlugin.childDir[6]).getAbsolutePath() + "/mapInfo.bin";
		this.creatMapInfo(TBTFileS, tagPairFileS, mapInfoFileS);
		return null;
	}

	private void creatMapInfo (String TBTFileS, String tagPairFileS, String mapInfoFileS) {
		TagsByTaxaByte tbt = new TagsByTaxaByte (TBTFileS, FilePacking.Byte);
		UTagPairs tp = new UTagPairs (tagPairFileS);
		UMapInfo umi = new UMapInfo (tbt, tp);
		umi.writeMapInfo(mapInfoFileS);
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
