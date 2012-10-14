/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.gbs.pipeline;

import java.awt.Frame;
import java.io.File;
import javax.swing.ImageIcon;
import net.maizegenetics.gbs.util.ArgsEngine;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import org.apache.log4j.Logger;
/**
 *
 * @author Fei Lu
 */
public class UCreatWorkingDirPlugin extends AbstractPlugin {
	private ArgsEngine engine = null;
	private Logger logger = Logger.getLogger(UCreatWorkingDirPlugin.class);
	private String parentDir=null;
	static String[] childDir = {"Illumina", "key", "tagCounts", "mergedTagCounts", "tagPair", "tagsByTaxa", "mapInfo", "hapMap"};

	public UCreatWorkingDirPlugin () {
        super(null, false);
    }

    public UCreatWorkingDirPlugin (Frame parentFrame) {
        super(parentFrame, false);
    }

    private void printUsage (){
        logger.info(
             "\n\nUsage is as follows:\n"
			+ " -w  Working directory to contain subdirectories\n"
        );
    }

	public DataSet performFunction (DataSet input) {
		File[] FileOfData = new File[childDir.length];
		for (int i = 0; i < childDir.length; i++) {
			FileOfData[i] = new File(parentDir, childDir[i]);
			FileOfData[i].mkdirs();
		}
		return null;
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
