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
public class UTagCountToTagPairPlugin extends AbstractPlugin {
	private ArgsEngine engine = null;
	private Logger logger = Logger.getLogger(UTagCountToTagPairPlugin.class);
	private String parentDir=null;
	private double etr = 0.03;

	public UTagCountToTagPairPlugin () {
        super(null, false);
    }

    public UTagCountToTagPairPlugin (Frame parentFrame) {
        super(parentFrame, false);
    }

    private void printUsage (){
        logger.info(
             "\n\nUsage is as follows:\n"
			+ " -w  Working directory to contain subdirectories\n"
			+ " -e  Error tolerance rate in the network filter. Default: 0.03\n"
        );
    }

	public DataSet performFunction(DataSet input) {
		File pd = new File (parentDir);
		String mergedTagCountOfAllS = new File (pd, UCreatWorkingDirPlugin.childDir[3]).getAbsolutePath() + "/mergedAll.cnt";
		String tagPairS = new File (pd, UCreatWorkingDirPlugin.childDir[4]).getAbsolutePath() + "/tagPair.tps";
		UNetworkFilter unf = new UNetworkFilter (mergedTagCountOfAllS, etr, tagPairS);
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
			engine.add("-e", "--error-tolerance-rate", true);
            engine.parse(args);
        }
		if (engine.getBoolean("-w")) { parentDir = engine.getString("-w");}
        else{ printUsage(); throw new IllegalArgumentException("Please specify the working directory."); }
		if (engine.getBoolean("-e")) { etr = Double.parseDouble(engine.getString("-e"));}
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
