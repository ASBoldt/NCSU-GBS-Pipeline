/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.gbs.pipeline;

import java.awt.Frame;
import java.io.File;
import javax.swing.ImageIcon;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.gbs.util.ArgsEngine;
import org.apache.log4j.Logger;
/**
 *
 * @author Fei Lu
 */
public class UMapInfoToHapMapPlugin extends AbstractPlugin {
	private ArgsEngine engine = null;
	private Logger logger = Logger.getLogger(UTBTToMapInfoPlugin.class);
	private String parentDir=null;
	private double mnMAF = 0.05;
	private double mxMAF = 0.5;
	private double mnCall = 0;
	private double mxCall = 1;

	public UMapInfoToHapMapPlugin () {
		super(null, false);
	}

	public UMapInfoToHapMapPlugin (Frame parentFrame) {
		super(parentFrame, false);
    }

	private void printUsage (){
        logger.info(
             "\n\nUsage is as follows:\n"
			+ " -w      Working directory to contain subdirectories\n"
			+ " -mnMAF  Mimimum minor allele frequency. Default: 0.05\n"
			+ " -mxMAF  Maximum minor allele frequency. Default: 0.5\n"
			+ " -mnC    Minimum call rate"
			+ " -mxC    Maximum call rate. Default: 1\n"
        );
    }

	public DataSet performFunction(DataSet input) {
		File pd = new File (parentDir);
		String mapInfoFileS = new File (pd, UCreatWorkingDirPlugin.childDir[6]).getAbsolutePath() + "/mapInfo.bin";
		String hapMapDirS = new File (pd, UCreatWorkingDirPlugin.childDir[7]).getAbsolutePath();
		UMapInfo umi = new UMapInfo (mapInfoFileS);
		umi.writeHapAll(hapMapDirS, mnMAF, mxMAF, mnCall, mxCall);
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
			engine.add("-mnMAF", "--minimum-MAF", true);
			engine.add("-mxMAF", "--maximum-MAF", true);
			engine.add("-mnC", "--minimum-callRate", true);
			engine.add("-mxC", "--maximum-callRate", true);
            engine.parse(args);
        }
		if (engine.getBoolean("-w")) { parentDir = engine.getString("-w");}
        else{ printUsage(); throw new IllegalArgumentException("Please specify the working directory."); }
		if (engine.getBoolean("-mnMAF")) { mnMAF = Double.parseDouble(engine.getString("-mnMAF"));}
		if (engine.getBoolean("-mxMAF")) { mxMAF = Double.parseDouble(engine.getString("-mxMAF"));}
		if (engine.getBoolean("-mnC")) { mnCall = Double.parseDouble(engine.getString("-mnC"));}
		if (engine.getBoolean("-mxC")) { mxCall = Double.parseDouble(engine.getString("-mxC"));}
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
