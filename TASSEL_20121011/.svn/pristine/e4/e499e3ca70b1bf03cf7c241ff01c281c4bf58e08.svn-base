package net.maizegenetics.gwas.jointlinkage;

import java.awt.Frame;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.regex.Pattern;

import javax.swing.ImageIcon;
import javax.swing.JOptionPane;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;

import net.maizegenetics.baseplugins.FileLoadPlugin;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;

public class StepwiseJointLinkagePlugin extends AbstractPlugin {
	private static final Logger myLogger = Logger.getLogger(StepwiseJointLinkagePlugin.class);
	private boolean parametersSet = false;
	private String genotypeFilename;
	private String phenotypeFilename;
	private String modelFilename;
	private String stepFilename;
	private String scanFilename;
	int numberOfTraits;
	int firstTraitColumn;
	int maxMarkers;
	boolean multithreaded;
	double[][] limits;
	
	public StepwiseJointLinkagePlugin(Frame parentFrame) {
		super(parentFrame, false);
		
		if (isInteractive()) {
			String msg = "Stepwise Joint Linkage is not implemented for the GUI.";
			JOptionPane.showMessageDialog(getParentFrame(), msg, "Not Implemented", JOptionPane.ERROR_MESSAGE);
			myLogger.error(msg);
			return;
		} else {
			
		}
	}

//	public static void main(String[] args) {
//		BasicConfigurator.configure();
//		StepwiseJointLinkagePlugin stepper = new StepwiseJointLinkagePlugin(null);
//		stepper.setParameters(new String[]{"-f", "/Volumes/Macintosh HD 2/project/namgbs/data/stepwise_configuration.txt"});
//		stepper.performFunction(null);
//	}
	
	/**
	 * @param filename	the name of the configuration file
	 * @return	true, if there were no errors parsing the configuration file, false otherwise
	 */
	private boolean readConfigurationFile(String filename) {
		String input;
		String[] info;
		HashMap<String, String> params = new HashMap<String, String>();
		Pattern ws = Pattern.compile("\\s+");
		try {
			BufferedReader br = new BufferedReader(new FileReader(filename));
			while ((input = br.readLine()) != null) {
				if (!(input.startsWith("#") || input.trim().length() <= 1)) {
					info = ws.split(input, 2);
					if (info[1].startsWith("\"")) { //delete enclosing double or single quotes
						info[1] = info[1].substring(1);
						if (info[1].endsWith("\"")) info[1] = info[1].substring(0, info[1].length() - 1);
					}
					else if (info[1].startsWith("'")) {
						info[1] = info[1].substring(1);
						if (info[1].endsWith("'")) info[1] = info[1].substring(0, info[1].length() - 1);
					}
					
					params.put(info[0], info[1]);
				}
			}
			
			br.close();
		} catch (IOException e) {
			myLogger.error(e);
			return false;
		}

		genotypeFilename = params.get("genotypes");
		if (genotypeFilename == null) {
			String msg = "Improper configuration file: no genotype file specified.  The configuration file must be corrected before the analysis can be run.";
			myLogger.error(msg);
			return false;
		}
		phenotypeFilename = params.get("phenotypes");
		if (phenotypeFilename == null) {
			String msg = "Improper configuration file: no phenotype file specified.  The configuration file must be corrected before the analysis can be run.";
			myLogger.error(msg);
			return false;
		}
		modelFilename = params.get("output.model");
		if (modelFilename == null) {
			String msg = "Improper configuration file: no model output file specified.  The configuration file must be corrected before the analysis can be run.";
			myLogger.error(msg);
			return false;
		}
		stepFilename = params.get("output.step");
		if (stepFilename == null) {
			String msg = "Improper configuration file: no step output file specified.  The configuration file must be corrected before the analysis can be run.";
			myLogger.error(msg);
			return false;
		}
		scanFilename = params.get("output.scan");
		if (scanFilename == null) {
			String msg = "Improper configuration file: no scan output file specified.  The configuration file must be corrected before the analysis can be run.";
			myLogger.error(msg);
			return false;
		}
		
		numberOfTraits = Integer.parseInt(params.get("number.traits"));
		firstTraitColumn = Integer.parseInt(params.get("first.trait"));
		maxMarkers = Integer.parseInt(params.get("max.markers"));
		String mt = params.get("multithreaded");
		if (mt == null || mt.toLowerCase().startsWith("t")) multithreaded = true;
		else multithreaded = false;
		
		info = params.get("limits").split(";");
		if (info.length != numberOfTraits && info.length != 1) {
			String msg = "The number of enter and exit limits must equal either the number of traits or one. The configuration file must be corrected before the analysis can be run.";
			myLogger.error(msg);
			return false;
		} else {
			limits = new double[numberOfTraits][2];
			if (info.length == 1) {
				String[] lim = info[0].split(",");
				for (int t = 0; t < numberOfTraits; t++) {
					limits[t][0] = Double.parseDouble(lim[0]);
					limits[t][1] = Double.parseDouble(lim[1]);
				}
			} else {
				for (int t = 0; t < numberOfTraits; t++) {
					String[] lim = info[t].split(",");
					limits[t][0] = Double.parseDouble(lim[0]);
					limits[t][1] = Double.parseDouble(lim[1]);
				}
			}
		}
		
		return true;
	}
	
	@Override
	public DataSet performFunction(DataSet input) {
		if (isInteractive()) {
			String msg = "Stepwise Joint Linkage is not implemented for the GUI.";
			JOptionPane.showMessageDialog(getParentFrame(), msg, "Not Implemented", JOptionPane.ERROR_MESSAGE);
			myLogger.error(msg);
		} else if (parametersSet) {
			StepwiseModelFitterForNAMJointLinkage stepper = new StepwiseModelFitterForNAMJointLinkage();
			stepper.setGenotypeFilename(genotypeFilename);
			stepper.setPhenotypeFilename(phenotypeFilename);
			stepper.setModelFilename(modelFilename);
			stepper.setStepFilename(stepFilename);
			stepper.setScanFilename(scanFilename);
			stepper.setNumberOfTraits(numberOfTraits);
			stepper.setFirstTrait(firstTraitColumn);
			stepper.setMaxNumberOfMarkers(maxMarkers);
			stepper.setLimits(limits);
			stepper.setMultithreaded(multithreaded);
			
			stepper.runAnalysis();
		} else {
			String msg = "No configuration file or incorrect configuration file. A configuration file name must be supplied using the -f option.";
			myLogger.error(msg);
		}
		return null;
	}

	@Override
	public ImageIcon getIcon() {
		// not implemented at this time
		return null;
	}

	@Override
	public String getButtonName() {
		// not implemented at this time
		return null;
	}

	@Override
	public String getToolTipText() {
		// not implemented at this time
		return null;
	}

	@Override
	public void setParameters(String[] args) {
		int narg = args.length;
		parametersSet = true;
		for (int i = 0; i < narg - 1; i++) {
			if (args[i].equalsIgnoreCase("-f")) {
				parametersSet = readConfigurationFile(args[++i]);
			}
			else if (args[i].equals("-g") || args[i].equals("-genotypes")) {
				genotypeFilename = args[++i];
			}
			else if (args[i].equals("-p") || args[i].equals("-phenotypes")) {
				phenotypeFilename = args[++i];
			}
			else if (args[i].equals("-m") || args[i].equals("-model_output")) {
				modelFilename = args[++i];
			}
			else if (args[i].equals("-s") || args[i].equals("-step_output")) {
				stepFilename = args[++i];
			}
			else if (args[i].equals("-c") || args[i].equals("-scan_output")) {
				scanFilename = args[++i];
			}
			else if (args[i].equals("-n") || args[i].equals("-number_traits")) {
				numberOfTraits = Integer.parseInt(args[++i]);
			}
			else if (args[i].equals("-t") || args[i].equals("-first_trait")) {
				firstTraitColumn = Integer.parseInt(args[++i]);
			}
			else if (args[i].equals("-x") || args[i].equals("-max_markers")) {
				maxMarkers = Integer.parseInt(args[++i]);
			}
			else if (args[i].equals("-l") || args[i].equals("-limits")) {
				limits = new double[numberOfTraits][2];
				String[] limitPairs = args[++i].split(";");
				int n = limitPairs.length;
				if (n == 1) {
					String[] pair = limitPairs[0].split(",");
					double enter = Double.parseDouble(pair[0]);
					double exit = Double.parseDouble(pair[1]);
					for (int t = 0; t < numberOfTraits; t++) {
						limits[t][0] = enter;
						limits[t][1] = exit;
					}
				} else {
					for (int t = 0; t < n; t++) {
						String[] pair = limitPairs[t].split(",");
						limits[t][0] = Double.parseDouble(pair[0]);
						limits[t][1] = Double.parseDouble(pair[1]);
					}
				}
			}
			else if (args[i].equals("-h") || args[i].equals("-multithread")) {
				i++;
				String opt = args[i].toUpperCase();
				if (opt.startsWith("T")) multithreaded = true;
				else if (opt.startsWith("F")) multithreaded = false;
			}
			else if (args[i].equals("?")) printUsage();
		}
	}
	
	public void printUsage() {
		StringBuilder usage = new StringBuilder("The StepwiseJointLinkagePlugin requires the following parameters:\n");
		usage.append("-g or -genotypes : a file containing genotypes in Hapmap format");
		usage.append("-p or -phenotypes : a file containing phenotypes");
		usage.append("-m or -model_output : model output will be appended to this file");
		usage.append("-s or -step_output : step output will be appended to this file");
		usage.append("-c or -scan_output : scan output will be appended to this file");
		usage.append("-n or -number_traits : the number of traits to be analyzed. These must be in contiguous columns in the input file.");
		usage.append("-t or -first_trait : the first column containing trait data.");
		usage.append("The following parameters are optional:");
		usage.append("-x or -max_markers : the maximum number of markers to be fit to the model.");
		usage.append("-l or -limits : the limits to be used for each trait. Comma separated enter and exit limit must be listed for each trait. Limit pairs are separated by semicolons.");
		usage.append("-h or -multithread : true, if the plugin is to use multiple cpu's for the analysis, false otherwise.");
		usage.append("? : print the parameter list.");
		
		myLogger.info(usage);
		
	}
}
