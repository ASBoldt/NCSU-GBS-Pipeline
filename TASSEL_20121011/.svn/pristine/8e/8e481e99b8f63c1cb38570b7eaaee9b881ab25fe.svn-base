/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pipeline;
import net.maizegenetics.gbs.util.ArgsEngine;
import org.apache.log4j.Logger;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import java.awt.Frame;
import javax.swing.ImageIcon;


/**
 *
 * @author Qi
 */
public class MergeDuplicateSNP_vcf_Plugin extends AbstractPlugin{
    String inputfile;
    String outputfile;
    private static ArgsEngine myArgsEngine = null;
    int alleleskept=3;
    private static Logger myLogger=Logger.getLogger(MergeDuplicateSNP_vcf_Plugin.class);
    
    public MergeDuplicateSNP_vcf_Plugin(){
        super(null, false);
    }
    public MergeDuplicateSNP_vcf_Plugin(Frame parentFrame) {
        super(parentFrame, false);
    }
    
        public DataSet performFunction(DataSet input) {
        vcf myVCF = new vcf(inputfile, alleleskept);
        myVCF.remove_duplicate();
        myVCF.write_vcf(outputfile);
        return null;
    }
        
private void printUsage() {
        myLogger.info(
            "\nUsage is as follows:\n"
            + "-i       Input vcf file\n"
            + "-o       Output vcf file\n" 
            + "-ak      Maximum number of alleles that are kept for each marker across the population , default: 3\n"
            );
        }


    public void setParameters(String[] args) {
        if(args.length==0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
 
        if (myArgsEngine == null) {
            myArgsEngine = new ArgsEngine();
            myArgsEngine.add("-i", "--input-file", true);
            myArgsEngine.add("-o", "--output-file", true);
            myArgsEngine.add("-ak", "--alleleskept", true);
        }
        myArgsEngine.parse(args);

        if (myArgsEngine.getBoolean("-i")) {
            inputfile = myArgsEngine.getString("-i");
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify an input file.");
        }
        
        if (myArgsEngine.getBoolean("-o")) {
            outputfile = myArgsEngine.getString("-o");
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify an output file.");
        }
        
         if (myArgsEngine.getBoolean("-ak")) {alleleskept = Integer.parseInt(myArgsEngine.getString("-ak"));}
    }
    
    @Override
    public ImageIcon getIcon() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getButtonName() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getToolTipText() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
   

    
}
