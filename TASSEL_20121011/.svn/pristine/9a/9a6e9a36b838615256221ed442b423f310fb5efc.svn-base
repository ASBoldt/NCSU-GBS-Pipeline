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
public class filter_vcf_Plugin extends AbstractPlugin{
    String inputfile;
    String outputfile;
    private double minMAF=0.01;
    private double mxMAF=1.0;
    static double defaultMinPropTaxaWithLocus=0.1;
    int minTaxaWithLocus;
    double minPropTaxaWithLocus;
    private static ArgsEngine myArgsEngine = null;
    boolean biallelic = false;
    int alleleskept = 3;
    private static Logger myLogger=Logger.getLogger(filter_vcf_Plugin.class);
    
    public filter_vcf_Plugin(){
        super(null, false);
    }
    public filter_vcf_Plugin(Frame parentFrame) {
        super(parentFrame, false);
    }
    
    public DataSet performFunction(DataSet input) {
        vcf myVCF = new vcf(inputfile, alleleskept);
        minTaxaWithLocus=(int)Math.round(minPropTaxaWithLocus * (double)myVCF.TAXA.length);
        System.out.println(minPropTaxaWithLocus);
        //System.err.println();
        myVCF.filter_vcf(minTaxaWithLocus, minMAF, mxMAF, biallelic);
        myVCF.write_vcf(outputfile);
        return null;
    }
        
private void printUsage() {
        myLogger.info(
            "\nUsage is as follows:\n"
            + "-i       Input vcf file\n"
            + "-o       Output vcf file\n" 
            + "-mnMAF   Minimum minor allele frequency (default: "+minMAF+")\n"
            + "-mxMAF   Minimum minor allele frequency (default: "+mxMAF+")\n"
            + "-mnLCov  Minimum locus coverage (proportion of Taxa with a genotype) (default: "+defaultMinPropTaxaWithLocus+")\n"
            + "-ak      Maximum number of alleles that are kept for each marker across the population , default: 3\n"
            + "-bi      BiAlleleic only (default: all)"
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
            myArgsEngine.add("-mnMAF", "--minMinorAlleleFreq", true);
            myArgsEngine.add("-mxMAF", "--mxMinorAlleleFreq", true);
            myArgsEngine.add("-mnLCov", "--minLocusCov", true);
            myArgsEngine.add("-ak", "--AllelesKept", true);
            myArgsEngine.add("-bi", "--bialleleic", false);
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
        
        if (myArgsEngine.getBoolean("-mnMAF")) {minMAF = Double.parseDouble(myArgsEngine.getString("-mnMAF"));}
        if (myArgsEngine.getBoolean("-mxMAF")) {mxMAF = Double.parseDouble(myArgsEngine.getString("-mxMAF"));}
        if (myArgsEngine.getBoolean("-mnLCov")) {
            minPropTaxaWithLocus= Double.parseDouble(myArgsEngine.getString("-mnLCov"));
        }
        if (myArgsEngine.getBoolean("-bi")) {biallelic = true;} else {biallelic = false;}
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
