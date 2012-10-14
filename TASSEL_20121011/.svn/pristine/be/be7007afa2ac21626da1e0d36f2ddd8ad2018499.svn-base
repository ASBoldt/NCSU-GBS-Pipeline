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
import java.io.File;


/**
 *
 * @author Qi
 */
public class MergeDuplicateTaxa_vcf_Plugin extends AbstractPlugin{
    String inputfile;
    String outputfile;
    int startChr=1, endChr=10;
    private static ArgsEngine myArgsEngine = null;
    int alleleskept=3;
    private static Logger myLogger=Logger.getLogger(MergeDuplicateTaxa_vcf_Plugin.class);
    
    public MergeDuplicateTaxa_vcf_Plugin(){
        super(null, false);
    }
    public MergeDuplicateTaxa_vcf_Plugin(Frame parentFrame) {
        super(parentFrame, false);
    }
    
        public DataSet performFunction(DataSet input) {
        for (int chr=startChr; chr<=endChr; chr++)
        {
           String myinfile = inputfile;
           
           
    
           myinfile = myinfile.replace("+", ""+chr);
           File findFile = new File(myinfile);
           if (findFile.isFile())
           {
                String myoutfile = outputfile;
                myoutfile = myoutfile.replace("+", ""+chr);
                vcf myVCF = new vcf(myinfile, alleleskept);
                myVCF.merge_duplicate_taxa();
                myVCF.write_vcf(myoutfile);
           }
           else
           {
               
           }
            

        }
        

        return null;
    }
        
private void printUsage() {
        myLogger.info(
            "\nUsage is as follows:\n"
            + "-i       Input vcf file template. using '+' for wildcard\n"
            + "-o       Output vcf file. using '+' for wildcard\n" 
            + "-ak      Maximum number of alleles that are kept for each marker across the population , default: 3\n"
            + "-s       Start chr\n"
            + "-e       End Chr\n"
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
            myArgsEngine.add("-s", "--startChromosome", true);
            myArgsEngine.add("-e", "--endChromosome", true);
            
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
        
        if (myArgsEngine.getBoolean("-s"))            startChr = Integer.parseInt(myArgsEngine.getString("-s"));
        if (myArgsEngine.getBoolean("-e"))            endChr = Integer.parseInt(myArgsEngine.getString("-e"));
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
