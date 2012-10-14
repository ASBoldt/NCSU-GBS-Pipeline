/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pipeline;
import java.io.File;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByteFileMap;
import net.maizegenetics.gbs.util.ArgsEngine;
import org.apache.log4j.Logger;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import java.awt.Frame;
import javax.swing.ImageIcon;
import org.biojava3.core.util.ConcurrencyTools;
/**
 *
 * @author Qi
 */
public class tbt2vcfPlugin extends AbstractPlugin{
        static int minTaxaCnt=0;  // a tag must show up in GREATER THAN minTaxaCnt taxa to be included in a sequence alignment
    static int maxSize=200000;  //normally 200K;
    private double minF=-2.0, minMAF=0.0;
    private int minMAC=10;
    private double mxMAF = 1.0;
    boolean top2tags = false;
//    static boolean ignoreTriallelic=false;
    //private boolean inclRare=false;  // false = only call the two most common alleles at a site
    //private boolean inclGaps=false;  // false = ignore sites where the major or the 1st minor alleles are gaps
    //private boolean isUpdateTOPM=false;
    //private final static int maxSNPsPerLocus=64;
    //private final static int maxAlignmentSize=150;
    static double defaultMinPropTaxaWithLocus=0.0;
    TagsOnPhysicalMap theTOPM = null;
    TagsByTaxaByteFileMap theTBT=null;
    File inputFile=null;
    private String inTOPMFile=null;
    //private String outTOPMFile=null;
    String outputFilePrefix = null;
    String outHapMap=null;
    int startChr=0;
    int endChr=0;
    private static ArgsEngine myArgsEngine = null;
    int minTaxaWithLocus;
    private static Logger myLogger=Logger.getLogger(tbt2vcfPlugin.class);
    int alleleskept=3;
        
    public tbt2vcfPlugin(){
        super(null, false);
    }
    public tbt2vcfPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }
    
        public DataSet performFunction(DataSet input) {
        myLogger.info("Convert tbt to vcf file");
        theTOPM.sortTable(true);
        theTOPM.printRows(5, true, true);
        System.out.println(String.format("StartChr:%d EndChr:%d %n", startChr, endChr));
        for(int i = startChr; i <= endChr; i++) {
             System.out.println("Processing chromosome "+i+"..."); 
             String out=outHapMap+".c"+i;
             System.out.println("Creating Mutable Alignment");
             vcf theVCF=new vcf(theTBT, theTOPM, i, minTaxaWithLocus, minMAF, mxMAF, alleleskept, 0);
             theVCF.write_vcf(out);
             System.out.println("Finish chromosome "+i+"..."); 
        }
        ConcurrencyTools.shutdown();
        return null;
    }
        
private void printUsage() {
         myLogger.info(
            "\nUsage is as follows:\n"
            + "-i       Input .tbt file\n"
            + "-o       Output directory (default current directory)\n" 
            + "-m       TagsOnPhysicalMap file containing genomic position of tags\n" 
            + "-mnMAF   Minimum minor allele frequency (default: "+minMAF+")\n"
            + "-mnLCov  Minimum locus coverage (proportion of Taxa with a genotype) (default: "+defaultMinPropTaxaWithLocus+")\n"
            + "-ak      Maximum number of alleles that are kept for each marker across the population , default: 3\n"
            + "-s       Start chromosome\n"
            + "-e       End chromosome"
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
            myArgsEngine.add("-o", "--output-directory", true);
            myArgsEngine.add("-m", "--physical-map", true);
            myArgsEngine.add("-mnMAF", "--minMinorAlleleFreq", true);
            myArgsEngine.add("-mnLCov", "--minLocusCov", true);
            myArgsEngine.add("-ak", "--AllelesKept", true);
            myArgsEngine.add("-s", "--start-chromosome", true);
            myArgsEngine.add("-e", "--end-chromosome", true);
        }
        myArgsEngine.parse(args);

        if (myArgsEngine.getBoolean("-i")) {
            inputFile = new File(myArgsEngine.getString("-i"));
            if (!inputFile.exists()) {
                printUsage();
                throw new IllegalArgumentException("Can't find the input file.");
            }
            outputFilePrefix = inputFile.getParentFile().getName();
            theTBT=new TagsByTaxaByteFileMap(myArgsEngine.getString("-i"));
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify an input file.");
        }

        //Set output directory and use the input TagsByTaxa filename for the output file prefix
        if (myArgsEngine.getBoolean("-o")) {
            outHapMap = myArgsEngine.getString("-o")+File.separator+outputFilePrefix;
        } else {
            outHapMap = inputFile.getParent()+File.separator+outputFilePrefix;
        }

        if (myArgsEngine.getBoolean("-m")) {
             inTOPMFile=myArgsEngine.getString("-m");
             boolean loadBinary=(inTOPMFile.endsWith(".txt"))?false:true;
             theTOPM = new TagsOnPhysicalMap(inTOPMFile, loadBinary);
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify a physical map file.");
        }

        if (myArgsEngine.getBoolean("-mnMAF")) {minMAF = Double.parseDouble(myArgsEngine.getString("-mnMAF"));}

        minTaxaWithLocus=(int)Math.round(theTBT.getTaxaCount()*defaultMinPropTaxaWithLocus);
        if (myArgsEngine.getBoolean("-mnLCov")) {
            double minPropTaxaWithLocus= Double.parseDouble(myArgsEngine.getString("-mnLCov"));
            minTaxaWithLocus=(int)Math.round(theTBT.getTaxaCount()*minPropTaxaWithLocus);
        }
        if (myArgsEngine.getBoolean("-ak")) {alleleskept = Integer.parseInt(myArgsEngine.getString("-ak"));}
       
        if (myArgsEngine.getBoolean("-s")) {
            startChr = Integer.parseInt(myArgsEngine.getString("-s"));
            if(startChr==0){
                printUsage();
                throw new IllegalArgumentException("Error: start chromosome is 0.");
            }
        } else {
                printUsage();
                throw new IllegalArgumentException("Please specify start and end chromosome numbers.");
        }
        if (myArgsEngine.getBoolean("-e")) {
            endChr = Integer.parseInt(myArgsEngine.getString("-e"));
        } else {
                printUsage();
                throw new IllegalArgumentException("Please specify start and end chromosome numbers.");
        }
        if (endChr-startChr<0) {
            printUsage();
            throw new IllegalArgumentException("Error: The start chromosome is higher than the end chromosome.");
        }
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
