package net.maizegenetics.gbs.pipeline;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import net.maizegenetics.util.DirectoryCrawler;
import org.apache.commons.lang.StringUtils;

/** @author jvh39 
 *  Stores a list of parameters that need to be recorded for each GBS analysis run.  Includes
 *  methods to load & write configuration files from disk.
 */
public class Provenance implements Serializable{
    
    /**The name of the document storing provenance for this experiment.*/
    String provenanceFileName;

    /**The name of the experiment referred to by this provenance document (by default, the name of the directory containing .qseq files).*/
    String experimentName;

    /**The full path name of the directory containing .qseq files */
    File homeDirectory;

    /**Names of all .qseq files used in the run*/
    String[] qseqFileName;

    /**Name of file with barcodes*/
    String keyFileName;

    /**Barcode sequences*/
    String[] barcode;

    /**Taxon names corresponding to discovered barcodes.  This should be in register with {@link #barcode}.*/
    String[] taxonName;

    /**Name of enzyme used to prepare GBS library*/
    String enzyme;

    /**All reads found in a <b>single qseq file</b>.  This should be in register with {@link #qseqFileName}.*/
    long[] rawReads;

    /**Reads with a known barcode and cut site overhang in a <b>single qseq file</b>.  This should be in register with {@link #qseqFileName}.*/
    long[] goodReads;

    /**"Good" reads that can be assigned to a known taxon in a <b>single qseq file</b>.  This should be in register with {@link #qseqFileName}.*/
    long[] goodMatchedReads;

    /**All reads found in <b>all qseq files</b>.*/
    long totalRawReads;

    /**Reads with a known barcode and cut site overhang in <b>all qseq files</b>.*/
    long totalGoodReads;

    /**"Good" reads that can be assigned to a known taxon in <b>all qseq files</b>.*/
    long totalGoodMatchedReads;

    /**Minimum # times a processed tag should appear in all .qseq files.*/
    int minCountPerExperiment;

    /**Minimum # times a processed tag should appear in a given line.*/
    int minCountPerTaxon;

    /**Names of count files produced by {@link QseqToTagCountPlugin}.  This should be in register with {@link #qseqFileName}.*/
    String[] countFileName;

    /**Names of TBT files produced by {@link QseqToTBTPlugin}.  This should be in register with {@link #qseqFileName}.*/
    String[] TBTFileName;

    /**Names of HapMap files produced by {@link TagsToSNPByAlignmentMTPlugin}.*/
    String[] hapmapFileName;

    /**Name of merged count file containing tags from all qseq files.*/
    String mergedCountFileName;

    /**Name of merged count file containing records from all TBT files.*/
    String mergedTBTFileName;

    /**Name of SAM alignment file.*/
    String SAMFileName;

    /**Name of .topm alignment file*/
    String TOPMFileName;

    /**All processed tags above {@link #minCount} in <b>all qseq files</b>.*/
    long totalTagsAboveMinCount;

    /**************************************************
     * Filter Options                                 *
     **************************************************/

    /**Names of HapMap files <b>after</b> filtering.  This should be in register with {@link #qseqFileName}.*/
    String[] filteredHapMapFileName;

    /**Minimum fraction of taxa that must contain a given tag.*/
    float minTaxonCoverage;

    /**Minimum fraction of all tags that must be found in any one taxon.*/
    float minSiteCoverage;

    /**Minimum value of the inbreeding coefficient (F statistic) measured for a given line.*/
    float minF;

    float minMAF, maxMAF,minLD,maxLD;

    /**Regular expression that limits analysis to taxon names that match it.*/
    String popMask;

    /**Names of taxa that matched the population mask of the filter step.*/
    String[] taxaMatchingPopulationMask;

    /**"Maximum error" value set in BiparentalErrorCorrection*/
    float maxError;

    /**"Minimum distortion" value set in BiparentalErrorCorrection*/
    float minDistortion;

    String currLine="";

    /**@param directoryName The full path name of the directory containing qseq files and other project files.*/
    Provenance(String directoryName){
        File directory = new File(directoryName);
        if(directory.isDirectory()){
            homeDirectory = directory;
            experimentName = directory.getName();
            provenanceFileName = experimentName+"_provenance.txt";

            try {
                BufferedReader br = new BufferedReader(new FileReader(
                    homeDirectory+
                    File.separator+
                    experimentName+
                    "_provenance.txt"
                ));
            }catch (FileNotFoundException e){
                System.out.println("WARNING: Couldn't find provenance file.");
                return;
            }
        }else{
            System.out.println("WARNING: The name of the directory containing the provenance file isn't a valid directory name.");
        }
    }

    /**Writes provenance to disk (in experiment directory) as a serialized object.*/
    public void save(){
        try{
            ObjectOutputStream os = new ObjectOutputStream(
                    new FileOutputStream(
                        homeDirectory+
                        File.separator+
                        experimentName+
                        "_provenance.dat"
                    )
            );
            os.writeObject(this);
            os.flush();
            os.close();
        } catch(Exception e){
            System.out.println("Caught exception while saving Provenance object to file: "+e);
        }
    }

    /**Loads provenance from disk (in experiment directory) as a serialized object.*/
    public static Provenance load(String directoryName){
        File directory = new File(directoryName);
        Provenance returnValue=null;
        if(directory.isDirectory()){

            String provenanceFileName =
                    directory.getAbsolutePath()+
                    File.separator+
                    directory.getName()+
                    "_provenance.dat";


        try{
            ObjectInputStream oi = new ObjectInputStream(new FileInputStream(provenanceFileName));
            returnValue = (Provenance)oi.readObject();
            oi.close();
//            BufferedReader br = new BufferedReader(new FileReader(homeDirectory+File.separator+provenanceFileName));
//            br.readLine(); //qseq category header
//            int currEntry = 0;
//            while(!currLine.equals("[taxa]\n")){
//                currLine = br.readLine();
//                if(currLine.equals("[taxa]\n")){break;}
////                String a=br.readLine();
////                String[] b=a.split("\t");
////                String c = b[1];
//                qseqFileName[currEntry] = currLine.split("\t")[1];
//                countFileName[currEntry] = br.readLine().split("\t")[1];
//                TBTFileName[currEntry] = br.readLine().split("\t")[1];
//                rawReads[currEntry] = Integer.parseInt(br.readLine().split("\t")[1]);
//                goodReads[currEntry] = Integer.parseInt(br.readLine().split("\t")[1]);
//                goodMatchedReads[currEntry] = Integer.parseInt(br.readLine().split("\t")[1]);
//                currEntry++;
//            }
//            br.readLine();//Single space
//            br.readLine();//taxa category header
//
//            currEntry = 0;
//            while(!br.readLine().equals("[hapmap files]\n")){
//                barcode[currEntry] = br.readLine().split("\t")[1];
//                taxonName[currEntry] = br.readLine().split("\t")[1];
//                currEntry++;
//            }
//            br.readLine();//Single space
//            br.readLine();//hapmap category header
//
//            currEntry = 0;
//            while(!br.readLine().equals("[dataset]\n")){
//                hapmapFileName[currEntry] = br.readLine().split("\t")[1];
//                filteredHapMapFileName[currEntry] = br.readLine().split("\t")[1];
//            }
//            br.readLine();//Single space
//            br.readLine();//Dataset category header
//            while(br.readLine()!=null){
//                 keyFileName= br.readLine().split("\t")[1];
//                  enzyme = br.readLine().split("\t")[1];
//               totalRawReads = Integer.parseInt(br.readLine().split("\t")[1]);
//               totalTagsAboveMinCount = Integer.parseInt(br.readLine().split("\t")[1]);
//               totalGoodReads = Integer.parseInt(br.readLine().split("\t")[1]);
//               totalGoodMatchedReads = Integer.parseInt(br.readLine().split("\t")[1]);
//               mergedCountFileName= br.readLine().split("\t")[1];
//               SAMFileName= br.readLine().split("\t")[1];
//               TOPMFileName= br.readLine().split("\t")[1];
//               mergedTBTFileName= br.readLine().split("\t")[1];
//               minCountPerExperiment= Integer.parseInt(br.readLine().split("\t")[1]);
//               minCountPerTaxon = Integer.parseInt(br.readLine().split("\t")[1]);
//               minTaxonCoverage = Float.parseFloat(br.readLine().split("\t")[1]);
//               minSiteCoverage = Float.parseFloat(br.readLine().split("\t")[1]);
//               minF = Float.parseFloat(br.readLine().split("\t")[1]);
//               minMAF = Float.parseFloat(br.readLine().split("\t")[1]);
//               maxMAF = Float.parseFloat(br.readLine().split("\t")[1]);
//               minLD = Float.parseFloat(br.readLine().split("\t")[1]);
//               maxLD = Float.parseFloat(br.readLine().split("\t")[1]);
//               maxError = Float.parseFloat(br.readLine().split("\t")[1]);
//               minDistortion = Float.parseFloat(br.readLine().split("\t")[1]);
//               popMask = br.readLine().split("\t")[1];
//            }

        }catch(Exception e){
            System.out.println("Caught exception while reading provenance file: "+e);
        }
        }else{
            System.out.println("WARNING: The name of the directory containing the provenance file isn't a valid directory name.");
        }
                    return returnValue;

    }

    /**Writes provenance to disk (in experiment directory) as a human-readable text file.*/
    public void write(){
        try {
            String s = homeDirectory.getAbsolutePath()+File.separator+provenanceFileName;
            BufferedWriter bw = new BufferedWriter(new FileWriter(s));

            for (int i= 0; i < qseqFileName.length; i++) {
                bw.write(
                    StringUtils.join(
                        new String[]{
                            "qseq file "+(i+1)+":",
                            "file name\t"+qseqFileName[i],
                            "tag count file\t"+countFileName[i],
                            "tbt file\t"+TBTFileName[i],
                            "raw reads\t"+Long.toString(rawReads[i]),
                            "processed tags\t"+Long.toString(goodReads[i]),
                            "tags matching taxa\t"+Long.toString(goodMatchedReads[i]),
                        }, "\n"
                    )+"\n\n"
                );
            }

            bw.write("[taxa]\nbarcode\ttaxon name\n");
            for (int i= 0; i < barcode.length; i++) {
                bw.write(barcode[i]+"\t"+taxonName[i]+"\n");
            }

            bw.write("\n[hapmap files]\n");
            for (int i= 0; i < hapmapFileName.length; i++) {
                bw.write(
                    StringUtils.join(
                            new String[]{
                                "raw file\t"+hapmapFileName[i],
                                "filtered file\t"+filteredHapMapFileName[i],
                            }, "\n"
                        )+"\n"
                );
            }

            bw.write("\n[dataset]\n");
            bw.write(
                    StringUtils.join(
                        new String[]{
                            "barcode key file\t"+keyFileName,
                            "enzyme\t"+enzyme,
                            "total raw reads\t"+totalRawReads+"",
                            "total processed tags above minimum count\t"+totalTagsAboveMinCount+"",
                            "total processed tags\t"+totalGoodReads+"",
                            "total tags matching taxa\t"+totalGoodMatchedReads+"",
                            "merged count file\t"+mergedCountFileName,
                            "sam file\t"+SAMFileName,
                            "topm file\t"+TOPMFileName,
                            "merged tbt file\t"+mergedTBTFileName,
                            "minimum unique tag count per experiment\t"+minCountPerExperiment+"",
                            "minimum unique tag count per taxon\t"+minCountPerTaxon+"",
                            "minimum taxon coverage\t"+minTaxonCoverage+"",
                            "minimum site coverage\t"+minSiteCoverage+"",
                            "minimum F\t"+minF+"",
                            "minimum MAF\t"+minMAF+"",
                            "maximum MAF\t"+maxMAF+"",
                            "minimum LD\t"+minLD+"",
                            "maximum LD\t"+maxLD+"",
                            "maximum error\t"+maxError+"",
                            "minimum distortion\t"+minDistortion+"",
                            "population mask\t"+popMask,
                        },"\n"
                    )+"\n"
            );
            bw.flush();
            bw.close();
        }catch (Exception e){
            System.out.println("WARNING: Caught exception while writing provenance file: "+e);
            return;
        }


    }
}
