package net.maizegenetics.gbs.pipeline;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.TreeSet;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByte;
import net.maizegenetics.gbs.tagdist.TagsByTaxaUtils;
import net.maizegenetics.genome.BaseEncoder;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.Locus;
import net.maizegenetics.pal.alignment.SimpleAlignment;
import net.maizegenetics.pal.datatype.IUPACNucleotides;
import net.maizegenetics.util.DirectoryCrawler;


/** @author jvh39 */
public class JVHPipeline {
    public static String[] sitesUniqueInDiscovery;
       
    public static void main(String[] args){
//        TagsOnPhysicalMap topm =  new TagsOnPhysicalMap("/media/jvh053111/all_sorghum_20120326/merged_032712.production.topm", true);
//        Alignment hapmap = ImportUtils.readFromHapmap("/media/jvh053111/all_sorghum_20120326/filt/all_sorghum_20120326_scv10mF8maf002.c1.hmp.txt");
//        TagsOnPhysicalMap.filter(topm, hapmap).writeTextFile(new File("/home/jvh39/Desktop/topm_filtered_040912.txt"));
        findProductionTOPM();
    }
    
    static int totalTaxa=0;
    static String[][] taxaNames;
    static DataInputStream[] tbts;
    static int[][] outputIndex;
    static boolean[] isOpen;

     public static void mergeTBTFiles() {
        String[] infiles = DirectoryCrawler.listFileNames(".*tbt.byte.*", "/media/SSD/mergetest");
        ArrayList<long[]> tags = new ArrayList<long[]>();
        TreeSet<String> allTaxa = new TreeSet<String>();
        int tagLengthInLong = 0;
        TagsOnPhysicalMap refTags = new TagsOnPhysicalMap("/media/jvh053111/all_zea_20120109/allZea_mappedonly_20120115.topm", true);
        tbts = new DataInputStream[infiles.length];
        int filesOpen=0;
        int[] numTags = new int[infiles.length];
        int[] tagsRead = new int[infiles.length];
        taxaNames = new String[infiles.length][];
        outputIndex=new int[infiles.length][];

        try {
            
            //Read headers
            for (int i = 0; i < tbts.length; i++) {
                tbts[i] = new DataInputStream(new BufferedInputStream(new FileInputStream(infiles[i]), 65536));

                numTags[i] = tbts[i].readInt();
                tagLengthInLong = tbts[i].readInt();
                int numTaxa = tbts[i].readInt();
                totalTaxa += numTaxa;
                
                taxaNames[i]=new String[numTaxa];
                for (int j = 0; j < numTaxa; j++) {
                    taxaNames[i][j]=tbts[i].readUTF();
                    allTaxa.add(taxaNames[i][j]);
                }
            }

            //Make lookup tables for input & output names
            for (int i = 0; i < taxaNames.length; i++) {
                String[] list = taxaNames[i];
                outputIndex[i]= MergeTagsByTaxaFilesPlugin.taxaRedirect(list, allTaxa.toArray(new String[allTaxa.size()]));
                for (int j = 0; j < list.length; j++) {
                }
            }

            filesOpen=infiles.length;
            byte[] outputDist = new byte[totalTaxa];
            byte[] currDist=null;
            long[] outputTag = null;
            long[][] inputTag = new long[infiles.length][];
            
            //Create array of #files X #taxa per file
            byte[][] inputDist = new byte[infiles.length][];
            for (int i = 0; i < taxaNames.length; i++) {
                inputDist[i]=new byte[taxaNames[i].length];
            }
            isOpen = new boolean[infiles.length];
            Arrays.fill(isOpen, true);
            
            for (int tag = 0; tag < refTags.tagNum; tag++){
                outputTag = refTags.getTag(tag);
//                System.out.println("Output tag is " + BaseEncoder.getSequenceFromLong(outputTag));
                //Loop through input files
                for (int i = 0; i < infiles.length; i++){

                    if (!isOpen[i]) continue; //Skip closed files
                    
                    //If queue is empty, load next tag.
                    if(inputTag[i]==null){
                        inputTag[i] = readTag(tbts[i]);
                        byte tagLength = tbts[i].readByte(); //Read through tag length
                        currDist = readDist(i);

                        System.out.println("Input tag "+i+" is: "+BaseEncoder.getSequenceFromLong(inputTag[i])+tagLength);
                        tagsRead[i]++;
                        if(tagsRead[i]==numTags[i]-1) isOpen[i]=false;  //Mark finished files
                    }
                    
                    //If input tag matches output tag, use its distribution and go to the next tag.
                    //If not, use a blank distribution and save the actual one for another round.
                    if (inputTag[i] == outputTag) {
                        System.out.println("Input tag "+i+" matches output tag.");
                        inputDist[i] = currDist;
                        inputTag[i] = null;
                    } else {
                        Arrays.fill(inputDist[i], (byte) 0);
                    }
               }
                //Now that all files have been touched, write out the merged distribution.
                for (int i = 0; i < inputDist.length; i++) {
                    for (int j = 0; j < inputDist[i].length; j++) {
                        //Increment output taxon that corresponds to input taxon
                        outputDist[outputIndex[i][j]] += inputDist[i][j]; 
                    }
                }
//                System.out.print(BaseEncoder.getSequenceFromLong(outputTag)+"\t");
//                System.out.print(refTags.getTagLength(tag) +"\t");
//                for (int i = 0; i < outputDist.length; i++) System.out.print(outputDist[i] +"\t");
//                System.out.println();
            }
            
        } catch (Exception e) {
            System.out.println("Caught exception while reading TBT files: " + e);
        }
    }
    
//    public static boolean allFilesOpen(){
//        for(boolean flag: isOpen){
//            if(flag==false) return false;
//        }
//        return true;
//    }
    public static long[] readTag(DataInputStream stream){
        long[] result=new long[2];

        try{
                result[0] = stream.readLong();
                result[1] = stream.readLong();
        } catch (Exception e) {
            System.out.println("Caught exception while reading tag: " + e);
        }
        return result;
    }
    
    /**Return the output distribution for the stream with the given index*/
    public static byte[] readDist(int file){
        byte[] result=new byte[taxaNames[file].length];
            
            for (int taxon = 0; taxon < taxaNames[file].length; taxon++) {
                //Loop over taxa in input file.  Read byte from file.
                //Put into output taxon corresponding to current input taxon
                try{
                        result[taxon] = tbts[file].readByte(); 
                } catch (Exception e) {
                    System.out.println("Caught exception while reading taxon distribution for file "+file+": " + e);
                }
            }
        return result;
    }
    
//                for (int j = 0; j < tagLengthInLong; j++) {
//                    currTag[j]=tbts[i].readLong();
//                }

//    public static void writeNetCDF(String filename){
//        TagsByTaxaByte tbt = new TagsByTaxaByte(filename, FilePacking.Byte);
//        int minCount=0;
//        int tagNum=tbt.getTagCount();
//        int taxaNum=tbt.getTaxaCount();
//        int tagLengthInLong=tbt.getTagSizeInLong();
//        String[] taxonNames=tbt.getTaxaNames();
//        int maxNameLength=0;
//        
//        for (String name: taxonNames) {
//            if(name.length() > maxNameLength) maxNameLength=name.length();
//        }
//
//        try {
//            //Open a new file: at this point it is in "define mode"
//           NetcdfFileWriteable output = NetcdfFileWriteable.createNew(filename+".nc");
//            
//           //Define dimensions.  There's an extra dimension for "taxon name character" because
//           //NetCDF can't store strings, so you have to store them as arrays of characters
//           Dimension tagAxis = output.addUnlimitedDimension("tag_number");
//           Dimension taxonAxis = output.addUnlimitedDimension("taxon_name");
//           Dimension taxonNameCharAxis=output.addUnlimitedDimension("taxon_name_char"); 
//           Dimension tagLengthInLongAxis = output.addDimension("tag_length_in_long", 2);
//           
//           //Define types
//           output.addVariable("taxon", DataType.STRING, new Dimension[]{taxonAxis, taxonNameCharAxis});
//           output.addVariable("tag", DataType.LONG, new Dimension[]{tagAxis, tagLengthInLongAxis});
//           output.addVariable("count", DataType.BYTE, new Dimension[]{tagAxis, taxonAxis});
//           
//           output.create(); //This ends "define mode" so I can start writing.
//           
//           //Prepare arrays for output
//           ArrayChar.D2 taxonArray = new ArrayChar.D2(maxNameLength, taxaNum);
//           ArrayByte.D2 countArray= new ArrayByte.D2(taxaNum, tagNum);
//           ArrayLong.D2 tagArray = new ArrayLong.D2(tagNum, tagLengthInLong);
//           
//           //Load arrays with data from TBT file.
//            for (int i = 0; i < taxaNum; i++) {
//                
//                //Store taxon name
//                char[] taxonName=taxonNames[i].toCharArray();
//                for (int j = 0; j < maxNameLength; j++) {
//                    taxonArray.set(j,i, taxonName[j]);
//                }
//                
//                for (int j = 0; j < tagNum; j++) {
//                    //Store tag sequences on first iteration of loop
//                    if(j==0){
//                        for (int k = 0; k < tagLengthInLong; k++) {
//                            tagArray.set(i,k, tbt.tags[k][i]);
//                        }
//                    }
//                    
//                    //Store tag distribution on every iteration of loop
//                    countArray.set(i,j, tbt.tagDist[i][j]);
//                }
//            }
//            
//            output.write("taxon", taxonArray);
//            output.write("tag", tagArray);
//            output.write("count", countArray);
//            
//        } catch (Exception e) {
//            System.out.println("Caught exception while writing netCDF file:"+e);
//        }
//        
//    }
//    
//    
//    public static void readHDF5IJV(String filename){
//        int minCount=0;
//        
//        IHDF5Reader reader = HDF5Factory.open(filename+".h5");
//        
//        String[] taxonNames=reader.readStringArray("taxon_names");
//        long[][] tagSequences=reader.readLongMatrix("tag_sequences");
//        int[] tagIndices= reader.readIntArray("tag_index");
//        int[] taxonIndices= reader.readIntArray("taxon_index");
//        byte[] counts= reader.readByteArray("count");
//        reader.close();
//        
//        
//        TagsByTaxaByte tbt = new TagsByTaxaByte(filename, FilePacking.Byte);
//        tbt.initMatrices(taxonNames.length, tagSequences[0].length);
//        
//        tbt.taxaNames=taxonNames;
//        tbt.tags=tagSequences;
//        
//        for (int i = 0; i < tagIndices.length; i++) {
//            tbt.setReadCountForTagTaxon(tagIndices[i], taxonIndices[i], counts[i]);
//        }
//        
//        tbt.writeDistFile(new File(filename+".output"), FilePacking.Byte, 0);
//    }
//    
//    public static void writeHDF5IJV(String filename){
//        TagsByTaxaByte tbt = new TagsByTaxaByte(filename, FilePacking.Byte);
//        int minCount=0;
//        int tagNum=tbt.getTagCount();
//        int taxaNum=tbt.getTaxaCount();
//        
//        System.out.println("Measuring matrix sparsity.");
//        int records=0;
//        for (int i = 0; i < tagNum; i++) {
//            if(i%1000000==0) System.out.println("Read "+i+" tags.");
//            for (int j = 0; j < taxaNum; j++) {
//                if( tbt.tagDist[j][i] > minCount) records++;
//            }
//        }
//        
//        int[] tagIndices= new int[records];
//        int[] taxonIndices= new int[records];
//        byte[] counts=new byte[records];
//        
//        System.out.println("Consolidating non-zero records.");
//        int currRecord=0;
//        for (int i = 0; i < tagNum; i++) {
//            if(i%1000000==0) System.out.println("Read "+i+" tags.");
//            for (int j = 0; j < taxaNum; j++) {
//                if( tbt.tagDist[j][i] > minCount){
//                    tagIndices[currRecord]=i;
//                    taxonIndices[currRecord]=j;
//                    counts[currRecord]=tbt.tagDist[j][i];
//                    currRecord++;
//                }
//            }
//        }
//
//        IHDF5Writer writer = HDF5Factory.open(filename+".h5");
//        
//        writer.writeStringArray("taxon_names", tbt.getTaxaNames());
//        writer.writeLongMatrix("tag_sequences", tbt.tags);
//        writer.writeIntArray("tag_index", tagIndices);
//        writer.writeIntArray("taxon_index", taxonIndices);
//        writer.writeByteArray("count", counts);
//        writer.close();
//    }
    
    public static void findProductionTOPM(){
        File[] files=DirectoryCrawler.listFiles(".*topm","/media/jvh053111");
        for(File file: files){
            System.out.print("File "+file.getName()+":  ");
            try {
                DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(file), 65536));
                
                //Read header
                int tagNum = dis.readInt();
                int tagLengthInLong = dis.readInt();
                int maxVariants = dis.readInt();
                
                //Read body
                boolean variantsFound=false;
                for (int i = 0; i < 1000; i++) {

                    long[] tags=new long[tagLengthInLong];
                    for (int j = 0; j < tagLengthInLong; j++) {
                        tags[j] = dis.readLong();
                    }
                    int tagLength = dis.readByte();
                    int multimaps = dis.readByte();
                    int chromosome = dis.readInt();
                    int strand = dis.readByte();
                    int startPosition = dis.readInt();
                    int endPosition = dis.readInt();
                    int divergence = dis.readByte();

                    byte[] variantPosOff=new byte[maxVariants];
                    byte[] variantDef=new byte[maxVariants];
                    for (int j = 0; j < maxVariants; j++) {
                        variantPosOff[j] = dis.readByte();
                        variantDef[j] = dis.readByte();
                        if(variantPosOff[j] != Byte.MIN_VALUE){
                            variantsFound=true;
                            break;
                        }
                    }
                    byte dcoP = dis.readByte();
                    byte mapP = dis.readByte();

                    if(variantsFound==true){
                        System.out.print("Variants found.");
                        break;
                    }
                }
                System.out.println();
                dis.close();
            } catch (Exception e) {
                System.out.println("Caught exception while counting tags in TOPM file.");
            }
        }
    }

    public static void testFastqToTagCountPlugin(){

        String[] args={
            "-i", "/mnt/nextgen/GBSFlowCells/81546ABXX/fastq",
            "-k", "/home/jvh39/Desktop/test.key",
            "-e", "ApekI",
            "-s","50000000",
            "-o", "/home/jvh39/Desktop",
        };

        FastqToTagCountPlugin testPlugin = new FastqToTagCountPlugin();
        testPlugin.setParameters(args);
        testPlugin.performFunction(null);
    }

    public static void testRawReadsToHapmapPlugin(){

        String[] args={
            "-i", "/mnt/nextgen/GBSFlowCells/D0MY9ACXX",
            "-k", "/media/jvh053111/yunbi_lines_032012/yunbi_032012.key",
            "-e", "ApekI",
            "-m","/mnt/nextgen/Zea/build20120110/topm/zea20120110c510.prod510.topm",
            "-o", "/media/jvh053111/yunbi_lines_032012",
        };

        RawReadsToHapMapPlugin testPlugin = new RawReadsToHapMapPlugin();
        testPlugin.setParameters(args);
        testPlugin.performFunction(null);
    }

    public static void measureTOPMFiles(){
        
        File[] files=DirectoryCrawler.listFiles(".*topm.bin.*","/media");
        for(File file: files){
            System.out.println("File "+file.getName());
            try {
                DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(file), 65536));
                int count=dis.readInt();
                System.out.println(count+" tags in "+file.getName()+".");
                dis.close();
            } catch (Exception e) {
                System.out.println("Caught exception while counting tags in TOPM file.");
            }
        }
    }
   
    public static void report(String directory){
        String[] hapmapFileNames = DirectoryCrawler.listFileNames(".*hmp.txt", directory);///mnt/nextgen/Zea/hapmap/build112311/unfilt");
        for (String name : hapmapFileNames) {
            
            //Use object factory in SimpleAlignment to get an object that implements Alignment
            Alignment a=ImportUtils.readFromHapmap(name);
            SimpleAlignment file=SimpleAlignment.getInstance(a);
            a=null;
            System.gc();
            
            int taxonCount=file.getIdGroup().getIdCount();
            int siteCount=file.getSiteCount();
            
            //Taxon properties
            int[] taxonCalls=new int[taxonCount];
            int[] taxonHetCalls=new int[taxonCount];
            
            //Histogram properties
            int numBins=20;
            int binWidth=100/numBins;
            int[] bins=new int[numBins];
            int[] hetBins=new int[numBins];
            
            System.out.println(
                name+"\n"+
                "Total taxa:\t"+file.getIdGroup().getIdCount()+"\n"+
                "Total SNPs:\t"+file.getSiteCount()+"\n"
            );
            
            System.out.println(
                    "Site\t"
                    + "Calls\t"
                    + "Het Calls\t"
                    + "Call Rate\t"
                    + "Het Call Rate\t"
                    +"MAF"
            );
                    
            System.out.println("Site Properties");
            for (int site = 0; site < siteCount; site++){
                int calls=0, hetCalls=0;
                //Skip uncalled bases.  Increment "calls" for all called bases and "hetCalls" for hets.
                for (int taxon = 0; taxon < taxonCount; taxon++) {
                    char base = (char)file.getBase(taxon, site);
                    switch(base){
                        case 'N': 
                            break;
                        case 'A':  case 'C':  case 'G':  case 'T':  
                            calls++; 
                            taxonCalls[taxon]++;
                            break;
                        default: 
                            calls++; 
                            hetCalls++; 
                            taxonCalls[taxon]++;
                            taxonHetCalls[taxon]++;
                            break;
                    }
                }
                
                //Calculate call rates (site coverage) as #calls/#taxa
                float callRate = ((float)calls/(float)taxonCount);
                float hetCallRate = ((float)hetCalls/(float)taxonCount);

                //Add SNP to an appropriate chart bin based on its call rate
                double percentile=callRate*100;
                for (int bin = 0; bin < bins.length; bin++) {
                    int upperBound = ((bin+1)*binWidth);
                    int lowerBound = (bin*binWidth);
                    if(percentile < upperBound && percentile > lowerBound)  bins[bin]+=calls;
                }
                
                System.out.println(
                    site+"\t"+
                    calls+"\t"+
                    hetCalls+"\t"+
                    callRate+"\t"+
                    hetCallRate+"\t"+
                    file.getMinorAlleleFrequency(site)
                );
            }
            System.out.println();
            
            
            //Print bar chart data for sites
            System.out.println("Site Coverage:");
            System.out.println("Call Rate(%)\t"+"SNPs");
            for (int bin = 0; bin < bins.length; bin++) {
                int callRatePct =(bin*binWidth);
                System.out.println(callRatePct+"\t"+bins[bin]);
                bins[bin]=0;
            }
            
            //Print taxon coverage stats
            System.out.println("Taxon Properties");
            System.out.println(
                    "Name\t"
                    +"Calls\t"
                    +"Het Calls\t"
                    +"Call Rate\t"
                    +"Het Call Rate\t"
            );
            
            for (int i = 0; i < taxonCount; i++) {
                double callRate=taxonCalls[i]/siteCount;
                double hetCallRate=taxonHetCalls[i]/siteCount;
                System.out.println(
                        file.getTaxaName(i)+"\t"
                        +taxonCalls[i]+"\t"
                        +taxonHetCalls[i]+"\t"
                        +callRate+"\t"
                        +hetCallRate+"\t"
                );
                
                //Add SNP to an appropriate chart bin based on its call rate
                double percentile=callRate*100;
                for (int bin = 0; bin < bins.length; bin++) {
                    int upperBound = ((bin+1)*binWidth);
                    int lowerBound = (bin*binWidth);
                    if(percentile < upperBound && percentile > lowerBound)  bins[bin]+=taxonCalls[i];
                }
            }
            
            //Print bar chart data for taxa
            System.out.println("Taxon Coverage:");
            System.out.println("Call Rate(%)\t"+"Taxa");
            for (int bin = 0; bin < bins.length; bin++){
                int callRatePct =(bin*binWidth);
                System.out.println(callRatePct+"\t"+bins[bin]);
            }
            System.gc();
        }
    }
    
    public static void testMergeTagsByTaxaFilesPlugin(){
        String[] args = new String[]{
            "-i","/media/jvh053111/all_zea_20120109/tbt/test",
            "-o","/media/SSD/merged.tbt.byte",
        };
        MergeTagsByTaxaFilesPlugin testPlugin=new MergeTagsByTaxaFilesPlugin();
        testPlugin.setParameters(args);
        testPlugin.performFunction(null);
    }
    
    public static void testSNPMergePlugin(){
        MergeDuplicateSNPsPlugin testPlugin = new MergeDuplicateSNPsPlugin();
        testPlugin.setParameters(new String[]{
            "-hmp","/media/jvh053111/qseq_to_hapmap_test/10247896_B08G4ABXX_s_1.hmp.txt",
            "-o","/media/jvh053111/qseq_to_hapmap_test/merged.c+.hmp.txt",
            "-s","1",
            "-callHets",
            "-misMat","0.2",
//            "-kpUnmergDups",
            "-e","10"
        });
        testPlugin.performFunction(null);
    }
    public static void testTOPMWithPhysicalPositions(){
        QseqToHapMapPlugin testPlugin = new QseqToHapMapPlugin();
//        testPlugin.setParameters(new String[]{
//            "-i","/media/Data/Paola/test",
//            "-k","/media/Data/Paola/key_gbs_grape.tsv",
//            "-e","apeki",
//            "-o","/media/Data/Paola/test",
//            "-m","/home/jvh39/topm_position_test.topm"
//        });
   
        testPlugin.setParameters(new String[]{
            "-i","/media/jvh053111/Ames",
            "-k","/media/NAM/all_maize_20110422.key",
            "-e","apeki",
            "-o","/media/jvh053111/qseq_to_hapmap_test/110411",
            "-m","/media/jvh053111/qseq_to_hapmap_test/mergedNAM282Ames_072011.prod.topm.bin",
            "-d","1"
//            "-m","/media/jvh053111/qseq_to_hapmap_test/chr10.topm.txt"
        });
   
       testPlugin.performFunction(null);
    }
    
        public static void runMergeDuplicateSNPsPlugin() {

           for(String dirName: new String[]{"6","7","8","9"}){
            
           String workDir = "/media/SSD/mnF0."+dirName+"/", 
                   inputName="SSD",
                   suffix=".c+.hmp.txt",
                   outputName=inputName+"_mergedSNPs";

           String[] args = new String[] {
                "-hmp", workDir+inputName+suffix, // Input HapMap file; use a plus sign (+) as a wild card character to specify multiple chromosome numbers
                "-o",   workDir+outputName+suffix, // Output HapMap file
                "-misMat", "0.15", // Threshold genotypic mismatch rate above which the duplicate SNPs won't be merged (default: 0.05)
                "-callHets",       // When two genotypes disagree at a SNP, call it a heterozygote (default: off = set to missing)
    //          "-kpUnmergDups",   // When two duplicate SNPs were not merged (different alleles or too many mismatches), keep them (default: off = delete them)
                "-s", "10",        // Start chromosome (default 1)
                "-e", "10"         // End chromosome (default 10)
            };

            MergeDuplicateSNPsPlugin plugin = new MergeDuplicateSNPsPlugin();
            plugin.setParameters(args);
            plugin.performFunction(null);
    }
        }
     public static void filterMergeImpute10KMaizeTaxa() {
        String rootIn="/media/jvh053111/all_zea_20120109/hapmap/SSD";
        String rootOut="/media/jvh053111/all_zea_20120109/hapmap/allZea20120110";
        int sC=3;
        int eC=3;
        String[] args = new String[] {
            "-hmp", rootIn+".c+.hmp.txt",
            "-o", rootOut+"_scv10mF8maf002.c+.hmp.txt",
            "-mnTCov", "0.01",
            "-mnF",".8",
            "-mnMAF","0.002",
//            "-hLD",
            "-mnSCov", "0.10", // Minimum locus coverage (proportion of Taxa)
            "-sC", ""+sC,  // Start chromosome
            "-eC", ""+eC // End chromosome
        };

//        GBSHapMapFiltersPlugin testClass = new GBSHapMapFiltersPlugin();
//        testClass.setParameters(args);
//        testClass.performFunction(null);
        
        args = new String[] {
            "-hmp", rootOut+"_scv10mF8maf002.c+.hmp.txt", // Input HapMap file; use a plus sign (+) as a wild card character to specify multiple chromosome numbers
            "-o", rootOut+"_scv10mF8maf002_mgs.c+.hmp.txt", // Output HapMap file
            "-misMat", "0.15", // Threshold genotypic mismatch rate above which the duplicate SNPs won't be merged (default: 0.05)
            "-callHets",       // When two genotypes disagree at a SNP, call it a heterozygote (default: off = set to missing)
//          "-kpUnmergDups",   // When two duplicate SNPs were not merged (different alleles or too many mismatches), keep them (default: off = delete them)
            "-s", ""+sC,        // Start chromosome (default 1)
            "-e", ""+eC         // End chromosome (default 10)
        };

//        MergeDuplicateSNPsPlugin plugin = new MergeDuplicateSNPsPlugin();
//        plugin.setParameters(args);
//        plugin.performFunction(null);
            
        args = new String[] {
            "-hmp", rootOut+"_scv10mF8maf002_mgs.c+.hmp.txt",
            "-o", rootOut+"_scv10mF8maf002_mgs_E1pLD5kpUn.c+.hmp.txt",
            "-oB", "/Users/edbuckler/SolexaAnal/GBS/test/errorBin.txt",
            "-oE", "/Users/edbuckler/SolexaAnal/GBS/test/errorBySNP.txt",
            "-popM","Z[0-9]{3}",
            "-sC",""+sC,
            "-eC",""+eC,
            "-mxE","0.01",
            "-mnD","2.0",
            "-mnPLD", "0.5",
            "-kpUT"
            };
        BiParentalErrorCorrectionPlugin bpec = new BiParentalErrorCorrectionPlugin();
        bpec.setParameters(args);
        bpec.performFunction(null);
        args = new String[] {
            "-hmp", rootOut+"_scv10mF8maf002_mgs_E1pLD5kpUn.c+.hmp.txt",
            "-o", rootOut+"_scv10mF8maf002_mgs_E1pLD5kpUn_mgNoHet.c+.hmp.txt",
            "-xHets",
            "-hetFreq", "0.76",
            "-sC", ""+sC,  // Start chromosome
            "-eC", ""+eC // End chromosome
        };
        MergeIdenticalTaxaPlugin mitp = new MergeIdenticalTaxaPlugin();
        mitp.setParameters(args);
        mitp.performFunction(null);

        args = new String[] {
            "-hmp", rootOut+"_scv10mF8maf002_mgs_E1pLD5kpUn_mgNoHet.c+.hmp.txt",
            "-o", rootOut+"_scv10mF8maf002_mgs_E1pLD5kpUn_mgNoHet_imp.c+.hmp.txt",
            "-sC", ""+sC,  // Start chromosome
            "-eC", ""+eC // End chromosome
        };
        FastImputationBitFixedWindow.main(args);
    }

     public static void imputeAllZea() {
        int sC=7;
        int eC=7;
        String[] args = new String[] {
            "-hmp", "/mnt/nextgen/Zea/hapmap/build111123/mergedTaxa/allZea_112311_SNPmerge15_cov10_fT1E1pLD_mergedTaxa_c+.hmp.txt",
            "-o", "/media/jvh053111/all_zea_tagsbytaxa_110411/imputed/allZea_112311_SNPmerge15_cov10_fT1E1pLD_mergedTaxa_imputed_c+.hmp.txt",
            "-sC", ""+sC,  // Start chromosome
            "-eC", ""+eC // End chromosome
        };
        FastImputationBitFixedWindow.main(args);
    }

    
    public static void createTOPMWithPhysicalPositions(){
        TagsOnPhysicalMap exampleTOPM=new TagsOnPhysicalMap("/media/Data/Paola/paola_grape_tags_072711.topm.bin", true);
        for (int i = 0; i < exampleTOPM.getTagCount(); i++) {
            int currChr=exampleTOPM.getChromosome(i);
                if(currChr==Integer.MIN_VALUE) continue;
            byte currOffset = (byte)(Math.random()*64);
            for (int j = 0; j < 4; j++) {
                byte base = randomBase();
                exampleTOPM.setVariantDef(i,j, base);
                exampleTOPM.setVariantPosOff(i,j, currOffset);
//                System.out.println("Def: "+exampleTOPM.getVariantDef(i, j));
//                System.out.println("Pos: "+exampleTOPM.getVariantPosOff(i, j));
            }
        }
        exampleTOPM.writeBinaryFile(new File("/home/jvh39/topm_position_test.topm"));
    }
    
       public static void runTagsToSNPByAlignmentPlugin() {
        String[] args = new String[] {
            "-i", "/media/SSD/all_Zea_high_coverage_121511.tbt.byte",
//            "-i", "/media/jvh053111/all_zea_tagsbytaxa_110411/chunk1/434GFAAXX_s_7.tbt.byte",
            "-o", "/media/SSD",
//            "-o","/media/Data",
            "-m", "/media/Data/mergedNam282Ames_072011_mappedOnly.topm.bin",
            "-y",
            "-mnF",   "0.8",
            "-mnMAF", "0.2",
            "-mnMAC", "99999",
            "-mnLCov","0.1", // Minimum locus coverage (proportion of Taxa)
//            "-callBiSNPsWGap",
//            "-inclGaps",  // Include sites where major or minor allele is a GAP
            "-s", "1",  // Start chromosome
            "-e", "10"  // End chromosome
        };
       
  
        TagsToSNPByAlignmentPlugin testClass = new TagsToSNPByAlignmentPlugin();
        testClass.setParameters(args);
        testClass.performFunction(null);
}
       
       public static void runTagsToSNPByAlignmentMTPlugin() {
        String[] args = new String[] {
            "-i", "/media/SSD/allZea_111611.tbt.filteredbytopm.byte",
//            "-i", "/media/Data/new_snp_caller/434GFAAXX_s_7.tbt.bin",
//            "-o", "/media/SSD",
            "-o","/media/Data/new_snp_caller/compare_old_snp_caller",
            "-m", "/media/Data/mergedNam282Ames_072011_mappedOnly.topm.bin",
            "-y",
            "-mnF",   "0.9",
            "-mnMAF", "0.005",
            "-mnMAC", "20",
            "-mnLCov","0.10", // Minimum locus coverage (proportion of Taxa)
//            "-inclGaps",  // Include sites where major or minor allele is a GAP
            "-s", "10",  // Start chromosome
            "-e", "10"  // End chromosome
        };

  
        TagsToSNPByAlignmentMTPlugin testClass = new TagsToSNPByAlignmentMTPlugin();
        testClass.setParameters(args);
        testClass.performFunction(null);
    }
       
    public static byte randomBase(){
        byte[] bases = new byte[]{(byte)'A',(byte)'C',(byte)'G',(byte)'T'};
        return bases[(int)(Math.random()*bases.length)];
    }

    public static void searchVariants(TagsOnPhysicalMap theTOPM){
                                 
         for (int tag = 0; tag < theTOPM.getTagCount(); tag++){
             int i = theTOPM.getReadIndexForPositionIndex(tag);
             int startPos=theTOPM.getStartPosition(i);

             for (int k = 0; k < theTOPM.maxVariants; k++) {
                int offset=theTOPM.getVariantPosOff(i, k);
                int position = startPos+offset;
                if(searchArray(sitesUniqueInDiscovery, position)){
                    System.out.println("at tag "+i+", variant "+k+".");
//                    theTOPM.printRow(i);
                    System.out.println(BaseEncoder.getSequenceFromLong(theTOPM.getTag(i)));
                }
            }
         }
    }
    
    public static boolean searchArray(String[] refPos, int queryPos){
        for (int i = 0; i < refPos.length; i++) {
            if(Integer.parseInt(refPos[i])==queryPos){
              System.out.println("Found discovery pipeline site: "+queryPos);
                return true;
            }
        }
        return false;
    }
    
    /**Read a file into a string array.*/
    public static String[] slurp(String reference) {
        ArrayList<String> result=new ArrayList<String>();
        String currLine;
        try {
            BufferedReader br = new BufferedReader(new FileReader(reference));
            while ((currLine = br.readLine()) != null) {
                result.add(currLine);
            }
        } catch (Exception e) {
            System.out.println("File not found.");
            return null;
        }
        return result.toArray(new String[result.size()]);
    }
     
        public static void getHapMapReport() {
            for (int fLevel = 6; fLevel < 10; fLevel++) {
                
        String infileName = "/media/SSD/mnF0."+fLevel+"/SSD_mergedSNPs_filt_mergedtaxa.c10.hmp.txt";
//            String infileName="/media/Data/mnF_benchmark/test.hmp.txt";
        //Use object factory in SimpleAlignment to get an object that implements Alignment
        Alignment a=ImportUtils.readFromHapmap(infileName);
        SimpleAlignment file=SimpleAlignment.getInstance(a);
        a=null;
        System.gc();

        int taxonCount=file.getIdGroup().getIdCount(), nonBlankTaxaCount=0;
        for (int i = 0; i < taxonCount; i++){
            String name=file.getIdGroup().getIdentifier(i).getName();
                if(!name.equalsIgnoreCase("blank")) nonBlankTaxaCount++; 
        }

        int siteCount=file.getSiteCount();

        //Taxon properties
        int[] taxonCalls=new int[taxonCount];
        int[] taxonHetCalls=new int[taxonCount];

        //Histogram properties
        int numBins=20;
        int binWidth=100/numBins;
        int[] bins=new int[numBins];
        int[] hetBins=new int[numBins];

        System.out.println(
            infileName+"\n"+
            "Total taxa:\t"+file.getIdGroup().getIdCount()+"\n"+
            "Total non-\"blank\" taxa:\t"+nonBlankTaxaCount+"\n"+
            "Total SNPs:\t"+file.getSiteCount()+"\n"
        );

//        System.out.println(
//                "Site\t"
//                + "Calls\t"
//                + "Het Calls\t"
//                + "Call Rate\t"
//                + "Het Call Rate\t"
//                +"MAF"
//        );

        System.out.println("Site Properties");
        for (int site = 0; site < siteCount; site++){
            int calls=0, hetCalls=0;
            //Skip uncalled bases.  Increment "calls" for all called bases and "hetCalls" for hets.
            for (int taxon = 0; taxon < taxonCount; taxon++) {
                if(file.getIdGroup().getIdentifier(taxon).getName().equalsIgnoreCase("blank")) continue;    //Skip blank samples

                char base = (char)file.getBase(taxon, site);
                switch(base){
                    case 'N':
                        break;
                    case 'A':  case 'C':  case 'G':  case 'T':
                        calls++;
                        taxonCalls[taxon]++;
                        break;
                    default:
                        calls++;
                        hetCalls++;
                        taxonCalls[taxon]++;
                        taxonHetCalls[taxon]++;
                        break;
                }
            }

            //Calculate call rates (site coverage) as #calls/#taxa
            float callRate = ((float)calls/(float)nonBlankTaxaCount);
            float hetCallRate = ((float)hetCalls/(float)nonBlankTaxaCount);

            //Add SNP to an appropriate chart bin based on its call rate
            double percentile=callRate*100;
            for (int bin = 0; bin < bins.length; bin++) {
                int upperBound = ((bin+1)*binWidth);
                int lowerBound = (bin*binWidth);
                if(percentile <= upperBound && percentile > lowerBound)  bins[bin]+=1;
            }

//            System.out.println(
//                site+"\t"+
//                calls+"\t"+
//                hetCalls+"\t"+
//                callRate+"\t"+
//                hetCallRate+"\t"+
//                file.getMinorAlleleFrequency(site)
//            );
        }
        System.out.println();


        //Print bar chart data for sites
        System.out.println("Site Coverage:");
        System.out.println("%Taxa containing site:\t"+"SNPs");
        for (int bin = 0; bin < bins.length; bin++) {
            int callRatePct =(bin*binWidth);
            System.out.println(callRatePct+"\t"+bins[bin]);
            bins[bin]=0; //Re-zero bins after they are printed
        }

        //Print taxon coverage stats
        System.out.println("Taxon Properties");
        System.out.println(
                "Name\t"
                +"Calls\t"
                +"Het Calls\t"
                +"Call Rate\t"
                +"Het Call Rate\t"
        );

        for (int i = 0; i < taxonCount; i++) {
            if(file.getIdGroup().getIdentifier(i).getName().equalsIgnoreCase("blank")) continue;    //Skip blank samples
            double callRate=(float)taxonCalls[i]/(float)siteCount;
            double hetCallRate=(float)taxonHetCalls[i]/(float)siteCount;
//            System.out.println(
//                    file.getTaxaName(i)+"\t"
//                    +taxonCalls[i]+"\t"
//                    +taxonHetCalls[i]+"\t"
//                    +callRate+"\t"
//                    +hetCallRate+"\t"
//            );

            //Add SNP to an appropriate chart bin based on its call rate
            double percentile=callRate*100;
            for (int bin = 0; bin < bins.length; bin++) {
                int upperBound = ((bin+1)*binWidth);
                int lowerBound = (bin*binWidth);
                if(percentile <= upperBound && percentile > lowerBound)  bins[bin]+=1;
            }
        }

        //Print bar chart data for taxa
        System.out.println("Taxon Coverage:");
        System.out.println("%Sites found in taxon:\t"+"Taxa");
        for (int bin = 0; bin < bins.length; bin++){
            int callRatePct =(bin*binWidth);
            System.out.println(callRatePct+"\t"+bins[bin]);
        }
        System.gc();
    }
                    }

    }

