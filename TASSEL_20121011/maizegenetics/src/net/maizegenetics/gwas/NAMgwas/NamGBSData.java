
package net.maizegenetics.gwas.NAMgwas;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Pattern;

import net.maizegenetics.gwas.Imputation.GbsGenotypeData;
import net.maizegenetics.gwas.Imputation.ProbabilityObservationGivenTrueState;
import net.maizegenetics.gwas.Imputation.Recombination;
import net.maizegenetics.gwas.Imputation.SnpInformation;
import net.maizegenetics.gwas.Imputation.TestTransitionMatrix;
import net.maizegenetics.gwas.Imputation.TransitionMatrix;
import net.maizegenetics.gwas.Imputation.ViterbiAlgorithm;
import net.maizegenetics.gwas.jointlinkage.SNP;


import cern.jet.stat.Probability;


public class NamGBSData {
	
	//data pipeline for translating gbs data in parental allele calls
	//step 1. create a NamGBSData instance for an input file. This reads snps into a byte array then imputes AB calls using imputeABH().
	//step 2. impute hets
	//step 3. create a SNP data set
	//step 4. concatenate SNP data sets
	
	int chr = 0;
	ArrayList<String> taxanames;
	ArrayList<String[]> snplabels;
	String[][] snpLabels;
	final String[] headerForLabels = new String[11];
	int[] snpPositions;
	byte[][] snps;
	byte[][] imputedSnps;
	byte[] b73Allele;
	byte[] nonb73Allele;
	int numberOfTaxa, numberOfSnps, numberOfPopulations;
	String inputFilename;
	String hdfFilename = "";
	BufferedWriter hapout;
	boolean writeAsHdf = true;
	
	final HashMap<String, Byte> snpToByteMap = new HashMap<String, Byte>();;
	HashMap<String, ArrayList<Integer>> popmap; //which taxa are in each population
	HashMap<String, boolean[]> segregatingSites;
	
	static String[] byteToSnp = new String[]{"A","C","G","T","R","Y","S","W","K","M","+","0","H","B","N","-","X"};

	static final byte A = (byte) 0;
	static final byte C = (byte) 1;
	static final byte G =(byte) 2;
	static final byte T =(byte) 3;
	static final byte R =(byte) 4;		//A/G
	static final byte Y =(byte) 5; 	//C/T
	static final byte S =(byte) 6; 	//G/C
	static final byte W =(byte) 7;		//A/T
	static final byte K =(byte) 8;		//G/T
	static final byte M =(byte) 9;		//A/C	
	static final byte plus =(byte) 10;	
	static final byte zero =(byte) 11;	//+/-
	static final byte H =(byte) 12;	
	static final byte B =(byte) 13;	
	static final byte N =(byte) 14;	
	static final byte minus =(byte) 15;
	static final byte X = (byte) 16;
	
	final byte[][] hetcodes = new byte[][]{{A,M,R,W},{M,C,S,Y},{R,S,G,K},{W,Y,K,T}}; 

	private void loadSnpToByte() {
		snpToByteMap.put("A", A);
		snpToByteMap.put("C", C);
		snpToByteMap.put("G", G);
		snpToByteMap.put("T", T);
		snpToByteMap.put("R", R);
		snpToByteMap.put("Y", Y);
		snpToByteMap.put("S", S);
		snpToByteMap.put("W", W);
		snpToByteMap.put("K", K);
		snpToByteMap.put("M", M);
		snpToByteMap.put("+", plus);
		snpToByteMap.put("0", zero);
		snpToByteMap.put("H", H);
		snpToByteMap.put("B", B);
		snpToByteMap.put("N", N);
		snpToByteMap.put("-", minus);
	}
	
	public NamGBSData() {
		loadSnpToByte();
	}

	public void imputeParents(String filename) {
		popmap = new HashMap<String, ArrayList<Integer>>();
		segregatingSites = new HashMap<String, boolean[]>();
		taxanames = new ArrayList<String>();
		inputFilename = filename;
		ArrayList<Integer> columns = new ArrayList<Integer>();
		String[] info;
		String input;
		Pattern tab = Pattern.compile("\t");
		Pattern slash = Pattern.compile("/");
		int count;
		File inputFile = new File(filename);
		File dir = inputFile.getParentFile();
		System.out.println("Reading the snp file...");
		try {
			BufferedReader br = new BufferedReader(new FileReader(filename));
			
			//count the number of lines in the file and get the chromosome number
			count = 1;
			br.readLine(); //the header
			info = tab.split(br.readLine());
			chr = Integer.parseInt(info[2]);
			while (br.readLine() != null) count++;
			br.close();
			numberOfSnps = count;
			
			//open the file again and read the taxanames. 
			br = new BufferedReader(new FileReader(filename));
			info = tab.split(br.readLine());
			
			for (int i = 0; i < 11; i++) headerForLabels[i] = new String(info[i]);
			count = 0;
			for (String str : info) {
				if (str.startsWith("Z0")) {
					taxanames.add(new String(str));
					columns.add(count);
				}
				count++;
			}
			numberOfTaxa = taxanames.size();
			
			//parse the population number from the taxa names, create a list of taxa column numbers for each population
			getpops();
			
			//load snp calls into a byte array
			snps = new byte[numberOfTaxa][numberOfSnps];
			b73Allele = new byte[numberOfSnps];
			nonb73Allele = new byte[numberOfSnps];
			snpPositions = new int[numberOfSnps];
//			snplabels = new ArrayList<String[]>(numberOfSnps);
			snpLabels = new String[numberOfSnps][];

			int incount = 0;
			while ((input = br.readLine()) != null) {
				info = tab.split(input);
				String[] alleles = slash.split(info[1]);
				b73Allele[incount] = snpToByteMap.get(alleles[0]);
				try {nonb73Allele[incount] = snpToByteMap.get(alleles[1]);} catch(Exception e) {
					System.out.println("Error in second allele: " + input);
				}
				
				snpPositions[incount] = Integer.parseInt(info[3]);
				count = 0;
				String[] label = new String[11];
				for (int i = 0; i < 11; i++) {
					label[i] = new String(info[i]);
				}
				snpLabels[incount] = label;
				for (Integer col:columns) {
					snps[count++][incount] = snpToByteMap.get(info[col]);
				}
				incount++;
			}
			
			br.close();
			
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}

		System.out.println("Imputing snps...");
		imputeABHv2();
	
	}
	
	public String imputeParentsWithoutDependentSnps(String filename) {
		popmap = new HashMap<String, ArrayList<Integer>>();
		segregatingSites = new HashMap<String, boolean[]>();
		taxanames = new ArrayList<String>();
		inputFilename = filename;
		snplabels = new ArrayList<String[]>();
		ArrayList<Integer> columns = new ArrayList<Integer>();
		String[] info;
		String input;
		Pattern tab = Pattern.compile("\t");
		Pattern slash = Pattern.compile("/");
		int count;
		File inputFile = new File(filename);
		File dir = inputFile.getParentFile();
		System.out.println("Reading the snp file...");
		try {
			BufferedReader br = new BufferedReader(new FileReader(filename));
			
			//open the file and read the taxanames, keeping track of those that start with Z0, i.e. the NAM lines 
			br = new BufferedReader(new FileReader(filename));
			info = tab.split(br.readLine());
			for (int i = 0; i < 11; i++) headerForLabels[i] = info[i];
			count = 0;
			for (String str : info) {
				if (str.startsWith("Z0")) {
					taxanames.add(new String(str));
					columns.add(count);
				}
				count++;
			}
			numberOfTaxa = taxanames.size();
			
			//parse the population number from the taxa names, create a list of taxa column numbers for each population
			getpops();
			
			//load snp calls into a byte array
			LinkedList<byte[]> genoList = new LinkedList<byte[]>();
			byte[] prevsnp = null;
			int prevpos = -1;
			
			//what chromosome is this?
			input = br.readLine();
			info = tab.split(input, 5);
			
			while ((input = br.readLine()) != null) {
				info = tab.split(input);
				String[] label = new String[11];
				for (int i = 0; i < 11; i++) label[i] = new String(info[i]);
				byte[] snpvals = new byte[numberOfTaxa];
				count = 0;
				for (Integer col:columns) {
					snpvals[count++] = snpToByteMap.get(info[col]);
				}
				
				//is this snp from the same tag as the previous snp?
				int curpos = Integer.parseInt(label[3]);
				if (prevsnp == null) { //no previous snp, keep this one
					snplabels.add(label);
					genoList.add(snpvals);
					prevsnp = snpvals;
					prevpos = curpos;
				} else {
					if (curpos - prevpos < 128) {
						int bothNonMissCount = 0;
						int oneNonMissCount = 0;
						for (int i = 0; i < numberOfTaxa; i++) {
							if (snpvals[i] != N && prevsnp[i] != N) bothNonMissCount++;
							else if (snpvals[i] != N || prevsnp[i] != N) oneNonMissCount++;
						}
						//if one is present the other should be if they are from the same tag
						if (oneNonMissCount > 0.7 * ((double) bothNonMissCount)) { //probably not from the same tag
							snplabels.add(label);
							genoList.add(snpvals);
							prevsnp = snpvals;
							prevpos = curpos;
						}
					} else { //it is not, add it
						snplabels.add(label);
						genoList.add(snpvals);
						prevsnp = snpvals;
						prevpos = curpos;
					}
				}
			}
			
			numberOfSnps = genoList.size();
			snps = new byte[numberOfTaxa][numberOfSnps];
			b73Allele = new byte[numberOfSnps];
			nonb73Allele = new byte[numberOfSnps];
			snpPositions = new int[numberOfSnps];
			
			int snpCount = 0;
			for (byte[] snpval:genoList) {
				String[] label = snplabels.get(snpCount);
				String[] alleles = slash.split(label[1]);
				b73Allele[snpCount] = snpToByteMap.get(alleles[0]);
				try {nonb73Allele[snpCount] = snpToByteMap.get(alleles[1]);} catch(Exception e) {
					System.out.println("Error in second allele: " + input);
				}
				snpPositions[snpCount] = Integer.parseInt(label[3]);
				
				int taxaCount = 0;
				for (byte val:snpval) {
					snps[taxaCount][snpCount] = val; 
					taxaCount++;
				}
				snpCount++;
			}
			
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
		System.out.println("Imputing snps...");
		snpLabels = new String[numberOfSnps][];
		for (int i = 0; i < numberOfSnps; i++) {
			snpLabels[i] = snplabels.get(i);
		}
		return imputeABHv2();

	}
	
	private void getpops() {
		int count = 0;
		for (String taxon:taxanames) {
			String pop = taxon.substring(0,4);
			ArrayList<Integer> poplist = popmap.get(pop);
			if (poplist == null) {
				poplist = new ArrayList<Integer>();
				popmap.put(pop, poplist);
				segregatingSites.put(pop, new boolean[numberOfSnps]);
			}
			poplist.add(count++);
		}
		numberOfPopulations = popmap.size();
		
	}
	
	public String imputeABHv2() {
		
		imputedSnps = new byte[numberOfTaxa][numberOfSnps];
		Set<String> popnames = popmap.keySet();
		boolean[] segByPop = new boolean[numberOfSnps];
		//for each snp
		System.out.println("Making initial ABH calls.");
		for (int s = 0; s < numberOfSnps; s++) {
			byte[] snp = new byte[numberOfTaxa];
			for (int t = 0; t < numberOfTaxa; t++) {
				snp[t] = snps[t][s];
				imputedSnps[t][s] = N;
			}

			//which pops are segregating?
			HashMap<String, Boolean> segMap = whichPopsAreSegregating(popnames, snp);
			int numberSegregating = 0;
			for (Boolean b:segMap.values()) {
				if (b) numberSegregating++;
			}
			if (numberSegregating > 0) segByPop[s] = true; 
			
			//decide which allele is B73, which non B73
			//if all pops are segregating leave B73, nonB73 call as is
			//otherwise the major allele in the non-segregating populations is B73
			if (numberSegregating < 25 && numberSegregating > 0) {
				int[] nucleotideCounts = new int[16];

				//count the nucleotides in the non-segregating populations
				for (String pop:popnames) {
					if (!segMap.get(pop)) {
						ArrayList<Integer> taxa = popmap.get(pop);
						for (Integer t:taxa) {
							nucleotideCounts[snp[t]]++;
						}
					}
				}

				//find the major allele
				int major = 0;
				for (int i = 1; i < 4; i++) {
					if (nucleotideCounts[i] > nucleotideCounts[major]) {
						major = i;
					}
				}

				//if the major allele count > 1 then set that to B73, otherwise leave B73/nonB73 as is
				if (nucleotideCounts[major] > 1) {
					b73Allele[s] = (byte) major;

					//find the other allele in the segregating pops and set that to nonB73
					nucleotideCounts = new int[16];

					//first count nucleotides in the segregating populations
					for (String pop:popnames) {
						if (segMap.get(pop)) {
							ArrayList<Integer> taxa = popmap.get(pop);
							for (Integer t:taxa) {
								nucleotideCounts[snp[t]]++;
							}
						}
					}

					//find the allele with the highest count other than the B73 allele (major) and make that the nonB73 allele
					int other;
					if (major == 0) {
						other = 1;
					} else {
						other = 0;
					}

					for (int i = other + 1; i < 4; i++) {
						if (i != major && nucleotideCounts[i] > nucleotideCounts[other]) {
							other = i;
						}
					}

					if (nucleotideCounts[major] == 0 || nucleotideCounts[other] == 0) {
						System.err.println("Unexpected nucleotide counts at snp " + s +": B73 allele count = " + nucleotideCounts[major] + ", nonB73 allele count = " + nucleotideCounts[other] );
					} else {
						b73Allele[s] = (byte) major;
						nonb73Allele[s] = (byte) other;
					}
				}


			}
			
			//impute this snp to ABH
			if (numberSegregating > 0) {
				byte Aallele = b73Allele[s];
				byte Ballele = nonb73Allele[s];
				byte het = hetcodes[Aallele][Ballele];
				for (String pop:popnames) {
					if (segMap.get(pop)) {
						ArrayList<Integer> taxa = popmap.get(pop);
						for (Integer t:taxa) {
							byte val = snp[t];
							if (val == Aallele) imputedSnps[t][s] = A;
							else if (val == Ballele) imputedSnps[t][s] = B;
							else if (val == het) imputedSnps[t][s] = H;
						}
					}
					
				}
			}
			
			
		}

		int numberOfSnpsSegregating = 0;
		for (int s = 0; s < numberOfSnps; s++) if (segByPop[s]) numberOfSnpsSegregating++;
		System.out.println("number of snps segregating = " + numberOfSnpsSegregating );
		
		//once all SNPs are imputed, check to see if each SNP agrees with its neighbors
		System.out.println("checking LD");
		checkLD();
		
		System.out.println("Writing to file");
		return writeABHv2();
		
	}
	
	private String writeABHv2() {
		int ndx = inputFilename.lastIndexOf('.');
		String outputFilename = inputFilename.substring(0, ndx) + ".abhv2" + inputFilename.substring(ndx);
		
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outputFilename));
			bw.write("rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode");
			for (String taxon:taxanames) bw.write("\t" + taxon);
			bw.newLine();
			for (int s = 0; s < numberOfSnps; s++) {
				String[] label = snpLabels[s];
				bw.write(label[0]);
				for (int i = 1; i < 11; i++) {
					bw.write("\t");
					bw.write(label[i]);
				}
				for (int t = 0; t < numberOfTaxa; t++) {
					bw.write("\t");
					bw.write(byteToSnp[imputedSnps[t][s]]);
				}
				bw.newLine();
			}
			bw.close();
			return outputFilename;
		} catch (IOException e) {
			e.printStackTrace();
			return null;
		}
		
	}
	
	private HashMap<String, Boolean> whichPopsAreSegregating(Set<String> popnames, byte[] snp) {
		HashMap<String, Boolean> segMap = new HashMap<String, Boolean>();
		boolean[] seg = new boolean[popnames.size()];
		for (String pop:popnames) {
			int[] nucleotideCount = new int[16];
			ArrayList<Integer> taxa = popmap.get(pop);
			for (Integer t:taxa) {
				nucleotideCount[snp[t]]++;
			}
			
			//do at least two nucleotides (A,C,G, or T) have a count of 2 or more
			int minCount = 0;
			for (int i = 0; i < 4; i++) {
				if (nucleotideCount[i] >= 2) minCount++;
			}
			
			if (minCount >= 2) segMap.put(pop, true);
			else segMap.put(pop, false);
		}
		return segMap;
	}
	
	private void checkLD() {
		int checkInterval = 10;
		double minMatch = 0.8;
		
		for (int s = 0; s < numberOfSnps; s++) {
			//find the taxa for which the snp has been genotyped
			ArrayList<byte[]> taxaList = new ArrayList<byte[]>();
			for (int t = 0; t < numberOfTaxa; t++) {
				if (imputedSnps[t][s] != N) {
					taxaList.add(imputedSnps[t]);
				}
			}
			
			//if the number of taxa > 3
			if (taxaList.size() > 3) {
				
				int total = 0;
				int match = 0;
				//for each taxon 
				for (byte[] taxon:taxaList) {
					//count the snps within the check interval and the number that match the snp
					//check to the left
					int snpnum = s - 1;
					int limit = checkInterval;
					while (snpnum >= 0 && total <= limit) {
						if (taxon[snpnum] != N) {
							total++;
							if (taxon[snpnum] == taxon[s]) {
								match++;
							}
						} 
						snpnum--;
					}
					
					//check to the right
					snpnum = s + 1;
					limit = total + checkInterval;
					while (snpnum < numberOfSnps && total <= limit) {
						if (taxon[snpnum] != N) {
							total++;
							if (taxon[snpnum] == taxon[s]) {
								match++;
							}
						} 
						snpnum++;
					}
					
				}
				
				double overallMatch = ((double) match) / ((double) total);
				
				if (overallMatch < (1 - minMatch)) {	//if the overall match is less than minMatch switch A to B and B to A
					for (byte[] taxon:taxaList) {
						if (taxon[s] == A) {
							taxon[s] = B;
						} else if (taxon[s] == B) {
							taxon[s] = A;
						}
					}
				} else if (overallMatch < minMatch) {	//if the overall match is between 1-minMatch and minMatch, set the snp to N for all taxa
					for (byte[] taxon:taxaList) {
						taxon[s] = N;
					}
				}
				
				
			}
			
		}
	}
	
	
	
	public void imputeABH() {
		imputedSnps = new byte[numberOfTaxa][numberOfSnps];
		
		//count alleles by population. 
		//consider all populations with total allele count >= 15.
		int minCallsPerPopulation = 15;
		
		//a pop is non-segregating if there is only one allele or if the p(sample drawn from a binary dist. with p = 0.5)<.01
		//if all non segregating pops have the same allele then call that B73
		//use the value 16(X) to indicated that a snp has not been imputed.
		
		//diagnostic variables
		HashMap<String, Integer> popCountMap = new HashMap<String, Integer>();
		for (String pop : popmap.keySet()) {
			popCountMap.put(pop, 0);
		}
		
		//for each snp
		System.out.println("Making initial ABH calls.");
		for (int s = 0; s < numberOfSnps; s++) {
			//count alleles in each population
			byte[][] popAlleles = new byte[numberOfPopulations][2];
			int popcount = 0;
			ArrayList<Integer> segList = new ArrayList<Integer>();
			for (Map.Entry<String, ArrayList<Integer>> ent:popmap.entrySet()) {
				boolean[] segsites = segregatingSites.get(ent.getKey());
				int[] allelecounts = new int[16];
				for (Integer taxon : ent.getValue()) {
					allelecounts[snps[taxon][s]]++;
					imputedSnps[taxon][s] = N; //initialize imputed snps to N (missing)
				}
				
				//determine the major and minor allele (A,C,G,or T)
				int high = allelecounts[0];
				int next = 0;
				byte major = A;
				byte minor = X;
				
				for (int i = 1; i < 4; i++) {
					if (allelecounts[i] > high) {
						next = high;
						high = allelecounts[i];
						minor = major;
						major = (byte) i;
					} else if (allelecounts[i] > next) {
						next = allelecounts[i];
						minor = (byte) i;
					}
				}
				
				//is this segregating?
				if (high + next < minCallsPerPopulation) { //not enough samples
					popAlleles[popcount][0] = X;
					popAlleles[popcount][1] = X;
					segsites[s] = false;
				} else {
					double pval = Probability.binomial(next, high + next, 0.5);
					if (pval < .01) {  //not segregating
						popAlleles[popcount][0] = major;
						popAlleles[popcount][1] = X;
						segsites[s] = false;
					} else { //segregating
						popAlleles[popcount][0] = major;
						popAlleles[popcount][1] = minor;
						segsites[s] = true;
					}
				}
				
				popcount++;
			}
			
			//can B73 be called?
			byte B73 = X;
			byte notB73 = X;
			boolean consistent = true;
			
			//is this site segregating or not callable in all populations?
			//if so, add to allPopSegList and do not call for now by setting consistent = false
			int nonsegpops = 0;
			int notcallable = 0;
			for (int p = 0; p < numberOfPopulations; p++) {
				if (popAlleles[p][0] == X && popAlleles[p][1] == X) notcallable++;
				if (popAlleles[p][0] == X || popAlleles[p][1] == X) nonsegpops++;
			}
			
			if (nonsegpops == 0 && notcallable < numberOfPopulations) {
				consistent = false;
			}
			
			//is the B73 allele the same in all non-seg populations. That is all non-seg populations should have the same allele.
			for (int p = 0; p < numberOfPopulations; p++) if (consistent){
				if (popAlleles[p][1] == X && popAlleles[p][0] != X) {
					if (B73 == X) B73 = popAlleles[p][0];
					else if (B73 != popAlleles[p][0]) consistent = false;
				}
			}
			
			if (consistent && B73 != X) { //find nonB73
				for (byte[] alleles : popAlleles) {
					if (alleles[0] == X) {
						//do nothing
					} else if (alleles[0] != B73) {
						if (notB73 == X) notB73 = alleles[0];
						else if (notB73 != alleles[0]) consistent = false;
					} else if (alleles[1] != X) {
						if (notB73 == X) notB73 = alleles[1];
						else if (notB73 != alleles[1]) consistent = false;
					}
				}
			}
			
			//make ABH calls only in segregating populations
			if (consistent && notB73 != X && B73 != X) {
				b73Allele[s] = B73;
				nonb73Allele[s] = notB73;
				byte het = hetcodes[B73][notB73];
				for (String pop : popmap.keySet()) {
					if (segregatingSites.get(pop)[s]) {
						popCountMap.put(pop, popCountMap.get(pop) + 1);
						ArrayList<Integer> taxaList = popmap.get(pop);
						for (Integer t:taxaList) {
							if (snps[t][s] == B73) imputedSnps[t][s] = A;
							else if (snps[t][s] == notB73) imputedSnps[t][s] = B;
							else if (snps[t][s] == het) imputedSnps[t][s] = H;
						}
					}
				}
			}
			
		} //end snp loop: for (int s = 0; s < numberOfSnps; s++)
		
		//output population counts
//		System.out.println("Initial set segregating sites\nAll,Selected");
//		for (String pop:segregatingSites.keySet()) {
//			int nseg = 0;
//			for (boolean b:segregatingSites.get(pop)) if (b) nseg++;
//			System.out.println(pop + "," + nseg + "," + popCountMap.get(pop));
//		}
		
		//check snps against surrounding haplotypes
		System.out.println("Checking snps against enclosing haplotypes.");
		
		//which sites are segregating in at least one population
		boolean[] checkSite = new boolean[numberOfSnps];
		for (int s = 0; s < numberOfSnps; s++) checkSite[s] = false;
		
		for (boolean[] segsites:segregatingSites.values()) {
			for (int s = 0; s < numberOfSnps; s++) {
				checkSite[s] = checkSite[s] || segsites[s];
			}
		}
		
		//set non-segregating sites to N in imputedSnps
		for (int s = 0; s < numberOfSnps; s++) {
			if (!checkSite[s]) {
				for (int t = 0; t < numberOfTaxa; t++) imputedSnps[t][s] = N; 
			}
		}
		
		//generate some stats at this point
//		System.out.println("pop,A,B,N,total");
//		for (String pop : popmap.keySet()) {
//			ArrayList<Integer> taxaList = popmap.get(pop);
//			int Acount = 0;
//			int Bcount = 0;
//			int Ncount = 0;
//			int total = 0;
//			for (Integer t : taxaList) {
//				for (int s = 0; s < numberOfSnps; s++) {
//					if (imputedSnps[t][s] == A) Acount++;
//					else if (imputedSnps[t][s] == B) Bcount++;
//					else if (imputedSnps[t][s] == N) Ncount++;
//					total++;
//				}
//			}
//			System.out.println(pop + "," + Acount + "," + Bcount + "," + Ncount + "," + total);
//		}
//		System.exit(0);
		
		//this checks all SNPs against the surrounding haplotype. It checks SNPs in the input data against the haplotypes in imputed SNPs
		//that have been called based on segregation ratios in individual populations. ImputedSNPs is a more stringent dataset at this point.
		//All SNPs are tested against it. Those that are consistent with the surrounding haplotype are added to ImputedSNPs.
		imputeSnpsBasedOnHaplotypes(0.85);
		
		//generate some stats at this point
		System.out.println("pop,A,B,H,N,total");
		for (String pop : popmap.keySet()) {
			ArrayList<Integer> taxaList = popmap.get(pop);
			int Acount = 0;
			int Bcount = 0;
			int Ncount = 0;
			int Hcount = 0;
			int total = 0;
			for (Integer t : taxaList) {
				for (int s = 0; s < numberOfSnps; s++) {
					if (imputedSnps[t][s] == A) Acount++;
					else if (imputedSnps[t][s] == B) Bcount++;
					else if (imputedSnps[t][s] == N) Ncount++;
					else if (imputedSnps[t][s] == H) Hcount++;
					total++;
				}
			}
			System.out.println(pop + "," + Acount + "," + Bcount + "," + Hcount + "," + Ncount + "," + total);
		}

//		System.exit(0);
		
		//for individual taxa, look for probable errors
		//change to missing if an A has hapsize B's on both sides change to B, likewise A's & H's
		//change on 10/17/2011: leave H calls unchanged. The previous version had the effect of setting all H's to 0.
		System.out.println("Identifying errors.");
		int hapsize = 6;
		for (String popname : popmap.keySet()) {
			ArrayList<Integer> taxaList = popmap.get(popname);
			for (Integer t : taxaList) {
				byte[] taxonSnps = imputedSnps[t];

				//count non-missing snps
				int nonmissCount = 0;
				for (byte snp:taxonSnps) {
					if (snp != N) nonmissCount++;
				}
				if (nonmissCount < 10) {
					System.out.println(nonmissCount + " non-missing data points for " + taxanames.get(t) );
					for (int s = 0; s < numberOfSnps; s++) imputedSnps[t][s] = N;

				} else {
					//create an index of non-missing counts
					int[] nonmissIndex = new int[nonmissCount];
					int count = 0;
					for (int s = 0; s < numberOfSnps; s++) {
						if (taxonSnps[s] != N) nonmissIndex[count++] = s;
					}

					//count number of adjacent matching snps from left to right
					int[] leftcount = new int[nonmissCount];
					leftcount[0] = 1;
					for (int s = 1; s < nonmissCount; s++) {
						if (taxonSnps[nonmissIndex[s]] == taxonSnps[nonmissIndex[s-1]]) leftcount[s] = leftcount[s-1] + 1;
						else leftcount[s] = 1;
					}

					//count number of adjacent matching snps from right to left
					int[] rightcount = new int[nonmissCount];
					rightcount[nonmissCount - 1] = 1;
					for (int s = nonmissCount - 2; s >= 0; s--) {
						if (taxonSnps[nonmissIndex[s]] == taxonSnps[nonmissIndex[s+1]]) rightcount[s] = rightcount[s+1] + 1;
						else rightcount[s] = 1;
					}

					
					
					
					//if left and right counts = one and left - 1 >= hapsize and right + 1 >= hapsize then set to N
					//also at position s, if leftcount[s-1] = s and rightcount[s+1] >= hapsize, then set to N
					//similarly if leftcount[s-1] >= hapsize and rightcount[s+1] == nmcount - s, then set to N
					//leave all Hs (hets) as is
					for (int s = 0; s < nonmissCount; s++) {
						if (taxonSnps[nonmissIndex[s]] != H) {
							if (leftcount[s] == 1 && rightcount[s] == 1) {
								int sIndex = nonmissIndex[s];
								if (s == 0) {
									if (rightcount[s+1] >= hapsize) taxonSnps[sIndex] = N;
								} else if (s <= hapsize) {
									if (rightcount[s+1] >= hapsize && leftcount[s - 1] == s) taxonSnps[sIndex] = N;
								} else if (s == nonmissCount - 1) {
									if (leftcount[s-1] >= hapsize) taxonSnps[sIndex] = N;
								} else if (s >= nonmissCount - hapsize - 1) {
									if (leftcount[s-1] >= hapsize && rightcount[s+1] == nonmissCount - s) taxonSnps[sIndex] = N;
								} else {
									if (leftcount[s-1] >= hapsize && rightcount[s+1] >= hapsize) taxonSnps[sIndex] = N;
								}
							}

							//if snp pairs are within 128 bp treat as a single snp
							if (s < nonmissCount - 1 && leftcount[s] == 1 && rightcount[s] == 2 && leftcount[s+1] == 2 && rightcount[s+1] == 1) {
								int sIndexLeft = nonmissIndex[s];
								int sIndexRight = nonmissIndex[s+1];
								if (snpPositions[sIndexRight] - snpPositions[sIndexLeft] < 128) { //the snps are within 128 bp of each other
									if (s == 0) {
										if (rightcount[s+2] >= hapsize) {
											taxonSnps[sIndexRight] = N;
											taxonSnps[sIndexLeft] = N;
										}
									} else if (s <= hapsize) {
										if (rightcount[s+2] >= hapsize && leftcount[s - 1] == s) {
											taxonSnps[sIndexRight] = N;
											taxonSnps[sIndexLeft] = N;
										}
									} else if (s == nonmissCount - 2) {
										if (leftcount[s-1] >= hapsize) {
											taxonSnps[sIndexRight] = N;
											taxonSnps[sIndexLeft] = N;
										}
									} else if (s >= nonmissCount - hapsize - 2) {
										if (leftcount[s-1] >= hapsize && rightcount[s+2] == nonmissCount - s - 1) {
											taxonSnps[sIndexRight] = N;
											taxonSnps[sIndexLeft] = N;
										}
									} else {
										if (leftcount[s-1] >= hapsize && rightcount[s+2] >= hapsize) {
											taxonSnps[sIndexRight] = N;
											taxonSnps[sIndexLeft] = N;
										}
									}
								}
								s++;
							}
						}
					}
				}
			}
		}

	}
	
	
	/**
	 * This function finds the proportion of snps five snps on either side of snpIndex that match the A/B call of snpIndex. 
	 * Fewer than five may be used at the ends of the chromosomes.
	 * @param snpIndex	the index of the snp being tested
	 * @param popnames	names of populations
	 * @return for each population, and int array containing the number of A matches, the number of B matches and the total count
	 */
	public int[][] testEnclosingHaplotype(int snpIndex, List<String> popnames) {
		//test this snp against enclosing haplotypes for each population
		int[][] nmatch = new int[popnames.size()][3];
		int sitecount;
		int hapsize = 5;
		byte b73 = b73Allele[snpIndex];
		byte nonb73 = nonb73Allele[snpIndex];
		int popcount = 0;
		for (String popname:popnames) { //iterate through populations
			ArrayList<Integer> taxaList = popmap.get(popname);
			int matchCount = 0;
			int AmatchCount = 0;
			int BmatchCount = 0;
			int totalCount = 0; 
			for (Integer t:taxaList) { //iterate through taxa in one population
				byte snp = snps[t][snpIndex];
				byte snptype;
				if (snp == b73) snptype = A;
				else if (snp == nonb73) snptype = B;
				else snptype = N;
				if (snptype != N) {
					//count matches to the left
					sitecount = 0;
					int testsite = snpIndex - 1;
					while (sitecount < hapsize && testsite >= 0) {
						byte testsnp = imputedSnps[t][testsite--];
//						if (testsnp == A || testsnp == B) {
//							totalCount++;
//							sitecount++;
//							if (testsnp == snptype) matchCount++;
//						}
						if (testsnp == A) {
							totalCount++;
							sitecount++;
							if (testsnp == snptype) {
								matchCount++;
								AmatchCount++;
							}
						}
						else if (testsnp == B) {
							totalCount++;
							sitecount++;
							if (testsnp == snptype) {
								matchCount++;
								BmatchCount++;
							}
						}
					}
					
					//count matches to the right
					sitecount = 0;
					testsite = snpIndex + 1;
					while (sitecount < hapsize && testsite < numberOfSnps) {
						byte testsnp = imputedSnps[t][testsite++];
//						if (testsnp == A || testsnp == B) {
//							totalCount++;
//							sitecount++;
//							if (testsnp == snptype) matchCount++;
//						}
						if (testsnp == A) {
							totalCount++;
							sitecount++;
							if (testsnp == snptype) {
								matchCount++;
								AmatchCount++;
							}
						}
						else if (testsnp == B) {
							totalCount++;
							sitecount++;
							if (testsnp == snptype) {
								matchCount++;
								BmatchCount++;
							}
						}
					}
				}
			}
			
//			if (totalCount == 0) pmatch[popcount] = Double.NaN;
//			else pmatch[popcount] = ((double) matchCount) / ((double) totalCount);
			
			nmatch[popcount][0] = AmatchCount;
			nmatch[popcount][1] = BmatchCount;
			nmatch[popcount][2] = totalCount;
			
			popcount++;
		}
		
		return nmatch;
	}
	
	/*
	public void imputeSnpBasedOnHaplotype(int snpIndex, double limit, List<String> popnames) {
		double[] pmatch = testEnclosingHaplotype(snpIndex, popnames);
		int overcount = 0;
		int undercount = 0;
		for (double p:pmatch) {
			if (p > limit) overcount++;
			else if (p < (1 - limit)) undercount++;
		}
		
		if (overcount > 0 && undercount == 0) { //impute snp values in the pops with overcount>0 and set to missing in the rest
			imputeSnpValuesInSegregatingPopulations(snpIndex, pmatch, limit);
			
		} else if (overcount == 0 && undercount > 0){  //switch b73/nonb73 and impute snp values
			byte temp = b73Allele[snpIndex];
			b73Allele[snpIndex] = nonb73Allele[snpIndex];
			nonb73Allele[snpIndex] = temp;
			imputeSnpValuesInSegregatingPopulations(snpIndex, pmatch, limit);
		} else { //set all snp values to missing
			for (int t = 0; t < numberOfTaxa; t++) imputedSnps[t][snpIndex] = N;
		}
	}
	*/
	
	/**
	 * Infers whether called SNPs are either A type, B type, or cannot be called reliably. 
	 * The function sets values to A/B/N in the array imputedSnps. Only existing calls are imputed as A/B. Missing data is not imputed.
	 * @param limit the proportion of SNPs that must match for the SNP to be considered okay
	 */
	public void imputeSnpsBasedOnHaplotypes (double limit) {
		byte[][] impute2 = new byte[numberOfTaxa][numberOfSnps];
		List<String> popnames = new ArrayList<String>(popmap.keySet());
		Collections.sort(popnames);
		for (int s = 0; s < numberOfSnps; s++) {
			//pmatch contains the proportion of neigboring SNPs with A/B matching the SNP at position s in each population
			//overcount contains the number of populations in which this SNP matches the surrounding haplotype
			//undercount contains the number of populations in which this SNP would match the surrounding haplotype if the SNP call were changed,
			//indicating that the A/B designation of the SNP should be changed
			int[][] nMatch = testEnclosingHaplotype(s, popnames);
			int n = nMatch.length;
			double[] pmatch = new double[n];
			for (int i = 0; i < n; i++) {
				if (nMatch[i][2] == 0) pmatch[i] = Double.NaN;
				else pmatch[i] = ((double)(nMatch[i][0] + nMatch[i][1])) / ((double) nMatch[i][2]);
			}
			
			int overcount = 0;
			int undercount = 0;
			for (double p:pmatch) {
				if ( !Double.isNaN(p) && p > limit ) overcount++;
				else if ( !Double.isNaN(p) && p < (1 - limit) ) undercount++;
			}
			
			//if both overcount > 0 and undercount > 0, there is a problem and the snp should not be used
			if (overcount > 0 && undercount == 0 || overcount == 0 && undercount > 0) { //impute snp values

				if (overcount == 0 && undercount > 0){  //switch b73/nonb73 
					byte temp = b73Allele[s];
					b73Allele[s] = nonb73Allele[s];
					nonb73Allele[s] = temp;
				}
				
				//impute snp values for populations in which the SNP matches the haplotype
				//and there are at least two B calls and two A calls to ensure that the population is segregating
				int popcount = 0;
				byte b73 = b73Allele[s];
				byte nonb73 = nonb73Allele[s];
				byte het = hetcodes[b73][nonb73];
				for (String popname:popnames) {
					double pval = pmatch[popcount];
					int Acount = nMatch[popcount][0];
					int Bcount = nMatch[popcount][1];
					int mincount = Math.max(2, (int)(0.05 * (Acount + Bcount) ));
					if ((pval > limit || pval < 1 - limit) && Acount >= mincount && Bcount >= mincount){
						for (Integer t:popmap.get(popname)) {
							if (snps[t][s] == b73) {
								impute2[t][s] = A;
							} else if (snps[t][s] == nonb73) {
								impute2[t][s] = B;
							} else if (snps[t][s] == het) {	//added het calls on 10/17/11
								impute2[t][s] = H;
							} else {
								impute2[t][s] = N;
							}
						}
					} else {
						for (Integer t:popmap.get(popname)) {
							impute2[t][s] = N;
						}
					}
					popcount++;
				}

			} else { //set all snp values to missing
				for (int t = 0; t < numberOfTaxa; t++) impute2[t][s] = N;
			}
			
			
		}
		
		imputedSnps = impute2;
	}
	
	public void imputeSnpValuesInSegregatingPopulations(int snpIndex, double[] pmatch, List<String> popnames, double limit) {
		byte b73 = b73Allele[snpIndex];
		byte nonb73 = nonb73Allele[snpIndex];
		byte het = getHetCode(b73, nonb73);
		int popcount = 0;
		for (String popname:popnames) {
			double pval = pmatch[popcount];
			if (pval > limit || pval < 1 - limit){
				for (Integer t:popmap.get(popname)) {
					if (snps[t][snpIndex] == b73) {
						imputedSnps[t][snpIndex] = A;
					} else if (snps[t][snpIndex] == nonb73) {
						imputedSnps[t][snpIndex] = B;
					} else {
						imputedSnps[t][snpIndex] = N;
					}
				}
			} else {
				for (Integer t:popmap.get(popname)) {
					imputedSnps[t][snpIndex] = N;
				}
			}
			popcount++;
		}
	}
	
	public static void callHets(byte[] snps) {
		int sizecount = 1;
		int segcount = 0;
		int hetstart = 0;
		int segstart = 0;
		byte prevbyte = -1;
		int hetlimit = 12;  //the largest number of contiguous calls of the same type allowed to call a het segment
		int minsegments = 1;  //the minimum number of contiguous het segments required to call a het
		int numberOfSnps = snps.length;
		
		for (int s = 0; s < numberOfSnps; s++) {
			byte thisbyte = snps[s];
			if (thisbyte != N) {
				if ( thisbyte == prevbyte) {
					sizecount++;
				} else {
					if (sizecount > hetlimit) {  //reset everything
						if (segcount >= minsegments && segstart != 0) {  //set bytes between start and current segment (exclusive) to H except at start of chromosome
							for (int h = hetstart; h < segstart; h++) {
								if (snps[h] != N) snps[h] = H;
							}
						}
						sizecount = 1;
						segcount = 0;
						hetstart = s;
						segstart = s;
						prevbyte = thisbyte;
					} else {  //count as a het segment, reset segstart, and keep going
						sizecount = 1;
						segcount++;
						segstart = s;
						prevbyte = thisbyte;
					}
				}
			}
		}
		
		//test for stuff at final snp
		if (sizecount > hetlimit) {
			if (segcount >= minsegments) {  //set bytes between start and current segment (exclusive) to H
				for (int h = hetstart; h < segstart; h++) {
					if (snps[h] != N) snps[h] = H;
				}
			}
		} else {
			if (segcount >= minsegments) {  //set bytes between start and this site (inclusive) to H. If only the final segment is short, it is not set to H.
				for (int h = hetstart; h < numberOfSnps; h++) {
					if (snps[h] != N) snps[h] = H;
				}
			}
		}
		
	}
	
	public static void callHetsFromABHFile(String inputfile, String outputfile) {
		Pattern tab = Pattern.compile("\t");
		HashMap<String, Byte> snpToByte = new HashMap<String, Byte>();
		snpToByte.put("A", A);
		snpToByte.put("H", H);
		snpToByte.put("B", B);
		snpToByte.put("N", N);
		ArrayList<String> labelList = new ArrayList<String>();
		byte snps[][];
		
		System.out.println("Reading " + inputfile + ", writing " + outputfile);
		try {
			File outFile = new File(outputfile);
//			if (outFile.exists()) {
//				String msg = "Output file already exists and will be overwritten.";
//				JOptionPane.showMessageDialog(null, msg, "Output Exists", JOptionPane.ERROR_MESSAGE);
////				return;
//			}
			BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));
			
			//count lines in the input file
			BufferedReader br = new BufferedReader(new FileReader(inputfile));
			int linecount = 0;
			while (br.readLine() != null) linecount++;
			br.close();
			int numberOfSnps = linecount - 1;
			
			String input;
			String[] info;
			br = new BufferedReader(new FileReader(inputfile));
			input = br.readLine();
			info = tab.split(input);
			int numberOfTaxa = info.length - 11;
			snps = new byte[numberOfTaxa][numberOfSnps];
			bw.write(input);
			bw.write("\n");
			
			int count = 0;
			while ((input = br.readLine()) != null) {
//				if (count % 1000 == 0) System.out.println("Processing line " + count);
				info = tab.split(input);
				StringBuilder sb = new StringBuilder(info[0]);
				for (int i = 1; i < 11; i++) sb.append("\t").append(info[i]);
				labelList.add(sb.toString());

				for (int t = 0; t < numberOfTaxa; t++) {
					snps[t][count] = snpToByte.get(info[t+11]);
				}
				count++;
			}
			br.close();
			
			for (int t = 0; t < numberOfTaxa; t++) {
				callHets(snps[t]);
			}
			
			//write to output
			for (int s = 0; s < numberOfSnps; s++) {
				bw.write(labelList.get(s));
				for (int t = 0; t < numberOfTaxa; t++) {
					bw.write("\t");
					bw.write(byteToSnp[snps[t][s]]);
				}
				bw.write("\n");
			}
			bw.close();
			
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
		System.out.println("Finished calling hets.");
	}
	
	public byte getHetCode(byte a1, byte a2) {
		if (a1 < R & a2 < R) return hetcodes[a1][a2];
		if (a1 == a2) return a1;
		if ((a1 == plus && a2 == minus) || (a1 == plus && a2 == minus)) return zero;
		return N;
	}
	
	public void writeTaxa(String filenamebase) {
		
		for (String popname:popmap.keySet()) {
			ArrayList<Integer> taxaList = popmap.get(popname);
			String filename;
			if (filenamebase == null) filename ="C:/users/peter/temp/gbs/" + popname + "_chr" + chr + "_imputedsnps.txt";
			else filename = filenamebase + "_chr" + chr + "_" + popname + ".txt";
			try {
				BufferedWriter bw = new BufferedWriter(new FileWriter(filename));
				bw.write("rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode");
				for (Integer t:taxaList) bw.write("\t" + taxanames.get(t));
				bw.write("\n");
				for (int i = 0; i < numberOfSnps; i++) {
					Iterator<Integer> it = taxaList.iterator();
					boolean allN = true;
					while (allN && it.hasNext()) {
						if (imputedSnps[it.next()][i] != N) allN = false;
					}
					if (!allN) {
						String[] label = snplabels.get(i);
						bw.write(label[0]);
						for (int j = 1; j < 11; j++) {
							bw.write("\t");
							bw.write(label[j]);
						}
						for (Integer t:taxaList) {
							bw.write("\t");
							bw.write(byteToSnp[imputedSnps[t][i]]);
						}
						bw.write("\n");
					}
				}
				bw.close();
			} catch (IOException e) {
				e.printStackTrace();
				System.exit(-1);
			}
		}
		
		System.out.println("finished writing files");
	}
	
	public static void countClasses() {
		String filename = "C:/Projects/NAM/namgbs/allfusion_110401.LD/allLD/allfusion_110401.lh.ld.c4.dedupe.txt";
		Pattern tab = Pattern.compile("\t");
		String input;
		String[] info;
		int[][] popcounts = new int[27][2];
		int[] pop;
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(filename));
			info = tab.split(br.readLine());
			int n = info.length;
			pop = new int[n - 11];
			for (int i = 11; i < n; i++) {
				pop[i - 11] = Integer.parseInt(info[i].substring(2,4));
			}
			
			while ((input = br.readLine()) != null) {
				info = tab.split(input);
				for (int i = 11; i < n; i++) {
					int taxon = i - 11;
					int popnum = pop[taxon];
					popcounts[popnum][1]++;
					if (!info[taxon].equals("N")) popcounts[popnum][0]++;
				}
			}
			br.close();
			
			for (int i = 0; i < 27; i++) {
				int totalcount = popcounts[i][1];
				int nonNcount = popcounts[i][0];
				double prop;
				if (totalcount != 0) prop = ((double) nonNcount) / ((double) totalcount);
				else prop = Double.NaN;
				
				System.out.println("pop " + i + ": " + nonNcount + ", " + totalcount + ", " + prop);
			}
			
		} catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
	}
	
	public static void imputeForHapmap(String inputfile, String outputfile) {
		String input;
		String[] info;
		Pattern tab = Pattern.compile("\t");
		String fill = "\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
		
		byte A = (byte) 0;
		byte H = (byte) 1;
		byte B = (byte) 2;
		byte N = (byte) 3;
		
		HashMap<String, Byte> string2byte = new HashMap<String, Byte>();
		string2byte.put("A", A);
		string2byte.put("H", H);
		string2byte.put("B", B);
		string2byte.put("N", N);
		
		String[] b2s = new String[] {"A", "M", "C", "N"};
		
		try{
			BufferedReader br = new BufferedReader(new FileReader(inputfile));
			BufferedWriter bw = new BufferedWriter(new FileWriter(outputfile));
			
			//count lines
			int linecount = 0;
			while (br.readLine() != null) linecount++;
			br.close();
			int numberOfSnps = linecount - 1;
			
			br = new BufferedReader(new FileReader(inputfile));
			info = tab.split(br.readLine());
			int numberOfTaxa = info.length - 4;
			byte[][] snps = new byte[numberOfTaxa][numberOfSnps];
			ArrayList<String> snplabel = new ArrayList<String>();
			bw.write("rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcentre\tprotLSID\tassayLSID\tpanelLSID\tQCcode");
			for (int i = 4; i < info.length; i++) bw.write("\t" + info[i]);
			bw.write("\n");
			
			System.out.println("reading data");
			linecount = 0;
			while ((input = br.readLine()) != null) {
				info = tab.split(input);
				snplabel.add(info[0] + "\t" + info[1] + "\t" + info[2] + "\t" + info[3]);
				for (int t = 0; t < numberOfTaxa; t++) {
					snps[t][linecount] = string2byte.get(info[t + 4]);
				}
				linecount++;
			}
			
			//impute snps
			System.out.println("imputing snps");
			for (int t = 0; t < numberOfTaxa; t++) {
				byte prevsnp = -1;
				int prevpos = -1;
				byte[] asnp = snps[t];
				for (int s = 0; s < numberOfSnps; s++) {
					if (asnp[s] != N) {
						if (prevsnp == -1) {
							for (int i = 0; i < s; i++) asnp[i] = asnp[s];
						} else if(asnp[s] == prevsnp) {
							for (int i = prevpos; i < s; i++) asnp[i] = prevsnp;
						}
						prevsnp = asnp[s];
						prevpos = s;
					}
				}
				if (prevpos > 0) {
					for (int i = prevpos; i < numberOfSnps; i++) asnp[i] = prevsnp;
				}
			}
			
			System.out.println("writing output");
			for (int s = 0; s < numberOfSnps; s++) {
				bw.write(snplabel.get(s));
				bw.write(fill);
				for (int t = 0; t < numberOfTaxa; t++) {
					bw.write("\t");
					bw.write(b2s[snps[t][s]]);
				}
				bw.write("\n");
			}
			
			
			br.close();
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

	public static void imputeNumericGenotype(String inputfile, String outputfile) {
		String input;
		String[] info;
		Pattern tab = Pattern.compile("\t");
		
		byte A = (byte) 0;
		byte H = (byte) 1;
		byte B = (byte) 2;
		byte N = (byte) 3;
		
		HashMap<String, Byte> string2byte = new HashMap<String, Byte>();
		string2byte.put("A", A);
		string2byte.put("H", H);
		string2byte.put("B", B);
		string2byte.put("N", N);
		
		String[] b2s = new String[] {"0", "1", "2", "N"};
		
		try{
			BufferedReader br = new BufferedReader(new FileReader(inputfile));
			BufferedWriter bw = new BufferedWriter(new FileWriter(outputfile));
			
			//count lines
			int linecount = 0;
			while (br.readLine() != null) linecount++;
			br.close();
			int numberOfSnps = linecount - 1;
			
			br = new BufferedReader(new FileReader(inputfile));
			info = tab.split(br.readLine());
			int numberOfTaxa = info.length - 4;
			byte[][] snps = new byte[numberOfTaxa][numberOfSnps];
			String[] taxa = new String[numberOfTaxa];
			for (int t = 0; t < numberOfTaxa; t++) taxa[t] = info[t + 4];
			ArrayList<String> snplabel = new ArrayList<String>();
			bw.write("<Numeric>\n");
			bw.write("<Marker>");
			for (int s = 0; s < numberOfSnps; s++) {
				bw.write("\t");
				bw.write("s" + s);
			}
			bw.write("\n");
			
			System.out.println("reading data");
			linecount = 0;
			while ((input = br.readLine()) != null) {
				info = tab.split(input);
				snplabel.add(info[0] + "\t" + info[1] + "\t" + info[2] + "\t" + info[3]);
				for (int t = 0; t < numberOfTaxa; t++) {
					snps[t][linecount] = string2byte.get(info[t + 4]);
				}
				linecount++;
			}
			
			//impute snps
			System.out.println("imputing snps");
			for (int t = 0; t < numberOfTaxa; t++) {
				byte prevsnp = -1;
				int prevpos = -1;
				byte[] asnp = snps[t];
				for (int s = 0; s < numberOfSnps; s++) {
					if (asnp[s] != N) {
						if (prevsnp == -1) {
							for (int i = 0; i < s; i++) asnp[i] = asnp[s];
						} else if(asnp[s] == prevsnp) {
							for (int i = prevpos; i < s; i++) asnp[i] = prevsnp;
						}
						prevsnp = asnp[s];
						prevpos = s;
					}
				}
				if (prevpos > 0) {
					for (int i = prevpos; i < numberOfSnps; i++) asnp[i] = prevsnp;
				}
			}
			
			System.out.println("writing output");
			
			for (int t = 0; t < numberOfTaxa; t++) {
				bw.write(taxa[t]);
				for (int s = 0; s < numberOfSnps; s++) {
					bw.write("\t");
					byte val = snps[t][s];
					bw.write(b2s[val]);
				}
				bw.write("\n");
			}
			
			
			br.close();
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

	/**
	 * compares gbs calls to 1106 marker genotyping array calls
	 * @param inputfile a hapmap formatted file of genotypes for a chromosome with hets imputed
	 * @return a HashMap with taxa names as keys and matchcount, total nam homozygote count, hetcount, totalcount as values
	 */
	public static HashMap<String, int[]> compareToNAMarray(String inputfile, int chr) {
		String arrayFileName = "C:/Projects/NAM/data/markergenotypes062508.txt";
		String mapFileName = "C:/Projects/NAM/data/markers061208agpv2.txt";
		int[] agpMarkerPositions = null;
		int[] namMarkerNumber = null;
		int numberOfMarkers = 0;
		HashMap <String, int[]> genotypeMap = new HashMap <String, int[]>();
		HashMap <String, int[]> taxaCountMap = new HashMap <String, int[]>();
		
		String input;
		String[] data;
		Pattern tab = Pattern.compile("\t");
		
		int startMarker = 0;
		
		HashMap <String, Integer> string2int = new HashMap <String, Integer>();
		string2int.put("0.0", 0);
		string2int.put("0.5", 1);
		string2int.put("1.0", 2);
		string2int.put("1.5", 3);
		string2int.put("2.0", 4);
		
		//read in the map file
		try { 
			
			//count the number of markers in this chromosome
			BufferedReader br = new BufferedReader(new FileReader(mapFileName));
			br.readLine();
			int snpcount = 0;
			int linecount = 0;
			input = br.readLine();
			data = tab.split(input);
			int chrnum = Integer.parseInt(data[1]);
			while (chrnum < chr) {
				input = br.readLine();
				data = tab.split(input);
				chrnum = Integer.parseInt(data[1]);
				linecount++;
			}
			startMarker = linecount;
			
			while (chrnum == chr) {
				snpcount++;
				input = br.readLine();
				if (input != null) {
					data = tab.split(input);
					chrnum = Integer.parseInt(data[1]);
				} else {
					chrnum = -1;
				}
				linecount++;
			}
			numberOfMarkers = snpcount;
			br.close();
			
			//dimension the marker position array and read the agpv2 positions into it
			agpMarkerPositions = new int[numberOfMarkers];
			namMarkerNumber = new int[numberOfMarkers];
			br = new BufferedReader(new FileReader(mapFileName));
			br.readLine();
			snpcount = 0;
			input = br.readLine();
			data = tab.split(input);
			chrnum = Integer.parseInt(data[1]);
			while (chrnum < chr) {
				input = br.readLine();
				data = tab.split(input);
				chrnum = Integer.parseInt(data[1]);
			}
			while (chrnum == chr) {
				agpMarkerPositions[snpcount] = Integer.parseInt(data[5]);
				namMarkerNumber[snpcount++] = Integer.parseInt(data[3]) - 1;
				input = br.readLine();
				if (input != null) {
					data = tab.split(input);
					chrnum = Integer.parseInt(data[1]);
				} else {
					chrnum = -1;
				}
			}
			br.close();
			
			
		} catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
		
		//read in the markergenotype data for this chromosome

		try { 
			BufferedReader br = new BufferedReader(new FileReader(arrayFileName));
			br.readLine();
			br.readLine(); //discard two header rows
			
			while ((input = br.readLine()) != null) {
				data = tab.split(input);
				String taxon = data[0];
				int[] geno = new int[1106];
				for (int m = 0; m < 1106; m++) {
					geno[m] = string2int.get(data[m + 5]);
				}
				genotypeMap.put(taxon, geno);
			}
			
			br.close();
		} catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
		
		//read the snp data and get counts

		try { 
			BufferedReader br = new BufferedReader(new FileReader(inputfile));
			input = br.readLine();
			data = tab.split(input);
			int numberOfTaxa = data.length - 11;
			String[] taxanames = new String[numberOfTaxa];
			for (int t = 0; t < numberOfTaxa; t++) {
				String taxon = data[t + 11];
				taxanames[t] = taxon.substring(0, taxon.indexOf(':'));
				taxaCountMap.put(taxanames[t], new int[4]);
			}
			
			while ((input = br.readLine()) != null) {
				data = tab.split(input);
				//get the position of this snp and find the flanking nam markers
				int snpPosition = Integer.parseInt(data[3]);
				int snpIndex = Arrays.binarySearch(agpMarkerPositions, snpPosition);
				int[] flankingMarkers = new int[2];
				if (snpIndex >= 0) {
					flankingMarkers[0] = snpIndex;
					flankingMarkers[1] = snpIndex;
				} else {
					int insertionPoint = - (snpIndex + 1);
					if (insertionPoint == 0) {
						flankingMarkers[0] = 0;
						flankingMarkers[1] = 0;
					} else if (insertionPoint == numberOfMarkers) {
						flankingMarkers[0] = numberOfMarkers - 1;
						flankingMarkers[1] = numberOfMarkers - 1;
					} else {
						flankingMarkers[0] = insertionPoint - 1;
						flankingMarkers[1] = insertionPoint;
					}
				}
				
				for (int t = 0; t < numberOfTaxa; t++) {
					String geno = data[t + 11];
					int[] namgeno = genotypeMap.get(taxanames[t]);
					if (namgeno != null) {
						int leftmarker = namgeno[namMarkerNumber[flankingMarkers[0]]];
						int rightmarker = namgeno[namMarkerNumber[flankingMarkers[1]]];
						int markersum = leftmarker + rightmarker;
						int[] theseCounts = taxaCountMap.get(taxanames[t]);
						if (geno.equals("A")) {
							theseCounts[3]++;
							if (markersum == 0) {
								theseCounts[0]++;
								theseCounts[1]++;
							} else if (markersum == 8) {
								theseCounts[1]++;
							}
						} else if (geno.equals("B")) {
							theseCounts[3]++;
							if (markersum == 8) {
								theseCounts[0]++;
								theseCounts[1]++;
							} else if (markersum == 0) {
								theseCounts[1]++;
							}
						} else if (geno.equals("H")) {
							theseCounts[2]++;
							theseCounts[3]++;
						} 
					}
				}
			}
			
			br.close();
		} catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
		
		//code to output results to a file
//		try {
//			BufferedWriter bw = new BufferedWriter(new FileWriter("c:/users/peter/temp/gbs/matchcounts.chr" + chr + ".txt"));
//			bw.write("taxon\tmatchAB\ttotalAB\thets\ttotal\n");
//			for (String t : taxaCountMap.keySet()) {
//				bw.write(t);
//				int counts[] = taxaCountMap.get(t);
//				for (int i = 0; i < 4; i++) {
//					bw.write("\t" + counts[i]);
//				}
//				bw.newLine();
//			}
//			bw.close();
//		} catch(IOException e) {
//			
//		}
		return taxaCountMap;
	}
	
	public static void countMatches() {
		int pops[] = new int[]{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,19,20,21,22,23,24,25,26};
		for (int pop:pops) {
			try {
				BufferedWriter bw = new BufferedWriter(new FileWriter("C:/Projects/NAM/namgbs/datasets/cov20/matchcounts.pop" + pop + ".txt"));
				bw.write("chromosome\ttaxon\tmatchAB\ttotalAB\tHcount\ttotal\tprmatch\tprhet\n");

				for (int chr = 1; chr <= 10; chr++) {
					String znum;
					if (pop < 10) znum = "Z00" + pop;
					else znum = "Z0" + pop;
					String filename = "C:/Projects/NAM/namgbs/datasets/cov20/cov20_chr" + chr + "_" + znum + ".hets.txt";
					HashMap<String, int[]> taxonCounts = compareToNAMarray(filename, chr);
					for (String t : taxonCounts.keySet()) {
						bw.write(chr + "\t" + t);
						int counts[] = taxonCounts.get(t);
						for (int i = 0; i < 4; i++) {
							bw.write("\t" + counts[i]);
						}
						if (counts[1] > 0) bw.write("\t" + (((double) counts[0])/counts[1]));
						else bw.write("\t");
						if (counts[3] > 0) bw.write("\t" + (((double) counts[2])/counts[3]));
						else bw.write("\t");
						bw.newLine();
					}

				}
				bw.close();
			} catch(IOException e) {
				e.printStackTrace();
				System.exit(-1);
			}
			System.out.println("Finished summarizing population " + pop + ".");
		}
	}
	
	/**
	 * creates a dataset of snps imputed at a given interval, e.g. every 0.2cM
	 * this uses files ending in .het.txt
	 * @param chr	the chromosome for which snps are to be generated
	 * @param outputFileName	the name of the output file of snps
	 * @param useArrayData	if true, uses the NAM 1106 array based markers, otherwise NAM gbs data
	 */
	public static void generateSNPdataset(int chr, String outputFileName, boolean useArrayData) {
		generateSNPdataset(chr, outputFileName, useArrayData, true);
	}
	
	/**
	 * @param chr	creates a dataset of snps imputed at a given interval, e.g. every 0.2cM
	 * @param outputFileName	the name of the output file of snps
	 * @param useArrayData	if true, uses the NAM 1106 array based markers, otherwise NAM gbs data
	 * @param useImputedHets if true, uses .het.txt files, otherwise uses .txt files which do not have hets imputed
	 */
	public static void generateSNPdataset(int chr, String outputFileName, boolean useArrayData, boolean useImputedHets) {
		
		HashMap<String,Integer> taxaMap = new HashMap<String,Integer>();
		HashMap<SNP, Integer> snpMap = new HashMap<SNP, Integer>();
		String input;
		String[] info;
		Pattern tab = Pattern.compile("\t");
		int numberOfTaxa;
		int numberOfSnps;
		ArrayList<String> taxaList = new ArrayList<String>();
		AGPMap agpmap = new AGPMap(true);

		//read in taxaList
		try {
			BufferedReader br = new BufferedReader(new FileReader("/Volumes/Macintosh HD 2/data/namgbs/namgbstaxalist.txt"));
			while ((input = br.readLine()) != null) {
				taxaList.add(input);
			}
			br.close();
		} catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		numberOfTaxa = taxaList.size();
		System.out.println("Number of taxa with genotypes: " + numberOfTaxa);
		
		//get a list of snp files
		final boolean usehets = useImputedHets;
		File snpdir = new File("/Volumes/Macintosh HD 2/data/namgbs/cov20/");
		final String prefix = "cov20_chr" + chr + "_Z";
		File[] snpfiles = snpdir.listFiles(new FilenameFilter() {

			@Override
			public boolean accept(File dir, String filename) {
				if (filename.startsWith(prefix)) {
					if (usehets) {
						if (filename.endsWith(".hets12.txt")) return true;
						else return false;
					} else {
						if (filename.contains("hets")) return false;
						else return true;
					}
				}
				return false;
			}
		});
		
		System.out.println("\nReading snp files.");

		ArrayList<SNP> snpList;
		if (useArrayData) {
			snpList = getNamArraySnps(chr, taxaList, agpmap);
		} else {
			snpList = getNamgbsSnps(taxaList, snpfiles);
		} 
		
		numberOfSnps = snpList.size();
		
		System.out.println("Number of snps = " + numberOfSnps);
//		for (int i = 0; i < 10; i++) System.gc();
//		total = Runtime.getRuntime().totalMemory();
//		free = Runtime.getRuntime().freeMemory();
//		used = total - free;
//		System.out.println("snplist memory: " + "total = " + total + ", free = " + free + ", used = " + used);
		
		//find the start and end of each chromosome rounded to 0.2 cM
		double startgenpos = agpmap.getCmFromPosition(chr, snpList.get(0).physicalPos);
		double endgenpos = agpmap.getCmFromPosition(chr, snpList.get(snpList.size() - 1).physicalPos);
		
		//create a list of imputed snps at 0.2 cM between these limits (inclusive)
		ArrayList<SNP> imputedSnpList = new ArrayList<SNP>();
		double interval = 0.2;
		
		double start = Math.ceil(startgenpos / interval) * interval;
		double end = Math.floor(endgenpos / interval) * interval;
		String snpPrefix;
		if (chr < 10) {
			snpPrefix = "P0" + chr + "_";
		} else {
			snpPrefix = "P10_";
		}
		
		int count = 0;
		for (double p = start; p <= end; p += 0.2) {
			int pos = agpmap.getPositionFromCm(chr, p);
			String snpname = snpPrefix + pos;
			imputedSnpList.add(new SNP(snpname, "A/B", chr, (float) p, pos, null, count++));
		}
		
		int numberOfImputedSnps = imputedSnpList.size();
		
//		for (int i = 0; i < 10; i++) System.gc();
//		total = Runtime.getRuntime().totalMemory();
//		free = Runtime.getRuntime().freeMemory();
//		used = total - free;
//		System.out.println("after loading all snps memory: " + "total = " + total + ", free = " + free + ", used = " + used);

		//impute snps
		System.out.println("Creating set of imputed snps.");
		
		float[][] imputedSnps = new float[numberOfImputedSnps][numberOfTaxa];
		for (int t = 0; t < numberOfTaxa; t++) {
			for (int s = 0; s < numberOfImputedSnps; s++) {
				SNP snp = imputedSnpList.get(s);
				int snpIndex = Collections.binarySearch(snpList, snp, new Comparator<SNP>(){

					@Override
					public int compare(SNP snp1, SNP snp2) {
						if (snp1.physicalPos > snp2.physicalPos) return 1;
						if (snp1.physicalPos < snp2.physicalPos) return -1;
						return 0;
					}
					
				});
				int endndx = snpList.size() - 1;
				int rightIndex = snpIndex;
				int leftIndex = snpIndex;
				if (snpIndex < 0) {
					rightIndex = -(snpIndex + 1);
					leftIndex = rightIndex - 1;
					if (leftIndex < 0) leftIndex = 0;
					if (rightIndex > endndx) rightIndex = endndx;
				}
				
				
				while (snpList.get(rightIndex).getFloatScore(t) == 3 && rightIndex < endndx) {
					rightIndex++;
				}
				while (snpList.get(leftIndex).getFloatScore(t) == 3 && leftIndex > 0) {
					leftIndex--;
				}
				
				SNP leftsnp = snpList.get(leftIndex);
				SNP rightsnp = snpList.get(rightIndex);
				float leftscore = leftsnp.getFloatScore(t);
				float rightscore = rightsnp.getFloatScore(t);
				if (rightscore == 3 || Math.abs(leftscore - rightscore) < 1e-10) imputedSnps[s][t] = leftsnp.getFloatScore(t);
				else if (leftsnp.getFloatScore(t) == 3) imputedSnps[s][t] = rightsnp.getFloatScore(t);
				else {
					float pd = ((float)(snp.physicalPos - leftsnp.physicalPos)) / ((float)(rightsnp.physicalPos - leftsnp.physicalPos));
					imputedSnps[s][t] = leftsnp.getFloatScore(t) * (1 - pd) + rightsnp.getFloatScore(t) * pd;
				}

			}
		}
		
		//write all imputed snp data to a single file
		System.out.println("Writing output to " + outputFileName);
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outputFileName));
			StringBuilder sb = new StringBuilder("Snp\tallele\tchr\tpos\tcm");
			for (int t = 0; t < numberOfTaxa; t++) {
				sb.append("\t").append(taxaList.get(t));
			}
			bw.write(sb.toString());
			bw.newLine();
			
			for (int s = 0; s < numberOfImputedSnps; s++) {
				SNP snp = imputedSnpList.get(s);
				sb = new StringBuilder(snp.name);
				sb.append("\t").append(snp.allele);
				sb.append("\t").append(snp.chr);
				sb.append("\t").append(snp.physicalPos);
				sb.append("\t").append(snp.geneticPos);
				for (int t = 0; t < numberOfTaxa; t++) {
					sb.append("\t");
					sb.append(imputedSnps[s][t]);
				}
				bw.write(sb.toString());
				bw.newLine();
			}
			bw.close();
		} catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		System.out.println("Finished.");
	}
	
	/**
	 * For the a chromosome, writes all SNP data as A or B for all populations, no hets, no missing data imputed
	 * Taxa with problems will be excluded (uses namgbstaxalist)
	 * @param chr	the chromosome
	 * @param outputFileName	snp data is written to this file, the file will be overwritten if it already exists
	 */
	public static void generateSNPdataAllSNPs(int chr, String outputFileName) {
		HashMap<String,Integer> taxaMap = new HashMap<String,Integer>();
		HashMap<SNP, Integer> snpMap = new HashMap<SNP, Integer>();
		String input;
		String[] info;
		Pattern tab = Pattern.compile("\t");
		int numberOfTaxa;
		int numberOfSnps;
		ArrayList<String> taxaList = new ArrayList<String>();
		AGPMap agpmap = new AGPMap(true);

		//read in taxaList
		try {
			BufferedReader br = new BufferedReader(new FileReader("/Volumes/Macintosh HD 2/data/namgbs/namgbstaxalist.txt"));
			while ((input = br.readLine()) != null) {
				taxaList.add(input);
			}
			br.close();
		} catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		numberOfTaxa = taxaList.size();
		System.out.println("Number of taxa with genotypes: " + numberOfTaxa);
		
		//get a list of snp files
		final boolean usehets = false;
		File snpdir = new File("/Volumes/Macintosh HD 2/data/namgbs/cov20");
		final String prefix = "cov20_chr" + chr + "_Z";
		File[] snpfiles = snpdir.listFiles(new FilenameFilter() {

			@Override
			public boolean accept(File dir, String filename) {
				if (filename.startsWith(prefix)) {
					if (usehets) {
						if (filename.endsWith(".hets20.txt")) return true;
						else return false;
					} else {
						if (filename.endsWith(".hets20.txt")) return false;
						else return true;
					}
				}
				return false;
			}
		});
		
		System.out.println("\nReading snp files.");

		ArrayList<SNP> snpList = getNamgbsSnps(taxaList, snpfiles);
		
		numberOfSnps = snpList.size();
		
		System.out.println("Number of snps = " + numberOfSnps);
		
		//write snps to the output file
		System.out.println("Writing SNPs.");
		
		//write all imputed snp data to a single file
		System.out.println("Writing output to " + outputFileName);
		char[] byte2char = new char[]{'A','H','B','N'};
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outputFileName));
			StringBuilder sb = new StringBuilder("Snp\tallele\tchr\tpos\tcm");
			for (int t = 0; t < numberOfTaxa; t++) {
				sb.append("\t").append(taxaList.get(t));
			}
			bw.write(sb.toString());
			bw.newLine();
			
			for (int s = 0; s < numberOfSnps; s++) {
				SNP snp = snpList.get(s);
				sb = new StringBuilder(snp.name);
				sb.append("\t").append(snp.allele);
				sb.append("\t").append(snp.chr);
				sb.append("\t").append(snp.physicalPos);
				sb.append("\t").append(snp.geneticPos);
				for (int t = 0; t < numberOfTaxa; t++) {
					sb.append("\t");
					sb.append(byte2char[snp.bytescore[t]]);
				}
				bw.write(sb.toString());
				bw.newLine();
			}
			bw.close();
		} catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		System.out.println("Finished.");
	}
	
	public static ArrayList<SNP> getNamgbsSnps(ArrayList<String> taxaList, File[] snpfiles) {
		HashMap<String,Byte> snp2byte = new HashMap<String,Byte>();
		snp2byte.put("A", (byte) 0);
		snp2byte.put("H", (byte) 1);
		snp2byte.put("B", (byte) 2);
		snp2byte.put("N", (byte) 3);
		String input;
		String[] info;
		Pattern tab = Pattern.compile("\t");
		int numberOfTaxa = taxaList.size();
		
		//read the snpfiles once to create a set of snps and a set of taxa
		HashMap<SNP, byte[]> snpMap = new HashMap<SNP, byte[]>();
		TreeSet<String> taxaSet = new TreeSet<String>();
		
		for (File file:snpfiles) {

			try{
				BufferedReader br = new BufferedReader(new FileReader(file));
				info = tab.split(br.readLine());
				int ntaxa = info.length - 11;
				String[] taxanames = new String[ntaxa];
				for (int t = 0; t < ntaxa; t++) {
					taxanames[t] = info[t + 11].substring(0, info[t+11].indexOf(':'));
					taxaSet.add(taxanames[t]);
				}
				
				while ((input = br.readLine()) != null) {
					info = tab.split(input);
					SNP snp = new SNP(info[0], info[1], Integer.parseInt(info[2]), 0, Integer.parseInt(info[3]), null, 0);
					byte[] geno = snpMap.get(snp);
					if (geno == null) {
						geno = new byte[numberOfTaxa];
						for (int t = 0; t < numberOfTaxa; t++) geno[t] = 3;
						snpMap.put(snp, geno);
					}
					
					//update geno with this snp information
					for (int t = 0; t < ntaxa; t++) {
						int ndxTaxon = Collections.binarySearch(taxaList, taxanames[t]);
						if (ndxTaxon > -1) {
							geno[ndxTaxon] = snp2byte.get(info[t + 11]);
						}
					}					
				}
				br.close();
			} catch(IOException e) {
				e.printStackTrace();
				System.exit(-1);
			}
		}
		
		//load the genotypes into the snps
		ArrayList<SNP> snpList = new ArrayList<SNP>(snpMap.keySet());
		Collections.sort(snpList);
		for (SNP snp:snpList) snp.bytescore = snpMap.get(snp);
		return snpList;
		
	}
	
	public static ArrayList<SNP> getNamArraySnps(int chr, ArrayList<String> taxaList, AGPMap agpmap) {
		File namdataFile = new File("/Volumes/Macintosh HD 2/data/namgbs/ImputedMarkerGenotypes_flowering_traits_092909.txt");
		ArrayList<SNP> snpList = new ArrayList<SNP>();
		int nmarkers = agpmap.marker.length;
		int markerNumber = 0;
		while (agpmap.markerChromosome[markerNumber] < chr) markerNumber++;
		
		int count = 0;
		int numberOfTaxa = taxaList.size();
		int firstMarkerColumn = markerNumber + 5;
		
		while (markerNumber < nmarkers && agpmap.markerChromosome[markerNumber] == chr) {
			int physicalPos = agpmap.getPositionFromCm(chr, agpmap.markercm[markerNumber]);
			float[] score = new float[numberOfTaxa]; 
			snpList.add(new SNP(agpmap.marker[markerNumber], "", chr, (float) agpmap.markercm[markerNumber], physicalPos, score, count++));
			markerNumber++;
		}
		
		int numberOfSnps = snpList.size();
		String input;
		String[] info;
		Pattern tab = Pattern.compile("\t");
		try {
			BufferedReader br = new BufferedReader(new FileReader(namdataFile));
			br.readLine();
			while((input = br.readLine()) != null) {
				info = tab.split(input);
				int ndx = taxaList.indexOf(info[0]);
				if (ndx >= 0) {
					for (SNP snp:snpList) {
						snp.score[ndx] = Float.parseFloat(info[firstMarkerColumn + snp.index]);
					}
				}
			}
			
			br.close();
		} catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
		return snpList;
	}
	
	public static void createImputedSNPDatasets() {
		//create imputed snps from nam array
		for (int c = 1; c <= 10; c++) {
			String filename = "/Volumes/Macintosh HD 2/data/namgbs/cov20/imputedFromArraySnpsChr" + c + ".txt";
			generateSNPdataset(c, filename, true);
		}

		//create imputed snps without imputing hets first
//		for (int c = 1; c <= 10; c++) {
//			String filename = "/Volumes/Macintosh HD 2/data/namgbs/cov20/imputedWithoutHetsSnpsChr" + c + ".txt";
//			generateSNPdataset(c, filename, false, false);
//		}
		
		//create imputed snps from imputed hets
//		for (int c = 1; c <= 10; c++) {
//			String filename = "/Volumes/Macintosh HD 2/data/namgbs/cov20/imputedFromHets12SnpsChr" + c + ".txt";
//			generateSNPdataset(c, filename, false, true);
//		}
		
		
	}
	
	public static void concatenateSomeFiles() {
		File[] input = new File[10];
		
		for (int chr = 1; chr <= 10; chr++) {
//			input[chr - 1] = new File("C:/Projects/NAM/namgbs/datasets/cov20/imputedSnpsChr" + chr + ".txt");
//			input[chr - 1] = new File("C:/Projects/NAM/namgbs/datasets/cov20/imputedFromArraySnpsChr" + chr + ".txt");
//			input[chr - 1] = new File("C:/Projects/NAM/namgbs/datasets/cov20/imputedWithoutHetsSnpsChr" + chr + ".txt");
//			input[chr - 1] = new File("C:/Projects/NAM/namgbs/datasets/cov20/imputedFromHetsSnpsChr" + chr + ".txt");
//			input[chr - 1] = new File("/Volumes/Macintosh HD 2/data/namgbs/cov20/imputedFromHets12SnpsChr" + chr + ".txt");
			input[chr - 1] = new File("/Volumes/Macintosh HD 2/data/namgbs/cov20/imputedFromArraySnpsChr" + chr + ".txt");
		}
//		File output = new File("C:/Projects/NAM/namgbs/datasets/cov20/imputedSnpsAllChr.txt");
//		File output = new File("C:/Projects/NAM/namgbs/datasets/cov20/imputedFromArraySnpsAllChr.txt");
//		File output = new File("C:/Projects/NAM/namgbs/datasets/cov20/imputedWithoutHetsSnpsAllChr.txt");
//		File output = new File("C:/Projects/NAM/namgbs/datasets/cov20/imputedFromHetsSnpsAllChr.txt");
//		File output = new File("/Volumes/Macintosh HD 2/data/namgbs/cov20/imputedFromHets12SnpsAllChr.txt");
		File output = new File("/Volumes/Macintosh HD 2/data/namgbs/cov20/imputedFromArraySnpsAllChr.txt");
		concatenateOutput(input, output);
	}
	
	public static void concatenateOutput(File[] inputFiles, File outputFile) {
		//assumes that the first line consists of column headers and that the headers are the same for all files
		int numberOfInputFiles = inputFiles.length;
		byte lf = (byte) '\n';
		
		try {
			FileOutputStream fos = new FileOutputStream(outputFile);
			FileChannel outfc = fos.getChannel();
			
			System.out.println("Processing " + inputFiles[0].getPath());
			FileInputStream fis = new FileInputStream(inputFiles[0]);
			FileChannel infc = fis.getChannel();
			int bufferCapacity = 100000;
			ByteBuffer bb = ByteBuffer.allocate(bufferCapacity);
			bb.clear();
			
			while (infc.read(bb) > 0) {
				bb.flip();
				outfc.write(bb);
				bb.clear();
			}
			infc.close();
			
			for (int f = 1; f < numberOfInputFiles; f++) {
				System.out.println("Processing " + inputFiles[f].getPath());
				fis = new FileInputStream(inputFiles[f]);
				infc = fis.getChannel();
				bb.clear();
				int bytesread = infc.read(bb);
				bb.flip();
				byte b = bb.get();
				while(b != lf) {
					b = bb.get();
				}
				outfc.write(bb);
				bb.clear();
				while (infc.read(bb) > 0) {
					bb.flip();
					outfc.write(bb);
					bb.clear();
				}
				infc.close();
			}

			outfc.close();
			
		} catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	public static void createListOfTaxa() {
		HashMap<String,Integer> taxaMap = new HashMap<String,Integer>();
		HashMap<SNP, Integer> snpMap = new HashMap<SNP, Integer>();
		String input;
		String[] info;
		Pattern tab = Pattern.compile("\t");
		int numberOfTaxa;
		int numberOfSnps;
		ArrayList<String> taxaList = new ArrayList<String>();
		AGPMap agpmap = new AGPMap(true);

		//read in taxanames and delete the poor quality lines
		String taxaNameFile = "C:/Projects/NAM/NAM_map_and_genos-090921/RILs_for_NAM_Map_20071102.txt";
		try {
			BufferedReader br = new BufferedReader(new FileReader(taxaNameFile));
			info = tab.split(br.readLine());
			for (String t:info) taxaList.add(t);
			br.close();
			
			System.out.println("Taxa names from NAM list: " + taxaList.size());
			String qualityFile = "C:/Projects/NAM/namgbs/datasets/cov20/sample quality issues.txt";
			br = new BufferedReader(new FileReader(qualityFile));
			br.readLine();
			while ((input = br.readLine()) != null) {
				info = tab.split(input);
				taxaList.remove(info[0]);
			}
			System.out.println("Taxa names after removing samples with problems: " + taxaList.size());
			
			br.close();
		} catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
		//sort the taxa
		Collections.sort(taxaList);
		
		for (int i = 0; i < 10; i++) System.gc();
		long total = Runtime.getRuntime().totalMemory();
		long free = Runtime.getRuntime().freeMemory();
		long used = total - free;
		System.out.println("taxalist memory: " + "total = " + total + ", free = " + free + ", used = " + used);
		
		//get a list of snp files
		final boolean usehets = false;
		File snpdir = new File("C:/Projects/NAM/namgbs/datasets/cov20/");
		final String prefix = "cov20_chr1_Z";
		File[] snpfiles = snpdir.listFiles(new FilenameFilter() {

			@Override
			public boolean accept(File dir, String filename) {
				if (filename.startsWith(prefix)) {
					if (usehets) {
						if (filename.endsWith(".hets.txt")) return true;
						else return false;
					} else {
						if (filename.endsWith(".hets.txt")) return false;
						else return true;
					}
				}
				return false;
			}
		});
		
		//read taxa from the snp files and delete any taxa without genotypes from the taxalist
		int n = taxaList.size();
		boolean[] taxaHaveGenotypes = new boolean[n];
		for (int i = 0 ; i < n; i++) taxaHaveGenotypes[i] = false;
		for (File file:snpfiles) {
			try {
				BufferedReader br = new BufferedReader(new FileReader(file));
				info = tab.split(br.readLine());
				int ntaxa = info.length - 11;
				for (int t = 0; t < ntaxa; t++) {
					String taxon = info[t + 11].substring(0, info[t+11].indexOf(':'));
					int ndx = taxaList.indexOf(taxon);
					if (ndx >= 0) taxaHaveGenotypes[ndx] = true;
				}
				
			} catch(IOException e) {
				e.printStackTrace();
				System.exit(-1);
			}
		}
		for (int i = n - 1; i >= 0 ; i--) {
			if (!taxaHaveGenotypes[i]) taxaList.remove(i);
		}
		
		numberOfTaxa = taxaList.size();
		
		//write taxa to file for comparison with nam array data
		String outFile = "c:/users/peter/temp/namgbstaxalist.txt";
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));
			for (String taxon:taxaList) bw.write(taxon + "\n");
			bw.close();
		} catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
		System.out.println("Taxa written to " + outFile);
	}
	
	public void imputeUsingAllData(int chromosome, int pop) {
		
		class SnpInfo {
			String name;
			boolean isgbs;
			int chr;
			int pos;
			int ndx;
			byte val;
		}

		String gbsGenotypeFilename = "/Volumes/Macintosh HD 2/data/namgbs/genos_20111003/NAM_20111003_ABHv2_chr" + chromosome + ".txt";
		String arrayGenotypeFilename = "/Volumes/Macintosh HD 2/data/namgbs/nam1106/ICIMdataForPopulation" + pop + ".txt";
		String arrayPositionsFilename = "/Volumes/Macintosh HD 2/data/namgbs/markers061208agpv2.txt";
		String taxalistFilename = "/Volumes/Macintosh HD 2/data/namgbs/namgbstaxalist.txt";
		Pattern tab = Pattern.compile("\t");
		
		HashMap<String, int[]> arrayMap = new HashMap<String, int[]>();
		
		//read the taxalist
		String znum;
		if (pop < 10) znum = "Z00" + pop;
		else znum = "Z0" + pop;
		System.out.println("Reading taxa list and array positions...");
		LinkedList<String> taxaList = new LinkedList<String>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(taxalistFilename));
			String input;
			while ((input = br.readLine()) != null) {
				if (input.startsWith(znum)) taxaList.add(input);
  			}
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
		
		//read array marker physical positions
		try {
			BufferedReader br = new BufferedReader(new FileReader(arrayPositionsFilename));
			String input = br.readLine();
			while ((input = br.readLine()) != null) {
				String[] info = tab.split(input);
				if (info.length > 6 && info[6].equals("estimated")) {
					//do not use this marker
				} else {
					int chr = Integer.parseInt(info[1]);
					int pos = Integer.parseInt(info[5]);
					arrayMap.put(info[0], new int[]{chr,pos});
				}
  			}
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
		//read gbs data
		System.out.println("Reading gbs data...");
		byte[][] gbsSnps = new byte[0][0];
		HashMap<String, Integer> gbsTaxaMap = new HashMap<String, Integer>();
		LinkedList<SnpInfo> gbsSnpList = new LinkedList<SnpInfo>();

		try {
			//count rows and columns
			BufferedReader br = new BufferedReader(new FileReader(gbsGenotypeFilename));
			String input = br.readLine();
			String[] info = tab.split(input);
			int ntaxa = info.length - 11;
			int nsnps = 0;
			while ((input = br.readLine()) != null) nsnps++;
			br.close();
			
			//read data into array
			gbsSnps = new byte[ntaxa][nsnps];
			br = new BufferedReader(new FileReader(gbsGenotypeFilename));
			input = br.readLine();
			info = tab.split(input);
			for (int t = 0; t < ntaxa; t++) {
				String name = info[t + 11];
				name = name.substring(0, name.indexOf(':'));
				gbsTaxaMap.put(name, t);
			}
			int count = 0;
			while ((input = br.readLine()) != null) {
				info = tab.split(input);
				SnpInfo snp = new SnpInfo();
				snp.name = info[0];
				snp.chr = Integer.parseInt(info[2]);
				snp.pos = Integer.parseInt(info[3]);
				snp.ndx = count;
				snp.isgbs = true;
				gbsSnpList.add(snp);
				for (int t = 0; t < ntaxa; t++) {
					gbsSnps[t][count] = snpToByteMap.get(info[t + 11]);
				}
				count++;
			}
			br.close();

		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
		//read array data
		System.out.println("Reading array data...");

		HashMap<String, Integer> arrayTaxaMap = new HashMap<String, Integer>();
		LinkedList<SnpInfo> arraySnpList = new LinkedList<SnpInfo>();
		ArrayList<byte[]> arraySnps = new ArrayList<byte[]>();
		
		try {
			//count rows and columns
			BufferedReader br = new BufferedReader(new FileReader(arrayGenotypeFilename));
			String input = br.readLine();
			String[] info = tab.split(input);
			int ntaxa = info.length;
			int nsnps = 0;
			while ((input = br.readLine()) != null) nsnps++;
			br.close();
			
			//read data into array
			br = new BufferedReader(new FileReader(arrayGenotypeFilename));
			input = br.readLine();
			String[] arrayTaxa = tab.split(input);
			int tcount = 0;
			for (String taxon : arrayTaxa) arrayTaxaMap.put(taxon, tcount++);
			
			int count = 0;
			while ((input = br.readLine()) != null) {
				info = tab.split(input);
				String name = info[0];
				int[] posinfo = arrayMap.get(name);
				if (posinfo != null && posinfo[0] == chromosome) {
					SnpInfo snp = new SnpInfo();
					snp.name = name;
					snp.ndx = count++;
					snp.chr = posinfo[0];
					snp.pos = posinfo[1];
					snp.isgbs = false;
					arraySnpList.add(snp);
					byte[] snps = new byte[ntaxa];
					for (int t = 0; t < ntaxa; t++) {
						snps[t] = Byte.parseByte(info[t + 1]);
					}
					arraySnps.add(snps);
				}
				
			}
			br.close();

		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
		//create a combined snp list sorted in order
		ArrayList<SnpInfo> snpList = new ArrayList<SnpInfo>();
		snpList.addAll(gbsSnpList);
//		snpList.addAll(arraySnpList);
		
		//sort
		Collections.sort(snpList, new Comparator<SnpInfo>() {

			@Override
			public int compare(SnpInfo s1, SnpInfo s2) {
				return s1.pos - s2.pos;
			}
		});
		
		//the imputed snps array
		int ntaxa = taxaList.size();
		numberOfSnps = snpList.size();
		byte[][] imputedSnps = new byte[ntaxa][];
		
		//for each taxon
		//create a vector of non-missing snp calls and an index into snpList for those calls
		//impute using the Viterbi algorithm
		//write to an imputed snps array
		
		//counts to keep track of agreement between observations and imputed states. First dimension is observation, second is state.
		int[][] NAMmatchCounts;
		int[][] GBSmatchCounts;

		//probabilities
		
		System.out.println("Imputing snps...");
		double[][] gbsProbObsGivenState = new double[][]{{0.995,0.2,0.002},{0.003,.6,0.003},{0.002,0.2,0.995}};
		double[][] arrayProbObsGivenState = new double[][]{{0.9998,0.0001,0.0001},{0.0001,0.9998,0.0001},{0.0001,0.0001,0.9998}};
		ProbabilityObservationGivenTrueState pObsGivenTrue= new ProbabilityObservationGivenTrueState();
		pObsGivenTrue.setGbsProbabilityMatrix(gbsProbObsGivenState);
		pObsGivenTrue.setArrayProbabilityMatrix(arrayProbObsGivenState);
		
		double[] pstate = new double[]{.485, .03, .485};
		Recombination recomb = new Recombination();
		recomb.readRates("/Volumes/Macintosh HD 2/data/namgbs/cov20/cov20_chr_chr" + chromosome + ".recombinationrate.het12.win1.txt");

		for (int iter = 0; iter < 8; iter++) {
			NAMmatchCounts = new int[3][3];
			GBSmatchCounts = new int[3][3];
			int[][] transitionCounts = new int[3][3]; //first dimension is state1, second is state2
			int taxaCount = 0;
			TransitionMatrix tm = new TestTransitionMatrix(.0003);
			int nNotMissing = 0;
			for (String taxon:taxaList) {
				for (int j = 0; j < gbsSnps[0].length; j++) {
					if (gbsSnps[gbsTaxaMap.get(taxon)][j] != N) nNotMissing++;
				}
			}
			System.out.println("non missing count = " + nNotMissing);
			
			for (String taxon:taxaList) {
				LinkedList<SnpInfo> nonMissing = new LinkedList<SnpInfo>();
				Integer taxonIndex = gbsTaxaMap.get(taxon);
				if (taxonIndex == null) break;
				byte[] gbsCalls = gbsSnps[gbsTaxaMap.get(taxon)];
				int arrayTaxonIndex = arrayTaxaMap.get(taxon);
//				int snpCount = 0;
//				for (SnpInfo snp:snpList) {
				for (int k = 0; k < numberOfSnps; k++) {
					SnpInfo snp = snpList.get(k);
					byte val;
					if (snp.isgbs) {
						val = gbsCalls[snp.ndx];
						if (val == N) val = -1;
						else if (val == A) val = 0;
						else if (val == B) val = 2;
						else val = 1;
					} else {
						val = arraySnps.get(snp.ndx)[arrayTaxonIndex];
					}
					if (val > -1) {
						SnpInfo s = new SnpInfo();
						s.pos = snp.pos;
						s.name = snp.name;
						s.isgbs = snp.isgbs;
						s.val = val;
//						s.ndx = snpCount;
						s.ndx = k;
						nonMissing.add(s);
					}
//					snpCount++;
				}

				int n = nonMissing.size();
				byte[] geno = new byte[n];
				boolean[] isgbs = new boolean[n];
				int[] snpPosition = new int[n];

				int snpCount = 0;
				for (SnpInfo snp:nonMissing) {
					geno[snpCount] = snp.val;
					isgbs[snpCount] = snp.isgbs;
					snpPosition[snpCount++] = snp.pos;
				}

//				if (taxon.equals("Z001E0032")) {
//					System.out.println("Z001E0032");
//				}

				recomb.setEvaluationPositions(snpPosition);
				pObsGivenTrue.setIsGBS(isgbs);
				ViterbiAlgorithm  va = new ViterbiAlgorithm(geno, recomb, pObsGivenTrue, pstate);
//				ViterbiAlgorithm  va = new ViterbiAlgorithm(geno, tm, pObsGivenTrue, pstate);
				va.calculate();
				byte[] genostates = va.getMostProbableStateSequence();
				
				numberOfSnps = snpList.size();
				byte[] imputedsnp = new byte[numberOfSnps];
				
				//count changes
				snpCount = 0;
				int[][] matches = new int[3][3];
				for (SnpInfo snp:nonMissing) {
					int state = genostates[snpCount];
					int obs = geno[snpCount];
					if (snp.isgbs) {
						GBSmatchCounts[obs][state]++;
					} else {
						NAMmatchCounts[obs][state]++;
					}
					matches[obs][state]++;
					snpCount++;
				}
//				System.out.println(taxon + ": " + matches[0][1] + "\t" + matches[0][2] + "\t" + matches[1][0] + "\t" + matches[1][2] + "\t" + matches[2][0] + "\t" + matches[2][1]+ "\t" + matches[0][0] + "\t" + matches[1][1] + "\t" + matches[2][2]);
				n = genostates.length;
				for (int i = 1; i < n; i++) {
					transitionCounts[genostates[i-1]][genostates[i]]++;
				}
				
				//initialize imputedsnps
				for (int s = 0; s < numberOfSnps; s++) {
					imputedsnp[s] = -1;
				}

				snpCount = 0;
				for (SnpInfo s:nonMissing) {
					imputedsnp[s.ndx] = genostates[snpCount];
					snpCount++;
				}

				imputedSnps[taxaCount] = imputedsnp;
				taxaCount++;
			}

			int nm = 0;
			for (int i = 0; i < imputedSnps.length; i++) {
				for (int j = 0; j < imputedSnps[0].length; j++) {
					if (imputedSnps[i][j] != -1) nm++;
				}
			}
			System.out.println("nonmissing in imputedSnps = " + nm);
			
			//summarize parameters
			System.out.println("GBS counts");
			System.out.println(GBSmatchCounts[0][0] + "\t" + GBSmatchCounts[0][1] + "\t" + GBSmatchCounts[0][2]);
			System.out.println(GBSmatchCounts[1][0] + "\t" + GBSmatchCounts[1][1] + "\t" + GBSmatchCounts[1][2]);
			System.out.println(GBSmatchCounts[2][0] + "\t" + GBSmatchCounts[2][1] + "\t" + GBSmatchCounts[2][2]);
			System.out.println("NAM counts");
			System.out.println(NAMmatchCounts[0][0] + "\t" + NAMmatchCounts[0][1] + "\t" + NAMmatchCounts[0][2]);
			System.out.println(NAMmatchCounts[1][0] + "\t" + NAMmatchCounts[1][1] + "\t" + NAMmatchCounts[1][2]);
			System.out.println(NAMmatchCounts[2][0] + "\t" + NAMmatchCounts[2][1] + "\t" + NAMmatchCounts[2][2]);
			System.out.println();
			
			int[] observed = new int[3];
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					observed[i] += GBSmatchCounts[i][j]; 
				}
			}

//			for (int i = 0; i < 3; i++) {
//				for (int j = 0; j < 3; j++) {
//					trueGivenObserved[i][j] = ((double) GBSmatchCounts[j][i]) / ((double) observed[j]); 
//				}
//			}
			
			//calculate p(state)
			int[] nstate = new int[3];
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					nstate[j] += GBSmatchCounts[i][j];
					nstate[j] += NAMmatchCounts[i][j];
				}
			}
			
			double total = nstate[0] + nstate[1] + nstate[2];
			for (int i = 0; i < 3; i++) {
				pstate[i] = ((double) nstate[i]) / total;
			}
			
			double[] gbsColumnTotal = new double[]{0,0,0};
			double[] arrayColumnTotal = new double[]{0,0,0};
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					gbsColumnTotal[j] += GBSmatchCounts[i][j];
					arrayColumnTotal[j] += NAMmatchCounts[i][j];
				}
			}
			
			//calculate p(obs|state) for GBS and array, obs row, state column 
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					gbsProbObsGivenState[i][j] = GBSmatchCounts[i][j] / gbsColumnTotal[j];
					arrayProbObsGivenState[i][j] = NAMmatchCounts[i][j] / arrayColumnTotal[j];
				}
			}
			pObsGivenTrue.setArrayProbabilityMatrix(arrayProbObsGivenState);
			pObsGivenTrue.setGbsProbabilityMatrix(gbsProbObsGivenState);
			
			//set markers with excess hets to missing
			int nsnps = imputedSnps[0].length;
			int markercount = 0;
			for (int s = 0; s < nsnps; s++) {
				int nhets = 0;
				int ntotal = 0;
				for (int t = 0; t < ntaxa; t++) {
					if (imputedSnps[t][s] != -1) {
						ntotal++;
						if (imputedSnps[t][s] == 1) nhets++;
					}
				}
				double phet = ((double) nhets) / ((double) ntotal);
				if (phet > .5) {
					SnpInfo snpinfo = snpList.get(s);
//					System.out.println("ndx = " + snpinfo.ndx);
					if (snpinfo.isgbs) {
						for (String taxon:taxaList) {
//							System.out.print("t = " + t + "before = " + gbsSnps[t][snpinfo.ndx]);
							gbsSnps[gbsTaxaMap.get(taxon)][snpinfo.ndx] = N;
//							System.out.println("after = " + gbsSnps[t][snpinfo.ndx]);
						}
						markercount++;
					} else {
						byte[] snpdata = arraySnps.get(snpinfo.ndx);
						for (int t = 0; t < ntaxa; t++) {
							snpdata[t] = -1;
						}
					}
				}
			}
			System.out.println("markers set to missing = " + markercount);
			
			//calculate transition probabilities
			

		}
		
		//write the output
		numberOfTaxa = taxaList.size();
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter("/Volumes/Macintosh HD 2/data/namgbs/genos_20111003/genos_20111003_chr" + chromosome +"_" + znum + ".hmm_abhv2.txt"));
			bw.write("rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode");
			for (String taxon:taxaList) {
				bw.write("\t");
				bw.write(taxon);
			}
			bw.newLine();
			
			int snpCount = 0;
			for (SnpInfo snp:snpList) {
				bw.write(snp.name);
				bw.write("\tNA\t");
				bw.write(Integer.toString(snp.chr));
				bw.write("\t");
				bw.write(Integer.toString(snp.pos));
				bw.write("\tNA\tNA\tNA\tNA\tNA\tNA\tNA");
				for (int t = 0; t < numberOfTaxa; t++) {
					bw.write("\t");
					int val;
					try {
						val = imputedSnps[t][snpCount];
					} catch (Exception e) {
						val = -1;
					}
					bw.write(Integer.toString(val));
				}
				snpCount++;
				bw.newLine();
			}
			
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		//write nam array data
//		try {
//			BufferedWriter bw = new BufferedWriter(new FileWriter("/Volumes/Macintosh HD 2/data/namgbs/genos_20111003/genos_20111003_chr1_" + znum + ".array.txt"));
//			bw.write("rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode");
//			for (String taxon:taxaList) {
//				bw.write("\t");
//				bw.write(taxon);
//			}
//			bw.newLine();
//			
//			int snpCount = 0;
//			for (SnpInfo snp:arraySnpList) {
//				bw.write(snp.name);
//				bw.write("\tNA\t");
//				bw.write(Integer.toString(snp.chr));
//				bw.write("\t");
//				bw.write(Integer.toString(snp.pos));
//				bw.write("\tNA\tNA\tNA\tNA\tNA\tNA\tNA");
//				byte[] snpvals = arraySnps.get(snp.ndx);
//				for (String taxon:taxaList) {
//					bw.write("\t");
//					int val = snpvals[arrayTaxaMap.get(taxon)];
//					bw.write(Integer.toString(val));
//				}
//				snpCount++;
//				bw.newLine();
//			}
//			
//			bw.close();
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
		
		//write nam gbs data
//		try {
//			BufferedWriter bw = new BufferedWriter(new FileWriter("/Volumes/Macintosh HD 2/data/namgbs/genos_20111003/genos_20111003_chr1_Z001.gbs.txt"));
//			bw.write("rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode");
//			for (String taxon:taxaList) {
//				bw.write("\t");
//				bw.write(taxon);
//			}
//			bw.newLine();
//			
//			int snpCount = 0;
//			for (SnpInfo snp:gbsSnpList) {
//				bw.write(snp.name);
//				bw.write("\tNA\t");
//				bw.write(Integer.toString(snp.chr));
//				bw.write("\t");
//				bw.write(Integer.toString(snp.pos));
//				bw.write("\tNA\tNA\tNA\tNA\tNA\tNA\tNA");
//				for (String taxon:taxaList) {
//					int ndx = gbsTaxaMap.get(taxon);
//					bw.write("\t");
//					byte snpval = gbsSnps[ndx][snp.ndx];
//					int val;
//					if (snpval == A) val = 0;
//					else if (snpval == H) val = 1;
//					else if (snpval == B) val = 2;
//					else val = -1;
//					bw.write(Integer.toString(val));
//				}
//				snpCount++;
//				bw.newLine();
//			}
//			
//			bw.close();
//		} catch (IOException e) {
//			e.printStackTrace();
//		}


		
		//write nam gbs data
	}
	
	public static void countOfftypesByLine(int chr) {
		System.out.println("Counting offtypes by line.");
		
		Pattern tab = Pattern.compile("\t");
		String filename = "/Volumes/Macintosh HD 2/data/namgbs/genos_1217/NAM282_20111217_scv10mF8maf002_mgs_E1pLD5_chr" + chr + ".hmp.txt";
		ArrayList<String> taxanames = new ArrayList<String>();
		ArrayList<Integer> colList = new ArrayList<Integer>();
		ArrayList<byte[]> snps = new ArrayList<byte[]>();
		
		//read the data
		try{
			BufferedReader br = new BufferedReader(new FileReader(filename));
			String input = br.readLine();
			String[] info = tab.split(input);
			for (int i = 11; i < info.length; i++) {
				if (info[i].startsWith("Z0")) {
					taxanames.add(info[i]);
					colList.add(i);
				}
			}
			int ntaxa = taxanames.size();
			
			while ((input = br.readLine()) != null) {
				info = tab.split(input);
				int taxacount = 0;
				byte[] snp = new byte[ntaxa];
				snps.add(snp);
				for (Integer col: colList) {
					String strval = info[col];
					if (strval.equals("A")) snp[taxacount] = 0;
					else if (strval.equals("C")) snp[taxacount] = 1;
					else if (strval.equals("G")) snp[taxacount] = 2;
					else if (strval.equals("T")) snp[taxacount] = 3;
					else snp[taxacount] = 4;
					taxacount++;
				}
			}
			br.close();
			
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
		
		//generate popmap
		HashMap<String, ArrayList<Integer>> popmap = new HashMap<String, ArrayList<Integer>>();
		int ntaxa = taxanames.size();
		int taxacount = 0;
		for (String name:taxanames) {
			String popname = name.substring(0,4);
			ArrayList<Integer> taxaList = popmap.get(popname);
			if (taxaList == null) {
				taxaList = new ArrayList<Integer>();
				popmap.put(popname, taxaList);
			}
			taxaList.add(taxacount++);
		}
		
		//within a population, find the (nearly) monomorphic loci
		//count the number of offtype allele calls for each line at the nearly monomorphic loci
		int[] offtypesByTaxon = new int[ntaxa];
		for (String popname:popmap.keySet()) {
			ArrayList<Integer> taxonIndex = popmap.get(popname);
			for (byte[] snp :snps) {
				//count the alleles
				int[] counts = new int[5];
				for (Integer t:taxonIndex) {
					counts[snp[t]]++;
				}
				
				int major = 0;
				int minor = 1;
				for (int i = 1; i < 4; i++) {
					if (counts[i] > counts[major]) {
						minor = major;
						major = i;
					} else if (counts[i] > counts[minor]){
						minor = i;
					}
				}
				
				int mincount = counts[minor];
				int minallele = minor;
				int total = counts[major] + counts[minor];
				if (mincount > 0 && Probability.binomial(mincount, total, 0.5) < .000001) {
					for (Integer t:taxonIndex) {
						if (snp[t] == minallele) offtypesByTaxon[t]++;
					}
				}
			}
		}
		
		//write the results to a file
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter("/Volumes/Macintosh HD 2/data/namgbs/genos_1217/NAM282_20111217_scv10mF8maf002_mgs_E1pLD5_chr" + chr + "offtype.count.hmp.txt"));
			bw.write("taxon\tofftypes\n");
			for (int t = 0; t < ntaxa; t++) {
				bw.write(taxanames.get(t));
				bw.write("\t" + offtypesByTaxon[t] + "\n");
			}
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.out.println("File error, did not write output.");
		}
		System.out.println("Finished counting offtypes.");
	}
	
	public static void calculateObservedGenotypeFrequenciesGBS() {
		File[] genotypefiles;
		int[] genoCounts = new int[3]; //A = 0, H = 1, B = 2
		File dir = new File("/Volumes/Macintosh HD 2/data/namgbs/genos_20111003");
		Pattern tab = Pattern.compile("\t");
		
		genotypefiles = dir.listFiles(new FilenameFilter(){

			@Override
			public boolean accept(File dir, String filename) {
				if (filename.startsWith("genos_20111003__chr") && filename.endsWith(".txt")) return true;
				return false;
			}
		});
		
		for (File file:genotypefiles) {
			System.out.println("Reading " + file.getName() + " ...");
			try {
				BufferedReader br = new BufferedReader(new FileReader(file));
				String input = br.readLine();
				while ((input = br.readLine()) != null) {
					String[] info = tab.split(input);
					int n = info.length;
					for (int i = 11; i < n; i++) {
						String val = info[i];
						if (val.equals("A")) genoCounts[0]++;
						else if (val.equals("H")) genoCounts[1]++;
						else if (val.equals("B")) genoCounts[2]++;
					}
				}
				br.close();
			} catch (IOException e) {
				e.printStackTrace();
				System.exit(-1);
			}
		}
		
		
		System.out.println("Counts: A = " + genoCounts[0] + ", H = " + genoCounts[1] + ", B = " + genoCounts[2]);
		double sumCounts = genoCounts[0] + genoCounts[1] + genoCounts[2];
		double freqA = ((double) genoCounts[0]) / sumCounts;
		double freqH = ((double) genoCounts[1]) / sumCounts;
		double freqB = ((double) genoCounts[2]) / sumCounts;
		System.out.println("Frequencies: A = " + freqA + ", H = " + freqH + ", B = " + freqB);
		System.out.println();
	}
	
	public static void calculateObservedGenotypeFrequenciesNAM() {
		File[] genotypefiles;
		int[] genoCounts = new int[3]; //A = 0, H = 1, B = 2
		File dir = new File("/Volumes/Macintosh HD 2/data/namgbs/nam1106");
		Pattern tab = Pattern.compile("\t");
		
		genotypefiles = dir.listFiles(new FilenameFilter(){

			@Override
			public boolean accept(File dir, String filename) {
				if (filename.startsWith("ICIMdataForPopulation") && filename.endsWith(".txt")) return true;
				return false;
			}
		});
		
		for (File file:genotypefiles) {
			System.out.println("Reading " + file.getName() + " ...");
			try {
				BufferedReader br = new BufferedReader(new FileReader(file));
				String input = br.readLine();
				while ((input = br.readLine()) != null) {
					String[] info = tab.split(input);
					int n = info.length;
					for (int i = 1; i < n; i++) {
						String val = info[i];
						if (val.equals("0")) genoCounts[0]++;
						else if (val.equals("1")) genoCounts[1]++;
						else if (val.equals("2")) genoCounts[2]++;
					}
				}
				br.close();
			} catch (IOException e) {
				e.printStackTrace();
				System.exit(-1);
			}
		}
		
		
		System.out.println("Counts: A = " + genoCounts[0] + ", H = " + genoCounts[1] + ", B = " + genoCounts[2]);
		double sumCounts = genoCounts[0] + genoCounts[1] + genoCounts[2];
		double freqA = ((double) genoCounts[0]) / sumCounts;
		double freqH = ((double) genoCounts[1]) / sumCounts;
		double freqB = ((double) genoCounts[2]) / sumCounts;
		System.out.println("Frequencies: A = " + freqA + ", H = " + freqH + ", B = " + freqB);
		System.out.println();
	}

	public static void countTransitionTypesNAM() {
		File[] genotypefiles;
		int[][] transitionCounts = new int[3][3];
		File dir = new File("/Volumes/Macintosh HD 2/data/namgbs/nam1106");
		Pattern tab = Pattern.compile("\t");
		
		genotypefiles = dir.listFiles(new FilenameFilter(){

			@Override
			public boolean accept(File dir, String filename) {
				if (filename.startsWith("ICIMdataForPopulation") && filename.endsWith(".txt")) return true;
				return false;
			}
		});
		
		for (File file:genotypefiles) {
			System.out.println("Reading " + file.getName() + " ...");
			try {
				BufferedReader br = new BufferedReader(new FileReader(file));
				String input = br.readLine();
				input = br.readLine();
				String[] info = tab.split(input);
				int nsnps = info.length - 1;
				int[] prev = new int[nsnps];
				for (int i = 0; i < nsnps; i++) prev[i] = Integer.parseInt(info[i + 1]);
				while ((input = br.readLine()) != null) {
					info = tab.split(input);
					for (int i = 0; i < nsnps; i++) {
						int val = Integer.parseInt(info[i + 1]);
						if (val > -1) {
							if (prev[i] > -1) {
								transitionCounts[prev[i]][val]++;
							}
							prev[i] = val;
						}
					}
				}
				br.close();
			} catch (IOException e) {
				e.printStackTrace();
				System.exit(-1);
			}
		}
		
		
		System.out.println("Counts: ");
		System.out.println(transitionCounts[0][0] + " " + transitionCounts[0][1] + " " + transitionCounts[0][2]);
		System.out.println(transitionCounts[1][0] + " " + transitionCounts[1][1] + " " + transitionCounts[1][2]);
		System.out.println(transitionCounts[2][0] + " " + transitionCounts[2][1] + " " + transitionCounts[2][2]);
		System.out.println();
	}
	
	public void setHdfFilename(String filename) {
		hdfFilename = filename;
	}
	
	public void setWriteHdf(boolean writeAsHdf) {
		this.writeAsHdf = writeAsHdf;
	}
}