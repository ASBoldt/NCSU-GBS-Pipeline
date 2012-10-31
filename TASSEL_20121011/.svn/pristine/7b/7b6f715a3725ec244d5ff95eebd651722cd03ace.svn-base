package net.maizegenetics.gwas.Imputation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Pattern;

import cern.jet.stat.Probability;

public class BiparentalGbsData {
	final Pattern tab = Pattern.compile("\t");
	ArrayList<String> taxaList;
	ArrayList<SnpInformation> snpList;
	int numberOfTaxa;
	int numberOfRils;
	int numberOfSnps;
	byte[][] imputedSnps; //-1 = missing, 0 = A, 1 = H, 2 = B
	byte[][] viterbiSnps;
	int[][] alleleCounts; //counts of A,C,G,T for each snp (rils only)
	byte[] parent1Genotype;
	byte[] parent2Genotype;
	int chromosome;
	String header;
	ArrayList<String> snpLabels;
	
	final static String[] byteToSnp = new String[]{"A","C","G","T","R","Y","S","W","K","M","+","0","H","B","N","-","X"};
	
	public BiparentalGbsData() {
		
	}
	
	public void imputeParents(String filename) {
		imputeBiparentalData(filename);
	}
	
	public void imputeBiparentalData(String filename) {
		String parent1 = "Al237";
		String parent2 = "L53";
		readMaronData(filename);
			
		//create merged parents
		//identify loci segregating 1:1 not 3:1 in progeny
		countAlleles();
		parent1Genotype = new byte[numberOfSnps];
		parent2Genotype = new byte[numberOfSnps];
		boolean[] isPolymorphic = new boolean[numberOfSnps];
		int p1a = 0;
		int p1b = 1;
		int p2a = numberOfTaxa - 2;
		int p2b = p2a + 1;
		
		int snpcount = 0;
		for (SnpInformation snp:snpList) {
			//merge parents
			byte val1a = snp.values[p1a];
			byte val1b = snp.values[p1b];
			byte val2a = snp.values[p2a];
			byte val2b = snp.values[p2b];
			
			if (val1a < 4 && val1b < 4 && val1a == val1b) parent1Genotype[snpcount] = val1a;
			else if (val1a < 4 && val1b == 14) parent1Genotype[snpcount] = val1a;
			else if (val1a == 14 && val1b < 4) parent1Genotype[snpcount] = val1b;
			else parent1Genotype[snpcount] = 14;
			if (val2a < 4 && val2b < 4 && val2a == val2b) parent2Genotype[snpcount] = val2a;
			else if (val2a < 4 && val2b == 14) parent2Genotype[snpcount] = val2a;
			else if (val2a == 14 && val2b < 4) parent2Genotype[snpcount] = val2b;
			else parent2Genotype[snpcount] = 14;
						
			//find major and minor alleles
			int major = 0;
			int minor = 0;
			int[] counts = alleleCounts[snpcount];
			for (int i = 1; i < 4; i++) {
				if (counts[i] > counts[major]) {
					minor = major;
					major = i;
				} else if (counts[i] > counts[minor]) {
					minor = i;
				}
			}
			
			//is the site segregating 1:1 (arbitrary cutoff)
			//are the parents segregating
			
			int total = counts[major] = counts[minor];
			double maf;
			if (minor == major) maf = 0;
			else maf = ((double) counts[minor]) / ((double) total);
//			if (total >= 40 && maf > 0.4 && parent1Genotype[snpcount] != 14 && parent2Genotype[snpcount] != 14 && parent2Genotype[snpcount] != parent1Genotype[snpcount]) isPolymorphic[snpcount] = true;
			if (total>=40 && maf>0.05 && maf < 0.3 && parent1Genotype[snpcount] != 14) isPolymorphic[snpcount] = true;
			else isPolymorphic[snpcount] = false;
			
			
			snpcount++;
		}
		
		//impute polymorphic snps, based on genotype of parent 1 (Al237)
		int numberOfPolymorphicSnps = 0;
		for (boolean b:isPolymorphic) if (b) numberOfPolymorphicSnps++;
		imputedSnps = new byte[numberOfRils][numberOfPolymorphicSnps];
		int[] snpPositions = new int[numberOfPolymorphicSnps];
		int polycount = 0;
		ArrayList<String> imputedSnpLabels = new ArrayList<String>(); 
		for (int s = 0; s < numberOfSnps; s++) {
			if (isPolymorphic[s]) {
				imputedSnpLabels.add(snpLabels.get(s));
				SnpInformation snp = snpList.get(s);
				snpPositions[polycount] = snp.position;
				for (int r = 0; r < numberOfTaxa -4; r++) {
					byte geno = snp.values[r + 2];
					if (geno == parent1Genotype[s]) {
						imputedSnps[r][polycount] = 0;
					} else if(geno < 4) {
						imputedSnps[r][polycount] = 2;
					} else if (geno < 10) {
						imputedSnps[r][polycount] = 1;
					} else {
						imputedSnps[r][polycount] = -1;
					}
				}
				polycount++;
			}
		}
		
		//check polymorphic SNPs for LD with imputed snps
		checkLD();
		
		//use Viterbi algorithm to impute states
		//get snp positions
		
		
		VariableAverageTransitionMatrix vatm = new VariableAverageTransitionMatrix(chromosome, snpPositions);
		vatm.calculateAverageDistance(imputedSnps);

		ProbabilityObsGivenStateBiparental probObsGivenState = new ProbabilityObsGivenStateBiparental();
		double[] pState = getStateProbabilities(imputedSnps);
		viterbiSnps = new byte[imputedSnps.length][imputedSnps[0].length];
		for (int t = 0; t < imputedSnps.length; t++) {
			for (int s = 0; s < imputedSnps[0].length; s++) {
				viterbiSnps[t][s] = -1;
			}
		}
		
		//for each taxon generate an index of nonmissing snps into imputedSnps
		int[][] nonMissingIndex = new int[numberOfRils][];
		int[] numberOfNonMissingSnps = new int[numberOfRils];
		int taxonCount = 0;
		for (byte[] row:imputedSnps) {
			int notmissing = 0;
			for (byte cell:row) {
				if (cell > -1) notmissing++;
			}
			numberOfNonMissingSnps[taxonCount] = notmissing;
			int[] index = new int[notmissing];
			nonMissingIndex[taxonCount] = index;
			int cellcount = 0;
			int nmcount = 0;
			for (byte cell:row) {
				if (cell > -1) {
					index[nmcount++] = cellcount;
				}
				cellcount++;
			}
			taxonCount++;
		}
		
		//run E-M iterations to convergence, 10 iterations is generally more than enough
		int nImputedSnps = imputedSnps[0].length;
		for (int iter = 0; iter < 10; iter++) {
			for (int t = 0; t < numberOfRils; t++) {
				byte[] obs = new byte[numberOfNonMissingSnps[t]];
				for (int i = 0; i < numberOfNonMissingSnps[t]; i++) {
					obs[i] = imputedSnps[t][nonMissingIndex[t][i]];
				}
				vatm.setPositionIndex(nonMissingIndex[t]);
				ViterbiAlgorithm va = new ViterbiAlgorithm(obs, vatm, probObsGivenState, pState);
				va.calculate();
				byte[] states = va.getMostProbableStateSequence();
				for (int i = 0; i < numberOfNonMissingSnps[t]; i++) {
					viterbiSnps[t][nonMissingIndex[t][i]] = states[i];
				}
			}
			
			//update probability matrices
			vatm.calculateAverageTransitionMatrix(viterbiSnps);
			pState = getStateProbabilities(viterbiSnps);
			
			int[][] stateCounts = new int[3][3];
			int[] rowCounts = new int[3];
			for (int t = 0; t < numberOfRils; t++) {
				for (int s = 0; s < nImputedSnps; s++) {
					if (viterbiSnps[t][s] > -1 && imputedSnps[t][s] > -1) {
						stateCounts[viterbiSnps[t][s]][imputedSnps[t][s]]++;
						rowCounts[viterbiSnps[t][s]]++;
					}
				}
			}
			
			System.out.println("Iteration " + iter + " state counts (rows are states, columns are observations)");
			for (int[] row:stateCounts) {
				for (int cell:row) System.out.print(cell +" ");
				System.out.println();
			}
			System.out.println();
			
			double[][] obsGivenState = new double[3][3];
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					obsGivenState[i][j] = ((double) stateCounts[i][j]) / ((double) rowCounts[i]);
				}
			}
			
			probObsGivenState.setProbabilityMatrix(obsGivenState);
		}
		
		//write results to a file
		String outfile = filename.substring(0, filename.lastIndexOf('.')) + ".imputedHMM.L53het.txt";
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
			bw.write(header);
			for (int i = 0; i < numberOfRils; i++) {
				bw.write("\t");
				bw.write(taxaList.get(i + 2));
			}
			bw.newLine();
			for (int i = 0; i < nImputedSnps; i++) {
				//does the snp have any data?
				int ril = 0;
				while (ril < numberOfRils && viterbiSnps[ril][i] == -1) ril++;
				if (ril < numberOfRils) { //if so, output to file
					bw.write(imputedSnpLabels.get(i));
					for (int j = 0; j < numberOfRils; j++) {
						bw.write("\t");
						bw.write(Byte.toString(viterbiSnps[j][i]));
					}
					bw.newLine();
				}
			}
			bw.close();
		} catch(IOException e) {
			e.printStackTrace();
		}
	}
	
	private double[] getStateProbabilities(byte[][] states) {
		int[] counts = new int[3];
		for (byte[] row:states) {
			for (byte cell:row) {
				if (cell > -1) counts[cell]++;
			}
		}
		
		double total = 0;
		for (int cnt:counts) total += cnt;
		
		double[] probs = new double[3];
		for (int i = 0; i < 3; i++) probs[i] = ((double) counts[i]) / total;
		return probs;
		
	}
	
	private void checkLD() {
		int window = 5;
		int nrils = imputedSnps.length;
		int nsnps = imputedSnps[0].length;
		double minSame = 0.8; //at least 80% of the neighboring snps must match
		double maxSame = 1 - minSame;  //if most are different, then parents can be swichted to make it acceptable
		System.out.println("Beginning check LD");

		//use framework imputedSnps for all comparisons and do not replace with new imputation until all have been checked 
		byte[][] tempImputed = new byte[nrils][nsnps];
		for (int i = 0; i < nrils; i++) {
			for (int j = 0; j < nsnps; j++) tempImputed[i][j] = -1;
		}

		for (int s = 0; s < nsnps; s++) {

				int match = 0;
				int total = 0;
				byte[] thisSnp = null;

				// if framework snp, then use values from imputedSnps, otherwise guess parentage from nucleotide data and check LD
				thisSnp = new byte[nrils];
				for (int r = 0; r < nrils; r++) {
					thisSnp[r] = imputedSnps[r][s];
				}

				for (int r = 0; r < nrils; r++) {
					if (thisSnp[r] == 0 || thisSnp[r] == 2) {
						//check previous non-missing, non-het snps
						int nNotMissing = 0;
						int comp = s - 1;
						while (comp >= 0 && nNotMissing < window) {
							byte compSnp = imputedSnps[r][comp];
							if (compSnp == 0 || compSnp == 2) {
								total++;
								nNotMissing++;
								if (compSnp == thisSnp[r]) match++;
							}
							comp--;
						}

						//check next non-missing, non-het snps
						nNotMissing = 0;
						comp = s + 1;
						while (comp < nsnps && nNotMissing < window) {
							byte compSnp = imputedSnps[r][comp];
							if (compSnp == 0 || compSnp == 2) {
								total++;
								nNotMissing++;
								if (compSnp == thisSnp[r]) match++;
							}
							comp++;
						}

					}
				}

				double psame = ((double) match) / ((double) total);
				if (psame > minSame) {
					for (int r = 0; r < nrils; r++) {
						tempImputed[r][s] = thisSnp[r];
					}
				} else if (psame < maxSame) {
					for (int r = 0; r < nrils; r++) {
						//switch parents
						byte val = thisSnp[r];
						if (val == 0) val = 2;
						else if (val == 2) val = 0;
						tempImputed[r][s] = val;
					}
				} 
			
		}

		//finished, set imputed to temp imputed
		imputedSnps = tempImputed;
		tempImputed = null;
	}
	
	private byte[] translateNucleotidesToParentCalls(byte[] nucleotideCalls, byte[] parents) {
		int n = nucleotideCalls.length;
		byte[] parentCalls = new byte[n - 4];
		
		for (int i = 0; i < n - 4; i++) {
			byte call = nucleotideCalls[i + 2];
			if (call == parents[0]) {
				parentCalls[i] = 0;
			} else if (call == parents[1]) {
				parentCalls[i] = 2;
			} else if (call > 3 && call < 10) {
				parentCalls[i] = 1;
			}
		}
		
		return parentCalls;
	}
	
	//returns null if the snp is not polymorphic, otherwise sets major allele to parent[0] and minor allele to parent[1]
	private byte[] setParents(int[] counts) {
		byte major = 0;
		byte minor = 0;
		for (byte i = 1; i < 4; i++) {
			if (counts[i] > counts[major]) {
				minor = major;
				major = i;
			} else if (counts[i] > counts[minor]) {
				minor = i;
			}
		}
		
		//criteria for being polymorphic
		boolean isPoly = false;
		if (counts[minor] > 4) {
			double maf = ((double) minor) / ((double) (major + minor));
			if (maf > 0.2) isPoly = true;
		}
		if (isPoly) return new byte[]{major, minor};
		else return null;
	}
	
	//count alleles in RILs not parents
	private void countAlleles() {
		alleleCounts = new int[numberOfSnps][4];
		int snpCount = 0;
		for (SnpInformation snp:snpList) {
			for (int t = 2; t < numberOfTaxa - 2; t++) {
				byte val = snp.values[t];
				if (val < 4) alleleCounts[snpCount][val]++;
			}
			snpCount++;
		}
	}
	
	public void readMaronData(String filename) {
		taxaList = new ArrayList<String>();
		snpList = new ArrayList<SnpInformation>();
		snpLabels = new ArrayList<String>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(filename));
			String input = br.readLine();
			String[] info = tab.split(input);
			numberOfTaxa = info.length - 11;
			numberOfRils = numberOfTaxa - 4;
			StringBuilder sb = new StringBuilder(info[0]);
			for (int i = 1; i < 11; i++) sb.append("\t").append(info[i]);
			header = new String(sb);
			for (int i = 0; i < numberOfTaxa; i++) {
				String longname = info[i + 11];
				taxaList.add(new String(longname.substring(0, longname.indexOf(":"))));
			}
			SnpInformation prevSnp = null;
			while ((input = br.readLine()) != null) {
				info = tab.split(input);
				SnpInformation snp = new SnpInformation();
				snp.name = new String(info[0]);
				snp.allele = new String(info[1]);
				snp.chromosome = Integer.parseInt(info[2]);
				snp.position = Integer.parseInt(info[3]);
				snp.values = new byte[numberOfTaxa];
				for (int t = 0; t < numberOfTaxa; t++) snp.values[t] = stringToByte(info[t + 11]);
				
				//check for duplicates
				if (prevSnp != null && prevSnp.name.equals(snp.name) && prevSnp.allele.equals(snp.allele)) {
					for (int t = 0; t < numberOfTaxa; t++) {
						if (prevSnp.values[t] == 14) prevSnp.values[t] = snp.values[t];
						else if (snp.values[t] != 14 && prevSnp.values[t] != snp.values[t]) prevSnp.values[t] = 14; 
					}
				} else {
					snpList.add(snp);
					
					sb = new StringBuilder(info[0]);
					for (int i = 1; i < 11; i++) sb.append("\t").append(info[i]);
					snpLabels.add(new String(sb));
					prevSnp = snp;
				}
				
			}
			br.close();
			numberOfSnps = snpList.size();
			
		} catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
		chromosome = snpList.get(0).chromosome;
	}
	
	public byte stringToByte(String geno) {
		if (geno.equals("N")) return 14;
		if (geno.equals("A")) return 0;
		if (geno.equals("C")) return 1;
		if (geno.equals("G")) return 2;
		if (geno.equals("T")) return 3;
		if (geno.equals("R")) return 4;
		if (geno.equals("Y")) return 5;
		if (geno.equals("S")) return 6;
		if (geno.equals("W")) return 7;
		if (geno.equals("K")) return 8;
		if (geno.equals("M")) return 9;
		return 14;
	}
	
	
	
	public void writeImputedData(String filename) {
		String[] geno = new String[]{"N", "A", "M", "C"};
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(filename));
			
			//header
			bw.write("rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode");
			for (int t = 0; t < numberOfTaxa - 4; t++) {
				bw.write("\t");
				bw.write(taxaList.get(t + 2));
			}
			
			//alleles coded as A,C,M where A=parent1, C=parent2, M=het
			int snpcount = 0;
			for (SnpInformation snp:snpList) {
				bw.write(snp.name);
				bw.write("\t");
				bw.write(snp.allele);
				bw.write("\t");
				bw.write(snp.chromosome);
				bw.write("\t");
				bw.write(snp.position);
				bw.write("\tNA\tNA\tNA\tNA\tNA\tNA\tNA");
				for (int t = 0; t < numberOfTaxa - 4; t++) {
					bw.write("\t");
					bw.write(geno[imputedSnps[t][snpcount] + 1]);
				}
				snpcount++;
			}
				
			
			bw.close();
		} catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
}
