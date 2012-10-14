 package net.maizegenetics.gwas.Imputation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.regex.Pattern;

import net.maizegenetics.gwas.NAMgwas.AGPMap;

import cern.jet.stat.Gamma;

public class NamGbsImputation {
	private int[] pops = null;
	private int[] subpops = null;
	
	public void imputeUsingViterbi(String filename) {
		
		//get gbs data
		System.out.println("Reading gbs data...");
//		String gbsGenotypeFilename = "/Volumes/Macintosh HD 2/data/namgbs/genos_20111003/NAM_20111003_ABHv2_chr" + chromosome + ".txt";
//		String gbsGenotypeFilename = "/Volumes/Macintosh HD 2/data/namgbs/genos_1123/NAM_20111123_SNPmerge15_cov10_fT1E1pLD_mergedTaxa_chr" + chromosome + ".hmp.hmm.txt";
		GbsGenotypeData gbsGenotype = new GbsGenotypeData(filename, false);
		
		//get the taxalist
		String taxalistFilename = "/Volumes/Macintosh HD 2/data/namgbs/namgbstaxalist.txt";
		System.out.println("Reading taxa list...");
		LinkedList<String> taxaList = new LinkedList<String>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(taxalistFilename));
			String input;
			while ((input = br.readLine()) != null) {
				if (input.startsWith("Z0")) taxaList.add(input);
  			}
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}

		//find the intersection of the taxa list and the taxa in the genotype data
		ArrayList<String> gbsTaxa = gbsGenotype.getTaxaList();
		Collections.sort(gbsTaxa);
		int ntaxa = taxaList.size();
		for (int t = ntaxa - 1; t >= 0; t--) {
			String taxon = taxaList.get(t);
			if (Collections.binarySearch(gbsTaxa, taxon) < 0) {
				taxaList.remove(t);
			}
		}
		getSubPops(taxaList);
		
		
		//create an array for the imputed Snps
		ntaxa = taxaList.size();
		int nsnps = gbsGenotype.getNumberOfSnps();
		ArrayList<SnpInformation> gbsSnpList = gbsGenotype.getSnpList();
		
		//counts to keep track of agreement between observations and imputed states. First dimension is observation, second is state.
//		int[][] NAMmatchCounts;
		int[][] GBSmatchCounts;

		//use EM + Viterbi to impute genotype
		System.out.println("Imputing snps from " + filename);
		double[][] gbsProbObsGivenState = new double[][]{{0.98,0.45,0.02},{0.01,.1,0.01},{0.01,0.45,0.97}};
		double[][] arrayProbObsGivenState = new double[][]{{0.9998,0.0001,0.0001},{0.0001,0.9998,0.0001},{0.0001,0.0001,0.9998}};
		ProbabilityObservationGivenTrueState pObsGivenTrue= new ProbabilityObservationGivenTrueState();
//		ProbObsGivenState pObsGivenTrue = new ProbObsGivenState();
		pObsGivenTrue.setGbsProbabilityMatrix(gbsProbObsGivenState);
		pObsGivenTrue.setArrayProbabilityMatrix(arrayProbObsGivenState);
		System.out.println("Setting initial probabilities...");
		pObsGivenTrue.setProbabilityOfAnObservation(gbsGenotype.getStateProbabilities());
		System.out.println("Setting probability monomorphic...");
		
		setProbabilityMonomorphicBySubpop(gbsGenotype, taxaList, pObsGivenTrue, 0.1, .002);
//		setProbabilityMonomorphicBySubpop5state(gbsGenotype, taxaList, pObsGivenTrue, 0.1, .002);

		double[] pstate = new double[]{.485, .03, .485};
		
		int[] positions = new int[nsnps];
		int count = 0;
		for (SnpInformation snp:gbsSnpList) {
			positions[count++] = snp.position;
		}
		
		byte[][] imputedStates = new byte[ntaxa][];
		
		int taxacount = 0;
		for (String taxon:taxaList) {
			imputedStates[taxacount++] = gbsGenotype.getSnpsForTaxon(taxon);
		}
		
//		TestTransitionMatrix vatm = new TestTransitionMatrix(.002);
		VariableAverageTransitionMatrix vatm = new VariableAverageTransitionMatrix(gbsGenotype.getChromosome(), positions);
//		VariableAverageTransitionMatrix vatm = new VariableTransitionMatrix5State(gbsGenotype.getChromosome(), positions);
		vatm.calculateAverageDistance(imputedStates);
		
		for (int iter = 0; iter < 10; iter++) {
			System.out.println("Starting iteration " + iter);
			
			taxacount = 0;
			GBSmatchCounts = new int[3][3];

			for (String taxon:taxaList) {
				byte[] snpcalls = gbsGenotype.getSnpsForTaxon(taxon);

				//extract nonmissing snps 
				LinkedList<int[]> nonmissingList = new LinkedList<int[]>();
				for (int s = 0; s < nsnps; s++) {
					if (snpcalls[s] != -1) {
						nonmissingList.add(new int[]{s, snpcalls[s], gbsSnpList.get(s).position});
					}
				}

				int numberNotMissing = nonmissingList.size();

				if (numberNotMissing < 50) {
//					System.out.println(taxon + " has " + numberNotMissing + " non-missing values and will not be imputed.");
				}
				else {
					byte[] geno = new byte[numberNotMissing];
					int[] snpPosition = new int[numberNotMissing];
					boolean[] isgbs = new boolean[numberNotMissing];
					int[] obsIndex = new int[numberNotMissing];

					count = 0;
					for (int[] snp : nonmissingList) {
						obsIndex[count] = snp[0];
						geno[count] = (byte) snp[1];
						snpPosition[count] = snp[2];
						isgbs[count++] = true;
					}

					//use the Viterbi algorithm to calculate the most probable state (maximization step)
					//				recomb.setEvaluationPositions(snpPosition);
					pObsGivenTrue.setIsGBS(isgbs);
					pObsGivenTrue.setSubpop(pops[taxacount], subpops[taxacount]);
					pObsGivenTrue.setObservationIndex(obsIndex);
					vatm.setPositionIndex(obsIndex);
					ViterbiAlgorithm  va = new ViterbiAlgorithm(geno, vatm, pObsGivenTrue, pstate);
					va.calculate();
					byte[] genostates = va.getMostProbableStateSequence();

					//update imputedStates
					count = 0;
					for (int[] snp : nonmissingList) {
						imputedStates[taxacount][snp[0]] = genostates[count++];
					}

					//count the changes
					for (int i = 0; i < numberNotMissing; i++) {
						GBSmatchCounts[geno[i]][genostates[i]]++;
					}

				}
				taxacount++;
			}
			
			//report the changes
			System.out.println("GBS match counts, row = observed, column = imputed:");
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					System.out.print(GBSmatchCounts[i][j] + "\t");
				}
				System.out.println();
			}
			System.out.println();
			
			//recalculate the probability matrices (expectation step)
			//calculate p(state)
			int[] nstate = new int[3];
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					nstate[j] += GBSmatchCounts[i][j];
//					nstate[j] += NAMmatchCounts[i][j];
				}
			}
			
			double total = nstate[0] + nstate[1] + nstate[2];
			for (int i = 0; i < 3; i++) {
				pstate[i] = ((double) nstate[i]) / total;
			}

			//calculate p(obs|state) for GBS and array, obs row, state column 
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					gbsProbObsGivenState[i][j] = ((double) GBSmatchCounts[i][j]) / ((double) nstate[j]);
//					arrayProbObsGivenState[i][j] = NAMmatchCounts[i][j] / arrayColumnTotal[j];
				}
			}
			
//			System.out.println("gbsProbObsGivenState:");
//			for (double[] row : gbsProbObsGivenState) {
//				for (double col:row) {
//					System.out.print(col +"  ");
//				}
//				System.out.println();
//			}
			
//			pObsGivenTrue.setArrayProbabilityMatrix(arrayProbObsGivenState);
			pObsGivenTrue.setGbsProbabilityMatrix(gbsProbObsGivenState);

			//calculate the transition matrix
			double transitionRatio = vatm.calculateAverageTransitionMatrix(imputedStates);
			
			if (iter == 9) {
				//percent het
				System.out.println("percent het\ttransition ratio");
				System.out.println(pstate[1] + "\t" + transitionRatio);
				
				//het/hom transition ratio
			}

		}
		
		//write the results to a file
		try {
			int lastPeriod = filename.lastIndexOf('.');
			String outFilename = filename.substring(0, lastPeriod) + ".imputeHMM" + filename.substring(lastPeriod);
			BufferedWriter bw = new BufferedWriter(new FileWriter(outFilename));
			bw.write("rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode");
			for (String taxon:taxaList) {
				bw.write("\t");
				bw.write(taxon);
			}
			bw.newLine();
			int snpCount = 0;
			for (SnpInformation snp:gbsGenotype.getSnpList()) {
				bw.write(snp.name);
				bw.write("\tNA\t");
				bw.write(Integer.toString(snp.chromosome));
				bw.write("\t");
				bw.write(Integer.toString(snp.position));
				bw.write("\tNA\tNA\tNA\tNA\tNA\tNA\tNA");
				for (int t = 0; t < ntaxa; t++) {
					bw.write("\t");
					bw.write(Byte.toString(imputedStates[t][snpCount]));
				}
				snpCount++;
				bw.newLine();
			}
			
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
	}
	
	public void imputeUsingFiveStateViterbi(String filename) {
		int numberOfStates = 5;
		//get gbs data
		System.out.println("Reading gbs data...");
		GbsGenotypeData gbsGenotype = new GbsGenotypeData(filename, false);
		
		//get the taxalist
		String taxalistFilename = "/Volumes/Macintosh HD 2/data/namgbs/namgbstaxalist.txt";
		System.out.println("Reading taxa list...");
		LinkedList<String> taxaList = new LinkedList<String>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(taxalistFilename));
			String input;
			while ((input = br.readLine()) != null) {
				if (input.startsWith("Z0")) taxaList.add(input);
  			}
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}

		//find the intersection of the taxa list and the taxa in the genotype data
		ArrayList<String> gbsTaxa = gbsGenotype.getTaxaList();
		Collections.sort(gbsTaxa);
		int ntaxa = taxaList.size();
		for (int t = ntaxa - 1; t >= 0; t--) {
			String taxon = taxaList.get(t);
			if (Collections.binarySearch(gbsTaxa, taxon) < 0) {
				taxaList.remove(t);
			}
		}
		getSubPops(taxaList);
		
		
		//create an array for the imputed Snps
		ntaxa = taxaList.size();
		int nsnps = gbsGenotype.getNumberOfSnps();
		ArrayList<SnpInformation> gbsSnpList = gbsGenotype.getSnpList();
		
		//counts to keep track of agreement between observations and imputed states. First dimension is observation, second is state.
		int[][] GBSmatchCounts;

		//use EM + Viterbi to impute genotype
		System.out.println("Imputing snps from " + filename);
		
		
		
		double[][] obsGivenState = new double[][]{{0.998,0.6,0.4,0.2,0.001},{0.001,0.2,0.2,0.2,0.001},{0.001,0.2,0.4,0.6,0.998}};
//		double[][] obsGivenState = new double[][]{{0.98,0.45,0.02},{0.01,.1,0.01},{0.01,0.45,0.97}};
		ProbObsGivenState pObsGivenTrue= new ProbObsGivenState();
		pObsGivenTrue.setGbsProbabilityMatrix(obsGivenState);

		System.out.println("Setting probability monomorphic...");
		
		setProbabilityMonomorphicBySubpop5state(gbsGenotype, taxaList, pObsGivenTrue, 0.1, .002);
		double[] pstate = new double[]{.485,.001,.0028,.001,.485};
//		double[] pstate = new double[]{.485,.003,.485};
		
		int[] positions = new int[nsnps];
		int count = 0;
		for (SnpInformation snp:gbsSnpList) {
			positions[count++] = snp.position;
		}
		
		byte[][] imputedStates = new byte[ntaxa][];
		
		int taxacount = 0;
		for (String taxon:taxaList) {
			imputedStates[taxacount++] = gbsGenotype.getSnpsForTaxon(taxon);
		}
		
		VariableTransitionMatrix5State vtm = new VariableTransitionMatrix5State(gbsGenotype.getChromosome(), positions);
		double[][] initialTransition = new double[][]{
				{0.9995,0.00006,0.00013,0.00006,0.00025},
				{0.00045,0.999,0.00005,0.00005,0.00045},
				{0.00045,0.00005,0.999,0.00005,0.00045},
				{0.00045,0.00005,0.00005,0.999,0.00045},
				{0.00025,0.00006,0.00013,0.00006,0.9995}};
//		double[][] initialTransition = new double[][]{
//				{0.9995,0.00025,0.00025},
//				{0.0005,0.999,0.0005},
//				{0.00025,0.00025,0.9995}};
		vtm.setAverageTransitionMatrix(initialTransition);
		vtm.calculateAverageDistance(imputedStates);
		
		for (int iter = 0; iter < 15; iter++) {
			System.out.println("Starting iteration " + iter);
			
			taxacount = 0;
			GBSmatchCounts = new int[3][numberOfStates];

			for (String taxon:taxaList) {
				byte[] snpcalls = gbsGenotype.getSnpsForTaxon(taxon);

				//extract nonmissing snps 
				LinkedList<int[]> nonmissingList = new LinkedList<int[]>();
				for (int s = 0; s < nsnps; s++) {
					if (snpcalls[s] != -1) {
						nonmissingList.add(new int[]{s, snpcalls[s], gbsSnpList.get(s).position});
					}
				}

				int numberNotMissing = nonmissingList.size();

				if (numberNotMissing < 50) {
//					System.out.println(taxon + " has " + numberNotMissing + " non-missing values and will not be imputed.");
				} else {
					byte[] geno = new byte[numberNotMissing];
					int[] snpPosition = new int[numberNotMissing];
					boolean[] isgbs = new boolean[numberNotMissing];
					int[] obsIndex = new int[numberNotMissing];

					count = 0;
					for (int[] snp : nonmissingList) {
						obsIndex[count] = snp[0];
						geno[count] = (byte) snp[1];
						snpPosition[count] = snp[2];
						isgbs[count++] = true;
					}

					//use the Viterbi algorithm to calculate the most probable state (maximization step)
					pObsGivenTrue.setSubpop(pops[taxacount], subpops[taxacount]);
					pObsGivenTrue.setObservationIndex(obsIndex);
					pObsGivenTrue.setIsGBS(isgbs);
					vtm.setPositionIndex(obsIndex);
					ViterbiAlgorithm  va = new ViterbiAlgorithm(geno, vtm, pObsGivenTrue, pstate);
					va.calculate();
					byte[] genostates = va.getMostProbableStateSequence();

					//update imputedStates
					count = 0;
					for (int[] snp : nonmissingList) {
						imputedStates[taxacount][snp[0]] = genostates[count++];
					}

					//count the changes
					for (int i = 0; i < numberNotMissing; i++) {
						GBSmatchCounts[geno[i]][genostates[i]]++;
					}
				}
				taxacount++;
			}
			
			//report the changes
			System.out.println("GBS match counts, row = observed, column = imputed:");
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < numberOfStates; j++) {
					System.out.print(GBSmatchCounts[i][j] + "\t");
				}
				System.out.println();
			}
			System.out.println();
			
			//recalculate the probability matrices (expectation step)
			//calculate p(state)
			int[] nstate = new int[numberOfStates];
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < numberOfStates; j++) {
					nstate[j] += GBSmatchCounts[i][j];
				}
			}
			
			double total = 0;
			for (int i = 0; i < numberOfStates; i++) total += nstate[i];
			for (int i = 0; i < numberOfStates; i++) {
				pstate[i] = ((double) nstate[i]) / total;
			}

			//calculate p(obs|state) for GBS and array, obs row, state column 
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < numberOfStates; j++) {
					obsGivenState[i][j] = ((double) GBSmatchCounts[i][j]) / ((double) nstate[j]);
				}
			}
			
			//calculate the transition matrix
			double transitionRatio = vtm.calculateAverageTransitionMatrix(imputedStates);
			
		}
		
		//write the results to a file
		try {
			int lastPeriod = filename.lastIndexOf('.');
			String outFilename = filename.substring(0, lastPeriod) + ".impute5stateHMM" + filename.substring(lastPeriod);
			BufferedWriter bw = new BufferedWriter(new FileWriter(outFilename));
			bw.write("rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode");
			for (String taxon:taxaList) {
				bw.write("\t");
				bw.write(taxon);
			}
			bw.newLine();
			int snpCount = 0;
			for (SnpInformation snp:gbsGenotype.getSnpList()) {
				bw.write(snp.name);
				bw.write("\tNA\t");
				bw.write(Integer.toString(snp.chromosome));
				bw.write("\t");
				bw.write(Integer.toString(snp.position));
				bw.write("\tNA\tNA\tNA\tNA\tNA\tNA\tNA");
				for (int t = 0; t < ntaxa; t++) {
					bw.write("\t");
					bw.write(Byte.toString(imputedStates[t][snpCount]));
				}
				snpCount++;
				bw.newLine();
			}
			
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
		//write the results converted to 3 state
		try {
			int lastPeriod = filename.lastIndexOf('.');
			String outFilename = filename.substring(0, lastPeriod) + ".impute5to3stateHMM" + filename.substring(lastPeriod);
			BufferedWriter bw = new BufferedWriter(new FileWriter(outFilename));
			bw.write("rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode");
			for (String taxon:taxaList) {
				bw.write("\t");
				bw.write(taxon);
			}
			bw.newLine();
			int snpCount = 0;
			for (SnpInformation snp:gbsGenotype.getSnpList()) {
				bw.write(snp.name);
				bw.write("\tNA\t");
				bw.write(Integer.toString(snp.chromosome));
				bw.write("\t");
				bw.write(Integer.toString(snp.position));
				bw.write("\tNA\tNA\tNA\tNA\tNA\tNA\tNA");
				for (int t = 0; t < ntaxa; t++) {
					bw.write("\t");
					String strval;
					switch (imputedStates[t][snpCount]) {
					case 0:
						strval = "0";
						break;
					case 1:
						strval = "1";
						break;
					case 2:
						strval = "1";
						break;
					case 3:
						strval = "1";
						break;
					case 4:
						strval = "2";
						break;
					default:
						strval = "-1";
						break;
					}
					bw.write(strval);
				}
				snpCount++;
				bw.newLine();
			}
			
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
	}
	
	public void monomorphismCheck(byte[][] states, ArrayList<SnpInformation> snpList, ArrayList<String> taxaList) {
		//find the subpopulations, entries 1-67, 68-134, 135-200
		ArrayList<String> sortedTaxa = new ArrayList<String>(taxaList);
		Collections.sort(sortedTaxa);
		
		//within subpops, check percent het in a window around each snp position across all taxa in the subpop
		
		//within each subpopulation, if a snp is in a region that is ~50% het and lacks evidence that it is polymorphic, and has the B73 allele, delete it
		
		
	}
	
	public void getSubPops(LinkedList<String> taxaList) {
		int ntaxa = taxaList.size();
		pops = new int[ntaxa];
		subpops = new int[ntaxa];
		int taxacount = 0;
		HashMap<String, Integer> popmap = new HashMap<String, Integer>();
		popmap.put("Z001", 0);
		popmap.put("Z002", 1);
		popmap.put("Z003", 2);
		popmap.put("Z004", 3);
		popmap.put("Z005", 4);
		popmap.put("Z006", 5);
		popmap.put("Z007", 6);
		popmap.put("Z008", 7);
		popmap.put("Z009", 8);
		popmap.put("Z010", 9);
		popmap.put("Z011", 10);
		popmap.put("Z012", 11);
		popmap.put("Z013", 12);
		popmap.put("Z014", 13);
		popmap.put("Z015", 14);
		popmap.put("Z016", 15);
		popmap.put("Z018", 16);
		popmap.put("Z019", 17);
		popmap.put("Z020", 18);
		popmap.put("Z021", 19);
		popmap.put("Z022", 20);
		popmap.put("Z023", 21);
		popmap.put("Z024", 22);
		popmap.put("Z025", 23);
		popmap.put("Z026", 24);
		
		for (String taxon:taxaList) {
			pops[taxacount] = popmap.get(taxon.substring(0,4));
			int ent = Integer.parseInt(taxon.substring(6));
			if (ent <= 67) subpops[taxacount] = 0;
			else if (ent <= 134) subpops[taxacount] = 1;
			else subpops[taxacount] = 2;
			taxacount++;
		}
	}
	
	public void setProbabilityMonomorphicBySubpopX(GbsGenotypeData gbs, LinkedList<String> taxaList, ProbabilityObservationGivenTrueState pObsGivenState) {
		//P(mono|obs) = P(obs|mono)P(mono)/[P(obs|mono)P(mono) + P(obs|!mono)P(!mono)]
		//nA = number of states = 0, nT = total number of states, B(n_success, n_trials, probability of a success) = binomial distribution
		//p(mono) = B(nA, nT, 0) / [B(nA, nT, 0.999) + B(nA, nT, 0.5)]
		
		//create a taxa index that indexes the gbs snps on taxalist
		int ntaxa = taxaList.size();
		int[] taxaIndex = new int[ntaxa];
		HashMap<String, Integer> gbsTaxaMap = gbs.getTaxaMap();
		int taxaCount = 0;
		for (String taxon:taxaList) taxaIndex[taxaCount++] = gbsTaxaMap.get(taxon);
		
		int nsnps = gbs.getNumberOfSnps();
		ArrayList<SnpInformation> snpList = gbs.getSnpList();
		double[][][] prob = new double[25][3][nsnps];
		for (SnpInformation snp:snpList) {
			int[][][] stateCountBySubpop = new int[25][3][3];
			int taxacount = 0;
			byte[] snpvalues = snp.values;
			for (int t = 0; t < ntaxa; t++) {
				int ndx = taxaIndex[t];
				if (snpvalues[ndx] != -1) stateCountBySubpop[pops[t]][subpops[t]][snpvalues[ndx]]++;
				taxacount++;
			}
			for (int i = 0; i < 25; i++) {
				for (int j = 0; j < 3; j++) {
					int nt = 0;
					int ns = stateCountBySubpop[i][j][0];
					for (int cnt : stateCountBySubpop[i][j]) nt += cnt;
					if (nt == 0) {
						prob[i][j][snp.ndx] = Double.NaN;
					} else {
						double prob0 = binomialProbability(ns, nt, 0.999);
						double prob5 = binomialProbability(ns, nt, 0.5);
						prob[i][j][snp.ndx] = prob0 / (prob0 + prob5);
					}
				}
			}
		}
		
		pObsGivenState.setProbMonomorphic(prob);
	}
	
	public void setProbabilityMonomorphicBySubpop(GbsGenotypeData gbs, LinkedList<String> taxaList, ProbabilityObservationGivenTrueState pObsGivenState, double pmono, double pError) {
		//P(mono|obs) = P(obs|mono)P(mono)/[P(obs|mono)P(mono) + P(obs|!mono)P(!mono)]
		
		
		//create a taxa index that indexes the gbs snps on taxalist
		int ntaxa = taxaList.size();
		int[] taxaIndex = new int[ntaxa];
		HashMap<String, Integer> gbsTaxaMap = gbs.getTaxaMap();
		int taxaCount = 0;
		for (String taxon:taxaList) taxaIndex[taxaCount++] = gbsTaxaMap.get(taxon);
		
		int nsnps = gbs.getNumberOfSnps();
		ArrayList<SnpInformation> snpList = gbs.getSnpList();
		double[][][] prob = new double[25][3][nsnps];
		for (SnpInformation snp:snpList) {
			int[][][] stateCountBySubpop = new int[25][3][3];
			int taxacount = 0;
			byte[] snpvalues = snp.values;
			for (int t = 0; t < ntaxa; t++) {
				int ndx = taxaIndex[t];
				if (snpvalues[ndx] != -1) stateCountBySubpop[pops[t]][subpops[t]][snpvalues[ndx]]++;
				taxacount++;
			}
			for (int i = 0; i < 25; i++) {
				for (int j = 0; j < 3; j++) {
					int nt = 0;
					int ns = stateCountBySubpop[i][j][0];
					for (int cnt : stateCountBySubpop[i][j]) nt += cnt;
//					int nerr = nt - Math.max(ns, stateCountBySubpop[i][j][2]);
					if (nt == 0) {
						prob[i][j][snp.ndx] = Double.NaN;
					} else {
						double pNotMono = binomialProbability(ns, nt, 0.5);
						double pMono = binomialProbability(ns, nt, 1 - pError);
						
						prob[i][j][snp.ndx] = pMono * pmono /(pMono * pmono + pNotMono * (1 - pmono));
					}
				}
			}
		}
		
		pObsGivenState.setProbMonomorphic(prob);
	}
	
	public void setProbabilityMonomorphicBySubpop5state(GbsGenotypeData gbs, LinkedList<String> taxaList, ProbObsGivenState pObsGivenState, double pmono, double pError) {
		//P(mono|obs) = P(obs|mono)P(mono)/[P(obs|mono)P(mono) + P(obs|!mono)P(!mono)]
		
		
		//create a taxa index that indexes the gbs snps on taxalist
		int ntaxa = taxaList.size();
		int[] taxaIndex = new int[ntaxa];
		HashMap<String, Integer> gbsTaxaMap = gbs.getTaxaMap();
		int taxaCount = 0;
		for (String taxon:taxaList) taxaIndex[taxaCount++] = gbsTaxaMap.get(taxon);
		
		int nsnps = gbs.getNumberOfSnps();
		ArrayList<SnpInformation> snpList = gbs.getSnpList();
		double[][][] prob = new double[25][3][nsnps];
		for (SnpInformation snp:snpList) {
			int[][][] stateCountBySubpop = new int[25][3][3];
			int taxacount = 0;
			byte[] snpvalues = snp.values;
			for (int t = 0; t < ntaxa; t++) {
				int ndx = taxaIndex[t];
				if (snpvalues[ndx] != -1) stateCountBySubpop[pops[t]][subpops[t]][snpvalues[ndx]]++;
				taxacount++;
			}
			for (int i = 0; i < 25; i++) {
				for (int j = 0; j < 3; j++) {
					int nt = 0;
					int ns = stateCountBySubpop[i][j][0];
					for (int cnt : stateCountBySubpop[i][j]) nt += cnt;
//					int nerr = nt - Math.max(ns, stateCountBySubpop[i][j][2]);
					if (nt == 0) {
						prob[i][j][snp.ndx] = Double.NaN;
					} else {
						double p1to3 = binomialProbability(ns, nt, 0.25);
						double p1to1 = binomialProbability(ns, nt, 0.5);
						double p3to1 = binomialProbability(ns, nt, 0.75);
						double pMono = binomialProbability(ns, nt, 1 - pError);
						
						prob[i][j][snp.ndx] = pMono * pmono /(pMono * pmono + p1to3 * 0.25 * (1 - pmono) + p1to1 * 0.5* (1 - pmono) + p3to1 * 0.25* (1 - pmono));
						
					}
				}
			}
		}
		
		pObsGivenState.setProbMonomorphic(prob);
	}
	
	public static double binomialProbability(double nSuccess, double nTrials, double p) {
		double prob;
		if (nSuccess == 0) {
			prob = Math.pow(1 - p, nTrials);
		} else {
			double logprob = Gamma.logGamma(nSuccess + nTrials) - Gamma.logGamma(nSuccess) - Gamma.logGamma(nTrials) + nSuccess*Math.log(p) + (nTrials - nSuccess)*Math.log(1 - p);
			prob = Math.exp(logprob);
		}
		return prob;
	}
	
	public static void imputeLinkageMarkers(double interval) {
		int chromosome = 1;
		String snpFilename = "/Volumes/Macintosh HD 2/data/namgbs/genos_1217/NAM282_20111217_scv10mF8maf002_mgs_E1pLD5_chr" + chromosome + ".hmp.abhv2.impute5to3stateHMM.txt";
		String outFilename = "/Volumes/Macintosh HD 2/data/namgbs/genos_1217/NAM282_20111217_scv10mF8maf002_mgs_E1pLD5_hmp_" + interval + "cmsnps.txt";
		Pattern tab = Pattern.compile("\t");
		AGPMap agpmap = new AGPMap(true);
		
		int ntaxa = 0;
		BufferedWriter bw = null;
		try {
			BufferedReader br = new BufferedReader(new FileReader(snpFilename));
			bw = new BufferedWriter(new FileWriter(outFilename));
			String header = br.readLine(); 
			String[] info = tab.split(header);
			int ncol = info.length;
			ntaxa = ncol - 11;
			
			bw.write("Snp\tallele\tchr\tpos\tcm");
			for (int t = 11; t < ncol; t++) {
				bw.write("\t");
				bw.write(info[t]);
			}
			bw.newLine();
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
		//impute data for each chromosome
		for (int chr = 1; chr <=1; chr++) {
			System.out.println("Imputing data for chromosome " + chr);
			snpFilename = "/Volumes/Macintosh HD 2/data/namgbs/genos_1217/NAM282_20111217_scv10mF8maf002_mgs_E1pLD5_chr" + chr + ".hmp.abhv2.impute5to3stateHMM.txt";
			int nsnps = 0;
			
			//count the number of snps
			try {
				BufferedReader br = new BufferedReader(new FileReader(snpFilename));
				br.readLine();
				while (br.readLine() != null) nsnps++;
				br.close();
			} catch(IOException e) {
				e.printStackTrace();
				System.exit(-1);
			}
			
			//read the snp data
			byte[][] snps = new byte[ntaxa][nsnps];
			int[] pos = new int[nsnps];
			
			try {
				BufferedReader br = new BufferedReader(new FileReader(snpFilename));
				String input = br.readLine();
				int snpcount = 0;
				while ((input = br.readLine()) != null) {
					String[] info = tab.split(input);
					pos[snpcount] = Integer.parseInt(info[3]);
					for (int t = 0; t < ntaxa; t++) snps[t][snpcount] = Byte.parseByte(info[t + 11]);
					snpcount++;
				}
				br.close();
			} catch (IOException e) {
				e.printStackTrace();
				System.exit(-1);
			}

			double startgenpos = agpmap.getCmFromPosition(chr, pos[0]);
			//round up to nearest interval
			startgenpos = ((double) (Math.ceil(startgenpos / interval))) * interval;
			
			double endgenpos = agpmap.getCmFromPosition(chr, pos[nsnps - 1]);
			//round down to nearest interval
			endgenpos = ((double)(Math.floor(endgenpos / interval))) * interval;
			
			int leftflank = 0;
			int rightflank = 0;
			
			String chrstr = Integer.toString(chr);
			try {
				for (double curpos = startgenpos; curpos <= endgenpos; curpos += interval) {
					int physpos = agpmap.getPositionFromCm(chr, curpos);
					String physposString = Integer.toString(physpos);
					String genpos = Double.toString(curpos);
					bw.write("S_");
					bw.write(physposString);
					bw.write("\timputed\t");
					bw.write(chrstr);
					bw.write("\t");
					bw.write(physposString);
					bw.write("\t");
					bw.write(genpos);
					
					while (physpos > pos[rightflank]) rightflank++;
					leftflank = rightflank - 1;
					for (int t = 0; t < ntaxa; t++) {
						bw.write("\t");
						int leftndx = leftflank;
						int rightndx = rightflank;
						while (snps[t][leftndx] == -1 && leftndx > 0) leftndx--;
						while (snps[t][rightndx] == -1 && rightndx < nsnps - 1) rightndx++;
						if (snps[t][leftndx] == snps[t][rightndx]) bw.write(Byte.toString(snps[t][leftndx]));
						else if (snps[t][leftndx] == -1) bw.write(Byte.toString(snps[t][rightndx]));
						else if (snps[t][rightndx] == -1) bw.write(Byte.toString(snps[t][leftndx]));
						else {
							double pd = ((double)(physpos - pos[leftndx])) / ((double)(pos[rightndx] - pos[leftndx]));
							double val = (1.0 - pd) * snps[t][leftndx] + pd * snps[t][rightndx];
							bw.write(Double.toString(val));
						}
//						else bw.write("-");   //use this for R/qtl type output
					}
					bw.newLine();
				}
				bw.close();
			} catch (IOException e) {
				e.printStackTrace();
				System.exit(-1);
			}
			
		}
		
		try {
			bw.close();
		} catch(IOException e) {
			e.printStackTrace();
		}
		
		System.out.println("Finished imputing markers.");
	}
}
