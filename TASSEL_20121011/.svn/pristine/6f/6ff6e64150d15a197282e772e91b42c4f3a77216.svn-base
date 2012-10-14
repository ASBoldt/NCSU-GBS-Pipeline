package net.maizegenetics.gwas.Imputation;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;

/**
 * @author pbradbury
 *
 */
public class VariableAverageTransitionMatrix implements TransitionMatrix {
	private double[][] rates;
	private int numberOfRates;
//	private double[][] averageTransitionMatrix = new double[][]{{.98, .004, .016},{.05, .9, .05},{.016, .004, .98}};
	private double[][] averageTransitionMatrix = new double[][]{{.9995, .00025, .00025},{.00050, .999, .00050},{.00025, .00025, .9995}};
	private final int[] positions;
	private double averageRecombinationRate;
	private double relativeRecombination;
	private double avgDistance;
	private int[] positionIndex;
	
	/**
	 * @param chromosome	the chromosome that this matrix represents
	 * @param positions		the snp positions on this chromosome
	 */
	public VariableAverageTransitionMatrix(int chromosome, int[] positions) {
		this.positions = positions;
//		String rateFilename = "/Volumes/Macintosh HD 2/data/namgbs/cov20/cov20_chr_chr" + chromosome + ".recombinationrate.het12.win1.txt";
//		readRates(rateFilename);
	}
	
	@Override
	public double getTransitionProbability(int state1, int state2) {
		return averageTransitionMatrix[state1][state2];
	}

	@Override
	public double getLnTransitionProbability(int state1, int state2) {
		return Math.log(getTransitionProbability(state1, state2));
	}

	@Override
	public double getLnTransitionProbability(int state1, int state2, int node) {
		int ndx2 = positionIndex[node];
		int ndx1 = positionIndex[node - 1];
		double distance = positions[ndx2] - positions[ndx1];
		distance = Math.max(distance, 10);
		double m = getTransitionProbability(state1, state2) * distance/avgDistance;
		double c = (1 - Math.exp(-2*m)) / 2;
		return Math.log(c);
	}
	
	@Override
	public void setObservation(int observation) {
		
		int thisPosition = positions[observation];
//		double dblPos = (double) thisPosition;
//		int ndx = Arrays.binarySearch(rates[0], dblPos);
//		if (ndx < 0) {
//			ndx = -(ndx + 1);
//			if (ndx == 0) relativeRecombination = rates[1][0] / averageRecombinationRate;
//			else if (ndx == numberOfRates) relativeRecombination = rates[1][numberOfRates - 1] / averageRecombinationRate;
//			else {
//				double pd = (dblPos - rates[0][ndx - 1])/(rates[0][ndx] - rates[0][ndx - 1]);
//				double thisRate = pd * rates[1][ndx] + (1 - pd) * rates[1][ndx - 1];
//				relativeRecombination = thisRate / averageRecombinationRate;
//			}
//		} else {
//			relativeRecombination = rates[1][ndx] / averageRecombinationRate;
//		}

	}

	@Override
	public int getNumberOfStates() {
		return 3;
	}

	/**
	 * @param filename	the name of the file containing recombination rates at positions along the chromosome
	 */
	public void readRates(String filename) {
		int count = 0;
		double sumRate = 0;
		try {
			BufferedReader br = new BufferedReader(new FileReader(filename));
			 //count the lines
			numberOfRates = 0;
			while (br.readLine() != null) numberOfRates++;
			br.close();
			
			rates = new double[2][numberOfRates - 1];
			br = new BufferedReader(new FileReader(filename));
			br.readLine(); //the header
			String input;
			while ((input = br.readLine()) != null) {
				String[] info = input.split("\t");
				rates[0][count] = Double.parseDouble(info[1]);
				rates[1][count] = Double.parseDouble(info[2]);
				sumRate += rates[1][count];
				count++;
			}
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
			
		}
		
		averageRecombinationRate = sumRate / count;
	}

	/**
	 * This function calculates a matrix of transition rates between marker states using all the transitions in the states array.
	 * @param states	an array of marker states.
	 */
	public double calculateAverageTransitionMatrix(byte[][] states) {
		int[][] transitionCounts = new int [3][3];
		
		for (byte[] taxonStates:states) {
			byte prevState = -1;
			for (byte state:taxonStates) {
				if (prevState == -1) {
					prevState = state;
				} else if (state != -1) {
					transitionCounts[prevState][state]++;
					prevState = state;
				}
			}
		}
		
		int[] rowSums = new int[3];
		int[] colSums = new int[3];
		int total = 0;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				rowSums[i] += transitionCounts[i][j];
				colSums[j] += transitionCounts[i][j];
				total++;
			}
		}
		
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				double top = transitionCounts[i][j];
				double bottom = rowSums[i];
				averageTransitionMatrix[i][j] = top / bottom;
 			}
		}
		System.out.println("Average transition matrix:");
		for (double[] row :averageTransitionMatrix) {
			for (double col:row) {
				System.out.print(col + "  ");
			}
			System.out.println();
		}
		System.out.println("Finished calculating average transition matrix.");
		
		int chrLength = positions[positions.length - 1] - positions[0];
		int ntaxa = states.length;
		avgDistance = ((double) chrLength) * ((double) ntaxa) / ((double) total);
		
		return (averageTransitionMatrix[1][0] + averageTransitionMatrix[1][2]) / (averageTransitionMatrix[0][2] + averageTransitionMatrix[2][0]);
	}
	
	/**
	 * Sets the index of taxon snps into the overall chromosome snps. This index is used by the ln transition probability function that takes node as an argument.
	 * @param index	the index of non missing snp positions for a taxon
	 */
	public void setPositionIndex(int[] index) {
		positionIndex = index;
	}
	
	/**
	 * This function calculates the average distance between non-missing loci.
	 * @param states	a matrix of marker states
	 */
	public void calculateAverageDistance(byte[][] states) {
		int totalCount = 0;
		double sumLengths = 0;
		int nsnps = states[0].length;
		for (byte[] taxonStates:states) {
			byte prevState = taxonStates[0];
			int prevSnp = 0;
			for (int s = 1; s < nsnps; s++) {
				byte state = taxonStates[s];
				if (prevState == -1) {
					prevState = state;
					prevSnp = s;
				} else if (state != -1) {
					totalCount++;
					sumLengths += positions[s] - positions[prevSnp];
					prevState = state;
					prevSnp = s;
				}
			}
		}
		avgDistance = sumLengths / totalCount;
	}
	
	public void setAverageDistance(double avgDistance) {
		this.avgDistance = avgDistance;
	}
	
	public void setAverageTransitionMatrix(double[][] atm) {
		this.averageTransitionMatrix = atm;
	}
}
