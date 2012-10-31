package net.maizegenetics.gwas.Imputation;

public class VariableTransitionMatrix5State extends
		VariableAverageTransitionMatrix {
	
	private final int[] positions;
	int numberOfStates = 5;	
	
	public VariableTransitionMatrix5State(int chromosome, int[] positions) {
		super(chromosome, positions);
		this.positions = positions;
	}

	/**
	 * This function calculates a matrix of transition rates between marker states using all the transitions in the states array.
	 * @param states	an array of marker states.
	 */
	public double calculateAverageTransitionMatrix(byte[][] states) {
		int[][] transitionCounts = new int [numberOfStates][numberOfStates];
		double[][] atm = new double[numberOfStates][numberOfStates];
		
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
		
		int[] rowSums = new int[numberOfStates];
		int[] colSums = new int[numberOfStates];
		int total = 0;
		for (int i = 0; i < numberOfStates; i++) {
			for (int j = 0; j < numberOfStates; j++) {
				rowSums[i] += transitionCounts[i][j];
				colSums[j] += transitionCounts[i][j];
				total++;
			}
		}
		
		for (int i = 0; i < numberOfStates; i++) {
			for (int j = 0; j < numberOfStates; j++) {
				double top = transitionCounts[i][j];
				double bottom = rowSums[i];
				atm[i][j] = top / bottom;
 			}
		}
		
		super.setAverageTransitionMatrix(atm);
		
		System.out.println("Average transition matrix:");
		for (double[] row :atm) {
			for (double col:row) {
				System.out.print(col + "  ");
			}
			System.out.println();
		}
		System.out.println("Finished calculating average transition matrix.");
		
		int chrLength = positions[positions.length - 1] - positions[0];
		int ntaxa = states.length;
		double avgdist = ((double) chrLength) * ((double) ntaxa) / ((double) total);
		super.setAverageDistance(avgdist);
		
		return 1;
	}
	
	@Override
	public int getNumberOfStates() {
		return numberOfStates;
	}


}
