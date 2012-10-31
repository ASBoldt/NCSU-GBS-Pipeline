package net.maizegenetics.gwas.Imputation;

public class ProbabilityObservationGivenTrueState implements
		ConditionalProbability {
	boolean[] isgbs;
	double[][] arrayProbabilityMatrix = new double[][]{{0.9998,0.0001,0.0001},{0.0001,0.9998,0.0001},{0.0001,0.0001,0.9998}}; 
	double[][] gbsProbabilityMatrix = new double[][]{{0.98,0.2,0.01},{0.01,.6,0.01},{0.01,0.2,0.98}};
	
	double[][][] probMonomorphic = null; //dimensions are pop, subpop, snp (observation number)
	double[][] gbsMonomorphicMatrix = null;
	int pop;
	int subpop;
	int[] observationIndex; //the index of the current set of observations into the overall snps
	
	@Override
	public double probabilityAgivenB(int stateA, int stateB, int observation) {
		if(isgbs[observation]) {
			if (gbsMonomorphicMatrix != null) {
				double pmono = probMonomorphic[pop][subpop][observationIndex[observation]];
				if (Double.isNaN(pmono)) {
//					System.out.println("pmono = NaN at obs " + observation + " for pop " + pop + " subpop " + subpop);
				}
				return (1 - pmono) * gbsProbabilityMatrix[stateA][stateB] + pmono * gbsMonomorphicMatrix[stateA][stateB];
//				double ifmono;
//				if (stateA == 0) ifmono = 1;
//				else ifmono = 0;
//				return (1 - pmono) * gbsProbabilityMatrix[stateA][stateB] + pmono * ifmono;
				
			}
			return gbsProbabilityMatrix[stateA][stateB];
		}
		
		return arrayProbabilityMatrix[stateA][stateB];
	}

	@Override
	public double lnProbabilityAgivenB(int stateA, int stateB, int observation) {
		return Math.log(probabilityAgivenB(stateA, stateB, observation));
	}

	public void setIsGBS(boolean[] isgbs) {
		this.isgbs = isgbs;
	}

	public void setArrayProbabilityMatrix(double[][] arrayProbabilityMatrix) {
		this.arrayProbabilityMatrix = arrayProbabilityMatrix;
	}

	public void setGbsProbabilityMatrix(double[][] gbsProbabilityMatrix) {
		this.gbsProbabilityMatrix = gbsProbabilityMatrix;
	}
	
	public void setProbMonomorphic(double[][][] prob) {
		probMonomorphic = prob;
	}

	public void setProbabilityOfAnObservation(double[] pobs) {
		gbsMonomorphicMatrix = new double[3][3];
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				gbsMonomorphicMatrix[i][j] = pobs[i];
			}
		}
	}
	
	public void setSubpop(int pop, int subpop) {
		this.pop = pop;
		this.subpop = subpop;
	}
	
	public void setObservationIndex(int[] obsIndex) {
		this.observationIndex = obsIndex;
	}
}

