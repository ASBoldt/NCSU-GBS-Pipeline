package net.maizegenetics.gwas.Imputation;

public class ProbabilityObsGivenStateBiparental implements
		ConditionalProbability {
	
	//B in rows, A in columns
//	double[][] probabilityObsGivenState = new double[][]{{.998, .001,.001},{.4,.2,.4},{.001,.001,.998}};
//	double[][] probabilityObsGivenState = new double[][]{{1, 0, 0},{0, 1, 0},{0,0,1}};
	double[][] probabilityObsGivenState = new double[][]{{.998, .001,.0001},{.25,.5,.25},{.001,.001,.998}};
	
	@Override
	public double probabilityAgivenB(int stateA, int stateB, int observation) {
		return probabilityObsGivenState[stateB][stateA];
	}

	@Override
	public double lnProbabilityAgivenB(int stateA, int stateB, int observation) {
		return Math.log(probabilityAgivenB(stateA, stateB, observation));
	}

	public void setProbabilityMatrix(double[][] probabilities) {
		probabilityObsGivenState = probabilities;
	}
}
