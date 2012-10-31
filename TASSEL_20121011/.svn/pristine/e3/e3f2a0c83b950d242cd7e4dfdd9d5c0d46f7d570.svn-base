package net.maizegenetics.gwas.Imputation;

public class TestTransitionMatrix implements TransitionMatrix {
	double recombinationRate;
	
	public TestTransitionMatrix(double recombination) {
		recombinationRate = recombination;
	}
	
	@Override
	public int getNumberOfStates() {
		return 3;
	}

	@Override
	public double getTransitionProbability(int state1, int state2) {
//		if (state1 == state2) return (1 - 2*recombinationRate);
//		else return 2 * recombinationRate;
		double prob = 0;
		switch (state1) {
		case 0:
			switch (state2) {
			case 0:
				prob = 1 - 2 * recombinationRate;
				break;
			case 1:
				prob = 0.3 * recombinationRate;
				break;
			case 2:
				prob = 1.7 * recombinationRate;
				break;
			}
			break;
		case 1:
			switch (state2) {
			case 0:
				prob = 5 * recombinationRate;
				break;
			case 1:
				prob = 1 - 10 * recombinationRate;
				break;
			case 2:
				prob = 5 * recombinationRate;
				break;
			}
			break;
		case 2:
			switch (state2) {
			case 0:
				prob = 1.7 * recombinationRate;
				break;
			case 1:
				prob = 0.3 * recombinationRate;
				break;
			case 2:
				prob = 1 - 2 * recombinationRate;
				break;
			}
			break;
		}
		return prob;
	}

	@Override
	public double getLnTransitionProbability(int state1, int state2) {
		return Math.log(getTransitionProbability(state1, state2));
	}

	@Override
	public void setObservation(int observation) {
		// do nothing
		
	}

	@Override
	public double getLnTransitionProbability(int state1, int state2, int node) {
		return getLnTransitionProbability(state1, state2);
	}

}
