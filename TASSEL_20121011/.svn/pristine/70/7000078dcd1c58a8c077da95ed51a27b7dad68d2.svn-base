package net.maizegenetics.gwas.Imputation;

public interface TransitionMatrix {
	
	/**
	 * @param state1	the state of the a snp
	 * @param state2	the state of the next snp
	 * @return	the probability of a transition from state1 to state2
	 */
	double getTransitionProbability(int state1, int state2);
	
	/**
	 * @param state1	the state of the a locus
	 * @param state2	the state of the next locus
	 * @return	the natural log of the probability of a transition from state1 to state2
	 */
	double getLnTransitionProbability(int state1, int state2);
	
	/**
	 * @param state1	the state of the a locus
	 * @param state2	the state of the next locus
	 * @param node	the node or number of locus 2 (the locus with state2) in the sequence
	 * @return	the natural log of the probability of a transition from state1 to state2 at this node
	 */
	double getLnTransitionProbability(int state1, int state2, int node);
	
	/**
	 * @param observation	the observation at which the transition is to be evaluated. This will be the relative position of the state2 locus.
	 */
	void setObservation(int observation);
	
	/**
	 * @return	the number of possible states for a locus, e.g. three for a diploid, bi-allelic locus
	 */
	int getNumberOfStates();
}
