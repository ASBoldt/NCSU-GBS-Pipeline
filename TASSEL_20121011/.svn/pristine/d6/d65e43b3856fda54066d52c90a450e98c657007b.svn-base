package net.maizegenetics.gwas.Imputation;

public class ProbObsGivenState extends ProbabilityObservationGivenTrueState {

	@Override
	public double probabilityAgivenB(int stateA, int stateB, int observation) {
		double pmono = probMonomorphic[pop][subpop][observationIndex[observation]];
		double ifpoly = gbsProbabilityMatrix[stateA][stateB];
		double ifmono;
		if (stateA == 0) ifmono = 1;
		else ifmono = 0;
		return ifpoly * (1 - pmono) + ifmono * pmono;
		
		
		//X in {A,B,H}, Y in {A,B,H}
		//P(obs=X|state=Y) = P(obs=X|state=Y, poly)P(poly) + P(obs=X|state=Y, monoA)P(monoA) [note that P(monoB = 0]]
		//P(obs=X|state=A, monoA) = 1; P(obs=X|state in {B,H}, monoA) = 0
		//P(obs=X|state=A, poly), use obsGivenTrue
		//P(obs=X|state=B, poly), use obsGivenTrue
		//P(obs=A|state=1A:3B, poly) = 0.25 * (1 - obsGivenTrue[1][1]) 
		//P(obs=A|state=1A:1B, poly) = 0.5 * (1 - obsGivenTrue[1][1]) 
		//P(obs=A|state=3A:1B, poly) = 0.75 * (1 - obsGivenTrue[1][1])
		//similarly for obs=B
		//P(obs=H|state=H, poly) = obsGivenTrue[1][1]
		//P(obs=A|state=H, poly) = P(obs=A|state=1A:3B)P(state=1A:3B| H,poly) + P(obs=A|state=1A:1B)P(state=1A:1B, poly) + P(obs=A|state=3A:1B)P(state=3A:1B, poly)
		// = 0.25*P(state=1A:3B| H,poly) + 0.5*P(state=1A:1B, poly) + 0.75*P(state=3A:1B, poly)
		//similarly for obs=B
		//the difficulty is estimating P(state = 1A:3B | H, poly), etc
	}

}
