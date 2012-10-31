package net.maizegenetics.gwas.Imputation;

import java.util.ArrayList;
import java.util.HashMap;

public interface GenotypeData {
	/**
	 * @return the chromosome for data will be returned
	 */
	int getChromosome();
	
	/**
	 * set the chromosome for which data will be returned
	 */
	void setChromosome(int chromosome);
	
	/**
	 * @return a List of SnpInfo objects representing the snps in this data set
	 */
	ArrayList<SnpInformation> getSnpList();
	
	/**
	 * @return a List of taxa names of the taxa contained in this data set
	 */
	ArrayList<String> getTaxaList();
	
	/**
	 * @return	a mapping of taxa names to the their order in the snp.values array
	 */
	HashMap<String, Integer> getTaxaMap();
	
	/**
	 * @param taxon the taxon for which the call is to be returned
	 * @param snp the snp for which the call is to be returned
	 * @return the snp call
	 */
	byte getCall(String taxon, SnpInformation snp);
	
	/**
	 * @param taxaIndex	the index in the underlying array of the taxon
	 * @param snpIndex	the index in the underlying array of the snp
	 * @return the snp call
	 */
	byte getCall(int taxaIndex, int snpIndex);
	
	/**
	 * @param taxon the taxon name
	 * @return an array of snp calls for this taxon in the same order as the snps returned by getSnpList()
	 */
	byte[] getSnpsForTaxon(String taxon);
	
	/**
	 * @param taxonIndex the index of the taxon, i.e. its order in taxaList
	 * @return an array of snp calls for this taxon in the same order as the snps returned by getSnpList()
	 */
	byte[] getSnpsForTaxon(int taxonIndex);
	
	/**
	 * @param snp a SnpInfo object representing a snp
	 * @return an array of snp calls for this snp, in the same order as the taxa returned by getTaxaList()
	 */
	byte[] getSnpsForSnp(SnpInformation snp);
	
	/**
	 * @return the number of taxa in the data set
	 */
	int getNumberOfTaxa();
	
	/**
	 * @return	the number of snps in the data set
	 */
	int getNumberOfSnps();
	
	/**
	 * @return	the probabilities or frequencies associated with the snp states in this dataset
	 */
	public double[] getStateProbabilities();
}
