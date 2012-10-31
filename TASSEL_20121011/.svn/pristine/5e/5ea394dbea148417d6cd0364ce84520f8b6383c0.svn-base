package net.maizegenetics.gwas.Imputation;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.regex.Pattern;

public class GbsGenotypeData implements GenotypeData {
	private int chromosome;
	private ArrayList<SnpInformation> mySnpList = new ArrayList<SnpInformation>();
	private ArrayList<String> myTaxaList;
	private HashMap<String, Integer> myTaxaMap = new HashMap<String, Integer>();
	
	private int numberOfTaxa;
	private int numberOfSnps;
	
	private final Pattern tab = Pattern.compile("\t");
	
	/** 
	 * @param filename	the name of the file containing GBS data in Hapmap format
	 */
	public GbsGenotypeData(String filename) {
		this(filename, false);
	}
	
	/**
	 * @param filename	the name of the file containing GBS data in Hapmap format
	 * @param setHetsToMissing	true, if hets should be set to missing, false otherwise
	 */
	public GbsGenotypeData(String filename, boolean setHetsToMissing) {
		//read gbs data
		System.out.println("Reading gbs data...");
		HashMap<String, Integer> gbsTaxaMap = new HashMap<String, Integer>();

		try {
			//read taxa and snp names
			BufferedReader br = new BufferedReader(new FileReader(filename));
			String input = br.readLine();
			String[] info = tab.split(input);
			int infolen = info.length;
			for (int i = 11; i < infolen; i++) {
				String taxaname = info[i];
				if (taxaname.startsWith("Z0")) {
					int colonpos = taxaname.indexOf(':');
					if (colonpos > 0) taxaname = new String(taxaname.substring(0, colonpos));
					gbsTaxaMap.put(taxaname, i);
				}
			}
			
			myTaxaList = new ArrayList<String>(gbsTaxaMap.keySet());
			Collections.sort(myTaxaList);
			int count = 0;
			for (String taxon:myTaxaList) {
				myTaxaMap.put(taxon,count++);
			}
			numberOfTaxa = myTaxaList.size();
			
			//an index of the saved taxa into the genotype data
			int[] taxonIndex = new int[numberOfTaxa];
			count = 0;
			for (String taxon:myTaxaList) {
				taxonIndex[count++] = gbsTaxaMap.get(taxon);
			}
			
			count = 0;
			while ((input = br.readLine()) != null) {
				info = tab.split(input);
				SnpInformation snp = new SnpInformation();
				snp.name = new String(info[0]);
				snp.chromosome = Integer.parseInt(info[2]);
				snp.position = Integer.parseInt(info[3]);
				snp.ndx = count++;
				snp.dataset = this;
				byte[] calls = new byte[numberOfTaxa];
				for (int t = 0; t < numberOfTaxa; t++) {
					byte byteval;
					String strval = info[taxonIndex[t]];
					if (strval.equals("A")) byteval = 0;
					else if (strval.equals("H")) {
						if (setHetsToMissing) {
							byteval = -1;
						} else {
							byteval = 1;
						}
					}
					else if (strval.equals("B")) byteval = 2;
					else byteval = -1;
					calls[t] = byteval;
				}
				snp.values = calls;
				mySnpList.add(snp);
			}
			br.close();
			numberOfSnps = mySnpList.size();

		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
		chromosome = mySnpList.get(0).chromosome;
	}
	
	@Override
	public int getChromosome() {
		return chromosome;
	}

	@Override
	public void setChromosome(int chromosome) {
		this.chromosome = chromosome;
	}

	@Override
	public ArrayList<SnpInformation> getSnpList() {
		return mySnpList;
	}

	@Override
	public ArrayList<String> getTaxaList() {
		return myTaxaList;
	}

	@Override
	public byte getCall(String taxon, SnpInformation snp) {
		return snp.values[myTaxaMap.get(taxon)];
	}

	@Override
	public byte getCall(int taxaIndex, int snpIndex) {
		return mySnpList.get(snpIndex).values[taxaIndex];
	}

	@Override
	public byte[] getSnpsForTaxon(String taxon) {
		int ndx = myTaxaMap.get(taxon);
		byte[] calls = new byte[numberOfSnps];
		for (int s = 0; s < numberOfSnps; s++) {
			calls[s] = mySnpList.get(s).values[ndx];
		}
		return calls;
	}

	@Override
	public byte[] getSnpsForTaxon(int taxonIndex) {
		byte[] calls = new byte[numberOfSnps];
		for (int s = 0; s < numberOfSnps; s++) {
			calls[s] = mySnpList.get(s).values[taxonIndex];
		}
		return calls;
	}

	@Override
	public byte[] getSnpsForSnp(SnpInformation snp) {
		return snp.values;
	}

	@Override
	public int getNumberOfTaxa() {
		return numberOfTaxa;
	}

	@Override
	public int getNumberOfSnps() {
		return numberOfSnps;
	}

	@Override
	public double[] getStateProbabilities() {
		double[] prob = new double[3];
		int[] stateCounts = new int[3];
		
		for (SnpInformation snp:mySnpList) {
			for (byte val:snp.values) {
				if (val != -1) stateCounts[val]++;
			}
		}
		
		int total = 0;
		for (int count:stateCounts) total += count;
		
		for (int i = 0; i < 3; i++) prob[i] = ((double) stateCounts[i]) / ((double) total);
				
		return prob;
	}

	@Override
	public HashMap<String, Integer> getTaxaMap() {
		return myTaxaMap;
	}
}
