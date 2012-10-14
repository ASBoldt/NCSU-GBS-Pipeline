package net.maizegenetics.gwas.Imputation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.TreeSet;
import java.util.regex.Pattern;

public class ArrayGenotypeData implements GenotypeData {
	private int chromosome;
	private int numberOfTaxa;
	private int numberOfSnps;
	private ArrayList<SnpInformation> mySnpList;
	private ArrayList<String> myTaxaList;
	private HashMap<String, Integer> myTaxaMap;
	private int[] chrEnd;
	
	private final Pattern tab = Pattern.compile("\t");
	
	public ArrayGenotypeData() {
		
	}
	
	private void readGenotypes() {
		String filename = "/Volumes/Macintosh HD 2/data/namgbs/markers061208agpv2.txt";
		myTaxaList = new ArrayList<String>();
		mySnpList = new ArrayList<SnpInformation>();
		chrEnd = new int[11];
		chrEnd[0] = -1;
		try {
			BufferedReader br = new BufferedReader(new FileReader(filename));
			String input = br.readLine();
			String[] info = tab.split(input);
			int infolen = info.length;
			for (int i = 11; i < infolen; i++) {
				myTaxaList.add(info[i]);
			}
			numberOfTaxa = infolen - 11;
			int snpcount = 0;
			while ((input = br.readLine()) != null) {
				SnpInformation snp = new SnpInformation();
				snp.chromosome = Integer.parseInt(info[2]);
				chrEnd[snp.chromosome] = snpcount++;
				snp.position = Integer.parseInt(info[3]);
				snp.dataset = this;
				byte[] snpvals = new byte[numberOfTaxa];
				for (int t = 0; t < numberOfTaxa; t++) snpvals[t] = Byte.parseByte(info[t + 11]);
				snp.values = snpvals;
				mySnpList.add(snp);
			}
			numberOfSnps = mySnpList.size();
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
		myTaxaMap = new HashMap<String, Integer>();
		int count = 0;
		for (String taxon : myTaxaList) {
			myTaxaMap.put(taxon, count++);
		}
		
	}
	
	@Override
	public int getChromosome() {
		return chromosome;
	}

	@Override
	public void setChromosome(int chromosome) {
		this.chromosome = chromosome;
		numberOfSnps = chrEnd[chromosome] - chrEnd[chromosome - 1];
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
			calls[s] = mySnpList.get(s + chrEnd[chromosome - 1] + 1).values[ndx];
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

	public static void createArrayGenotypeFile() {
		//read array marker physical positions
		final String mapFilename = "/Volumes/Macintosh HD 2/data/namgbs/markers061208agpv2.txt";
		HashMap<String, int[]> arrayMap = new HashMap<String, int[]>();
		
		Pattern tab = Pattern.compile("\t");
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(mapFilename));
			String input = br.readLine();
			while ((input = br.readLine()) != null) {
				String[] info = tab.split(input);
				if (info.length > 6 && info[6].equals("estimated")) {
					//do not use this marker
				} else {
					int chr = Integer.parseInt(info[1]);
					int pos = Integer.parseInt(info[5]);
					arrayMap.put(info[0], new int[]{chr,pos});
				}
  			}
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}

		//read genotypes for each pop
		ArrayList<String[]> taxanameList = new ArrayList<String[]>();
		ArrayList<HashMap<String, byte[]>> genotypeList = new ArrayList<HashMap<String, byte[]>>();
		int[] pop = new int[]{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,19,20,21,22,23,24,25,26};
		for (int p = 0; p < 25; p++) {
			String filename = "/Volumes/Macintosh HD 2/data/namgbs/nam1106/ICIMdataForPopulation" + pop[p] + ".txt";
			String input;
			String[] info;
			HashMap<String, byte[]> popgeno = new HashMap<String, byte[]>();
			genotypeList.add(popgeno);
			try {
				BufferedReader br = new BufferedReader(new FileReader(filename));
				input = br.readLine();
				String[] taxanames = tab.split(input);
				int ntaxa = taxanames.length;
				taxanameList.add(taxanames);
				while ((input = br.readLine()) != null) {
					info = tab.split(input);
					byte[] geno = new byte[ntaxa];
					for (int t = 0; t < ntaxa; t++) {
						geno[t] = (byte) Integer.parseInt(info[t+1]);
					}
					popgeno.put(info[0], geno);
				}
				
				br.close();
			} catch (IOException e) {
				e.printStackTrace();
				System.exit(-1);
			}
		}
		
		//create a single array with all genotypes for all pops
		//make a set of distinct markers
		TreeSet<String> markerSet = new TreeSet<String>();
		for (HashMap<String, byte[]> geno : genotypeList) {
			markerSet.addAll(geno.keySet());
		}
		
		//write the data in hapmap format
		String header = "rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode";
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter("/Volumes/Macintosh HD 2/data/namgbs/nam1106/ArrayDataNAM.hmp.txt"));
			bw.write(header);
			ArrayList<Integer> taxaCounts = new ArrayList<Integer>();
			for (String[] names:taxanameList) {
				taxaCounts.add(names.length);
				for (String taxon:names) {
					bw.write("\t");
					bw.write(taxon);
				}
			}
			bw.newLine();
			
			for (String markerName : markerSet) {
				int[] positionInfo = arrayMap.get(markerName);
				bw.write(markerName);
				bw.write("\t0/2\t");
				bw.write(Integer.toString(positionInfo[0]));
				bw.write("\t");
				bw.write(Integer.toString(positionInfo[1]));
				for (int p = 0; p < 25; p++) {
					byte[] genovals = genotypeList.get(p).get(markerName);
					
					int n = taxaCounts.get(p);
					if (genovals == null) {
						for (int i = 0; i < n; i++) {
							bw.write("\t-1");
						}
					} else {
						if (n != genovals.length) System.err.println("genovals != n at marker " + markerName + " in pop " + p);

						for (int i = 0; i < n; i++) {
							bw.write("\t");
							bw.write(Byte.toString(genovals[i]));
						}
					}
					
				}
			}
			
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
	}

	@Override
	public HashMap<String, Integer> getTaxaMap() {
		return myTaxaMap;
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
	public byte[] getSnpsForTaxon(int taxonIndex) {
		// TODO Auto-generated method stub
		return null;
	}
}
