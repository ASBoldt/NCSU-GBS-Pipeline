package net.maizegenetics.gwas.jointlinkage;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.regex.Pattern;

public class SNPdataNAMGBS implements SNPdata {
	final private ArrayList<SNP> snpList = new ArrayList<SNP>();
	private String phenotypeName;
	private int phenotypeColumn;
	private String snpType = "NAM GBS SNPs";
//	final private String inputGenotypes = "/Volumes/Macintosh HD 2/data/namgbs/cov20/imputedFromArraySnpsAllChr.txt";
//	final private String inputGenotypes = "C:/Projects/NAM/namgbs/datasets/cov20/imputedWithoutHetsSnpsAllChr.txt";
//	final private String inputGenotypes = "/Users/pbradbury/Documents/namgbs/snpdata/cov20/imputedFromHetsSnpsAllChr.txt";
//	final private String inputGenotypes = "/Volumes/Macintosh HD 2/data/namgbs/cov20/imputedFromHets12SnpsAllChr.txt";
//	final private String inputPhenotypes = "/Volumes/Macintosh HD 2/data/namgbs/ImputedMarkerGenotypes_flowering_traits_092909.txt";
//	final private String inputPhenotypes = "/Volumes/Macintosh HD 2/data/lipka/BLUPs_Three_Carotenoids.txt";
//	final private String inputPhenotypes = "/Volumes/Macintosh HD 2/data/lipka/BLUPs_Remaining_Carotenoids.txt";
//	final private String inputPhenotypes = "/Volumes/Macintosh HD 2/data/lipka/BLUPs_Tocochromanols.txt";
	
	private int startSNP = 0;
	private int endSNP = 0;
	private int currentSNP = 0;
	private int numberOfSNPs = 0;
	private double[] phenotype;
	private String[] populations;
	private Pattern tab = Pattern.compile("\t");

	public SNPdataNAMGBS(int phenotypeCol, String inputGenotypes, String inputPhenotypes) {
		phenotypeColumn = phenotypeCol;
		
		//get a list of taxa with both genotypes and phenotypes, sort it
		HashMap<String, Integer> taxaMap = new HashMap<String, Integer>();
		HashMap<String, Double> phenoMap = new HashMap<String, Double>();
		String input;
		String[] info;
		int numberOfTaxa;
		try {
			BufferedReader br = new BufferedReader(new FileReader(inputGenotypes));
			info = tab.split(br.readLine());
			int n = info.length;
			for (int i = 5; i < n; i++) taxaMap.put(info[i], i);
			br.close();
			
			br = new BufferedReader(new FileReader(inputPhenotypes));
			info = tab.split(br.readLine());
			phenotypeName = info[phenotypeCol];
			while ((input = br.readLine()) != null) {
				info = tab.split(input);
				if (phenotypeColumn < info.length) {
					String ph = info[phenotypeColumn];
					String taxon = info[0];
					if (!ph.startsWith("-99") && taxaMap.containsKey(taxon)) {
						Double phenovalue;
						try {
							phenovalue = Double.valueOf(ph);
							phenoMap.put(taxon, phenovalue);
						} catch(Exception e) {

						}
					}
				}
			}
				
			br.close();
			
		} catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
		//remove taxa with no phenotype
		ArrayList<String> taxaWithoutPhenotypes = new ArrayList<String>();
		for (String t:taxaMap.keySet()) {
			if (phenoMap.get(t) == null) taxaWithoutPhenotypes.add(t);
		}
		for (String t:taxaWithoutPhenotypes) taxaMap.remove(t);
		
		numberOfTaxa = taxaMap.size();
	
		//import genotypes with values in taxa order
		int[] taxaIndex = new int[numberOfTaxa];
		ArrayList<String> taxaList = new ArrayList<String>(taxaMap.keySet());
		Collections.sort(taxaList);
		for (int i = 0; i < numberOfTaxa; i++) {
			taxaIndex[i] = taxaMap.get(taxaList.get(i));
		}
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(inputGenotypes));
			br.readLine();
			int count = 0;
			while((input = br.readLine()) != null) {
				info = tab.split(input);
				float[] geno = new float[numberOfTaxa];
				for (int t = 0; t < numberOfTaxa; t++) {
					geno[t] = Float.parseFloat(info[taxaIndex[t]]);
				}
				String name = info[0];
				String allele = info[1];
				int chr = Integer.parseInt(info[2]);
				int phypos = Integer.parseInt(info[3]);
				float genpos = Float.parseFloat(info[4]);
				snpList.add(new SNP(name, allele, chr, genpos, phypos, geno, count++));
			}
			br.close();
		} catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
		numberOfSNPs = snpList.size();
		endSNP = numberOfSNPs - 1;
		
		//import phenotypes with values in taxa order
		phenotype = new double[numberOfTaxa];
		for (int t = 0; t < numberOfTaxa; t++) {
			phenotype[t] = phenoMap.get(taxaList.get(t));
		}
		
		//populate populations
		populations = new String[numberOfTaxa];
		for (int t = 0; t < numberOfTaxa; t++) {
			populations[t] = taxaList.get(t).substring(0, 4);
		}
	}
	
	@Override
	public boolean hasNext() {
		return currentSNP < endSNP;
	}

	@Override
	public synchronized SNP next() {
		if (currentSNP == endSNP) return null;
		return snpList.get(++currentSNP);
	}

	@Override
	public synchronized SNP getSnp(int index) {
		return snpList.get(index);
	}

	@Override
	public int getNumberOfSNPs() {
		return numberOfSNPs;
	}

	@Override
	public void setStartingSNP(int start) {
		startSNP = start;
	}

	@Override
	public void setEndingSNP(int end) {
		endSNP = end;
	}

	@Override
	public int getChromosome(int index) {
		return snpList.get(index).chr;
	}

	@Override
	public int getPosition(int index) {
		return snpList.get(index).physicalPos;
	}

	@Override
	public void resetSNPs() {
		currentSNP = startSNP - 1;
	}

	@Override
	public double[] getPhenotype() {
		return phenotype;
	}

	@Override
	public String[] getPopulations() {
		return populations;
	}

	@Override
	public String getPhenotypeName() {
		return phenotypeName;
	}
	
}
