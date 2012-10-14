package net.maizegenetics.gwas.Imputation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.ListIterator;
import java.util.TreeSet;
import java.util.regex.Pattern;

public class Recombination implements TransitionMatrix {
//	private int chromosome;
	enum MapFunction {Haldane, Kosambi};
	
	private double[][] rates;
	private int numberOfPoints = 100; //calculate rates at this many points
	private double window = .1; //the window for each point is this fraction of the chromosome
	private int[] popnums = new int[]{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,19,20,21,22,23,24,25,26};
	private double[] cumulativeRecombinations;
	private double[][] cumulativeTransitions;
	private MapFunction myMapFunction = MapFunction.Haldane;
	
	ArrayList<Integer> posList;
	
	private static final String[] byte2String = new String[]{"A", "H", "B", "N"};
	private HashMap<String, Byte> string2byteMap = new HashMap<String, Byte>();
	private static final byte A = 0;
	private static final byte H = 1;
	private static final byte B = 2;
	private static final byte N = 3;
	private double[][] transitionProbability;
	private int[] evaluationPositions;
	double recombinationRate;
	
	public Recombination() {
		string2byteMap.put("A", A);
		string2byteMap.put("H", H);
		string2byteMap.put("B", B);
		string2byteMap.put("N", N);
	}
	
	public void countCumulativeRecombinants() {
		class snpinfo {
			final int pos;
			final String geno;

			snpinfo(int pos, String geno) {
				this.pos = pos;
				this.geno = geno;
			}
		}

		Pattern tab = Pattern.compile("\t");
		String basename = "/Volumes/Macintosh HD 2/data/namgbs/cov20/cov20_chr";
		String readEnding = ".hets12.txt";
		String writeEnding = ".hets12.recomb.window.1.txt";
		
		
		//find all of the snp positions
		System.out.println("Reading snp positions...");
		for (int chr = 1; chr <= 10; chr++) {
			TreeSet<Integer> posSet = new TreeSet<Integer>();
			for (int pop = 0; pop < 25; pop++) {
				String znum;
				if (popnums[pop] < 10) znum = "Z00" + popnums[pop];
				else znum = "Z0" + popnums[pop];
				String filename = basename + chr + "_" + znum + readEnding;

				try {
					BufferedReader br = new BufferedReader(new FileReader(filename));
					String input = br.readLine();
					String[] info = tab.split(input);
					while((input = br.readLine()) != null) {
						info = tab.split(input, 5);
						posSet.add(Integer.decode(info[3]));
					}
					br.close();
				} catch (IOException e) {
					e.printStackTrace();
					System.exit(-1);
				}
			}

			posList = new ArrayList<Integer>(posSet);
			cumulativeRecombinations = new double[posList.size()];

			//for each population get a list of position, genotype pairs for all non-N snps
			int numberOfTaxa = 0;
			for (int pop = 0; pop < 25; pop++) {

				String znum;
				if (popnums[pop] < 10) znum = "Z00" + popnums[pop];
				else znum = "Z0" + popnums[pop];
				System.out.println("Reading snps and counting recombinants for chr " + chr + ", znum = " + znum + " ...");
				String filename = basename + chr + "_" + znum + readEnding;
				ArrayList<ArrayList<snpinfo>> taxainfo = new ArrayList<ArrayList<snpinfo>>();

				try {
					BufferedReader br = new BufferedReader(new FileReader(filename));
					String input = br.readLine();
					String[] info = tab.split(input);
					int ntaxa = info.length - 11;
					for (int t = 0; t < ntaxa; t++) taxainfo.add(new ArrayList<snpinfo>());
					while((input = br.readLine()) != null) {
						info = tab.split(input);
						int pos = Integer.parseInt(info[3]);
						for (int t = 0; t < ntaxa; t++) {
							String val = info[t + 11];
							if (!val.equals("N")) taxainfo.get(t).add(new snpinfo(pos, val));
						}
					}
					br.close();
				} catch (IOException e) {
					e.printStackTrace();
					System.exit(-1);
				}

				//use taxa info to add counts for each taxon to total recombinations
				numberOfTaxa += taxainfo.size();
				Integer snpPos;
				int snpIndex = 0;
				int numberOfPositions = posList.size();
				for (ArrayList<snpinfo> taxon:taxainfo) {
					int recombinants = 0;
//					int[][] transitions = new int[3][3];
					ListIterator<Integer> it = posList.listIterator();
					snpIndex = 0;
					snpPos = it.next();
					snpinfo prevSnp = null;
					for (snpinfo thisSnp : taxon) {
						if (prevSnp == null || prevSnp.geno.equals(thisSnp.geno)) {
							//iterate through posList until pos is less than or equal to this position
							//add number of recombinants to each pos 
							while (snpPos > -1 && snpPos <= thisSnp.pos) {
								cumulativeRecombinations[snpIndex] += recombinants;
//								for (int i = 0; i < 3; i++) {
//									for (int j = 0; j < 3; j++) {
//										cumulativeTransitions[i][j] += transitions[i][j];
//									}
//								}
								if (it.hasNext()) {
									snpPos = it.next();
									snpIndex++; 
								} else snpPos = -1;
							}
						} else {
							//iterate through posList until pos is less than or equal to this position
							//add a number number of recombinants plus a number between 0 and 1 proportional to distance at each position
							double dist = thisSnp.pos - prevSnp.pos;
							while (snpPos > -1 && snpPos <= thisSnp.pos) {
								double pd = snpPos - prevSnp.pos;
								cumulativeRecombinations[snpIndex] += recombinants + pd / dist;
								if (it.hasNext()) {
									snpPos = it.next();
									snpIndex++; 
								} else snpPos = -1;
							}
							recombinants++;
						}
						prevSnp = thisSnp;
					}
					//fill in values for the remaining positions at the end of the chromosome
					while (snpIndex < numberOfPositions) {
						cumulativeRecombinations[snpIndex++] += recombinants;
					}
				}

			}
			calculateRates(numberOfTaxa);
			String writeFilename = basename + "_chr" + chr + ".recombinationrate.het12.win1.txt";
			writeRates(writeFilename, chr);
		}

	}
	
	public void calculateRates(int numberOfTaxa) {
		//physical length is first pos - last pos
		//this needs to calculate slope by finding position, recombination pairs in each interval across all populations
		int numberOfPos = posList.size();
		int lastpos = posList.get(posList.size() - 1);
		int physlength = lastpos - posList.get(0);
		int firstPosition = (int) (physlength * window / 2);
		int lastPosition = lastpos - firstPosition;
		int interval = (lastPosition - firstPosition) / (numberOfPoints - 1);
//		rates = new double[numberOfPoints][2];
		rates = new double[2][numberOfPoints];
		for (int i = 0; i < numberOfPoints; i++) {
			int thisPos = firstPosition + i * interval;
			double thisSlope = slope(thisPos - firstPosition, thisPos + firstPosition);
			rates[0][i] = thisPos;
			rates[1][i] = thisSlope / numberOfTaxa;
		}
		
	}
	
	private double slope(int start, int end) {
		//this needs to calculate slope by finding position, recombination pairs in each interval across all populations
		//sums of squares and products
		double sspos = 0;
		double sumpos = 0;
		double sumrecomb = 0;
		double sumprod = 0;
		int datacnt = 0;
		
		int ndx = Collections.binarySearch(posList, start);
		if (ndx < 0) ndx = -(ndx + 1);
		
		double thisPos = posList.get(ndx);
		double cr = cumulativeRecombinations[ndx];
		while (thisPos <= end) {
			double squarePos = thisPos *thisPos;
			sspos += thisPos * thisPos;
			sumpos += thisPos;
			sumrecomb += cr;
			sumprod += thisPos * cr;
			datacnt++;
			if (ndx == posList.size() - 1) break;
			thisPos = posList.get(++ndx);
			cr = cumulativeRecombinations[ndx];
		}


		double n = (double) datacnt;
		return (sumprod - sumpos * sumrecomb / n) / (sspos - sumpos * sumpos / n);
	}
	
	public void writeRates(String filename, int chromosome) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(filename));
			bw.write("chromosome\tposition\trate\n");
			for (int i = 0; i < numberOfPoints; i++) {
				bw.write(Integer.toString(chromosome));
				bw.write("\t");
				bw.write(rates[0][i] + "\t");
				bw.write(rates[1][i] + "\n");
			}
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void readRates(String filename) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(filename));
			 //count the lines
			int nlines = 0;
			while (br.readLine() != null) nlines++;
			br.close();
			
			rates = new double[2][nlines - 1];
			br = new BufferedReader(new FileReader(filename));
			br.readLine(); //the header
			String input;
			int count = 0;
			while ((input = br.readLine()) != null) {
				String[] info = input.split("\t");
				rates[0][count] = Double.parseDouble(info[1]);
				rates[1][count] = Double.parseDouble(info[2]);
				count++;
			}
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
			
		}
	}
	
	public double getRateAtAPosition(double position) {
		int ndx = Arrays.binarySearch(rates[0], position);
		int n = rates[0].length;
		if (ndx >= n) {
			return rates[1][n - 1];
		}
		if (ndx >= 0) {
			return rates[1][ndx];
		}
		
		ndx = -(ndx + 1);
		if (ndx == 0) {
			return rates[1][0];
		}
		if (ndx >= n) {
			return rates[1][n - 1];
		}
		
		double partialDistance = (position - rates[0][ndx - 1]) / (rates[0][ndx] - rates[0][ndx - 1]);
		return rates[1][ndx - 1] + partialDistance * (rates[1][ndx] - rates[1][ndx - 1]);
	}
	
	public double recombinationProbability(double start, double end) {
		double midpoint = (start + end) / 2;
		double rate = getRateAtAPosition(midpoint);
		return rate * (end - start + 1);
	}

	@Override
	public double getTransitionProbability(int state1, int state2) {
		if (state1 == 0) {
			if (state2 == 0) {
				return 1 - recombinationRate;
			} else if (state2 == 1) {
				return 0.15 * recombinationRate;
			} else {
				return 0.85 * recombinationRate;
			}
		} else if (state1 == 1) {
			if (state2 == 0) {
				return 2.5 * recombinationRate;
			} else if (state2 == 1) {
				return 1 - 5 * recombinationRate;
			} else {
				return 2.5 * recombinationRate;
			}
		} else {
			if (state2 == 0) {
				return 0.85 * recombinationRate;
			} else if (state2 == 1) {
				return 0.15 * recombinationRate;
			} else {
				return 1 - recombinationRate;
			}
		}
//		return transitionProbability[state1][state2];
	}

	public void setObservation(int observation) {
		recombinationRate = recombinationProbability(evaluationPositions[observation - 1], evaluationPositions[observation]);
		//apply the map function
		if (myMapFunction == MapFunction.Haldane) {
			recombinationRate = (1 - Math.exp(-2*recombinationRate))/2;
		} else if (myMapFunction == MapFunction.Kosambi) {
			double val = Math.exp(4 * recombinationRate);
			recombinationRate = (val - 1) / 2 / (val + 1);
		}
	}
	
	@Override
	public double getLnTransitionProbability(int state1, int state2) {
		return Math.log(getTransitionProbability(state1, state2));
	}

	public void setEvaluationPositions(int[] positions) {
		evaluationPositions = positions;
	}

	@Override
	public int getNumberOfStates() {
		return 3;
	}
	
	public void setTransitionProbability(double[][] transitionProbability) {
		this.transitionProbability = transitionProbability;
	}
	
	public double getChrLength() {
		int n = rates[0].length;
		return rates[0][n-1];
	}
	
	public double getRecombinationRate(int observation) {
		setObservation(observation);
		return recombinationRate;
	}

	@Override
	public double getLnTransitionProbability(int state1, int state2, int node) {
		setObservation(node);
		return getLnTransitionProbability(state1, state2);
	}
}
