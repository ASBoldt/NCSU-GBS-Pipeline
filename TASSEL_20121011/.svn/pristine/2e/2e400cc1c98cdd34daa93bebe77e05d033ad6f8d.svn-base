package net.maizegenetics.gwas.jointlinkage;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import net.maizegenetics.jGLiM.BasicShuffler;
import net.maizegenetics.jGLiM.LinearModelUtils;
import net.maizegenetics.jGLiM.dm.FactorModelEffect;
import net.maizegenetics.jGLiM.dm.ModelEffect;
import net.maizegenetics.jGLiM.dm.ModelEffectUtils;
import net.maizegenetics.jGLiM.dm.NestedCovariateModelEffect;
import net.maizegenetics.jGLiM.dm.SweepFastLinearModel;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory;

public class PermutationTesterForNAMJointLinkage {
	double[] data;
	ArrayList<double[]> permutedData = new ArrayList<double[]>();
	String[] pops;
	SweepFastLinearModel sflm;
	ArrayList<ModelEffect> model = new ArrayList<ModelEffect>();
	FactorModelEffect popEffect;
	SNPdata snpdata;
	int[][] snpIndices;
	ExecutorService exec;
	NamMap nammap;
	int numberOfPermutations = 10;
	String outputFile = "/Users/peter/temp/permutationTestOutput.txt";
	String genotypeFilename = "";
	String phenotypeFilename = "";
	
	public PermutationTesterForNAMJointLinkage(int numberOfPermutations) {
		this.numberOfPermutations = numberOfPermutations;
		loadData();
		runPermutationTest_noMissingData();

	}
	
	private void loadData() {
		System.out.println("loading data...");
		int phenoColumn = 3;
		snpdata = new SNPdataNAMGBS(phenoColumn, genotypeFilename, phenotypeFilename);
//		snpdata = new SNPdataFromArray();
		data = snpdata.getPhenotype();
		pops = snpdata.getPopulations();
	}
	
	public void runPermutationTest_noMissingData() {
		System.out.println("Running permutations...");

		int numberOfObs = data.length;
		double totalSS = 0;
		for (int i = 0; i < numberOfObs; i++) totalSS += data[i] * data[i];
		
		ArrayList<double[]> baseModelSSdf = new ArrayList<double[]>();
		
		int[] mean = new int[numberOfObs];
		FactorModelEffect meanME = new FactorModelEffect(mean, false, "mean");
		FactorModelEffect popME = new FactorModelEffect(ModelEffectUtils.getIntegerLevels(pops), true);
		ArrayList<ModelEffect> effectList = new ArrayList<ModelEffect>();
		effectList.add(meanME);
		effectList.add(popME);
		DoubleMatrix[][] Xmatrices = new DoubleMatrix[1][3];
		Xmatrices[0][0] = meanME.getX();
		Xmatrices[0][1] = popME.getX();
		
		//permute data
		double[] minP = new double[numberOfPermutations];
		for (int p = 0; p < numberOfPermutations; p++) {
			minP[p] = 1;
			double[] pdata = new double[numberOfObs];
			System.arraycopy(data, 0, pdata, 0, numberOfObs);
			BasicShuffler.shuffle(pdata);
			permutedData.add(pdata);
			
			SweepFastLinearModel sflm = new SweepFastLinearModel(effectList, pdata);
			baseModelSSdf.add(sflm.getFullModelSSdf());
		}
		
		for (int m = 0; m < snpdata.getNumberOfSNPs(); m++) {
			if (m % 50 == 0) System.out.println("Testing marker " + m);
			SNP snp = snpdata.getSnp(m);
			float[] snpscore = snp.score;
			double[] dscore = new double[numberOfObs];
			for (int i = 0; i < numberOfObs; i++) {
				dscore[i] = snpscore[i];
			}
			
			NestedCovariateModelEffect ncme = new NestedCovariateModelEffect(dscore, popME);
			
			Xmatrices[0][2] = ncme.getX();
			DoubleMatrix X = DoubleMatrixFactory.DEFAULT.compose(Xmatrices);
			
			effectList.add(ncme);
			sflm = new SweepFastLinearModel(effectList, data);
			double[] modelSSdf = sflm.getFullModelSSdf();
			DoubleMatrix G = sflm.getInverseOfXtX();
			
			for (int p = 0; p < numberOfPermutations; p++) {
				double[] pdata = permutedData.get(p);
				DoubleMatrix y = DoubleMatrixFactory.DEFAULT.make(numberOfObs, 1, pdata);
				DoubleMatrix Xty = X.crossproduct(y);
				double[] reducedSSdf = baseModelSSdf.get(p);
				double fullSS = Xty.crossproduct(G.mult(Xty)).get(0, 0);
				double fulldf = modelSSdf[1];
				double markerSS = fullSS - reducedSSdf[0];
				double markerdf = fulldf - reducedSSdf[1];
				double errorSS = totalSS - fullSS;
				double errordf = numberOfObs - fulldf;
				double F = markerSS / markerdf / errorSS * errordf;
				double probF;
				try {
					probF = LinearModelUtils.Ftest(F, markerdf, errordf);
				} catch(Exception e) {
					probF = 1;
				}
				
				minP[p] = Math.min(minP[p], probF);
				
			}
			
			effectList.remove(ncme);
		}
		
		writeResultsToFile(minP);
		System.out.println("Finished writing results to " + outputFile);

	}
	
	public void writeResultsToFile(double[] minP) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile));
			for (double p : minP) bw.write(p + "\n");
			bw.close();
		} catch(IOException e) {
			e.printStackTrace();
		}
	}
	
	public BufferedWriter openlog() {
		BufferedWriter bw = null;
		try {
			bw = new BufferedWriter(new FileWriter("c:/users/peter/temp/permutation_log.txt"));
		} catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return bw;
	}
	
	public void writeToLog(String output, BufferedWriter bw) {
		try {
			bw.write(output);
			bw.newLine();
		} catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	public void closewriter(BufferedWriter bw) {
		try {
			bw.close();
		} catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
}
