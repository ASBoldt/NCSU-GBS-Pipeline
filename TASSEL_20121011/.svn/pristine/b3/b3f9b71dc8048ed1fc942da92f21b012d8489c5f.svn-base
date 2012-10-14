/*
 * Input ReadsByTaxa for counting (Taxa Number, total read number in each taxa, number of each read in each taxa) and sequences
 * Input CommonReadsOnMap for getting position
 * Output Fasta sequences for Blast after setting up a cutoff repetitive number
 * Output normalized number of each read in each taxa(by numner of each read / total read number in each taxa)
 * Output reads counts distribution (optional)
 *
 * Input genotype data of IBM population
 * Output genotype value table and phenotype value table for tassel analysis (this is limited by tassel's small throughput, so this step becomes optional)
 * Output genotype value table and phenotype value table for SAS analysis
 * Output genotype and phenotype value (2 into 1) for SAS analysis
 * 
 */

package net.maizegenetics.genome.GBS.repeat;

import net.maizegenetics.genome.GBS.GBSPipelineFei;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import net.maizegenetics.genome.BaseEncoder;
/**
 *
 * @author fl262
 */
public class setQuery {
	int totalReadNum;
	String[] taxaName;
	int repeatCutoff = 500;
	ReadWCountPosition[] RWCP;
	static String desDir = GBSPipelineFei.parentDir.replaceFirst("\\w+/$", "repeat/query/");
	static String freDesDir = desDir.replaceFirst("\\w+/$", "4QTL/");
	static String forTasselDir = freDesDir + "R_4Tassel/";
	static String forSASDir = freDesDir + "R_4SAS/";
	static File sourceGenoFile644 = new File(GBSPipelineFei.parentDir + GBSPipelineFei.childDir[11]).listFiles()[0];
	static File sourceGenoFile55k = new File("E:/Database/InfoFile/SNP55hapmap05162010_del191.txt");
	static File sourceGenoFileMerge = new File("E:/repeat/4QTL/R_4SAS/marker_selection/644plus55k.txt");

	public setQuery () {
		String dir, sourceFile, sourceFile1, sourceFile2;
		dir = GBSPipelineFei.parentDir+GBSPipelineFei.childDir[9] + "/";
		sourceFile1 = dir + new File(dir).list()[0];
		dir = GBSPipelineFei.parentDir+GBSPipelineFei.childDir[13] + "/";
		sourceFile2 = dir + new File(dir).list()[0];
		getCount(sourceFile1, sourceFile2);
		dir = GBSPipelineFei.parentDir+GBSPipelineFei.childDir[10] + "/";
		sourceFile = dir + new File(dir).list()[0];
		getPosition(sourceFile);
		Arrays.sort(RWCP);
		writeQuery(desDir, repeatCutoff);
		//writeTotalCountDis (); //optional choice
	}

	public void writeTable4QTL_SAS (String forSASDir, int hitOfCutoff) {
		GeneticInfoArray gif = new GeneticInfoArray(sourceGenoFileMerge);
		gif.setGenoValue();
		new File (forSASDir).mkdir();
		String genoPhenoFile = forSASDir + "genoPhenofile.txt";
		try {
			PrintWriter pw = new PrintWriter(genoPhenoFile);
			int effectTaxaNum = 0;
			int phenoNum = totalReadNum - hitOfCutoff;
			ArrayList<Integer> phenoIndex = new ArrayList<Integer>();
			ArrayList<Integer> genoIndex = new ArrayList<Integer>();
			for (int i = 0; i < taxaName.length; i++ ) {
				for (int j = 0; j < gif.taxaNum; j++) {
					if (taxaName[i].equals(gif.taxaName[j])) {
						phenoIndex.add(i);
						genoIndex.add(j);
						effectTaxaNum++;
						break;
					}
				}
			}
			Integer[] phenoIndexArray = phenoIndex.toArray(new Integer[phenoIndex.size()]);
			Integer[] genoIndexArray = genoIndex.toArray(new Integer[genoIndex.size()]);
			//transposition 1
			pw.printf("G/P"+"\t");
			for (int i = 0; i < effectTaxaNum - 1; i++) {
				pw.printf(taxaName[phenoIndexArray[i]]+"\t");
			}
			pw.printf(taxaName[phenoIndexArray[effectTaxaNum-1]]+"\n");
			for (int i = 0; i < gif.markerNum; i++) {
				pw.printf(gif.markers[i].markerName.replace(".", "_")+"\t"); //SAS doesn't accept header name with "."
				for (int j = 0; j < effectTaxaNum - 1; j++) {
					if (gif.genoValue[i][genoIndexArray[j]] == -1) {
						pw.printf("1"+"\t");   //for simple linear regression, it is "."; for GLMselect, it is "1".  Missing Data ("NN");
					}
					else if (gif.genoValue[i][genoIndexArray[j]] == -2) {
						pw.printf("."+"\t");
					}
					else {
						pw.printf(gif.genoValue[i][genoIndexArray[j]]+"\t");
					}
				}
				if (gif.genoValue[i][effectTaxaNum - 1] == -1) {
					pw.printf(1+"\n");
				}
				else {
					pw.printf(gif.genoValue[i][effectTaxaNum - 1]+"\n");
				}
			}
			for (int i = 0; i < phenoNum; i++) {
				int j = i + 1;
				pw.printf(j+"\t");
				for (j = 0; j < effectTaxaNum - 1; j++) {
					pw.printf(RWCP[i+hitOfCutoff].normalCountInTaxa[phenoIndexArray[j]]+"\t");
				}
				pw.printf(RWCP[i+hitOfCutoff].normalCountInTaxa[phenoIndexArray[effectTaxaNum - 1]]+"\n");
			}


/*
			//transposition 2
			pw.printf("Taxa"+"\t");
			for (int i = 0; i < gif.markerNum; i++) {
				pw.printf(gif.markers[i].markerName+"\t");
			}
			for (int i = 0; i < phenoNum-1; i++) {
				int j = i + 1;
				pw.print(j); pw.print("\t");
			}
			pw.println(phenoNum);
			for (int i = 0; i < phenoIndexArray.length; i++) {
				pw.print(taxaName[phenoIndexArray[i]]+"\t");
				for (int j = 0; j < gif.markerNum; j++) {
					if (gif.genoValue[j][genoIndexArray[i]] == -1) {
						pw.printf(1+"\t");
					}
					else {
						pw.printf(gif.genoValue[j][genoIndexArray[i]]+"\t");
					}
				}
				for (int j = 0; j < phenoNum-1; j++) {
					pw.print(RWCP[j+hitOfCutoff].freRelaCount[phenoIndexArray[i]]);
					pw.print("\t");
				}
				pw.println(RWCP[totalReadNum-1].freRelaCount[phenoIndexArray[i]]);
			}
 *
 */
			pw.flush();
			pw.close();
		}
		catch (Exception e) {
			System.out.println (genoPhenoFile);
		}


	}

	public void writeTable4QTL_Tassel (String forTasselDir, int hitOfCutoff) {
		GeneticInfoArray gif = new GeneticInfoArray(sourceGenoFile644);
		gif.setGenoValue4Tassel();
		new File (forTasselDir).mkdir();
		String genoFile = forTasselDir + "genofile.txt";
		String phenoFile = forTasselDir + "phenofile.txt";
		try {
			PrintWriter pwG = new PrintWriter(genoFile);
			PrintWriter pwP = new PrintWriter(phenoFile);
			int effectTaxaNum = 0;
			int phenoNum = totalReadNum - hitOfCutoff;
			ArrayList<Integer> phenoIndex = new ArrayList<Integer>();
			ArrayList<Integer> genoIndex = new ArrayList<Integer>();

			for (int i = 0; i < taxaName.length; i++ ) {
				for (int j = 0; j < gif.taxaNum; j++) {
					if (taxaName[i].equals(gif.taxaName[j])) {
						phenoIndex.add(i);
						genoIndex.add(j);
						effectTaxaNum++;
						break;
					}
				}
			}
			Integer[] pIndexArray = phenoIndex.toArray(new Integer[phenoIndex.size()]);
			Integer[] genoIndexArray = genoIndex.toArray(new Integer[genoIndex.size()]);
			pwG.println (genoIndexArray.length+"\t"+gif.markerNum+"\t"+1);
			for (int i = 0; i < gif.markerNum-1; i++) {
				pwG.printf(gif.markers[i].markerName+"\t");
			}
			pwG.println (gif.markers[gif.markerNum-1].markerName);
			for (int i = 0; i < genoIndexArray.length; i++) {
				pwG.print (gif.taxaName[genoIndexArray[i]]+"\t");
				for (int j = 0; j < gif.markerNum-1; j++) {
					pwG.printf(gif.genoValue4Tassel[j][genoIndexArray[i]]+"\t");
				}
				pwG.println(gif.genoValue4Tassel[gif.markerNum-1][genoIndexArray[i]]);
			}
			pwG.flush();
			pwG.close();
			pwP.println(pIndexArray.length+"\t"+phenoNum+"\t"+1);
			for (int i = 0; i < phenoNum-1; i++) {
				int j = i + 1;
				pwP.print(j); pwP.print("\t");
			}
			pwP.println(phenoNum);
			for (int i = 0; i < pIndexArray.length; i++) {
				pwP.print(taxaName[pIndexArray[i]]+"\t");
				for (int j = 0; j < phenoNum-1; j++) {
					pwP.print(RWCP[j+hitOfCutoff].normalCountInTaxa[pIndexArray[i]]);
					pwP.print("\t");
				}
				pwP.println(RWCP[totalReadNum-1].normalCountInTaxa[pIndexArray[i]]);
			}
			pwP.flush();
			pwG.close();
		}
		catch (Exception e) {
			System.out.println (genoFile+" or "+phenoFile);
		}
	}

	
	public void writeFreRelaCount (String freDesDir, int hitOfCutoff, boolean ifbinary) {
		new File(freDesDir).mkdir();
		if (!ifbinary) {
			String distributionFile = freDesDir+"freRelaCount.txt";
			try {
				PrintWriter pw = new PrintWriter(distributionFile);
				pw.println(totalReadNum-hitOfCutoff+"\t"+taxaName.length);
				for (int i = 0; i < taxaName.length - 1; i++) {
					pw.printf (taxaName[i]+"\t");
				}
				pw.println(taxaName[taxaName.length-1]);
				for (int i = hitOfCutoff; i < totalReadNum; i++) {
					for (int j = 0; j < taxaName.length - 1; j++) {
						pw.printf (RWCP[i].normalCountInTaxa[j]+"\t");
					}
					pw.println(RWCP[i].normalCountInTaxa[taxaName.length-1]);
				}
				pw.flush();
				pw.close();
			}
			catch (Exception e) {
				System.out.println (distributionFile);
			}
		}
		else {
			String distributionFile = freDesDir+"freRelaCount_bina.txt";
			try {
				DataOutputStream dos = new DataOutputStream(new BufferedOutputStream (new FileOutputStream(distributionFile),65536));
				dos.writeInt(totalReadNum-hitOfCutoff);
				dos.writeInt(taxaName.length);
				for (int i = 0; i < taxaName.length; i++) {
					dos.writeUTF(taxaName[i]);
				}
				for (int i = hitOfCutoff; i < totalReadNum; i++) {
					for (int j = 0; j < taxaName.length; j++) {
						dos.writeFloat(RWCP[i].normalCountInTaxa[j]);
					}
				}
				dos.flush();
				dos.close();
			}
			catch (Exception e) {
				System.out.println (distributionFile);
			}
		}
		//writeTable4QTL_Tassel(forTasselDir, hitOfCutoff); //optional choice
		writeTable4QTL_SAS (forSASDir, hitOfCutoff);
	}
	public void writeTotalCountDis () {
		String distributionFile = desDir+"totalCountDis.txt";
		int span = 2;
		int n = RWCP[totalReadNum-1].totalCount/span+1;
		int[] disA = new int[n];
		for (int i = 0; i < totalReadNum; i++) {
			n = RWCP[i].totalCount / span;
			disA[n] = disA[n]+1;
		}
		try {
			PrintWriter pw = new PrintWriter(distributionFile);
			for (int i = 0; i < disA.length; i++) {
				pw.println(disA[i]);
			}
			pw.flush();
			System.out.println("Done!");
		}
		catch (Exception e) {
			System.out.println (distributionFile);
		}
	}

	public void writeQuery (String desDir, int repeatCutoff) {
		int hit;
		File fn = new File(desDir);
		fn.mkdir();
		String desFile = desDir.replaceFirst("((\\w+)/$)", "$1$2\\.txt");
		ReadWCountPosition key = new ReadWCountPosition();
		key.setTotalCount(repeatCutoff);
		hit = Arrays.binarySearch(RWCP, key);
		if (hit < 0) {
			hit = -hit - 1;
		}
		else {
			while (RWCP[hit].totalCount == repeatCutoff) {
				hit--;
			}
			hit++;
		}
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(desFile), 65536);
			int j = 0;
			for (int i = hit; i < totalReadNum; i++) {
				j++;
				bw.write(RWCP[i].setTitle(j));
				bw.newLine();
				bw.write(BaseEncoder.getSequenceFromLong(RWCP[i].haplotype));
				bw.newLine();
			}
			int i = totalReadNum-hit;
			System.out.println ("There are "+i+" repeats in total.");
			bw.flush();
			bw.close();
		}
		catch (Exception e) {
			System.out.println (desFile);
		}
		writeFreRelaCount(freDesDir, hit, false);
	}
	
	public void setInfoArray (int size) {
		RWCP = new ReadWCountPosition[size];
	}

	public void getCount (String sourceFile1, String sourceFile2) {
		try {
			DataInputStream dis  = new DataInputStream (new BufferedInputStream (new FileInputStream(sourceFile1), 65536));
			int taxaNum = dis.readInt();
			totalReadNum = dis.readInt();
			setInfoArray (totalReadNum);
			taxaName = new String[taxaNum];
			for (int i = 0; i < taxaNum; i++) {
				taxaName[i] = dis.readUTF();
			}
			for (int i = 0; i < totalReadNum; i++) {
				long haplotype[] = new long[2];
				int totalCount = 0;
				byte[] countInTaxa = new byte[taxaNum];
				haplotype[0] = dis.readLong();
				haplotype[1] = dis.readLong();
				for (int j = 0; j < taxaNum; j++) {
					countInTaxa[j] = dis.readByte();
					totalCount = totalCount + countInTaxa[j];
				}
				RWCP[i] = new ReadWCountPosition();
				RWCP[i].setLongSeq(haplotype);
				RWCP[i].setTotalCount(totalCount);
				RWCP[i].setCountInTaxa(countInTaxa);
			}
			dis.close();
			dis  = new DataInputStream (new BufferedInputStream (new FileInputStream(sourceFile2), 65536));
			dis.readInt(); dis.readInt();
			for (int i = 0; i < taxaNum; i++) {
				dis.readUTF();
			}
			for (int i = 0; i < taxaNum; i++) {
				dis.readInt();
			}
			for (int i = 0; i < totalReadNum; i++) {
				float normalTotalCount = 0;
				float[] normalCountInTaxa = new float[taxaNum];
				dis.readLong(); dis.readLong();
				for (int j = 0; j < taxaNum; j++) {
					normalCountInTaxa[j] = dis.readFloat();
					normalTotalCount = normalTotalCount + normalCountInTaxa[j];
				}

				RWCP[i].setNormalTotalCount(normalTotalCount);
				RWCP[i].setNormalCountInTaxa(normalCountInTaxa);
			}
			dis.close();
		}
		catch (Exception e) {

		}
	}


	public void getPosition (String sourceFile) {
		try {
			DataInputStream dis = new DataInputStream(new BufferedInputStream (new FileInputStream(sourceFile),65536));
			for(int i=0; i < totalReadNum; i++) {
				long[] haplotype = new long[2];
                haplotype[0] = dis.readLong();
                haplotype[1] = dis.readLong();
                byte chrB=dis.readByte();
                byte strand=dis.readByte();
                int positionMin=dis.readInt();
                int positionMax=dis.readInt();
                short nextCutDistance=dis.readShort();
                byte divergence=dis.readByte();
                float dcoP=dis.readFloat();
                float mapP=dis.readFloat();
                byte multimaps=dis.readByte();
				RWCP[i].setPosition(chrB, strand, positionMin, positionMax, nextCutDistance, divergence, totalReadNum, totalReadNum, multimaps);
            }
			dis.close();
		}
		catch (Exception e) {
			System.out.println("Reading or writing error: " + sourceFile);
		}
	}

	public static void main (String args[]) {
		setQuery sq = new setQuery();
	}

}
