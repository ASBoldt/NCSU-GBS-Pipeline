/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.gbs.pav;

import cern.jet.random.Binomial;
import edu.cornell.lassp.houle.RngPack.RandomJava;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.TreeSet;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByte;
import net.maizegenetics.gbs.util.OpenBitSet;
import net.maizegenetics.gbs.util.SBitAlignment;
import net.maizegenetics.genome.BaseEncoder;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.ids.IdGroup;

/**
 * FamilyCode holds 64 families at most
 * Tag presence in one family is 0.01, p-value in one family is 0.05, p-value of joint linkage is 0.05, which can be filtered later
 * @author Fei Lu
 */
public class TagAgainstAnchorNAM {
	int ancChrCnt=10;
	int threadNumPerCore = 2;
	int threadNum;
	TagsByTaxaByte[] theTBTs;
	SBitAlignment[] theAnchorChrInBits;
	int[][] family2TBTIndex;
	int[][] tbt2anchorRedirect;
	Task[] jobs;
	int familyGroupNum = 2;

	public TagAgainstAnchorNAM (String tagsByTaxaFileS, String anchorFileS, String pedigreeFileS, String outfileS) {
		System.out.println("Program started at " + this.getCurrentTime());
		System.out.println("maxMemory: " + Runtime.getRuntime().maxMemory()/1024/1024 + "Mb; freeMemory: " + Runtime.getRuntime().freeMemory()/1024/1024 + "Mb");

		this.loadAnchorMaps(anchorFileS);

		this.calculateThreadNum();

		this.cutTBT(tagsByTaxaFileS);

		this.redirect(pedigreeFileS, true);

        double tMinP = 0.01; //Minimum tag presence in one family
        double pOne = 0.05; //p-value in one family to test segregation
        double pJoint = 0.05; //p-value of joint linkage in multiple family
		this.MTMapping(true, tMinP, pOne, pJoint);

		this.writeMappingResult(outfileS);
	}

	public static void main(String[] args) {
		String tagsByTaxaFileS = args[0];
		String anchorFileS = args[1];
		String pedigreeFileS = args[2];
		String outfileS = args[3];
		new TagAgainstAnchorNAM (tagsByTaxaFileS, anchorFileS, pedigreeFileS, outfileS);
	}

	private void writeMappingResult (String outfileS) {
		try {
			BufferedWriter bw = new BufferedWriter (new FileWriter (outfileS), 65536);
			String header = "TestTag\tTagTaxaCount\tFamilyGroupNum\t";
			for (int i = 0; i < familyGroupNum; i++) {
				String sufix = "-G"+String.valueOf(i+1);
				header = header + "LDChr" + sufix + "\t";
				header = header + "LDSite" + sufix + "\t";
				header = header + "LDPos" + sufix + "\t";
				header = header + "BinomP" + sufix + "\t";
				header = header + "ChrSig" + sufix + "\t";
				header = header + "Likelyhood" + sufix + "\t";
				header = header + "FamilyNum" + sufix + "\t";
				header = header + "FamilyCode" + sufix + "\t";
			}
			bw.write(header);
			bw.newLine();
			for (int i = 0; i < jobs.length; i++) {
				String[] result = jobs[i].getResult();
				for (int j = 0; j < result.length; j++) {
					bw.write(result[j]);
				}

				System.out.println("Mapping result of slice "+ (i + 1) + " with " + result.length + " tags is written");
			}
			bw.flush();
			bw.close();
			System.out.println("Mapping result written at " + this.getCurrentTime());
		}
		catch (Exception e) {
			System.out.println(e.toString());
		}
	}

	public void MTMapping(boolean ifMT, double tMinP, double pOne, double pJoint) {
		System.out.println ("Mapping started at " + this.getCurrentTime());
		long startTime = System.currentTimeMillis();
		jobs = new Task[threadNum];
		Thread[] mts = new Thread[threadNum];
		if (ifMT) {
			for (int i = 0; i < threadNum; i++) {
				jobs[i] = new Task (theTBTs[i], i, tMinP, pOne, pJoint);
				mts[i] = new Thread(jobs[i]);
				mts[i].start();
			}
			for (int i = 0; i < threadNum; i++) {
				try {
					mts[i].join();
				}
				catch (Exception e) {
					System.out.println(e.toString());
				}
			}
		}
		else {
			for (int i = 0; i < threadNum; i++) {
				jobs[i] = new Task (theTBTs[i], i, tMinP, pOne, pJoint);
				jobs[i].run();
			}
		}
		
		long endTime = System.currentTimeMillis();
		double secs = (double)(endTime - startTime)*0.001;
		System.out.println("Mapping finished in " + secs + " seconds!");
		System.out.println("Mapping finished at " + this.getCurrentTime());
	}

	private void redirect (String pedigreeFileS, boolean ifHMM) {
		Pedigree ped = new Pedigree (pedigreeFileS);
		String[] NAMfamily = ped.getFamilyStartWith("NAM_", 100);
		String[][] samples = ped.getSampleOfFamilies(NAMfamily, 100);
		if (ifHMM) {
			for (int i = 0; i < samples.length; i++) {
				for (int j = 0; j < samples[i].length; j++) {
					int index = samples[i][j].indexOf(":");
					samples[i][j] = samples[i][j].substring(0, index);
				}
			}
		}
		for (int i = 0; i < samples.length; i++) {
			Arrays.sort(samples[i]);
		}
		tbt2anchorRedirect = new int[NAMfamily.length][theTBTs[0].getTaxaNames().length];
		IdGroup anchorIdGroup=theAnchorChrInBits[0].getIdGroup();
		String[] taxaNameAnchor = new String[anchorIdGroup.getIdCount()];
		int[] taxaNameIndexAnchor = new int[anchorIdGroup.getIdCount()];
		int[] familyMark = new int[taxaNameAnchor.length];
		if (ifHMM) {
			for (int i = 0; i < taxaNameAnchor.length; i++) {
				taxaNameAnchor[i] = anchorIdGroup.getIdentifier(i).getFullName();
				String[] temp = taxaNameAnchor[i].split(":");
				taxaNameAnchor[i] = temp[0];
				taxaNameIndexAnchor[i] = i;
			}
		}
		else {
			for (int i = 0; i < taxaNameAnchor.length; i++) {
				taxaNameAnchor[i] = anchorIdGroup.getIdentifier(i).getFullName();
				String[] temp = taxaNameAnchor[i].split(":");
				taxaNameAnchor[i] = taxaNameAnchor[i].replaceFirst(":"+temp[temp.length-1], ""); //need to change after the new build is done
				taxaNameIndexAnchor[i] = i;
			}
		}
		

		for (int i = 0; i < taxaNameAnchor.length-1; i++) {
			for (int j = i + 1; j < taxaNameAnchor.length; j++) {
				if (taxaNameAnchor[i].compareTo(taxaNameAnchor[j])>0) {
					String temp;
					temp = taxaNameAnchor[i]; taxaNameAnchor[i] = taxaNameAnchor[j]; taxaNameAnchor[j] = temp;
					int tem;
					tem = taxaNameIndexAnchor[i]; taxaNameIndexAnchor[i] = taxaNameIndexAnchor[j]; taxaNameIndexAnchor[j] = tem;
				}
			}
		}

		for (int i = 0; i < taxaNameAnchor.length; i++) {
			familyMark[i] = -1;
			int hit = -1;
			for (int j = 0; j < samples.length; j++) {
				hit = Arrays.binarySearch(samples[j], taxaNameAnchor[i]);
				if (hit >= 0 ) {
					familyMark[i] = j;
					break;
				}
			}
			//System.out.println(taxaNameAnchor[i]+"\t"+hit);
		}

		for (int i = 0; i < theTBTs[0].getTaxaNames().length; i++) {
			String[] temp = theTBTs[0].getTaxaName(i).split(":");
			String query;
			if (ifHMM) {
				query = temp[0];
			}
			else {
				query = theTBTs[0].getTaxaName(i).replaceFirst(":"+temp[temp.length-1], ""); //need to change after the new build is done
			}
			int hit = Arrays.binarySearch(taxaNameAnchor, query);
			if (hit >= 0 && familyMark[hit] != -1) {
				tbt2anchorRedirect[familyMark[hit]][i] = taxaNameIndexAnchor[hit];
				for (int j = 0; j < tbt2anchorRedirect.length; j++) {
					if (j == familyMark[hit]) continue;
					tbt2anchorRedirect[j][i] = -1;
				}
			}
			else {
				for (int j = 0; j < tbt2anchorRedirect.length; j++) {
					tbt2anchorRedirect[j][i] = -1;
				}
			}
		}

		family2TBTIndex = new int[tbt2anchorRedirect.length][];
		for (int i = 0; i < tbt2anchorRedirect.length; i++) {
			ArrayList<Integer> temp = new ArrayList();
			for (int j = 0; j < tbt2anchorRedirect[i].length; j++) {
				if (tbt2anchorRedirect[i][j]>-1) temp.add(j);
			}
			family2TBTIndex[i] = new int[temp.size()];
			for (int j = 0; j < family2TBTIndex[i].length; j++) {
				family2TBTIndex[i][j] = temp.get(j);
			}
		}
	}

	private void cutTBT (String tagsByTaxaFileS) {
		theTBTs = new PAVUtils().sliceTBT(tagsByTaxaFileS, threadNum);
		System.out.println("TBT pieces cut at " + this.getCurrentTime());
	}

	private void calculateThreadNum () {
		int numOfProcessors = Runtime.getRuntime().availableProcessors();
		threadNum = numOfProcessors * threadNumPerCore;
		System.out.println("This node has " + numOfProcessors + " processors");
		System.out.println("TBT will be mapped by " + threadNum + " tasks");
		System.out.println("Each core runs " + threadNumPerCore + " tasks, or " + threadNumPerCore + " threads");
	}

	private void loadAnchorMaps(String anchorFileName) {
		System.out.println("Alignment Loaded:"+anchorFileName);
		if(anchorFileName.contains("+")||anchorFileName.contains("=")) {
			theAnchorChrInBits=new SBitAlignment[ancChrCnt];
			LoadChr[] lc = new LoadChr[ancChrCnt];
			Thread[] lt = new Thread[ancChrCnt];
            for (int i = 0; i < ancChrCnt; i++) {
                System.out.println("anchorFileName "+anchorFileName);
                String file=anchorFileName.replace("+",""+(i+1));
                file=file.replace("=",""+(i+1));
                System.out.println("Reading:"+file);
				lc[i] = new LoadChr(file,""+(i+1));
				lt[i] = new Thread(lc[i]);
                lt[i].start();
            }
			for (int i = 0; i < lt.length; i++) {
				try {
					lt[i].join();
				}
				catch (Exception e) {
					System.out.println(e.toString());
					System.exit(1);
				}
			}
			for (int i = 0; i < lc.length; i++) {
				theAnchorChrInBits[i] = lc[i].getBitAlignment();
				System.out.printf("Chr %s Sites %d Taxa %d %n", theAnchorChrInBits[i].getLocus(0),
                        theAnchorChrInBits[i].getSiteCount(),theAnchorChrInBits[i].getSequenceCount());
			}
			System.out.println("maxMemory: " + Runtime.getRuntime().maxMemory()/1024/1024 + "Mb; freeMemory: " + Runtime.getRuntime().freeMemory()/1024/1024 + "Mb");
        } else {
            System.out.println ("Can't find the anchor map files");
			System.exit(0);
        }
		System.out.println("Anchor maps load at " + this.getCurrentTime());
    }

	class Task implements Runnable {
		TagsByTaxaByte tbt;
		int taskID;
		String[] result;
        double tMinP;
		double pOne;
        double pJoint;

		Task (TagsByTaxaByte tbt, int index, double tMinP, double pOne, double pJoint) {
			this.tbt = tbt;
			this.taskID = index + 1;
            this.tMinP = tMinP;
            this.pOne = pOne;
            this.pJoint = pJoint;
		}

		public String[] getResult () {
			return result;
		}

		private int[] getMinPresentFamily (int tagIndex, double minPresentInFamily) {
			ArrayList<Integer> overMinPresentFamilyList = new ArrayList();
			for (int i = 0; i < family2TBTIndex.length; i++) {
				int cnt = 0;
				for (int j = 0; j < family2TBTIndex[i].length; j++) {
					if (tbt.getReadCountForTagTaxon(tagIndex, family2TBTIndex[i][j]) != 0) cnt++;
				}
				double rate = (double)cnt/family2TBTIndex[i].length;
				if (rate > minPresentInFamily) overMinPresentFamilyList.add(i);
			}
			if (overMinPresentFamilyList.isEmpty()) return null;
			int[] overMinPresentFamily = new int[overMinPresentFamilyList.size()];
			for (int i = 0; i < overMinPresentFamily.length; i++) {
				overMinPresentFamily[i] = overMinPresentFamilyList.get(i);
			}
			return overMinPresentFamily;
		}

		public ArrayList getSegregationFamilyGroups (int tagIndex, int[] minPresentFamilies, int step, double pThresh) {
			ArrayList<Integer> segregationFamilyList = new ArrayList();
			ArrayList<Integer> segregationChromList = new ArrayList();
			for (int i = 0; i < minPresentFamilies.length; i++) {
				ScanChromosome[] scanOnChr = new ScanChromosome[theAnchorChrInBits.length];
				int[] testFamilies = new int[1];
				testFamilies[0] = minPresentFamilies[i];
				double[][] theResults=new double[theAnchorChrInBits.length][];
				for (int j = 0; j < theAnchorChrInBits.length; j++) {
					scanOnChr[j] = new ScanChromosome(tbt, tagIndex, theAnchorChrInBits[j], testFamilies, j, theResults, pThresh, step);
					scanOnChr[j].scan();
				}
				int bestChr = -1;
				double bestP = pThresh;
				for (int j = 0; j < theAnchorChrInBits.length; j++) {
					if (theResults[j] == null) continue;
					if (theResults[j][3] < bestP) {
						bestP = theResults[j][3]; //or use the countSig to decide if the tag belong to a chromosome
						bestChr = (int)theResults[j][0];
					}
				}
				if (bestChr != -1) {
					segregationFamilyList.add(minPresentFamilies[i]);
					segregationChromList.add(bestChr);
				}
			}
			if (segregationFamilyList.isEmpty()) return null;
			TreeSet<Integer> bestChrSet = new TreeSet();
			for (int i = 0; i < segregationChromList.size(); i++) {
				bestChrSet.add(segregationChromList.get(i));
			}
			Integer[] bestChrArray = bestChrSet.toArray(new Integer[bestChrSet.size()]);
			int[] bestChrCount = new int[bestChrArray.length];
			for (int i = 0; i < bestChrArray.length; i++) {
				bestChrCount[i] = 0;
				for (int j = 0; j < segregationChromList.size(); j++) {
					if (bestChrArray[i] == segregationChromList.get(j)) {
						bestChrCount[i]++;
					}
				}
			}
			for (int i = 0; i < bestChrCount.length - 1; i++) {
				for (int j = i + 1; j < bestChrCount.length; j++) {
					if (bestChrCount[i] < bestChrCount[j]) {
						int temp;
						temp = bestChrCount[i]; bestChrCount[i] = bestChrCount[j]; bestChrCount[j] = temp;
						temp = bestChrArray[i]; bestChrArray[i] = bestChrArray[j]; bestChrArray[j] = temp;
					}
				}
			}
			int[][] segregationFamilyGroups = new int[bestChrArray.length][];
			for (int i = 0; i < segregationFamilyGroups.length; i++) {
				segregationFamilyGroups[i] = new int[bestChrCount[i]];
				int cnt = 0;
				for (int j = 0; j < segregationChromList.size(); j++) {
					if (bestChrArray[i] == segregationChromList.get(j)) {
						segregationFamilyGroups[i][cnt] = segregationFamilyList.get(j);
						cnt++;
					}
				}
			}
			int[] bestChrs = new int[bestChrArray.length];
			for (int i = 0; i < bestChrArray.length; i++) bestChrs[i] = bestChrArray[i];
			ArrayList<Object> segregationList = new ArrayList();
			segregationList.add(segregationFamilyGroups);
			segregationList.add(bestChrs);
			return segregationList;
		}

		public String getOutputRecord (int tagIndex, double[][] theResults) {
			int GNum = 0;
			for (int i = 0; i < theResults.length; i++) {
				if (theResults[i] != null) GNum++;
			}
			StringBuilder sb = new StringBuilder();
			sb.append(BaseEncoder.getSequenceFromLong(tbt.getTag(tagIndex))).append("\t");
			sb.append(String.valueOf(tbt.getNumberOfTaxaWithTag(tagIndex))).append("\t");
			sb.append(String.valueOf(GNum)).append("\t");
	//chr, bestSite, refAlignment.getPositionInLocus(bestSite), bestP, countSig, (double)-1, bestNewTestFamilies.length, (double)this.getFamilyCode(bestNewTestFamilies)}
			for (int i = 0; i < theResults.length; i++) {
				if (theResults[i] == null) continue;
				sb.append(String.valueOf((int)theResults[i][0])).append("\t");
				sb.append(String.valueOf((int)theResults[i][1])).append("\t");
				sb.append(String.valueOf((int)theResults[i][2])).append("\t");
				sb.append(String.valueOf(theResults[i][3])).append("\t");
				sb.append(String.valueOf((int)theResults[i][4])).append("\t");
				sb.append(String.valueOf(theResults[i][5])).append("\t");
				sb.append(String.valueOf((int)theResults[i][6])).append("\t");
				sb.append(String.valueOf((long)theResults[i][7])).append("\t");
			}
			sb.append("\n");
			return sb.toString();
		}

		public void run() {
            long starttime = System.currentTimeMillis();
			ArrayList<String> resultList = new ArrayList();
            int tagInFamilycnt = 0;
			for (int i = 0; i < tbt.getTagCount(); i++) {
				int[] minPresentFamilies = getMinPresentFamily(i, tMinP);
                if(minPresentFamilies == null) continue;
                tagInFamilycnt++;
				ArrayList segregationList = this.getSegregationFamilyGroups(i, minPresentFamilies, 10, pOne);
				if (segregationList == null) continue;
				int[][]	segregationFamilyGroups = (int[][])segregationList.get(0);
				int[] bestChrs = (int[])segregationList.get(1);
				int minGnum = segregationFamilyGroups.length;
				if (familyGroupNum < minGnum) minGnum = familyGroupNum;
				double[][] theResults = new double[minGnum][];
				boolean nullFlag = true;
				for (int j = 0; j < minGnum; j++) { //it should be bestChrs.length, only caculate 2 loci for now
					ScanChromosome scanOnBestChrom = new ScanChromosome(tbt, i, theAnchorChrInBits[bestChrs[j]-1], segregationFamilyGroups[j], j, theResults, pJoint, 1);
					scanOnBestChrom.scan();
					if (theResults[j] == null) continue;
					nullFlag = false;
				}
				if (nullFlag) continue;
				resultList.add(this.getOutputRecord(i, theResults));
				
            }
            System.out.println("Number of tested tags in families is " + String.valueOf(tagInFamilycnt));
			result = resultList.toArray(new String[resultList.size()]);
			double secs = (double)(System.currentTimeMillis() - starttime) * 0.001;
			System.out.println("Task " + taskID + " with " + tbt.getTagCount() + " tags is finished in " + secs + " seconds");
		}

	}

	private class ScanChromosome {
		TagsByTaxaByte tbt;
		int tagIndex;
        SBitAlignment refAlignment;
		int[] testFamilies;
		int resultIndex;
		double[][] resultReport;
		double sigThreshold;
		int step = 1;
        int minSigPos=-1, maxSigPos=-1;

		public ScanChromosome (TagsByTaxaByte tbt, int tagIndex, SBitAlignment refAlignment, int[] testFamilies, int resultIndex, double[][] resultReport, double sigThreshold, int step) {
			this.tbt = tbt;
			this.tagIndex = tagIndex;
			this.refAlignment = refAlignment;
			this.testFamilies = testFamilies;
			this.resultIndex = resultIndex;
			this.resultReport = resultReport;
			this.sigThreshold = sigThreshold;
			this.step = step;
		}

		public void scan () {
			int bestSite = -1, countSig = 0;
			double bestP = 1;
			int[] bestNewTestFamilies = null;
			for (int i = 0; i < refAlignment.getSiteCount(); i+=step) {
				OpenBitSet dA1 = refAlignment.getSiteBitsNoClone(i, 0);
                OpenBitSet dA2 = refAlignment.getSiteBitsNoClone(i, 1);
				int[] newTestFamilies = this.getNewTestFamiliesBySegregationOfAnchor(dA1, dA2);
				if (newTestFamilies == null) continue;
				long[] tagDistBits = this.getTagDistInBits(newTestFamilies);
				OpenBitSet dT = new OpenBitSet(tagDistBits, tagDistBits.length);
				long[] testFamiliesBits = this.getFamiliesMarkInBits(newTestFamilies);
				OpenBitSet dF = new OpenBitSet(testFamiliesBits, testFamiliesBits.length);
				double currP = this.testSiteBio(dA1, dA2, dT, dF);
				if (currP < bestP) {
					bestP = currP;
					bestSite = i;
					bestNewTestFamilies = newTestFamilies.clone();
				}
				if (currP < this.sigThreshold) {
					countSig++;
                    if(minSigPos==-1) minSigPos=refAlignment.getPositionInLocus(i);
                    maxSigPos=refAlignment.getPositionInLocus(i);
				}
			}
			if (bestP > this.sigThreshold) {
				resultReport[resultIndex] = null;
			}
			else {
				int chr=Integer.parseInt(refAlignment.getLocus(bestSite).getChromosomeName());
				double[] result={chr, bestSite, refAlignment.getPositionInLocus(bestSite), bestP, countSig, (double)-1, bestNewTestFamilies.length, (double)this.getFamilyCode(bestNewTestFamilies)};
				resultReport[resultIndex]=result;
			}
		}

		public double testSiteBio (OpenBitSet dA1, OpenBitSet dA2, OpenBitSet dT, OpenBitSet dF) {
			double result = 1, minorP = 1;
			int cntIA1 = 0, cntIA2 = 0, cntIA1T = 0, cntIA2T = 0, cntIA1F = 0, cntIA2F = 0;
			cntIA1T = (int)OpenBitSet.intersectionCount(dT, dA1);
			cntIA2T = (int)OpenBitSet.intersectionCount(dT, dA2);
			cntIA1F = (int)OpenBitSet.intersectionCount(dF, dA1);
			cntIA2F = (int)OpenBitSet.intersectionCount(dF, dA2);
			int sumTag = cntIA1T + cntIA2T;
			if (sumTag==0) return result;
			int sumAnchor = cntIA1F + cntIA2F;
			if (cntIA1T < cntIA2T) {
				minorP = (double) cntIA1F/sumAnchor;
			}
			else {
				minorP = (double) cntIA2F/sumAnchor;
			}
			Binomial binomFunc=new Binomial(5, 0.5, new RandomJava());
			binomFunc.setNandP(sumTag, minorP);
			try {
				result = (cntIA1T<cntIA2T)?binomFunc.cdf(cntIA1T):binomFunc.cdf(cntIA2T);
			}
			catch (Exception e) {
				System.out.println(e.toString());
			}
			return result;
		}

		public double testSiteSimple (OpenBitSet dA1, OpenBitSet dA2, OpenBitSet dT, OpenBitSet dF) {
			double result=1;
			int cntIA1 = 0, cntIA2 = 0;
			cntIA1 = (int)OpenBitSet.intersectionCount(dT, dA1);
			cntIA2 = (int)OpenBitSet.intersectionCount(dT, dA2);

			int sum = cntIA1 + cntIA2;
			if (sum > 0) {
				if (cntIA1 < cntIA2) result = (double)cntIA1/sum;
				else result = (double)cntIA2/sum;
			}
/*
			if (result < 0.0001) {
				int cA1 = (int)dA1.lastCalculatedCardinality();
				int cA2 = (int)dA2.lastCalculatedCardinality();
				int cT = (int)dT.lastCalculatedCardinality();
				int cF = (int)dF.lastCalculatedCardinality();
				int cntIA1F = (int)OpenBitSet.intersectionCount(dF, dA1);
				int cntIA2F = (int)OpenBitSet.intersectionCount(dF, dA2);
				System.out.println(cntIA1 + "\t" + cntIA2 + "\t" + sum +"\t"+cA1+"\t"+cA2+"\t"+cT+"\t"+cF+"\t"+cntIA1F+"\t"+cntIA2F);
			}
 *
 */
			return result;
		}

		public int[] getNewTestFamiliesBySegregationOfAnchor (OpenBitSet dA1, OpenBitSet dA2) {
			ArrayList<Integer> familyList = new ArrayList();
			for (int i = 0; i < testFamilies.length; i++) {
				long[] testFamilyBits = this.getFamilyMarkInBits(testFamilies[i]);
				OpenBitSet familyMark = new OpenBitSet(testFamilyBits, testFamilyBits.length);
				if (this.checkAnchorSegregationInFamilies(familyMark, dA1, dA2, family2TBTIndex[testFamilies[i]].length)) familyList.add(testFamilies[i]);
			}
			if (familyList.isEmpty()) return null;
			int[] newTestFamilies = new int[familyList.size()];
			for (int i = 0; i < newTestFamilies.length; i++) {
				newTestFamilies[i] = familyList.get(i);
			}
			return newTestFamilies;
		}

		public long[] getFamilyMarkInBits (int familyIndex) {
			long[] mark = new long[(refAlignment.getSequenceCount()/64) + 1];
			for (int i = 0; i < family2TBTIndex[familyIndex].length; i++) {
				int index = tbt2anchorRedirect[familyIndex][family2TBTIndex[familyIndex][i]]/64;
				int offset = tbt2anchorRedirect[familyIndex][family2TBTIndex[familyIndex][i]]%64;
				mark[index]=mark[index]|(1L<<offset);
			}
			return mark;
		}

		public long[] getFamiliesMarkInBits (int[] testFamilies) {
			long[] mark = new long[(refAlignment.getSequenceCount()/64) + 1];
			for (int i = 0; i < testFamilies.length; i++) {
				for (int j = 0; j < family2TBTIndex[testFamilies[i]].length; j++) {
					int index = tbt2anchorRedirect[testFamilies[i]][family2TBTIndex[testFamilies[i]][j]]/64;
					int offset = tbt2anchorRedirect[testFamilies[i]][family2TBTIndex[testFamilies[i]][j]]%64;
					mark[index]=mark[index]|(1L<<offset);
				}
			}
			return mark;
		}

		public long[] getTagDistInBits (int[] newTestFamilies) {
			int lgPerSite = (refAlignment.getSequenceCount() / 64) + 1;
			long[] seq = new long[lgPerSite];
			for (int i = 0; i < newTestFamilies.length; i++) {
				for (int j = 0; j < family2TBTIndex[newTestFamilies[i]].length; j++) {
					int index = tbt2anchorRedirect[newTestFamilies[i]][family2TBTIndex[newTestFamilies[i]][j]]/64;
					int offset = tbt2anchorRedirect[newTestFamilies[i]][family2TBTIndex[newTestFamilies[i]][j]]%64;
					if (tbt.getReadCountForTagTaxon(tagIndex, family2TBTIndex[newTestFamilies[i]][j])>0) {  //reference alleles
						seq[index]=seq[index]|(1L<<offset);
					}
				}
			}
			return seq;
		}

		public boolean checkAnchorSegregationInFamilies (OpenBitSet familiesMark, OpenBitSet dA1, OpenBitSet dA2, int familySize) {
			int cntA1 = (int)OpenBitSet.intersectionCount(familiesMark, dA1);
			int cntA2 = (int)OpenBitSet.intersectionCount(familiesMark, dA2);
			if (cntA1 > cntA2) {
				int temp = cntA1; cntA1 = cntA2; cntA2 = temp;
			}
			int sum = cntA1 + cntA2;
			double rate = (double)cntA1 / sum;
			if (rate < 0.2) return false; //check segregation
			//if ((double)sum/familySize < 0.8) return false; //check missing data
			return true;
		}

		public long getFamilyCode (int[] newTestFamilies) {
			long code = 0;
			for (int i = 0; i < newTestFamilies.length; i++) {
				code = code|(1L<<(newTestFamilies[i]));
			}
			return code;
		}

		public int[] getFamiliesFromCode (long code) {
			ArrayList<Integer> list = new ArrayList();
			for (int i = 0; i < 64; i++) {
				long mark = 1L << i;
				if ((code & mark) != 0) list.add(i);
			}
			int[] families = new int[list.size()];
			for (int i = 0; i < families.length; i++) {
				families[i] = list.get(i);
			}
			return families;
		}
    }

	class LoadChr implements Runnable {
		String filename;
		String chrom;
		Alignment a;
		SBitAlignment chrInBits;

		LoadChr (String filename, String chrom) {
			 this.filename = filename;
			 this.chrom = chrom;
		}

		public void run() {
			System.out.println("Chromosome " + chrom + " start loading at " + getCurrentTime());
			a = ImportUtils.readFromHapmap(filename, chrom);
			chrInBits = new SBitAlignment(a);
			System.out.println("Chromosome " + chrInBits.getLocusName(0) + " is loaded " + getCurrentTime());
		}

		public SBitAlignment getBitAlignment () {
			return chrInBits;
		}
	}

	private String getCurrentTime() {
		SimpleDateFormat sd = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        Date date1 = new Date();
        String str1 = sd.format(date1);
        return str1;
	}
}
