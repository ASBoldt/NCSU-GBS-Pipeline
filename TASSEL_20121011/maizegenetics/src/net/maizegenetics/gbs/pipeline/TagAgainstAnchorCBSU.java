/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.gbs.pipeline;

import cern.jet.random.Binomial;
import edu.cornell.lassp.houle.RngPack.RandomJava;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import net.maizegenetics.gbs.pav.PAVUtils;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByte;
import net.maizegenetics.gbs.util.SBitAlignment;
import net.maizegenetics.gbs.util.OpenBitSet;
import net.maizegenetics.genome.BaseEncoder;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.ids.IdGroup;

/**
 * Read TBT slices from a directory, process one after another, each slice is cut into pieces, each piece has a thread, results is written to a dir one by one
 * Memory need: 10 chrs (7G) + TBT slice (depends)
 * @author Fei Lu
 */
public class TagAgainstAnchorCBSU {
	int ancChrCnt=10;
	int threadNumPerCore = 32; //key to performance, 32 or 64 tasks have the best performance, almost linear to the number of cores.
	int threadNum;
	File[] TBTFiles;
	SBitAlignment[] theAnchorChrInBits;
	TagsByTaxaByte[] subTBTs;
	int[][] tbt2anchRedirect;
	Task[] jobs;

	public TagAgainstAnchorCBSU (String tagsByTaxaDirS, String anchorFileS, String outfileDirS) {
		System.out.println("Program started at " + this.getCurrentTime());
		System.out.println("maxMemory: " + Runtime.getRuntime().maxMemory()/1024/1024 + "Mb; freeMemory: " + Runtime.getRuntime().freeMemory()/1024/1024 + "Mb");

		this.loadAnchorMaps(anchorFileS);
		
		this.calculateThreadNum();

		this.getTBTFileNameInDir(tagsByTaxaDirS);

		this.loopMapping(outfileDirS);

	}

	private void loopMapping (String outfileDirS) {
		Arrays.sort(TBTFiles);
		for (int i = 0; i < TBTFiles.length; i++) {
			subTBTs = null;
			System.gc();
			this.cutTBT(TBTFiles[i].getAbsolutePath());
			if (i == 0) this.redirect();
			this.MTmapping();
			String outfileS = (outfileDirS + "/" + TBTFiles[i].getName()).replaceAll("tbt.byte", "mapping.txt");
			this.writeMappingResult(outfileS);
		}
	}

	private void getTBTFileNameInDir (String tagsByTaxaDirS) {
		TBTFiles = new File (tagsByTaxaDirS).listFiles();
		for (int i = 0; i < TBTFiles.length; i++) {
			System.out.println(TBTFiles[i].getAbsolutePath() + " will be processed");
		}
	}

	private void writeMappingResult (String outfileS) {
		try {
			BufferedWriter bw = new BufferedWriter (new FileWriter (outfileS), 65536);
			String header="TestTag\tTestTagNum\tBlastChr\tBlastPos\trefDiv\tLDChr\tLDSite\tLDPos\tBinomP\tSigTests\tTagTaxaCnt\tChrSig\tLRatioB:2\tLRatioB:M\tSiteOnBestChrThanNextChr\tMinSigPos\tMaxSigPos";
			bw.write(header);
			bw.newLine();
			for (int i = 0; i < jobs.length; i++) {
				String[] result = jobs[i].getResult();
				for (int j = 0; j < result.length; j++) {
					bw.write(result[j]);
				}
			}
			bw.flush();
			bw.close();
			System.out.println("Mapping result written in " + outfileS + " at " + this.getCurrentTime());
		}
		catch (Exception e) {
			System.out.println(e.toString());
		}
	}

	private void MTmapping () {
		System.out.println ("Mapping started at " + this.getCurrentTime());
		long startTime = System.currentTimeMillis();
		jobs = new Task[threadNum];
		Thread[] mts = new Thread[threadNum];
		for (int i = 0; i < threadNum; i++) {
			jobs[i] = new Task (subTBTs[i], i);
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
		long endTime = System.currentTimeMillis();
		double secs = (double)(endTime - startTime)*0.001;
		System.out.println("Mapping finished in " + secs + " seconds!");
		System.out.println("Mapping finished at " + this.getCurrentTime());
	}

	class Task implements Runnable {
		TagsByTaxaByte tbt;
		int taskID;
		String[] result;
		Task (TagsByTaxaByte tbt, int index) {
			this.tbt = tbt;
			this.taskID = index + 1;
		}

		public String[] getResult () {
			return result;
		}

		public void run() {
            long starttime = System.currentTimeMillis();
			ArrayList<String> resultList = new ArrayList();
			for (int i = 0; i < tbt.getTagCount(); i++) {
                if(tbt.getNumberOfTaxaWithTag(i)<30) continue;
                double[] bestR={-1,-1,-1, 1, -1};
                int blastChr=-1, blastPos=-1, bestAlignment=-1, refDiv=-1;
                long[] testTag = tbt.getTag(i);
                long[] testTagDist;
                double[][] theResults=new double[theAnchorChrInBits.length][];
                ScanChromosomeMTR[] scanOnChr = new ScanChromosomeMTR[theAnchorChrInBits.length];
                for (int cn=0; cn<theAnchorChrInBits.length; cn++) {
                    testTagDist=getTagsInBits(tbt,i,tbt2anchRedirect[cn],theAnchorChrInBits[cn].getSequenceCount());
                    scanOnChr[cn] = new ScanChromosomeMTR(theAnchorChrInBits[cn],testTagDist, cn, theResults,0.000001,1, blastPos);
                    scanOnChr[cn].scan();
                }
                int countRealSig=0;
                double[] pRank=new double[theAnchorChrInBits.length];
                for (int cn = 0; cn < theAnchorChrInBits.length; cn++) {
                    double[] r=theResults[cn];
                    pRank[cn]=r[3];
                    if(r[3]<bestR[3]) {bestR=r.clone(); bestAlignment=cn;}
                    if(r[3]<0.000001) countRealSig++;
                }
                Arrays.sort(pRank);
                double[][] bestResWithNewThreshold=new double[1][];
                testTagDist=getTagsInBits(tbt,i,tbt2anchRedirect[bestAlignment],theAnchorChrInBits[bestAlignment].getSequenceCount());
                ScanChromosomeMTR bestChrNewThres = new ScanChromosomeMTR(theAnchorChrInBits[bestAlignment],testTagDist, 0, bestResWithNewThreshold,pRank[1],1, blastPos);
                bestChrNewThres.scan();
                int countOfSitesBetterThanNextBestChr=(int)bestResWithNewThreshold[0][4];
                String s=String.format("%s %d %d %d %d %d %d %d %g %d %d %d %g %g %d %d %d %n",BaseEncoder.getSequenceFromLong(testTag),i,blastChr, blastPos, refDiv,
                        (int)bestR[0],(int)bestR[1],(int)bestR[2],bestR[3], (int)bestR[4], tbt.getNumberOfTaxaWithTag(i),
                        countRealSig, Math.log10(pRank[1]/pRank[0]), Math.log10(pRank[theAnchorChrInBits.length/2]/pRank[0]), countOfSitesBetterThanNextBestChr,
                        bestChrNewThres.minSigPos, bestChrNewThres.maxSigPos);
                resultList.add(s);
            }
			result = resultList.toArray(new String[resultList.size()]);
			double secs = (double)(System.currentTimeMillis() - starttime) * 0.001;
			System.out.println("Task " + taskID + " with " + tbt.getTagCount() + " tags is finished in " + secs + " seconds");
		}

	}

	private class ScanChromosomeMTR {
        SBitAlignment refAlignment;
        OpenBitSet obsymj;
        Binomial binomFunc=new Binomial(5, 0.5, new RandomJava());
        int resultIndex;
        double[][] resultReport;
        double sigThreshold;
        int minSigPos=-1, maxSigPos=-1;
        int step=1;
        double bionomialThreshold=0.2;
        int blockPosition=-1000;  //position to block from testing.
        int blockWindow=100;

        public ScanChromosomeMTR(SBitAlignment refAlignment, long[] ymj, int resultIndex, double[][] resultReport, double sigThreshold, int step, int blockPosition) {
            this.refAlignment=refAlignment;
            obsymj=new OpenBitSet(ymj,ymj.length);
            this.resultIndex=resultIndex;
            this.resultReport=resultReport;
            this.sigThreshold=sigThreshold;
            this.blockPosition=blockPosition;
            this.step=step;
        }

        public void scan() {
            //in the future this may call other scanOnChr to split the effort up more
            long tests=0;
            int bestSite=-1, countSig=0;
            double bestP=2;
            for (int i = 0; i < refAlignment.getSiteCount(); i+=step) {
             //   if(i%10000==0) System.out.println("scanChromosome chr:"+refAlignment.getChr()+"+site:"+i+" with Chr:"+refAlignment.getChr()+" test:"+tests);
                if(Math.abs(refAlignment.getPositionInLocus(i)-blockPosition)<blockWindow) continue;
                OpenBitSet mj=refAlignment.getSiteBitsNoClone(i, 0);
                OpenBitSet mn=refAlignment.getSiteBitsNoClone(i, 1);
                if(mn.lastCalculatedCardinality()>4) {
                    double p=TagCallerAgainstAnchorMT.testSites(obsymj, mj, mn, binomFunc);
                    if(p<bestP) {bestP=p; bestSite=i;}
                    if(p<sigThreshold) {
                        countSig++;
                        if(minSigPos==-1) minSigPos=refAlignment.getPositionInLocus(i);
                        maxSigPos=refAlignment.getPositionInLocus(i);
                    }
                }
                tests++;
            }
            int chr=Integer.parseInt(refAlignment.getLocus(bestSite).getChromosomeName());
            double[] result={chr, bestSite, refAlignment.getPositionInLocus(bestSite),bestP, countSig};
          //  if(bestP<0.000001) System.out.println(Arrays.toString(result));
            resultReport[resultIndex]=result;
        }

    }

	private long[] getTagsInBits(TagsByTaxaByte aTBT, int tagIndex, int[] reDirect, int anchorTaxa) {
        int lgPerSite = (anchorTaxa / 64) + 1;
        long[] seq = new long[lgPerSite];
        for (int j = 0; j < aTBT.getTaxaCount(); j++) {
            if(reDirect[j]<0) continue;
            int index=reDirect[j]/64;
            int offset=reDirect[j]%64;
            if (aTBT.getReadCountForTagTaxon(tagIndex, j)>0) {  //reference alleles
                seq[index]=seq[index]|(1L<<offset);
            }
        }
        return seq;
	}

	private String getCurrentTime() {
		SimpleDateFormat sd = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        Date date1 = new Date();
        String str1 = sd.format(date1);
        return str1;
	}

	private void redirect () {
		tbt2anchRedirect = new int[ancChrCnt][subTBTs[0].getTaxaNames().length];
        for(int c=0; c<ancChrCnt; c++) {
            for (int t = 0; t < subTBTs[0].getTaxaNames().length; t++) {
                IdGroup anchorIdGroup=theAnchorChrInBits[c].getIdGroup();
                tbt2anchRedirect[c][t] = anchorIdGroup.whichIdNumber(subTBTs[0].getTaxaNames()[t]);
            }
        }
		System.out.println("Taxa name redirected at " + this.getCurrentTime());
	}

	private void cutTBT (String tagsByTaxaFileS) {
		subTBTs = new PAVUtils().sliceTBT(tagsByTaxaFileS, threadNum);
		System.out.println("TBT pieces of " + tagsByTaxaFileS + " are cut at " + this.getCurrentTime());
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
            Alignment[] a = ImportUtils.readFromHapmap(anchorFileName).getAlignments();
            theAnchorChrInBits=new SBitAlignment[a.length];
            for (int i = 0; i < a.length; i++) {
                theAnchorChrInBits[i] = new SBitAlignment(a[i]);
    //            System.out.printf("Chr %d Sites %d Taxa %d %n", theAnchorChrInBits[i].getChr(),
    //                    theAnchorChrInBits[i].getNumSites(),theAnchorChrInBits[i].getNumTaxa());
                System.out.printf("Chr %s Sites %d Taxa %d %n", theAnchorChrInBits[i].getLocus(0),
                        theAnchorChrInBits[i].getSiteCount(),theAnchorChrInBits[i].getSequenceCount());
            }
        }
		System.out.println("Anchor maps load at " + this.getCurrentTime());
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

	public static void main (String[] args) {
//		String tagsByTaxaDirS = "M:/slice/";
//		String anchorFileS = "M:/pav/imp" + "/Zea20120110_scv10mF8maf002_mgs_E1pLD5kpUn.imp95_1024.c+.hmp.txt";
//		String outfileDirS = "M:/Result/";
		String tagsByTaxaDirS = args[0];
		String anchorFileS = args[1] + "/Zea20120110_scv10mF8maf002_mgs_E1pLD5kpUn.imp95_1024.c+.hmp.txt";
		String outfileDirS = args[2];
		new TagAgainstAnchorCBSU (tagsByTaxaDirS, anchorFileS, outfileDirS);
	}
} 