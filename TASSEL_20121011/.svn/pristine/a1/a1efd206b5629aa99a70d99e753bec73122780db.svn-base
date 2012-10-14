/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.genome.GBS.switchgrass;

import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import net.maizegenetics.genome.BaseEncoder;

/**
 *
 * @author Fei Lu
 */
public class NetworkFilter {
    TagAlignment[] tps;

	public NetworkFilter () {
	}

    public NetworkFilter (ReadCounts rc, double errorRate) {
        this.getTagPair(rc, errorRate);
		this.reciprocate();
    }

	public NetworkFilter (ReadCounts rc, double errorRate, String outfileS) {
		this(rc, errorRate);
		this.writeTagPair(rc, outfileS);
	}

	private final void getAlignment (ReadCounts rc){// this is only for the paper to show how effective the filter is, filter doesn't work here
		TagpairFinder tf = new TagpairFinder (rc);
		ArrayList<Float> freList = new ArrayList();
        ArrayList<TagAlignment> tpList = new ArrayList();
		int total = 3000000;
		String outfileS = "M:/dis.txt";
		for (int i = 0; i < rc.getSize(); i++) {
			long[] queryLongSeq = new long[2];
			queryLongSeq[0] = rc.haplotype[0][i];
			queryLongSeq[1] = rc.haplotype[1][i];
			ArrayList<Integer> hitIndex = tf.findOneMismatch(queryLongSeq);
			if (hitIndex.isEmpty()) {
				continue;
			}
			Integer[] hitIndexArray = hitIndex.toArray(new Integer[hitIndex.size()]);
			for (int j = 0; j < hitIndexArray.length; j++) {
				int n1 = rc.hapcount[i];
				int n2 = rc.hapcount[hitIndexArray[j]];
				int sum = n1 + n2;
				if (sum < 30) continue;
				float f = (float)n1/(float)sum;
				freList.add(f);
			}
		}
		Float[] fre = freList.toArray(new Float[freList.size()]);
		System.out.println(fre.length);
		int[] count = new int[100];
		for (int i = 0; i < count.length; i++) count[i] = 0;
		for (int i = 0; i < fre.length; i++) {
			int index = (int)(fre[i]*100);
			count[index]++;
		}
		float[] disFraction = new float[count.length];
		int sum = 0;
		for (int i = 0; i < count.length; i++) {
			sum += count[i];
		}
		for (int i = 0; i < disFraction.length; i++) {
			disFraction[i] = (float) count[i] / (float) sum;
		}
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outfileS), 65536);
			for (int i = 0; i < disFraction.length; i++) {
				bw.write(String.valueOf(disFraction[i]));
				bw.newLine();
			}
			bw.flush();
			bw.close();
		}
		catch (Exception e) {

		}
		
	}

	public final void getTagPair (ReadCounts rc, double errorRate) {
        TagpairFinder tf = new TagpairFinder (rc);
        ArrayList<TagAlignment> tpList = new ArrayList();
		int cnt1 = 0, cnt2 = 0;
		for (int i = 0; i < rc.getSize(); i++) {
			long[] queryLongSeq = new long[2];
			queryLongSeq[0] = rc.haplotype[0][i];
			queryLongSeq[1] = rc.haplotype[1][i];
			ArrayList<Integer> hitIndex = tf.findOneMismatch(queryLongSeq);
			if (hitIndex.isEmpty()) continue;
			Integer[] hitIndexArray = hitIndex.toArray(new Integer[hitIndex.size()]);
			if (hitIndexArray.length == 1) {
				TagAlignment cl = new TagAlignment(i, hitIndexArray[0]);
				tpList.add(cl);
				cnt1++;
			}
			else {
				int maxCountHitIndex = this.errorFilter2(rc, errorRate, i, hitIndexArray);
				if (maxCountHitIndex == -1) continue;
				TagAlignment cl = new TagAlignment(i, maxCountHitIndex);
				tpList.add(cl);
				cnt2++;
			}
		}
		tps = tpList.toArray(new TagAlignment[tpList.size()]);
		System.out.println(tps.length + " TagAlignments found. " + cnt1 + " are size = 1, " + cnt2 + " are error curated from size > 1");
    }

    public int errorFilter1 (ReadCounts rc, double errorRate, int queryIndex, Integer[] hitIndexArray) {//(sum of other count)/ < error, error is greater than 0.01 based on binomial
        int[] counts = new int[hitIndexArray.length+1];
		int queryCount = rc.getReadCount(queryIndex);
		counts[0] = queryCount;
		for (int i = 0; i < hitIndexArray.length; i++) {
			counts[i+1] = rc.getReadCount(hitIndexArray[i]);
		}
		Arrays.sort(counts);
		int maxIndex = -1;
		int hit = Arrays.binarySearch(counts, queryCount);
		if (hit == (counts.length - 1) || hit == (counts.length - 2)) {
			if (queryCount == counts[counts.length-3]) return -1;
			int otherCount = 0;
			for (int i = 0; i < counts.length - 2; i++) {
				otherCount += counts[i];
			}
			double pError = (double)otherCount/(double)(otherCount+queryCount);
			if (pError < errorRate) {
				int max = 0;
				for (int i = 0; i < hitIndexArray.length; i++) {
					int cnt = rc.getReadCount(hitIndexArray[i]);
					if (cnt > max) {
						max = cnt;
						maxIndex = hitIndexArray[i];
					}
				}
				return maxIndex;
			}
		}
		return maxIndex;
    }

	public int errorFilter2 (ReadCounts rc, double errorRate, int queryIndex, Integer[] hitIndexArray) {// each other count < erroRate, errorrate is 0.01 or 0.02
        int[] counts = new int[hitIndexArray.length+1];
		int queryCount = rc.getReadCount(queryIndex);
		counts[0] = queryCount;
		for (int i = 0; i < hitIndexArray.length; i++) {
			counts[i+1] = rc.getReadCount(hitIndexArray[i]);
		}
		Arrays.sort(counts);
		int maxIndex = -1;
		int hit = Arrays.binarySearch(counts, queryCount);
		if (hit == (counts.length - 1) || hit == (counts.length - 2)) {
			if (queryCount == counts[counts.length-3]) return -1;
			for (int i = 0; i < counts.length - 2; i++) {
				if ((double)counts[i]/(double)(counts[i]+queryCount) > errorRate) return -1;
			}
			int max = 0;
			for (int i = 0; i < hitIndexArray.length; i++) {
				int cnt = rc.getReadCount(hitIndexArray[i]);
				if (cnt > max) {
					max = cnt;
					maxIndex = hitIndexArray[i];
				}
			}
			
		}
		return maxIndex;
    }

	public void reciprocate () {
		for (int i = 0; i < tps.length; i++) {
			tps[i].swapQueryHit();
		}
		Arrays.sort(tps);
		ArrayList<TagAlignment> tpList = new ArrayList();
		TagAlignment currentTagPair = tps[0];
		int count = 1;
		for (int i = 1; i < tps.length; i++) {
			if (tps[i].equals(currentTagPair)) {
				count++;
			}
			else {
				if (count == 2) {
					tpList.add(currentTagPair);
				}
				currentTagPair = tps[i];
				count = 1;
			}
		}
		if (count == 2) tpList.add(currentTagPair);
		TagAlignment[] tpArray = tpList.toArray(new TagAlignment[tpList.size()]);
		tps = tpArray;
		System.out.println(tps.length + " reciprocal TagPairs found");
	}

	public void writeTagAlignment (ReadCounts rc, String outfileS) {//used to output the cluster, for data visilization
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outfileS), 65536);
			bw.write("queryIndex\thitOIndex\tqueryCount\thitCount\tclusterSize");
			bw.newLine();
			TagpairFinder tf = new TagpairFinder (rc);
			for (int i = 0; i < rc.getSize(); i++) {
				long[] queryLongSeq = new long[2];
				queryLongSeq[0] = rc.haplotype[0][i];
				queryLongSeq[1] = rc.haplotype[1][i];
				ArrayList<Integer> hitIndex = tf.findOneMismatch(queryLongSeq);
				if (hitIndex.isEmpty()) continue;
				Integer[] hitIndexArray = hitIndex.toArray(new Integer[hitIndex.size()]);
				for (int j = 0; j < hitIndexArray.length; j++) {
					bw.write(String.valueOf(i)+"\t"+String.valueOf(hitIndexArray[j])+"\t");
					bw.write(String.valueOf(rc.hapcount[i])+"\t"+String.valueOf(rc.hapcount[hitIndexArray[j]])+"\t");
					bw.write(String.valueOf(hitIndexArray.length+1));
					//bw.write("\t"+BaseEncoder.getSequenceFromLong(queryLongSeq));
					bw.newLine();
				}
			}
			bw.flush();
			bw.close();
			System.out.println("Initial tagAlignment is written");
		}
		catch (Exception e) {
			System.out.println(e.toString());
		}
	}

	public void writeTagPair (ReadCounts rc, String outfileS) {
		try {
			DataOutputStream  dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfileS), 65536));
			for (int i = 0; i < tps.length; i++) {
				tps[i].writeTagPair(rc, dos, 2*i);
			}
			dos.flush();
			dos.close();
		}
		catch (Exception e) {
			System.out.println("Error occurred while writing TagPair file " + outfileS + " " + e.toString());
		}
	}

	public void screenPrintTps () {
		int total = 1000;
		for (int i = 0; i < total; i++) {
			System.out.println(tps[i].queryIndex + "\t" + tps[i].hitIndex);
		}
	}

	public void screenPrint (ReadCounts rc) {
		int total = 10;
		if (tps.length < total) total = tps.length;
		for (int i = 0; i < total; i++) {
			long[] lseq = new long[2];
			lseq[0] = rc.haplotype[0][tps[i].queryIndex];
			lseq[1] = rc.haplotype[1][tps[i].queryIndex];
			System.out.println(BaseEncoder.getSequenceFromLong(lseq));
			lseq[0] = rc.haplotype[0][tps[i].hitIndex];
			lseq[1] = rc.haplotype[1][tps[i].hitIndex];
			System.out.println(BaseEncoder.getSequenceFromLong(lseq));
			System.out.println();
		}
	}

    class TagAlignment implements Comparable <TagAlignment> {
        int queryIndex;
        int hitIndex;

        TagAlignment (int queryIndex, int hitIndex) {
            this.queryIndex = queryIndex;
            this.hitIndex = hitIndex;
        }

		void swapQueryHit () {
			if (queryIndex > hitIndex) {
				int mid = queryIndex;
				queryIndex = hitIndex;
				hitIndex = mid;
			}
		}

		void writeTagPair (ReadCounts rc, DataOutputStream dos, int order) {
			try {
				dos.writeLong(rc.haplotype[0][queryIndex]);
				dos.writeLong(rc.haplotype[1][queryIndex]);
				dos.writeInt(order);
				dos.writeLong(rc.haplotype[0][hitIndex]);
				dos.writeLong(rc.haplotype[1][hitIndex]);
				dos.writeInt(order+1);
			}
			catch (Exception e) {
				System.out.println(e.toString());
			}
		}

		public boolean equals (TagAlignment o) {
			if (queryIndex == o.queryIndex && hitIndex == o.hitIndex) return true;
			return false;
		}

		public int compareTo(TagAlignment o) {
			if (queryIndex < o.queryIndex) return -1;
			else if (queryIndex > o.queryIndex) return 1;
			else {
				if (hitIndex < o.hitIndex) return -1;
				else if (hitIndex > o.hitIndex) return 1;
				else return 0;
			}
		}

    }
}
