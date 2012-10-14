/*
record mapping result and counting number for a repetitive read
 */
package net.maizegenetics.genome.GBS.repeat;

/**
 *
 * @author fl262
 */
public class ReadWCountPosition implements Comparable <ReadWCountPosition> {
	long[] haplotype = new long[2];
	int totalCount;
	byte countInTaxa[];
	float normalTotalCount;
	float normalCountInTaxa[];
	byte chrB;
	byte strand;
	int positionMin;
	int positionMax;
	short nextCutDistance;
	byte divergence;
	float dcoP;
	float mapP;
	byte multimaps;

	public ReadWCountPosition() {}

	
	public void setLongSeq(long[] haplotype) {
		for (int i = 0; i < haplotype.length; i++) {
			this.haplotype[i] = haplotype[i];
		}
	}
	public void setTotalCount (int totalCount) {
		this.totalCount = totalCount;
	}
	public void setCountInTaxa (byte[] countInTaxa) {

		this.countInTaxa = countInTaxa.clone();
	}
	public void setNormalTotalCount(float normalTotalCount) {
		this.normalTotalCount = normalTotalCount;
	}

	public void setNormalCountInTaxa (float[] normalCountInTaxa) {
		this.normalCountInTaxa = normalCountInTaxa.clone();
	}
	

	public void setPosition(byte chrB, byte strand, int positionMin, int positionMax, short nextCutDistance, byte divergence, float dcoP, float mapP, byte multimaps) {
		this.chrB = chrB;
		this.strand = strand;
		this.positionMin = positionMin;
		this.positionMax = positionMax;
		this.nextCutDistance = nextCutDistance;
		this.divergence = divergence;
		this.dcoP = dcoP;
		this.mapP = mapP;
		this.multimaps = multimaps;
	}

	public String setTitle (int j) {
		StringBuilder sb = new StringBuilder (">");
		sb = sb.append(j).append("|").append(totalCount).append("|").append(chrB).append("|").append(strand).append("|").append(positionMin).append("|").append(positionMax).append("|").append(nextCutDistance).append("|").append(divergence).append("|").append(dcoP).append("|").append(mapP).append("|").append(multimaps);
		return sb.toString();
	}

	public int compareTo (ReadWCountPosition o) {
		if (totalCount < o.totalCount) return -1;
		else if (totalCount > o.totalCount) return 1;
		else {
			if (chrB < o.chrB) return -1;
			else if (chrB > o.chrB) return 1;
			else {
				if (positionMin < o.positionMin) return -1;
				else if (positionMin > o.positionMin)  return 1;
				else return 0;
			}
		}
	}
}
