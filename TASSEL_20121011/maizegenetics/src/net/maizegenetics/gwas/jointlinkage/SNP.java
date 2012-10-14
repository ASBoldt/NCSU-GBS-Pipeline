package net.maizegenetics.gwas.jointlinkage;

public class SNP implements Comparable<SNP> {
	public float[] score;
	public byte[] bytescore;
	public final int chr;
	public final float geneticPos;
	public final String name;
	public final int physicalPos;
	public final String allele;
	
	public int index = 0;
	
	public SNP(int chr, float geneticPosition, float[] score, int index) {
		this("", "", chr, geneticPosition, 0, score, index);
	}
	
	public SNP(int chr, float geneticPosition, float[] score, String name, int index) {
		this(name, "", chr, geneticPosition, 0, score, index);
	}
	
	public SNP(String name, String allele, int chr, float geneticPosition, int physicalPosition, float[] score, int index) {
		this.chr = chr;
		geneticPos = geneticPosition;
		this.score = score;
		this.name = name;
		this.index = index;
		physicalPos = physicalPosition;
		this.allele = allele;
	}

	public float getFloatScore(int index) {
		if (score != null) return score[index];
		if (bytescore != null) return bytescore[index];
		return Float.NaN;
	}
	
	public double[] getScoresAsDoubles() {
		int n = score.length;
		double[] result = new double[n];
		for (int i = 0; i < n; i++) result[i] = score[i];
		return result;
	}
	
	@Override
	public int compareTo(SNP snp) {
		if (chr > snp.chr) return 1;
		if (chr < snp.chr) return -1;
		if (geneticPos > snp.geneticPos) return 1;
		if (geneticPos < snp.geneticPos) return -1;
		if (physicalPos > snp.physicalPos) return 1;
		if (physicalPos < snp.physicalPos) return -1;
		int comp = allele.compareTo(snp.allele);
		if (comp != 0) return comp;
		return name.compareTo(snp.name);
	}

	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof SNP)) return false;
		if (compareTo((SNP) obj) == 0) return true;
		return false;
	}

	@Override
	public int hashCode() {
		return name.hashCode();
	}

}
