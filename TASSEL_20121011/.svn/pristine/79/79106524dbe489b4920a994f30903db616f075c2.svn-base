package net.maizegenetics.gwas.Imputation;

public class SnpInformation implements Comparable<SnpInformation> {
	public String name;
	public String allele;
	public int chromosome;
	public int position;
	public int ndx;
	public byte val;
	public byte[] values;
	public GenotypeData dataset;
	
	@Override
	public int compareTo(SnpInformation otherSnp) {
		if (chromosome > otherSnp.chromosome) return 1;
		if (chromosome < otherSnp.chromosome) return -1;
		if (position > otherSnp.position) return 1;
		if (position < otherSnp.position) return -1;
		if (ndx > otherSnp.ndx) return 1;
		if (ndx < otherSnp.ndx) return -1;
		return name.compareTo(otherSnp.name);
	}
}
