/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS.switchgrass;

import java.io.File;
import java.util.ArrayList;
import java.util.Random;
import net.maizegenetics.genome.GBS.ReadsByTaxaPlasmid;
/**
 *
 * @author fl262
 */
public class FilterReadsByTaxa {
	ReadsByTaxaPlasmid inRbt;
	ReadsByTaxaPlasmid outRbt;
	int minCountOfAllTaxa = 6; //filter 0;
	int maxCountOfAllTaxa = 60;
	int minCountOfEachTaxa = 0;
	public FilterReadsByTaxa (String infile, String outfile, int filterType) {
		inRbt = new ReadsByTaxaPlasmid(infile, true);
		doFilter (outfile, filterType);
	}
	public void doFilter (String outfile, int filterType) {
		ArrayList<Integer> indexList = new ArrayList();
		for (int i = 0; i < inRbt.haplotypeNum; i++) {
			int[] count = inRbt.hapDist[i];
			if (passFilter(count, filterType, true)) {
				indexList.add(i);
			}
		}
		Integer[] index = indexList.toArray(new Integer[indexList.size()]);
		System.out.println (index.length);
		long[][] reads = new long[2][index.length];
		int[][] readDist = new int[index.length][inRbt.taxaNum];
		for (int i = 0; i < index.length; i++) {
			reads[0][i] = inRbt.haplotype[0][index[i]];
			reads[1][i] = inRbt.haplotype[1][index[i]];
			readDist[i] = inRbt.hapDist[index[i]];
		}
		String[] namesForTaxa = inRbt.taxaNames;
		outRbt = new ReadsByTaxaPlasmid(reads, readDist, namesForTaxa);
		outRbt.writeDistFile(new File(outfile), true, 5);
	}
	public boolean passFilter (int[] count, int filterType, boolean ifrandom) {
		if (filterType == 0) {
			int total = 0;
			for (int i = 0; i < count.length; i++ )	{
				total = total + count[i];
			}
			Random ran = new Random();
			if (total > minCountOfAllTaxa) {
				if (ifrandom) {
					if (ran.nextFloat() < 0.01) {
						return true;
					}
					else {
						return false;
					}	
				}
				return true;
			}
			return false;
		}
		else if (filterType == 1) {
			int total = 0;
			for (int i = 0; i < count.length; i++) {
				total = total + count[i];
			}
			Random ran = new Random();
			if (total > minCountOfAllTaxa && total < maxCountOfAllTaxa) {
				if (ifrandom) {
					if (ran.nextFloat() < 0.01) {
						return true;
					}
					else {
						return false;
					}
				}
				return true;
			}
			return false;
		}
		else {
			return false;
		}
	}
}
