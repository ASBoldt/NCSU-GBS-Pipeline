/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS.switchgrass;

/**
 *
 * @author fl262
 */
public class switchRun {
	public switchRun (int type) {
		if (type == 0) {
			step0 ();
		}
		else if (type == 1) {
			step1 ();
		}
		else if (type == 2) {
			step2 ();
		}
		else if (type == 3) {
			step3 ();
		}
	}
	public void step0 () {
		String infile = "M:/GBS/ReadsByTaxaMin/ReadsByTaxaMin.txt";
		String outfile = "M:/switch/FilteredReadByTaxa/ReadsByTaxaFil.txt";
		new FilterReadsByTaxa (infile, outfile, 1);
	}
	public void step1 () {
		String infile = "M:/switch/FilteredReadByTaxa/ReadsByTaxaFil.txt";
		String outfile = "M:/switch/FilteredReadByTaxa/nonbirary.txt";
		new mkTxtFile(infile, outfile);
	}
	public void step2 () {
		String infile = "M:/GBS/ReadsByTaxaMin/ReadsByTaxaMin.txt";
		String outfile = "M:/list.txt";
		new linePersent (infile, outfile);
	}
	public void step3 () {
		String infile = "M:/switch/ReadByTaxaAll/ReadsByTaxa.txt";
		String outfile = "M:/switch/ReadByTaxaBit/TagByTaxaBitAll.txt";
		new intToBit (infile, outfile);
	}
	public static void main (String[] args) {
		new switchRun(3);
	}
}
