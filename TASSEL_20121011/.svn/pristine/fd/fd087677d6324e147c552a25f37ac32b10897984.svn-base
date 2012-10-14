/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS.repeat;

/**
 *
 * @author fl262
 */
public class pipepineGBSRepeat {
	public static void main (String[] args) {
		mAndSGeneticMap masgm = new mAndSGeneticMap(); //select markers, set genotype value for them
		setQuery sq = new setQuery(); //select repetitive reads, set phenotype value for them

		//***************************//
		//*****Run SAS Glmselect*****//
		//***************************//


		processGlmSAS pgs = new processGlmSAS();//filtered SAS result and run BLAST againt GenBank nr database
		Blastn2Table b2t = new Blastn2Table(Blastn2Table.sourceFileS, Blastn2Table.desFileS);//convert Blastn format to table
		filterBlastRec fbr = new filterBlastRec(); //filter the blast table result by various filters
		annotateByKeyword abkw = new annotateByKeyword(); //select blast recs and Glm rec by keywords and draw pictures using visualizeGlm, also do the GUI
		visualizeGlm vg = new visualizeGlm();//draw a picture for Glmselect
	}
}
