/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.gbs.pipeline;

/**
 *
 * @author Fei Lu
 */
public class UNEAK {
	String workingDirS;

	public UNEAK (String workingDirS) {
		this.workingDirS = workingDirS;
		this.pipe();
	}

	public void pipe () {
		this.creatDir(workingDirS);
		this.creatTagCounts(workingDirS, false);
		this.mergeTaxaAndTagCount(workingDirS);
		this.getTagPairs(workingDirS);
		this.creatTBT(workingDirS);
		this.creatMapInfo();
		this.creatHapMap();
	}

	public void creatHapMap () {
		UMapInfoToHapMapPlugin umithm = new UMapInfoToHapMapPlugin();
		String arguments = "-w " + workingDirS + " -mnMAF 0.05 -mxMAF 0.5 -mnC 0 -mxC 1";
		String[] args = arguments.split(" ");
		umithm.setParameters(args);
		umithm.performFunction(null);
	}

	public void creatMapInfo () {
		UTBTToMapInfoPlugin uttmi = new UTBTToMapInfoPlugin();
		String arguments = "-w " + workingDirS;
		String[] args = arguments.split(" ");
		uttmi.setParameters(args);
		uttmi.performFunction(null);
	}

	public void creatTBT (String workingDirS) {
		UTagPairToTBTPlugin utptt = new UTagPairToTBTPlugin();
		String arguments = "-w " + workingDirS;
		String[] args = arguments.split(" ");
		utptt.setParameters(args);
		utptt.performFunction(null);
	}

	public void getTagPairs (String workingDirS) {
		UTagCountToTagPairPlugin utpbn = new UTagCountToTagPairPlugin ();
		String arguments = "-w " + workingDirS + " -e 0.03";
		String[] args = arguments.split(" ");
		utpbn.setParameters(args);
		utpbn.performFunction(null);
	}

	public void mergeTaxaAndTagCount (String workingDirS) {
		UMergeTaxaTagCountPlugin umttc = new UMergeTaxaTagCountPlugin ();
		String arguments = "-w " + workingDirS + " -t y -c 1";
		String[] args = arguments.split(" ");
		umttc.setParameters(args);
		umttc.performFunction(null);
	}

	public void creatTagCounts (String workingDirS, boolean ifQseq) {
		String arguments = "-w " + workingDirS + " -e PstI";
		String[] args = arguments.split(" ");
		if (ifQseq) {
			UQseqToTagCountPlugin uqtc = new UQseqToTagCountPlugin();
			uqtc.setParameters(args);
			uqtc.performFunction(null);
		}
		else {
			UFastqToTagCountPlugin uftc = new UFastqToTagCountPlugin();
			uftc.setParameters(args);
			uftc.performFunction(null);
		}
	}

	public void creatDir (String workingDirS) {
		UCreatWorkingDirPlugin cwd = new UCreatWorkingDirPlugin();
		String arguments = "-w " + workingDirS;
		String[] args = arguments.split(" ");
		cwd.setParameters(args);
		cwd.performFunction(null);
	}

	public static void main (String[] args) {
		String workingDirS = "M:/UNEAK/";
		new UNEAK (workingDirS);
	}
}
