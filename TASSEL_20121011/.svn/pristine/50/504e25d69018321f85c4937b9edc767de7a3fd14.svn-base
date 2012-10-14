/*
 filter the Blast record table
 */

package net.maizegenetics.genome.GBS.repeat;

import java.io.File;
import java.util.ArrayList;

/**
 *
 * @author fl262
 */
public class filterBlastRec {
	static String sourceFileS = Blastn2Table.desFileS;
	static String desFilesS = processGlmSAS.blastDirS + "filtered_oout.txt";
	float efilter = (float)1e-15;
	blastRec[] br;
	blastRecArray brad;
	public filterBlastRec () {
		getAndFilterData();
		writeData();
	}
	public void getAndFilterData() {
		blastRecArray bras = new blastRecArray (new File(sourceFileS));
		ArrayList<blastRec> filteredRec = new ArrayList();
		for (int i = 0; i < bras.recNum; i++) {
			if (Float.valueOf(bras.recs[i].item[5]) > efilter) {
				continue;
			}
			if (bras.recs[i].item[2].contains("clone")) {
				continue;
			}
			if (bras.recs[i].item[2].contains("locusiuiiuh")) {
				continue;
			}
			if (bras.recs[i].item[2].contains("chromosome")) {
				continue;
			}
			filteredRec.add(bras.recs[i]);
		}
		br = filteredRec.toArray(new blastRec[filteredRec.size()]);
		brad = new blastRecArray(br, bras.header);
	}
	public void writeData () {
		brad.writeArray(desFilesS);
	}
	public static void main (String[] args) {
		new filterBlastRec ();
	}
}
