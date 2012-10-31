/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.gbs.pav;

import java.io.File;
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByte;
import net.maizegenetics.genome.BaseEncoder;
import net.maizegenetics.gbs.util.OpenBitSet;

/**
 *
 * @author Fei Lu
 */
public class Test {

	public Test () {}

	public void makeSubTBTFromTagCount() {
		String tagCountFileS = "E:/Research/plastid/mito.cnt.txt";
		String inputTBTFileS = "Z:/tbt_20.byte";
		String outputTBTFileS = "M:/out.txt";
		TagCounts tc = new TagCounts (tagCountFileS, FilePacking.Text);
		PAVUtils uti = new PAVUtils();
		uti.writeSubTBTTxtByTagCount(inputTBTFileS, outputTBTFileS, tc);

	}

	public void makeSubTBTFromTagPGMap () {
		String tagPGMapFileS = "M:/pav/PhyGenMapping/pav.pgmap.txt";
		String inputTBTFileS = "Z:/tbt_20.byte";
		String outputTBTFileS = "M:/out.txt";
		TagPGMap pgmap = new TagPGMap (tagPGMapFileS, false);
		PAVUtils uti = new PAVUtils();
		uti.writeSubTBTTxtByTagPGMap(inputTBTFileS, outputTBTFileS, pgmap);
	}

	//**output first n row of TBT*/
	public void makeNrowOfTBT () {
		String inTBTFileS = "M:/tbt_test.byte";
		//String inTBTFileS = "N:/Zea/build20120110/mergedtbt/merged_chr1-4.tbt.byte";
		String outTBTFileS = "M:/pav/mergedTBT/tbt_test_50.byte";
		int nRow = 50;
		PAVUtils uti = new PAVUtils();
		uti.writeSubTBTByNumber(inTBTFileS, outTBTFileS, nRow);
	}

	public void convertTxtTBT () {
		String tbtInputFileS = "M:/test55againstAnchor/fake2.tbt.txt";
		String outFileS = "M:/test55againstAnchor/fake2.tbt.byte";
		TagsByTaxaByte tbt = new TagsByTaxaByte (tbtInputFileS,FilePacking.Text);
		tbt.writeDistFile(new File(outFileS), FilePacking.Byte, 0);
	}

	public void getTagCountOfTBT () {
		String tbtInputFileS = "N:/Zea/build20120110/mergedtbt/zea20120110c510.tbt.byte";
		TagsByTaxaByte tbt = new TagsByTaxaByte (tbtInputFileS,FilePacking.Byte);
	}

	//**output text format of TBT*/
	public void txtTBT () {
		String tbtInputFileS = "M:/pav/mergedTBT/tbt_test.byte";
		String outputFileS = "M:/tbt.txt";
		TagsByTaxaByte tbt = new TagsByTaxaByte (tbtInputFileS,FilePacking.Byte);
		tbt.writeDistFile(new File(outputFileS), FilePacking.Text, 1);

	}

	public void txtTagCount () {
		String inTagCount = "M:/RCG-51_C05F2ACXX_merged.cnt";
		String outTagCount = "M:/RCG-51_C05F2ACXX_merged.txt";
		TagCounts tc = new TagCounts (inTagCount, FilePacking.Byte);
		tc.writeTagCountFile(outTagCount, FilePacking.Text, 0);

	}

	//**calculate the minimum count of TBT*/
	public void getMinimumCountOfTBT () {
		String infileS = "M:/pav/tbtByLane/10225395_61VBRAAXX_s_1.tbt.byte";
		TagsByTaxaByte tbt = new TagsByTaxaByte (infileS,FilePacking.Byte);
		int min = 500;
		for (int i = 0; i < tbt.getTagCount(); i++) {
			if (tbt.getReadCount(i) < min) {
				min = tbt.getReadCount(i);
			}
		}
		System.out.println(min);
	}

	//**calculate the minimum count of TagCount*/
	public void getMinimumCountOfTagCount () {
		String infileS = "M:/pav/tags/merged.cnt";
		TagCounts tc = new TagCounts (infileS, FilePacking.Byte);
		int min = 50;
		for (int i = 0; i < tc.getSize(); i++) {
			if (tc.getReadCount(i) < min) {
				min = tc.getReadCount(i);
			}
		}
		System.out.println(min);
	}

	//**convert first n row of TBT to TagCount*/
	public void convertNrowTBTtoTagCount () {
		String tbtFileS = "M:/pav/mergedTBT/tbt_20_unmap_copy.byte";
		//String tbtFileS = "M:/tbts/434GFAAXX_s_2.tbt.byte";
		String tagCountFileS = "M:/merged_test.cnt";
		int nRow = 300000;
		PAVUtils uti = new PAVUtils ();
		uti.convertNrowTBTToTagCount(tbtFileS, tagCountFileS, nRow);
	}

	public void getCountByID () {
		String infileS = "M:/merged_20.cnt";
		TagCounts tc = new TagCounts (infileS, FilePacking.Byte);
		//int[] ID = {88635, 155467, 323427, 348021, 446194, 2968881, 3538578, 4533150, 4533151, 5566251, 7114118, 7114119, 7114120, 8748806, 8822261};
		int[] ID = {3538578, 4533150, 4533151, 5566251, 5574240, 8748806};
		for (int i = 0; i < ID.length; i++) {
			System.out.println (ID[i] + "\t" + BaseEncoder.getSequenceFromLong(tc.getTag(ID[i] -1)) + "\t" + tc.getReadCount(ID[i]-1));
		}
	}

	public void seqDiff () {
		long targetLong = BaseEncoder.getLongFromSeq("CAGAGATAAAAAAGAAGAAAAAAAACTTCGCT");
		long queryLong = BaseEncoder.getLongFromSeq("CAGCGATTAAAAAGAAGAAAAAAAAAAAAAAA");

            int c=BaseEncoder.seqDifferencesForSubset(targetLong,  queryLong, 8, 1);
            System.out.println(c);

	}

	public void testBits () {
		long[] a = new long[2];
		a[0] = 4;
		long[] b = new long[2];
		b[0] = 7;
		OpenBitSet sa = new OpenBitSet (a, 2);
		OpenBitSet sb = new OpenBitSet (b, 2);
		System.out.println(OpenBitSet.intersectionCount(sa, sb));
	}

	public void BitRoundTrip () {
		int nseq = 10000;
		long[][] tags = new long[2][nseq];
		for (int i = 0; i < tags.length; i++) {
			for (int j = 0; j < tags[0].length; j++) {
				tags[i][j] = Math.round(Math.random()*Double.MAX_VALUE);
			}
		}
		int noMatch = 0;
		for (int i = 0; i < tags[0].length; i++) {
			long[] temp = new long[2];
			for (int j = 0; j < temp.length; j++) {
				temp[j] = tags[j][0];
			}
			String seq = BaseEncoder.getSequenceFromLong(temp);
			long[] tem = BaseEncoder.getLongArrayFromSeq(seq);
			for (int j = 0; j < tem.length; j++) {
				if (tem[j] != temp[j]) {
					noMatch++;
					break;
				}
			}
		}
		System.out.println((double)noMatch/nseq);
	}

	public void BitRoundTrip2 () {
		int nseq = 10000;
		String[] seq = new String[nseq];
		for (int i = 0; i < nseq; i++) {
			seq[i] = "";
			for (int j = 0; j < 64; j++) {
				double r = Math.random();
				if (r < 0.25) seq[i] = seq[i] + "A";
				else if (r < 0.50) seq[i] = seq[i] + "T";
				else if (r < 0.75) seq[i] = seq[i] + "G";
				else if (r < 1.00) seq[i] = seq[i] + "C";
			}
		}
		int same = 0;
		for (int i = 0; i < seq.length; i++) {
			long[] temp = BaseEncoder.getLongArrayFromSeq(seq[i]);
			String t = BaseEncoder.getSequenceFromLong(temp);
			if (t.equals(seq[i])) {
				same++;
			}
			else {
				System.out.println("No!");
			}
		}
		System.out.println((double)same/nseq);
	}
}
