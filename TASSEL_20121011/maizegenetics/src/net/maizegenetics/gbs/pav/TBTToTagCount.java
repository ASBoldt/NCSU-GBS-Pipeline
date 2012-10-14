/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.gbs.pav;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import net.maizegenetics.gbs.tagdist.TagCountMutable;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByte;

/**
 *
 * @author Fei Lu
 */
public class TBTToTagCount {

	public TBTToTagCount (File tbtFolder, File tagCountFolder, File tempFolder) {
		//this.deleteTagCountFiles(tempFolder);
		this.parse(tbtFolder, tempFolder);
		this.convert(tempFolder, tagCountFolder, 2);
	}

	private void deleteTagCountFiles (File tagCountFolder) {
		File[] fs = tagCountFolder.listFiles();
		for (int i = 0; i < fs.length; i++) {
			fs[i].delete();
		}
	}

	public void parse (File tbtFolder, File tempFolder) {
		tempFolder.mkdir();
		File[] tbtLaneFiles = tbtFolder.listFiles();
		for (int i = 0; i < tbtLaneFiles.length; i++) {
			TagsByTaxaByte tbt = new TagsByTaxaByte (tbtLaneFiles[i].getAbsolutePath(),FilePacking.Byte);
			String[] outfileSs = new String[tbt.getTaxaCount()];
			ParseOut[] pos = new ParseOut[outfileSs.length];
			for (int j = 0; j < outfileSs.length; j++) {
				String name = tbt.getTaxaName(j);
				name = name.replaceAll("_", "-");
				name = name.replaceAll(":", "_");
				Pattern p = Pattern.compile("[^a-zA-Z0-9()_]");
				Matcher m = p.matcher(name);
				name = m.replaceAll("-");
				outfileSs[j] = tempFolder.getAbsolutePath() + "/" +name+".tem";
				pos[j] = new ParseOut (outfileSs[j]);
			}
			for (int j = 0; j < tbt.getTagCount(); j++) {
				long[] seq = tbt.getTag(j);
				byte[] dis = tbt.getTaxaReadCountsForTag(j);
				for (int k = 0; k < outfileSs.length; k++) {
					if (dis[k] == 0) continue;
					pos[k].writeTempFile(seq, (byte)tbt.getTagLength(j), dis[k]);
				}
			}
			for (int j = 0; j < outfileSs.length; j++) {
				pos[j].closeOut();
			}
			System.gc();
		}
	}

	public void convert (File tempFolder, File tagCountFolder, int tagLengthInLong) {
		File[] tempFiles = tempFolder.listFiles();
		for (int i = 0; i < tempFiles.length; i++) {
			String outfileS = tagCountFolder.getAbsolutePath() + "/" + tempFiles[i].getName().replaceFirst("tem", "cnt");
			System.out.println(outfileS);
			ParseIn pi = new ParseIn (tempFiles[i], tagLengthInLong);
			TagCountMutable tcm = pi.getTagCountMutable();
			tcm.writeTagCountFile(outfileS, FilePacking.Byte, 1);
			tempFiles[i].delete();
		}
		tempFolder.delete();
	}

	private class ParseIn {
		DataInputStream dis;
		int tagNum;
		int tagLengthInLong;
		ParseIn (File infile, int tagLengthInLong) {
			try {
				this.dis = new DataInputStream (new BufferedInputStream (new FileInputStream(infile), 65536));
			}
			catch (Exception e) {
				System.out.println(e.toString());
			}
			tagNum = (int)infile.length()/(8*tagLengthInLong+2);
			this.tagLengthInLong = tagLengthInLong;
		}

		TagCountMutable getTagCountMutable () {
			TagCountMutable tc = new TagCountMutable (tagLengthInLong, tagNum);
			try {
				long[] tag = new long[tagLengthInLong];
				for (int i = 0; i < tagNum; i++) {
					for (int j = 0; j < tagLengthInLong; j++) {
						tag[j] = dis.readLong();
					}
					byte tagLength = dis.readByte();
					tc.addReadCount(tag, tagLength, dis.readByte());
				}
				dis.close();
			}
			catch (Exception e) {
				System.out.println(e.toString());
			}
			tc.collapseCounts();
			return tc;
		}
	}

	private class ParseOut {
		DataOutputStream dos;

		ParseOut (String outfileS) {
			try {
				this.dos = new DataOutputStream (new BufferedOutputStream (new FileOutputStream(outfileS), 65536));
			}
			catch (Exception e) {
				System.out.println(e.toString());
			}
		}

		void writeTempFile (long[] seq, byte tagLength, byte count) {
			try {
				for (int i = 0; i < seq.length; i++) {
                    dos.writeLong(seq[i]);
                }
				dos.writeByte(tagLength);
				dos.writeByte(count);
			}
			catch (Exception e) {
				System.out.println(e.toString());
			}
		}

		void closeOut () {
			try {
				dos.flush();
				dos.close();
			}
			catch (Exception e) {
				System.out.println(e.toString());
			}
		}
	}
}
