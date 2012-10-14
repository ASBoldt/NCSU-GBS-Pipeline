/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS.plasmid;

import java.io.*;

/**
 *
 * @author fl262
 */
public class showFile {
	String parentDir = plastidPipe.parentDir.replaceFirst("/$", "_Visual/");
	String[] childDir = plastidPipe.childDir.clone();
	int totalLine = 3000;
	public showFile () {
		setupDirectory();
		for (int i = 0; i < childDir.length; i++) {
				String dirS = childDir[i];
				String sourceDirS = plastidPipe.parentDir + childDir[i] + "/";
				String desDirS = parentDir + childDir[i] + "/";
				if (dirS.equals("Qseq")) {
					copyPasteDir(sourceDirS, desDirS);
				} else if (dirS.equals("Key")) {
					copyPasteDir(sourceDirS, desDirS);
				} else if (dirS.equals("ParsedTaxaDir")) {
					formatDir(sourceDirS, desDirS);
				} else if (dirS.equals("CollapsedTags")) {
					formatDir(sourceDirS, desDirS);
				} else if (dirS.equals("CombinedTags")) {
					formatDir(sourceDirS, desDirS);
				} else if (dirS.equals("CombinedRepTags")) {
					formatDir(sourceDirS, desDirS);
				} else if (dirS.equals("TBTtxt")) {
					copyPasteDir(sourceDirS, desDirS);
				} else {
					System.out.println("There must be missing methods for source directories");
					System.exit(1);
				}
			}
	}
	public void formatDir (String sourceDirS, String desDirS) {
		String[] files = new File(sourceDirS).list();
		for (int i = 0; i < files.length; i++) {
			String sourceFile = sourceDirS + files[i];
			String desFile = desDirS + files[i];
			format(sourceFile, desFile);
		}
		System.gc();
	}
	public void format (String sourceFileS, String desFileS) {
		File fn = new File(sourceFileS);
		int rows = (int) (fn.length() / (baseEncoder.byteLen+4));
		int total;
		if (rows < totalLine) {
			total = rows;
		} else {
			total = totalLine;
		}
		try {
			DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(sourceFileS), 65536));
			PrintWriter pw = new PrintWriter(desFileS);
			for (int i = 0; i < total; i++) {
				byte[] tag = new byte[baseEncoder.byteLen];
				dis.read(tag, 0, baseEncoder.byteLen);
				String seq = baseEncoder.getSeqFromByteArray(tag);
				int counting = dis.readInt();
				pw.println(seq + "\t" + counting);
			}
			pw.flush();
			pw.close();
			dis.close();
		} catch (Exception e) {
			System.out.println("Reading or writing error: " + sourceFileS + " / " + desFileS);
		}
	}
	public void copyPasteDir (String sourceDirS, String desDirS) {
		String[] files = new File(sourceDirS).list();
		for (int i = 0; i < files.length; i++) {
			String sourceFile = sourceDirS + files[i];
			String desFile = desDirS + files[i];
			simpleCopyPaste(sourceFile, desFile);
		}
		System.gc();
	}
	private void simpleCopyPaste(String sourceFileS, String desFileS) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(sourceFileS), 65536);
			BufferedWriter bw = new BufferedWriter(new FileWriter(desFileS), 65536);
			String strLine;
			for (int i = 0; i < totalLine; i++) {
				strLine = br.readLine();
				if (strLine == null) {
					break;
				} else {
					bw.write(strLine);
					bw.newLine();
				}
			}
			bw.flush();
			bw.close();
			br.close();
		} catch (Exception e) {
			System.out.println("Reading or writing error: " + sourceFileS + " / " + desFileS);
		}
	}
	private void forReadsCounts(String sourceFile, String desFile) {
		
	}
	private void setupDirectory() {
		File[] FileOfData = new File[childDir.length];
		for (int i = 0; i < childDir.length; i++) {
			String path = parentDir + childDir[i] + "/";
			FileOfData[i] = new File(path);
			FileOfData[i].mkdirs();
		}
	}
	public static void main (String[] args) {
		new showFile ();
	}
}
