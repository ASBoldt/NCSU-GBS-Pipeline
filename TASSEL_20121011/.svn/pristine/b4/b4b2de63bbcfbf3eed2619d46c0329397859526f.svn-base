/*
 This program convert GenBank files in a folder to one fasta file and then make it as a blast library
 This request NCBI BLAST+ applications
 */
package net.maizegenetics.genome.GBS.repeat;

import java.io.*;

/**
 *
 * @author fl262
 */
public class GenBank2Fasta {

	static String sourceDirS = "E:/Database/GenBank";
	static File desFile = new File("E:/Database/GenBank/GenBank_Plant_nr.lib");

	public GenBank2Fasta() {
		desFile.delete();
		File sourceDir = new File(sourceDirS);
		File[] gbFileArray = sourceDir.listFiles();
		for (int i = 0; i < gbFileArray.length; i++) {
			if (gbFileArray[i].toString().endsWith("seq")) {
				readAndWrite(gbFileArray[i], desFile);
			}
		}
	}

	public void readAndWrite(File sourceFile, File desFile) {
		try {
			BufferedReader br = new BufferedReader (new FileReader(sourceFile), 65536);
			BufferedWriter bw = new BufferedWriter(new FileWriter(desFile, true), 65536);
			String temp = br.readLine();
			String def = "", acc = "", org = "";
			while (temp != null) {
				if (temp.startsWith("DEFINITION")) {
					def = "";
					acc = "";
					org = "";
					StringBuilder sb = new StringBuilder(temp.substring(12));
					temp = br.readLine();
					while(!temp.startsWith("ACCESSION")) {
						sb.append(" ").append(temp.substring(12));
						temp = br.readLine();
					}
					def = sb.toString();
					acc = temp.substring(12);
				}
				else if (temp.startsWith("  ORGANISM")) {
					org = temp.substring(12);
				}
				else if (temp.startsWith("ORIGIN")) {
					StringBuilder sb = new StringBuilder(">");
					sb.append(acc).append("|").append(def).append("|").append(org);
					bw.write(sb.toString());
					bw.newLine();
					temp = br.readLine();
					while(!temp.startsWith("//")) {
						temp = temp.replaceAll("\\d| ", "");
						bw.write(temp);
						bw.newLine();
						temp = br.readLine();
					}
				}
				temp = br.readLine();
			}
			bw.flush();
			bw.close();
			br.close();
		}
		catch (Exception e) {
			System.out.println(sourceFile.toString());
		}
	}

	public static void main(String args[]) throws IOException {
		GenBank2Fasta gb2f = new GenBank2Fasta();
		String cmd = "Makeblastdb -in " + desFile + " -dbtype nucl";
		Runtime.getRuntime().exec(cmd);
	}
}
