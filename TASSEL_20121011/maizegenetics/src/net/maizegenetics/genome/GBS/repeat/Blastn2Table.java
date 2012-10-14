/*
convert BLAST+ blastn result format to table format
 */

package net.maizegenetics.genome.GBS.repeat;

import java.io.*;

/**
 *
 * @author fl262
 */
public class Blastn2Table {
	String[] header = {"Query_ID", "Query_Lenth", "Subject_ID", "Subject_Lenth", "Score",
						"E-Value", "Identities", "Identities_P", "Gap", "Gap_P",
						"Frame", "Query_Origin", "Query_End", "Subject_Origin", "Subject_End"};
	static String sourceFileS = processGlmSAS.blastReS;
	static String desFileS = processGlmSAS.blastDirS + "oout.txt";
	
	public Blastn2Table (String infile, String outfile) {
		readAndWriteFile(new File(infile), new File(outfile));
	}


	public void readAndWriteFile(File sourceFile, File desFile) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(sourceFile), 65536);
			BufferedWriter bw = new BufferedWriter(new FileWriter(desFile), 65536);
			bw.write(getLine(header));
			bw.newLine();
			String temp = br.readLine();
			String[] item = new String[header.length];
			while (temp != null) {
				if (temp.startsWith("Query=")) {
					item[0] = temp.substring(8);
					temp = br.readLine();
					while (!temp.startsWith("Length")) {
						item[0] = item[0] + temp;
						temp = br.readLine();
					}
					item[1] = temp.substring(7);
				}
				else if (temp.startsWith(">")) {
					item[2] = temp.substring(2);
					temp = br.readLine();
					while (!temp.startsWith("Length")) {
						item[2] = item[2] + temp;
						temp = br.readLine();
					}
					item[3] = temp.substring(7);
				}
				else if (temp.startsWith(" Score")) {
					String[] temps = temp.split(",");
					item[4] = temps[0].split(" +")[3];
					item[5] = temps[1].split(" +")[3];
					temp = br.readLine();
					temps = temp.split(",");
					String[] tem = temps[0].split(" +");
					item[6] = tem[3];
					item[7] = tem[4].replaceAll("\\(|\\)|%", "");
					tem = temps[1].split(" +");
					item[8] = tem[3];
					item[9] = tem[4].replaceAll("\\(|\\)|%", "");
					temp = br.readLine();
					item[10] = temp.substring(13);
				}
				else if (temp.startsWith("Query ")) {
					String[] temps = temp.split(" +");
					item[11] = temps[1];
					item[12] = temps[3];
					temp = br.readLine();
					temp = br.readLine();
					temps = temp.split(" +");
					item[13] = temps[1];
					item[14] = temps[3];
					temp = br.readLine();
					temp = br.readLine();
					while (temp.startsWith("Query ")) {
						temps = temp.split(" +");
						item[12] = temps[3];
						temp = br.readLine();
						temp = br.readLine();
						temps = temp.split(" +");
						item[14] = temps[3];
						temp = br.readLine();
						temp = br.readLine();
					}
					bw.write(getLine(item));
					bw.newLine();
				}
				temp = br.readLine();
			}
			bw.flush();
			bw.close();
			br.close();
		}
		catch (Exception e) {
			System.out.println("Exception occurred in " + sourceFile.toString() + " " + desFile.toString());
		}
		
	}
	public String getLine(String[] item) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < item.length - 1; i++) {
			sb.append(item[i]).append("\t");
		}
		sb.append(item[item.length-1]);
		return sb.toString();
	}

	public static void main(String[] args) {
		Blastn2Table b2t = new Blastn2Table(sourceFileS, desFileS);
	}
}
